"""
Functions for Rapids SingleCell Harmony batch correction
"""

import logging
import numpy as np
import pandas as pd
from ..preprocessing import io

# Optional imports for Rapids SingleCell Harmony
try:
    import anndata
    import rapids_singlecell as rsc
    RAPIDS_AVAILABLE = True
except ImportError:
    RAPIDS_AVAILABLE = False

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Configuration
N_PCA_COMPONENTS = 300


def apply_harmonyrsc_correction(dframe_path, batch_key, n_clusters, output_path):
    """Perform Rapids SingleCell Harmony batch correction on feature data.

    Parameters
    ----------
    dframe_path : str
        Path to the input parquet file containing metadata and features.
    batch_key : str
        Column name in metadata that identifies the batch information.
    n_clusters : int
        Number of clusters for Harmony correction (default: 300).
    output_path : str
        Path where the corrected data will be saved as a parquet file.

    Returns
    -------
    None
        The corrected data is saved to the specified output path.

    Notes
    -----
    The function performs the following steps:
    1. Loads input parquet and splits into metadata and features
    2. Handles special case for high-cardinality batch keys like Metadata_Well
    3. Creates AnnData object and performs PCA reduction
    4. Applies Rapids SingleCell Harmony correction
    5. Saves harmonized features with original metadata to output path

    This implementation is significantly faster than the standard harmonypy
    implementation, typically completing in minutes instead of hours.

    Requires rapids-singlecell and anndata to be installed.
    """

    if not RAPIDS_AVAILABLE:
        raise ImportError(
            "Rapids SingleCell Harmony requires 'rapids-singlecell' and 'anndata' packages. "
            "Install them with: pixi add rapids-singlecell anndata"
        )

    logger.info(f"Loading data from {dframe_path}")
    df = pd.read_parquet(dframe_path)

    # Handle high-cardinality batch keys like Metadata_Well
    if batch_key == "Metadata_Well":
        logger.warning("Running harmony on Metadata_Well may not be meaningful due to high cardinality.")
        logger.info("Grouping wells by plate size to reduce cardinality")
        # Group the wells into groups based on plate size to remove the chance of fitting incorrect trends
        # Separate 384 and 1536 well plates, i.e. source 1 and 9
        batch_key = "Metadata_Well_Grouped"
        df[batch_key] = df.apply(
            lambda row: f"{row['Metadata_Well']}_1536"
            if row['Metadata_Source'] in ['source_1', 'source_9']
            else f"{row['Metadata_Well']}_384",
            axis=1
        )

    # Split features and metadata
    meta = df[[c for c in df.columns if c.startswith("Metadata_")]].copy()

    # Convert string columns to categorical for memory efficiency
    for col in meta.columns:
        if meta[col].dtype == 'object':
            meta[col] = meta[col].astype('category')

    # Create meaningful index from Source, Plate, and Well columns
    meta.index = (meta['Metadata_Source'].astype(str) + ':' +
                  meta['Metadata_Plate'].astype(str) + ':' +
                  meta['Metadata_Well'].astype(str))

    # Extract feature data
    feats = df[[c for c in df.columns if not c.startswith("Metadata_")]].values

    # Handle NaN values
    feats = np.nan_to_num(feats, nan=0.0)

    logger.info(f"Original shape: {feats.shape}")

    # Create AnnData object with explicit string index
    adata = anndata.AnnData(X=feats, obs=meta)

    logger.info(f"Running PCA with {N_PCA_COMPONENTS} components")
    # Run PCA first
    rsc.tl.pca(adata, n_comps=N_PCA_COMPONENTS)
    logger.info(f"PCA complete: {adata.obsm['X_pca'].shape}")

    logger.info(f"Running Rapids SingleCell Harmony on key: {batch_key} with {n_clusters} clusters")

    # Run harmony on PCA coordinates (use GEMM to avoid CUDA alignment issues)
    rsc.pp.harmony_integrate(
        adata,
        key=batch_key,
        use_gemm=True,
        n_clusters=int(n_clusters),
        max_iter_harmony=20
    )

    logger.info(f"Harmony complete: {adata.obsm['X_pca_harmony'].shape}")

    # Create DataFrame with harmonized features
    harmony_features = adata.obsm['X_pca_harmony']
    harmony_cols = [f"harmony_{i+1}" for i in range(harmony_features.shape[1])]
    harmony_df = pd.DataFrame(harmony_features, columns=harmony_cols, index=adata.obs.index)

    # Add original metadata columns
    for col in adata.obs.columns:
        harmony_df[col] = adata.obs[col].values

    # Save to parquet
    logger.info(f"Saving harmonized data to {output_path}")
    harmony_df.to_parquet(output_path)
    logger.info(f"Rapids SingleCell Harmony correction complete. Output saved to {output_path}")