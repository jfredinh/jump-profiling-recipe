#!/usr/bin/env python
"""
Process JUMP profiles: join AP/MAP scores, filter plates, and convert to AnnData.

Usage:
    pixi run python process_profiles.py
"""

import pandas as pd
import anndata as ad
import numpy as np
from pathlib import Path


def load_and_join_ap(profiles_path: str, ap_path: str, map_path: str) -> pd.DataFrame:
    """Load profiles and join with average precision (AP) and mean average precision (MAP) scores."""
    print("Loading profiles...")
    profiles = pd.read_parquet(profiles_path)
    print(f"  Profiles shape: {profiles.shape}")

    # --- Join AP scores (per-well) ---
    print("\nLoading AP scores...")
    ap = pd.read_parquet(ap_path)
    print(f"  AP shape: {ap.shape}")

    # Join keys and AP-specific columns
    join_keys = ["Metadata_Source", "Metadata_Plate", "Metadata_Well"]
    ap_cols = ["n_pos_pairs", "n_total_pairs", "average_precision"]

    # Select and rename AP columns with Metadata_ prefix
    ap_subset = ap[join_keys + ap_cols].copy()
    ap_subset = ap_subset.rename(columns={
        "n_pos_pairs": "Metadata_n_pos_pairs",
        "n_total_pairs": "Metadata_n_total_pairs",
        "average_precision": "Metadata_average_precision"
    })

    # Left join AP
    print("Joining AP scores...")
    merged = profiles.merge(ap_subset, on=join_keys, how="left")
    print(f"  Merged shape: {merged.shape}")
    print(f"  Rows with AP scores: {merged['Metadata_average_precision'].notna().sum()}")
    print(f"  Rows without AP scores: {merged['Metadata_average_precision'].isna().sum()}")

    # --- Join MAP scores (per-perturbation) ---
    print("\nLoading MAP scores...")
    map_df = pd.read_parquet(map_path)
    print(f"  MAP shape: {map_df.shape}")

    # MAP is indexed by Metadata_JCP2022, reset index to make it a column
    map_df = map_df.reset_index()

    # Select and rename MAP columns with Metadata_ prefix (exclude 'indices' column)
    map_cols = ["mean_average_precision", "p_value", "corrected_p_value", "below_p", "below_corrected_p"]
    map_subset = map_df[["Metadata_JCP2022"] + map_cols].copy()
    map_subset = map_subset.rename(columns={
        "mean_average_precision": "Metadata_mean_average_precision",
        "p_value": "Metadata_map_p_value",
        "corrected_p_value": "Metadata_map_corrected_p_value",
        "below_p": "Metadata_map_below_p",
        "below_corrected_p": "Metadata_map_below_corrected_p"
    })

    # Left join MAP
    print("Joining MAP scores...")
    merged = merged.merge(map_subset, on="Metadata_JCP2022", how="left")
    print(f"  Merged shape: {merged.shape}")
    print(f"  Rows with MAP scores: {merged['Metadata_mean_average_precision'].notna().sum()}")
    print(f"  Rows without MAP scores: {merged['Metadata_mean_average_precision'].isna().sum()}")

    return merged


def filter_plates(df: pd.DataFrame, exclude_plate_types: list[str]) -> pd.DataFrame:
    """Filter out specified plate types."""
    print(f"\nFiltering out plate types: {exclude_plate_types}")
    print(f"  PlateTypes before: {df['Metadata_PlateType'].unique().tolist()}")

    mask = ~df["Metadata_PlateType"].isin(exclude_plate_types)
    filtered = df[mask].copy()

    print(f"  PlateTypes after: {filtered['Metadata_PlateType'].unique().tolist()}")
    print(f"  Rows removed: {len(df) - len(filtered)}")
    print(f"  Rows remaining: {len(filtered)}")

    return filtered


def convert_to_anndata(df: pd.DataFrame) -> ad.AnnData:
    """Convert DataFrame to AnnData format."""
    print("\nConverting to AnnData...")

    # Separate metadata and feature columns
    meta_cols = [c for c in df.columns if c.startswith("Metadata_")]
    feat_cols = [c for c in df.columns if not c.startswith("Metadata_")]

    print(f"  Metadata columns: {len(meta_cols)}")
    print(f"  Feature columns: {len(feat_cols)}")

    # Create AnnData object
    X = df[feat_cols].values.astype(np.float32)
    obs = df[meta_cols].copy()
    obs.index = obs.index.astype(str)

    # Convert boolean/mixed columns to strings for h5ad compatibility
    for col in obs.columns:
        if obs[col].dtype == bool or obs[col].dtype == 'boolean':
            obs[col] = obs[col].astype(str)
        # Handle object columns that contain booleans (mixed with NaN)
        elif obs[col].dtype == 'object':
            sample = obs[col].dropna().head(1)
            if len(sample) > 0 and isinstance(sample.iloc[0], (bool, np.bool_)):
                obs[col] = obs[col].map(lambda x: str(x) if pd.notna(x) else "NaN")

    var = pd.DataFrame(index=feat_cols)
    var["feature_name"] = feat_cols

    adata = ad.AnnData(X=X, obs=obs, var=var)
    print(f"  AnnData shape: {adata.shape}")

    return adata


def main():
    # Configuration
    input_dir = Path("outputs/orf-crispr-PT")
    profiles_path = input_dir / "profiles_var_mad_int_featselect_harmonyrsc_nc300_pca300_reconstructed.parquet"
    ap_path = input_dir / "metrics/profiles_var_mad_int_featselect_harmonyrsc_nc300_pca300_ap_negcon.parquet"
    map_path = input_dir / "metrics/profiles_var_mad_int_featselect_harmonyrsc_nc300_pca300_map_negcon.parquet"
    exclude_plate_types = ["TARGET2"]

    # Output paths
    output_parquet = input_dir / "profiles_var_mad_int_featselect_harmonyrsc_nc300_pca300_reconstructed_with_metrics_filtered.parquet"
    output_h5ad = input_dir / "profiles_var_mad_int_featselect_harmonyrsc_nc300_pca300_reconstructed_with_metrics_filtered.h5ad"

    print("=" * 60)
    print("JUMP Profiles Processing Pipeline")
    print("=" * 60)

    # Step 1: Load and join AP/MAP scores
    df = load_and_join_ap(str(profiles_path), str(ap_path), str(map_path))

    # Step 2: Filter out unwanted plate types
    df = filter_plates(df, exclude_plate_types)

    # Step 3: Save filtered parquet
    print(f"\nSaving parquet: {output_parquet}")
    df.to_parquet(output_parquet)

    # Step 4: Convert to AnnData and save
    adata = convert_to_anndata(df)
    print(f"Saving h5ad: {output_h5ad}")
    adata.write_h5ad(output_h5ad)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Final shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"Plate types: {df['Metadata_PlateType'].value_counts().to_dict()}")
    print(f"\nOutputs:")
    print(f"  Parquet: {output_parquet}")
    print(f"  H5AD:    {output_h5ad}")


if __name__ == "__main__":
    main()
