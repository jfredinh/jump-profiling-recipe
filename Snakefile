wildcard_constraints:
    pipeline=r"[_a-zA-Z.~0-9\-]*",
    scenario=r"[_a-zA-Z.~0-9\-]*",


import jump_profiling_recipe.correct as correct
import jump_profiling_recipe.preprocessing as pp


include: "rules/sphering.smk"
include: "rules/map.smk"


rule all:
    input:
        f"outputs/{config['scenario']}/reformat.done",
        # Include harmonyrsc sweep and summary if config specifies the parameter lists
        f"outputs/{config['scenario']}/harmonyrsc_sweep.done" if config.get("harmonyrsc_n_clusters_list") else [],
        f"outputs/{config['scenario']}/metrics/harmonyrsc_sweep_summary.csv" if config.get("harmonyrsc_n_clusters_list") else [],
        # ap_negcon_path=f"outputs/{config['scenario']}/metrics/{config['pipeline']}_ap_negcon.parquet",
        # map_negcon_path=f"outputs/{config['scenario']}/metrics/{config['pipeline']}_map_negcon.parquet",
        # ap_nonrep_path=f"outputs/{config['scenario']}/metrics/{config['pipeline']}_ap_nonrep.parquet",
        # map_nonrep_path=f"outputs/{config['scenario']}/metrics/{config['pipeline']}_map_nonrep.parquet",

rule reformat:
    input:
        f"outputs/{config['scenario']}/{config['pipeline']}.parquet",
    output:
        touch("outputs/{scenario}/reformat.done"),
    params:
        profile_dir=lambda w: f"outputs/{w.scenario}/",
        meta_col_new=config.get("meta_col_new", None),
    run:
        if "meta_col_new" in config:
            correct.format_check.run_format_check(params.profile_dir, params.meta_col_new)
        else:
            correct.format_check.run_format_check(params.profile_dir)

rule write_parquet:
    output:
        "outputs/{scenario}/profiles.parquet",
    params:
        existing_profile_file=config.get("existing_profile_file", None),
    benchmark:
        "benchmarks/{scenario}/write_parquet.txt"
    run:
        if "existing_profile_file" in config:
            shell("mkdir -p $(dirname {output}) && cp {input} {output}".format(input=params.existing_profile_file, output=output))
        else:
            pp.io.write_parquet(
                config["sources"],
                config["plate_types"],
                output[0],
                profile_type=config.get("profile_type"),
                search_additional_metadata=config.get("search_additional_metadata", False)
            )


rule compute_norm_stats:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/norm_stats/{pipeline}.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_normstats.txt"
    params:
        use_negcon=config["use_mad_negcon"],
    run:
        pp.stats.compute_norm_stats(*input, *output, **params)


rule select_variant_feats:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
        "outputs/{scenario}/norm_stats/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_var.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_var.txt"
    run:
        pp.stats.select_variant_features(*input, *output)


rule mad_normalize:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
        "outputs/{scenario}/norm_stats/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_mad.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_mad.txt"
    run:
        pp.normalize.mad(*input, *output)


rule INT:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_int.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_int.txt"
    run:
        pp.transform.rank_int(*input, *output)


rule well_correct:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_wellpos.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_wellpos.txt"
    run:
        correct.corrections.subtract_well_mean(*input, *output)


rule cc_regress:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_cc.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_cc.txt"
    params:
        cc_path=config.get("cc_path"),
    run:
        correct.corrections.regress_out_cell_counts_parallel(
            *input, *output, params.cc_path
        )


rule remove_outliers:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_outlier.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_outlier.txt"
    run:
        pp.clean.remove_outliers(*input, *output)


rule drop_na_rows:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_dropna.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_dropna.txt"
    params:
        na_threshold=config.get("na_threshold", 0.1),
        max_rows_to_drop=config.get("max_rows_to_drop", 100),
    run:
        correct.corrections.remove_na_rows(
            *input,
            *output,
            na_threshold=params.na_threshold,
            max_rows_to_drop=params.max_rows_to_drop
        )


rule annotate_genes:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_annotated.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_annotated.txt"
    params:
        df_gene_path="inputs/metadata/crispr.csv.gz",
        df_chrom_path="inputs/metadata/gene_chromosome_map.tsv",
    run:
        correct.corrections.annotate_dataframe(
            *input, *output, params.df_gene_path, params.df_chrom_path
        )


rule pca_transform:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_PCA.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_PCA.txt"
    run:
        correct.corrections.transform_data(*input, *output)


rule correct_arm:
    input:
        "outputs/{scenario}/{pipeline}_annotated.parquet",
    output:
        "outputs/{scenario}/{pipeline}_corrected.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_corrected.txt"
    params:
        gene_expression_path="inputs/metadata/Recursion_U2OS_expression_data.csv.gz",
    run:
        correct.corrections.arm_correction(*input, *output, params.gene_expression_path)


rule featselect:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_featselect.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_featselect.txt"
    params:
        keep_image_features=config["keep_image_features"],
    run:
        pp.select_features(*input, *output, *params)


rule harmony:
    input:
        "outputs/{scenario}/{pipeline}.parquet",
    output:
        "outputs/{scenario}/{pipeline}_harmony.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_harmony.txt"
    params:
        batch_key=config["batch_key"],
        use_gpu=config.get("use_gpu", False),
    run:
        correct.apply_harmony_correction(input[0], params.batch_key, params.use_gpu, output[0])


rule setup_harmonyrsc:
    output:
        directory("resources/harmonyrsc"),
        "resources/harmonyrsc/harmonyrsc.py",
        "resources/harmonyrsc/pyproject.toml"
    params:
        repo_url="https://github.com/shntnu/harmonyrsc.git",
        commit_hash="main"
    shell:
        """
        git clone {params.repo_url} resources/harmonyrsc
        cd resources/harmonyrsc
        git checkout {params.commit_hash}
        pixi install
        """


rule harmonyrsc:
    input:
        profiles="outputs/{scenario}/{pipeline}.parquet",
        harmonyrsc_setup="resources/harmonyrsc/harmonyrsc.py"
    output:
        "outputs/{scenario}/{pipeline}_harmonyrsc_nc{n_clusters}_pca{n_pca}.parquet",
    benchmark:
        "benchmarks/{scenario}/{pipeline}_harmonyrsc_nc{n_clusters}_pca{n_pca}.txt"
    wildcard_constraints:
        n_clusters=r"\d+",
        n_pca=r"\d+"
    params:
        batch_key=config["batch_key"],
        harmonyrsc_dir="resources/harmonyrsc"
    shell:
        """
        cd {params.harmonyrsc_dir} && \
        PIXI_PROJECT_MANIFEST="" pixi run -e default python harmonyrsc.py \
        ../../{input.profiles} ../../{output} {params.batch_key} {wildcards.n_clusters} {wildcards.n_pca}
        """


# Rule to run all harmonyrsc parameter combinations with MAP calculations
rule harmonyrsc_sweep:
    input:
        profiles=expand(
            "outputs/{scenario}/{pipeline}_harmonyrsc_nc{n_clusters}_pca{n_pca}.parquet",
            scenario=config["scenario"],
            pipeline=config.get("harmonyrsc_input_pipeline", "profiles_var_mad_int_featselect"),
            n_clusters=config.get("harmonyrsc_n_clusters_list", [100, 200, 300]),
            n_pca=config.get("harmonyrsc_n_pca_list", [50, 100, 200, 300])
        ),
        maps=expand(
            "outputs/{scenario}/metrics/{pipeline}_harmonyrsc_nc{n_clusters}_pca{n_pca}_map_negcon.parquet",
            scenario=config["scenario"],
            pipeline=config.get("harmonyrsc_input_pipeline", "profiles_var_mad_int_featselect"),
            n_clusters=config.get("harmonyrsc_n_clusters_list", [100, 200, 300]),
            n_pca=config.get("harmonyrsc_n_pca_list", [50, 100, 200, 300])
        )
    output:
        touch("outputs/{scenario}/harmonyrsc_sweep.done")


# Rule to create summary heatmap of MAP scores across parameter combinations
rule harmonyrsc_sweep_summary:
    input:
        maps=expand(
            "outputs/{scenario}/metrics/{pipeline}_harmonyrsc_nc{n_clusters}_pca{n_pca}_map_negcon.parquet",
            scenario=config["scenario"],
            pipeline=config.get("harmonyrsc_input_pipeline", "profiles_var_mad_int_featselect"),
            n_clusters=config.get("harmonyrsc_n_clusters_list", [100, 200, 300]),
            n_pca=config.get("harmonyrsc_n_pca_list", [50, 100, 200, 300])
        ),
        sweep_done="outputs/{scenario}/harmonyrsc_sweep.done"
    output:
        heatmap_map="outputs/{scenario}/metrics/harmonyrsc_sweep_heatmap_map.png",
        heatmap_pct="outputs/{scenario}/metrics/harmonyrsc_sweep_heatmap_pct_significant.png",
        heatmap_pct_fdr10="outputs/{scenario}/metrics/harmonyrsc_sweep_heatmap_pct_significant_fdr10.png",
        summary_csv="outputs/{scenario}/metrics/harmonyrsc_sweep_summary.csv"
    params:
        metrics_dir=lambda w: f"outputs/{w.scenario}/metrics",
        threshold=config.get("map_params", {}).get("threshold", 0.05),
        harmonyrsc_dir="resources/harmonyrsc"
    shell:
        """
        cd {params.harmonyrsc_dir} && \
        PIXI_PROJECT_MANIFEST="" pixi run -e default python -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import glob
import os

# Use absolute path
base_dir = os.path.dirname(os.path.dirname(os.getcwd()))
metrics_dir = os.path.join(base_dir, '{params.metrics_dir}')
threshold = {params.threshold}

print(f'Looking for files in: {{metrics_dir}}')

# Find all MAP files
map_files = glob.glob(os.path.join(metrics_dir, '*_harmonyrsc_nc*_pca*_map_negcon.parquet'))
print(f'Found {{len(map_files)}} MAP files')

results = []
for map_file in map_files:
    match = re.search(r'_nc([0-9]+)_pca([0-9]+)_map_negcon[.]parquet', map_file)
    if match:
        n_clusters = int(match.group(1))
        n_pca = int(match.group(2))

        map_df = pd.read_parquet(map_file)
        mean_norm_map = map_df['mean_average_precision'].mean()
        n_significant = (map_df['corrected_p_value'] < threshold).sum()
        n_significant_fdr10 = (map_df['corrected_p_value'] < 0.1).sum()
        n_total = len(map_df)
        pct_significant = (n_significant / n_total) * 100 if n_total > 0 else 0
        pct_significant_fdr10 = (n_significant_fdr10 / n_total) * 100 if n_total > 0 else 0

        results.append({{
            'n_clusters': n_clusters,
            'n_pca': n_pca,
            'mean_normalized_map': mean_norm_map,
            'n_significant': n_significant,
            'n_significant_fdr10': n_significant_fdr10,
            'n_total': n_total,
            'pct_significant': pct_significant,
            'pct_significant_fdr10': pct_significant_fdr10,
            'p_value_threshold': threshold
        }})

print(f'Processed {{len(results)}} results')

summary_df = pd.DataFrame(results)
summary_df = summary_df.sort_values(['n_pca', 'n_clusters'])
output_csv = os.path.join(base_dir, '{output.summary_csv}')
summary_df.to_csv(output_csv, index=False)
print(f'Saved summary to: {{output_csv}}')

# Heatmap for normalized MAP
heatmap_map_data = summary_df.pivot(index='n_clusters', columns='n_pca', values='mean_normalized_map')
heatmap_map_data = heatmap_map_data.sort_index(ascending=True)
heatmap_map_data = heatmap_map_data.reindex(sorted(heatmap_map_data.columns), axis=1)

plt.figure(figsize=(12, 10))
vmin = heatmap_map_data.values.min()
vmax = heatmap_map_data.values.max()
sns.heatmap(heatmap_map_data, annot=True, fmt='.4f', cmap='viridis',
            vmin=vmin, vmax=vmax, cbar_kws={{'label': 'Mean Normalized MAP'}})
plt.title('HarmonyRSC Parameter Sweep: Mean Normalized MAP')
plt.xlabel('Number of PCA Components')
plt.ylabel('Number of Clusters')
plt.tight_layout()
output_heatmap_map = os.path.join(base_dir, '{output.heatmap_map}')
plt.savefig(output_heatmap_map, dpi=150)
plt.close()
print(f'Saved MAP heatmap to: {{output_heatmap_map}}')

# Heatmap for % significant
heatmap_pct_data = summary_df.pivot(index='n_clusters', columns='n_pca', values='pct_significant')
heatmap_pct_data = heatmap_pct_data.sort_index(ascending=True)
heatmap_pct_data = heatmap_pct_data.reindex(sorted(heatmap_pct_data.columns), axis=1)

plt.figure(figsize=(12, 10))
vmin = heatmap_pct_data.values.min()
vmax = heatmap_pct_data.values.max()
sns.heatmap(heatmap_pct_data, annot=True, fmt='.1f', cmap='viridis',
            vmin=vmin, vmax=vmax, cbar_kws={{'label': f'% Significant (p < {{threshold}})'}})
plt.title(f'HarmonyRSC Parameter Sweep: % Compounds with Corrected p-value < {{threshold}}')
plt.xlabel('Number of PCA Components')
plt.ylabel('Number of Clusters')
plt.tight_layout()
output_heatmap_pct = os.path.join(base_dir, '{output.heatmap_pct}')
plt.savefig(output_heatmap_pct, dpi=150)
plt.close()
print(f'Saved pct significant heatmap to: {{output_heatmap_pct}}')

# Heatmap for % significant with FDR 0.1
heatmap_pct_fdr10_data = summary_df.pivot(index='n_clusters', columns='n_pca', values='pct_significant_fdr10')
heatmap_pct_fdr10_data = heatmap_pct_fdr10_data.sort_index(ascending=True)
heatmap_pct_fdr10_data = heatmap_pct_fdr10_data.reindex(sorted(heatmap_pct_fdr10_data.columns), axis=1)

plt.figure(figsize=(12, 10))
vmin = heatmap_pct_fdr10_data.values.min()
vmax = heatmap_pct_fdr10_data.values.max()
sns.heatmap(heatmap_pct_fdr10_data, annot=True, fmt='.1f', cmap='viridis',
            vmin=vmin, vmax=vmax, cbar_kws={{'label': '% Significant (FDR < 0.1)'}})
plt.title('HarmonyRSC Parameter Sweep: % Compounds with Corrected p-value < 0.1')
plt.xlabel('Number of PCA Components')
plt.ylabel('Number of Clusters')
plt.tight_layout()
output_heatmap_pct_fdr10 = os.path.join(base_dir, '{output.heatmap_pct_fdr10}')
plt.savefig(output_heatmap_pct_fdr10, dpi=150)
plt.close()
print(f'Saved FDR 0.1 heatmap to: {{output_heatmap_pct_fdr10}}')
"
        """
