# nf-index
 
Nextflow (DSL2) pipelines that build a **consensus index of regulatory elements** from per-sample calls. They project every sample onto that index as element × sample matrices, normalize the signal, and run follow-up analyses (variance partition, core sets, differential accessibility).
 
We have used it to build two kinds of index:
- **DNase I hypersensitive sites (DHSs)**, from [hotspot3](https://github.com/vierstralab/hotspot3) peak calls;
- **consensus footprints**, from per-sample footprint calls.
The index-building step works on any per-sample interval calls in the same BED layout (see [The samples file](#the-samples-file)). The rest of this README says "DHS" for an index element, because that is how the code and output files name them (`dhs_id`, `autosomal_dhs`, `masterlist_DHSs_*`).
 
This is the second stage of the processing pipeline. Per-sample inputs (filtered CRAMs and [hotspot3](https://github.com/vierstralab/hotspot3) peaks, density, fit stats and p-values) come from [`nf-aligning`](https://github.com/vierstralab/nf-aligning).
 
```
per-sample calls (DHS peaks or footprints)
   └─► build_masterlist.nf ─► index.anndata.zarr            (consensus index + binary presence matrix)
          └─► generate_matrices.nf ─► index+matrices.anndata.zarr   (counts, binary, density, p-values, background)
                 └─► normalize_signal.nf ─► index+matrices+normalized.anndata.zarr  (lowess + DESeq2 VST + variance partition)
                        ├─► core_sets.nf             (group-specific "core" DHS sets)
                        └─► differential_deseq.nf    (DESeq2 differential tests)
```
 
Each stage reads the AnnData (zarr) written by the previous one, so run them in order from the **same launch directory**. `outdir` defaults to `$launchDir/output`, and later stages look for `output/zarr/...` there.
 
## Workflows at a glance
 
Every workflow and entry point, what it needs, and what it writes. Run each with `nextflow run <script> [-entry <entry>] -profile <profile> -resume [--param value ...]`. Details are in [Workflows](#workflows); the samples-file columns are described in [The samples file](#the-samples-file).
 
| Script / entry | Purpose | Inputs it needs | Main output |
|---|---|---|---|
| `build_masterlist.nf` | Build the consensus index | `--samples_file` with `sample_id` + peak/footprint-file column (`--index_peaks_column`, default `peaks_file_0.001fdr`); params `chrom_sizes`, `genome_fasta`, `mappable_not_blacklist_file`, `gencode`, `repeats` | `zarr/index.anndata.zarr` |
| `build_masterlist.nf -entry annotateMasterlist` | Re-annotate an existing index | `--index_dir` (a previous build's `outdir`, containing `samples_order.txt` and `unfiltered_masterlists/`); same annotation params | `masterlist_DHSs_all_chunks.Altius.annotated.bed` |
| `generate_matrices.nf` | Element × sample matrices | `--index_anndata` (default `output/zarr/index.anndata.zarr`); `--samples_file` with `sample_id`, peak column (`--matrix_peaks_column`, default `peaks_file_0.01fdr`), `cram_file` [`cram_index`], `normalized_density_bw`, `hotspot3_fit_stats_file`, `hotspot3_pvals_parquet`, plus any covariates needed for normalization; params `genome_fasta`, `nuclear_chrom_sizes` | `zarr/index+matrices.anndata.zarr` |
| `generate_matrices.nf -entry extractDensityAtSummit` | Density matrix only | `--index_anndata`; `--samples_file` with `sample_id`, `normalized_density_bw` | `raw_matrices/matrix.density.npy` |
| `normalize_signal.nf` | Lowess + DESeq2 VST + variance partition | `--matrices_anndata` (default `output/zarr/index+matrices.anndata.zarr`) with `binary` and `--normalization_layer` layers, and `obs` columns used by `--vst_design_formula` and `--variance_partition_formula`; params `dhs_mask_name`, `chunk_size` | `index+matrices+normalized.anndata.zarr` |
| `normalize_signal.nf -entry existingModel` | Normalize new samples with a previous model | Same as above, plus `--template_anndata` (a previously normalized anndata) | `index+matrices+normalized.anndata.zarr` |
| `normalize_signal.nf -entry existingModelToScaleFactors` | Lowess scale factors only, from a previous model | `--matrices_anndata`, `--template_anndata`, `--normalization_layer` | `index+matrices+normalized.anndata.zarr` (scale factors only) |
| `variance_partition.nf` | Variance partition only | `--norm_matrix` (`.npy`), `--samples_order`, `--masterlist`, `--samples_meta`; params `variance_partition_formula`, `chunk_size` | `variance_partition.vp_annotated.sorted.bed` |
| `core_sets.nf` | Group-specific core sets | `--matrices_anndata` with `binary` and `neglog10_pvals` layers; `--grouping_column` (an `obs` column); params `core_set_correction_method`, `core_set_n_shuffles` | `core_sets/` |
| `differential_deseq.nf -entry diffDeseq` | Differential tests, chunked | `--dds` (DESeq2 dataset `.RDS`), `--total_chunks` | Concatenated result tables in `outdir` |
| `differential_deseq.nf -entry diffDeseqSpecies` | Differential tests per dataset | `--deseq_meta` (TSV with columns `prefix`, `dataset_anndata`) | `deseq_results.all_cell_types.tsv` |
 
Typical full run (from one launch directory):
 
```bash
nextflow run build_masterlist.nf  -profile Altius,new_cluster -resume --samples_file samples.tsv
nextflow run generate_matrices.nf -profile Altius,new_cluster -resume --samples_file samples.tsv
nextflow run normalize_signal.nf  -profile Altius,new_cluster -resume
nextflow run core_sets.nf         -profile Altius,new_cluster -resume --grouping_column <column>
```
 
## Requirements
 
- [Nextflow](https://www.nextflow.io/) (the conda env pins 23.10.1)
- [conda/mamba](https://conda.io): every process runs in the conda env given by `params.conda`
- A SLURM cluster or another executor configured as a profile in `nextflow.config`
Create the environment once:
 
```bash
mamba env create -n super-index -f environment.yml
```
 
Then set `conda` in `params.config` to that environment's path (for example `~/miniconda3/envs/super-index`). The environment installs `genome-tools` and `nmf-tools` from GitHub with pip.
 
Some processes call command-line tools that **are not in `environment.yml`** and must also be on `PATH` inside that env:
 
- `hotspot3-pvals`, from [hotspot3](https://github.com/vierstralab/hotspot3) (`pip install hotspot3`), used by `generate_matrices.nf` to extract per-DHS p-values from each sample's `pvals.parquet` (marked experimental in hotspot3).
- `normalize-matrix` and `qc-normalization`, used by `normalize_signal.nf`. Install the package that provides them into the same env.
- `bin/index_scripts/writeRowsPerChunkForMatrices` is a prebuilt binary. If it doesn't run on your system, rebuild it from `writeRowsPerChunkForMatrices.cpp` (`g++ -O3 -o writeRowsPerChunkForMatrices writeRowsPerChunkForMatrices.cpp`).
## Configuration
 
- `nextflow.config` sets the executor profiles and includes `params.config`.
  - `Altius`: SLURM, `conda.enabled = true`, and the memory labels `medmem`, `highmem`, `bigmem` (50 CPUs, 1.9 TB on the `bigmem` queue) and `bg_params`.
  - `new_cluster` / `old_cluster`: only override the queue. Combine them with a base profile, e.g. `-profile Altius,new_cluster`.
  - To run elsewhere, add your own profile with the same labels.
- `params.config` holds the parameters for **human GRCh38** (identical to `params.human.config`). `params.mouse.config` has the mm10 equivalents. Copy it over `params.config` or override individual params on the command line.
Parameters can be set in `params.config` or on the command line (`--name value`).
 
### Common parameters
 
| Parameter | Default | Used by | Description |
|---|---|---|---|
| `conda` | path to `super-index` env | all | Conda env used by every process |
| `outdir` | `$launchDir/output` | all | Output directory |
| `samples_file` | — | all | Tab-separated sample metadata, one row per sample, with a `sample_id` column |
| `chrom_sizes`, `genome_fasta` | GRCh38 no-alts | build, generate | Reference genome |
| `nuclear_chrom_sizes` | GRCh38 nuclear | generate | Passed to `hotspot3-pvals` |
| `mappable_not_blacklist_file` | K36 mappable minus blacklist | build | Used for GC-content annotation |
| `gencode`, `repeats` | GENCODE v49 GFF3, RepeatMasker | build | DHS annotation |
| `index_peaks_column` | `peaks_file_0.001fdr` | build | Column of `samples_file` holding the peak files used to **build** the index |
| `matrix_peaks_column` | `peaks_file_0.01fdr` | generate | Column holding the peak files used for the **binary matrix** |
| `index_anndata` | `${outdir}/zarr/index.anndata.zarr` | generate, offsets | Index produced by `build_masterlist.nf` |
| `matrices_anndata` | `${outdir}/zarr/index+matrices.anndata.zarr` | normalize, core sets | Output of `generate_matrices.nf` |
| `normalization_layer` | `counts` | normalize | Layer to normalize (`counts` or `bg_corrected_agg_cutcounts`) |
| `dhs_mask_name` | `autosomal_dhs` | normalize | Boolean `var` column selecting which DHSs are normalized |
| `vst_design_formula` | `~ 0 + interaction(core_annotation, extended_annotation, drop = TRUE)` | normalize | DESeq2 design for VST. The named columns must exist in sample metadata |
| `variance_partition_formula` | `~ nuclear_dup_rate + SPOT3_score + nuclear_reads + alignment_quality` | normalize | variancePartition formula. The named columns must exist in sample metadata |
| `chunk_size` | `5000` | normalize | DHSs per variance-partition job |
| `core_set_correction_method` | `cauchy` | core sets | P-value combination: `cauchy` or `stouffer` |
| `core_set_n_shuffles` | `1000` | core sets | Number of shuffles |
| `grouping_column` | *(none)* | core sets | **Required**: sample-metadata column that defines the groups |
 
## The samples file
 
A single tab-separated file with a header, one row per sample. Different workflows read different columns. Extra columns are allowed, and they are carried into the AnnData `obs` table by `generate_matrices.nf`. That table is where the normalization formulas look for their covariates.
 
| Column | Needed by | Content (as produced by nf-aligning) |
|---|---|---|
| `sample_id` | all | Unique sample ID. Must match the 4th column of that sample's peak file ([hotspot3](https://github.com/vierstralab/hotspot3) writes the sample ID there) |
| `peaks_file_0.001fdr` (or `--index_peaks_column`) | build | `<sample_id>/fdr0.001/<sample_id>.peaks.fdr0.001.bed.gz` |
| `peaks_file_0.01fdr` (or `--matrix_peaks_column`) | generate | `<sample_id>/fdr0.01/<sample_id>.peaks.fdr0.01.bed.gz` |
| `cram_file` | generate | `<sample_id>/<sample_id>.filtered.cram` |
| `cram_index` | generate (optional) | Defaults to `<cram_file>.crai` |
| `normalized_density_bw` | generate | `<sample_id>/<sample_id>.normalized_density.bw` |
| `hotspot3_fit_stats_file` | generate | `<sample_id>/<sample_id>.fit_stats.bed.gz` (tabix-indexed; older [hotspot3](https://github.com/vierstralab/hotspot3) versions name it `.fit_stats.tsv.gz`) |
| `hotspot3_pvals_parquet` | generate | `<sample_id>/<sample_id>.pvals.parquet` |
| covariates (e.g. `core_annotation`, `extended_annotation`, `nuclear_dup_rate`, `SPOT3_score`, `nuclear_reads`, `alignment_quality`) | normalize | Whatever columns `vst_design_formula` and `variance_partition_formula` use. `nuclear_reads` is also required when `normalization_layer = bg_corrected_agg_cutcounts` |
| your grouping column | core sets | Any categorical column passed as `--grouping_column` |
 
Example (paths shortened):
 
```tsv
sample_id	peaks_file_0.001fdr	peaks_file_0.01fdr	cram_file	normalized_density_bw	hotspot3_fit_stats_file	hotspot3_pvals_parquet
AG81965	out/AG81965/fdr0.001/AG81965.peaks.fdr0.001.bed.gz	out/AG81965/fdr0.01/AG81965.peaks.fdr0.01.bed.gz	out/AG81965/AG81965.filtered.cram	out/AG81965/AG81965.normalized_density.bw	out/AG81965/AG81965.fit_stats.bed.gz	out/AG81965/AG81965.pvals.parquet
```
 
Peak (or footprint) files may be gzipped or plain. Only their first 7 columns are used, in this order: `chrom, start, end, sample_id, score, density_summit, wavelet_summit`. hotspot3 peak files already have this layout. Footprint calls need to be written the same way, with the `sample_id` in column 4.
 
## Workflows
 
### 1. `build_masterlist.nf`: build the DHS index
 
```bash
nextflow run build_masterlist.nf -profile Altius,new_cluster -resume \
  --samples_file samples.tsv \
  --index_peaks_column peaks_file_0.001fdr
```
 
Builds the index with the approach of [Meuleman et al., 2020](https://www.nature.com/articles/s41586-020-2559-3), using the R implementation by W. Meuleman and E. Rynes in `bin/index_scripts/`:
 
1. **collate_and_chunk**: concatenates all peak files, sorts them, and splits the genome into independent chunks (`chunk_bed.awk`).
2. **process_chunk**: builds DHSs within each chunk (`code_build.R`).
3. **resolve_overlaps**: produces the `all`, `nonovl_core` and `nonovl_any` DHS sets per chunk (`code_overlap.R`).
4. **merge_chunks**: concatenates the chunks into genome-wide masterlists (`code_gen_masterlist.sh`).
5. **write_rows → collect_chunks → convert_to_numpy**: builds the binary DHS × sample presence matrix.
6. **annotate_masterlist**: adds GC content, GENCODE and RepeatMasker annotations.
7. **convert_index_to_anndata**: writes the AnnData.
Outputs in `outdir`:
 
- `samples_order.txt`: sample order (the row order of `samples_file`) used for all matrix columns
- `unfiltered_masterlists/masterlist_DHSs_Altius_*`: masterlists (`all`, `nonovl_core`, `nonovl_any`) as `_chunkIDs.bed`, `.bed12`, `.bb`
- `masterlist_DHSs_all_chunks.Altius.annotated.bed`: annotated masterlist
- `raw_matrices/binary.index.raw.matrix.{mtx.gz,npy}`
- `zarr/index.anndata.zarr`: AnnData with `obs` = samples, `var` = DHSs (with an `autosomal_dhs` mask and annotations), `X` = binary peak matrix, plus `obs['n_peaks']`
Extra entry: `-entry annotateMasterlist --index_dir <outdir>` re-runs only the annotation step on an existing build.
 
### 2. `generate_matrices.nf`: DHS × sample matrices
 
```bash
nextflow run generate_matrices.nf -profile Altius,new_cluster -resume \
  --samples_file samples.tsv \
  --matrix_peaks_column peaks_file_0.01fdr
```
 
Reads the index from `--index_anndata` (default `output/zarr/index.anndata.zarr`, so run it from the same launch directory as step 1). For each sample it builds one column of each matrix:
 
| Layer | Process | How |
|---|---|---|
| `binary` | `project_peak_calls` | DHS overlaps a peak by ≥ 50% of either (`bedmap --fraction-either 0.5`) |
| `counts` | `count_tags` | `featureCounts` reads per DHS (`-p` for paired data). CRAMs are converted to BAM with `genome_fasta` |
| `neglog10_pvals` | `extract_max_pval` | `hotspot3-pvals` from the p-value parquet |
| `density` | `extract_max_density` | Max normalized density at the DHS summit |
| `mean_bg_agg_cutcounts` | `extract_bg_mean` | [hotspot3](https://github.com/vierstralab/hotspot3) background fit at the DHS summit, from the fit-stats file |
 
The columns are assembled into `raw_matrices/matrix.<layer>.npy` and added as layers. `samples_file` is merged into `obs`; values in `samples_file` win when a column already exists.
 
Output: `zarr/index+matrices.anndata.zarr`.
 
The `samples_file` here can list the same samples as the index, or a different set (for example new samples projected onto an existing index). If the sample IDs differ from the index, `obs` and `X` are rebuilt from `samples_file`.
 
Extra entry: `-entry extractDensityAtSummit` builds only the density matrix, from `normalized_density_bw`.
 
### 3. `normalize_signal.nf`: normalization
 
```bash
nextflow run normalize_signal.nf -profile Altius,new_cluster -resume \
  --matrices_anndata output/zarr/index+matrices.anndata.zarr
```
 
For the DHSs selected by `dhs_mask_name`:
 
1. **normalize_matrix**: lowess normalization of `normalization_layer` against a reference (`normalize-matrix --normalization_type lowess`), plus QC plots.
2. **prepare_deseq_dataset → deseq2_vst**: DESeq2 dataset with lowess scale factors, then a variance-stabilizing transform with `vst_design_formula`.
3. **variancePartition**: `variancePartition` on the VST matrix, in chunks of `chunk_size` DHSs.
Outputs:
 
- `index+matrices+normalized.anndata.zarr`: new layers (VST matrix and lowess scale factors), lowess params and the VST dispersion function stored in `uns`/`varm`, variance-partition results in `varm`, and the formulas recorded in `uns`
- `normalization/`, `params/`, `qc/`, `deseq2/`, `<dhs_mask_name>.vp_annotated.sorted.bed`
Applying an existing model to new samples:
 
- `-entry existingModel --template_anndata <normalized anndata>`: re-uses the lowess params and VST dispersion function from a previous run, then runs VST and variance partition.
- `-entry existingModelToScaleFactors --template_anndata <...>`: lowess scale factors only.
### 4. `variance_partition.nf` (standalone)
 
Normally called from `normalize_signal.nf`. To run it on its own:
 
```bash
nextflow run variance_partition.nf -profile Altius -resume \
  --norm_matrix <vst.npy> --samples_order <samples_order.txt> \
  --masterlist <masterlist.bed> --samples_meta <samples_meta.tsv>
```
 
### 5. `core_sets.nf`: group-specific core DHS sets
 
```bash
nextflow run core_sets.nf -profile Altius,new_cluster -resume \
  --matrices_anndata output/zarr/index+matrices.anndata.zarr \
  --grouping_column <metadata column>
```
 
For each value of `grouping_column`, combines the per-sample `neglog10_pvals` of that group's samples (`cauchy` or `stouffer`, with `core_set_n_shuffles` shuffles) and calls core DHSs at FDR 0.1, 0.05, 0.01, 0.001 and 0.0001. Requires the `binary` and `neglog10_pvals` layers.
 
Outputs go to `core_sets/`: `core_set_mapping.<column>.tsv` (group value → numeric `path_safe_id`), `<fdr>/<path_safe_id>.fdr<fdr>.*`, and `core_sets.meta.tsv`.
 
### 6. `differential_deseq.nf`: differential accessibility
 
- `-entry diffDeseq --dds <DESeqDataSet.RDS> --total_chunks 500`: runs `bin/r_scripts/differential_deseq.R` in chunks and concatenates the results.
- `-entry diffDeseqSpecies --deseq_meta <tsv with prefix, dataset_anndata>`: runs `differential_deseq.species.R` once per row.
Both use a hard-coded conda env (`/home/sabramov/miniconda3/envs/r-jupyter/`), and `diffDeseq` defaults `--dds` to a project-specific path. Set both before using them.
 
### Experimental: `extract_data_with_offsets.nf`
 
Extracts density and background matrices at jittered DHS positions (`--jitter_matrix`). The `extractDataWithOffsets` sub-workflow is meant to be called from other code. The default `workflow {}` in this file calls it with the wrong number of arguments and will not run as is.
 
### Not maintained
 
`defunc/` (NMF, track averaging, peak filtering, signal propagation) and `bin/defunc/` are kept for reference only.
 
## Notes
 
- The `Altius` profile uses `errorStrategy 'ignore'` for failures other than out-of-memory errors (exit codes 137–143 are retried up to 3 times). A run can therefore finish with some samples missing. Check `.nextflow.log` for ignored tasks, and check that every sample in `samples_file` appears in `samples_order.txt` and the matrices.
- Column order in every matrix follows `samples_order.txt` (row order of `samples_file`). Row order follows the masterlist.
- `build_masterlist.nf` and `normalize_signal.nf` use the `highmem` and `bigmem` labels, which need hundreds of GB to ~2 TB of RAM for large (10k+ sample) indices.
## Authors
 
Jeff Vierstra, Sergey Abramov. Index construction code by Wouter Meuleman and Eric Rynes.
 
