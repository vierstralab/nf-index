# AGENTS.md — nf-index
 
Guidance for AI coding agents (Claude Code and others) working in this repository. `README.md` is the human documentation; this file is the operational summary for running, debugging and modifying the pipelines.
 
## What this repo is
 
The second stage of the paper's processing. It takes per-sample [hotspot3](https://github.com/vierstralab/hotspot3) outputs and filtered CRAMs from [`nf-aligning`](https://github.com/vierstralab/nf-aligning), builds a consensus index of regulatory elements (Meuleman et al. 2020 method), projects all samples onto it, and normalizes the signal. It has been used for both **DHS** indices (from hotspot3 peaks) and **consensus footprint** indices (from per-sample footprint calls). The code and outputs always say "DHS" (`dhs_id`, `autosomal_dhs`, `masterlist_DHSs_*`), whatever the element type. Input BEDs must have `chrom, start, end, sample_id, score, density_summit, wavelet_summit` as their first 7 columns. Everything is passed between stages as **zarr-backed AnnData**:
 
| Order | Script | Reads | Writes (under `outdir`, default `$launchDir/output`) |
|---|---|---|---|
| 1 | `build_masterlist.nf` | `samples_file` (peak column) | `zarr/index.anndata.zarr` |
| 2 | `generate_matrices.nf` | `index_anndata` + `samples_file` | `zarr/index+matrices.anndata.zarr` |
| 3 | `normalize_signal.nf` | `matrices_anndata` | `index+matrices+normalized.anndata.zarr` (note: in `outdir`, not `outdir/zarr`) |
| 4a | `core_sets.nf` | `matrices_anndata`, `--grouping_column` | `core_sets/` |
| 4b | `differential_deseq.nf` | a DESeq2 `.RDS` | DE tables |
 
The default anndata paths are `${outdir}/zarr/...`, so **launch stages 1–3 from the same directory**. Otherwise pass `--index_anndata` / `--matrices_anndata` explicitly.
 
## How to run
 
```bash
nextflow run build_masterlist.nf  -profile Altius,new_cluster -resume --samples_file S.tsv
nextflow run generate_matrices.nf -profile Altius,new_cluster -resume --samples_file S.tsv
nextflow run normalize_signal.nf  -profile Altius,new_cluster -resume
nextflow run core_sets.nf         -profile Altius,new_cluster -resume --grouping_column <col>
```
 
Other entries (`-entry`):
 
- `build_masterlist.nf`: `annotateMasterlist` (re-annotate an existing build; uses `--index_dir`).
- `generate_matrices.nf`: `extractDensityAtSummit`.
- `normalize_signal.nf`: `existingModel` and `existingModelToScaleFactors`. Both need `--template_anndata` pointing at a previously normalized anndata.
- `differential_deseq.nf`: `diffDeseq` and `diffDeseqSpecies`.
`variance_partition.nf` has a standalone `workflow {}` but is normally called by `normalize_signal.nf`.
 
Do not use:
 
- `extract_data_with_offsets.nf` main workflow: it calls `extractDataWithOffsets` with 3 arguments while that workflow `take:`s 1, so it errors.
- `generate_matrices.nf` `extractBGParams`: references the undefined process `extract_bg_params`.
- Anything in `defunc/` or `bin/defunc/`.
## The samples file (TSV, header, one row per sample)
 
- `sample_id`: required everywhere. It **must equal column 4 of that sample's peak file**; the index scripts use column 4 as `sampleID`, and [hotspot3](https://github.com/vierstralab/hotspot3) writes the sample ID there. If a peak file's column 4 isn't in `samples_order.txt`, `write_rows` fails with `Sample name "…" was not found in the list of sample IDs`. Because that failure is ignored under the `Altius` profile, the binary matrix can end up incomplete.
- Build: the column named by `index_peaks_column` (default `peaks_file_0.001fdr`).
- Generate: `matrix_peaks_column` (default `peaks_file_0.01fdr`), `cram_file`, optional `cram_index`, `normalized_density_bw`, `hotspot3_fit_stats_file`, `hotspot3_pvals_parquet`.
- Normalize: any covariates named in `vst_design_formula` (default `core_annotation`, `extended_annotation`) and `variance_partition_formula` (default `nuclear_dup_rate`, `SPOT3_score`, `nuclear_reads`, `alignment_quality`). These reach the anndata `obs` only through the `samples_file` given to `generate_matrices.nf`, so include them there.
Matrix column order = row order of `samples_file` at build time (`samples_order.txt`). Matrix row order = masterlist order.
 
## Environment
 
- Every process uses `conda params.conda`. `params.config` points to an Altius-specific env path. Elsewhere, create the env from `environment.yml` and set `--conda`.
- Note that `generate_matrices.nf`, `normalize_signal.nf` and `extract_data_with_offsets.nf` declare `params.conda = "$moduleDir/environment.yml"` in the script. Values from `params.config` and the command line take precedence over that.
- Not in `environment.yml`, but required on `PATH`:
  - `hotspot3-pvals` (pip `hotspot3`), for `extract_max_pval`.
  - `normalize-matrix` and `qc-normalization`, for `normalize_matrix`.
  - `featureCounts` (subread, which is in the env).
  - The prebuilt binary `bin/index_scripts/writeRowsPerChunkForMatrices`; rebuild it from the `.cpp` if it fails with exec or GLIBC errors.
- `genome_tools` (pip, from GitHub) supplies `read_zarr_backed`, `add_matrices_to_anndata` and related helpers used by `bin/convert_to_anndata/*.py`.
- `params.config` == `params.human.config` (GRCh38). Use `params.mouse.config` for mm10. All paths are absolute paths on the Altius `/net/seq` filesystem.
## Resource and failure behaviour
 
- Labels: `medmem` (50 GB × attempt²), `highmem` (450 GB × attempt), `bigmem` (50 CPUs, 1.9 TB, `bigmem` queue), `bg_params`. `normalize_matrix`, `prepare_deseq_dataset` and `deseq2_vst` are `bigmem`.
- In the `Altius` profile, exit codes 137–143 are retried up to 3 times. **Any other failure is `ignore`d**, and `annotate_masterlist` is always `ignore`. A run can therefore "complete" with missing samples or matrices. After each stage check `.nextflow.log` for `IGNORED`/`FAILED`, and check that the output anndata has the expected `n_obs` and layers:
```python
from genome_tools.data.anndata import read_zarr_backed
a = read_zarr_backed("output/zarr/index+matrices.anndata.zarr")
print(a.shape, list(a.layers))
```
 
## Code map
 
- `*.nf` at the top level: workflows. `converters.nf` holds the shared AnnData read/write processes (`convert_index_to_anndata`, `add_matrices_to_anndata`, `add_normalized_matrices_to_anndata`, `extract_meta_from_anndata`, `convert_to_numpy`).
- `bin/index_scripts/`: Meuleman/Rynes DHS index code (R + awk + C++). Treat it as vendored and avoid editing it.
- `bin/helpers/`: per-sample vector → `.npy` and matrix assembly helpers.
- `bin/convert_to_anndata/`: AnnData I/O.
- `bin/r_scripts/`: DESeq2 VST, variancePartition, differential tests.
- `bin/core_sets/`: core-set statistics.
- `bin/annotations/`: GC, GENCODE and repeat annotation shell scripts.
## Working rules for agents
 
- Don't start full runs unless asked. They submit large SLURM jobs (`queueSize = 800`, TB-scale memory). Test with a few samples, a scratch launch directory and `-resume`.
- Never delete `work/` or existing `output/zarr/*` unless told to. Later stages and `-resume` depend on them.
- Before changing a process, grep for its uses. `extract_max_density`, `extract_bg_mean`, `restore_masterlist_order` and `generate_matrix` from `generate_matrices.nf` are reused by `extract_data_with_offsets.nf`, and the `converters.nf` processes are used everywhere.
- Process inputs and outputs are positional tuples, and channel shapes are documented in inline comments (`// sample_id, suffix, np_array`). Keep those comments accurate when you change a tuple.
- Update `README.md` when you change entry points, samples-file columns, params or output names.
- Check syntax and params without running: `nextflow config` (resolved params) and `nextflow run <file>.nf -preview`.
 
