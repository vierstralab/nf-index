import argparse
import pandas as pd
from genome_tools.data.anndata import (
    read_zarr_backed,
    add_matrices_to_anndata,
    get_matrices_mapping_by_types,
)
import anndata as ad
from scipy import sparse as sp


def main(adata: ad.AnnData, matrices):
    matrices_types = [
        'counts', 'binary', 'density',
        'mean_bg_agg_cutcounts', 'neglog10_pvals'
    ] + list(adata.layers.keys())
    matrices_mapping = get_matrices_mapping_by_types(matrices, matrices_types)

    add_matrices_to_anndata(adata, matrices_mapping)

    return adata


def merge_samples_meta(adata_obj, samples_meta_file):
    samples_meta = pd.read_table(samples_meta_file, low_memory=False).set_index('sample_id')
    samples_meta[
        samples_meta.select_dtypes(include=['object', 'O', 'category']).columns
    ] = samples_meta.select_dtypes(include=['object', 'O', 'category']).fillna('None')

    missing_indices = adata_obj.obs.index.difference(samples_meta.index)
    if len(missing_indices) > 0:
        print(f"Warning: Missing indices: {missing_indices}")
        print("Creating a new anndata object (dropping X and obs)")
        return ad.AnnData(X=None, obs=samples_meta, var=adata_obj.var)

    samples_meta = samples_meta.loc[adata_obj.obs.index, :]
    overlap = samples_meta.columns.intersection(adata_obj.obs.columns)

    for col in overlap:
        if not samples_meta[col].equals(adata_obj.obs[col]):
            print(f"Column '{col}' differs between samples_meta and adata.obs. Using samples meta as refeence")

    adata_obj.obs = samples_meta.drop(columns=overlap).join(adata_obj.obs)
    return adata_obj


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add matrices to an existing AnnData object as layers.")
    parser.add_argument("input", help="Input Zarr file")
    parser.add_argument("output", help="Output Zarr file")
    parser.add_argument("matrices", nargs="+", help="Matrices to add to the AnnData layers")
    parser.add_argument(
        "--samples_meta",
        default=None,
        help="Optional samples metadata table (tsv with a sample_id column) to merge into adata.obs"
    )
    args = parser.parse_args()

    adata_obj = read_zarr_backed(args.input)

    if args.samples_meta is not None:
        adata_obj = merge_samples_meta(adata_obj, args.samples_meta)

    main(adata_obj, args.matrices).write_zarr(args.output)
