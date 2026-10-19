#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { extract_meta_from_anndata } from "./converters"
include { extract_max_density; extract_bg_mean; restore_masterlist_order; generate_matrix } from "./generate_matrices"

params.conda = "$moduleDir/environment.yml"


process apply_jitter {
    conda params.conda
    tag "${sample_id}"

    input:
        tuple path(masterlist), val(sample_id), path(jitter_matrix)

    output:
        tuple val(sample_id), path(name), path(orig_idx)

    script:
    name = "${sample_id}.jittered_masterlist.bed"
    orig_idx = "${sample_id}.jittered_masterlist.orig_idx.txt"
    """
    python3 $moduleDir/bin/helpers/apply_jitter.py \
        ${masterlist} \
        ${sample_id} \
        ${jitter_matrix} \
        ${name} \
        ${orig_idx}
    """
}


// Extracts per-sample summit-based matrices, restricted to (and jittered around)
// only the DHSs that were sampled for that sample/epoch — see vinson's sample_examples.py.
workflow extractDataWithOffsets {
    take:
        index_anndata   // channel: path to the reference index anndata
        samples_meta    // channel: tuple(sample_id, hotspot3_fit_stats_file, normalized_density_bw)
        jitter_matrix   // channel: path to per-sample jitter offsets .npz, keyed by sample_id
    main:
        meta = index_anndata
            | extract_meta_from_anndata // masterlist, saf_masterlist, samples_order, samples_meta_file

        masterlist = meta
            | map(it -> it[0])
            | first()

        n_dhs = masterlist.countLines().toInteger()

        // sample_id, jittered_masterlist (sampled DHSs only), orig_idx
        jittered_masterlists = masterlist
            | combine(samples_meta.map(it -> it[0])) // masterlist, sample_id
            | combine(jitter_matrix) // masterlist, sample_id, jitter_matrix
            | apply_jitter

        samples_order = meta
            | map(it -> it[2])
            | first()

        // sample_id, bg_params_tabix, density_bw, jittered_masterlist, orig_idx
        data = samples_meta
            | join(jittered_masterlists)

        orig_idx = data
            | map(it -> tuple(it[0], it[4])) // sample_id, orig_idx

        density_cols = data
            | map(it -> tuple(it[3], it[0], it[2])) // masterlist, sample_id, density_bw
            | extract_max_density

        out = data
            | map(it -> tuple(it[3], it[0], it[1])) // masterlist, sample_id, bg_params_tabix
            | extract_bg_mean // sample_id, suffix, np_array
            | mix(density_cols)
            | combine(orig_idx, by: 0) // sample_id, suffix, np_array, orig_idx
            | map(it -> tuple(it[1], it[2], it[3])) // suffix, np_array, orig_idx
            | combine(n_dhs) // suffix, np_array, orig_idx, n_dhs
            | restore_masterlist_order
            | combine(samples_order.countLines().toInteger())
            | map(it -> tuple(groupKey(it[0], it[2]), it[1]))
            | groupTuple(by: 0) // suffix, columns
            | combine(samples_order) // suffix, columns, samples_order
            | generate_matrix
    emit:
        out
}


workflow extractDataWithOffsetsFromParams {
    println "Extracting data with offsets from ${params.index_anndata}"

    index_anndata = Channel.of(params.index_anndata)

    samples_meta = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.hotspot3_fit_stats_file),
            file(row.normalized_density_bw),
        ))

    jitter_matrix = Channel.fromPath(params.jitter_matrix)

    out = extractDataWithOffsets(index_anndata, samples_meta, jitter_matrix)

}

// main workflow
workflow {
    extractDataWithOffsetsFromParams()
}
