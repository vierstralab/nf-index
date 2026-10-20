#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { extract_meta_from_anndata } from "./converters"
include { extract_max_density; extract_bg_mean; restore_masterlist_order; generate_matrix } from "./generate_matrices"

params.conda = "$moduleDir/environment.yml"


process apply_jitter {
    conda params.conda
    tag "${sample_id}"

    input:
        tuple path(masterlist), val(sample_id), val(jitter_id), path(jitter_matrix)

    output:
        tuple val(sample_id), val(prefix), path(name), path(orig_idx)

    script:
    prefix = "${sample_id}-${jitter_id}"
    name = "${prefix}.jittered_masterlist.bed"
    orig_idx = "${prefix}.jittered_masterlist.orig_idx.txt"
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
        data // channel: jitter_id, jitter_matrix, anndata
    main:
        meta = data
            | map(it -> it[2])
            | unique() // anndata should be the same for all entries
            | extract_meta_from_anndata // masterlist, saf_masterlist, samples_order, samples_meta_file

        masterlist = meta
            | map(it -> it[0])
            | first()
        
        n_dhs = masterlist.countLines().toInteger()
        
        samples_order = meta
            | map(it -> it[2])
            | first()
        
        samples_meta = meta 
            | map(it -> it[3])
            | splitCsv(header: true, sep: '\t')
            | map(row -> tuple(
                row.sample_id,
                file(row.hotspot3_fit_stats_file),
                file(row.normalized_density_bw),
            ))

        // sample_id, jittered_masterlist (sampled DHSs only), orig_idx
        jittered_masterlists = masterlist
            | combine(samples_meta.map(it -> it[0])) // masterlist, sample_id
            | combine(data.map(it -> tuple(it[0], it[1]))) // masterlist, sample_id, jitter_id, jitter_matrix
            | apply_jitter // sample_id, joint_id, jittered_masterlist, orig_idx


        
        merged_data = samples_meta // sample_id, bg_params_tabix, density_bw
            | combine(jittered_masterlists, by: 0) // sample_id, bg_params_tabix, density_bw, joint_id, jittered_masterlist, orig_idx

        orig_idx = merged_data
            | map(it -> tuple(it[3], it[5])) // joint_id, orig_idx

        density_cols = merged_data
            | map(it -> tuple(it[4], it[3], it[2])) // masterlist, sample_id, density_bw
            | extract_max_density

        bg_cols = merged_data
            | map(it -> tuple(it[4], it[3], it[1])) // masterlist, sample_id, bg_params_tabix
            | extract_bg_mean // sample_id, suffix, np_array

        out = bg_cols
            | mix(density_cols)
            | combine(orig_idx, by: 0) // sample_id, suffix, np_array, orig_idx
            | combine(n_dhs) // sample_id, suffix, np_array, orig_idx, n_dhs
            | restore_masterlist_order // sample_id, suffix, column
            | combine(samples_order.countLines().toInteger()) 
            | map(
                it -> tuple(
                    groupKey(tuple(it[1], it[0].split('-')[1]),
                    it[3]),
                    it[2]
                )
            )
            | combine(samples_order) // suffix, columns, samples_order
            | generate_matrix

    emit:
        out
}


// main workflow
workflow {
    println "Extracting data with offsets from ${params.index_anndata}"

    index_anndata = Channel.of(params.index_anndata)

    samples_meta = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.hotspot3_fit_stats_file),
            file(row.normalized_density_bw),
        ))

    jitter_matrix = Channel.of(
        tuple(
            "custom_jitter",
            file(params.jitter_matrix)
        )
    )

    out = extractDataWithOffsets(index_anndata, samples_meta, jitter_matrix)

}
