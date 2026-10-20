#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { add_matrices_to_anndata; convert_to_numpy; extract_meta_from_anndata } from "./converters"

params.conda = "$moduleDir/environment.yml"


process convert_regions_to_summits {
    conda params.conda

    input:
        path masterlist

    output:
        tuple path(summit_masterlist), path(orig_idx)

    script:
    summit_masterlist = "masterlist.summits.bed"
    orig_idx = "masterlist.summits.orig_idx.txt"
    """
    awk -F'\t' -v OFS='\t' \
        '{ print \$1, \$5, \$5 + 1, \$4, NR }' ${masterlist} \
        | sort-bed - > ${summit_masterlist}

    awk '{print \$NF - 1}' ${summit_masterlist} > ${orig_idx}
    """
}

// Scatters a np_array (values for a subset or all of the masterlist, in orig_idx's order)
// back into a full-length n_dhs array, filling positions absent from orig_idx with 0.
process restore_masterlist_order {

    conda params.conda
    tag "${prefix}"

    input:
        tuple val(prefix), path(np_array, name: "input/*"), path(inverse_argsort)

    output:
        tuple val(prefix), path(name)

    script:
    name = np_array.name
    """
    python3 $moduleDir/bin/helpers/restore_masterlist_order.py \
        ${np_array} \
        ${orig_idx} \
        ${n_dhs} \
        ${name}
    """
}


process project_peak_calls {
    conda params.conda
    tag "${sample_id}"
    scratch true

    input:
		tuple path(masterlist), val(sample_id), path(peaks_file)
    
    output:
        tuple val(sample_id), val(suffix), path(name)
    
    script:
    suffix = 'binary'
    name = "${sample_id}.${suffix}.npy"
    """
    zcat ${peaks_file} \
        | grep -v "^#" \
        | bedmap --fraction-either 0.5 \
        --indicator \
        --sweep-all \
        ${masterlist} \
        - > tmp.txt

    python3 $moduleDir/bin/helpers/txt_to_np.py tmp.txt ${name} --dtype bool
    """
}


process extract_max_density {
    conda params.conda
    tag "${sample_id}"
    scratch true

    input:
        tuple path(masterlist), val(sample_id), path(density_bw)
    
    output:
        tuple val(sample_id), val(suffix), path(name)
    
    script:
    suffix = "density"
    name = "${sample_id}.${suffix}.npy"
    """
    bigWigToBedGraph ${density_bw} tmp.bg

    cat tmp.bg \
        | awk -v OFS='\t' '{print \$1,\$2,\$3,"${sample_id}",\$4}' \
        | bedmap --sweep-all \
            --delim "\t" \
            --max ${masterlist} - \
        | sed 's/\\<NAN\\>/0/g' \
        > tmp.txt

    python3 $moduleDir/bin/helpers/txt_to_np.py tmp.txt ${name} --dtype float
    """
}

process count_tags {
	tag "${sample_id}"
	conda "${params.conda}"
	scratch true

	input:
		tuple path(saf), val(sample_id), path(bam_file), path(bam_file_index)

	output:
		tuple val(sample_id), val(suffix), path(name)

	script:
    suffix = "counts"
    name = "${sample_id}.${suffix}.npy"
    ext = bam_file.extension
	"""
    if [ ${ext} != 'bam' ]; then 
        samtools view -bh \
            --reference ${params.genome_fasta} \
            ${bam_file} > align.bam
        samtools index align.bam
    else
        ln -s ${bam_file} align.bam
        ln -s ${bam_file_index} align.bam.bai
    fi

    count=`samtools stats align.bam \
        | grep "^SN" \
        | grep "mapped and paired" \
        | cut -f3` && [[ \$count -gt 0 ]] && tag="-p" || tag=""

	featureCounts -a ${saf} -O -o counts.txt -F SAF \$tag align.bam
	cat counts.txt \
        | awk 'NR > 2 {print \$(NF)}' > tmp.txt
        
    python3 $moduleDir/bin/helpers/txt_to_np.py tmp.txt ${name} --dtype int
	"""
}

process extract_max_pval {
    tag "${sample_id}"
    conda params.conda

    input:
        tuple path(masterlist), val(sample_id), path(pvals_parquet)

    output:
        tuple val(sample_id), val(suffix), path(name)

    script:
    suffix = 'neglog10_pvals'
    name = "${sample_id}.${suffix}.npy"
    """
    hotspot3-pvals \
        ${pvals_parquet} \
        ${masterlist} \
        ${name} \
        --chrom_sizes ${params.nuclear_chrom_sizes} \
        --format npy
    """
}

// returns a np array of shape n_dhs x 2
process extract_bg_mean {
    conda params.conda
    tag "${sample_id}"
    label "bg_params"

    input:
        tuple path(masterlist), val(sample_id), path(bg_params_tabix)
    
    output:
        tuple val(sample_id), val(suffix), path(name)
    
    script:
    suffix = "mean_bg_agg_cutcounts"
    name = "${sample_id}.${suffix}.npy"
    """
    echo -e "chrom_masterlist\tstart_masterlist\tend_masterlist\tdhs_id\tsummit\t\$(head -1 <(zcat ${bg_params_tabix}))" > tmp.bed

    zcat ${bg_params_tabix} \
        | grep -v "^#" \
        | grep segment \
        | bedtools intersect \
            -loj \
            -a ${masterlist} \
            -b stdin \
            -sorted \
            >> tmp.bed

    python3 $moduleDir/bin/helpers/extract_bg_params.py tmp.bed ${name}
    """
}

// merge a np array of shape n_dhs x 2
process generate_matrix {
	publishDir "${params.outdir}/raw_matrices", pattern: "${name}"
	
	label "highmem"
    conda "${params.conda}"
	scratch true
    tag "${prefix}"

	input:
		tuple val(prefix), path(files), path(samples_order)

	output:
		tuple val(prefix), path(name)

	script:
    name = "matrix.${prefix}.npy"
	"""
    python3 $moduleDir/bin/helpers/matrix_from_vectors.py \
        ${prefix} \
        ${samples_order} \
        ${name} \
        --input_ext npy
	"""
}

 workflow getSummitBasedColumns {
    take:
        data // sample_id, bg_params_tabix, density_bw, masterlist, orig_idx
        n_dhs
    main:
        orig_idx = data
            | map(it -> tuple(it[0], it[4])) // sample_id, orig_idx

        density_cols = data
            | map(it -> tuple(it[3], it[0], it[2])) // masterlist, sample_id, density_bw
            | extract_max_density

        summit_based_cols = data
            | map(it -> tuple(it[3], it[0], it[1])) // masterlist, sample_id, bg_params_tabix
            | extract_bg_mean // sample_id, suffix, np_array
            | mix(density_cols)
            | combine(orig_idx, by: 0) // sample_id, suffix, np_array, orig_idx
            | map(it -> tuple(it[1], it[2], it[3])) // suffix, np_array, orig_idx
            | combine(n_dhs) // suffix, np_array, orig_idx, n_dhs
            | restore_masterlist_order

    emit:
        summit_based_cols
 }


workflow generateMatrices {
    take:
        data // masterlist, samples_order, saf_masterlist, id, bam, bam_index, peaks, density, fit_stats, hotspot3_pvals_parquet
    main:

        binary_cols = data
            | map(it -> tuple(it[0], it[3], it[6]))
            | project_peak_calls
        
        count_cols = data
            | map(it -> tuple(it[2], it[3], it[4], it[5]))
			| count_tags
        
        max_pvals = data
            | map(it -> tuple(it[0], it[3], it[9]))
            | extract_max_pval

        // Summit based extraction
        masterlist = data
            | map(it -> it[0])
            | first()

        n_dhs = masterlist.countLines().toInteger()

        summits_masterlist = masterlist
            | convert_regions_to_summits // summits_masterlist, orig_idx

        samples_order = data
            | map(it -> it[1])
            | first()

        samples_masterlists = data
            | map(it -> it[3])
            | combine(summits_masterlist)

        summit_based_cols = getSummitBasedColumns(
            data
                | map(it -> tuple(it[3], it[8], it[7])) // sample_id, bg_params_tabix, density_bw
                | join(samples_masterlists), // sample_id, bg_params_tabix, density_bw, summits_masterlist, orig_idx
            n_dhs
        )

        out = binary_cols
            | mix(count_cols)
            | mix(max_pvals) // sample_id, suffix, np_array
            | map(it -> tuple(it[1], it[2])) // suffix, np_array
            | mix(summit_based_cols)
            | combine(samples_order.countLines().toInteger())
            | map(it -> tuple(groupKey(it[0], it[2]), it[1]))
            | groupTuple(by: 0) // suffix, columns
            | combine(samples_order) // suffix, columns, samples_order
            | generate_matrix // matrix_name, matrix
    emit:
        out // matrix_name, matrix
}


workflow generateMatricesFromAnndata {
    take:
        bams_hotspots
    main:
        print "Looking for anndata input (generated by build_masterlist.nf) in - params.index_anndata = ${params.index_anndata}"

        index_anndata = Channel.of(params.index_anndata)

        matrices = index_anndata
            | extract_meta_from_anndata // masterlist , saf_masterlist, samples_order, samples_meta
            | map(it -> tuple(it[0], it[2], it[1])) // masterlist, samples_order, saf_masterlist
            | combine(bams_hotspots)
            | generateMatrices
            | map(it -> it[1])
            | collect(sort: true, flat: true)

        add_matrices_to_anndata(
            index_anndata,
            matrices
        )
    emit:
        add_matrices_to_anndata.out
}

// main workflow 
workflow {
    bams_hotspots = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.cram_file),
            file(row?.cram_index ?: "${row.cram_file}.crai"),
            file(row[params.matrix_peaks_column]),
            file(row.normalized_density_bw),
            file(row.hotspot3_fit_stats_file),
            file(row.hotspot3_pvals_parquet)
        ))
        | generateMatricesFromAnndata
}


workflow extractDensityAtSummit {
    println "Extracting normalized density at summits from ${params.index_anndata}"
    data = Channel.of(params.index_anndata)
        | extract_meta_from_anndata // masterlist , saf_masterlist, samples_order, samples_meta
    
    masterlist = data
        | map(it -> it[0])

    n_dhs = masterlist.countLines().toInteger()

    summits_masterlist = masterlist
        | convert_regions_to_summits // summits_masterlist, orig_idx

    samples_meta = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.normalized_density_bw)
        ))
        | combine(summits_masterlist.map(it -> it[0])) // summits_masterlist
        | map(it -> tuple(it[2], it[0], it[1]))
        | extract_max_density // sample_id, suffix, np_array
        | map(it -> tuple(it[1], it[2])) // suffix, np_array
        | combine(summits_masterlist.map(it -> it[1])) // suffix, np_array, orig_idx
        | combine(n_dhs) // suffix, np_array, orig_idx, n_dhs
        | restore_masterlist_order
        | groupTuple()
        | combine(data.map(it -> it[2])) // suffix, columns, samples_order
        | generate_matrix
}


////////////////////// DEFUNC  //////////////////////
workflow extractBGParams {
    println "Extracting bg params at ${params.reference_bed}"
    reference_bed = Channel.fromPath(params.reference_bed)
    samples_meta = Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(
            row.sample_id,
            file(row.peak_stats),
            file("${row.peak_stats}.tbi")
        ))
        | combine(reference_bed)
        | map(it -> tuple(it[3], it[0], it[1], it[2]))
        | extract_bg_params
}
