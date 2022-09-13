#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process count_tags {
	tag "${indiv_id}"
	scratch true
	conda params.conda

	input:
		tuple val(indiv_id), val(bam_file), val(peaks_file)

	output:
		tuple val(indiv_id), file("${prefix}.counts.txt"), file("${prefix}.bin.txt")

	script:
	prefix = "${indiv_id}"
	"""
	awk '{print \$1}' ${params.chrom_sizes} > chroms.txt
	while read chrom; do
		echo \$chrom
		cat ${params.index_file} | grep "\$chrom\t" > chrom_index.txt
		if [ -s chrom_index.txt ]; then
			bedtools intersect -sorted -g ${params.chrom_sizes} -c -a chrom_index.txt -b ${bam_file} | awk '{print \$(NF)}' >> ${prefix}.counts.txt
		fi
	done < chroms.txt
	bedmap --indicator ${params.index_file} ${peaks_file} > ${prefix}.bin.txt
	"""
}

process generate_count_matrix {

	publishDir params.outdir + '/index'
	conda params.conda

	input:
		tuple val(indiv_ids), file(count_files), file(bin_files)

	output:
		tuple file("matrix.all.signal.txt.gz"), file("matrix.all.peaks.txt.gz"), file("indivs_order.txt")

	script:
	indiv_ids_join = indiv_ids.join("\t")
	"""
	echo "${indiv_ids_join}" > indivs_order.txt
	paste - ${count_files} | cut -c2- | gzip -c > matrix.all.signal.txt.gz
	paste - ${bin_files} | cut -c2- | gzip -c > matrix.all.peaks.txt.gz
	"""
}


workflow generateMatrix {
	take:
		BAMS_HOTSPOTS
	main:
		COUNT_FILES = count_tags(BAMS_HOTSPOTS)
		generate_count_matrix(COUNT_FILES.groupTuple(by: 0))
	emit:
		generate_count_matrix.out
}
workflow {
	BAMS_HOTSPOTS = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map{ row -> tuple(row.ag_id, row.bam_file, row.hotspots_file) }
	generateMatrix(BAMS_HOTSPOTS)
}