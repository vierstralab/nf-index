#!/usr/bin/env nextflow

params.samples_file = '/net/seq/data/projects/regulotyping-h.CD3+/metadata.txt'
params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'
params.genome_fasta_file = "/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa"

//Build cell-type specific indicies as well as a combined index
params.build_ct_index = 1 

params.outdir='output'

nuclear_chroms="$params.genome" + ".nuclear.txt"
chrom_sizes="$params.genome"  + ".chrom_sizes"
chrom_sizes_bed="$params.genome"  + ".chrom_sizes.bed"
mappable="$params.genome" + ".K76.mappable_only.bed"
centers="$params.genome" + ".K76.center_sites.n100.nuclear.starch"

// Read samples file
Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.indiv_id, row.cell_type, row.bam_file ) }
	.groupTuple(by: [0, 1])
	.map{ it -> tuple(it[0], it[1], it[2].join(" ")) }
	.set{ SAMPLES_AGGREGATIONS_MERGE }

process merge_bamfiles {
	tag "${indiv_id}:${cell_type}"

	publishDir params.outdir + '/merged', mode: 'symlink' 

        module "samtools/1.3"

	cpus 2

	input:
	set val(indiv_id), val(cell_type), val(bam_files) from SAMPLES_AGGREGATIONS_MERGE

	output:
	set val(indiv_id), val(cell_type), file('*.bam'), file('*.bam.bai') into BAMS_HOTSPOTS

	script:
	"""
	samtools merge -f -@${task.cpus} --reference ${params.genome_fasta_file} ${indiv_id}_${cell_type}.bam ${bam_files}
	samtools index ${indiv_id}_${cell_type}.bam
	"""
}	

process call_hotspots {
	tag "${indiv_id}:${cell_type}"

	// only publish varw_peaks and hotspots
	publishDir params.outdir + '/hotspots', mode: 'copy'

	module "modwt/1.0:kentutil/302:bedops/2.4.35-typical:bedtools/2.25.0:hotspot2/2.1.2:samtools/1.3"

	input:
	file 'nuclear_chroms.txt' from file("${nuclear_chroms}")
	file 'mappable.bed' from file("${mappable}")
	file 'chrom_sizes.txt' from file("${chrom_sizes}")
	file 'chrom_sizes.bed' from file("${chrom_sizes_bed}")
	file 'centers.starch' from file("${centers}")

	set val(indiv_id), val(cell_type), file(bam_file), file(bam_index_file) from BAMS_HOTSPOTS

	output:
	set val(indiv_id), val(cell_type), file(bam_file), file(bam_index_file), file("${indiv_id}_${cell_type}.varw_peaks.fdr0.001.starch") into PEAKS
	file("${indiv_id}_${cell_type}.hotspots.fdr*.starch")
	file("${indiv_id}_${cell_type}.SPOT.fdr0.05.txt")
	file("${indiv_id}_${cell_type}.normalized.density.starch")
	file("${indiv_id}_${cell_type}.normalized.density.bw")

	script:
	"""
	TMPDIR=\$(mktemp -d)

	echo "Temporary directory =  \${TMPDIR}"
	
	samtools view -H ${bam_file} > \${TMPDIR}/header.txt

	cat nuclear_chroms.txt \
	| xargs samtools view -b ${bam_file} \
	| samtools reheader \${TMPDIR}/header.txt - \
	> \${TMPDIR}/nuclear.bam

	hotspot2.sh -F 0.05 -f 0.05 -p varWidth_20_${indiv_id}_${cell_type} \
		-M mappable.bed \
		-c chrom_sizes.bed \
		-C centers.starch \
		\${TMPDIR}/nuclear.bam \
		peaks

	cd peaks

	hsmerge.sh -f 0.001 nuclear.allcalls.starch nuclear.hotspots.fdr0.001.starch

	rm -f nuclear.varw_peaks.*

	density-peaks.bash \
		\${TMPDIR} \
		varWidth_20_${indiv_id}_${cell_type} \
		nuclear.cutcounts.starch \
		nuclear.hotspots.fdr0.001.starch \
		../chrom_sizes.bed \
		nuclear.varw_density.fdr0.001.starch \
		nuclear.varw_peaks.fdr0.001.starch \
		\$(cat nuclear.cleavage.total)

	cp nuclear.varw_peaks.fdr0.001.starch ../${indiv_id}_${cell_type}.varw_peaks.fdr0.001.starch
	cp nuclear.hotspots.fdr0.05.starch ../${indiv_id}_${cell_type}.hotspots.fdr0.05.starch
	cp nuclear.hotspots.fdr0.001.starch ../${indiv_id}_${cell_type}.hotspots.fdr0.001.starch
	cp nuclear.SPOT.fdr0.05.txt ../${indiv_id}_${cell_type}.SPOT.fdr0.05.txt

	tagcounts=\$(samtools view -c \${TMPDIR}/nuclear.bam)
	echo "tagcounts = \${tagcounts}"

	unstarch nuclear.density.starch \
		| awk -v tagcount=\${tagcounts} \
		      -v scale=1000000 \
			  -v OFS="\\t" \
			  '{ z=\$5; n=(z/tagcount)*scale; print \$1, \$2, \$3, \$4, n }' \
		| starch - \
	> ../${indiv_id}_${cell_type}.normalized.density.starch

	unstarch ../${indiv_id}_${cell_type}.normalized.density.starch \
		| awk -v OFS="\\t" \
			  -v bin=20 \
			  'BEGIN { \
				  chr=""; \
				} \
				{ \
					if( \$5 != 0 && \$2 != 0 ){ \
						if( chr=="" ) { \
							chr=\$1; \
							print "variableStep chrom=" chr " span=" bin; \
						} \
						if( \$1==chr ){ \
							print \$2, \$5; \
						} else { \
							chr=\$1; \
							print "variableStep chrom=" chr " span=" bin; \
							print \$2, \$5; \
						} \
					} \
				}' \
	> \${TMPDIR}/normalized.density.wig

	wigToBigWig -clip \${TMPDIR}/normalized.density.wig ../chrom_sizes.txt ../${indiv_id}_${cell_type}.normalized.density.bw

	rm -r \${TMPDIR}
	"""
}

PEAKS.into{PEAK_LIST;PEAK_FILES}

// PEAK_FILES
// 	.map{ it -> tuple(it[1], it[3]) }
// 	.groupTuple(by: 0)
// 	.tap{PEAK_FILES_BY_CELLTYPE}
// 	.map{ it -> it[1].flatten() }
// 	.set{PEAK_FILES_ALL}

PEAK_FILES
	.map{ it -> tuple(it[1], it[4]) }
	.groupTuple(by: 0)
	.tap{PEAK_FILES_BY_CELLTYPE}
	.map{ it -> tuple("all", it[1]) }
	.groupTuple(by: 0)
	.map{ it -> tuple(it[0], it[1].flatten())}
	.set{PEAK_FILES_ALL}

// Include cell-type specific indices or just one large index covering all samples
PEAK_INDEX_FILES = params.build_ct_index ? PEAK_FILES_ALL.concat(PEAK_FILES_BY_CELLTYPE) : PEAK_FILES_ALL

process build_index {
	tag "${cell_type}"
	
	module "R/4.0.5:bedops/2.4.35-typical:kentutil/388"
	// R module should have caTools package installed

	publishDir params.outdir + '/index', mode: 'copy' 

	input:
	set val(cell_type), file('*') from PEAK_INDEX_FILES
	file chrom_sizes from file("${chrom_sizes_bed}")

	output:
	file "masterlist*"
	set val(cell_type), file("masterlist_DHSs_*_nonovl_any_chunkIDs.bed") into INDEX_FILES, INDEX_FILES_FOR_ANNOTATION

	script:
	"""
	ls *.varw_peaks.fdr0.001.starch > filelist.txt

	/home/jvierstra/.local/src/Index/run_sequential.sh \
		\$(pwd) \
		${chrom_sizes} \
		filelist.txt \
		${cell_type}

	"""
}

PEAK_LIST
	.map{ it -> tuple(it[1], it[1], it[0], it[2], it[3], it[4])}
	.tap{ PEAK_LIST_BY_CELLTYPE }
	.map{ it -> tuple("all", it[1], it[2], it[3], it[4], it[5])}
	.set{ PEAK_LIST_ALL }

PEAK_LIST_COMBINED = params.build_ct_index ? PEAK_LIST_ALL.concat(PEAK_LIST_BY_CELLTYPE) : PEAK_LIST_ALL

process count_tags {
	tag "${index_id}:${indiv_id}:${cell_type}"
	
	conda '/home/jvierstra/.local/miniconda3/envs/py3.9_default'
	module "bedops/2.4.35-typical"

	publishDir params.outdir + '/counts', mode: 'symlink'

	input:
	set val(index_id), val(cell_type), val(indiv_id), file(bam_file), file(bam_index_file), file(peaks_file), file(index_file) from PEAK_LIST_COMBINED.combine(INDEX_FILES, by: 0)

	output:
	set val(index_id), val(cell_type), val(indiv_id), file("${index_id}_${indiv_id}_${cell_type}.counts.txt"), file("${index_id}_${indiv_id}_${cell_type}.bin.txt") into COUNTS_FILES

	script:
	"""
	count_tags.py ${bam_file} < ${index_file} > ${index_id}_${indiv_id}_${cell_type}.counts.txt
	
	bedmap --indicator ${index_file} ${peaks_file} > ${index_id}_${indiv_id}_${cell_type}.bin.txt
	"""
}

process generate_count_matrix {
	tag "${cell_type}"
	publishDir params.outdir + '/index', mode: 'copy' 

	input:
	set val(index_id), val(cell_types), val(indiv_ids), file(count_files), file(bin_files), file(index_file) from COUNTS_FILES.groupTuple(by: 0).combine(INDEX_FILES_FOR_ANNOTATION, by: 0)

	output:
	file "matrix_*.txt.gz"

	script:
	col_names = [indiv_ids, cell_types].transpose().collect{e -> e.join(">") }
	"""
	echo -n "region_id" > header.txt
	echo -e "\\t${col_names.join("\t")}" >> header.txt

	cat header.txt <(cut -f4 ${index_file} | paste - ${count_files}) | gzip -c >  matrix_counts.${index_id}.txt.gz
	cat header.txt <(cut -f4 ${index_file} | paste - ${bin_files}) | gzip -c >  matrix_bin.${index_id}.txt.gz
	"""
}

//COUNTS_FILES.combine(INDEX_FILES_FOR_ANNOTATION, by: 0).groupTuple(by: 0).println()
