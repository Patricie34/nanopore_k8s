process DELLY_SVs {
	tag "DELLY_SVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Delly/SVs/", mode:'copy'
	container "dellytools/delly:latest"
 label "small_process"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.delly.bcf")
	
	script:
	"""
 delly lr -y ont  -o ${name}.delly.bcf -g ${sample.ref} ${bam}
	"""
} 


process DELLY_merge_SVs {
	tag "DELLY_merge_SVs  using $task.cpus CPUs $task.memory"
	container "dellytools/delly:latest"
 label "small_process"

	input:
 path bcfs_filter

	output:
	path "delly.merged.bcf"
	
	script:
	"""
 delly merge -n 250000000 -o delly.merged.bcf ${bcfs_filter}
	"""
} 

process DELLY_genotype_SVs {
	tag "DELLY_genotype_SVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Delly/SVs/", mode:'copy'
	container "dellytools/delly:latest"
 label "small_process"

	input:
 tuple val(name), val(sample), path(bam), path(bai), path(bcf), path(bcfMerged)

	output:
	tuple val(name), val(sample), path("${name}.geno.bcf")
	
	script:
	""" 
 delly call -g ${sample.ref} -v ${bcfMerged} -o ${name}.geno.bcf -x ${sample.nonMappableRepeats} ${bam}
	"""
} 


process DELLY_merge_genotyped {
	tag "DELLY_merge_genotyped using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/Delly/", mode:'copy'
	container "staphb/bcftools:1.18"
 label "small_process"

	input:
 path genotyped_bcfs

	output:
	path "delly.geno.merged.bcf"
	
	script:
	"""
	# List of files
 file_list=($genotyped_bcfs)

 # Loop over files
 for file in "\${file_list[@]}"; do
     echo "\$file"
					bcftools index \$file
 done
	bcftools merge -m id -O b -o delly.geno.merged.bcf $genotyped_bcfs
	"""
} 

process DELLY_filter_merge_genotyped {
	tag "DELLY_filter_merge_genotyped using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/Delly/", mode:'copy'
	container "quay.io/biocontainers/cyvcf2:0.30.25--py38he403ad2_0"
 label "small_process"

	input:
 path merged_genotyped_bcfs

	output:
	path "filtered.geno.merged.tsv"
	
	script:
	"""
	python $params.FilterSingletons -v $merged_genotyped_bcfs > filtered.geno.merged.tsv
	"""
} 

process DELLY_get_genesBNDs {
	tag "DELLY_get_genesBNDs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Delly/SVs/", mode:'copy'
	container "staphb/bcftools:1.18"
 label "small_process"

	input:
 tuple val(name), val(sample), path(bcf)

	output:
	tuple val(name), val(sample), path("${name}.input.tsv")
	
	script:
	""" 
 bcftools query -f "%CHROM\\t%POS\\t%CHROM\\t%INFO/END\\t%ID\\n" ${bcf} | grep -v "BND" > ${name}.sv.tsv 
	bcftools query -f "%CHROM\\t%POS\\t%INFO/CHR2\\t%INFO/POS2\\t%ID\\n" ${bcf} | grep "BND" >> ${name}.sv.tsv
	cat ${name}.sv.tsv | cut -f 1,2,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Left";}' > ${name}.input.tsv
	cat ${name}.sv.tsv | cut -f 3,4,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Right";}' >> ${name}.input.tsv
	"""
} 

process DELLY_annot_genesBNDs {
	tag "DELLY_annot_genesBNDs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Delly/SVs/", mode:'copy'
	container "trausch/alfred:latest"
 label "small_process"

	input:
 tuple val(name), val(sample), path(input_tsv)

	output:
	tuple val(name), val(sample), path("${name}.sv.gene.tsv")
	
	script:
	""" 
	alfred annotate -d 3000 -g $sample.alfredAnnotateBed -o ${name}.sv.gene.tsv ${input_tsv}
	"""
} 

process DELLY_CNVs {
	tag "DELLY_CNVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Delly/CNVs/", mode:'copy'
	container "dellytools/delly:latest"
 label "small_process"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.cov.gz")
	path "*"
	
	script:
	"""
	# -m parameter should not be hardcoded
	delly cnv -w 25000 -i 25000 -j 25000 -a -g ${sample.ref} -m ${sample.nonMappableRepeatsTobias} -c ${name}.cov.gz -o ${name}.cnv.bcf ${bam}
	"""
}

	process DELLY_CNV_PLOT {
	tag "DELLY_CNVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Delly/CNVs/", mode:'copy'
 label "small_process"

	input:
	tuple val(name), val(sample), path(covfile)

 output:
 path "*"
	
	script:
	"""
	Rscript ${params.CnvPlotTobias} $covfile
	"""
} 