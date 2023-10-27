
process SVIM{
	tag "Variant calling using SVIM on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	label "medium_mem"
	label "medium_cpus"

	input:
	tuple val(name), val(sample), path(bam), path(bai)


	output:
	tuple val(name), path ("*")
	tuple val(name), val(sample), path("${name}.variants.vcf")
	
	script:
	"""
 find /tmp/ -name "*.bam" -exec ls -lh {} \\; >> bamSize.log
	svim alignment ./ ${bam} ${params.GrCh38ref} --minimum_depth 1 --read_names --all_bnds
	mv variants.vcf ${name}.variants.vcf
	"""
} 

process SNIFFLES{
	tag "Variant calling using SNIFFLES on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	label "medium_mem"
	label "medium_cpus"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), path ("*")
	tuple val(name), val(sample), path("${name}.variants.vcf")
	
	script:
	"""
 find /tmp/ -name "*.bam" -exec ls -lh {} \\; >> bamSize.log
	sniffles --input ${bam} --vcf ${sample.name}.vcf --reference ${params.GrCh38ref} --mosaic
	"""
} 