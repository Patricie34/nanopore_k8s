process MINIMAP2 {
	tag "Mapping on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${sample.name}/nano/mapped/", mode:'copy'
	label "xl_cpu"
	label "xxl_mem"

	input:
	tuple val(name), path(reads), val(sample)

	output:
	tuple val(sample.name), val(sample), path("${name}*.bam")

	script:
	"""
	echo MINIMAP2 on $sample.name
	minimap2 -ax map-ont -t ${task.cpus} ${sample.ref} $reads > ${name}.sam
	samtools view -bS ${name}.sam > ${name}.remapped.bam
	"""
}

process BAM_INDEX_SORT {
	tag "Mapping on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${sample.name}/nano/mapped/", mode:'copy'
	label "l_mem"

	input:
	tuple val(name), val(sample), path(bam)

	output:
	tuple val(sample.name), val(sample), path("${name}.remapped.sorted.bam"), path("${name}.remapped.sorted.bam.bai")

	script:
	"""
	echo BAM_INDEX_SORT on $name
	samtools sort -o ${name}.remapped.sorted.bam ${bam}
	samtools index ${name}.remapped.sorted.bam ${name}.remapped.sorted.bam.bai	
	"""
} 

process BAM2FASTQ {
	tag "BAM2FASTQ on $sample.name using $task.cpus CPUs $task.memory"
	//publishDir  "${params.outDir}/${sample.name}/nano/mapped/", mode:'copy'
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(sample.name), val(sample), path("${name}.fastq")

	script:
	"""
	echo BAM2FASTQ on $sample.name
	samtools fastq $bam > ${name}.fastq
	"""
}