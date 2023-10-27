process MINIMAP2 {
	tag "Mapping on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${sample.name}/nano/mapped/", mode:'copy'
	label "medium_cpus"

	input:
	tuple val(name), path(reads), val(sample)

	output:
	tuple val(sample.name), val(sample), path("${name}.sorted.bam"), path("${name}.sorted.bam.bai")

	script:
	"""
	minimap2 --MD -a -t ${task.cpus} ${sample.ref} $reads > ${name}.sam
	samtools view -bS ${name}.sam > ${name}.bam
	samtools sort -o ${name}.sorted.bam ${name}.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bam.bai	
	"""
} 