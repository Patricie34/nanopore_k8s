
process COLLECT_FQs {
	tag "Collecting Fastq on $sample.name using $task.cpus CPUs $task.memory"
	label "small_process"

	input:
	val sample

	output:
		tuple val(sample.name), path("basecalled/pass/*.gz")

	when:
		sample.type == 'fastq'

	script:
	"""
	echo "Collecting fastqs for $sample.name"
	mkdir -p basecalled/pass/
	find ${sample.path} -type f -name '*.fastq.gz' | xargs -I % cp % ./basecalled/pass/
	"""
} 

process COLLECT_BAMs {
	tag "Collecting bam on $sample.name using $task.cpus CPUs $task.memory"
	label "small_process"

	input:
	tuple val(name), val(sample)

	output:
	tuple val(sample.name), val(sample), path("${sample.name}.sorted.bam"), path("${sample.name}.sorted.bam.bai")
	
	when:
		sample.type == 'bam'

	script:
	"""
	echo "Collecting bams for $sample.name"
	find ${sample.path} -type f -name '*.ba*' | xargs -I % cp % .
	"""
} 