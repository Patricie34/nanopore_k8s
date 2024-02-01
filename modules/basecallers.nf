process GUPPY {
	tag "Basecalling on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${sample.name}/nano/", mode:'copy'
	container "cerit.io/docker/genomicpariscentre/guppy-gpu"
	accelerator 1, type: 'nvidia.com/gpu'
	label "medium_mem"

	input:
	val sample

	output:
	tuple val(sample.name), path("basecalled/pass/*.gz")
	
	when:
	sample.type == 'fast5'

	script:
	"""
	/opt/ont/guppy/bin/guppy_basecaller -i ${sample.path} -s ./basecalled -c ${sample.config} --compress_fastq --recursive --num_callers ${task.cpus} -x auto --progress_stats_frequency 999 --verbose_logs --trace_category_logs --fast5_out
	"""
} 