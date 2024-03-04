process GUPPY {
	tag "Basecalling on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${sample.name}/nano/", mode:'copy'
	container "cerit.io/docker/genomicpariscentre/guppy-gpu"
	accelerator 1, type: 'nvidia.com/gpu'
	label "l_cpu"
	label "l_mem"

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


process DORADO {
	tag "DORADO on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${sample.name}/nano/mapped", mode:'copy'
 container 'ontresearch/dorado:sha087b7b8d8fc047f531926ba064c2f2503fe9a25a'
	accelerator 1, type: 'nvidia.com/mig-1g.10gb'
	//'nvidia.com/gpu'
	label "l_cpu"
	label "l_mem"
	debug true

	input:
	val sample

	output:
	tuple val(sample.name), path("${sample.name}.dorado.bam")
	
	when:
	sample.type == 'raw'

	script:
	"""
	echo DORADO on $sample.name
	# Conditional block to set model based on the value of var
	if [[ ${sample.config} == "sup" ]]; then
					model="dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
	elif [[ ${sample.config} == "hac" ]]; then
					model="dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
	else
					echo "Unknown value for model"
					exit
	fi
 echo "Model: \$model"

	echo downloading model..
	dorado download --model \$model
	
	echo basecalling...
	dorado basecaller \$model ${sample.path} --reference ${sample.ref} --recursive > ${sample.name}.dorado.bam
		"""
}