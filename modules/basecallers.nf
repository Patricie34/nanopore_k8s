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
	publishDir  "${params.outDir}/${sample.name}/nano/", mode:'copy'
 container 'ontresearch/dorado:sha087b7b8d8fc047f531926ba064c2f2503fe9a25a'
	accelerator 1, type: 'nvidia.com/gpu'
	label "l_cpu"
	label "l_mem"

	input:
	val sample

	output:
	tuple val(sample.name), path("basecalled/pass/*.gz")
	
	script:
	"""
	sleep infinity
	dorado download --model dna_r10.4.1_e8.2_400bps_hac@v4.3.0
	dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 /storage/share/110000-MED/110323-imgg/111323-plevova/01.NanoBreak/data/raw_data/BRNO1837/ --emit-fastq > test.fastq
	"""
} 