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
	tag "DORADO on $sample.name using $task.cpus CPUs $task.memory and GPU $sample.gpu"
	// publishDir  "${params.outDir}/${sample.name}/nano/mapped", mode:'copy'
 container 'ontresearch/dorado:sha087b7b8d8fc047f531926ba064c2f2503fe9a25a'
	accelerator 1, type: "${sample.gpu}"
	label "xl_cpu"
	label "xxxl_mem"

	input:
	tuple val(name), val(sample)

	output:
	tuple val(sample.name), path("${sample.name}.dorado.bam")
	
	when:
	sample.type == 'raw'

	script:
	"""
	echo DORADO on $sample.name
	# Conditional block to set model based on the value of var
if [[ ${sample.basecall_qual} == "sup" ]]; then
    if [[ ${sample.chemistry} == "v10" ]]; then
        model="dna_r9.4.1_e8_sup@v3.6"
    elif [[ ${sample.chemistry} == "v14" ]]; then
        model="dna_r10.4.1_e8.2_400bps_sup@v4.3.0"
    else
        echo "Unknown value for chemistry"
        exit
    fi
elif [[ ${sample.basecall_qual} == "hac" ]]; then
    if [[ ${sample.chemistry} == "v10" ]]; then
        model="dna_r9.4.1_e8_hac@v3.3"
    elif [[ ${sample.chemistry} == "v14" ]]; then
        model="dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
    else
        echo "Unknown value for chemistry"
        exit
    fi
else
    echo "Unknown value for model"
    exit
fi

 echo "Model: \$model"
	echo downloading model..
	dorado download --model \$model
	
	echo basecalling...
	dorado basecaller \$model ${sample.path} --reference ${sample.ref} --verbose --recursive > ${sample.name}.dorado.bam
	"""
}