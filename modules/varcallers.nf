process SVIM {
	tag "Variant calling using SVIM on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/svim", mode:'copy'
 container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.svim.bcf")
	path "*"
	
	script:
	"""
	echo SVIM on ${name}
	svim alignment ./ ${bam} ${sample.ref} --minimum_depth 2 --read_names --sample ${name}.Sniffles
	#--all_bnds
	mv variants.vcf ${name}.svim.bcf
	"""
} 



process PEPPER_DEEPVARIANT {
	tag "Variant calling using PEPPER_DEEPVARIANT on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/DeepVariant/", mode:'copy'
	accelerator 1, type: 'nvidia.com/gpu'
 container "kishwars/pepper_deepvariant:r0.8-gpu"
	label "s_cpu"
	label "l_mem"

 input:
	tuple val(name), val(sample), path(bam), path(bai)

 output:
	tuple val(name), path ("*")
	tuple val(name), val(sample), path("${name}.variants.vcf")
	
	when:
	name == 'BRNO2641'

	script:
 def model = sample.config == 'dna_r9.4.1_450bps_sup.cfg' ? 'ont_r9_guppy5_sup' : 'ont_r10_q20'

	"""
	echo PEPPER_DEEPVARIANT on ${name}
 run_pepper_margin_deepvariant call_variant -b $bam -f ${sample.ref} -o ./ -p ${name} -s ${name} -t ${task.cpus} --${model} --phased_output
	"""
} 

process DEEPVARIANT {
	tag "Variant calling using DEEPVARIANT on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/DeepVariant/", mode:'copy'
	container "google/deepvariant:1.6.0-gpu"
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), path ("*")
	tuple val(name), val(sample), path("${name}.variants.vcf")

	when:
	sample.config == 'dna_r10.4.1_e8.2_400bps_sup.cfg'	
	
	script:
	"""
	echo DEEPVARIANT on ${name}
	samtools view -b -f 0x900 -q 10 > primary_reads.bam
	/opt/deepvariant/bin/run_deepvariant --model_type=ONT_R104 --ref=${sample.ref} --reads=primary_reads.bam --output_vcf=${name}.DeepVariant --num_shards=${task.cpus}
	"""
} 


process SNIFFLES2 {
	tag "SNIFFLES2 on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/sniffles", mode:'copy'
	container 'quay.io/biocontainers/sniffles:2.2--pyhdfd78af_0'
	label "m_cpu"
	label "m_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${sample.name}.sniffles.vcf")
	
	script:
	"""
	echo sniffles on ${name}
	sniffles --input ${bam} --reference ${sample.ref} --tandem-repeats ${sample.nonMappableRepeats} --minsupport 2 --output-rnames -t ${task.cpus} --vcf ${sample.name}.sniffles.vcf 
	"""
} 