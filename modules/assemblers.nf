
process FLYE {
	container 'staphb/flye'
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'

	input:
	tuple val(name), path (fastq)

	output:
	path "Flye_output/*"

	script:
	"""
	flye --nano-hq ${fastq} --threads ${task.cpus} -m 1 -o ./Flye_output
	"""
}


process ASSEMBLY_PREFILTER {
	tag "Prefiltering BNDs on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/SHASTA/", mode:'copy'

	input:
	tuple val(name), path(bam), path(bai), path(vars)

	output:
	tuple val(name), path("${name}.BND.fasta")
	path 'readnames.txt'
	path 'filtered.bam'


	script:
	"""
 cat ${vars} | grep -oP '(?<=READS\\=)([\\w\\-]*,)*[\\w\\-]*' | awk  -v FS="," -v OFS="\n" '{\$1=\$1;print}' > readnames.txt
	samtools view ${bam} -b -N readnames.txt > filtered.bam
	samtools fasta filtered.bam > ${name}.BND.fasta
	"""
}


process SHASTA {
	publishDir  "${params.outDir}/${name}/nano/SHASTA/", mode:'copy'
	tag "Shasta on $name using $task.cpus CPUs $task.memory"

	input:
	tuple val(name), path(fasta)
	//tuple val(name), path(bam)


	output:
	path "ShastaRun/*"
	tuple val(name), path("ShastaRun/Assembly.fasta")
 
	script:
	"""
	#shasta --help
	#gunzip -d --force *.fastq.gz
	#cat *.fastq* | sed -n '1~4s/^@/>/p;2~4p' > all.fasta
	#cat *.fastq* | samtools fasta - > all.fasta
	#cat all.fasta | wc
#--memoryMode filesystem --memoryBacking disk
	echo Running_Shasta
	cat $fasta | wc -l
	shasta --input $fasta --config Nanopore-May2022 --threads ${task.cpus} --Reads.minReadLength 5000
	"""
}


process MUMMER {
 container 'staphb/mummer'
	//container 'quay.io/biocontainers/mummer'
	publishDir  "${params.outDir}/${name}/nano/Dnadiff", mode:'copy'
	tag "Mummer on $name using $task.cpus CPUs $task.memory"

	input:
	tuple val(name), path(assembly)

	output:
	path "*"

	script:
	"""
	#dnadiff -p dnadiff -t $task.cpus ${params.GrCh38ref} ${assembly}
	#	sleep infinity

	nucmer -t $task.cpus ${params.GrCh38ref} ${assembly}
#run-mummer3 -t $task.cpus ${params.GrCh38ref} ${assembly}
#mummerplot out.delta
	"""
}