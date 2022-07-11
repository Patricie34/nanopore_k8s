nextflow.enable.dsl = 2

process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"

	script:
	"""
	guppy_basecaller -i $params.data -s ${launchDir}/basecalled/ --flowcell ${params.flowcell} --kit ${params.kit} -x auto 
	cat ${launchDir}/basecalled/*.fastq >  ${launchDir}/basecalled/all_guppy.fastq
	"""
} 

process MAPPING {
	tag "Mapping on $name using $task.cpus CPUs $task.memory"
	publishDir  "${launchDir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
 	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bai")
	script:
	"""
	minimap2 --MD -a /mnt/shared/999993-Bioda/data/ssd_3/references/homsap/GRCh37-p13/seq/ $reads > ${name}.mapped.sam	\
	samtools view –bS ${name}.mapped.sam > ${name}.mapped.bam
	samtools sort -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai	
	"""
} 



workflow {
println("${params.data}")

GUPPY_BASECALL()
rawfastq = channel
			.fromPath("${launchDir}/basecalled/*.fastq", checkIfExists: true)
			.map({ file -> [file.getSimpleName(), file]})
			.view()

MAPPING(rawfastq)
}
