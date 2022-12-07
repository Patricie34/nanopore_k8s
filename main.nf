run = "${params.data}".split("/")
run = run[run.size()-1]
resultsDir = "${params.outDir}/${run}"



process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/", mode:'copy'

	output:
	path "basecalled/pass/*.fastq"

	script:
	"""
	guppy_basecaller -i ${params.data}/ -s ./basecalled --flowcell ${params.flowcell} --kit ${params.kit} --recursive --num_callers 4 -x auto
	"""
} 

process MAPPING {
	tag "Mapping on $name using $task.cpus CPUs $task.memory"
	publishDir  "${resultsDir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bam.bai")
	
	script:
	"""
	minimap2 --MD -a ${params.ref} $reads > ${name}.mapped.sam
	samtools view -bS ${name}.mapped.sam > ${name}.mapped.bam
	samtools sort -o ${name}.sorted.bam ${name}.mapped.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bam.bai	
	"""
} 

process SVIM{
	tag "Variant calling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${resultsDir}/VarCal/", mode:'copy'

	input:
	tuple val(name), path(mapped)
	tuple val(name), path(bai)

	output:
	path '*'
	
	script:
	"""
	svim alignment ./ ${mapped} ${params.ref} --minimum_depth 1 --read_names --all_bnds
	python ${params.calcdist} variants.vcf > Distanced_${name}.tsv
	"""
	

} 



process QUALIMAP {
 container 'pegi3s/qualimap'
	publishDir  "${resultsDir}/", mode:'copy'
	
 input:
  tuple val(name), path (bams)

	output:
	path "qualimap/*"

 script:
   """
   qualimap bamqc -bam ${bams} -nw 5000 -nt 14 --java-mem-size=30G -c -outdir ./qualimap
   """
}

workflow {
println("${params.data}")
println("Run name: ${run}")


rawfastq	= GUPPY_BASECALL()
mapped		= MAPPING(rawfastq.map({ file -> [run, file]}))
mapped[0].view()

SVIM(mapped[0],mapped[1])
QUALIMAP(mapped[0])
}

