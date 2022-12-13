run = "${params.data}".split("/")
run = run[run.size()-1]
resultsDir = "${params.outDir}/${run}/nano"



process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${resultsDir}/", mode:'copy'
	accelerator 1, type: 'nvidia.com/mig-1g.10gb'

	output:
	path "basecalled/pass/*.fastq"

	script:
	"""
	guppy_basecaller -i ${params.data}/ -s ./basecalled --flowcell ${params.flowcell} --kit ${params.kit} --recursive --num_callers $task.cpus -x auto
	"""
} 

process MAPPING {
    label "big_cpu"


	tag "Mapping on $name using $task.cpus CPUs $task.memory"
	publishDir  "${resultsDir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bam.bai")
	
	script:
	"""
	minimap2 --MD -a -t 4*${task.cpus} ${params.ref} $reads > ${name}.mapped.sam
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
	echo $PATH
	svim alignment ./ ${mapped} ${params.ref} --minimum_depth 1 --read_names --all_bnds
echo Computing distance...
	python ${params.ComputeDistance} variants.vcf Distanced_${name}
	"""
} 



process QUALIMAP {
	container 'pegi3s/qualimap'
	publishDir  "${resultsDir}/", mode:'copy'
		    label "big_mem"

	
	input:
	 tuple val(name), path (bams)

	output:
	path "qualimap/*"

	 script:
	  """
	 qualimap bamqc -bam ${bams} -nw 5000 -nt 14 --java-mem-size=60G -c -outdir ./qualimap
	 """
}

process FLYE {
	container 'staphb/flye'
	publishDir  "${resultsDir}/", mode:'copy'
		    label "big_mem"


	input:
	tuple val(name), path (fastq)

	output:
	path "Flye_output/*"

	script:
	"""
	flye --nano-hq ${fastq} --threads 16 -m 1 -o ./Flye_output
	"""
}

process SHASTA {
 container 'cimendes/shasta:0.6.0'
	publishDir  "${resultsDir}/SHASTA", mode:'copy'

	input:
	path(fastq)
	tuple val(name), path(reads)


	output:
	path "*"

	script:
	"""
	shasta --help
	#cat *.fastq | sed -n '1~4s/^@/>/p;2~4p' > all.fasta
	#shasta --input all.fasta --MarkerGraph.minCoverage 1 --config Nanopore-May2022
	"""
}

process CNVKIT {
 container 'etal/cnvkit'
	publishDir  "${resultsDir}/CNVkit", mode:'copy'

	input:
	tuple val(name), path(mapped)
	tuple val(name), path(bai)

	output:
	path "*"

	script:
	"""
	cnvkit.py batch ${mapped} -n -m wgs -f ${params.ref} --annotate ${params.refFlat} -p $task.cpus --target-avg-size 30000
	cat ${name}.sorted.cns | grep -v GL000 > ${name}.cns
	cat ${name}.sorted.cnr | grep -v GL000 > ${name}.cnr
	cnvkit.py scatter -s ${name}.cn{s,r} -o ${name}.scatter.svg
	cnvkit.py diagram -s ${name}.cn{s,r}
	"""
}


workflow {
println("${params.data}")
println("Run name: ${run}")


rawfastq	= GUPPY_BASECALL()
//FLYE(rawfastq.map(rawfastq)
mapped		= MAPPING(rawfastq.map({ file -> [run, file]}))
mapped[0].view()

SVIM(mapped[0],mapped[1])
CNVKIT(mapped[0],mapped[1])
QUALIMAP(mapped[0])
rawfastq.collect().view()
SHASTA(rawfastq.collect(),mapped[0])
}

