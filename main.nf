
run = "${params.data}".split("/")
println(run)
run = run[run.size()-2]
println(run)


launchDir = "${launchDir}/${run}"

process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${launchDir}/", mode:'copy'

	output:
	path "basecalled/pass/*.fastq"

	script:
	"""
	guppy_basecaller -i ${params.data}/fast5_pass -s ./basecalled --flowcell ${params.flowcell} --kit ${params.kit} -x auto 
	"""
} 

process MAPPING {
	tag "Mapping on $name using $task.cpus CPUs $task.memory"
	publishDir  "${launchDir}/mapped/", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bam.bai")
	tuple val(name), path("${name}.mapped.sam")
	
	script:
	"""
	minimap2 --MD -a /mnt/shared/999993-Bioda/data/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa $reads > ${name}.mapped.sam
	samtools view -bS ${name}.mapped.sam > ${name}.mapped.bam
	samtools sort -o ${name}.sorted.bam ${name}.mapped.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bam.bai	
	"""
} 


process SNIFF {
	tag "Variant calling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${launchDir}/VarCal/", mode:'copy'

	input:
	tuple val(name), path(mapped)
	tuple val(name), path(bai)

	output:
	tuple val(name), path("${name}.vcf")
	
	script:
	"""
	sniffles --input ${name}.sorted.bam --vcf ${run}.vcf
	"""
} 

process QUALIMAP {
        container 'pegi3s/qualimap'
	publishDir  "${launchDir}/", mode:'copy'
	
        input:
        tuple val(name), path (bams)

	output:
	path "qualimap/*"

        script:
        """
        qualimap bamqc -bam ${bams} -nw 5000 -nt 14 --java-mem-size=3G -c -outdir ./qualimap

        """
}

workflow {
println("${params.data}")

rawfastq	= GUPPY_BASECALL()
mapped		= MAPPING(rawfastq.map({ file -> [run, file]}))
mapped[0].view()
		SNIFF(mapped[0],mapped[1])
QUALIMAP(mapped[0])
}

