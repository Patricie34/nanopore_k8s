process QUALIMAP {
	tag "Qualimap $name using $task.cpus CPUs $task.memory"
	container 'pegi3s/qualimap'
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'
	label "s_cpu"
	label "xxl_mem"
	
	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	path "qualimap/*"

	script:
	"""
	echo QUALIMAP on $name
	qualimap bamqc -bam ${bam} -nw 5000 -nt 14 --java-mem-size=60G -c -outdir ./qualimap
	"""
}