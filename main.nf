nextflow.enable.dsl = 2
launchDir = "${launchDir}/${params.datain.replaceAll(".*/","")}"


process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/basecalled/", mode:'copy'
	
	input:
	tuple val(name), path(f5)
	
	output:
	path '*'

	script:
	"""
	guppy_basecaller -i 
	"""
}


 
workflow {
 rawfast5 = channel.fromPath("${params.data}/*.fast5", checkIfExists: true)
	
rawfast5.view()   
}
