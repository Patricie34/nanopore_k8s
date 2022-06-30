nextflow.enable.dsl = 2

process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"

	script:
	"""
	guppy_basecaller -i $params.data -s ${launchDir}/basecalled/ --flowcell FLO-FLG001 --kit SQK-LSK110 -x auto 
	"""
}


 
workflow {
println("${params.data}")
rawfast5 = channel
		.fromPath("${params.data}/*.fast5", checkIfExists: true)
		.map({ file -> [file, file.getSimpleName() ]})
		.view()

GUPPY_BASECALL()
}
