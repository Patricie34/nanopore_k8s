nextflow.enable.dsl = 2
k8s.storageClaimName = 'pvc-janek-storage-elixir1-cerit-sc-cz'

process GUPPY_BASECALL {
	containerOptions '--gpus all'
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"

	script:
	"""
	guppy_basecaller -i $params.data -s ${launchDir}/basecalled/ --flowcell ${params.flowcell} --kit ${params.kit} -x auto 
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
