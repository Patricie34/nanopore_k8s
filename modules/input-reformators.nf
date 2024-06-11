process REFORMAT_SAMPLE {
	tag "Reformating $sample.name using $task.cpus CPUs $task.memory"
	label "s_cpu"
	label "s_mem"

	input:
	val sample

	output:
	tuple val(sample.name),val(sample.reference), val(sample)

 """ """ //this is not an error
}


process REFORMAT_PARAMS {
	tag "Reformating $references.refname using $task.cpus CPUs $task.memory"
	label "s_cpu"
	label "s_mem"
 
	input:
	val references

	output:
	tuple val(references), val(references.refname)

	""" """ //this is not an error
}
