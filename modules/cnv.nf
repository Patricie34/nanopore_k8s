process CNVkit {
	tag "CNVkit $name using $task.cpus CPUs $task.memory"
	container 'etal/cnvkit:0.9.10'
	publishDir "${params.outDir}/${name}/nano/CNVkit", mode:'copy'
	label "l_mem"
    label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	path "*"

	script:
	"""
	echo CNVKIT on $name
	find /tmp/ -name "*.bam" -exec ls -lh {} \\; >> bamSize.log
	cnvkit.py batch ${bam} -n -m wgs -f ${sample.ref} --annotate ${sample.refFlat} -p $task.cpus
	cat ${name}.remapped.sorted.cns | grep -v GL000 > ${name}.cns
	cat ${name}.remapped.sorted.cnr | grep -v GL000 > ${name}.cnr

	#cnvkit.py segment ${name}.cnr -p $task.cpus -m hmm -o ${name}.cns

	cnvkit.py scatter -s ${name}.cn{s,r} --y-max=4 --y-min=-4 -o ${name}.scatter.pdf
	for Chr in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y'
		do
    	echo "\$Chr"
    	cnvkit.py scatter -s ${name}.cn{s,r} -c \$Chr --y-max=3 --y-min=-3 -o ${name}.chr\${Chr}.pdf
	done
	cnvkit.py diagram -s ${name}.cn{s,r} -o ${name}.diagram.pdf
	"""
}