
process CIRCOS{
	tag "Creating circos on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Circos/", mode:'copy'

	input:
	tuple val(name), path(vcfs), path(cnv_sorted), path(vars_edited)

	output:
	path("*.html")
	
	script:
	"""

	cat $vars_edited | awk -v FS="\\t" '(\$3!~"GL|MT|KI") && (\$9!~"GL|MT|KI") {print \$2,\$3,\$4,\$9,\$10,\$15}' > vars_slimmed.tsv
	head -n 20 vars_slimmed.tsv
	Rscript --vanilla ${params.Circos} $cnv_sorted vars_slimmed.tsv $name ${params.GrCh38lens}
	"""
} 

process HEATMAP{
	tag "Creating heatmap on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Heatmap/", mode:'copy'
	container "rnakato/juicer:1.6.2"
	label "small_process"

	input:
	tuple val(name), path(deduplicatedTsv)

	output:
	tuple val(name), path("*.hic")
	
	script:
	"""
	cat $deduplicatedTsv| awk '{printf "%s chr%s %s %s %s chr%s %s %s\\n", 0,\$3,\$4,0,12,\$9,\$10,1;}' FS='\\t' > ${name}_preHiC.txt
	cat ${name}_preHiC.txt | awk '{if (\$2 <= \$6) print \$0; else print \$5,\$6,\$7,\$8,\$1,\$2,\$3,\$4}' OFS=" " > ${name}_preHic_orderedChroms.txt
	cat ${name}_preHic_orderedChroms.txt| sort -k2,2d -k6,6d -k3,3n -k7,7n > ${name}_preHiC_sorted.txt
	#java -Xms512m -Xmx32384m -jar
	java -Xmx2g -jar /opt/juicer/scripts/common/juicer_tools.jar pre ${name}_preHiC_sorted.txt ${name}.hic hg38 -n
	"""
} 