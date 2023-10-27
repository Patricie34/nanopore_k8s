
process TAG_UNIQUE_VARS{
	tag "Tagging unique vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'

	input:
 tuple val(name), val(sample), path(vcf), path(vcf_filter_with)

	output:
	tuple val(name), val(sample), path("*UniqueTag.vcf"), val(sample)
	
	script:
	"""
	bedtools intersect -a ${vcf} -b ${vcf_filter_with} -v | awk '{print \$3}' > uniqueIDs.txt
 python ${params.TagUniqes} ${vcf} uniqueIDs.txt ${name}
	"""
} 


process ANNOTATE{
	tag "Annotating vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
        container "registry.gitlab.ics.muni.cz:443/450402/btk_k8s:21"
	label "big_mem"

	input:
 tuple val(name),val(sample), path(vcf)

	output:
	tuple val(name), val(sample), path("${name}.UniqueTag.Annotated.vcf")

	script:
	"""
	vep -i ${vcf} --cache --cache_version 90 --dir_cache ${params.GrCh38vep} --fasta ${params.GrCh38ref} --merged --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.UniqueTag.Annotated.vcf
	"""
} 


process EDITVCF {
	tag "Post process VCF on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'

	input:
	tuple val(name), path(inputvcf), path(bam) ,path(bai)

	output:
	tuple val(name), path("Dedup_1000distDistanced_${name}.tsv"), val(sample)
	path('*')
		
	script:
	"""
	#cat ${inputvcf} | grep BND | grep -vE 'GL000|MT|KI27' | awk -v OFS="\\t" '{print \$1,\$2-1,\$2,\$3}' | sort -k1,1 -k2,2V > BNDlocs.bed
	#bedtools coverage causes problems with strange chroms, remove first line
	#bedtools coverage -sorted -abam BNDlocs.bed -b ${bam} -d -split -g genome.txt > Coverage.txt
	#bedtools map -c 4 -a Coverage.txt -b $params.GrCh38cytomap -o concat > CoverageCytomap.txt
	#bedtools map -c 4 -a CoverageCytomap.txt -b $params.GrCh38CNVdb -o collapse,count > CovCytoCNVdb.txt


	cat ${inputvcf} | grep BND | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") && (\$5!~"GL|MT|KI") {\$1=\$1;print \$1,\$2-1,\$2,\$3}' | sort -k1,1 -k2,2V | tail -n +2 > BNDlocs.bed 
	samtools view -H $bam | grep @SQ | sed 's/@SQ\tSN:\\|LN://g' > genome.txt 
	head -n 20 BNDlocs.bed
	head -n 20 genome.txt
	bedtools coverage -sorted -abam BNDlocs.bed -b ${bam} -d -split -g genome.txt \
	| bedtools map -c 4 -a - -b $params.GrCh38cytomap -o concat \
	| bedtools map -c 4,4,5 -a - -b $params.GrCh38CNVdb -o collapse,count,collapse \
	> CovCytoCNVdb.txt

	python ${params.ComputeDistance} ${inputvcf} CovCytoCNVdb.txt Distanced_${name}
	"""
} 