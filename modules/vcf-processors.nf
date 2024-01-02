
process TAG_UNIQUE_VARS{
	tag "Tagging unique vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	label "big_mem"

	input:
 tuple val(name), val(sample), path(vcf), path(vcf_filter_with)

	output:
	tuple val(name), val(sample), path("*UniqueTag.vcf")
	
	script:
	"""
	bedtools intersect -a ${vcf} -b ${vcf_filter_with} -v | awk '{print \$3}' > uniqueIDs.txt
 python ${params.TagUniqes} ${vcf} uniqueIDs.txt ${name}
	"""
} 


process BCF2TVC {
	tag "Annotating vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/SVs/", mode:'copy'
 label "small_process"

	input:
 tuple val(name),val(sample), path(bcf)

	output:
	path "${name}.Delly.Genotyped.vcf"

	script:
	"""
	bcftools view ${bcf} > ${name}.Delly.Genotyped.vcf
	"""
} 

process ANNOTATE{
	tag "Annotating vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
 //container "registry.gitlab.ics.muni.cz:443/450402/btk_k8s:21"
 container "ensemblorg/ensembl-vep:release_110.1"
 label "big_mem"

	input:
 tuple val(name),val(sample), path(vcf)

	output:
	tuple val(name), val(sample), path("${name}.UniqueTag.Annotated.vcf")

	script:
	"""
	echo annotating $name
	vep -i ${vcf} --cache --cache_version 95 --dir_cache ${params.GrCh38vep} --buffer_size 2000 --fasta ${params.GrCh38ref} --merged --fork ${task.cpus} --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.UniqueTag.Annotated.vcf
	"""
} 

process FILTER_SVIM {
	tag "FILTER_SVIM vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'

	input:
	tuple val(name),val(sample), path(vcf)

	output:
	tuple val(name), val(sample), path("${name}.UniqueTag.Annotated.Sup2.vcf")

	script:
	"""
	echo FILTER_SVIM $name
	cat $vcf | grep -v 'SUPPORT=1' > ${name}.UniqueTag.Annotated.Sup2.vcf
	"""
} 


process PEPPER_PREFILT {
	tag "PEPPER_PREFILT on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/mapped/", mode:'copy'
 container "quay.io/biocontainers/picard:2.27.0--hdfd78af_0"

 input:
	tuple val(name), val(sample), path(bam), path(bai)

 output:
	tuple val(name), val(sample), path("${name}.first.md.bam"), path("${name}.first.md.bam.bai")
	
	when:
	name == 'BRNO2641'

	script:
	"""
	echo PEPPER_PREFILT on ${name}
	sleep infinity
	picard MarkDuplicates I=$bam M=${name}.MD.metrics.txt O=${name}.first.md.bam REMOVE_DUPLICATES=true
	samtools index ${name}.first.md.bam
	"""
} 

process PHASE_WHATSHAP {
	tag "Phasing vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
 container "quay.io/biocontainers/whatshap:2.1--py38h2494328_0"

	input:
 tuple val(name), val(sample), path(vcf), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.phased.vcf")

	script:
	"""
	echo phasing $name
 whatshap phase -o ${name}.phased.vcf --reference=${params.GrCh38ref} ${vcf} ${bam}
	"""
} 

process CALC_COVERAGE {
tag "CALC_COVERAGE VCF on $name using $task.cpus CPUs $task.memory"
publishDir  "${params.outDir}/${name}/nano/VarCal/debug", mode:'copy'

input:
tuple val(name), val(sample), path(inputvcf), path(bam), path(bai)

output:
tuple val(name), path("${name}.CovCytoCNVdb.txt")

script:
"""
echo CALC_COVERAGE $name
samtools view -H $bam | grep @SQ | sed 's/@SQ\tSN:\\|LN://g' > genome.txt 

cat ${inputvcf} | grep svim.BND | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") && (\$5!~"GL|MT|KI") {\$1=\$1;print \$1,\$2-1,\$2,\$3,\$5}' | sort -k1,1 -k2,2V | tail -n +2 > fromBNDs.bed 
cat	fromBNDs.bed | awk '{print \$5}' | grep -oE '[A-Za-z0-9]*:[0-9]*' | sed 's/:/\\t/' | awk -v OFS="\\t" '{print \$1,\$2-1,\$2}' > toBNDs.temp
#add IDs to join in the end
paste toBNDs.temp <(cat fromBNDs.bed | cut -f 4) > toBNDs.bed
bedtools map -c 4 -a <(cat toBNDs.bed | sort -k1,1 -k2,2V ) -b $params.GrCh38cytomap -o concat -g genome.txt > toBNDsCyto.temp

bedtools coverage -sorted -abam fromBNDs.bed -b ${bam} -d -split -g genome.txt > Coverage.txt
bedtools map -c 4 -a Coverage.txt -b $params.GrCh38cytomap -o concat -g genome.txt > CoverageCytomap.txt
bedtools map -c 4,4,5 -a CoverageCytomap.txt -b $params.GrCh38CNVdb -o collapse,count,collapse -g genome.txt > from.CovCytoCNVdb.txt

join -1 4 -2 4 -t \$'\\t' <(sort -k4 from.CovCytoCNVdb.txt) <(sort -k4 toBNDsCyto.temp) | awk -v FS="\\t" -v OFS="\\t" '{print \$2,\$3,\$4,\$1,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15}' > ${name}.CovCytoCNVdb.txt

#bedtools map -c 4 -a <(cat sortedindex.txt| cut -f2- ) -b $params.GrCh38cytomap -o concat -g genome.txt > toBNDsCyto.temp
#sort the toBNDs in order to work with bedtools
#awk '{print NR "\\t" \$0}' toBNDs.bed | sort -k2,2 -k3,3V > sortedindex.txt 
#bedtools coverage -sorted -abam <(cat sortedindex.txt| cut -f2- ) -b ${bam} -d -split -g genome.txt NOT YET
#revert sort
#awk 'NR==FNR{a[\$0];next} (FNR in a)' <(cut -f1 sortedindex.txt) toBNDsCyto.temp > toBNDsCyto.bed
#join -1 4 -2 4 -t \$'\\t' from.CovCytoCNVdb.txt toBNDsCyto.temp > joined.txt
#get rid of column 5, no meaning, coming from coverage
#paste from.CovCytoCNVdb.txt toBNDsCyto.bed | cut -f 1-4,6- > ${name}.CovCytoCNVdb.txt

 """
} 

process EDITVCF {
tag "EDITVCF VCF on $name using $task.cpus CPUs $task.memory"
publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'

input:
tuple val(name), val(sample), path(inputvcf), path(bam), path(bai), path(covFile)

output:
tuple val(name), val(sample), path("Dedup.1000dfilt.${name}.tsv")
path('*.tsv')

script:
"""
echo python parsing $name
python ${params.ComputeDistance} ${inputvcf} ${covFile} ${name}
"""
} 