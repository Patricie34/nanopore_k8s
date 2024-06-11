
process SURVIVOR {
	tag "SURVIVOR VCF on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(name), val(sample), path(vcfs)


	output:
	tuple val(name), val(sample), path("${name}.SURVIVOR.BNDs.bcf")

	script:
	//bcftools view -i 'ID~".*BND.*"' $vcf1 -o vcf1.BNDs.vcf
	//bcftools view -i 'ID~".*BND.*"' $vcf2 -o vcf2.BNDs.vcf
	//bcftools view -i 'ID~".*BND.*"' $vcf3 -o vcf3.BNDs.vcf
	//ls *.BNDs.vcf > sample_files
	"""
	echo SURVIVOR $name

	ls *.vcf > sample_files
	cat sample_files
	SURVIVOR merge sample_files 100 1 1 1 0 30 ${name}.SURVIVOR.BNDs.vcf
	bcftools view ${name}.SURVIVOR.BNDs.vcf -O b > ${name}.SURVIVOR.BNDs.bcf

	rm *.vcf
	"""
} 

process SURVIVOR_INFLATE_BNDs {
	tag "SURVIVOR_INFLATE_BNDs on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/survivor", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(name), val(sample), path(vcf), val(varcaller)

	output:
	tuple val(name), val(sample), path("${name}.${varcaller}.inflatedBNDs.vcf")

	script:
	"""
	echo SURVIVOR_INFLATE_BNDs
	bcftools view -i 'ID~".*BND.*"' $vcf > tmp.vcf
	python ${params.ExploadBNDs} -v tmp.vcf -o ${name}.${varcaller}.inflatedBNDs.vcf
	"""
} 

process APPEND_SUPP_FIELD {
	tag "APPEND_SUPP_FIELD on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(name), val(sample), path(inputBcf)

	output:
	tuple val(name), val(sample), path("${name}.SURVIVOR.merged.BNDs.SUPP.vcf")
	tuple val(name), val(sample), path("${name}.temp.vcf")


	script:
	"""
	echo APPEND_SUPP_FIELD $name
	python ${params.SurvivorPutSuppToSamples} -v ${inputBcf} -o ${name}.temp.vcf
	nsamples="\$(bcftools query -l ${inputBcf} | wc -w)"	

	for ((i = 1; i <= \$nsamples; i++)); do
		echo $name >> rename.txt
	done

	#cannot be converted to .bcf because same names cannot be, .vcf is fine. dntknw why
	#bcftools view ${name}.temp.vcf -O b > ${name}.temp.bcf 
	bcftools reheader --samples rename.txt ${name}.temp.vcf -o ${name}.SURVIVOR.merged.BNDs.SUPP.vcf
	"""
} 

process SURVIVOR_INTERSECT_SAMPLES {
	tag "SURVIVOR_INTERSECT_SAMPLES using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/survivor", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "l_mem"

	input:
		path vcfs

	output:
	path "SURVIVOR.samples.intersect.bcf"

	script:
	"""
	echo SURVIVOR_INTERSECT_SAMPLES
	ls *.vcf > sample_files
	SURVIVOR merge sample_files 10 1 1 1 0 30 temp.vcf
	bcftools view temp.vcf -O b > SURVIVOR.samples.intersect.bcf

	rm *.vcf
	"""
} 

process SURVIVOR_FILTER_SINGLETONS {
	tag "SURVIVOR_FILTER_SINGLETONS on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(name), val(sample), path(samplevcf), path(multibcf)

	output:
	tuple val(name), val(sample), path("${name}.SURVIVOR.singletons.vcf")
	tuple val(name), val(sample), path("merged.bcf")


	script:
	"""
	echo SURVIVOR_FILTER_SINGLETONS $name
	bcftools view $samplevcf -O b > sample.bcf
	for bcf_file in *.bcf; do
		echo "Sorting \$bcf_file..."
		bcftools sort -o "\${bcf_file}.sorted.bcf" "\$bcf_file"
	echo "Indexing \$\${bcf_file}.sorted.bcf..."
			bcftools index "\${bcf_file}.sorted.bcf"
	done

	bcftools merge *.sorted.bcf -O b -o merged.bcf
	python ${params.FilterSingletonsSurvivor} -v merged.bcf -n $name -o SURVIVOR.singletons.vcf
	bcftools view --samples ${name} -O v SURVIVOR.singletons.vcf > ${name}.SURVIVOR.singletons.vcf
	"""
} 


process ANNOTATE_SURVIROR {
	tag "ANNOTATE_SURVIROR on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/btk_k8s:25"


	input:
	tuple val(name),val(sample), path(vcf)

	output:
	tuple val(name), val(sample), path("${name}.SURVIVOR.Annotated.vcf")

	script:
	"""
	echo ANNOTATE_SURVIROR $name
	bcftools view -i 'ID~".*BND.*"' ${vcf} > ${name}.BND.vcf
	bcftools sort ${name}.BND.vcf -o ${name}.BND.sorted.vcf

	vep -i ${name}.BND.sorted.vcf --cache --cache_version 95 --dir_cache ${sample.vep} --buffer_size 100 --fasta ${sample.ref} --merged --fork ${task.cpus} --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.SURVIVOR.Annotated.vcf
	"""
} 


process CALC_COVERAGE_SURVIVOR {
	tag "CALC_COVERAGE VCF on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor/debug", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:58"
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(inputvcf), path(bam), path(bai)

	output:
	tuple val(name), path("${name}.CovCytoCNVdb.txt")

	script:
	"""
	echo CALC_COVERAGE_SURVIVOR $name
	samtools view -H $bam | grep @SQ | sed 's/@SQ\tSN:\\|LN://g' | sort -k1,1 -k2,2V > genome.txt 

	bcftools view ${inputvcf} | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") && (\$5!~"GL|MT|KI") && (\$3 ~ "BND") {\$1=\$1;print \$1,\$2-1,\$2,\$3,\$5}' | sort -k1,1 -k2,2V > fromBNDs.bed

	bedtools coverage -sorted -abam fromBNDs.bed -b ${bam} -d -split -g genome.txt > Coverage.txt
	bedtools map -c 4 -a Coverage.txt -b $sample.cytoMap -o concat -g genome.txt > CoverageCytomap.txt
	bedtools map -c 4,4,5 -a CoverageCytomap.txt -b $sample.CNVdb -o collapse,count,collapse -g genome.txt > ${name}.CovCytoCNVdb.txt
	"""
} 

process PARSE_SURVIVOR {
	tag "PARSE_SURVIVOR VCF on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor/", mode:'copy'
	container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:60"
	label "s_cpu"
	label "m_mem"

	input:
	tuple val(name), val(sample), path(inputvcf), path(covfile)

	output:
	tuple val(name), path("${name}.survivor.tsv")

	script:
	"""
	echo PARSE_SURVIVOR $name
	python ${params.ParseSurvivor} -v $inputvcf -c $covfile -o ${name}.survivor.tsv
	"""
}