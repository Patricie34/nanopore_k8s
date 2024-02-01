
process TAG_UNIQUE_VARS {
	tag "Tagging unique vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	label "big_mem"

	input:
 tuple val(name), val(sample), path(vcf), path(vcf_filter_with)

	output:
	tuple val(name), val(sample), path("*UniqueTag.vcf")
	
	script:
	"""
	echo TAG_UNIQUE_VARS on $name
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

process ANNOTATE {
	tag "Annotating vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
 container "ensemblorg/ensembl-vep:release_110.1"
 label "big_mem"

	input:
 tuple val(name),val(sample), path(vcf)

	output:
	tuple val(name), val(sample), path("${name}.UniqueTag.Annotated.vcf")

	script:
	"""
	echo ANNOTATE $name
	vep -i ${vcf} --cache --cache_version 95 --dir_cache ${params.GrCh38vep} --buffer_size 2000 --fasta ${params.GrCh38ref} --merged --fork ${task.cpus} --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.UniqueTag.Annotated.vcf
	"""
} 


process ANNOTATE_SURVIROR {
	tag "ANNOTATE_SURVIROR on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/survivor", mode:'copy'
 container "registry.gitlab.ics.muni.cz:443/450402/btk_k8s:25"
 // label "big_mem"
 memory '128 GB'

	input:
 tuple val(name),val(sample), path(vcf)

	output:
	tuple val(name), val(sample), path("${name}.SURVIVOR.Annotated.vcf")

	script:
	"""
	echo ANNOTATE_SURVIROR $name
	bcftools view -i 'ID~".*BND.*"' ${vcf} > ${name}.BND.vcf
	bcftools sort ${name}.BND.vcf -o ${name}.BND.sorted.vcf

	vep -i ${name}.BND.sorted.vcf --cache --cache_version 95 --dir_cache ${params.GrCh38vep} --buffer_size 100 --fasta ${params.GrCh38ref} --merged --fork ${task.cpus} --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.SURVIVOR.Annotated.vcf
	"""
} 

process FILTER_SVIM {
	tag "FILTER_SVIM vcf vars on $name using $task.cpus CPUs $task.memory"
	//publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
 label "small_process"

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
	publishDir  "${params.outDir}/${name}/nano/mapped/pepper_prefilt", mode:'copy'
	container "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
 label "small_process"

 input:
	tuple val(name), val(sample), path(bam), path(bai)

 output:
	tuple val(name), val(sample), path("${name}.sambamba.sorted.dedup.bam"), path("${name}.sambamba.sorted.dedup.bam.bai")
	
	when:
	name == 'BRNO2641'

	script:
	"""
	echo PEPPER_PREFILT on ${name}
	sambamba markdup -r -t ${task.cpus} $bam ${name}.first.md.bam
	sambamba sort -t ${task.cpus} -o ${name}.sambamba.sorted.dedup.bam ${name}.first.md.bam
	sambamba index -t ${task.cpus} ${name}.sambamba.sorted.dedup.bam
	"""
} 


process PEPPER_UMIREADS {
	tag "PEPPER_UMIREADS on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/mapped/pepper_prefilt", mode:'copy'
 container "registry.gitlab.ics.muni.cz:443/450402/qc_cmbg:26"
 label "small_process"

 input:
	tuple val(name), val(sample), path(bam), path(bai)

 output:
	tuple val(name), val(sample), path("${name}.uniqueNames.bam"), path("${name}.uniqueNames.bam.bai")
	
	when:
	name == 'BRNO2641'

	script:
	"""
	echo PEPPER_UMIREADS on $name
 samtools view -h ${bam} | \
 awk -v output_bam=${name}.uniqueNames.bam \
     'BEGIN {OFS="\\t"} 
      { 
         if (\$1 ~ /^@/) {
             print \$0;
         } else {
             cmd = "openssl rand -hex 4";
             cmd | getline random_string;
             close(cmd);
             \$1 = \$1 "_" random_string;
             print \$0;
         }
      }' | \
 samtools view -b -o ${name}.uniqueNames.bam -
 samtools index ${name}.uniqueNames.bam ${name}.uniqueNames.bam.bai
	"""
} 

process PHASE_WHATSHAP {
	tag "Phasing vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
 container "registry.gitlab.ics.muni.cz:443/450402/qc_cmbg:26"

	input:
 tuple val(name), val(sample), path(vcf), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.phased.vcf")

	script:
	"""
	echo PHASE_WHATSHAP $name
 whatshap phase -o ${name}.phased.vcf --reference=${params.GrCh38ref} ${vcf} ${bam}
	"""
} 

process CALC_COVERAGE {
 tag "CALC_COVERAGE VCF on $name using $task.cpus CPUs $task.memory"
 //publishDir "${params.outDir}/${name}/nano/VarCal/debug", mode:'copy'

 input:
 tuple val(name), val(sample), path(inputvcf), path(bam), path(bai)

 output:
 tuple val(name), path("${name}.CovCytoCNVdb.txt")

 script:
 """
 echo CALC_COVERAGE $name
 samtools view -H $bam | grep @SQ | sed 's/@SQ\tSN:\\|LN://g' > genome.txt 

 cat ${inputvcf} | grep svim.BND | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") && (\$5!~"GL|MT|KI") {\$1=\$1;print \$1,\$2-1,\$2,\$3,\$5}' | sort -k1,1 -k2,2V | tail -n  +2 > fromBNDs.bed 
 cat	fromBNDs.bed | awk '{print \$5}' | grep -oE '[A-Za-z0-9]*:[0-9]*' | sed 's/:/\\t/' | awk -v OFS="\\t" '{print \$1,\$2-1,\$2}' > toBNDs.temp
 #add IDs to join in the end
 paste toBNDs.temp <(cat fromBNDs.bed | cut -f 4) > toBNDs.bed
 bedtools map -c 4 -a <(cat toBNDs.bed | sort -k1,1 -k2,2V ) -b $params.GrCh38cytomap -o concat -g genome.txt > toBNDsCyto.temp

 bedtools coverage -sorted -abam fromBNDs.bed -b ${bam} -d -split -g genome.txt > Coverage.txt
 bedtools map -c 4 -a Coverage.txt -b $params.GrCh38cytomap -o concat -g genome.txt > CoverageCytomap.txt
 bedtools map -c 4,4,5 -a CoverageCytomap.txt -b $params.GrCh38CNVdb -o collapse,count,collapse -g genome.txt > from.CovCytoCNVdb.txt

 join -1 4 -2 4 -t \$'\\t' <(sort -k4 from.CovCytoCNVdb.txt) <(sort -k4 toBNDsCyto.temp) | awk -v FS="\\t" -v OFS="\\t" '{print \$2,\$3,\$4,\$1,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15}' > ${name}.CovCytoCNVdb.txt
 """
} 



process CALC_COVERAGE_SURVIVOR {
 tag "CALC_COVERAGE VCF on $name using $task.cpus CPUs $task.memory"
 //publishDir  "${params.outDir}/${name}/nano/VarCal/survivor/debug", mode:'copy'

 input:
 tuple val(name), val(sample), path(inputvcf), path(bam), path(bai)

 output:
 tuple val(name), path("${name}.CovCytoCNVdb.txt")

 script:
 """
 echo CALC_COVERAGE_SURVIVOR $name
 samtools view -H $bam | grep @SQ | sed 's/@SQ\tSN:\\|LN://g' > genome.txt 

 cat ${inputvcf} | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") && (\$5!~"GL|MT|KI") && (\$3 ~ "BND") {\$1=\$1;print \$1,\$2-1,\$2,\$3,\$5}' | sort -k1,1 -k2,2V | tail -n +2 > fromBNDs.bed 
cat	fromBNDs.bed | awk '{print \$5}' | grep -oE '[A-Za-z0-9]*:[0-9]*' | sed 's/:/\\t/' | awk -v OFS="\\t" '{print \$1,\$2-1,\$2}' > toBNDs.temp
 #add IDs to join in the end
 paste toBNDs.temp <(cat fromBNDs.bed | cut -f 4) > toBNDs.bed
 bedtools map -c 4 -a <(cat toBNDs.bed | sort -k1,1 -k2,2V ) -b $params.GrCh38cytomap -o concat -g genome.txt > toBNDsCyto.temp

 bedtools coverage -sorted -abam fromBNDs.bed -b ${bam} -d -split -g genome.txt > Coverage.txt
 bedtools map -c 4 -a Coverage.txt -b $params.GrCh38cytomap -o concat -g genome.txt > CoverageCytomap.txt
 bedtools map -c 4,4,5 -a CoverageCytomap.txt -b $params.GrCh38CNVdb -o collapse,count,collapse -g genome.txt > from.CovCytoCNVdb.txt

 join -1 4 -2 4 -t \$'\\t' <(sort -k4 from.CovCytoCNVdb.txt) <(sort -k4 toBNDsCyto.temp) | awk -v FS="\\t" -v OFS="\\t" '{print \$2,\$3,\$4,\$1,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15}' > ${name}.CovCytoCNVdb.txt
 """
} 

process PARSE_SVIM_VCF {
 tag "PARSE_SVIM_VCF VCF on $name using $task.cpus CPUs $task.memory"
 publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
 label "small_process"

 input:
 tuple val(name), val(sample), path(inputvcf), path(bam), path(bai), path(covFile)

 output:
 tuple val(name), val(sample), path("Dedup.1000dfilt.${name}.tsv")
 path('*.tsv')

 script:
 """
 echo PARSE_SVIM_VCF $name
 python ${params.ComputeDistance} ${inputvcf} ${covFile} ${name}
 """
} 

process SURVIVOR {
 tag "SURVIVOR VCF on $name using $task.cpus CPUs $task.memory"
 //publishDir  "${params.outDir}/${name}/nano/VarCal/survivor", mode:'copy'
 container "registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:55"
 label "small_process"

 input:
 tuple val(name), val(sample), path(vcf1), path(vcf2)

 output:
 tuple val(name), val(sample), path("${name}.SURVIVOR.merged.vcf")

 script:
 """
 echo SURVIVOR $name
 ls *.vcf > sample_files
 SURVIVOR merge sample_files 100 1 1 1 0 30 ${name}.SURVIVOR.merged.vcf
 """
} 