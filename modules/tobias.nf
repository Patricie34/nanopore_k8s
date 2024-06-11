process DELLY_SVs {
	tag "DELLY_SVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/SVs/", mode:'copy'
	container "dellytools/delly:latest"
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.delly.bcf")
	
	script:
	"""
	delly lr -y ont  -o ${name}.delly.bcf -g ${sample.ref} ${bam}
	"""
}


process INDEX_BCF {
	tag "INDEX_BCF on $sample.name using $task.cpus CPUs $task.memory"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(name), val(sample), path(bcf)

	output:
	tuple val(name), val(sample), path("${sample.name}.Delly.Chr.bcf"), path("*.csi")
	
	script:
	"""
	bcftools view $bcf | awk '{ if(\$0 !~ /^#/) print "chr"\$0; else if(match(\$0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print \$0 }' | bcftools view - -o ${sample.name}.Delly.Chr.bcf
	bcftools index ${sample.name}.Delly.Chr.bcf
	"""
} 


process INDEX_POP {
	tag "INDEX_POP using $task.cpus CPUs $task.memory"
	label "s_cpu"
	label "s_mem"

	input:
	path popBcf

	output:
	tuple path(popBcf), path("*.csi")
	
	script:
	"""
	bcftools index $popBcf
	"""
} 



process SANSA_MARKDUP {
	tag "SANSA_MARKDUP using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/Delly/", mode:'copy'
	container "dellytools/sansa:v0.2.1"
	label "s_cpu"
	label "xl_mem"

	input:
	tuple val(popBcf), val(popCsi)

	output:
	tuple path("rmdup.bcf"), path("rmdup.bcf.csi")
	
	script:
	"""
	sansa markdup -o rmdup.bcf ${popBcf}
	"""
}

process VCFsPerSample {
	tag "VCFsPerSample on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${params.outDir}/${sample.name}/nano/VarCal/Delly/SVs/", mode:'copy'
	// container "broadinstitute/gatk:4.5.0.0"
	label "s_cpu"
	label "m_mem"
	
	input:
	tuple val(sample), path(mergedBcf), path(mergedBcfCsi)

	output:
	tuple val(sample), path("${sample.name}.markDup.Chr.bcf"), path("${sample.name}.markDup.Chr.bcf.csi")

	script:
	"""
	bcftools view -s ${sample.name} $mergedBcf | bcftools view -i 'GT="1/1" | GT="0/1"' -O b -o ${sample.name}.markDup.bcf
	bcftools view ${sample.name}.markDup.bcf| awk '{ if(\$0 !~ /^#/) print "chr"\$0; else if(match(\$0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print \$0 }' | bcftools view - -o ${sample.name}.markDup.Chr.bcf
	bcftools index ${sample.name}.markDup.Chr.bcf
	"""
}


process DELLY_merge_SVs {
	tag "DELLY_merge_SVs using $task.cpus CPUs $task.memory"
	container "dellytools/delly:latest"
	label "s_cpu"
	label "s_mem"

	input:
	path bcfs_filter

	output:
	path "*"
	
	script:
	"""
	delly merge -n 250000000 -o delly.merged.bcf ${bcfs_filter}
	"""
} 



process DELLY_genotype_SVs {
	tag "DELLY_genotype_SVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/SVs/", mode:'copy'
	container "dellytools/delly:latest"
	label "s_cpu"
	label "l_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai), path(bcf), path(bcfMerged)

	output:
	tuple val(name), val(sample), path("${name}.geno.bcf")
	
	script:
	""" 
	delly call -g ${sample.ref} -v ${bcfMerged} -o ${name}.geno.bcf -x ${sample.nonMappableRepeats} ${bam}
	"""
} 



process DELLY_rename_BCFs {
	tag "DELLY_rename_BCFs on $sample.name using $task.cpus CPUs $task.memory"
	// publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/SVs/", mode:'copy'
	container "staphb/bcftools:1.18"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(name), val(sample), path(bcf)

	output:
	tuple val(name), val(sample), path("${sample.name}.genotyped.renamed.bcf")
	
	script:
	""" 
	bcftools query -l $bcf | tr '\\n' ' ' > RenamePair.txt
	echo -n $sample.name >> RenamePair.txt
	cat RenamePair.txt
	bcftools reheader -s RenamePair.txt -o ${sample.name}.genotyped.renamed.bcf $bcf

	"""
} 

process DELLY_merge_genotyped {
	tag "DELLY_merge_genotyped using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/Delly/", mode:'copy'
	container "staphb/bcftools:1.18"
	label "s_cpu"
	label "s_mem"

	input:
	path genotyped_bcfs

	output:
	path "delly.geno.merged.bcf"
	
	script:
	"""
	# List of files
	file_list=($genotyped_bcfs)

	# Loop over files
	for file in "\${file_list[@]}"; do
    	echo "\$file"
		bcftools index \$file
	done
	bcftools merge -m id -O b -o delly.geno.merged.bcf $genotyped_bcfs
	"""
} 



process SANSA_FILTER_1000G {
	tag "SANSA_FILTER_1000G on $sample.name using $task.cpus CPUs $task.memory"
	// publishDir  "${params.outDir}/${sample.name}/nano/VarCal/Delly/SVs/", mode:'copy'
	container "dellytools/sansa:v0.2.1"
	label "s_cpu"
	label "xl_mem"

	input:
	tuple val(sample), path(bcf), path(csi)

	output:
	tuple val(sample), path("${sample.name}.1000Genomes.tsv")
	
	script:
	"""
	sansa compvcf -a /mnt/shared_resources/homo_sapiens/GRCh38/ONT/all.delly.hg38.1kGP.ont.bcf -e 0 -m 0 -p -n 250000000 $bcf -o $sample.name
	head -n 1 ${sample.name}.sv.classification > ${sample.name}.1000Genomes.tsv
	cat ${sample.name}.sv.classification | awk '\$2=="FP"' >> ${sample.name}.1000Genomes.tsv
	"""
} 

process PARSE_SANSA {
	tag "PARSE_SANSA on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDir}/${sample.name}/nano/VarCal/Delly/SVs/", mode:'copy'
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(sample), path(bcf), path(csi), path(tsv)

	output:
	tuple val(sample), path("${sample.name}.1000G.tsv")
	
	script:
	"""
	bcftools view $bcf > tmp.vcf
	python $params.Edit1kONT -v tmp.vcf -t $tsv > ${sample.name}.1000G.tsv
	"""
} 

process DELLY_filter_merge_genotyped {
	tag "DELLY_filter_merge_genotyped using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/Delly/", mode:'copy'
	container "quay.io/biocontainers/cyvcf2:0.30.25--py38he403ad2_0"
	label "s_cpu"
	label "s_mem"

	input:
	path merged_genotyped_bcfs

	output:
	path "filtered.geno.merged.tsv"
	
	script:
	"""
	python $params.FilterSingletons -v $merged_genotyped_bcfs > filtered.geno.merged.tsv
	"""
} 

process DELLY_get_genesBNDs {
	tag "DELLY_get_genesBNDs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/SVs/", mode:'copy'
	container "staphb/bcftools:1.18"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(name), val(sample), path(bcf)

	output:
	tuple val(name), val(sample), path("${name}.input.tsv")
	
	script:
	""" 
	bcftools query -f "%CHROM\\t%POS\\t%CHROM\\t%INFO/END\\t%ID\\n" ${bcf} | grep -v "BND" > ${name}.sv.tsv 
	bcftools query -f "%CHROM\\t%POS\\t%INFO/CHR2\\t%INFO/POS2\\t%ID\\n" ${bcf} | grep "BND" >> ${name}.sv.tsv
	cat ${name}.sv.tsv | cut -f 1,2,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Left";}' > ${name}.input.tsv
	cat ${name}.sv.tsv | cut -f 3,4,5 | sed 's/^chr//' | awk '{print \$1"\\t"(\$2-100)"\\t"(\$2+100)"\\t"\$3"Right";}' >> ${name}.input.tsv
	"""
} 

process DELLY_annot_genesBNDs {
	tag "DELLY_annot_genesBNDs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/SVs/", mode:'copy'
	container "trausch/alfred:latest"
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(name), val(sample), path(input_tsv)

	output:
	tuple val(name), val(sample), path("${name}.sv.gene.tsv")
	
	script:
	""" 
	alfred annotate -d 3000 -g $sample.alfredAnnotateBed -o ${name}.sv.gene.tsv ${input_tsv}
	"""
} 

process DELLY_CNVs {
	tag "DELLY_CNVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/CNVs/", mode:'copy'
	container "dellytools/delly:latest"
	label "s_cpu"
	label "m_mem"


	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.cov.gz")
	path "*"
	
	script:
	"""
	# -m parameter should not be hardcoded
	delly cnv -w 25000 -i 25000 -j 25000 -a -g ${sample.ref} -m ${sample.nonMappableRepeatsTobias} -c ${name}.cov.gz -o ${name}.cnv.bcf ${bam}
	"""
}

	process DELLY_CNV_PLOT {
	tag "DELLY_CNVs on $sample.name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/Delly/CNVs/", mode:'copy'
	label "s_cpu"
	label "s_mem"

	input:
	tuple val(name), val(sample), path(covfile)

	output:
	path "*"
	
	script:
	"""
	Rscript --vanilla Rscript ${params.CnvPlotTobias} $covfile
	"""
} 