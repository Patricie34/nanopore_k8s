// run = "${params.data}".split("/")
// run = run[run.size()-1]
// resultsDir = "${params.outDir}/${run}/nano"



process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'
	container "cerit.io/docker/genomicpariscentre/guppy-gpu"
	accelerator 1, type: 'nvidia.com/gpu'
			label "small_process"

	// accelerator 1, type: 'nvidia.com/mig-1g.10gb'


	input:
	tuple val(name), val(path), val(type), val(config)

	output:
	tuple val(name), path("basecalled/pass/*.gz")
	tuple val(name), path("*")


	when:
		type == 'fast5'

	script:
	"""
	#/opt/ont/guppy/bin/guppy_basecaller -i ${path}/ -s ./basecalled --flowcell ${params.flowcell} --kit ${params.kit} --compress_fastq --recursive --num_callers ${task.cpus} -x auto
	/opt/ont/guppy/bin/guppy_basecaller -i ${path} -s ./basecalled -c ${config} --compress_fastq --recursive --num_callers ${task.cpus} -x auto --progress_stats_frequency 600 --verbose_logs --trace_category_logs --fast5_out
	"""
} 


process COLLECT_BASECALLED {
	tag "Collecting Fastq on $name using $task.cpus CPUs $task.memory"
	//publishDir  "${params.outDir}/${name}/nano/", mode:'copy'

		label "small_process"


	input:
	tuple val(name), val(path), val(type), val(config)

	output:
	tuple val(name), path("basecalled/pass/*.fastq.gz")

	when:
		type == 'fastq'

	script:
	"""
	mkdir -p basecalled/pass/
	find ${path} -type f -name '*.fastq.gz' | xargs -I % cp % ./basecalled/pass/
	"""
} 

process COLLECT_MAPPED {
	tag "Collecting bam on $name using $task.cpus CPUs $task.memory"
	//publishDir  "${params.outDir}/${name}/mapped/", mode:'copy'
		label "small_process"


	input:
	tuple val(name), val(path), val(type), val(config)

	output:
	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bam.bai")

	when:
		type == 'bam'

	script:
	"""
	find ${path} -type f -name '*.ba*' | xargs -I % cp % .
	"""
} 

process MAPPING {
	tag "Mapping on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/mapped/", mode:'copy'

	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bam.bai")
	
	script:
	"""
	minimap2 --MD -a -t ${task.cpus} ${params.ref} $reads > ${name}.mapped.sam
	samtools view -bS ${name}.mapped.sam > ${name}.mapped.bam
	samtools sort -o ${name}.sorted.bam ${name}.mapped.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bam.bai	
	"""
} 

process SVIM{
	tag "Variant calling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'


	input:
	tuple val(name), path(mapped), path(bai)

	output:
	tuple val(name),path ("*")
	tuple val(name), path("${name}.variants.vcf")

	
	script:
	"""
	svim alignment ./ ${mapped} ${params.ref} --minimum_depth 1 --read_names --all_bnds
	cp variants.vcf ${name}.variants.vcf
	rm variants.vcf
	"""
} 

process TAG_UNIQUE_VARS{
	tag "Tagging unique vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	//label "small_process"


	input:
 tuple val(name), path(vcf), path(vcf_filter_with)

	output:
	tuple val(name), path("*UniqueTag.vcf")

	
	script:
	"""
	bedtools intersect -a ${vcf} -b ${vcf_filter_with} -v | awk '{print \$3}' > uniqueIDs.txt
 python ${params.TagUniqes} ${vcf} uniqueIDs.txt ${name}
	"""
} 


process ANNOTATE{
	tag "Annotating vcf vars on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	label "big_mem"


	input:
 tuple val(name), path(vcf)

	output:
	tuple val(name), path("${name}.UniqueTag.Annotated.vcf")

	
	script:
	"""
	vep -i ${vcf} --cache --cache_version 90 --dir_cache ${params.vep} --fasta ${params.ref} --merged --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.UniqueTag.Annotated.vcf
	"""
} 


process EDITVCF {
	tag "Post process VCF on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'

	input:
	tuple val(name), path(inputvcf), path(bam) ,path(bai)

	output:
	tuple val(name), path("Dedup_1000distDistanced_${name}.tsv")
	path('*')
		
	script:
	"""
	#cat ${inputvcf} | grep BND | grep -vE 'GL000|MT|KI27' | awk -v OFS="\\t" '{print \$1,\$2-1,\$2,\$3}' | sort -k1,1 -k2,2V > BNDlocs.bed
	#bedtools coverage causes problems with strange chroms, remove first line
	cat ${inputvcf} | grep BND | awk -v OFS="\\t" '(\$1!~"GL|MT|KI") && (\$5!~"GL|MT|KI") {\$1=\$1;print \$1,\$2-1,\$2,\$3}' | sort -k1,1 -k2,2V | tail -n +2 > BNDlocs.bed 
	samtools view -H $bam | grep @SQ | sed 's/@SQ\tSN:\\|LN://g' > genome.txt 
	head -n 20 BNDlocs.bed
	head -n 20 genome.txt
	bedtools coverage -sorted -abam BNDlocs.bed -b ${bam} -d -g genome.txt > Coverage.txt
	python ${params.ComputeDistance} ${inputvcf} Coverage.txt Distanced_${name}
	"""
} 


process CIRCOS{
	tag "Creating circos on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/Circos/", mode:'copy'

	input:
	tuple val(name),	path(vcfs), path(cnv_sorted),	path(vars_edited)


	output:
	path("*.html")
	
	script:
	"""
	echo $vcfs
	echo $cnv_sorted
	echo $vars_edited
	cat $vars_edited | awk -v FS="\\t" '(\$3!~"GL|MT|KI") && (\$9!~"GL|MT|KI") {print \$2,\$3,\$4,\$9,\$10,\$15}' > vars_slimmed.tsv
	head -n 20 vars_slimmed.tsv
	Rscript --vanilla ${params.circos} $cnv_sorted vars_slimmed.tsv $name ${params.grch38_lens}
	"""
} 


process QUALIMAP {
		tag "Qualimap $name using $task.cpus CPUs $task.memory"

	container 'pegi3s/qualimap'
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'
	label "big_mem"
	
	input:
	tuple val(name), path(bams), path(bais)

	output:
	path "qualimap/*"

	script:
	"""
	qualimap bamqc -bam ${bams} -nw 5000 -nt 14 --java-mem-size=60G -c -outdir ./qualimap
	"""
}

process FLYE {
	container 'staphb/flye'
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'
	label "big_mem"


	input:
	tuple val(name), path (fastq)

	output:
	path "Flye_output/*"

	script:
	"""
	flye --nano-hq ${fastq} --threads ${task.cpus} -m 1 -o ./Flye_output
	"""
}


process ASSEMBLY_PREFILTER {
		tag "Prefiltering BNDs on $name using $task.cpus CPUs $task.memory"

	publishDir  "${params.outDir}/${name}/nano/SHASTA/", mode:'copy'

	input:
	tuple val(name), path(bam), path(bai), path(vars)

	output:
	tuple val(name), path("${name}.BND.fasta")
	path 'readnames.txt'
	path 'filtered.bam'


	script:
	"""
 cat ${vars} | grep -oP '(?<=READS\\=)([\\w\\-]*,)*[\\w\\-]*' | awk  -v FS="," -v OFS="\n" '{\$1=\$1;print}' > readnames.txt
	samtools view ${bam} -b -N readnames.txt > filtered.bam
	samtools fasta filtered.bam > ${name}.BND.fasta
	"""
}


process SHASTA {
	publishDir  "${params.outDir}/${name}/nano/SHASTA/", mode:'copy'
	// label "biggest_mem"
	//label "big_cpus"
	tag "Shasta on $name using $task.cpus CPUs $task.memory"

	input:
	tuple val(name), path(fasta)
	//tuple val(name), path(bam)


	output:
	path "ShastaRun/*"
	tuple val(name), path("ShastaRun/Assembly.fasta")
 
	script:
	"""
	#shasta --help
	#gunzip -d --force *.fastq.gz
	#cat *.fastq* | sed -n '1~4s/^@/>/p;2~4p' > all.fasta
	#cat *.fastq* | samtools fasta - > all.fasta
	#cat all.fasta | wc
#--memoryMode filesystem --memoryBacking disk
	echo Running_Shasta
	cat $fasta | wc -l
	shasta --input $fasta --config Nanopore-May2022 --threads ${task.cpus} --Reads.minReadLength 5000
	"""
}

process CNVKIT {
	tag "CNVkit $name using $task.cpus CPUs $task.memory"
 container 'etal/cnvkit'
	publishDir  "${params.outDir}/${name}/nano/CNVkit", mode:'copy'

	input:
	tuple val(name), path(mapped), path(bai)

	output:
	path "*"
	tuple val(name), path("${name}.sorted.cns")

	script:
	"""
	cnvkit.py batch ${mapped} -n -m wgs -f ${params.ref} --annotate ${params.refFlat} -p $task.cpus --target-avg-size 30000
	cat ${name}.sorted.cns | grep -v GL000 > ${name}.cns
	cat ${name}.sorted.cnr | grep -v GL000 > ${name}.cnr
	cnvkit.py scatter -s ${name}.cn{s,r} -o ${name}.scatter.svg
	cnvkit.py diagram -s ${name}.cn{s,r}
	"""
}


process MUMMER {
 container 'staphb/mummer'
	//container 'quay.io/biocontainers/mummer'
	publishDir  "${params.outDir}/${name}/nano/Dnadiff", mode:'copy'
	label "biggest_mem"
	tag "Mummer on $name using $task.cpus CPUs $task.memory"

	input:
	tuple val(name), path(assembly)

	output:
	path "*"

	script:
	"""
	#dnadiff -p dnadiff -t $task.cpus ${params.ref} ${assembly}
	#	sleep infinity

	nucmer -t $task.cpus ${params.ref} ${assembly}
#run-mummer3 -t $task.cpus ${params.ref} ${assembly}
#mummerplot out.delta
	"""
}

workflow {
//println("${params.data}")
//println("Run name: ${run}")
runlist = channel.fromList(params.runs)
FQcalled = GUPPY_BASECALL(runlist)
FQcollected = COLLECT_BASECALLED(runlist)
FQfiles = FQcalled[0].mix(FQcollected)

Bamcollected = COLLECT_MAPPED(runlist)
Bams	= MAPPING(FQfiles)
mapped_bams = Bamcollected[0].mix(Bams[0]) // combine all bams from differents sources
mapped_bais = Bamcollected[1].mix(Bams[1])
mapped = mapped_bams.join(mapped_bais)

Vcfs = SVIM(mapped)
Vcf_paths = Vcfs[1].map({it -> [it[1]]})
Combined_collected_vcf = Vcfs[1].combine(Vcf_paths.collect().map({it -> [it]}))
Combined_filtered = Combined_collected_vcf.map({
	 row ->
            def name      = row[0]
												def vcf       = row[1]
            def filtered    = removeSame(vcf, row[2])
            [name,vcf, filtered]	
	})

 Tagged = TAG_UNIQUE_VARS(Combined_filtered)
 Annotated = ANNOTATE(Tagged)


// //  Annotated = ANNOTATE(Vcfs[1])
 Mapped_annot = Annotated.join(mapped) //.view() // join mapped with annotated, vcf fits bam for given sample
 Mapped_vcfs = mapped.join(Vcfs[1]) //.view() // join mapped with variants, vcf fits bam for given sample
 Editedvcfs = EDITVCF(Mapped_annot)
 Cnvs = CNVKIT(mapped)
 QUALIMAP(mapped)

// Prefiltered = ASSEMBLY_PREFILTER(Mapped_vcfs)
// Vcfs[1].view()
// Cnvs[0].view()
// Editedvcfs[0].view() 

Synchronized = Vcfs[0].join(Cnvs[1]).join(Editedvcfs[0])
// SynchronizedParsed = Synchronized.map({
// row -> def name = row[0]
// 							def vcfs = row[1]
// 							def cnvs = row[2]
// 							def editedVcf = row[3]
// 							[name,vcfs,cnvs,editedVcf]
// })
// SynchronizedParsed.view()
CIRCOS(Synchronized)
// CIRCOS(Synchronized[0],Synchronized[1][0], Synchronized[1][1], Synchronized[1][2])

// CIRCOS(Vcfs[0], Cnvs[1], Editedvcfs[0])

//Prefiltered = ASSEMBLY_PREFILTER(mapped_bams,mapped_bais,Vcfs[1])
//Assembly = SHASTA(Prefiltered[0])

//assembly = SHASTA(FQfiles)
//MUMMER(Assembly[1])

/////////////////////////////////
 //spolecne(FQfiles)
//rawfastq	= GUPPY_BASECALL()
//FLYE(rawfastq.map(rawfastq)
//mapped		= MAPPING(rawfastq.map({ file -> [run, file]}))
// rawfastq.collect().view()
// assembly = SHASTA(rawfastq.collect(),mapped[0])
}


def removeSame(nm,lst) {

	def list=[]
        for (int i = 0; i < lst.size(); i++) {
         if (lst[i] != nm){
	     list.add(lst[i])
       // lst.remove(i)
	  }
        }
        return(list)
  }
