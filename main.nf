// run = "${params.data}".split("/")
// run = run[run.size()-1]
// resultsDir = "${params.outDir}/${run}/nano"



process GUPPY_BASECALL {
	tag "Basecalling on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'
	accelerator 1, type: 'nvidia.com/mig-1g.10gb'

	input:
	tuple val(name), val(path), val(type)

	output:
	tuple val(name), path("basecalled/pass/*.gz")

	when:
		type == 'fast5'

	script:
	"""
	guppy_basecaller -i ${path}/ -s ./basecalled --flowcell ${params.flowcell} --kit ${params.kit} --compress_fastq --recursive --num_callers ${task.cpus} -x auto
	"""
} 


process COLLECT_BASECALLED {
	tag "Collecting Fastq on $name using $task.cpus CPUs $task.memory"
	//publishDir  "${params.outDir}/${name}/nano/", mode:'copy'

	input:
	tuple val(name), val(path), val(type)

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

	input:
	tuple val(name), val(path), val(type)

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
	tuple val(name), path(mapped)
	tuple val(name), path(bai)

	output:
	path '*'
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
	label "small_process"


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
	label "small_process"


	input:
 tuple val(name), path(vcf)

	output:
	tuple val(name), path("${name}.UniqueTag.Annotated.vcf")

	
	script:
	"""
	vep -i ${vcf} --cache --cache_version 95 --dir_cache ${params.vep} 	--fasta ${params.ref} --merged --offline --format vcf --vcf --everything --canonical --force_overwrite -o ${name}.UniqueTag.Annotated.vcf
	"""
} 


process EDITVCF {
	tag "Post process VCF on $name using $task.cpus CPUs $task.memory"
	publishDir  "${params.outDir}/${name}/nano/VarCal/", mode:'copy'
	label "small_process"

	input:
	tuple val(name), path(bam)
	tuple val(name), path(bai)
	tuple val(name), path(inputvcf)

	output:
	tuple val(name), path("*.tsv")
	
	script:
	"""
 cat ${inputvcf} | grep BND | awk -v OFS="\\t" '{print \$1,\$2-1,\$2,\$3}' > BNDlocs.bed
	bedtools coverage -abam BNDlocs.bed -b ${bam} -d > Coverage.txt
	python ${params.ComputeDistance} ${inputvcf} Coverage.txt Distanced_${name}
	"""
} 



process QUALIMAP {
	container 'pegi3s/qualimap'
	publishDir  "${params.outDir}/${name}/nano/", mode:'copy'
	label "big_mem"
	
	input:
	tuple val(name), path (bams)

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
	publishDir  "${params.outDir}/${name}/SHASTA/", mode:'copy'

	input:
	tuple val(name), path (fastq)
	tuple val(name), path (bam)


	output:
	path "*.txt"
	path "${name}.BND.fasta"

	script:
	"""
 
	samtools view ${bam} -b -N ${name}.BNDreadnames.txt | samtools fasta - > ${name}.BND.fasta
	"""
}


process SHASTA {
	publishDir  "${params.outDir}/${name}/nano/SHASTA/", mode:'copy'
	memory '300 GB'
	label "biggest_mem"
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
	shasta --input $fasta --config Nanopore-May2022 --threads ${task.cpus} --Reads.minReadLength 5000
	"""
}

process CNVKIT {
 container 'etal/cnvkit'
	publishDir  "${params.outDir}/${name}/nano/CNVkit", mode:'copy'

	input:
	tuple val(name), path(mapped)
	tuple val(name), path(bai)

	output:
	path "*"

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
	nucmer -t $task.cpus ${params.ref} ${assembly}
	run-mummer3 -t $task.cpus ${params.ref} ${assembly}
	mummerplot out.delta
	"""
}

workflow {
//println("${params.data}")
//println("Run name: ${run}")
runlist = channel.fromList(params.runs)
FQcalled = GUPPY_BASECALL(runlist)
FQcollected = COLLECT_BASECALLED(runlist)
FQfiles = FQcalled.mix(FQcollected)
Bamcollected = COLLECT_MAPPED(runlist)
Bamcollected[0].view()

Bams	= MAPPING(FQfiles)
Bams[0].view()
mapped_bams = Bamcollected[0].mix(Bams[0])
mapped_bais = Bamcollected[1].mix(Bams[1])

 Vcfs = SVIM(mapped_bams,mapped_bais)

 Vcf_paths = Vcfs[1].map({it -> [it[1]]})



// vcf_paths =  Vcfs[1].collect()
// 													.map({
//   															it -> [it[0]]
// 																	}).view()

 Combined_collected_vcf = Vcfs[1].combine(
																																													Vcf_paths.collect().map({it -> [it]})
																																													)

Combined_filtered = Combined_collected_vcf.map({
	 row ->
            def name      = row[0]
												def vcf       = row[1]
												// println(name)
            def filtered    = removeSame(vcf, row[2])
            [name,vcf, filtered]	
												// [name]
	})


Combined_filtered.view()

Tagged = TAG_UNIQUE_VARS(Combined_filtered)

Annotated = ANNOTATE(Tagged)


editedvcfs = EDITVCF(mapped_bams,mapped_bais,Annotated)
//CNVKIT(mapped[0],mapped[1])
//QUALIMAP(mapped[0])

//prefiltered = ASSEMBLY_PREFILTER(editedvcfs, mapped[0])

//assembly = SHASTA(prefiltered[1])

//assembly = SHASTA(FQfiles)
//MUMMER(assembly[1])

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
