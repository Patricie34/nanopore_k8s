include { REFORMAT_SAMPLE; REFORMAT_PARAMS } from "${params.projectDirectory}/modules/input-reformators"
include { GUPPY } from "${params.projectDirectory}/modules/basecallers"
include { COLLECT_FQs; COLLECT_BAMs } from "${params.projectDirectory}/modules/collectors"
include { MINIMAP2 } from "${params.projectDirectory}/modules/mappers"
include { QUALIMAP } from "${params.projectDirectory}/modules/qc"
include { SVIM } from "${params.projectDirectory}/modules/varcallers"
include { CNVkit } from "${params.projectDirectory}/modules/cnv"


workflow {
runlist = channel.fromList(params.samplesheet)
references = channel.fromList(params.ref_specific)
samplesReformated = REFORMAT_SAMPLE(runlist)
referencesReformated = REFORMAT_PARAMS(references)

samplesWithReferences = samplesReformated.combine(referencesReformated,by:1)
	.map({ sample -> [sample[1],sample[2]+(sample[3]) ]})
	// .view{"____________samplesWithReferences_____________: $it"}

FQcalled = GUPPY(runlist)[0]
FQcollected = COLLECT_FQs(runlist)
FQs = FQcalled.mix(FQcollected) //combine fq from different sources

//FQs.join(samplesWithReferences).view {"____________Combined FastQs_____________: $it"}
BAMmapped	= MINIMAP2(FQs.join(samplesWithReferences))
 BAMcollected = COLLECT_BAMs(samplesWithReferences)
	// .view{"____________Collected BAMS_____________: $it"}

BAMs = BAMcollected.mix(BAMmapped)
// .view{"____________All BAMS_____________: $it"}
// BAMs.view()

QUALIMAP(BAMs)

Vcfs = SVIM(BAMs)
Cnvs = CNVkit(BAMs)
Vcf_paths = Vcfs[1].map({it -> [it[2]]})
Combined_collected_vcf = Vcfs[1].combine(Vcf_paths.collect().map({it -> [it]})) //
Combined_filtered = Combined_collected_vcf.map({
	// get one vcf and all other vcfs
	row ->
	def name = row[0]
	def sample = row[1]
	def vcf = row[2]
	def filtered  = removeSame(vcf, row[2])
	[name,sample, vcf, filtered]	
	})

Combined_filtered.view()
//  Tagged = TAG_UNIQUE_VARS(Combined_filtered)
//  Annotated = ANNOTATE(Tagged)
//  BAMs_annot = Annotated.join(BAMs) //.view() // join BAMs with annotated, vcf fits bam for given sample
//  BAMs_vcfs = BAMs.join(Vcfs[1]) //.view() // join BAMs with variants, vcf fits bam for given sample
//  Editedvcfs = EDITVCF(BAMs_annot)
//  Cnvs = CNVKIT(BAMs)
//  QUALIMAP(BAMs)
// Synchronized = Vcfs[0].join(Cnvs[1]).join(Editedvcfs[0])
// CIRCOS(Synchronized)
// HEATMAP(Editedvcfs[0])

// Prefiltered = ASSEMBLY_PREFILTER(BAMs_BAMmapped,BAMs_bais,Vcfs[1])
// Assembly = SHASTA(Prefiltered[0])

// assembly = SHASTA(FQs)
// MUMMER(Assembly[1])

// ///////////////////////////////
//  spolecne(FQs)
// rawfastq	= GUPPY_BASECALL()
// FLYE(rawfastq.map(rawfastq)
// BAMs		= MAPPING(rawfastq.map({ file -> [run, file]}))
// rawfastq.collect().view()
// assembly = SHASTA(rawfastq.collect(),BAMs[0])
}


// def removeSame(nm,lst) {

// 	def list=[]
//         for (int i = 0; i < lst.size(); i++) {
//          if (lst[i] != nm){
// 	     list.add(lst[i])
//        // lst.remove(i)
// 	  }
//         }
//         return(list)
//   }


// 	case "$sample.reference" in
//   "GrCh38")
//     echo "GrCh38 reference selected."
//     	minimap2 --MD -a -t ${task.cpus} ${params.GrCh38ref} $reads > ${name}.BAMs.sam
//     ;;
//   "GrCh37")
//     echo "GrCh37 reference selected."
//    	minimap2 --MD -a -t ${task.cpus} ${params.GrCh38ref} $reads > ${name}.BAMs.sam
//     ;;
//   "T2T")
//     echo "T2T reference selected."
//     # Commands for Option C
//     ;;
//   *)
//     echo "Invalid option selected."
//     exit 1
//     ;;
// esac