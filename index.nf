include { REFORMAT_SAMPLE; REFORMAT_PARAMS } from "${params.projectDirectory}/modules/input-reformators"
include { GUPPY } from "${params.projectDirectory}/modules/basecallers"
include { COLLECT_FQs; COLLECT_BAMs } from "${params.projectDirectory}/modules/collectors"
include { MINIMAP2 } from "${params.projectDirectory}/modules/mappers"
include { QUALIMAP } from "${params.projectDirectory}/modules/qc"
include { SVIM; SNIFFLES2; PEPPER_DEEPVARIANT } from "${params.projectDirectory}/modules/varcallers"
include { CNVkit } from "${params.projectDirectory}/modules/cnv"
include {TAG_UNIQUE_VARS; ANNOTATE; EDITVCF; CALC_COVERAGE; PHASE_WHATSHAP; BCF2TVC; FILTER_SVIM; PEPPER_PREFILT} from  "${params.projectDirectory}/modules/vcf-processors"
include {DELLY_SVs; DELLY_merge_SVs; DELLY_genotype_SVs;DELLY_CNVs;DELLY_CNV_PLOT;DELLY_get_genesBNDs;DELLY_annot_genesBNDs} from  "${params.projectDirectory}/modules/tobias"



workflow {
	///  GATHERING INPUTS
//*****************************************************
	runlist = channel.fromList(params.samplesheet)
	references = channel.fromList(params.ref_specific)
	samplesReformated = REFORMAT_SAMPLE(runlist)
	referencesReformated = REFORMAT_PARAMS(references)
	samplesWithReferences = samplesReformated.combine(referencesReformated,by:1) //match parameters to individual samples
		.map({ sample -> [sample[1],sample[2]+(sample[3]) ]})
//*******************************************************
/// CALLING, MAPPING
//*******************************************************
	FQcalled = GUPPY(runlist)[0]
	FQcollected = COLLECT_FQs(runlist)
	FQs = FQcalled.mix(FQcollected) //combine fq from different sources
	BAMmapped	= MINIMAP2(FQs.join(samplesWithReferences))
	BAMcollected = COLLECT_BAMs(samplesWithReferences)
	BAMs = BAMcollected.mix(BAMmapped)
// .view{"____________All BAMS_____________: $it"}
// BAMs.view()
//*******************************************************
/// SVIM, CNVkit, QC...
//*******************************************************
	QUALIMAP(BAMs)
	Vcfs = SVIM(BAMs)
	Cnvs = CNVkit(BAMs)
	Vcf_paths = Vcfs[1].map({it -> [it[2]]})
//*******************************************************
/// ANNOTATING, CUSTOM VCFFILTERING...
//*******************************************************
	Combined_collected_vcf = Vcfs[1].combine(Vcf_paths.collect().map({it -> [it]})) //
	Combined_filtered = Combined_collected_vcf.map({
		row -> // get one vcf and all other vcfs
		def name = row[0]
		def sample = row[1]
		def vcf = row[2]
		def filtered  = removeSame(vcf, row[3])
		[name,sample, vcf, filtered]	
		})

		Tagged = TAG_UNIQUE_VARS(Combined_filtered)
		Annotated = ANNOTATE(Tagged)

		SVIM_Annotated_Sup2 = FILTER_SVIM(Annotated)

		BAMs_annot = SVIM_Annotated_Sup2.join(BAMs).map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1]
			def vcf = row[2]
			def bam  = row[4]
			def bai  = row[5]

			[name, sample, vcf, bam, bai]	
	})//.view{"____________Annotated+Bams_____________: $it"}// join BAMs with annotated, vcf fits bam for given sample
	CovFiles = CALC_COVERAGE(BAMs_annot)
	BAMs_annot_cov = BAMs_annot.join(CovFiles)


	Editedvcfs = EDITVCF(BAMs_annot_cov) // NOT ALL FINISHED

 Sniffles_vcfs = SNIFFLES2(BAMs)

	//*******************************************************
/// TOBIAS, DELLY...
//*******************************************************
Delly_CNVs = DELLY_CNVs(BAMs)
DELLY_CNV_PLOT(Delly_CNVs[0])

//Primary_bam = PEPPER_PREFILT(BAMs)
//DeepVarOutput = PEPPER_DEEPVARIANT(Primary_bam)
// DeepVarOutput[0].view()
// DeepVarOutput[1].view()

//Phased_VCFs = PHASE_WHATSHAP(BAMs_annot)
Delly_SVs = DELLY_SVs(BAMs)
Delly_BC_paths = Delly_SVs.map({it -> [it[2]]})
//Delly_BC_paths.view{"____________Delly_BC_paths_____________: $it"}
Delly_SVs_collected = Delly_BC_paths.collect()

// Delly_SVs_preFilter = Delly_SVs_combined.map({
// 		row -> // get one vcf and all other vcfs
// 		def name = row[0]
// 		def sample = row[1]
// 		def vcf = row[2]
// 		def filtered  = removeSame(vcf, row[3])
// 		[name,sample, vcf, filtered]	
// 		})//.view{"____________Delly_BC_paths_filtered_____________: $it"}
 Delly_merged_BCF = DELLY_merge_SVs(Delly_SVs_collected)

	Delly_SVs_combined = Delly_SVs.combine(Delly_merged_BCF).join(BAMs).map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1]
			def vcf = row[2]
			def mergedVcf = row[3]
			def bam  = row[5]
			def bai  = row[6]

			[name, sample, bam, bai, vcf,mergedVcf ]	
	})//.view{"____________Delly_BC_paths_bams_____________: $it"}

	//Delly_genotyped = DELLY_genotype_SVs(Delly_SVs_combined) // TOTO PAK NEJAK FILTROVAT?
//BCF2TVC(Delly_genotyped)
Delly_BNDs = DELLY_get_genesBNDs(Delly_SVs) //Delly_genotyped
DELLY_annot_genesBNDs(Delly_BNDs)

}

//*******************************************************
/// CUSTOM FILTERING SCRIPTS
//*******************************************************
def removeSame(nm,lst) {
	def list=[]
 for (int i = 0; i < lst.size(); i++) {
   if (lst[i] != nm){
	   list.add(lst[i])}
   }
 return(list)
}



//*******************************************************
/// OLD PARTS
//*******************************************************
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