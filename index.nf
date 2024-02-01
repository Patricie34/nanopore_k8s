include { QUALIMAP } from "${params.projectDirectory}/modules/qc"
include { CNVkit } from "${params.projectDirectory}/modules/cnv"

include { wf_process_till_bam } from "${params.projectDirectory}/subworkflows/wf_process_till_bam"
include { wf_svim } from "${params.projectDirectory}/subworkflows/wf_svim"
include { wf_delly } from "${params.projectDirectory}/subworkflows/wf_delly"
include { wf_survivor } from "${params.projectDirectory}/subworkflows/wf_survivor"



workflow {
	///  GATHERING INPUTS
	Runlist = channel.fromList(samplesheet)
	References = channel.fromList(reference_sample_specific)
 
	BAMs = wf_process_till_bam(Runlist, References)
 Svim_annot_vcfs = wf_svim(BAMs)
	wf_delly(BAMs)
	wf_survivor(BAMs, Svim_annot_vcfs)
 


	// 	QUALIMAP(BAMs)
	// Cnvs = CNVkit(BAMs)



	///  GATHERING INPUTS
//*****************************************************
	// runlist = channel.fromList(params.samplesheet)
	// references = channel.fromList(params.ref_specific)
	// samplesReformated = REFORMAT_SAMPLE(runlist)
	// referencesReformated = REFORMAT_PARAMS(references)
	// samplesWithReferences = samplesReformated.combine(referencesReformated,by:1) //match parameters to individual samples
	// 	.map({ sample -> [sample[1],sample[2]+(sample[3]) ]})
//*******************************************************
/// CALLING, MAPPING
//*******************************************************
	// FQcalled = GUPPY(runlist)[0]
	// FQcollected = COLLECT_FQs(runlist)
	// FQs = FQcalled.mix(FQcollected) //combine fq from different sources
	// BAMmapped	= MINIMAP2(FQs.join(samplesWithReferences))
	// BAMcollected = COLLECT_BAMs(samplesWithReferences)
	// BAMs = BAMcollected.mix(BAMmapped)
// .view{"____________All BAMS_____________: $it"}
// BAMs.view()
//*******************************************************
/// SVIM, CNVkit, QC...
//*******************************************************
// 	QUALIMAP(BAMs)
// 	Vcfs = SVIM(BAMs)
// 	Cnvs = CNVkit(BAMs)
// 	// Vcf_paths = Vcfs[1].map({it -> [it[2]]})
// //*******************************************************
// /// ANNOTATING, CUSTOM VCFFILTERING...
// //*******************************************************
// 	Combined_collected_vcf = Vcfs[1].combine(Vcf_paths.collect().map({it -> [it]})) //
// 	Combined_filtered = Combined_collected_vcf.map({
// 		row -> // get one vcf and all other vcfs
// 		def name = row[0]
// 		def sample = row[1]
// 		def vcf = row[2]
// 		def filtered  = removeSame(vcf, row[3])
// 		[name,sample, vcf, filtered]	
// 		})

// 		Tagged = TAG_UNIQUE_VARS(Combined_filtered)
// 		Annotated = ANNOTATE(Tagged)

// 		SVIM_Annotated_Sup2 = FILTER_SVIM(Annotated)

// 		BAMs_annot = SVIM_Annotated_Sup2.join(BAMs).map({
// 			// get rid of double sample info from 
// 			row ->
// 			def name = row[0]
// 			def sample = row[1]
// 			def vcf = row[2]
// 			def bam  = row[4]
// 			def bai  = row[5]

// 			[name, sample, vcf, bam, bai]	
// 	})//.view{"____________Annotated+Bams_____________: $it"}// join BAMs with annotated, vcf fits bam for given sample
// 	CovFiles = CALC_COVERAGE(BAMs_annot)
// 	BAMs_annot_cov = BAMs_annot.join(CovFiles)

// 	Editedvcfs = EDITVCF(BAMs_annot_cov) // NOT ALL FINISHED

//  Sniffles_vcfs = SNIFFLES2(BAMs)

//  joinedVCFs = Sniffles_vcfs.join(SVIM_Annotated_Sup2).map({
// 			// get rid of double sample info from 
// 			row ->
// 			def name = row[0]
// 			def sample = row[1]
// 			def vcf1 = row[2]
// 			def vcf2  = row[4]

// 			[name, sample, vcf1, vcf2]	
// 	})

		
// 	Survivor_vcfs = SURVIVOR(joinedVCFs)
// 	Survivor_annot = ANNOTATE_SURVIROR(Survivor_vcfs)

// BAMs_Survivor = Survivor_vcfs.join(BAMs).map({
// 			// get rid of double sample info from 
// 			row ->
// 			def name = row[0]
// 			def sample = row[1]
// 			def vcf = row[2]
// 			def bam  = row[4]
// 			def bai  = row[5]

// 			[name, sample, vcf, bam, bai]	
// 	})//.view{"____________BAMs_Survivor_____________: $it"}// join BAMs with annotated, vcf fits bam for given sample
		// CovFilesSurvivor = CALC_COVERAGE_SURVIVOR(BAMs_Survivor)


	//*******************************************************
/// TOBIAS, DELLY...
//*******************************************************
// Delly_CNVs = DELLY_CNVs(BAMs)
// DELLY_CNV_PLOT(Delly_CNVs[0])

// Primary_bam = PEPPER_PREFILT(BAMs)
// Bam_unique_reads = PEPPER_UMIREADS(BAMs)
//DeepVarOutput = PEPPER_DEEPVARIANT(Bam_unique_reads)

// DeepVarOutput[0].view()
// DeepVarOutput[1].view()

//Phased_VCFs = PHASE_WHATSHAP(BAMs_annot)


// Delly_SVs = DELLY_SVs(BAMs)
// Delly_BC_paths = Delly_SVs.map({it -> [it[2]]})
// //Delly_BC_paths.view{"____________Delly_BC_paths_____________: $it"}
// Delly_SVs_collected = Delly_BC_paths.collect()

// // Delly_SVs_preFilter = Delly_SVs_combined.map({
// // 		row -> // get one vcf and all other vcfs
// // 		def name = row[0]
// // 		def sample = row[1]
// // 		def vcf = row[2]
// // 		def filtered  = removeSame(vcf, row[3])
// // 		[name,sample, vcf, filtered]	
// // 		})//.view{"____________Delly_BC_paths_filtered_____________: $it"}
//  Delly_merged_BCF = DELLY_merge_SVs(Delly_SVs_collected)

// 	Delly_SVs_combined = Delly_SVs.combine(Delly_merged_BCF).join(BAMs).map({
// 			// get rid of double sample info from 
// 			row ->
// 			def name = row[0]
// 			def sample = row[1]
// 			def vcf = row[2]
// 			def mergedVcf = row[3]
// 			def bam  = row[5]
// 			def bai  = row[6]

// 			[name, sample, bam, bai, vcf,mergedVcf ]	
// 	})//.view{"____________Delly_BC_paths_bams_____________: $it"}

// 	Delly_genotyped = DELLY_genotype_SVs(Delly_SVs_combined) // TOTO PAK NEJAK FILTROVAT?
// 	Delly_genotyped_paths = Delly_genotyped.map({it -> [it[2]]})
//  Delly_genotyped_paths_collected = Delly_genotyped_paths.collect()

//  Delly_merged_genotyped = DELLY_merge_genotyped(Delly_genotyped_paths_collected)
// 	DELLY_filter_merge_genotyped(Delly_merged_genotyped)
// 	//then parse back to sample-specific

// //BCF2TVC(Delly_genotyped)
// Delly_BNDs = DELLY_get_genesBNDs(Delly_SVs) //Delly_genotyped
// DELLY_annot_genesBNDs(Delly_BNDs)

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