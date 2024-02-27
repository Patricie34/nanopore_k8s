include {DELLY_SVs; DELLY_merge_SVs; DELLY_genotype_SVs;DELLY_CNVs;DELLY_CNV_PLOT;DELLY_get_genesBNDs;DELLY_annot_genesBNDs;DELLY_merge_genotyped;DELLY_filter_merge_genotyped} from  "${params.projectDirectory}/modules/tobias"

workflow wf_delly {
 take:
 BAMs

	main:
 (Delly_CNVs, _) = DELLY_CNVs(BAMs)
 DELLY_CNV_PLOT(Delly_CNVs)
 Delly_SVs = DELLY_SVs(BAMs)
//  Delly_BC_paths = Delly_SVs.map({it -> [it[2]]})
//  Delly_SVs_collected = Delly_BC_paths.collect()
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

	emit:
	//Delly_merged_genotyped
	Delly_SVs
}