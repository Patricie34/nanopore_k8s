include {
	DELLY_SVs;
	DELLY_merge_SVs;
	DELLY_genotype_SVs;
	DELLY_CNVs;
	DELLY_CNV_PLOT;
	DELLY_get_genesBNDs;
	DELLY_annot_genesBNDs;
	DELLY_merge_genotyped;
	DELLY_filter_merge_genotyped;
	INDEX_BCF;
	SANSA_FILTER_1000G;
	SANSA_MARKDUP;
	INDEX_POP;
	DELLY_rename_BCFs;
	VCFsPerSample;
	PARSE_SANSA;
	} from  "${params.projectDirectory}/modules/tobias"

workflow wf_delly {
	take:
	BAMs

	main:
//  (Delly_CNVs, _) = DELLY_CNVs(BAMs)
//  DELLY_CNV_PLOT(Delly_CNVs)
	Delly_SVs = DELLY_SVs(BAMs)
	Delly_BC_paths = Delly_SVs.map({it -> [it[2]]})
	Delly_SVs_collected = Delly_BC_paths.collect()
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
	Delly_genotyped = DELLY_genotype_SVs(Delly_SVs_combined)
	Delly_genotyped_renamed = DELLY_rename_BCFs(Delly_genotyped)
	Delly_genotyped_paths = Delly_genotyped_renamed.map({it -> [it[2]]})
	Delly_genotyped_paths_collected = Delly_genotyped_paths.collect()
	Delly_merged_genotyped = DELLY_merge_genotyped(Delly_genotyped_paths_collected)
	Indexed_merged = INDEX_POP(Delly_merged_genotyped)
	MarkDup_Pop_Bcf = SANSA_MARKDUP(Indexed_merged)
	MarkDupBcfs = VCFsPerSample(Delly_SVs_combined.map{ [it[1]]}.combine(MarkDup_Pop_Bcf))//Parse multiBcf to sample-vcfs
	SansaTsvs = SANSA_FILTER_1000G(MarkDupBcfs)
	PARSE_SANSA(MarkDupBcfs.join(SansaTsvs))

	// emit:
	// Delly_SVs
}