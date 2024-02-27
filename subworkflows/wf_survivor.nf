include { SNIFFLES2; PEPPER_DEEPVARIANT } from "${params.projectDirectory}/modules/varcallers"
include { PEPPER_PREFILT; PEPPER_UMIREADS } from "${params.projectDirectory}/modules/vcf-processors"
include { ANNOTATE_SURVIROR; CALC_COVERAGE_SURVIVOR;  SURVIVOR; APPEND_SUPP_FIELD; SURVIVOR_INTERSECT_SAMPLES; SURVIVOR_FILTER_SINGLETONS; SURVIVOR_INFLATE_BNDs; PARSE_SURVIVOR} from "${params.projectDirectory}/modules/survivor"

workflow wf_survivor {
 take:
 BAMs
	Svim_vcfs
	Delly_vcfs

	main:
 Sniffles_vcfs = SNIFFLES2(BAMs)
	// Primary_bam = PEPPER_PREFILT(BAMs)
 // Bam_unique_reads = PEPPER_UMIREADS(BAMs)
	svim = Channel.of( 'svim' )
	sniffles = Channel.of( 'sniffles' )
	delly = Channel.of( 'delly' )


	Mixed_vcfs = Svim_vcfs.combine(svim).mix(Sniffles_vcfs.combine(sniffles), Delly_vcfs.combine(delly))

	Inflated_vcfs = SURVIVOR_INFLATE_BNDs(Mixed_vcfs)//.view{"Inflated_vcfs: $it"}
	Inflated_vcfs_perSample =	Inflated_vcfs.groupTuple().map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1][0]
			def vcfs = row[2]

			def vcf1 = row[2][0]
			def vcf2  = row[2][1]
			def vcf3  = row[2][2]

			[name, sample, vcfs]	
	})//.view{"Inflated_vcfs_perSample: $it"}


	All_joinedVCFs = Sniffles_vcfs.join(Svim_vcfs).join(Delly_vcfs).map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1]
			def vcf1 = row[2]
			def vcf2  = row[4]
			def vcf3  = row[6]

			[name, sample, vcf1, vcf2, vcf3]	
	})//.view{"All_joinedVCFs: $it"}
	// joinedVCFs = Sniffles_vcfs.join(Svim_vcfs).map({
	// 		// get rid of double sample info from 
	// 		row ->
	// 		def name = row[0]
	// 		def sample = row[1]
	// 		def vcf1 = row[2]
	// 		def vcf2  = row[4]

	// 		[name, sample, vcf1, vcf2]	
	// })


	Survivor_vcfs = SURVIVOR(Inflated_vcfs_perSample)//.view{"Survivor_vcfs: $it"}
	// Survivor_annot = ANNOTATE_SURVIROR(Survivor_vcfs)
	(Survivor_vcfs_supp_field,Survivor_vcfs_supp_field_nonRenamed) = APPEND_SUPP_FIELD(Survivor_vcfs)

	Survivor_collected_sample_paths = Survivor_vcfs_supp_field.map({it -> it[2]}).collect()
	Survivor_multiVcf = SURVIVOR_INTERSECT_SAMPLES(Survivor_collected_sample_paths)//.view{"SURVIVOR_INTERSECT_SAMPLES: $it"}
	Survivor_preFilter = Survivor_vcfs_supp_field_nonRenamed.combine(Survivor_multiVcf)
 (Singletons_Survivor, merged) = SURVIVOR_FILTER_SINGLETONS(Survivor_preFilter)



 BAMs_Survivor = Singletons_Survivor.join(BAMs).map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1]
			def vcf = row[2]
			def bam  = row[4]
			def bai  = row[5]

			[name, sample, vcf, bam, bai]	
	})//.view{"____________BAMs_Survivor_____________: $it"}// join BAMs with annotated, vcf fits bam for given sample
		
		CovFilesSurvivor = CALC_COVERAGE_SURVIVOR(BAMs_Survivor)//.view{"CovFilesSurvivor: $it"}
		PARSE_SURVIVOR(Singletons_Survivor.join(CovFilesSurvivor))//.view{"BAMs_Survivor+CovFilesSurvivor: $it"}

	// emit:
	// Survivor_annot
}