include { SNIFFLES2; PEPPER_DEEPVARIANT } from "${params.projectDirectory}/modules/varcallers"
include { ANNOTATE_SURVIROR; CALC_COVERAGE_SURVIVOR;  PEPPER_PREFILT; SURVIVOR; PEPPER_UMIREADS} from "${params.projectDirectory}/modules/vcf-processors"

workflow wf_survivor {
 take:
 BAMs
	Svim_vcfs

	main:
 Sniffles_vcfs = SNIFFLES2(BAMs)
	Primary_bam = PEPPER_PREFILT(BAMs)
 Bam_unique_reads = PEPPER_UMIREADS(BAMs)
	joinedVCFs = Sniffles_vcfs.join(Svim_vcfs).map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1]
			def vcf1 = row[2]
			def vcf2  = row[4]

			[name, sample, vcf1, vcf2]	
	})

	Survivor_vcfs = SURVIVOR(joinedVCFs)
	Survivor_annot = ANNOTATE_SURVIROR(Survivor_vcfs)

 BAMs_Survivor = Survivor_vcfs.join(BAMs).map({
			// get rid of double sample info from 
			row ->
			def name = row[0]
			def sample = row[1]
			def vcf = row[2]
			def bam  = row[4]
			def bai  = row[5]

			[name, sample, vcf, bam, bai]	
	})//.view{"____________BAMs_Survivor_____________: $it"}// join BAMs with annotated, vcf fits bam for given sample
		
		// CovFilesSurvivor = CALC_COVERAGE_SURVIVOR(BAMs_Survivor)

	emit:
	Survivor_annot
}