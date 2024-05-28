include { SVIM } from "${params.projectDirectory}/modules/varcallers"
include { TAG_UNIQUE_VARS; ANNOTATE; PARSE_SVIM_VCF; CALC_COVERAGE; FILTER_SVIM } from "${params.projectDirectory}/modules/vcf-processors"

workflow wf_svim {
 take:
 BAMs

	main:
	(SvimRawVcf, _) = SVIM(BAMs)
	Vcf_paths = SvimRawVcf.map({it -> [it[2]]})
	Combined_collected_vcf = SvimRawVcf.combine(Vcf_paths.collect().map({it -> [it]}))
	Combined_filtered = Combined_collected_vcf.map({
		row -> // get one vcf and all other vcfs
		def name = row[0]
		def sample = row[1]
		def vcf = row[2]
		def filtered  = removeSame(vcf, row[3])
		[name,sample, vcf, filtered]	
		})
		(Tagged, uniqueIds) = TAG_UNIQUE_VARS(Combined_filtered)
		Annotated = ANNOTATE(Tagged)
		//SVIM_Annotated_Sup2 = FILTER_SVIM(Annotated) //already min 2 supplementary reads
		BAMs_annot = Annotated.join(BAMs).map({
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
	Editedvcfs = PARSE_SVIM_VCF(BAMs_annot_cov) // NOT ALL FINISHED

	emit:
	SvimRawVcf
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
