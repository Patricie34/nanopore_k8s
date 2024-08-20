include { QUALIMAP } from "${params.projectDirectory}/modules/qc"
include { CNVkit } from "${params.projectDirectory}/modules/cnv"

include { wf_process_till_bam } from "${params.projectDirectory}/subworkflows/wf_process_till_bam"
include { wf_svim } from "${params.projectDirectory}/subworkflows/wf_svim"
include { wf_delly } from "${params.projectDirectory}/subworkflows/wf_delly"
include { wf_survivor } from "${params.projectDirectory}/subworkflows/wf_survivor"


workflow {
	///  GATHERING INPUTS
	Runlist = channel.fromList(params.samplesheet)
	References = channel.fromList(params.ref_specific)


	BAMs = wf_process_till_bam(Runlist, References)
	// Svim_annot_vcfs = wf_svim(BAMs)
	// Delly_vcfs = 
	wf_delly(BAMs)
	// wf_survivor(BAMs, Svim_annot_vcfs, Delly_vcfs)
	// CNVkit(BAMs)
	QUALIMAP(BAMs)

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