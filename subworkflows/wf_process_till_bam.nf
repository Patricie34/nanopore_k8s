
include { REFORMAT_SAMPLE; REFORMAT_PARAMS } from "${params.projectDirectory}/modules/input-reformators"
include { GUPPY; DORADO } from "${params.projectDirectory}/modules/basecallers"
include { COLLECT_FQs; COLLECT_BAMs } from "${params.projectDirectory}/modules/collectors"
include { MINIMAP2; BAM_INDEX_SORT } from "${params.projectDirectory}/modules/mappers"


workflow wf_process_till_bam {
 take:
	Runlist
	References

	main:
	SamplesReformated = REFORMAT_SAMPLE(Runlist)
	ReferencesReformated = REFORMAT_PARAMS(References)
	SamplesWithReferences = SamplesReformated.combine(ReferencesReformated,by:1) //match parameters to individual samples
		.map({ sample -> [sample[1],sample[2]+(sample[3]) ]})
/// CALLING, MAPPING
	//FQcalled = GUPPY(Runlist)
	BAMdorado = DORADO(Runlist)
	FQcollected = COLLECT_FQs(Runlist)
	BAMminimap2	= MINIMAP2(FQcollected.join(SamplesWithReferences))
	BAMSorted = BAM_INDEX_SORT(BAMminimap2.mix(BAMdorado))
	BAMcollectedSorted = COLLECT_BAMs(SamplesWithReferences)
	BAMs = BAMcollected.mix(BAMSorted)

	emit:
	BAMs
}