include { REFORMAT_SAMPLE; REFORMAT_PARAMS } from "${params.projectDirectory}/modules/input-reformators"
include { FAST5TOPOD5; DORADO } from "${params.projectDirectory}/modules/basecallers"
include { COLLECT_FQs; COLLECT_BAMs } from "${params.projectDirectory}/modules/collectors"
include { MINIMAP2; BAM_INDEX_SORT; BAM2FASTQ } from "${params.projectDirectory}/modules/mappers"


workflow wf_process_till_bam {
	take:
	Runlist
	References

	main:
	SamplesReformated = Runlist.map {[it.name, it.reference, it]}
	ReferencesReformated = References.map {[it, it.refname]}

	SamplesWithReferences = SamplesReformated.combine(ReferencesReformated,by:1) //match parameters to individual samples
		.map({ sample -> [sample[1],sample[2]+(sample[3]) ]})

// SamplesWithReferences.view()

	BAMcollectedSorted = SamplesWithReferences.map { row ->
			if (row[1].type == 'bam') {

        return [row[0], row[1],
            file("${row[1].path}/${row[0]}.remapped.sorted.bam", checkIfExists: true),
            file("${row[1].path}/${row[0]}.remapped.sorted.bam.bai", checkIfExists: true)
		]
    }
	}

	// BAMcollectedSorted = COLLECT_BAMs(SamplesWithReferences)
	// extractedFQs = BAM2FASTQ(BAMcollectedSorted)

	Pod5s = FAST5TOPOD5(SamplesWithReferences)//.view{"____________Pod5s_____________: $it"}
	BAMdorado = DORADO(SamplesWithReferences)
	FQcollected = COLLECT_FQs(Runlist)
	BAMminimap2	= MINIMAP2(FQcollected.join(SamplesWithReferences))
		// BAMminimap2	= MINIMAP2(extractedFQs.mix(FQcollected.join(SamplesWithReferences)))

	BAMSorted = BAM_INDEX_SORT(BAMminimap2.mix(BAMdorado))//.view{"____________Sorted_____________: $it"}

	BAMs = BAMcollectedSorted.mix(BAMSorted)//.view{"____________BAM_____________: $it"} //COMMENTED OUT TO MAKE DORADO INDEPENDENT


	emit:
	BAMs
}