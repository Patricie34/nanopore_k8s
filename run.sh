# Run the Nextflow process using kuberun
nextflow kuberun xsvato01/nanopore_k8s/main.nf -r main -head-image 'cerit.io/nextflow/nextflow:22.11.1' \
-resume -with-report -params-file samplesheet.json
