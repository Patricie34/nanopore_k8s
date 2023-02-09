nextflow kuberun xsvato01/nanopore_k8s/main.nf -r main -pod-image 'cerit.io/nextflow/nextflow:22.06.1' \
-c zaloha_nextflow.config -resume -with-report -params-file ./samplesheet.yaml --flowcell FLO-MIN106 --kit SQK-LSK110
