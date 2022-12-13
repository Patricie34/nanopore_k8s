nextflow kuberun xsvato01/nanopore_k8s/main.nf -r main -pod-image 'cerit.io/nextflow/nextflow:22.06.1' \
-c zaloha_nextflow.config -resume -process.echo --flowcell FLO-MIN106 --kit SQK-LSK110 \
 --data /mnt/shared/MedGen/ONTdata/Sabina_HMW/BRNO2013/

