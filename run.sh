nextflow kuberun xsvato01/nanopore_k8s/main.nf -r main -pod-image 'cerit.io/nextflow/nextflow:22.05.0' \
-resume -c zaloha_nextflow.config -process.echo --flowcell FLO-MIN106 --kit SQK-LSK110 \
 --data /mnt/shared/MedGen/ONTdata/Sabina_HMW/1727_x100_lig_r941/20220615_1152_MN16014_FAR41390_26e0e79d

