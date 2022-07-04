# GPU-accelerated ONT Nanopore K8S basecalling 
## Run basecalling 

Computations are running from bioit proxy (use ssh to connect). Analysis is then deployed to k8s on kuba cluster.

1. Create folders for results, project and temps
 add their location to local_nextflow.config
 ```
 k8s {
...
   launchDir = 'PATH_TO_YOUR_LAUNCHDIR'
   projectDir = 'PATH_TO_LOCAL_CODE_FROM_GITHUB_REPO'
   workDir = 'PATH_TO_YOUR_TEMP'
...
}
```
2. Enter flowcell and kit info with path for fast5 files 
Run the script as the example bellow (run.sh)
```
nextflow kuberun xsvato01/nanopore_k8s/main.nf -r main -pod-image 'cerit.io/nextflow/nextflow:22.05.0' -resume -c local_nextflow.config -process.echo //
        --flowcell FLO-FLG001 --kit SQK-LSK110 --data '/mnt/shared/MedGen/ONTdata/k8s_testing/20220614_1150_MN16014_ais387_bfa74e7a/fast5_pass'
```
## Using the GPU-accelerated container in other pipelines
To use the preinstalled Guppy-GPU container in eg. Nextflow pipelines, use the container from
```
'registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:latest'
```
In Nextflow, you can facilitate this container by specifying it in the .config file or the actual script:
```
process XXX {
      container = 'registry.gitlab.ics.muni.cz:443/450402/nanopore_k8s:latest'
      containerOptions '--gpus all'

    script:
      """
      guppy_basecaller ....
      """

}

```
