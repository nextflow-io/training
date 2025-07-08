# Nextflow Command Reference

This is a tutorial  [link][ref].

 [ref]: https://training.nextflow.io/latest/  "Tutorial"

This is a container reg (seqera)  [link][ref].

 [ref]: https://seqera.io/containers/ "Container"

Searched container example (MAC arm) [link][ref].

[ref]: https://seqera.io/containers/?packages=conda-forge::cowpy=1.1.5 "Container Search cowpy"

## Check the Nextflow Version
```sh
nextflow -v
```

# Run a Nextflow Pipeline 
```sh
nextflow run <File.nf>
```
# Resume a Previous Run
```sh
nextflow run hello-world.nf -resume
```
# View Nextflow Logs
```sh
nextflow log 
```

# Dry Run (Preview):
```sh
nextflow clean -before <runnamefromlog> -n
```
-n is for a dry run (shows what would be deleted).

Force Clean (Delete):
```sh
nextflow clean -before <runnamefromlog> -f
```
-f is for force (actually deletes files).

# Passing a variabe 
```sh
nextflow run hello-world.nf --greeting "Hek varlden"
```
---
# More options 
Trace Report
Creates a tabular text file with detailed info about every process:
```sh 
nextflow run pipeline.nf -with-trace
```

Execution Report (HTML)
Generates a report.html file with pipeline stats.
```sh
nextflow run pipeline.nf -with-report
```
Timeline Visualization (HTML)
Generates a timeline.html showing process flow and duration.
```sh
nextflow run pipeline.nf -with-timeline
```
Execution DAG
Creates a flow chart of your pipeline.
```sh
nextflow run pipeline.nf -with-dag flowchart.png
```
Input variables 
Calling an input varibale with the value, we can also use the default value
```sh
nextflow run hello-world.nf --greeting "He Ali" -with-report  
nextflow run hello-world.nf
```
Output 
```
nextflow run hello-channels.nf -ansi-log false 
```
Operators
```
 #useed for array, ref ex 03
 .flatten()

 #Used to read file csv, ref ex 03
 .splitCsv()

 # view before and after, ref ex 03
 .view()

 # Concatinate all the output to single file, ref ex 04
 .collect() 

 # Count the size of the input channel and assign it to the variable 
.size()
```
Workflow 
```
Pipeline run in parallel, till now we only check the single process workflow
1. How we will use channels to connect those steps 
2. Multiple tasks, which can collapse into a single process 
3. We also look the process with multiple inputs and multiple outputs 
```
MultiInput, ref ex 05
```
# overriding the default value from input variable 
nextflow run hello-workflow.nf --batch_name "demo-output"

```
Multiple Output, ref ex 06 
```
# Using the emit to access them in positional manner, like a decorator
```
Nextflow Module Ex 07 
```
Example NF Core Module repository 
include { sayHello } from './modules/sayHello.nf'

```
Container 
```
Community run container 

Run the singularity 

# Pull the singularity image 
apptainer pull cowpy.sif https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bde91ba24078478208035944326ee69a6063ae0f29523e017840e863c42e0b93/data

Run the docker

# Pull the image from seqera 
docker pull community.wave.seqera.io/library/cowpy:1.1.5--1a5414f41afdfd25

# run the image
docker run --rm community.wave.seqera.io/library/cowpy:1.1.5--1a5414f41afdfd25  cowpy

docker run -it community.wave.seqera.io/library/cowpy:1.1.5--1a5414f41afdfd25  /bin/bash

# mount a directory 
docker run --rm -it -v .:/data community.wave.seqera.io/library/cowpy:1.1.5--1a5414f41afdfd25  /bin/bash

# run this CLI inside the container 
cowpy "Hello Containers" -c tux
cowpy "Hello Containers" -c cheese 
cowpy "Hello Containers" -c dragonandcow
Run the podman
```