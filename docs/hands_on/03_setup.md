# Environment Setup

## Pipeline data

All the files needed for the hands-on activity are stored in the directory shown below:

```bash
tree /workspace/gitpod/hands-on
```

```bash
/workspace/gitpod/hands-on
hands-on
├── README.md
├── bin
│   └── gghist.R
├── data
│   ├── blacklist.bed
│   ├── genome.fa
│   ├── known_variants.vcf.gz
│   └── reads
│       ├── ENCSR000COQ1_1.fastq.gz
│       ├── ENCSR000COQ1_2.fastq.gz
│       ├── ENCSR000COQ2_1.fastq.gz
│       ├── ENCSR000COQ2_2.fastq.gz
│       ├── ENCSR000COR1_1.fastq.gz
│       ├── ENCSR000COR1_2.fastq.gz
│       ├── ENCSR000COR2_1.fastq.gz
│       ├── ENCSR000COR2_2.fastq.gz
│       ├── ENCSR000CPO1_1.fastq.gz
│       ├── ENCSR000CPO1_2.fastq.gz
│       ├── ENCSR000CPO2_1.fastq.gz
│       └── ENCSR000CPO2_2.fastq.gz
├── final_main.nf
└── nextflow.config
```

4 directories, 19 files

## Pulling the Docker image

Nextflow can pull Docker images at runtime, but let’s just download it manually to see how Docker works:

```bash
docker pull cbcrg/callings-with-gatk:latest
```

You should see the progress of the download:

```console
sha256:93910bf77bc197cb790eca776e42950bc8eff117bdc6e67157295e09c98fc381: Pulling from cbcrg/callings-with-gatk
915665fee719: Downloading [=============================================>     ] 47.08 MB/51.36 MB
f332de2321e6: Downloading [===========>                                       ] 41.96 MB/187.8 MB
1577a6dd9e43: Downloading [===============================>                   ] 46.72 MB/73.45 MB
7059d9bb5245: Waiting
71863f70269f: Waiting
ce2a2879246d: Waiting
e38ba5d5f9fb: Waiting
90158da87bb2: Waiting
```

and the following message when the pull is completed:

```console
Digest: sha256:93910bf77bc197cb790eca776e42950bc8eff117bdc6e67157295e09c98fc381
Status: Downloaded newer image for cbcrg/callings-with-gatk:latest
```

!!! note

    Ideally, you should have a light container image for every process, containing exclusively what is required for that task. This not only contributes to the image being light to pull, run and stop but also makes it easier to debug in case of problems. For this course specifically, we won't bother with that and use a single bulky container image (mentioned above) that has everything we need.

## Script permission

Make sure the following R script has execute permissions:

```bash
chmod +x /workspace/gitpod/hands-on/bin/gghist.R
```
