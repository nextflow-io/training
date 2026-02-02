# Partie 1 : Aperçu de la méthode et tests manuels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Il existe plusieurs méthodes valides pour traiter et analyser les données RNAseq en vrac.
Pour cette formation, nous suivons la méthode décrite [ici](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) par les Drs. Simon Andrews et Laura Biggins au [Babraham Institute](https://www.babraham.ac.uk/).

Notre objectif est de développer un workflow qui implémente les étapes de traitement suivantes : exécuter un contrôle qualité initial sur les lectures dans un échantillon RNAseq en vrac, couper les séquences d'adaptateurs des lectures, aligner les lectures sur un génome de référence et produire un rapport de contrôle qualité (QC) complet.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC :** Effectuer le QC sur les données de lecture avant la coupe en utilisant FastQC
- **TRIM_GALORE :** Couper les séquences d'adaptateurs et effectuer le QC après la coupe en utilisant Trim Galore (regroupe Cutadapt et FastQC)
- **HISAT2_ALIGN :** Aligner les lectures sur le génome de référence en utilisant Hisat2
- **MULTIQC :** Générer un rapport QC complet en utilisant MultiQC

Cependant, avant de nous lancer dans l'écriture de code de workflow, nous allons tester les commandes manuellement sur des données de test.
Les outils dont nous avons besoin ne sont pas installés dans l'environnement GitHub Codespaces, nous allons donc les utiliser via des conteneurs (voir [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Note"

     Assurez-vous d'être dans le répertoire `nf4-science/rnaseq`. La dernière partie du chemin affichée lorsque vous tapez `pwd` devrait être `rnaseq`.

---

## 1. QC initial et coupe des adaptateurs

Nous allons récupérer une image de conteneur qui a à la fois `fastqc` et `trim_galore` installés, la lancer en mode interactif et exécuter les commandes de coupe et de QC sur l'un des fichiers de données d'exemple.

### 1.1. Récupérer le conteneur

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Cela vous donne la sortie console suivante pendant que le système télécharge l'image :

??? success "Sortie de la commande"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. Lancer le conteneur en mode interactif

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Sortie de la commande"

    ```console

    ```
-->

Votre invite changera pour quelque chose comme `(base) root@b645838b3314:/tmp#`, ce qui indique que vous êtes maintenant à l'intérieur du conteneur.

La partie `-v ./data:/data` de la commande nous permettra d'accéder au contenu du répertoire `data/` depuis l'intérieur du conteneur.

```bash
ls /data/reads
```

??? success "Sortie de la commande"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. Exécuter la première commande `fastqc`

Exécutons `fastqc` pour collecter les métriques de contrôle qualité sur les données de lecture.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Sortie de la commande"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Cela devrait s'exécuter très rapidement.
Vous pouvez trouver les fichiers de sortie dans le même répertoire que les données originales :

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="Sortie"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Couper les séquences d'adaptateurs avec `trim_galore`

Maintenant, exécutons `trim_galore`, qui regroupe Cutadapt et FastQC, pour couper les séquences d'adaptateurs et collecter les métriques QC après la coupe.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

Le flag `--fastqc` fait en sorte que la commande exécute automatiquement une étape de collecte QC une fois la coupe terminée.

_La sortie est très verbeuse, ce qui suit est donc abrégé._

??? success "Sortie de la commande"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Vous pouvez trouver les fichiers de sortie dans le répertoire de travail :

```bash
ls ENCSR000COQ1_1*
```

```console title="Sortie"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Déplacer les fichiers de sortie vers le système de fichiers en dehors du conteneur

Tout ce qui reste à l'intérieur du conteneur sera inaccessible pour les travaux futurs, déplaçons donc ces fichiers vers un nouveau répertoire.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Quitter le conteneur

```bash
exit
```

---

## 2. Aligner les lectures sur le génome de référence

Nous allons récupérer une image de conteneur qui a `hisat2` installé, la lancer en mode interactif et exécuter la commande d'alignement pour aligner les données RNAseq sur un génome de référence.

### 2.1. Récupérer le conteneur `hisat2`

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Sortie de la commande"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. Lancer le conteneur `hisat2` en mode interactif

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

La commande est la même qu'avant, avec l'URI du conteneur pertinent échangée.

### 2.3. Créer les fichiers d'index du génome Hisat2

Hisat2 nécessite que la référence du génome soit fournie dans un format très spécifique, et ne peut pas simplement consommer le fichier FASTA `genome.fa` que nous fournissons, nous allons donc profiter de cette occasion pour créer les ressources pertinentes.

```bash
hisat2-build /data/genome.fa genome_index
```

La sortie est très verbeuse, ce qui suit est donc abrégé :

<!-- TODO: switch to full output -->

??? success "Sortie de la commande"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

Cela crée plusieurs fichiers d'index du génome, que vous pouvez trouver dans le répertoire de travail.

```bash
ls genome_index.*
```

```console title="Sortie"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Nous les utiliserons dans un instant, mais générons d'abord une archive tar compressée avec ces fichiers d'index du génome ; nous en aurons besoin plus tard et générer ces fichiers n'est généralement pas quelque chose que nous voulons faire dans le cadre d'un workflow.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

Cela stocke une archive tar `genome_index.tar.gz` contenant les fichiers d'index du génome dans le répertoire `data/` de notre système de fichiers, ce qui sera utile dans la Partie 2 de cette formation.

### 2.4. Exécuter la commande `hisat2`

Maintenant nous pouvons exécuter la commande d'alignement, qui effectue l'étape d'alignement avec `hisat2` puis redirige la sortie vers `samtools` pour écrire la sortie sous forme de fichier BAM.

L'entrée de données de lecture est le fichier `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que nous avons généré avec `trim_galore` dans l'étape précédente.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Sortie de la commande"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Cela s'exécute presque instantanément car c'est un fichier de test très petit.
À l'échelle réelle, cela pourrait prendre beaucoup plus de temps.

Une fois de plus, vous pouvez trouver les fichiers de sortie dans le répertoire de travail :

```bash
ls ENCSR000COQ1_1*
```

```console title="Sortie"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Déplacer les fichiers de sortie vers le système de fichiers en dehors du conteneur

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Quitter le conteneur

```bash
exit
```

---

## 3. Générer un rapport QC complet

Nous allons récupérer une image de conteneur qui a `multiqc` installé, la lancer en mode interactif et exécuter une commande de génération de rapport sur les fichiers de rapport FastQC avant/après.

### 3.1. Récupérer le conteneur `multiqc`

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Sortie de la commande"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. Lancer le conteneur `multiqc` en mode interactif

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. Exécuter la commande `multiqc`

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Sortie de la commande"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC est capable de rechercher dans les répertoires des rapports QC compatibles et agrégera tout ce qu'il trouve.

Ici, nous voyons que l'outil a trouvé les trois rapports QC que nous avons générés : le QC initial que nous avons fait avec `fastqc`, le rapport après coupe de `cutadapt` (fait via `trim_galore`) et le QC après alignement produit par `hisat2`.

Les fichiers de sortie sont une fois de plus dans le répertoire de travail :

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Sortie"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Déplacer les fichiers de sortie vers le système de fichiers en dehors du conteneur

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Quitter le conteneur

```bash
exit
```

---

### À retenir

Vous avez testé toutes les commandes individuelles de manière interactive dans les conteneurs pertinents.

### Et ensuite ?

Apprenez à encapsuler ces mêmes commandes dans un workflow multi-étapes qui utilise des conteneurs pour exécuter le travail.
