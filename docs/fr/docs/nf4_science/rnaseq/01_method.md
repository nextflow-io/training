# Partie 1 : Aperçu de la méthode

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Il existe plusieurs méthodes valides pour traiter et analyser les données RNAseq en vrac.
Pour cette formation, nous suivons la méthode décrite [ici](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) par les Drs. Simon Andrews et Laura Biggins au [Babraham Institute](https://www.babraham.ac.uk/).

Notre objectif est de développer un workflow qui implémente les étapes de traitement suivantes : exécuter un contrôle qualité initial sur les lectures dans un échantillon RNAseq en vrac, couper les séquences d'adaptateurs des lectures, aligner les lectures sur un génome de référence et produire un rapport de contrôle qualité (QC) complet.

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC :** Effectuer le QC sur les données de lecture avant la coupe en utilisant FastQC
- **TRIM_GALORE :** Couper les séquences d'adaptateurs et effectuer le QC après la coupe en utilisant Trim Galore (regroupe Cutadapt et FastQC)
- **HISAT2_ALIGN :** Aligner les lectures sur le génome de référence en utilisant Hisat2
- **MULTIQC :** Générer un rapport QC complet en utilisant MultiQC

### Méthodes

Nous allons vous montrer comment appliquer ces étapes de traitement en deux phases.
Nous commencerons par le **traitement d'un échantillon unique** qui exécute les outils de QC, de coupe et d'alignement sur un échantillon.
Ensuite, nous étendrons au **traitement multi-échantillons** qui exécute les mêmes outils sur plusieurs échantillons et génère un rapport de contrôle qualité agrégé.

Avant de nous lancer dans l'écriture de code de workflow pour l'une ou l'autre approche, nous allons tester les commandes manuellement sur des données de test.

### Jeu de données

Nous fournissons les données et ressources associées suivantes :

- **Données RNAseq** (`reads/`) : fichiers FASTQ de six échantillons, réduits à une petite région pour limiter la taille des fichiers. Chaque échantillon a des lectures paired-end (deux fichiers par échantillon), bien que nous commencions par travailler uniquement avec des lectures single-end.
- **Un génome de référence** (`genome.fa`) : une petite région du chromosome 20 humain (de hg19/b37).
- **Feuilles d'échantillons CSV** (`single-end.csv` et `paired-end.csv`) : fichiers listant les identifiants et chemins des fichiers de données d'exemple.

### Logiciels

Les quatre outils principaux impliqués sont [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) pour la collecte de métriques de contrôle qualité, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) pour la coupe des adaptateurs (regroupe Cutadapt et FastQC pour le QC après coupe), [HISAT2](http://daehwankimlab.github.io/hisat2/) pour l'alignement épissé sur un génome de référence, et [MultiQC](https://multiqc.info/) pour la génération de rapports QC agrégés.

Ces outils ne sont pas installés dans l'environnement GitHub Codespaces, nous allons donc les utiliser via des conteneurs récupérés via le service Seqera Containers (voir [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Astuce"

     Assurez-vous d'être dans le répertoire `nf4-science/rnaseq`. La dernière partie du chemin affichée lorsque vous tapez `pwd` devrait être `rnaseq`.

---

## 1. Traitement d'un échantillon unique

Dans cette section, nous testons les commandes qui traitent un seul échantillon RNAseq : contrôle qualité, coupe des adaptateurs et alignement sur un génome de référence.
Ce sont les commandes que nous encapsulerons dans un workflow Nextflow dans la Partie 2 de cette formation.

1. Exécuter le QC initial sur un fichier FASTQ en utilisant FastQC
2. Couper les séquences d'adaptateurs et exécuter le QC après coupe en utilisant Trim Galore
3. Aligner les lectures coupées sur le génome de référence en utilisant HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Nous commençons par tester ces commandes sur un seul échantillon.

### 1.1. QC et coupe des adaptateurs

Tout d'abord, nous voulons exécuter les commandes de QC et de coupe sur l'un des fichiers de données d'exemple.

#### 1.1.1. Récupérer le conteneur

Récupérons une image de conteneur qui a à la fois `fastqc` et `trim_galore` installés :

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

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

Si vous n'avez pas téléchargé cette image auparavant, cela peut prendre une minute.
Une fois terminé, vous avez une copie locale de l'image du conteneur.

#### 1.1.2. Lancer le conteneur en mode interactif

Pour exécuter le conteneur en mode interactif, utilisez `docker run` avec les flags `-it`.
L'option `-v ./data:/data` monte notre répertoire local `data/` afin que nous puissions accéder aux fichiers d'entrée depuis l'intérieur du conteneur.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Sortie de la commande"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Votre invite changera pour quelque chose comme `(base) root@b645838b3314:/tmp#`, ce qui indique que vous êtes maintenant à l'intérieur du conteneur.

Vérifiez que vous pouvez voir les fichiers de données de séquençage sous `/data/reads` :

```bash
ls /data/reads
```

??? abstract "Contenu du répertoire"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Avec cela, vous êtes prêt·e à essayer votre première commande.

#### 1.1.3. Exécuter la commande FastQC

La méthode référencée ci-dessus nous donne la ligne de commande pour exécuter le QC sur un seul fichier.
Nous devons seulement fournir le fichier d'entrée ; l'outil générera automatiquement les fichiers de sortie dans le même répertoire que les données originales.

Exécutez la commande `fastqc` sur un fichier de données :

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

??? abstract "Contenu du répertoire"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Vous devriez voir un rapport HTML et une archive ZIP contenant les métriques QC.
Cela complète le test de la première étape.

#### 1.1.4. Couper les séquences d'adaptateurs avec Trim Galore

Maintenant, exécutons `trim_galore`, qui regroupe Cutadapt et FastQC, pour couper les séquences d'adaptateurs et collecter les métriques QC après la coupe.
Comme indiqué ci-dessus, le logiciel est inclus dans le même conteneur, donc aucun changement nécessaire à ce niveau.

La commande est simple ; nous devons simplement ajouter le flag `--fastqc` pour exécuter automatiquement une étape de collecte QC une fois la coupe terminée.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Sortie de la commande"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

La sortie est très verbeuse, nous avons donc mis en évidence les lignes les plus pertinentes dans l'exemple ci-dessus.
Vous pouvez trouver les fichiers de sortie dans le répertoire de travail :

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contenu du répertoire"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Cela inclut les lectures coupées, le rapport de coupe et les fichiers QC après coupe.

#### 1.1.5. Déplacer les fichiers de sortie

Tout ce qui reste à l'intérieur du conteneur sera inaccessible pour les travaux futurs, nous devons donc déplacer ces fichiers vers un répertoire sur le système de fichiers monté.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Contenu du répertoire"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Les fichiers sont maintenant accessibles dans votre système de fichiers normal.

#### 1.1.6. Quitter le conteneur

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite devrait revenir à la normale ; cela complète le test des deux premières étapes.

### 1.2. Aligner les lectures sur le génome de référence

Ensuite, nous voulons exécuter la commande d'alignement pour aligner les lectures RNAseq coupées sur un génome de référence.

#### 1.2.1. Récupérer le conteneur

Récupérons une image de conteneur qui a `hisat2` et `samtools` installés :

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

Vous remarquerez que certaines couches affichent `Already exists` car elles sont partagées avec l'image du conteneur Trim Galore que nous avons récupérée précédemment.
En conséquence, ce téléchargement devrait être plus rapide que le premier.

#### 1.2.2. Lancer le conteneur en mode interactif

Lancez le conteneur en mode interactif, en utilisant la même approche qu'avant avec l'URI du conteneur pertinent échangée.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Votre invite changera à nouveau pour indiquer que vous êtes à l'intérieur du conteneur.

#### 1.2.3. Créer les fichiers d'index du génome

HISAT2 nécessite que la référence du génome soit fournie dans un format très spécifique, et ne peut pas simplement consommer le fichier FASTA `genome.fa` que nous fournissons, nous allons donc profiter de cette occasion pour créer les ressources pertinentes.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Sortie de la commande"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

La sortie est très verbeuse, nous avons donc mis en évidence quelques lignes pertinentes dans l'exemple ci-dessus.

Cela crée plusieurs fichiers d'index du génome, que vous pouvez trouver dans le répertoire de travail.

```bash
ls genome_index.*
```

??? abstract "Contenu du répertoire"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Nous aurons besoin de ces fichiers plus tard, et générer ces fichiers n'est généralement pas quelque chose que nous voulons faire dans le cadre d'un workflow, nous allons donc générer une archive tar compressée contenant les fichiers d'index du génome que nous pouvons facilement transmettre selon les besoins.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Sortie de la commande"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

Nous déplacerons l'archive tar `genome_index.tar.gz` résultante contenant les fichiers d'index du génome vers le répertoire `data/` de notre système de fichiers dans quelques minutes.
Cela sera utile dans la Partie 2 de cette formation.

#### 1.2.4. Exécuter la commande d'alignement

Maintenant nous pouvons exécuter la commande d'alignement, qui effectue l'étape d'alignement avec `hisat2` puis redirige la sortie vers `samtools` pour écrire la sortie sous forme de fichier BAM.

L'entrée de données de lecture est le fichier `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que nous avons généré avec `trim_galore` dans l'étape précédente.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
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

??? abstract "Contenu du répertoire"

    ```console title="Sortie"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

L'alignement a produit un fichier BAM et un fichier journal avec les statistiques d'alignement.

#### 1.2.5. Déplacer les fichiers de sortie

Comme précédemment, déplacez les fichiers de sortie vers un répertoire sur le système de fichiers monté afin qu'ils restent accessibles après avoir quitté le conteneur.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Avec cela fait, nous avons tout ce dont nous avons besoin.

#### 1.2.6. Quitter le conteneur

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite devrait revenir à la normale.
Cela conclut l'exécution de test du traitement d'un échantillon unique.

!!! example "Écrivez-le sous forme de workflow !"

    N'hésitez pas à passer directement à la [Partie 2](./02_single-sample.md) si vous souhaitez commencer à implémenter cette analyse sous forme de workflow Nextflow.
    Vous devrez simplement revenir pour terminer le deuxième cycle de tests avant de passer à la Partie 3.

---

## 2. Agrégation QC multi-échantillons

Les commandes que nous venons de tester traitent un échantillon à la fois.
En pratique, nous devons généralement traiter de nombreux échantillons puis agréger les résultats QC sur tous pour évaluer la qualité de l'ensemble de données global.

[MultiQC](https://multiqc.info/) est un outil qui recherche dans les répertoires des rapports QC de nombreux outils bioinformatiques courants et les agrège dans un seul rapport HTML complet.
Il peut reconnaître les sorties de FastQC, Cutadapt (via Trim Galore) et HISAT2, parmi beaucoup d'autres.

Ici, nous traitons deux échantillons supplémentaires avec les mêmes outils par échantillon, puis utilisons MultiQC pour agréger les rapports QC sur les trois échantillons.
Ce sont les commandes que nous encapsulerons dans un workflow Nextflow dans la Partie 3 de cette formation.

1. Exécuter le QC et la coupe sur des échantillons supplémentaires en utilisant Trim Galore
2. Exécuter l'alignement sur des échantillons supplémentaires en utilisant HISAT2
3. Agréger tous les rapports QC dans un rapport complet en utilisant MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC et coupe d'échantillons supplémentaires

Les commandes de QC et de coupe par échantillon sont identiques à ce que nous avons exécuté dans la section 1.1.
Nous avons déjà récupéré l'image du conteneur, nous pouvons donc la lancer directement.

#### 2.1.1. Lancer le conteneur

Nous avons déjà récupéré cette image de conteneur dans la section 1.1, nous pouvons donc la lancer directement :

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Votre invite change pour indiquer que vous êtes à l'intérieur du conteneur.

#### 2.1.2. Exécuter le QC et la coupe sur des échantillons supplémentaires

Exécutez FastQC et Trim Galore sur deux échantillons supplémentaires, l'un après l'autre.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Une fois terminé, vous devriez avoir les fichiers de sortie de Trim Galore pour les deux échantillons dans le répertoire de travail.

#### 2.1.3. Déplacer les fichiers de sortie

Déplacez les fichiers de sortie de Trim Galore vers le même répertoire que nous avons utilisé dans la section 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Contenu du répertoire"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

Les fichiers sont maintenant accessibles dans votre système de fichiers normal.

#### 2.1.4. Quitter le conteneur

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite devrait revenir à la normale.

### 2.2. Aligner des échantillons supplémentaires

Les commandes d'alignement sont identiques à ce que nous avons exécuté dans la section 1.2.
Nous devons extraire l'index du génome de l'archive tar que nous avons sauvegardée précédemment, car les fichiers d'index originaux ont été créés à l'intérieur d'un conteneur qui n'existe plus.

#### 2.2.1. Lancer le conteneur

Nous avons déjà récupéré cette image de conteneur dans la section 1.2, nous pouvons donc la lancer directement :

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Votre invite change pour indiquer que vous êtes à l'intérieur du conteneur.

#### 2.2.2. Extraire l'index du génome

Extrayez les fichiers d'index du génome de l'archive tar que nous avons sauvegardée sur le système de fichiers monté :

```bash
tar -xzf /data/genome_index.tar.gz
```

Cela restaure les fichiers `genome_index.*` dans le répertoire de travail.

#### 2.2.3. Exécuter l'alignement sur des échantillons supplémentaires

Exécutez l'alignement HISAT2 sur les deux échantillons nouvellement coupés, l'un après l'autre.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Une fois terminé, vous devriez avoir des fichiers BAM et des fichiers journaux pour les deux échantillons dans le répertoire de travail.

#### 2.2.4. Déplacer les fichiers de sortie

Déplacez les fichiers de sortie d'alignement vers le même répertoire que nous avons utilisé dans la section 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Contenu du répertoire"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Les fichiers sont maintenant accessibles dans votre système de fichiers normal.

#### 2.2.5. Quitter le conteneur

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite devrait revenir à la normale.

### 2.3. Générer un rapport QC complet

Maintenant que nous avons les sorties de QC, de coupe et d'alignement pour trois échantillons, nous pouvons utiliser MultiQC pour les agréger dans un seul rapport.
MultiQC recherche dans les répertoires des rapports QC compatibles et agrège tout ce qu'il trouve.

#### 2.3.1. Récupérer le conteneur

Récupérons une image de conteneur qui a `multiqc` installé :

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Sortie de la commande"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
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
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

Vous remarquerez que certaines couches affichent `Already exists` car elles sont partagées avec les images de conteneur que nous avons récupérées précédemment.
En conséquence, ce téléchargement devrait être plus rapide que les précédents.

#### 2.3.2. Lancer le conteneur en mode interactif

Lancez le conteneur en mode interactif avec le répertoire de données monté, comme précédemment.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Votre invite changera pour indiquer que vous êtes à l'intérieur du conteneur.

#### 2.3.3. Exécuter la commande MultiQC

Exécutez `multiqc`, en le pointant vers les répertoires où nous avons stocké les fichiers de sortie liés au QC pour les trois échantillons.
Le flag `-n` définit le nom du rapport de sortie.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Sortie de la commande"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

Ici, nous voyons que l'outil a trouvé les rapports QC pour les trois échantillons : le QC initial que nous avons fait avec `fastqc`, les rapports après coupe de `cutadapt` (via `trim_galore`) et les résumés d'alignement produits par `hisat2`.

Les fichiers de sortie sont dans le répertoire de travail :

```bash
ls all_samples_QC*
```

??? abstract "Contenu du répertoire"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

La sortie principale est le rapport `all_samples_QC.html`, accompagné d'un répertoire de données contenant les métriques sous-jacentes.

#### 2.3.4. Déplacer les fichiers de sortie

Déplacez le rapport et son répertoire de données vers le système de fichiers monté.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Les fichiers sont maintenant accessibles dans votre système de fichiers normal.

#### 2.3.5. Quitter le conteneur

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite devrait revenir à la normale.
Cela conclut le test de toutes les commandes de traitement RNAseq.

---

### À retenir

Vous savez comment exécuter les commandes FastQC, Trim Galore, HISAT2 et MultiQC dans leurs conteneurs respectifs, y compris comment traiter plusieurs échantillons et agréger les rapports QC.

### Et ensuite ?

Faites une pause, puis passez à la [Partie 2](./02_single-sample.md) pour apprendre à encapsuler ces mêmes commandes dans des workflows qui utilisent des conteneurs pour exécuter le travail.
