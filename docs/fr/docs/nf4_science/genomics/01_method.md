# Partie 1 : Présentation de la méthode et tests manuels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'appel de variants est une méthode d'analyse génomique qui vise à identifier les variations dans une séquence génomique par rapport à un génome de référence.
Nous allons utiliser ici des outils et des méthodes conçus pour appeler des variants germinaux courts, _c'est-à-dire_ des SNP et des indels, dans des données de séquençage de génome entier.

![Pipeline GATK](img/gatk-pipeline.png)

Un pipeline complet d'appel de variants implique généralement de nombreuses étapes, notamment l'alignement sur la référence (parfois appelé alignement du génome) ainsi que le filtrage et la priorisation des variants.
Pour des raisons de simplicité, dans ce cours nous allons nous concentrer uniquement sur la partie appel de variants.

### Méthodes

Nous allons vous montrer deux façons d'appliquer l'appel de variants à des échantillons de séquençage de génome entier pour identifier les SNP et les indels germinaux.
Nous commencerons d'abord par une **approche par échantillon** simple qui appelle les variants indépendamment pour chaque échantillon.
Ensuite, nous vous montrerons une **approche d'appel conjoint** plus sophistiquée qui analyse plusieurs échantillons ensemble, produisant des résultats plus précis et informatifs.

Avant de nous lancer dans l'écriture de code de workflow pour l'une ou l'autre approche, nous allons tester les commandes manuellement sur des données de test.

### Jeu de données

Nous fournissons les données et ressources associées suivantes :

- **Un génome de référence** constitué d'une petite région du chromosome 20 humain (issu de hg19/b37) et ses fichiers accessoires (index et dictionnaire de séquence).
- **Trois échantillons de séquençage de génome entier** correspondant à un trio familial (mère, père et fils), qui ont été réduits à une petite tranche de données sur le chromosome 20 pour maintenir de petites tailles de fichier.
  Il s'agit de données de séquençage Illumina en lectures courtes qui ont déjà été alignées sur le génome de référence, fournies au format [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, une version compressée de SAM, Sequence Alignment Map).
- **Une liste d'intervalles génomiques**, c'est-à-dire des coordonnées sur le génome où nos échantillons ont des données appropriées pour appeler des variants, fournie au format BED.

### Logiciels

Les deux outils principaux impliqués sont [Samtools](https://www.htslib.org/), une boîte à outils largement utilisée pour manipuler les fichiers d'alignement de séquences, et [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un ensemble d'outils pour la découverte de variants développé au Broad Institute.

Ces outils ne sont pas installés dans l'environnement GitHub Codespaces, nous allons donc les utiliser via des conteneurs (voir [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

    Assurez-vous d'être dans le répertoire `nf4-science/genomics` afin que la dernière partie du chemin affichée lorsque vous tapez `pwd` soit `genomics`.

---

## 1. Appel de variants par échantillon

L'appel de variants par échantillon traite chaque échantillon indépendamment : le programme d'appel de variants examine les données de séquençage pour un échantillon à la fois et identifie les positions où l'échantillon diffère de la référence.

Dans cette section, nous testons les deux commandes qui constituent l'approche d'appel de variants par échantillon : l'indexation d'un fichier BAM avec Samtools et l'appel de variants avec GATK HaplotypeCaller.
Ce sont les commandes que nous encapsulerons dans un workflow Nextflow dans la Partie 2 de ce cours.

1. Générer un fichier d'index pour un fichier d'entrée BAM en utilisant [Samtools](https://www.htslib.org/)
2. Exécuter GATK HaplotypeCaller sur le fichier BAM indexé pour générer des appels de variants par échantillon au format VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Nous commençons par tester les deux commandes sur un seul échantillon.

### 1.1. Indexer un fichier d'entrée BAM avec Samtools

Les fichiers d'index sont une caractéristique courante des formats de fichiers en bioinformatique ; ils contiennent des informations sur la structure du fichier principal qui permettent à des outils comme GATK d'accéder à un sous-ensemble des données sans avoir à lire l'intégralité du fichier.
Ceci est important en raison de la taille que peuvent atteindre ces fichiers.

Les fichiers BAM sont souvent fournis sans index, donc la première étape dans de nombreux workflows d'analyse consiste à en générer un en utilisant `samtools index`.

Nous allons télécharger un conteneur Samtools, le lancer en mode interactif et exécuter la commande `samtools index` sur l'un des fichiers BAM.

#### 1.1.1. Télécharger le conteneur Samtools

Exécutez la commande `docker pull` pour télécharger l'image du conteneur Samtools :

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Sortie de la commande"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

Si vous n'avez pas téléchargé cette image auparavant, cela peut prendre une minute. Une fois terminé, vous disposez d'une copie locale de l'image du conteneur.

#### 1.1.2. Lancer le conteneur Samtools en mode interactif

Pour exécuter le conteneur en mode interactif, utilisez `docker run` avec les options `-it`.
L'option `-v ./data:/data` monte le répertoire local `data` dans le conteneur afin que les outils puissent accéder aux fichiers d'entrée.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Votre invite de commande change pour quelque chose comme `(base) root@a1b2c3d4e5f6:/tmp#`, indiquant que vous êtes maintenant à l'intérieur du conteneur.
Les fichiers de données sont accessibles sous `/data`.

#### 1.1.3. Exécuter la commande d'indexation

La [documentation de Samtools](https://www.htslib.org/doc/samtools-index.html) nous donne la ligne de commande à exécuter pour indexer un fichier BAM.

Nous devons uniquement fournir le fichier d'entrée ; l'outil générera automatiquement un nom pour la sortie en ajoutant `.bai` au nom du fichier d'entrée.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Contenu du répertoire"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Vous devriez maintenant voir un fichier nommé `reads_mother.bam.bai` dans le même répertoire que le fichier BAM d'entrée d'origine.

#### 1.1.4. Quitter le conteneur Samtools

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite de commande devrait maintenant être revenue à ce qu'elle était avant de démarrer le conteneur.

### 1.2. Appeler des variants avec GATK HaplotypeCaller

Nous allons télécharger un conteneur GATK, le lancer en mode interactif et exécuter la commande `gatk HaplotypeCaller` sur le fichier BAM que nous venons d'indexer.

#### 1.2.1. Télécharger le conteneur GATK

Exécutez la commande `docker pull` pour télécharger l'image du conteneur GATK :

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Sortie de la commande"

    Certaines couches affichent `Already exists` car elles sont partagées avec l'image du conteneur Samtools que nous avons téléchargée précédemment.

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

Cela devrait être plus rapide que le premier téléchargement car les deux images de conteneur partagent la plupart de leurs couches.

#### 1.2.2. Lancer le conteneur GATK en mode interactif

Lancez le conteneur GATK en mode interactif avec le répertoire de données monté, exactement comme nous l'avons fait pour Samtools.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Votre invite de commande change pour indiquer que vous êtes maintenant à l'intérieur du conteneur GATK.

#### 1.2.3. Exécuter la commande d'appel de variants

La [documentation de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) nous donne la ligne de commande à exécuter pour effectuer l'appel de variants sur un fichier BAM.

Nous devons fournir le fichier d'entrée BAM (`-I`) ainsi que le génome de référence (`-R`), un nom pour le fichier de sortie (`-O`) et une liste d'intervalles génomiques à analyser (`-L`).

Cependant, nous n'avons pas besoin de spécifier le chemin vers le fichier d'index ; l'outil le recherchera automatiquement dans le même répertoire, en se basant sur la convention de nommage et de co-localisation établie.
Il en va de même pour les fichiers accessoires du génome de référence (fichiers d'index et de dictionnaire de séquence, `*.fai` et `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Sortie de la commande"

    L'outil produit une sortie de journalisation verbeuse. Les lignes surlignées confirment l'achèvement réussi.

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

Le fichier de sortie `reads_mother.vcf` est créé dans votre répertoire de travail à l'intérieur du conteneur, vous ne le verrez donc pas dans l'explorateur de fichiers VS Code à moins de modifier le chemin du fichier de sortie.
Cependant, c'est un petit fichier de test, vous pouvez donc utiliser `cat` pour l'ouvrir et voir son contenu.
Si vous remontez jusqu'au début du fichier, vous trouverez un en-tête composé de nombreuses lignes de métadonnées, suivi d'une liste d'appels de variants, un par ligne.

??? abstract "Contenu du fichier"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Chaque ligne décrit un variant possible identifié dans les données de séquençage de l'échantillon. Pour des conseils sur l'interprétation du format VCF, consultez [cet article utile](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Le fichier de sortie VCF est accompagné d'un fichier d'index nommé `reads_mother.vcf.idx` qui a été automatiquement créé par GATK.
Il a la même fonction que le fichier d'index BAM, permettre aux outils de rechercher et de récupérer des sous-ensembles de données sans charger l'intégralité du fichier.

#### 1.2.4. Quitter le conteneur GATK

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite de commande devrait être revenue à la normale.
Cela conclut le test d'appel de variants par échantillon.

---

## 2. Appel conjoint sur une cohorte

L'approche d'appel de variants que nous venons d'utiliser génère des appels de variants par échantillon.
C'est bien pour examiner les variants de chaque échantillon isolément, mais cela fournit des informations limitées.
Il est souvent plus intéressant d'examiner comment les appels de variants diffèrent entre plusieurs échantillons.
GATK propose une méthode alternative appelée appel de variants conjoint à cette fin.

L'appel de variants conjoint implique de générer un type spécial de sortie de variant appelé GVCF (pour Genomic VCF) pour chaque échantillon, puis de combiner les données GVCF de tous les échantillons et d'exécuter une analyse statistique de 'génotypage conjoint'.

![Analyse conjointe](img/joint-calling.png)

Ce qui est spécial dans le GVCF d'un échantillon, c'est qu'il contient des enregistrements résumant les statistiques des données de séquençage pour toutes les positions dans la zone ciblée du génome, et pas seulement les positions où le programme a trouvé des preuves de variation.
Ceci est essentiel pour le calcul du génotypage conjoint ([lecture complémentaire](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

Le GVCF est produit par GATK HaplotypeCaller, le même outil que nous venons de tester, avec un paramètre supplémentaire (`-ERC GVCF`).
La combinaison des GVCF est effectuée avec GATK GenomicsDBImport, qui combine les appels par échantillon dans un magasin de données (analogue à une base de données).
L'analyse de 'génotypage conjoint' proprement dite est ensuite effectuée avec GATK GenotypeGVCFs.

Ici, nous testons les commandes nécessaires pour générer des GVCF et exécuter le génotypage conjoint.
Ce sont les commandes que nous encapsulerons dans un workflow Nextflow dans la Partie 3 de ce cours.

1. Générer un fichier d'index pour chaque fichier d'entrée BAM en utilisant Samtools
2. Exécuter GATK HaplotypeCaller sur chaque fichier d'entrée BAM pour générer un GVCF d'appels de variants génomiques par échantillon
3. Collecter tous les GVCF et les combiner dans un magasin de données GenomicsDB
4. Exécuter le génotypage conjoint sur le magasin de données GVCF combiné pour produire un VCF au niveau de la cohorte

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Nous devons maintenant tester toutes ces commandes, en commençant par indexer les trois fichiers BAM.

### 2.1. Indexer les fichiers BAM pour les trois échantillons

Dans la première section ci-dessus, nous n'avons indexé qu'un seul fichier BAM.
Maintenant, nous devons indexer les trois échantillons pour que GATK HaplotypeCaller puisse les traiter.

#### 2.1.1. Lancer le conteneur Samtools en mode interactif

Nous avons déjà téléchargé l'image du conteneur Samtools, nous pouvons donc la lancer directement :

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Votre invite de commande change pour indiquer que vous êtes à l'intérieur du conteneur, avec le répertoire de données monté comme précédemment.

#### 2.1.2. Exécuter la commande d'indexation sur les trois échantillons

Exécutez la commande d'indexation sur chacun des trois fichiers BAM :

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "Contenu du répertoire"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Cela devrait produire les fichiers d'index dans le même répertoire que les fichiers BAM correspondants.

#### 2.1.3. Quitter le conteneur Samtools

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite de commande devrait être revenue à la normale.

### 2.2. Générer des GVCF pour les trois échantillons

Pour exécuter l'étape de génotypage conjoint, nous avons besoin de GVCF pour les trois échantillons.

#### 2.2.1. Lancer le conteneur GATK en mode interactif

Nous avons déjà téléchargé l'image du conteneur GATK précédemment, nous pouvons donc la lancer directement :

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Votre invite de commande change pour indiquer que vous êtes à l'intérieur du conteneur GATK.

#### 2.2.2. Exécuter la commande d'appel de variants avec l'option GVCF

Afin de produire un VCF génomique (GVCF), nous ajoutons l'option `-ERC GVCF` à la commande de base, ce qui active le mode GVCF de HaplotypeCaller.

Nous modifions également l'extension du fichier de sortie de `.vcf` en `.g.vcf`.
Ce n'est techniquement pas une exigence, mais c'est une convention fortement recommandée.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Sortie de la commande"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

Cela crée le fichier de sortie GVCF `reads_mother.g.vcf` dans le répertoire de travail actuel dans le conteneur.

Si vous utilisez `cat` pour voir son contenu, vous constaterez qu'il est beaucoup plus long que le VCF équivalent que nous avons généré dans la section 1. Vous ne pouvez même pas remonter jusqu'au début du fichier, et la plupart des lignes semblent très différentes de ce que nous avons vu dans le VCF.

??? abstract "Contenu du fichier"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Celles-ci représentent des régions non-variantes où le programme d'appel de variants n'a trouvé aucune preuve de variation, il a donc capturé quelques statistiques décrivant son niveau de confiance dans l'absence de variation.
Cela permet de distinguer entre deux situations très différentes : (1) il existe des données de bonne qualité montrant que l'échantillon est homozygote-référence, et (2) il n'y a pas suffisamment de bonnes données disponibles pour faire une détermination dans un sens ou dans l'autre.

Dans un GVCF, il y a généralement beaucoup de telles lignes non-variantes, avec un plus petit nombre d'enregistrements de variants dispersés parmi elles.
Essayez d'exécuter `head -176` sur le GVCF pour charger uniquement les 176 premières lignes du fichier afin de trouver un appel de variant réel.

??? abstract "Contenu du fichier"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

La deuxième ligne montre le premier enregistrement de variant dans le fichier, qui correspond au premier variant dans le fichier VCF que nous avons examiné précédemment.

Tout comme le VCF original, le fichier de sortie GVCF est également accompagné d'un fichier d'index, appelé `reads_mother.g.vcf.idx`.

#### 2.2.3. Répéter le processus sur les deux autres échantillons

Générez des GVCF pour les deux échantillons restants en exécutant les commandes ci-dessous, l'une après l'autre.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Une fois cela terminé, vous devriez avoir trois fichiers se terminant par `.g.vcf` dans votre répertoire actuel (un par échantillon) et leurs fichiers d'index respectifs se terminant par `.g.vcf.idx`.

Mais ne quittez pas le conteneur !
Nous allons utiliser le même conteneur dans l'étape suivante.

### 2.3. Exécuter le génotypage conjoint

Maintenant que nous avons tous les GVCF, nous pouvons essayer l'approche de génotypage conjoint pour générer des appels de variants pour une cohorte d'échantillons.
C'est une méthode en deux étapes qui consiste à combiner les données de tous les GVCF dans un magasin de données, puis à exécuter l'analyse de génotypage conjoint proprement dite pour générer le VCF final de variants appelés conjointement.

#### 2.3.1. Combiner tous les GVCF par échantillon

Cette première étape utilise un autre outil GATK, appelé GenomicsDBImport, pour combiner les données de tous les GVCF dans un magasin de données GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Sortie de la commande"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

La sortie de cette étape est en fait un répertoire contenant un ensemble de répertoires supplémentaires imbriqués contenant les données de variants combinées sous la forme de plusieurs fichiers différents.
Vous pouvez l'explorer mais vous verrez rapidement que ce format de magasin de données n'est pas destiné à être lu directement par des humains.

!!! note

    GATK inclut des outils qui permettent d'inspecter et d'extraire les données d'appel de variants du magasin de données selon les besoins.

#### 2.3.2. Exécuter l'analyse de génotypage conjoint proprement dite

Cette deuxième étape utilise encore un autre outil GATK, appelé GenotypeGVCFs, pour recalculer les statistiques de variants et les génotypes individuels à la lumière des données disponibles dans tous les échantillons de la cohorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Sortie de la commande"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

Cela crée le fichier de sortie VCF `family_trio.vcf` dans le répertoire de travail actuel dans le conteneur.
C'est un autre fichier raisonnablement petit, vous pouvez donc utiliser `cat` sur ce fichier pour voir son contenu, et remonter pour trouver les premières lignes de variants.

??? abstract "Contenu du fichier"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Cela ressemble au VCF que nous avons généré précédemment, sauf que cette fois nous avons des informations au niveau du génotype pour les trois échantillons.
Les trois dernières colonnes du fichier sont les blocs de génotype pour les échantillons, listés par ordre alphabétique.

Si nous examinons les génotypes appelés pour notre trio familial de test pour le tout premier variant, nous voyons que le père est hétérozygote-variant (`0/1`), et la mère et le fils sont tous deux homozygotes-variant (`1/1`).

C'est en fin de compte l'information que nous cherchons à extraire du jeu de données !

#### 2.3.3. Quitter le conteneur GATK

Pour quitter le conteneur, tapez `exit`.

```bash
exit
```

Votre invite de commande devrait être revenue à la normale.
Cela conclut le test manuel des commandes d'appel de variants.

---

### À retenir

Vous savez comment tester les commandes d'indexation Samtools et d'appel de variants GATK dans leurs conteneurs respectifs, y compris comment générer des GVCF et exécuter le génotypage conjoint sur plusieurs échantillons.

### Et ensuite ?

Apprenez à encapsuler ces mêmes commandes dans des workflows qui utilisent des conteneurs pour exécuter le travail.
