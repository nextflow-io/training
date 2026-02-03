# Partie 2 : Implémentation pour un seul échantillon

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette partie du cours, nous allons écrire le workflow le plus simple possible qui encapsule toutes les commandes que nous avons exécutées dans la Partie 1 pour automatiser leur exécution, et nous viserons simplement à traiter un échantillon à la fois.

Nous le ferons en trois étapes :

1. Écrire un workflow à une seule étape qui exécute l'étape initiale de QC
2. Ajouter le nettoyage des adaptateurs et le QC post-nettoyage
3. Ajouter l'alignement sur le génome de référence

!!! warning "Prérequis"

    Vous devez travailler la Partie 1 du cours avant de commencer cette leçon.
    Plus précisément, travailler les sections 2.1-3 crée le fichier d'index du génome (`data/genome_index.tar.gz`) requis pour l'étape d'alignement dans cette leçon.

---

## 1. Écrire un workflow à une seule étape qui exécute le QC initial

Commençons par écrire un workflow simple qui exécute l'outil FastQC sur un fichier FASTQ contenant des lectures RNAseq single-end.

Nous vous fournissons un fichier de workflow, `rnaseq.nf`, qui décrit les parties principales du workflow.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Instructions INCLUDE de modules

/*
 * Paramètres du pipeline
 */

// Entrée principale

workflow {

    // Créer le canal d'entrée

    // Appeler les processus

}
```

Gardez à l'esprit que ce code de workflow est correct mais il n'est pas fonctionnel ; son but est juste de servir de squelette que vous utiliserez pour écrire le workflow réel.

### 1.1. Créer un répertoire pour stocker les modules

Nous allons créer des modules autonomes pour chaque processus afin de faciliter leur gestion et leur réutilisation, créons donc un répertoire pour les stocker.

```bash
mkdir modules
```

### 1.2. Créer un module pour le processus de collecte des métriques QC

Créons un fichier de module appelé `modules/fastqc.nf` pour héberger le processus `FASTQC` :

```bash
touch modules/fastqc.nf
```

Ouvrez le fichier dans l'éditeur de code et copiez-y le code suivant :

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

Vous devriez reconnaître tous les éléments de ce que vous avez appris dans les Parties 1 et 2 de cette série de formation ; le seul changement notable est que cette fois nous utilisons `mode: symlink` pour la directive `publishDir`, et nous utilisons un paramètre pour définir le `publishDir`.

!!! note

    Même si les fichiers de données que nous utilisons ici sont très petits, en génomique ils peuvent devenir très volumineux. À des fins de démonstration dans l'environnement d'enseignement, nous utilisons le mode de publication 'symlink' pour éviter des copies de fichiers inutiles. Vous ne devriez pas faire cela dans vos workflows finaux, car vous perdrez les résultats lorsque vous nettoierez votre répertoire `work`.

### 1.3. Importer le module dans le fichier de workflow

Ajoutez l'instruction `include { FASTQC } from './modules/fastqc.nf'` au fichier `rnaseq.nf` :

```groovy title="rnaseq.nf" linenums="3"
// Instructions INCLUDE de modules
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Ajouter une déclaration d'entrée

Déclarez un paramètre d'entrée avec une valeur par défaut :

```groovy title="rnaseq.nf" linenums="10"
params {
    // Entrée principale
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Créer un canal d'entrée dans le bloc workflow

Utilisez une fabrique de canaux basique `.fromPath()` pour créer le canal d'entrée :

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Créer le canal d'entrée à partir d'un chemin de fichier
    read_ch = channel.fromPath(params.reads)

    // Appeler les processus

}
```

### 1.6. Appeler le processus `FASTQC` sur le canal d'entrée

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Créer le canal d'entrée à partir d'un chemin de fichier
    read_ch = channel.fromPath(params.reads)

    // Contrôle qualité initial
    FASTQC(read_ch)

}
```

### 1.7. Exécuter le workflow pour tester qu'il fonctionne

Nous pourrions utiliser le paramètre `--reads` pour spécifier une entrée depuis la ligne de commande, mais pendant le développement nous pouvons être paresseux et simplement utiliser le test par défaut que nous avons configuré.

```bash
nextflow run rnaseq.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Cela devrait s'exécuter très rapidement si vous avez travaillé la Partie 1 et avez déjà téléchargé le conteneur.
Si vous l'avez sautée, Nextflow téléchargera le conteneur pour vous ; vous n'avez rien à faire pour que cela se produise, mais vous devrez peut-être attendre jusqu'à une minute.

Vous pouvez trouver les sorties sous `results/fastqc` comme spécifié dans le processus `FASTQC` par la directive `publishDir`.

```bash
ls results/fastqc
```

```console title="Sortie"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Ajouter le nettoyage des adaptateurs et le contrôle qualité post-nettoyage

Nous allons utiliser le wrapper Trim_Galore, qui intègre Cutadapt pour le nettoyage lui-même et FastQC pour le contrôle qualité post-nettoyage.

### 2.1. Créer un module pour le processus de nettoyage et QC

Créons un fichier de module appelé `modules/trim_galore.nf` pour héberger le processus `TRIM_GALORE` :

```bash
touch modules/trim_galore.nf
```

Ouvrez le fichier dans l'éditeur de code et copiez-y le code suivant :

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Importer le module dans le fichier de workflow

Ajoutez l'instruction `include { TRIM_GALORE } from './modules/trim_galore.nf'` au fichier `rnaseq.nf` :

```groovy title="rnaseq.nf" linenums="3"
// Instructions INCLUDE de modules
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Appeler le processus sur le canal d'entrée

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Créer le canal d'entrée à partir d'un chemin de fichier
    read_ch = channel.fromPath(params.reads)

    // Contrôle qualité initial
    FASTQC(read_ch)

    // Nettoyage des adaptateurs et QC post-nettoyage
    TRIM_GALORE(read_ch)
}
```

### 2.4. Exécuter le workflow pour tester qu'il fonctionne

```bash
nextflow run rnaseq.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

Cela devrait également s'exécuter très rapidement, car nous travaillons sur un fichier d'entrée si petit.

Vous pouvez trouver les sorties sous `results/trimming` comme spécifié dans le processus `TRIM_GALORE` par la directive `publishDir`.

```bash
ls results/trimming
```

```console title="Sortie"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Aligner les lectures sur le génome de référence

Finalement, nous pouvons exécuter l'étape d'alignement du génome en utilisant Hisat2, qui émettra également des métriques de contrôle qualité de type FastQC.

### 3.1. Créer un module pour le processus HiSat2

Créons un fichier de module appelé `modules/hisat2_align.nf` pour héberger le processus `HISAT2_ALIGN` :

```bash
touch modules/hisat2_align.nf
```

Ouvrez le fichier dans l'éditeur de code et copiez-y le code suivant :

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Importer le module dans le fichier de workflow

Ajoutez l'instruction `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` au fichier `rnaseq.nf` :

```groovy title="rnaseq.nf" linenums="3"
// Instructions INCLUDE de modules
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Ajouter une déclaration de paramètre pour fournir l'index du génome

Déclarez un paramètre d'entrée avec une valeur par défaut :

```groovy title="rnaseq.nf" linenums="8"
params {
    // Entrée principale
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Archive du génome de référence
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Appeler le processus `HISAT2_ALIGN` sur les lectures nettoyées produites par `TRIM_GALORE`

Les lectures nettoyées sont dans le canal de sortie `TRIM_GALORE.out.trimmed_reads` produit par l'étape précédente.

De plus, nous utilisons `file (params.hisat2_index_zip)` pour fournir à l'outil Hisat2 l'archive tar compressée de l'index du génome.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Créer le canal d'entrée à partir d'un chemin de fichier
    read_ch = channel.fromPath(params.reads)

    // Contrôle qualité initial
    FASTQC(read_ch)

    // Nettoyage des adaptateurs et QC post-nettoyage
    TRIM_GALORE(read_ch)

    // Alignement sur un génome de référence
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Exécuter le workflow pour tester qu'il fonctionne

```bash
nextflow run rnaseq.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Vous pouvez trouver les sorties sous `results/align` comme spécifié dans le processus `HISAT2_ALIGN` par la directive `publishDir`.

```bash
ls results/align
```

```console title="Sortie"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Cela complète le traitement de base que nous devons appliquer à chaque échantillon.

_Nous ajouterons l'agrégation de rapports MultiQC dans la Partie 2, après avoir fait en sorte que le workflow accepte plusieurs échantillons à la fois._

---

### À retenir

Vous savez comment encapsuler toutes les étapes principales pour traiter des échantillons RNAseq single-end individuellement.

### Et ensuite ?

Apprenez comment modifier le workflow pour traiter plusieurs échantillons en parallèle, agréger les rapports QC sur toutes les étapes pour tous les échantillons, et permettre l'exécution du workflow sur des données RNAseq paired-end.
