# Partie 3 : Implémentation multi-échantillons en lecture appariée

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Précédemment, vous avez construit un pipeline d'appel de variants par échantillon qui traitait les données de chaque échantillon de manière indépendante.
Dans cette partie du cours, nous allons faire passer notre simple workflow au niveau supérieur en le transformant en un puissant outil d'automatisation par lots capable de traiter un nombre arbitraire d'échantillons.
Et pendant que nous y sommes, nous allons également le modifier pour qu'il accepte des données en lecture appariée, ce qui est plus courant dans les études récentes.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez terminé la [Partie 1 : Vue d'ensemble de la méthode](./01_method.md), la [Partie 2 : Implémentation mono-échantillon](./02_single-sample.md) et que vous disposez d'un pipeline `rnaseq.nf` fonctionnel avec des fichiers de modules complétés.

    Si vous n'avez pas terminé la Partie 2 ou souhaitez repartir de zéro pour cette partie, vous pouvez utiliser la solution de la Partie 2 comme point de départ.
    Exécutez ces commandes depuis le répertoire `nf4-science/rnaseq/` :

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Cela vous donne un workflow complet de traitement mono-échantillon.
    Vous pouvez tester qu'il s'exécute avec succès :

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Objectif

Dans cette partie du cours, nous allons étendre le workflow pour effectuer les opérations suivantes :

1. Lire les informations d'échantillons à partir d'un fichier CSV
2. Exécuter le QC, le trimming et l'alignement par échantillon sur tous les échantillons en parallèle
3. Agréger tous les rapports QC dans un rapport MultiQC complet

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Cela automatise les étapes de la deuxième section de la [Partie 1 : Vue d'ensemble de la méthode](./01_method.md#2-multi-sample-qc-aggregation), où vous avez exécuté ces commandes manuellement dans leurs conteneurs.

## Plan de la leçon

Nous avons divisé cela en trois étapes :

1. **Faire accepter au workflow plusieurs échantillons en entrée.**
   Cela couvre le passage d'un seul chemin de fichier à un fichier CSV d'échantillons, son analyse avec `splitCsv()`, et l'exécution de tous les processus existants sur plusieurs échantillons.
2. **Ajouter la génération de rapport QC complet.**
   Cela introduit l'opérateur `collect()` pour agréger les sorties entre échantillons, et ajoute un processus MultiQC pour produire un rapport combiné.
3. **Passer aux données RNAseq en lecture appariée.**
   Cela couvre l'adaptation des processus pour les entrées en lecture appariée (en utilisant des tuples), la création de modules en lecture appariée, et la configuration d'un profil de test séparé.

Cela implémente la méthode décrite dans la [Partie 1 : Vue d'ensemble de la méthode](./01_method.md) (deuxième section couvrant le cas d'usage multi-échantillons) et s'appuie directement sur le workflow produit par la Partie 2.

!!! tip "Astuce"

     Assurez-vous d'être dans le bon répertoire de travail :
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Faire accepter au workflow plusieurs échantillons en entrée

Pour exécuter sur plusieurs échantillons, nous devons modifier la façon dont nous gérons l'entrée : au lieu de fournir un seul chemin de fichier, nous allons lire les informations d'échantillons à partir d'un fichier CSV.

Nous fournissons un fichier CSV contenant les identifiants d'échantillons et les chemins de fichiers FASTQ dans le répertoire `data/`.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Ce fichier CSV inclut une ligne d'en-tête qui nomme les colonnes.

Notez qu'il s'agit toujours de données de lecture simple.

!!! warning "Avertissement"

    Les chemins de fichiers dans le CSV sont des chemins absolus qui doivent correspondre à votre environnement.
    Si vous n'exécutez pas cela dans l'environnement de formation que nous fournissons, vous devrez mettre à jour les chemins pour qu'ils correspondent à votre système.

### 1.1. Changer l'entrée principale pour qu'elle soit un CSV de chemins de fichiers dans le profil de test

Tout d'abord, nous devons mettre à jour le profil de test dans `nextflow.config` pour fournir le chemin du fichier CSV au lieu du chemin FASTQ unique.

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Ensuite, nous devrons mettre à jour la création du canal pour lire à partir de ce CSV.

### 1.2. Mettre à jour la fabrique de canal pour analyser l'entrée CSV

Nous devons charger le contenu du fichier dans le canal plutôt que simplement le chemin du fichier lui-même.

Nous pouvons le faire en utilisant le même modèle que nous avons utilisé dans la [Partie 2 de Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file) : appliquer l'opérateur [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) pour analyser le fichier, puis une opération `map` pour extraire le chemin du fichier FASTQ de chaque ligne.

=== "Après"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Créer le canal d'entrée à partir du contenu d'un fichier CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Créer le canal d'entrée à partir d'un chemin de fichier
        read_ch = channel.fromPath(params.input)
    ```

Une chose qui est nouvelle par rapport à ce que vous avez rencontré dans le cours Hello Nextflow est que ce CSV a une ligne d'en-tête, donc nous ajoutons `#!groovy header: true` à l'appel `splitCsv()`.
Cela nous permet de référencer les colonnes par nom dans l'opération `map` : `#!groovy row.fastq_path` extrait le chemin du fichier de la colonne `fastq_path` de chaque ligne.

La gestion de l'entrée est mise à jour et le workflow est prêt à être testé.

### 1.3. Exécuter le workflow

Le workflow lit maintenant les informations d'échantillons à partir d'un fichier CSV et traite tous les échantillons en parallèle.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Cette fois, chaque étape est exécutée 6 fois, une fois pour chaque échantillon dans le fichier CSV.

C'est tout ce qu'il a fallu pour que le workflow s'exécute sur plusieurs fichiers.
Nextflow gère tout le parallélisme pour nous.

### À retenir

Vous savez comment passer d'une entrée mono-fichier à une entrée multi-échantillons basée sur CSV que Nextflow traite en parallèle.

### Et ensuite ?

Ajouter une étape d'agrégation de rapport QC qui combine les métriques de tous les échantillons.

---

## 2. Agréger les métriques QC de pré-traitement dans un seul rapport MultiQC

Tout cela produit beaucoup de rapports QC, et nous ne voulons pas avoir à fouiller dans les rapports individuels.
C'est le moment idéal pour ajouter une étape d'agrégation de rapport MultiQC.

Rappelez-vous la commande `multiqc` de la [Partie 1](01_method.md) :

```bash
multiqc . -n <output_name>.html
```

La commande analyse le répertoire courant à la recherche de fichiers de sortie QC reconnus et les agrège dans un seul rapport HTML.
L'URI du conteneur était `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Nous devons configurer un paramètre supplémentaire, préparer les entrées, écrire le processus, le connecter et mettre à jour la gestion des sorties.

### 2.1. Configurer les entrées

Le processus MultiQC a besoin d'un paramètre de nom de rapport et des sorties QC collectées de toutes les étapes précédentes regroupées ensemble.

#### 2.1.1. Ajouter un paramètre `report_id`

Ajoutez un paramètre pour nommer le rapport de sortie.

=== "Après"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Entrée principale
        input: Path

        // Archive du génome de référence
        hisat2_index_zip: Path

        // ID de rapport
        report_id: String
    }
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Entrée principale
        input: Path

        // Archive du génome de référence
        hisat2_index_zip: Path
    }
    ```

Ajoutez la valeur par défaut de l'ID de rapport au profil de test :

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Ensuite, nous devrons préparer les entrées pour le processus MultiQC.

#### 2.1.2. Collecter et combiner les sorties QC des étapes précédentes

Nous devons donner au processus `MULTIQC` toutes les sorties liées au QC des étapes précédentes regroupées ensemble.

Pour cela, nous utilisons l'opérateur `.mix()`, qui agrège plusieurs canaux en un seul.
Nous partons de `channel.empty()` et y mélangeons tous les canaux de sortie que nous voulons combiner.
C'est plus propre que de chaîner `.mix()` directement sur l'un des canaux de sortie, car cela traite toutes les entrées de manière symétrique.

Dans notre workflow, les sorties liées au QC à agréger sont :

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Nous les mélangeons dans un seul canal, puis utilisons `.collect()` pour agréger les rapports de tous les échantillons en une seule liste.

Ajoutez ces lignes au corps du workflow après l'appel `HISAT2_ALIGN` :

=== "Après"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Alignement sur un génome de référence
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Génération de rapport QC complet
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="38"
        // Alignement sur un génome de référence
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

L'utilisation de variables intermédiaires rend chaque étape claire : `multiqc_files_ch` contient tous les fichiers QC individuels mélangés dans un canal, et `multiqc_files_list` est le paquet collecté prêt à être transmis à MultiQC.

### 2.2. Écrire le processus d'agrégation QC et l'appeler dans le workflow

Comme précédemment, nous devons compléter la définition du processus, importer le module et ajouter l'appel du processus.

#### 2.2.1. Compléter le module pour le processus d'agrégation QC

Ouvrez `modules/multiqc.nf` et examinez le squelette de la définition du processus.

Allez-y et complétez la définition du processus par vous-même en utilisant les informations fournies ci-dessus, puis vérifiez votre travail par rapport à la solution dans l'onglet "Après" ci-dessous.

=== "Avant"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Agréger les rapports QC avec MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Après"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Agréger les rapports QC avec MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Ce processus utilise `#!groovy path '*'` comme qualificateur d'entrée pour les fichiers QC.
Le joker `'*'` indique à Nextflow de placer tous les fichiers collectés dans le répertoire de travail sans nécessiter de noms spécifiques.
L'entrée `val output_name` est une chaîne qui contrôle le nom du fichier de rapport.

La commande `multiqc .` analyse le répertoire courant (où se trouvent tous les fichiers QC placés) et génère le rapport.

Une fois que vous avez terminé cela, le processus est prêt à être utilisé.

#### 2.2.2. Inclure le module

Ajoutez l'instruction d'importation à `rnaseq.nf` :

=== "Après"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Instructions INCLUDE de modules
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="3"
    // Instructions INCLUDE de modules
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Maintenant, ajoutez l'appel du processus au workflow.

#### 2.2.3. Ajouter l'appel du processus

Transmettez les fichiers QC collectés et l'ID de rapport au processus `MULTIQC` :

=== "Après"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

Le processus MultiQC est maintenant connecté au workflow.

### 2.3. Mettre à jour la gestion des sorties

Nous devons ajouter les sorties MultiQC à la déclaration de publication et configurer où elles vont.

#### 2.3.1. Ajouter des cibles de publication pour les sorties MultiQC

Ajoutez les sorties MultiQC à la section `publish:` :

=== "Après"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

Ensuite, nous devrons indiquer à Nextflow où placer ces sorties.

#### 2.3.2. Configurer les nouvelles cibles de sortie

Ajoutez des entrées pour les cibles MultiQC dans le bloc `output {}`, en les publiant dans un sous-répertoire `multiqc/` :

=== "Après"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Avant"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

La configuration de sortie est complète.

### 2.4. Exécuter le workflow

Nous utilisons `-resume` pour que les étapes de traitement précédentes soient mises en cache et que seule la nouvelle étape MultiQC s'exécute.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Un seul appel à MULTIQC a été ajouté après les appels de processus mis en cache.

Vous pouvez trouver les sorties MultiQC dans le répertoire des résultats.

```bash
tree -L 2 results/multiqc
```

```console title="Sortie"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Ce dernier fichier `all_single-end.html` est le rapport agrégé complet, commodément emballé dans un seul fichier HTML facile à consulter.

### À retenir

Vous savez comment collecter des sorties de plusieurs canaux, les regrouper avec `.mix()` et `.collect()`, et les transmettre à un processus d'agrégation.

### Et ensuite ?

Adapter le workflow pour gérer des données RNAseq en lecture appariée.

---

## 3. Activer le traitement de données RNAseq en lecture appariée

Actuellement, notre workflow ne peut gérer que des données RNAseq en lecture simple.
Il est de plus en plus courant de voir des données RNAseq en lecture appariée, donc nous voulons pouvoir les gérer.

Rendre le workflow complètement agnostique du type de données nécessiterait d'utiliser des fonctionnalités légèrement plus avancées du langage Nextflow, donc nous n'allons pas le faire ici, mais nous pouvons créer une version de traitement en lecture appariée pour démontrer ce qui doit être adapté.

### 3.1. Copier le workflow et mettre à jour les entrées

Nous commençons par copier le fichier de workflow en lecture simple et le mettre à jour pour les données en lecture appariée.

#### 3.1.1. Copier le fichier de workflow

Créez une copie du fichier de workflow à utiliser comme point de départ pour la version en lecture appariée.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Maintenant, mettez à jour les paramètres et la gestion des entrées dans le nouveau fichier.

#### 3.1.2. Ajouter un profil de test en lecture appariée

Nous fournissons un deuxième fichier CSV contenant les identifiants d'échantillons et les chemins de fichiers FASTQ appariés dans le répertoire `data/`.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Ajoutez un profil `test_pe` à `nextflow.config` qui pointe vers ce fichier et utilise un ID de rapport en lecture appariée.

=== "Après"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Avant"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

Le profil de test pour les données en lecture appariée est prêt.

#### 3.1.3. Mettre à jour la fabrique de canal

L'opérateur `.map()` doit maintenant récupérer les deux chemins de fichiers FASTQ et les retourner sous forme de liste.

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Créer le canal d'entrée à partir du contenu d'un fichier CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Créer le canal d'entrée à partir du contenu d'un fichier CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

La gestion des entrées est configurée pour les données en lecture appariée.

### 3.2. Adapter le module FASTQC pour les données en lecture appariée

Copiez le module pour créer une version en lecture appariée :

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

L'entrée du processus FASTQC n'a pas besoin de changer — lorsque Nextflow reçoit une liste de deux fichiers, il place les deux et `reads` se développe en les deux noms de fichiers.
Le seul changement nécessaire est dans le bloc de sortie : puisque nous obtenons maintenant deux rapports FastQC par échantillon, nous passons de modèles basés sur `simpleName` à des jokers.

=== "Après"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Avant"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Cela généralise le processus d'une manière qui le rend capable de gérer des données RNAseq en lecture simple ou appariée.

Mettez à jour l'importation dans `rnaseq_pe.nf` pour utiliser la version en lecture appariée :

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

Le module FASTQC et son importation sont mis à jour pour les données en lecture appariée.

### 3.3. Adapter le module TRIM_GALORE pour les données en lecture appariée

Copiez le module pour créer une version en lecture appariée :

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Ce module nécessite des changements plus substantiels :

- L'entrée passe d'un seul chemin à un tuple de deux chemins
- La commande ajoute le drapeau `--paired` et prend les deux fichiers de lecture
- La sortie change pour refléter les conventions de nommage en lecture appariée de Trim Galore, produisant des rapports FastQC séparés pour chaque fichier de lecture

=== "Après"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Avant"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

Mettez à jour l'importation dans `rnaseq_pe.nf` :

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Le module TRIM_GALORE et son importation sont mis à jour pour les données en lecture appariée.

### 3.4. Adapter le module HISAT2_ALIGN pour les données en lecture appariée

Copiez le module pour créer une version en lecture appariée :

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Ce module nécessite des changements similaires :

- L'entrée passe d'un seul chemin à un tuple de deux chemins
- La commande HISAT2 passe de `-U` (non apparié) aux arguments de lecture `-1` et `-2` (appariés)
- Toutes les utilisations de `reads.simpleName` changent en `read1.simpleName` puisque nous référençons maintenant un membre spécifique de la paire

=== "Après"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Avant"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

Mettez à jour l'importation dans `rnaseq_pe.nf` :

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Le module HISAT2_ALIGN et son importation sont mis à jour pour les données en lecture appariée.

### 3.5. Mettre à jour l'agrégation MultiQC pour les sorties en lecture appariée

Le processus `TRIM_GALORE` en lecture appariée produit maintenant deux canaux de rapports FastQC séparés (`fastqc_reports_1` et `fastqc_reports_2`) au lieu d'un seul.
Mettez à jour le bloc `.mix()` dans `rnaseq_pe.nf` pour inclure les deux :

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

L'agrégation MultiQC inclut maintenant les deux ensembles de rapports FastQC en lecture appariée.

### 3.6. Mettre à jour la gestion des sorties pour les sorties en lecture appariée

La section `publish:` et le bloc `output {}` doivent également refléter les deux canaux de rapports FastQC séparés du processus `TRIM_GALORE` en lecture appariée.

Mettez à jour la section `publish:` dans `rnaseq_pe.nf` :

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Mettez à jour les entrées correspondantes dans le bloc `output {}` :

=== "Après"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Avant"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

Le workflow en lecture appariée est maintenant entièrement mis à jour et prêt à être exécuté.

### 3.7. Exécuter le workflow

Nous n'utilisons pas `-resume` car cela ne serait pas mis en cache, et il y a deux fois plus de données à traiter qu'avant, mais cela devrait quand même se terminer en moins d'une minute.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Maintenant nous avons deux versions légèrement divergentes de notre workflow, une pour les données de lecture simple et une pour les données de lecture appariée.
L'étape logique suivante serait de faire en sorte que le workflow accepte l'un ou l'autre type de données à la volée, ce qui est hors du cadre de ce cours, mais nous pourrions aborder cela dans une suite.

---

### À retenir

Vous savez comment adapter un workflow mono-échantillon pour paralléliser le traitement de plusieurs échantillons, générer un rapport QC complet et adapter le workflow pour utiliser des données de lecture appariée.

### Et ensuite ?

Félicitez-vous chaleureusement ! Vous avez terminé le cours Nextflow pour RNAseq.

Rendez-vous au [résumé du cours](./next_steps.md) final pour passer en revue ce que vous avez appris et découvrir ce qui vient ensuite.
