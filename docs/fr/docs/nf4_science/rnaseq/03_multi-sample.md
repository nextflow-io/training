# Partie 3 : Implémentation multi-échantillons en lecture appariée

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette dernière partie du cours, nous allons faire passer notre simple workflow au niveau supérieur en le transformant en un puissant outil d'automatisation par lots capable de traiter un nombre arbitraire d'échantillons.
Et pendant que nous y sommes, nous allons également le modifier pour qu'il accepte des données en lecture appariée, ce qui est plus courant dans les études récentes.

Nous procéderons en trois étapes :

1. Faire accepter au workflow plusieurs échantillons en entrée et paralléliser l'exécution
2. Ajouter la génération de rapport QC complet
3. Passer aux données RNAseq en lecture appariée

---

## 1. Faire accepter au workflow plusieurs échantillons en entrée et paralléliser l'exécution

Nous allons devoir modifier la façon dont nous gérons l'entrée.

### 1.1. Changer l'entrée principale pour qu'elle soit un CSV de chemins de fichiers au lieu d'un seul fichier

Nous fournissons un fichier CSV contenant les identifiants d'échantillons et les chemins de fichiers FASTQ dans le répertoire `data/`.
Ce fichier CSV inclut une ligne d'en-tête.
Notez que les chemins de fichiers FASTQ sont des chemins absolus.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Renommons le paramètre d'entrée principal en `input_csv` et changeons la valeur par défaut pour qu'elle soit le chemin vers le fichier `single-end.csv`.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Entrée principale
    input_csv: Path = "data/single-end.csv"

    // Archive de l'index du génome de référence
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Mettre à jour la fabrique de canal d'entrée pour gérer un CSV en entrée

Nous allons vouloir charger le contenu du fichier dans le canal plutôt que simplement le chemin du fichier lui-même, donc nous utilisons l'opérateur `.splitCsv()` pour analyser le format CSV, puis l'opérateur `.map()` pour récupérer l'information spécifique que nous voulons (le chemin du fichier FASTQ).

```groovy title="rnaseq.nf" linenums="16"
    // Créer le canal d'entrée à partir du contenu d'un fichier CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Exécuter le workflow pour vérifier qu'il fonctionne

```bash
nextflow run rnaseq.nf
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

Cette fois, nous voyons que chaque étape est exécutée 6 fois, sur chacun des 6 fichiers de données que nous avons fournis.

C'est tout ce qu'il a fallu pour que le workflow s'exécute sur plusieurs fichiers !
Nextflow gère tout le parallélisme pour nous.

---

## 2. Agréger les métriques QC de pré-traitement dans un seul rapport MultiQC

Tout cela produit beaucoup de rapports QC, et nous ne voulons pas avoir à fouiller dans les rapports individuels.
C'est le moment idéal pour ajouter une étape d'agrégation de rapport MultiQC !

### 2.1. Créer un module pour le processus d'agrégation QC

Créons un fichier module appelé `modules/multiqc.nf` pour héberger le processus `MULTIQC` :

```bash
touch modules/multiqc.nf
```

Ouvrez le fichier dans l'éditeur de code et copiez-y le code suivant :

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. Importer le module dans le fichier de workflow

Ajoutez l'instruction `include { MULTIQC } from './modules/multiqc.nf'` au fichier `rnaseq.nf` :

```groovy title="rnaseq.nf" linenums="3"
// Instructions INCLUDE de modules
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Ajouter un paramètre `report_id` et lui donner une valeur par défaut appropriée

```groovy title="rnaseq.nf" linenums="9"
params {
    // Entrée principale
    input_csv: Path = "data/single-end.csv"

    // Archive de l'index du génome de référence
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID de rapport
    report_id: String = "all_single-end"
}
```

### 2.4. Appeler le processus sur les sorties des étapes précédentes

Nous devons donner au processus `MULTIQC` toutes les sorties liées au QC des étapes précédentes.

Pour cela, nous allons utiliser l'opérateur `.mix()`, qui agrège plusieurs canaux en un seul.

Si nous avions quatre processus appelés A, B, C et D avec chacun un simple canal `.out`, la syntaxe ressemblerait à ceci : `A.out.mix( B.out, C.out, D.out )`. Comme vous pouvez le voir, vous l'appliquez au premier des canaux que vous voulez combiner (peu importe lequel) et vous ajoutez simplement tous les autres, séparés par des virgules, dans les parenthèses qui suivent.

Dans le cas de notre workflow, nous avons les sorties suivantes à agréger :

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

L'exemple de syntaxe devient donc :

```groovy title="Application de .mix() dans l'appel MULTIQC"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

Cela collectera les rapports QC par échantillon.
Mais puisque nous voulons les agréger pour tous les échantillons, nous devons ajouter l'opérateur `collect()` afin de rassembler les rapports de tous les échantillons en un seul appel à `MULTIQC`.
Et nous devons également lui donner le paramètre `report_id`.

Cela nous donne ce qui suit :

```groovy title="L'appel MULTIQC complet" linenums="33"
    // Génération de rapport QC complet
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Dans le contexte du bloc de workflow complet, cela finit par ressembler à ceci :

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Créer le canal d'entrée à partir du contenu d'un fichier CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Contrôle qualité initial
    FASTQC(read_ch)

    // Coupe des adaptateurs et QC post-coupe
    TRIM_GALORE(read_ch)

    // Alignement sur un génome de référence
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Génération de rapport QC complet
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Exécuter le workflow pour vérifier qu'il fonctionne

```bash
nextflow run rnaseq.nf -resume
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

Cette fois, nous voyons un seul appel à MULTIQC ajouté après les appels de processus mis en cache :

Vous pouvez trouver les sorties sous `results/trimming` comme spécifié dans le processus `TRIM_GALORE` par la directive `publishDir`.

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

---

## 3. Activer le traitement de données RNAseq en lecture appariée

Actuellement, notre workflow ne peut gérer que des données RNAseq en lecture simple.
Il est de plus en plus courant de voir des données RNAseq en lecture appariée, donc nous voulons pouvoir les gérer.

Rendre le workflow complètement agnostique du type de données nécessiterait d'utiliser des fonctionnalités légèrement plus avancées du langage Nextflow, donc nous n'allons pas le faire ici, mais nous pouvons créer une version de traitement en lecture appariée pour démontrer ce qui doit être adapté.

### 3.1. Créer une copie du workflow appelée `rnaseq_pe.nf`

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Modifier le `input_csv` par défaut pour pointer vers les données en lecture appariée

Nous fournissons un deuxième fichier CSV contenant les identifiants d'échantillons et les chemins de fichiers FASTQ appariés dans le répertoire `data/`

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Changeons la valeur par défaut de `input_csv` pour qu'elle soit le chemin vers le fichier `paired-end.csv`.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Entrée principale
    input_csv: Path = "data/paired-end.csv"

    // Archive de l'index du génome de référence
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID de rapport
    report_id: String = "all_single-end"
}
```

### 3.3. Mettre à jour la fabrique de canal

Nous devons indiquer à l'opérateur `.map()` de récupérer maintenant les deux chemins de fichiers FASTQ.

Donc `row -> file(row.fastq_path)` devient `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Créer le canal d'entrée à partir du contenu d'un fichier CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Créer une version en lecture appariée du processus FASTQC

Créons une copie du module pour que nous puissions avoir les deux versions à portée de main.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Ouvrez le nouveau fichier module `fastqc_pe.nf` dans l'éditeur de code et effectuez les modifications de code suivantes :

- Changez `fastqc $reads` en `fastqc ${reads}` dans le bloc `script` (ligne 17) pour que l'entrée `reads` soit déballée, puisqu'elle est maintenant un tuple de deux chemins au lieu d'un seul chemin.
- Remplacez `${reads.simpleName}` par un joker (`*`) pour éviter d'avoir à gérer les fichiers de sortie individuellement.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Techniquement, cela généralise le processus `FASTQC` d'une manière qui le rend capable de gérer des données RNAseq en lecture simple ou appariée.

Enfin, mettez à jour l'instruction d'importation de module pour utiliser la version en lecture appariée du module.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Créer une version en lecture appariée du processus TRIM_GALORE

Créez une copie du module pour que nous puissions avoir les deux versions à portée de main.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Ouvrez le nouveau fichier module `trim_galore_pe.nf` dans l'éditeur de code et effectuez les modifications de code suivantes :

- Changez la déclaration d'entrée de `path reads` en `tuple path(read1), path(read2)`
- Mettez à jour la commande dans le bloc `script`, en remplaçant `$reads` par `--paired ${read1} ${read2}`
- Mettez à jour les déclarations de sortie pour refléter les fichiers ajoutés et les différentes conventions de nommage, en utilisant des jokers pour éviter d'avoir à tout lister.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Enfin, mettez à jour l'instruction d'importation de module pour utiliser la version en lecture appariée du module.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Mettre à jour l'appel au processus MULTIQC pour attendre deux rapports de TRIM_GALORE

Le processus `TRIM_GALORE` produit maintenant un canal de sortie supplémentaire, donc nous devons le fournir à MultiQC.

Remplacez `TRIM_GALORE.out.fastqc_reports,` par `TRIM_GALORE.out.fastqc_reports_1,` plus `TRIM_GALORE.out.fastqc_reports_2,` :

```groovy title="rnaseq_pe.nf" linenums="33"
    // Génération de rapport QC complet
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Pendant que nous sommes sur MultiQC, mettons également à jour la valeur par défaut du paramètre `report_id` de `"all_single-end"` à `"all_paired-end"`.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Entrée principale
    input_csv: Path = "data/paired-end.csv"

    // Archive de l'index du génome de référence
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID de rapport
    report_id: String = "all_paired-end"
}
```

### 3.7. Créer une version en lecture appariée du processus HISAT2_ALIGN

Créez une copie du module pour que nous puissions avoir les deux versions à portée de main.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Ouvrez le nouveau fichier module `hisat2_align_pe.nf` dans l'éditeur de code et effectuez les modifications de code suivantes :

- Changez la déclaration d'entrée de `path reads` en `tuple path(read1), path(read2)`
- Mettez à jour la commande dans le bloc `script`, en remplaçant `-U $reads` par `-1 ${read1} -2 ${read2}`
- Remplacez toutes les instances de `${reads.simpleName}` par `${read1.simpleName}` dans la commande du bloc `script` ainsi que dans les déclarations de sortie.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Enfin, mettez à jour l'instruction d'importation de module pour utiliser la version en lecture appariée du module.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Exécuter le workflow pour vérifier qu'il fonctionne

Nous n'utilisons pas `-resume` car cela ne serait pas mis en cache, et il y a deux fois plus de données à traiter qu'avant, mais cela devrait quand même se terminer en moins d'une minute.

```bash
nextflow run rnaseq_pe.nf
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

Et voilà ! Maintenant nous avons deux versions légèrement divergentes de notre workflow, une pour les données de lecture simple et une pour les données de lecture appariée.
L'étape logique suivante serait de faire en sorte que le workflow accepte l'un ou l'autre type de données à la volée, ce qui est hors du cadre de ce cours, mais nous pourrions aborder cela dans une suite.

---

### À retenir

Vous savez comment adapter un workflow mono-échantillon pour paralléliser le traitement de plusieurs échantillons, générer un rapport QC complet et adapter le workflow pour utiliser des données de lecture appariée si nécessaire.

### Et maintenant ?

Félicitations, vous avez terminé le mini-cours Nextflow pour RNAseq ! Célébrez votre succès et prenez une pause bien méritée !

Ensuite, nous vous demandons de compléter une très courte enquête sur votre expérience avec ce cours de formation, puis nous vous emmènerons vers une page avec des liens vers d'autres ressources de formation et des liens utiles.
