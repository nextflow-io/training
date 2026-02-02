# Partie 1 : Exécuter un pipeline de démonstration

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette première partie de la formation Hello nf-core, nous vous montrons comment trouver et essayer un pipeline nf-core, comprendre comment le code est organisé, et reconnaître en quoi il diffère du code Nextflow simple tel que présenté dans [Hello Nextflow](../hello_nextflow/index.md).

Nous allons utiliser un pipeline appelé nf-core/demo qui est maintenu par le projet nf-core dans le cadre de son inventaire de pipelines pour démontrer la structure du code et le fonctionnement des outils.

Assurez-vous que votre répertoire de travail est défini sur `hello-nf-core/` comme indiqué sur la page [Premiers pas](./00_orientation.md).

---

## 1. Trouver et récupérer le pipeline nf-core/demo

Commençons par localiser le pipeline nf-core/demo sur le site web du projet à l'adresse [nf-co.re](https://nf-co.re), qui centralise toutes les informations telles que : la documentation générale et les articles d'aide, la documentation pour chacun des pipelines, les articles de blog, les annonces d'événements, etc.

### 1.1. Trouver le pipeline sur le site web

Dans votre navigateur web, allez sur [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) et tapez `demo` dans la barre de recherche.

![résultats de recherche](./img/search-results.png)

Cliquez sur le nom du pipeline, `demo`, pour accéder à la page de documentation du pipeline.

Chaque pipeline publié dispose d'une page dédiée qui comprend les sections de documentation suivantes :

- **Introduction :** Une introduction et un aperçu du pipeline
- **Usage :** Descriptions de la manière d'exécuter le pipeline
- **Parameters :** Paramètres du pipeline regroupés avec descriptions
- **Output :** Descriptions et exemples des fichiers de sortie attendus
- **Results :** Exemples de fichiers de sortie générés à partir du jeu de données de test complet
- **Releases & Statistics :** Historique des versions du pipeline et statistiques

Chaque fois que vous envisagez d'adopter un nouveau pipeline, vous devez d'abord lire attentivement la documentation du pipeline pour comprendre ce qu'il fait et comment il doit être configuré avant de tenter de l'exécuter.

Jetez-y un coup d'œil maintenant et voyez si vous pouvez découvrir :

- Quels outils le pipeline exécutera (Vérifiez l'onglet : `Introduction`)
- Quelles entrées et paramètres le pipeline accepte ou requiert (Vérifiez l'onglet : `Parameters`)
- Quelles sont les sorties produites par le pipeline (Vérifiez l'onglet : `Output`)

#### 1.1.1. Aperçu du pipeline

L'onglet `Introduction` fournit un aperçu du pipeline, incluant une représentation visuelle (appelée carte de métro) et une liste des outils qui sont exécutés dans le cadre du pipeline.

![carte de métro du pipeline](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Exemple de ligne de commande

La documentation fournit également un exemple de fichier d'entrée (discuté plus en détail ci-dessous) et un exemple de ligne de commande.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Vous remarquerez que l'exemple de commande ne spécifie PAS de fichier de workflow, juste la référence au dépôt du pipeline, `nf-core/demo`.

Lorsqu'il est invoqué de cette manière, Nextflow supposera que le code est organisé d'une certaine manière.
Récupérons le code afin de pouvoir examiner cette structure.

### 1.2. Récupérer le code du pipeline

Une fois que nous avons déterminé que le pipeline semble convenir à nos besoins, essayons-le.
Heureusement, Nextflow facilite la récupération des pipelines à partir de dépôts correctement formatés sans avoir à télécharger quoi que ce soit manuellement.

Retournons au terminal et exécutons ce qui suit :

```bash
nextflow pull nf-core/demo
```

??? success "Sortie de la commande"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow effectue un `pull` du code du pipeline, c'est-à-dire qu'il télécharge le dépôt complet sur votre disque local.

Pour être clair, vous pouvez faire cela avec n'importe quel pipeline Nextflow qui est correctement configuré dans GitHub, pas seulement les pipelines nf-core.
Cependant, nf-core est la plus grande collection open-source de pipelines Nextflow.

Vous pouvez demander à Nextflow de vous fournir une liste des pipelines que vous avez récupérés de cette manière :

```bash
nextflow list
```

??? success "Sortie de la commande"

    ```console
    nf-core/demo
    ```

Vous remarquerez que les fichiers ne se trouvent pas dans votre répertoire de travail actuel.
Par défaut, Nextflow les enregistre dans `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Contenu du répertoire"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    Le chemin complet peut différer sur votre système si vous n'utilisez pas notre environnement de formation.

Nextflow garde intentionnellement le code source téléchargé 'à l'écart' sur le principe que ces pipelines doivent être utilisés davantage comme des bibliothèques que comme du code avec lequel vous interagiriez directement.

Cependant, pour les besoins de cette formation, nous voulons pouvoir explorer et voir ce qu'il contient.
Donc, pour faciliter cela, créons un lien symbolique vers cet emplacement depuis notre répertoire de travail actuel.

```bash
ln -s $NXF_HOME/assets pipelines
```

Cela crée un raccourci qui facilite l'exploration du code que nous venons de télécharger.

```bash
tree -L 2 pipelines
```

```console title="Contenu du répertoire"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Maintenant, nous pouvons plus facilement jeter un coup d'œil au code source si nécessaire.

Mais d'abord, essayons d'exécuter notre premier pipeline nf-core !

### À retenir

Vous savez maintenant comment trouver un pipeline via le site web nf-core et récupérer une copie locale du code source.

### Et ensuite ?

Apprenez comment essayer un pipeline nf-core avec un minimum d'effort.

---

## 2. Essayer le pipeline avec son profil de test

De manière pratique, chaque pipeline nf-core est fourni avec un profil de test.
Il s'agit d'un ensemble minimal de paramètres de configuration permettant au pipeline de s'exécuter en utilisant un petit jeu de données de test hébergé dans le dépôt [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
C'est un excellent moyen d'essayer rapidement un pipeline à petite échelle.

!!! note

    Le système de profils de configuration de Nextflow vous permet de basculer facilement entre différents moteurs de conteneurs ou environnements d'exécution.
    Pour plus de détails, consultez [Hello Nextflow Partie 6 : Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Examiner le profil de test

C'est une bonne pratique de vérifier ce que spécifie le profil de test d'un pipeline avant de l'exécuter.
Le profil `test` pour `nf-core/demo` se trouve dans le fichier de configuration `conf/test.config` et est présenté ci-dessous.

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Données d'entrée
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Vous remarquerez tout de suite que le bloc de commentaires en haut inclut un exemple d'utilisation montrant comment exécuter le pipeline avec ce profil de test.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Les seules choses que nous devons fournir sont ce qui est montré entre crochets dans l'exemple de commande : `<docker/singularity>` et `<OUTDIR>`.

Pour rappel, `<docker/singularity>` fait référence au choix du système de conteneurs. Tous les pipelines nf-core sont conçus pour être utilisables avec des conteneurs (Docker, Singularity, etc.) afin de garantir la reproductibilité et d'éliminer les problèmes d'installation de logiciels.
Nous devrons donc spécifier si nous voulons utiliser Docker ou Singularity pour tester le pipeline.

La partie `--outdir <OUTDIR>` fait référence au répertoire où Nextflow écrira les sorties du pipeline.
Nous devons lui fournir un nom, que nous pouvons simplement inventer.
S'il n'existe pas déjà, Nextflow le créera pour nous lors de l'exécution.

Passant à la section après le bloc de commentaires, le profil de test nous montre ce qui a été préconfiguré pour les tests : plus particulièrement, le paramètre `input` est déjà défini pour pointer vers un jeu de données de test, nous n'avons donc pas besoin de fournir nos propres données.
Si vous suivez le lien vers l'entrée préconfigurée, vous verrez qu'il s'agit d'un fichier csv contenant des identifiants d'échantillons et des chemins de fichiers pour plusieurs échantillons expérimentaux.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

C'est ce qu'on appelle une feuille d'échantillons (samplesheet), et c'est la forme d'entrée la plus courante pour les pipelines nf-core.

!!! note

    Ne vous inquiétez pas si vous n'êtes pas familier avec les formats et types de données, ce n'est pas important pour la suite.

Cela confirme donc que nous avons tout ce dont nous avons besoin pour essayer le pipeline.

### 2.2. Exécuter le pipeline

Décidons d'utiliser Docker pour le système de conteneurs et `demo-results` comme répertoire de sortie, et nous sommes prêts à exécuter la commande de test :

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Si votre sortie correspond à celle-ci, félicitations ! Vous venez d'exécuter votre premier pipeline nf-core.

Vous remarquerez qu'il y a beaucoup plus de sortie console que lorsque vous exécutez un pipeline Nextflow basique.
Il y a un en-tête qui inclut un résumé de la version du pipeline, des entrées et sorties, et quelques éléments de configuration.

!!! note

    Votre sortie affichera des horodatages, des noms d'exécution et des chemins de fichiers différents, mais la structure globale et l'exécution des processus devraient être similaires.

Passons maintenant à la sortie d'exécution, et jetons un coup d'œil aux lignes qui nous indiquent quels processus ont été exécutés :

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Cela nous indique que trois processus ont été exécutés, correspondant aux trois outils présentés dans la page de documentation du pipeline sur le site web nf-core : FASTQC, SEQTK_TRIM et MULTIQC.

Les noms complets des processus tels qu'affichés ici, comme `NFCORE_DEMO:DEMO:MULTIQC`, sont plus longs que ce que vous avez pu voir dans le matériel d'introduction Hello Nextflow.
Ils incluent les noms de leurs workflows parents et reflètent la modularité du code du pipeline.
Nous entrerons plus en détail à ce sujet dans un instant.

### 2.3. Examiner les sorties du pipeline

Enfin, jetons un coup d'œil au répertoire `demo-results` produit par le pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contenu du répertoire"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

Cela peut sembler beaucoup.
Pour en savoir plus sur les sorties du pipeline `nf-core/demo`, consultez sa [page de documentation](https://nf-co.re/demo/1.0.2/docs/output/).

À ce stade, ce qui est important d'observer est que les résultats sont organisés par module, et il y a en plus un répertoire appelé `pipeline_info` contenant divers rapports horodatés sur l'exécution du pipeline.

Par exemple, le fichier `execution_timeline_*` vous montre quels processus ont été exécutés, dans quel ordre et combien de temps ils ont pris pour s'exécuter :

![rapport de chronologie d'exécution](./img/execution_timeline.png)

!!! note

    Ici, les tâches n'ont pas été exécutées en parallèle car nous fonctionnons sur une machine minimaliste dans Github Codespaces.
    Pour voir ces tâches s'exécuter en parallèle, essayez d'augmenter l'allocation CPU de votre codespace et les limites de ressources dans la configuration de test.

Ces rapports sont générés automatiquement pour tous les pipelines nf-core.

### À retenir

Vous savez comment exécuter un pipeline nf-core en utilisant son profil de test intégré et où trouver ses sorties.

### Et ensuite ?

Apprenez comment le code du pipeline est organisé.

---

## 3. Examiner la structure du code du pipeline

Maintenant que nous avons exécuté avec succès le pipeline en tant qu'utilisateur·trices, changeons de perspective pour voir comment les pipelines nf-core sont structurés en interne.

Le projet nf-core applique des directives strictes sur la façon dont les pipelines sont structurés, et sur la façon dont le code est organisé, configuré et documenté.
Comprendre comment tout cela est organisé est la première étape vers le développement de vos propres pipelines compatibles nf-core, que nous aborderons dans la Partie 2 de ce cours.

Jetons un coup d'œil à la façon dont le code du pipeline est organisé dans le dépôt `nf-core/demo`, en utilisant le lien symbolique `pipelines` que nous avons créé précédemment.

Vous pouvez soit utiliser `tree` soit utiliser l'explorateur de fichiers pour trouver et ouvrir le répertoire `nf-core/demo`.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Contenu du répertoire"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

Il se passe beaucoup de choses là-dedans, nous allons donc aborder cela étape par étape.

Tout d'abord, notons qu'au niveau supérieur, vous pouvez trouver un fichier README avec des informations récapitulatives, ainsi que des fichiers accessoires qui résument les informations du projet telles que les licences, les directives de contribution, les citations et le code de conduite.
La documentation détaillée du pipeline se trouve dans le répertoire `docs`.
Tout ce contenu est utilisé pour générer les pages web sur le site web nf-core de manière programmatique, elles sont donc toujours à jour avec le code.

Maintenant, pour le reste, nous allons diviser notre exploration en trois étapes :

1. Composants du code du pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuration du pipeline
3. Entrées et validation

Commençons par les composants du code du pipeline.
Nous allons nous concentrer sur la hiérarchie des fichiers et l'organisation structurelle, plutôt que de plonger dans le code à l'intérieur de fichiers individuels.

### 3.1. Composants du code du pipeline

L'organisation standard du code d'un pipeline nf-core suit une structure modulaire conçue pour maximiser la réutilisation du code, comme introduit dans [Hello Modules](../hello_nextflow/04_hello_modules.md), Partie 4 du cours [Hello Nextflow](../hello_nextflow/index.md), bien que dans le véritable style nf-core, cela soit implémenté avec un peu de complexité supplémentaire.
Plus précisément, les pipelines nf-core font un usage abondant des subworkflows, c'est-à-dire des scripts de workflow qui sont importés par un workflow parent.

Cela peut sembler un peu abstrait, alors jetons un coup d'œil à la façon dont cela est utilisé en pratique dans le pipeline `nf-core/demo`.

!!! note

    Nous ne passerons pas en revue le code réel de la _façon_ dont ces composants modulaires sont connectés, car il y a une complexité supplémentaire associée à l'utilisation des subworkflows qui peut être déroutante, et comprendre cela n'est pas nécessaire à ce stade de la formation.
    Pour l'instant, nous allons nous concentrer sur l'organisation générale et la logique.

#### 3.1.1. Vue d'ensemble générale

Voici à quoi ressemblent les relations entre les composants de code pertinents pour le pipeline `nf-core/demo` :

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

Il y a un script dit _point d'entrée_ appelé `main.nf`, qui agit comme une enveloppe pour deux types de workflows imbriqués : le workflow contenant la logique d'analyse réelle, situé sous `workflows/` et appelé `demo.nf`, et un ensemble de workflows de gestion situés sous `subworkflows/`.
Le workflow `demo.nf` fait appel aux **modules** situés sous `modules/` ; ceux-ci contiennent les **processus** qui effectueront les étapes d'analyse réelles.

!!! note

    Les subworkflows ne sont pas limités aux fonctions de gestion, et ils peuvent utiliser des modules de processus.

    Le pipeline `nf-core/demo` présenté ici se trouve être du côté le plus simple du spectre, mais d'autres pipelines nf-core (tels que `nf-core/rnaseq`) utilisent des subworkflows qui sont impliqués dans l'analyse réelle.

Maintenant, passons en revue ces composants à tour de rôle.

#### 3.1.2. Le script de point d'entrée : `main.nf`

Le script `main.nf` est le point d'entrée à partir duquel Nextflow démarre lorsque nous exécutons `nextflow run nf-core/demo`.
Cela signifie que lorsque vous exécutez `nextflow run nf-core/demo` pour exécuter le pipeline, Nextflow trouve et exécute automatiquement le script `main.nf`.
Cela fonctionne pour tout pipeline Nextflow qui suit cette convention de nommage et cette structure, pas seulement les pipelines nf-core.

L'utilisation d'un script de point d'entrée facilite l'exécution de subworkflows de 'gestion' standardisés avant et après l'exécution du script d'analyse réel.
Nous examinerons ceux-ci après avoir passé en revue le workflow d'analyse réel et ses modules.

#### 3.1.3. Le script d'analyse : `workflows/demo.nf`

Le workflow `workflows/demo.nf` est l'endroit où la logique centrale du pipeline est stockée.
Il est structuré de la même manière qu'un workflow Nextflow normal, sauf qu'il est conçu pour être appelé depuis un workflow parent, ce qui nécessite quelques fonctionnalités supplémentaires.
Nous couvrirons les différences pertinentes dans la partie suivante de ce cours, lorsque nous aborderons la conversion du simple pipeline Hello de Hello Nextflow en une forme compatible nf-core.

Le workflow `demo.nf` fait appel aux **modules** situés sous `modules/`, que nous examinerons ensuite.

!!! note

    Certains workflows d'analyse nf-core affichent des niveaux supplémentaires d'imbrication en appelant des subworkflows de niveau inférieur.
    Ceci est principalement utilisé pour encapsuler deux modules ou plus qui sont couramment utilisés ensemble dans des segments de pipeline facilement réutilisables.
    Vous pouvez voir quelques exemples en parcourant les [subworkflows nf-core](https://nf-co.re/subworkflows/) disponibles sur le site web nf-core.

    Lorsque le script d'analyse utilise des subworkflows, ceux-ci sont stockés sous le répertoire `subworkflows/`.

#### 3.1.4. Les modules

Les modules sont l'endroit où réside le code des processus, comme décrit dans la [Partie 4 du cours de formation Hello Nextflow](../hello_nextflow/04_hello_modules.md).

Dans le projet nf-core, les modules sont organisés selon une structure imbriquée à plusieurs niveaux qui reflète à la fois leur origine et leur contenu.
Au niveau supérieur, les modules sont différenciés comme étant soit `nf-core` soit `local` (ne faisant pas partie du projet nf-core), puis placés dans un répertoire nommé d'après le(s) outil(s) qu'ils encapsulent.
Si l'outil appartient à une boîte à outils (c'est-à-dire un package contenant plusieurs outils), il y a un niveau de répertoire intermédiaire nommé d'après la boîte à outils.

Vous pouvez voir cela appliqué en pratique aux modules du pipeline `nf-core/demo` :

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Contenu du répertoire"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Ici, vous voyez que les modules `fastqc` et `multiqc` se situent au niveau supérieur au sein des modules `nf-core`, alors que le module `trim` se trouve sous la boîte à outils à laquelle il appartient, `seqtk`.
Dans ce cas, il n'y a pas de modules `local`.

Le fichier de code du module décrivant le processus s'appelle toujours `main.nf`, et est accompagné de tests et de fichiers `.yml` que nous ignorerons pour l'instant.

Pris ensemble, le workflow de point d'entrée, le workflow d'analyse et les modules sont suffisants pour exécuter les parties 'intéressantes' du pipeline.
Cependant, nous savons qu'il y a aussi des subworkflows de gestion là-dedans, alors regardons-les maintenant.

#### 3.1.5. Les subworkflows de gestion

Comme les modules, les subworkflows sont différenciés en répertoires `local` et `nf-core`, et chaque subworkflow a sa propre structure de répertoire imbriquée avec son propre script `main.nf`, ses tests et son fichier `.yml`.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Contenu du répertoire"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Comme noté ci-dessus, le pipeline `nf-core/demo` n'inclut aucun subworkflow spécifique à l'analyse, donc tous les subworkflows que nous voyons ici sont des workflows dits de 'gestion' ou 'utilitaires', comme indiqué par le préfixe `utils_` dans leurs noms.
Ces subworkflows sont ce qui produit le joli en-tête nf-core dans la sortie console, parmi d'autres fonctions accessoires.

!!! tip

    Mis à part leur modèle de nommage, une autre indication que ces subworkflows n'effectuent aucune fonction véritablement liée à l'analyse est qu'ils n'appellent aucun processus du tout.

Ceci complète le tour d'horizon des composants de code de base qui constituent le pipeline `nf-core/demo`.
Maintenant, jetons un coup d'œil aux éléments restants que vous devriez connaître un peu avant de plonger dans le développement : la configuration du pipeline et la validation des entrées.

### 3.2. Configuration du pipeline

Vous avez appris précédemment que Nextflow offre de nombreuses options pour configurer l'exécution du pipeline, que ce soit en termes d'entrées et de paramètres, de ressources de calcul et d'autres aspects d'orchestration.
Le projet nf-core applique des directives hautement standardisées pour la configuration du pipeline qui visent à s'appuyer sur les options de personnalisation flexibles de Nextflow d'une manière qui offre une plus grande cohérence et maintenabilité entre les pipelines.

Le fichier de configuration central `nextflow.config` est utilisé pour définir les valeurs par défaut des paramètres et d'autres options de configuration.
La majorité de ces options de configuration sont appliquées par défaut tandis que d'autres (par exemple, les profils de dépendances logicielles) sont incluses comme profils optionnels.

Il existe plusieurs fichiers de configuration supplémentaires qui sont stockés dans le dossier `conf` et qui peuvent être ajoutés à la configuration par défaut ou optionnellement comme profils :

- `base.config` : Un fichier de configuration 'de base', approprié pour une utilisation générale sur la plupart des environnements de calcul haute performance. Cela définit de larges catégories d'utilisation des ressources, par exemple, qui sont pratiques à appliquer aux modules.
- `modules.config` : Directives et arguments de modules supplémentaires.
- `test.config` : Un profil pour exécuter le pipeline avec des données de test minimales, que nous avons utilisé lorsque nous avons exécuté le pipeline de démonstration.
- `test_full.config` : Un profil pour exécuter le pipeline avec un jeu de données de test de taille complète.

Nous toucherons quelques-uns de ces fichiers plus tard dans le cours.

### 3.3. Entrées et validation

Comme nous l'avons noté précédemment, lorsque nous avons examiné le profil de test du pipeline `nf-core/demo`, il est conçu pour prendre comme entrée une feuille d'échantillons contenant des chemins de fichiers et des identifiants d'échantillons.
Les chemins de fichiers sont liés à de vraies données situées dans le dépôt `nf-core/test-datasets`.

Un exemple de feuille d'échantillons est également fourni dans le répertoire `assets`, bien que les chemins dans celui-ci ne soient pas réels.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Cette feuille d'échantillons particulière est assez simple, mais certains pipelines s'exécutent sur des feuilles d'échantillons plus complexes, avec beaucoup plus de métadonnées associées aux entrées principales.

Malheureusement, parce que ces fichiers peuvent être difficiles à vérifier à l'œil, un formatage incorrect des données d'entrée est une source très courante d'échecs de pipeline.
Un problème connexe est lorsque les paramètres sont fournis de manière incorrecte.

La solution à ces problèmes est d'exécuter des vérifications de validation automatisées sur tous les fichiers d'entrée pour s'assurer qu'ils contiennent les types d'informations attendus, correctement formatés, et sur les paramètres pour s'assurer qu'ils sont du type attendu.
C'est ce qu'on appelle la validation d'entrée, et devrait idéalement être effectuée _avant_ d'essayer d'exécuter un pipeline, plutôt que d'attendre que le pipeline échoue pour découvrir qu'il y avait un problème avec les entrées.

Tout comme pour la configuration, le projet nf-core a des opinions très arrêtées sur la validation des entrées, et recommande l'utilisation du [plugin nf-schema](https://nextflow-io.github.io/nf-schema/latest/), un plugin Nextflow qui fournit des capacités de validation complètes pour les pipelines Nextflow.

Nous couvrirons ce sujet plus en détail dans la Partie 5 de ce cours.
Pour l'instant, soyez simplement conscient qu'il existe deux fichiers JSON fournis à cet effet, `nextflow_schema.json` et `assets/schema_input.json`.

Le `nextflow_schema.json` est un fichier utilisé pour stocker des informations sur les paramètres du pipeline, y compris le type, la description et le texte d'aide dans un format lisible par machine.
Ceci est utilisé à diverses fins, notamment la validation automatisée des paramètres, la génération de texte d'aide et le rendu de formulaires de paramètres interactifs dans les interfaces utilisateur.

Le `schema_input.json` est un fichier utilisé pour définir la structure de la feuille d'échantillons d'entrée.
Chaque colonne peut avoir un type, un modèle, une description et un texte d'aide dans un format lisible par machine.
Le schéma est utilisé à diverses fins, notamment la validation automatisée et la fourniture de messages d'erreur utiles.

### À retenir

Vous savez quels sont les principaux composants d'un pipeline nf-core et comment le code est organisé ; où se trouvent les principaux éléments de configuration ; et vous êtes conscient·e de l'utilité de la validation des entrées.

### Et ensuite ?

Faites une pause ! C'était beaucoup. Lorsque vous vous sentez rafraîchi·e et prêt·e, passez à la section suivante pour appliquer ce que vous avez appris afin d'écrire un pipeline compatible nf-core.

!!! tip

    Si vous souhaitez apprendre à composer des workflows avec des subworkflows avant de passer à la partie suivante, consultez la [Quête Secondaire Workflows de Workflows](../side_quests/workflows_of_workflows.md).
