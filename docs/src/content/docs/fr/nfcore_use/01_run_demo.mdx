# Partie 1 : Exécuter un pipeline de démonstration

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette première partie du cours "Use nf-core", nous vous montrons comment trouver un pipeline nf-core et l'essayer en utilisant son profil de test intégré.

Nous allons utiliser un pipeline appelé nf-core/demo, maintenu par le projet nf-core dans le cadre de son inventaire de pipelines à des fins de démonstration et de formation.

Assurez-vous que votre répertoire de travail est défini sur `nfcore-use/` comme indiqué sur la page [Premiers pas](./00_orientation.md).

---

## 1. Trouver et récupérer le pipeline nf-core/demo

Commençons par localiser le pipeline nf-core/demo sur le site web du projet à l'adresse [nf-co.re](https://nf-co.re), qui centralise toutes les informations telles que : la documentation générale et les articles d'aide, la documentation pour chacun des pipelines, les articles de blog, les annonces d'événements, etc.

### 1.1. Trouver le pipeline sur le site web

Dans votre navigateur web, rendez-vous sur [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) et tapez `demo` dans la barre de recherche.

![résultats de recherche](./img/search-results.png)

Cliquez sur le nom du pipeline, `demo`, pour accéder à la page de documentation du pipeline.

Chaque pipeline publié dispose d'une page dédiée qui comprend les sections de documentation suivantes :

- **Introduction :** Une introduction et une vue d'ensemble du pipeline
- **Usage :** Des descriptions de la façon d'exécuter le pipeline
- **Parameters :** Les paramètres du pipeline regroupés avec leurs descriptions
- **Output :** Des descriptions et des exemples des fichiers de sortie attendus
- **Results :** Des exemples de fichiers de sortie générés à partir du jeu de données de test complet
- **Releases & Statistics :** L'historique des versions du pipeline et les statistiques

Chaque fois que vous envisagez d'adopter un nouveau pipeline, vous devriez d'abord lire attentivement sa documentation pour comprendre ce qu'il fait et comment il doit être configuré avant de tenter de l'exécuter.

Jetez-y un œil maintenant et voyez si vous pouvez découvrir :

- Quels outils le pipeline va exécuter (Consultez l'onglet : `Introduction`)
- Quelles entrées et quels paramètres le pipeline accepte ou requiert (Consultez l'onglet : `Parameters`)
- Quelles sont les sorties produites par le pipeline (Consultez l'onglet : `Output`)

#### 1.1.1. Vue d'ensemble du pipeline

L'onglet `Introduction` fournit une vue d'ensemble du pipeline, incluant une représentation visuelle (appelée carte de métro) et une liste des outils exécutés dans le cadre du pipeline.

![carte de métro du pipeline](./img/nf-core-demo-subway-cropped.png)

1. Contrôle qualité des lectures ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Rognage des adaptateurs et de la qualité ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Présentation du contrôle qualité des lectures brutes ([MULTIQC](http://multiqc.info/))
4. Génération d'un message textuel fantaisiste depuis une vache ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Exemple de ligne de commande

La documentation fournit également un exemple de fichier d'entrée (abordé plus loin) et un exemple de ligne de commande.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Vous remarquerez que l'exemple de commande ne spécifie PAS de fichier workflow, seulement la référence au dépôt du pipeline, `nf-core/demo`.

Lorsqu'il est invoqué de cette façon, Nextflow suppose que le code est organisé d'une certaine manière.
Récupérons le code afin de pouvoir examiner cette structure.

### 1.2. Récupérer le code du pipeline

Une fois que nous avons déterminé que le pipeline semble convenir à nos besoins, essayons-le.
Heureusement, Nextflow facilite la récupération de pipelines depuis des dépôts correctement formatés sans avoir à télécharger quoi que ce soit manuellement.

#### 1.2.1. Utiliser `nextflow pull`

Retournons au terminal et exécutons la commande suivante :

```bash
nextflow pull nf-core/demo
```

??? success "Sortie de la commande"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow effectue un `pull` du code du pipeline, ce qui signifie qu'il télécharge le dépôt complet sur votre disque local.

Pour être clair, vous pouvez faire cela avec n'importe quel pipeline Nextflow correctement configuré sur GitHub, pas seulement les pipelines nf-core.
Cependant, nf-core est la plus grande collection open-source de pipelines Nextflow.

#### 1.2.2. Utiliser `nextflow list`

Vous pouvez demander à Nextflow de vous fournir une liste des pipelines que vous avez récupérés de cette façon :

```bash
nextflow list
```

??? success "Sortie de la commande"

    ```console
    nf-core/demo
    ```

Vous pouvez essayer de récupérer quelques autres pipelines pour voir comment ils sont listés lorsque vous en avez plus d'un.

#### 1.2.3. Trouver où le pipeline a été téléchargé

Vous remarquerez que les fichiers ne se trouvent pas dans votre répertoire de travail actuel.
Par défaut, Nextflow sauvegarde les pipelines récupérés sous `$NXF_HOME/assets`.

Pour trouver où se trouve un pipeline spécifique, demandez directement à Nextflow :

```bash
nextflow info nf-core/demo
```

??? success "Sortie de la commande"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Info"

    Le chemin complet peut différer sur votre système si vous n'utilisez pas notre environnement de formation.

Nextflow conserve le code source téléchargé intentionnellement "à l'écart" selon le principe que ces pipelines doivent être utilisés davantage comme des bibliothèques que comme du code avec lequel vous interagiriez directement.

En coulisses, Nextflow stocke chaque pipeline récupéré en tant que dépôt git sous `$NXF_HOME/assets/.repos/`, et extrait le code pour chaque révision dans un sous-répertoire `clones/<commit>/`.
Comme `.repos` est un répertoire caché, un simple `tree -L 2 $NXF_HOME/assets/` semblera vide.

#### 1.2.4. Créer un lien symbolique pour accéder facilement au code source

Nous n'allons pas examiner le code en détail, mais jetons-y un rapide coup d'œil pour avoir une idée de l'organisation générale.

Pour faciliter la navigation dans le code source du pipeline, créez un lien symbolique pointant vers la copie extraite du pipeline :

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Cela crée un raccourci pour vous permettre d'explorer le code avec `tree -L 2 pipelines/nf-core/demo` ou d'ouvrir des fichiers directement.

#### 1.2.5. Vue d'ensemble de l'organisation du code

Vous pouvez utiliser `tree` ou l'explorateur de fichiers pour trouver et ouvrir le répertoire `nf-core/demo`.

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

    7 directories, 12 files
    ```

Comme vous pouvez le voir, il s'y passe beaucoup de choses, dont la plupart ne vous concernent pas.

En bref, notons qu'au niveau supérieur, vous pouvez trouver un fichier README avec des informations récapitulatives, ainsi que des fichiers accessoires qui résument les informations du projet telles que la licence, les directives de contribution, les citations et le code de conduite.
La documentation détaillée du pipeline se trouve dans le répertoire `docs`.
Tout ce contenu est utilisé pour générer les pages web du site nf-core de manière programmatique, de sorte qu'elles sont toujours à jour avec le code.

Pour le reste, nous pouvons distinguer trois groupes fonctionnels de fichiers de code :

1. Composants du code du pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuration du pipeline
3. Paramètres / entrées et validation du pipeline

Nous n'allons pas passer en revue les composants du code du pipeline dans cette partie du cours, mais nous aborderons des éléments de configuration et de validation qui vous seront probablement utiles en tant qu'utilisateur·trice final·e des pipelines nf-core.

!!! tip "Astuce"

    Vous pouvez également parcourir le code source de n'importe quel pipeline nf-core sur GitHub, par exemple [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Chaque pipeline nf-core suit la même organisation de répertoires, donc une fois que vous connaissez la structure, vous pouvez trouver les fichiers de configuration, les modules et les workflows de n'importe quel pipeline de la même façon.

Pour l'instant, passons à l'exécution du pipeline !

### À retenir

Vous savez maintenant comment trouver un pipeline via le site web nf-core et récupérer une copie locale du code source.

### Et ensuite ?

Apprenez comment essayer un pipeline nf-core avec un minimum d'effort.

---

## 2. Essayer le pipeline avec son profil de test

Chaque pipeline nf-core est livré avec un profil de test, ce qui est très pratique.
Il s'agit d'un ensemble minimal de paramètres de configuration permettant au pipeline de s'exécuter en utilisant un petit jeu de données de test hébergé dans le dépôt [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
C'est un excellent moyen d'essayer rapidement un pipeline à petite échelle.

!!! tip "Astuce"

    Le système de profils de configuration de Nextflow vous permet de basculer facilement entre différents moteurs de conteneurs ou environnements d'exécution.
    Pour plus de détails, consultez [Hello Nextflow Partie 6 : Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Examiner le profil de test

Il est recommandé de vérifier ce que spécifie le profil de test d'un pipeline avant de l'exécuter.
Le profil `test` pour `nf-core/demo` se trouve dans le fichier de configuration `conf/test.config`.
Vous pouvez le trouver localement dans le code source du pipeline téléchargé par `nextflow pull`, via le lien symbolique `pipelines` créé à la section 1.2.4 :

```bash
code pipelines/nf-core/demo/conf/test.config
```

Voici le contenu de ce fichier :

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
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Données d'entrée
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Vous remarquerez immédiatement que le bloc de commentaires en haut inclut un exemple d'utilisation montrant comment exécuter le pipeline avec ce profil de test.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Les seules choses que nous devons fournir sont ce qui est indiqué entre chevrons dans l'exemple de commande : `<docker/singularity>` et `<OUTDIR>`.

Pour rappel, `<docker/singularity>` fait référence au choix du système de conteneurs. Tous les pipelines nf-core sont conçus pour être utilisables avec des conteneurs (Docker, Singularity, etc.) afin de garantir la reproductibilité et d'éliminer les problèmes d'installation de logiciels.
Nous devrons donc spécifier si nous souhaitons utiliser Docker ou Singularity pour tester le pipeline.

La partie `--outdir <OUTDIR>` fait référence au répertoire dans lequel Nextflow écrira les sorties du pipeline.
Nous devons lui fournir un nom, que nous pouvons simplement inventer.
S'il n'existe pas déjà, Nextflow le créera pour nous au moment de l'exécution.

En passant à la section après le bloc de commentaires, le profil de test nous montre ce qui a été préconfiguré pour les tests : notamment, le paramètre `input` est déjà défini pour pointer vers un jeu de données de test, nous n'avons donc pas besoin de fournir nos propres données.
Si vous suivez le lien vers l'entrée préconfigurée, vous verrez qu'il s'agit d'un fichier CSV contenant des identifiants d'échantillons et des chemins de fichiers pour plusieurs échantillons expérimentaux.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

C'est ce qu'on appelle un samplesheet, et c'est la forme d'entrée la plus courante pour les pipelines nf-core.
Ne vous inquiétez pas si vous n'êtes pas familier·ère avec les formats et types de données, ce n'est pas important pour la suite.

Nous avons maintenant tout ce qu'il nous faut pour essayer le pipeline.

### 2.2. Exécuter le pipeline

Comme indiqué ci-dessus, nous pouvons utiliser l'exemple de commande de test presque tel quel ; nous devons simplement spécifier quel système de gestion de logiciels utiliser et comment nommer le répertoire de sortie.
Ici, nous utiliserons Docker pour le système de conteneurs et `demo-results`, respectivement.

Avec cela, nous pouvons exécuter la commande de test :

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Si votre sortie correspond à cela, félicitations ! Vous venez d'exécuter votre premier pipeline nf-core.

Vous remarquerez qu'il y a beaucoup plus de sortie console que lorsque vous exécutez un pipeline Nextflow basique.
Il y a un en-tête qui inclut un résumé de la version du pipeline, des entrées et sorties, et quelques éléments de configuration.

!!! info "Info"

    Votre sortie affichera des horodatages, des noms d'exécution et des chemins de fichiers différents, mais la structure générale et l'exécution des processus devraient être similaires.

Remarquez la ligne près du haut de la sortie :

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Cela vous indique quelle révision du pipeline a été utilisée.
Comme nous n'avons pas spécifié de version, Nextflow a utilisé le dernier commit sur `master`.
Pour des exécutions reproductibles, vous devriez fixer une version spécifique en utilisant le flag `-r` :

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

Cela garantit que le même code de pipeline est utilisé à chaque fois, indépendamment des nouveaux commits ou des nouvelles versions.
Pour cette formation, nous omettons `-r` par souci de simplicité, mais en production vous devriez toujours le spécifier.

Passons à la sortie d'exécution et examinons les lignes qui nous indiquent quels processus ont été exécutés :

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Cela nous indique que quatre processus ont été exécutés, correspondant aux quatre outils présentés sur la page de documentation du pipeline sur le site web nf-core : `FASTQC`, `SEQTK_TRIM`, `MULTIQC` et `COWPY`.

Les noms complets des processus tels qu'ils apparaissent ici, comme `NFCORE_DEMO:DEMO:MULTIQC`, sont plus longs que ce que vous avez peut-être vu dans le matériel d'introduction Hello Nextflow.
Ceux-ci incluent les noms de leurs workflows parents et reflètent la modularité du code du pipeline.
Si vous souhaitez apprendre à développer vous-même des pipelines de style nf-core, consultez le cours [Build with nf-core](../nfcore_build/index.md).

### 2.3. Examiner les sorties du pipeline

Enfin, jetons un coup d'œil au répertoire `demo-results` produit par le pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contenu du répertoire"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
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
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Cela peut sembler beaucoup.
Pour en savoir plus sur les sorties du pipeline `nf-core/demo`, consultez sa [page de documentation](https://nf-co.re/demo/1.2.0/docs/output/).

À ce stade, ce qu'il est important d'observer, c'est que les résultats sont organisés par module, et qu'il existe en outre un répertoire appelé `pipeline_info` contenant divers rapports horodatés sur l'exécution du pipeline.

Par exemple, le fichier `execution_timeline_*` vous montre quels processus ont été exécutés, dans quel ordre et combien de temps ils ont pris :

![rapport de chronologie d'exécution](./img/execution_timeline.png)

!!! info "Info"

    Ici, les tâches n'ont pas été exécutées en parallèle car nous fonctionnons sur une machine minimaliste dans Github Codespaces.
    Pour les voir s'exécuter en parallèle, essayez d'augmenter l'allocation CPU de votre codespace et les limites de ressources dans la configuration de test.

Ces rapports sont générés automatiquement pour tous les pipelines nf-core.

### À retenir

Vous savez comment exécuter un pipeline nf-core en utilisant son profil de test intégré et où trouver ses sorties.

### Et ensuite ?

Rendez-vous à la [Partie 2](./02_configure_execution.md), où vous apprendrez comment configurer l'exécution du pipeline.

---

## Résumé

Dans cette partie, vous avez appris à :

- Trouver et récupérer un pipeline nf-core et examiner sa structure de code
- Exécuter un pipeline en utilisant son profil de test intégré
