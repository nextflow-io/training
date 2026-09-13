---
title: Partie 2 : Configurer l'exécution du pipeline
description: Apprenez à configurer l'exécution d'un pipeline nf-core en définissant des paramètres, en comprenant la validation et en personnalisant l'allocation des ressources.
---

# Partie 2 : Configurer l'exécution du pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans la [Partie 1](./01_run_demo.md), vous avez trouvé et exécuté le pipeline nf-core/demo en utilisant son profil de test.
Nous allons maintenant voir comment configurer l'exécution du pipeline : définir des paramètres, comprendre la validation, et personnaliser l'allocation des ressources et les arguments des outils.

Comme expliqué dans [Hello Config](../hello_nextflow/06_hello_config.md), nous souhaitons pouvoir modifier les données sur lesquelles notre pipeline s'exécutera et la façon dont il s'exécutera, sans modifier le code du pipeline lui-même.
À cette fin, Nextflow prend en charge plusieurs façons de contrôler la configuration du pipeline, ce qui peut être un peu déroutant.

Le projet nf-core définit des conventions pour organiser les éléments de configuration, en distinguant deux types de configuration au niveau supérieur : les **paramètres du pipeline** et la **configuration** au sens strict.

- Les **paramètres du pipeline** (définis via le système `params`) incluent généralement des éléments tels que les fichiers d'entrée, les indicateurs de comportement des outils et les paramètres d'analyse.
- La **configuration** au sens strict désigne la logistique de la façon dont le pipeline est exécuté, c'est-à-dire l'executor, l'allocation des ressources de calcul, etc.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Commençons par les paramètres du pipeline, puis nous examinerons la configuration au sens strict.

---

## 1. Paramètres du pipeline

Pour tous les pipelines nf-core, vous pouvez obtenir la liste complète des paramètres du pipeline directement depuis la ligne de commande en utilisant l'indicateur `--help`, qui est lui-même un paramètre du pipeline.

### 1.1. Obtenir la liste des paramètres avec `--help`

Exécutez la commande d'aide pour le pipeline demo :

```bash
nextflow run nf-core/demo --help
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Comme vous pouvez le voir, la sortie regroupe les paramètres par catégories (options d'entrée/sortie, options du génome de référence, etc.) avec les types et les descriptions pour chacun.

Cette catégorisation est déterminée par un fichier de schéma, qui est présenté plus en détail ci-dessous.
Dans les pipelines Nextflow simples, `--help` ne fonctionne que si le·la développeur·se l'a implémenté manuellement.

!!! tip "Astuce"

    Utilisez `--help --show_hidden` pour voir les paramètres supplémentaires qui sont masqués par défaut, tels que `--publish_dir_mode` ou `--monochrome_logs`.

### 1.2. Définir les valeurs des paramètres

Comme expliqué dans [Hello Config](../hello_nextflow/06_hello_config.md), vous pouvez définir les valeurs des paramètres sur la ligne de commande avec `--param_name` ou regrouper un ensemble de paramètres dans un fichier YAML et le passer avec `-params-file`.
Les deux approches fonctionnent de la même façon avec les pipelines nf-core.

Par exemple, pour ignorer l'étape de rognage, nous souhaitons définir le paramètre booléen `skip_trim` à `true`.
Un fichier de paramètres appelé `my_params.yml` est fourni dans votre répertoire de travail avec cette valeur déjà définie :

```yaml title="my_params.yml"
skip_trim: true
```

Passez-le avec `-params-file` :

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


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
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
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

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Le processus `SEQTK_TRIM` n'apparaît plus dans la sortie.

!!! warning "Avertissement : Limitations importantes concernant les entrées de paramètres"

    **Définir des paramètres booléens sur la ligne de commande**

    À partir de la version 26.04 de Nextflow, toutes les valeurs fournies sur la ligne de commande sont typées comme des chaînes de caractères.
    Pour un paramètre booléen comme `skip_trim`, le passer comme un indicateur nu (`--skip_trim`) ou comme `--skip_trim true` est évalué comme la **chaîne** `"true"`, ce qui échoue à la validation du schéma :

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Pour définir un paramètre booléen à une valeur `true`/`false` authentique, utilisez un `-params-file` comme indiqué ci-dessus, ou définissez-le dans un fichier de configuration.
    Les paramètres de type chaîne, entier et chemin de fichier ne sont pas affectés et peuvent toujours être définis directement sur la ligne de commande.
    Ce cours utilise ce modèle tout au long pour les paramètres booléens.

    **Utiliser des fichiers de configuration personnalisés**

    Bien qu'il soit techniquement possible de définir des paramètres du pipeline dans un fichier de configuration personnalisé passé avec `-c`, cela peut ne pas remplacer les valeurs par défaut déjà définies dans le `nextflow.config` du pipeline, selon les règles de priorité de configuration de Nextflow.
    Utiliser `--param_name` sur la ligne de commande ou `-params-file` est plus fiable, car ces méthodes ont toujours la priorité.

    En règle générale : si un paramètre apparaît dans la sortie de `--help`, définissez-le via la ligne de commande ou un fichier de paramètres plutôt que via un fichier de configuration.

### 1.3. Validation des paramètres

Fait intéressant : la commande `--help` fonctionne pour tous les pipelines nf-core parce que le projet nf-core exige des développeur·ses de définir formellement tous les paramètres du pipeline dans un fichier de schéma JSON (`nextflow_schema.json`).
Ce schéma enregistre le type, la description, la valeur par défaut et le regroupement de chaque paramètre.

En plus d'alimenter la sortie de `--help`, le fichier de schéma permet également une validation automatisée au moment du lancement.
Cela signifie que Nextflow peut vérifier que chaque paramètre que vous passez existe et a reçu une valeur appropriée (du type approprié, dans la plage de valeurs autorisée, etc.).

Nous abordons cela plus en détail dans la [section sur la validation des entrées](../nfcore_build/04_input_validation.md), mais vous pouvez déjà le voir en action en fournissant au pipeline demo des entrées de paramètres invalides.

#### 1.3.1. Paramètres non reconnus

Essayez de passer un paramètre qui n'existe pas :

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

La sortie de la console inclut un avertissement :

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

Le pipeline s'exécute quand même, mais l'avertissement vous alerte immédiatement que `--foobar` n'est pas un paramètre reconnu.
Cela est destiné à attirer votre attention sur les fautes de frappe non bloquantes, comme l'utilisation de `--outDir` au lieu de `--outdir`, ce qui peut vous aider à éviter de perdre du temps et des ressources de calcul.

#### 1.3.2. Valeurs de paramètres invalides

La validation vérifie également les **valeurs** des paramètres.
Le paramètre `--skip_trim` est un indicateur booléen, donc passer une valeur de type chaîne provoque l'échec immédiat du pipeline :

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Le pipeline s'arrête avant l'exécution de tout processus, vous évitant ainsi une exécution échouée ou incorrecte.
Comme indiqué dans la [section 1.2](#12-set-parameter-values), les paramètres booléens doivent être définis à une valeur `true`/`false` authentique dans un fichier de paramètres plutôt que passés sur la ligne de commande, car les valeurs de la ligne de commande sont typées comme des chaînes de caractères.

### 1.4. Validation des entrées

La même logique de validation peut également être utilisée pour vérifier la validité des fichiers d'entrée.
Par exemple, si un pipeline attend une samplesheet comme entrée de données principale (ce qui est le cas de nombreux pipelines nf-core, voire de la plupart), le·la développeur·se peut fournir un schéma d'entrée (distinct du schéma des paramètres) décrivant comment le fichier d'entrée doit être structuré.

Ensuite, au moment de l'exécution, Nextflow peut vérifier que le fichier d'entrée fourni est valide.

Nous abordons également cela plus en détail dans la [section sur la validation des entrées](../nfcore_build/04_input_validation.md), mais vous pouvez déjà le voir en action en fournissant au pipeline demo une samplesheet d'entrée invalide.

Le pipeline `nf-core/demo` attend un fichier CSV avec les colonnes `sample`, `fastq_1` et `fastq_2`.
Cela est défini dans un fichier de schéma (`assets/schema_input.json`) qui spécifie la structure attendue, les types de colonnes et les contraintes.

??? abstract "Fichier de schéma pour les entrées"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Le schéma spécifie que `sample` et `fastq_1` sont obligatoires, tandis que `fastq_2` est optionnel (prenant en charge les données paired-end et single-end).
Les chemins de fichiers sont validés pour leur existence et leur extension.

Pour illustrer cela, nous fournissons une samplesheet malformée appelée `malformed_samplesheet.csv` dans votre répertoire de travail :

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Cette samplesheet est dépourvue de la colonne obligatoire `fastq_1` et contient un chemin de fichier inexistant dans `fastq_2`.

Exécutez le pipeline demo en utilisant `malformed_samplesheet.csv` comme entrée :

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Comme vous pouvez le voir, le pipeline échoue immédiatement et signale **toutes** les erreurs de validation en même temps.
nf-schema ne s'arrête pas à la première erreur — il collecte tous les problèmes et les liste ensemble, afin que vous puissiez tout corriger en une seule fois plutôt que de découvrir les problèmes un par un.

Chaque erreur identifie l'entrée et le champ exacts qui ont causé le problème, afin que vous puissiez corriger votre samplesheet puis relancer le pipeline avec la certitude qu'il ne va pas échouer à un moment ultérieur lorsque Nextflow tentera d'accéder au chemin du fichier.

Pour les développeur·ses, tout cela est abordé plus en détail dans la [Partie 4 de Build with nf-core](../nfcore_build/04_input_validation.md).

### À retenir

Vous savez comment obtenir la liste complète des paramètres d'un pipeline avec `--help`, les définir via la ligne de commande ou un fichier de paramètres, et comment Nextflow valide à la fois les valeurs des paramètres et les fichiers d'entrée par rapport aux schémas du pipeline.

### Et ensuite ?

Découvrez l'autre type de configuration : la façon dont le pipeline s'exécute, couvrant l'allocation des ressources et les arguments des outils.

---

## 2. Configuration

La configuration au sens strict contrôle **comment** le pipeline s'exécute : l'allocation des ressources, les arguments spécifiques aux outils, l'endroit où les tâches s'exécutent et le système d'empaquetage logiciel à utiliser.

Les pipelines nf-core incluent une configuration par défaut dans `nextflow.config` et le répertoire `conf/`.
Avant de remplacer quoi que ce soit, il est utile de savoir où se trouvent les valeurs par défaut.

### 2.1. Explorer les fichiers de configuration

Vous avez déjà vu dans la [Partie 1](./01_run_demo.md) que le code source du pipeline se trouve sous `$NXF_HOME/assets`.
En utilisant le lien symbolique `pipelines` que vous avez créé dans la [Partie 1](./01_run_demo.md), listez les fichiers de configuration pour voir ce qui est disponible :

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

Les fichiers de configuration les plus importants sont :

- **`conf/base.config`** : Définit les labels de ressources (`process_low`, `process_medium`, `process_high`) qui attribuent des CPUs, de la mémoire et du temps aux processus. Lorsque vous voyez un processus utiliser plus de ressources que prévu, c'est là que proviennent ces valeurs par défaut.
- **`conf/modules.config`** : Définit les arguments des outils par processus (`ext.args`) et les paramètres de publication des sorties (`publishDir`). Ouvrez ce fichier pour voir quels arguments chaque outil reçoit par défaut.
- **`conf/test.config`** : Le profil de test que vous avez utilisé dans la [Partie 1](./01_run_demo.md), qui limite les ressources via `resourceLimits` et définit une samplesheet de test. Activé avec `-profile test`.
  Il existe également un `conf/test_full.config` pour l'exécution avec un jeu de données de test de taille réelle, utile pour les benchmarks.

Le `nextflow.config` central charge tous les fichiers ci-dessus et définit les valeurs par défaut appropriées pour tout.

Si vous souhaitez modifier l'un des paramètres spécifiés dans ces fichiers, ne modifiez aucun d'entre eux directement.
Créez plutôt votre propre fichier de configuration et passez-le avec `-c`.
Les valeurs que vous spécifiez remplaceront les valeurs par défaut définies dans ces autres fichiers.

Essayons cela en pratique.

### 2.2. Personnaliser les ressources des processus et les arguments des outils

Les modules nf-core prennent en charge deux types courants de remplacement de configuration : **l'allocation des ressources** (CPUs, mémoire, temps) et les **arguments des outils** via `ext.args`.

De nombreux outils en ligne de commande ont des arguments qui ne sont pas suffisamment utilisés pour être exposés comme paramètres du pipeline.
La convention `ext.args` vous permet de passer ces arguments à l'outil sous-jacent via un fichier de configuration.

Le fichier `custom.config` fourni dans votre répertoire de travail illustre ces deux types de remplacement :

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

Le premier bloc remplace l'allocation des ressources de `FASTQC`.
Par défaut, `FASTQC` utilise le label `process_medium` de `base.config`, qui alloue 6 CPUs et 36 Go de mémoire ; ici nous le limitons à 2 CPUs et 4 Go.

Le second bloc passe un argument supplémentaire à `SEQTK_TRIM` via `ext.args`.
L'indicateur `-b 5` indique à `seqtk trimfq` de rogner 5 bases au début de chaque lecture en plus du rognage par qualité.

Exécutez le pipeline avec cette configuration :

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Sortie de la commande"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

L'indicateur `-c` ajoute votre configuration par-dessus la configuration intégrée du pipeline.

Pour vérifier que le remplacement de `ext.args` a bien pris effet, trouvez le hash du répertoire de travail de `SEQTK_TRIM` dans la sortie de l'exécution (par exemple `work/17/428668...`) et vérifiez le fichier `.command.sh` à l'intérieur :

```bash
cat work/17/428668/.command.sh
```

??? success "Sortie de la commande"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Vous devriez voir `-b 5` dans la commande `seqtk trimfq`.

Une chose importante à savoir sur `ext.args` : si un module a déjà une valeur par défaut définie, votre valeur la **remplacera complètement** plutôt que de s'y ajouter.
Par exemple, `FASTQC` a `ext.args = '--quiet'` défini par défaut dans `conf/modules.config` :

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Si vous définissez `ext.args = '--kmers 8'` pour `FASTQC`, l'indicateur `--quiet` ne sera plus appliqué.
Pour conserver les deux, définissez `ext.args = '--quiet --kmers 8'`.

Vous devriez toujours vérifier la configuration par défaut d'un module avant de remplacer `ext.args`.

### À retenir

Vous savez où se trouvent les valeurs par défaut de configuration des pipelines nf-core, et comment remplacer les allocations de ressources et les arguments des outils avec un fichier de configuration personnalisé.

### Et ensuite ?

Passez à la [Partie 3](./03_run_production_pipeline.md), où vous appliquerez ce que vous avez appris à un vrai pipeline de production.

---

## Résumé

Dans cette partie, vous avez appris à :

- Obtenir de l'aide, définir des paramètres et comprendre la validation des paramètres et des entrées
- Personnaliser l'allocation des ressources et les arguments des outils via des fichiers de configuration
