# Partie 4 : Créer un module nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette quatrième partie du cours de formation Hello nf-core, nous vous montrons comment créer un module nf-core en appliquant les conventions clés qui rendent les modules portables et maintenables.

Le projet nf-core fournit une commande (`nf-core modules create`) qui génère automatiquement des modèles de modules correctement structurés, similaire à ce que nous avons utilisé pour le workflow dans la Partie 2.
Cependant, à des fins pédagogiques, nous allons commencer par le faire manuellement : transformer le module local `cowpy` dans votre pipeline `core-hello` en un module de style nf-core étape par étape.
Après cela, nous vous montrerons comment utiliser la création de modules basée sur des modèles pour travailler plus efficacement à l'avenir.

??? info "Comment commencer à partir de cette section"

    Cette section suppose que vous avez terminé la [Partie 3 : Utiliser un module nf-core](./03_use_module.md) et que vous avez intégré le module `CAT_CAT` dans votre pipeline.

    Si vous n'avez pas terminé la Partie 3 ou souhaitez repartir de zéro pour cette partie, vous pouvez utiliser la solution `core-hello-part3` comme point de départ.
    Exécutez ces commandes depuis l'intérieur du répertoire `hello-nf-core/` :

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Cela vous donne un pipeline avec le module `CAT_CAT` déjà intégré.
    Vous pouvez tester qu'il s'exécute avec succès en exécutant la commande suivante :

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Transformer `cowpy` en un module nf-core

Dans cette section, nous appliquerons les conventions nf-core au module local `cowpy` dans votre pipeline `core-hello`, le transformant en un module qui suit les normes de la communauté nf-core.

Voici le code actuel du module de processus `cowpy` :

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Nous appliquerons les conventions nf-core suivantes de manière incrémentale :

1. **Mettre le nom du processus en majuscules `COWPY`** pour suivre la convention.
2. **Mettre à jour `COWPY` pour utiliser des tuples de métadonnées** afin de propager les métadonnées d'échantillon à travers le workflow.
3. **Centraliser la configuration des arguments de l'outil avec `ext.args`** pour augmenter la polyvalence du module tout en gardant l'interface minimale.
4. **Standardiser la dénomination des sorties avec `ext.prefix`** pour promouvoir la cohérence.
5. **Centraliser la configuration de publication** pour promouvoir la cohérence.

Après chaque étape, nous exécuterons le pipeline pour tester que tout fonctionne comme prévu.

!!! warning "Répertoire de travail"

    Assurez-vous d'être dans le répertoire `core-hello` (la racine de votre pipeline) pour toutes les modifications de fichiers et l'exécution de commandes dans cette section.

    ```bash
    cd core-hello
    ```

### 1.1. Mettre le nom du processus en majuscules

Il s'agit purement d'une convention stylistique (il n'y a pas de justification technique) mais comme c'est la norme pour les modules nf-core, conformons-nous.

Nous devons effectuer trois séries de modifications :

1. Mettre à jour le nom du processus dans le module
2. Mettre à jour l'instruction d'importation du module dans l'en-tête du workflow
3. Mettre à jour l'appel du processus et la déclaration emit dans le corps du workflow

Commençons !

#### 1.1.1. Mettre à jour le nom du processus dans le module

Ouvrez le fichier du module `cowpy.nf` (sous `core-hello/modules/local/`) et modifiez le nom du processus en majuscules :

=== "Après"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Avant"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

Dans ce cas, la mise en majuscules est complètement directe.

Si le nom du processus était composé de plusieurs mots, par exemple si nous avions un processus appelé MyCowpyTool à l'origine en camel case, la convention nf-core serait d'utiliser des underscores pour les séparer, donnant MY_COWPY_TOOL.

#### 1.1.2. Mettre à jour l'instruction d'importation du module

Les noms de processus sont sensibles à la casse, donc maintenant que nous avons changé le nom du processus, nous devons mettre à jour l'instruction d'importation du module en conséquence dans l'en-tête du workflow de `hello.nf` :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Nous pourrions utiliser un alias dans l'instruction d'importation pour éviter d'avoir à mettre à jour les appels au processus, mais cela irait quelque peu à l'encontre de l'objectif d'adopter la convention de mise en majuscules.

#### 1.1.3. Mettre à jour l'appel du processus et la déclaration emit

Maintenant, mettons à jour les deux références au processus dans le bloc workflow de `hello.nf` :

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // générer de l'art ASCII des salutations avec cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Rassembler et enregistrer les versions des logiciels
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // générer de l'art ASCII des salutations avec cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Rassembler et enregistrer les versions des logiciels
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Assurez-vous de faire les **deux** modifications, sinon vous obtiendrez une erreur lorsque vous exécuterez ceci.

#### 1.1.4. Exécuter le pipeline pour le tester

Exécutons le workflow pour tester que tout fonctionne correctement après ces modifications.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Très bien, cela fonctionne ! Maintenant, passons à des modifications plus substantielles.

### 1.2. Mettre à jour `COWPY` pour utiliser des tuples de métadonnées

Dans la version actuelle du pipeline `core-hello`, nous extrayons le fichier du tuple de sortie de `CAT_CAT` pour le passer à `COWPY`, comme indiqué dans la moitié supérieure du diagramme ci-dessous.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Il serait préférable que `COWPY` accepte directement des tuples de métadonnées, permettant aux métadonnées de circuler à travers le workflow, comme indiqué dans la moitié inférieure du diagramme.

À cette fin, nous devrons effectuer les modifications suivantes :

1. Mettre à jour les définitions d'entrée et de sortie
2. Mettre à jour l'appel du processus dans le workflow
3. Mettre à jour le bloc emit dans le workflow

Une fois que nous aurons fait tout cela, nous exécuterons le pipeline pour tester que tout fonctionne encore comme avant.

#### 1.2.1. Mettre à jour les définitions d'entrée et de sortie

Revenez au fichier du module `cowpy.nf` et modifiez-le pour accepter des tuples de métadonnées comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Avant"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Comme vous pouvez le voir, nous avons modifié à la fois l'**entrée principale** et la **sortie** en un tuple qui suit le modèle `tuple val(meta), path(input_file)` introduit dans la Partie 3 de cette formation.
Pour la sortie, nous en avons également profité pour ajouter `emit: cowpy_output` afin de donner un nom descriptif au canal de sortie.

Maintenant que nous avons modifié ce que le processus attend, nous devons mettre à jour ce que nous lui fournissons dans l'appel du processus.

#### 1.2.2. Mettre à jour l'appel du processus dans le workflow

La bonne nouvelle est que ce changement simplifiera l'appel du processus.
Maintenant que la sortie de `CAT_CAT` et l'entrée de `COWPY` ont la même 'forme', c'est-à-dire qu'elles consistent toutes deux en une structure `tuple val(meta), path(input_file)`, nous pouvons simplement les connecter directement au lieu d'avoir à extraire explicitement le fichier de la sortie du processus `CAT_CAT`.

Ouvrez le fichier workflow `hello.nf` (sous `core-hello/workflows/`) et mettez à jour l'appel à `COWPY` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // générer de l'art ASCII des salutations avec cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extraire le fichier du tuple puisque cowpy n'utilise pas encore les métadonnées
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // générer de l'art ASCII des salutations avec cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Nous appelons maintenant `COWPY` directement sur `CAT_CAT.out.file_out`.

Par conséquent, nous n'avons plus besoin de construire le canal `ch_for_cowpy`, donc cette ligne (et sa ligne de commentaire) peuvent être entièrement supprimées.

#### 1.2.3. Mettre à jour le bloc emit dans le workflow

Puisque `COWPY` émet maintenant une sortie nommée, `cowpy_output`, nous pouvons mettre à jour le bloc `emit:` du workflow `hello.nf` pour l'utiliser.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Ce n'est techniquement pas nécessaire, mais c'est une bonne pratique de faire référence aux sorties nommées chaque fois que possible.

#### 1.2.4. Exécuter le pipeline pour le tester

Exécutons le workflow pour tester que tout fonctionne correctement après ces modifications.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Le pipeline devrait s'exécuter avec succès, avec les métadonnées circulant maintenant de `CAT_CAT` à travers `COWPY`.

Cela complète ce que nous devions faire pour que `COWPY` gère les tuples de métadonnées.
Maintenant, voyons ce que nous pouvons faire d'autre pour tirer parti des modèles de modules nf-core.

### 1.3. Centraliser la configuration des arguments de l'outil avec `ext.args`

Dans son état actuel, le processus `COWPY` s'attend à recevoir une valeur pour le paramètre `character`.
Par conséquent, nous devons fournir une valeur à chaque fois que nous appelons le processus, même si nous serions satisfaits des valeurs par défaut définies par l'outil.
Pour `COWPY`, ce n'est certes pas un gros problème, mais pour des outils avec de nombreux paramètres optionnels, cela peut devenir assez lourd.

Le projet nf-core recommande d'utiliser une fonctionnalité de Nextflow appelée [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) pour gérer les arguments de l'outil plus commodément via des fichiers de configuration.

Au lieu de déclarer des entrées de processus pour chaque option de l'outil, vous écrivez le module pour référencer `ext.args` dans la construction de sa ligne de commande.
Ensuite, il suffit de configurer la variable `ext.args` pour contenir les arguments et les valeurs que vous souhaitez utiliser dans le fichier `modules.config`, qui consolide les détails de configuration pour tous les modules.
Nextflow ajoutera ces arguments avec leurs valeurs dans la ligne de commande de l'outil au moment de l'exécution.

Appliquons cette approche au module `COWPY`.
Nous allons devoir effectuer les modifications suivantes :

1. Mettre à jour le module `COWPY`
2. Configurer `ext.args` dans le fichier `modules.config`
3. Mettre à jour le workflow `hello.nf`

Une fois que nous aurons fait tout cela, nous exécuterons le pipeline pour tester que tout fonctionne encore comme avant.

#### 1.3.1. Mettre à jour le module `COWPY`

Faisons-le.
Ouvrez le fichier du module `cowpy.nf` (sous `core-hello/modules/local/`) et modifiez-le pour référencer `ext.args` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Avant"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Vous pouvez voir que nous avons fait trois modifications.

1. **Dans le bloc `input:`, nous avons supprimé l'entrée `val character`.**
   Dorénavant, nous fournirons cet argument via la configuration `ext.args` comme décrit ci-dessous.

2. **Dans le bloc `script:`, nous avons ajouté la ligne `def args = task.ext.args ?: ''`.**
   Cette ligne utilise l'opérateur `?:` pour déterminer la valeur de la variable `args` : le contenu de `task.ext.args` s'il n'est pas vide, ou une chaîne vide s'il l'est.
   Notez que bien que nous fassions généralement référence à `ext.args`, ce code doit référencer `task.ext.args` pour extraire la configuration `ext.args` au niveau du module.

3. **Dans la ligne de commande, nous avons remplacé `-c "$character"` par `$args`.**
   C'est ici que Nextflow injectera tous les arguments de l'outil définis dans `ext.args` dans le fichier `modules.config`.

Par conséquent, l'interface du module est maintenant plus simple : elle n'attend que les entrées essentielles de métadonnées et de fichiers.

!!! note

    L'opérateur `?:` est souvent appelé 'opérateur Elvis' car il ressemble à un visage d'Elvis Presley de côté, avec le caractère `?` symbolisant la vague dans ses cheveux.

#### 1.3.2. Configurer `ext.args` dans le fichier `modules.config`

Maintenant que nous avons retiré la déclaration `character` du module, nous devons l'ajouter à `ext.args` dans le fichier de configuration `modules.config`.

Plus précisément, nous allons ajouter ce petit morceau de code au bloc `process {}` :

```groovy title="Code à ajouter"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

La syntaxe `withName:` assigne cette configuration uniquement au processus `COWPY`, et `ext.args = { "-c ${params.character}" }` compose simplement une chaîne qui inclura la valeur du paramètre `character`.
Notez l'utilisation des accolades, qui indiquent à Nextflow d'évaluer la valeur du paramètre au moment de l'exécution.

C'est clair ? Ajoutons-le.

Ouvrez `conf/modules.config` et ajoutez le code de configuration à l'intérieur du bloc `process {}` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Avant"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Vous pouvez probablement imaginer avoir tous les modules d'un pipeline avec leur `ext.args` spécifié dans ce fichier, avec les avantages suivants :

- L'**interface du module reste simple** - Elle n'accepte que les entrées essentielles de métadonnées et de fichiers
- Le **pipeline expose toujours `params.character`** - Les utilisateur·trices final·es peuvent toujours le configurer comme avant
- Le **module est maintenant portable** - Il peut être réutilisé dans d'autres pipelines sans attendre un nom de paramètre spécifique
- La configuration est **centralisée** dans `modules.config`, gardant la logique du workflow propre

En utilisant le fichier `modules.config` comme l'endroit où tous les pipelines centralisent la configuration par module, nous rendons nos modules plus réutilisables dans différents pipelines.

#### 1.3.3. Mettre à jour le workflow `hello.nf`

Puisque le module `COWPY` ne nécessite plus le paramètre `character` comme entrée, nous devons mettre à jour l'appel du workflow en conséquence.

Ouvrez le fichier workflow `hello.nf` (sous `core-hello/workflows/`) et mettez à jour l'appel à `COWPY` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // générer de l'art ASCII des salutations avec cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // générer de l'art ASCII des salutations avec cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

Le code du workflow est maintenant plus propre : nous n'avons pas besoin de passer `params.character` directement au processus.
L'interface du module est maintenue minimale, la rendant plus portable, tandis que le pipeline fournit toujours l'option explicite via la configuration.

#### 1.3.4. Exécuter le pipeline pour le tester

Testons que le workflow fonctionne toujours comme prévu, en spécifiant un caractère différent pour vérifier que la configuration `ext.args` fonctionne.

Exécutez cette commande en utilisant `kosh`, l'une des options les plus... énigmatiques :

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Cela devrait s'exécuter avec succès comme précédemment.

Vérifions que la configuration `ext.args` a fonctionné en vérifiant la sortie.
Trouvez la sortie dans le navigateur de fichiers ou utilisez le hash de tâche (la partie `38/eb29ea` dans l'exemple ci-dessus) pour regarder le fichier de sortie :

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Sortie de la commande"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

Vous devriez voir l'art ASCII affiché avec le caractère `kosh`, confirmant que la configuration `ext.args` a fonctionné !

??? info "(Optionnel) Inspecter le fichier de commande"

    Si vous voulez voir exactement comment la configuration a été appliquée, vous pouvez inspecter le fichier `.command.sh` :

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Vous verrez la commande `cowpy` avec l'argument `-c kosh` :

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Cela montre que le fichier `.command.sh` a été généré correctement en fonction de la configuration `ext.args`.

Prenez un moment pour réfléchir à ce que nous avons accompli ici.
Cette approche maintient l'interface du module focalisée sur les données essentielles (fichiers, métadonnées et tous paramètres obligatoires par échantillon), tandis que les options qui contrôlent le comportement de l'outil sont gérées séparément via la configuration.

Cela peut sembler inutile pour un outil simple comme `cowpy`, mais cela peut faire une grande différence pour les outils d'analyse de données qui ont beaucoup d'arguments optionnels.

Pour résumer les avantages de cette approche :

- **Interface propre** : Le module se concentre sur les entrées de données essentielles (métadonnées et fichiers)
- **Flexibilité** : Les utilisateur·trices peuvent spécifier les arguments de l'outil via la configuration, y compris des valeurs spécifiques aux échantillons
- **Cohérence** : Tous les modules nf-core suivent ce modèle
- **Portabilité** : Les modules peuvent être réutilisés sans options d'outils codées en dur
- **Pas de changements de workflow** : L'ajout ou la modification d'options d'outils ne nécessite pas de mise à jour du code du workflow

!!! note

    Le système `ext.args` a des capacités supplémentaires puissantes non couvertes ici, y compris le changement dynamique des valeurs d'arguments en fonction des métadonnées. Consultez les [spécifications des modules nf-core](https://nf-co.re/docs/guidelines/components/modules) pour plus de détails.

### 1.4. Standardiser la dénomination des sorties avec `ext.prefix`

Maintenant que nous avons donné au processus `COWPY` accès au metamap, nous pouvons commencer à profiter d'un autre modèle utile de nf-core : nommer les fichiers de sortie en fonction des métadonnées.

Ici, nous allons utiliser une fonctionnalité de Nextflow appelée `ext.prefix` qui nous permettra de standardiser la dénomination des fichiers de sortie dans tous les modules en utilisant `meta.id` (l'identifiant inclus dans le metamap), tout en étant capable de configurer les modules individuellement si désiré.

Ce sera similaire à ce que nous avons fait avec `ext.args`, avec quelques différences que nous détaillerons au fur et à mesure.

Appliquons cette approche au module `COWPY`.
Nous allons devoir effectuer les modifications suivantes :

1. Mettre à jour le module `COWPY`
2. Configurer `ext.prefix` dans le fichier `modules.config`

(Aucune modification nécessaire du workflow.)

Une fois que nous aurons fait cela, nous exécuterons le pipeline pour tester que tout fonctionne encore comme avant.

#### 1.4.1. Mettre à jour le module `COWPY`

Ouvrez le fichier du module `cowpy.nf` (sous `core-hello/modules/local/`) et modifiez-le pour référencer `ext.prefix` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Avant"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Vous pouvez voir que nous avons fait trois modifications.

1. **Dans le bloc `script:`, nous avons ajouté la ligne `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Cette ligne utilise l'opérateur `?:` pour déterminer la valeur de la variable `prefix` : le contenu de `task.ext.prefix` s'il n'est pas vide, ou l'identifiant du metamap (`meta.id`) s'il l'est.
   Notez que bien que nous fassions généralement référence à `ext.prefix`, ce code doit référencer `task.ext.prefix` pour extraire la configuration `ext.prefix` au niveau du module.

2. **Dans la ligne de commande, nous avons remplacé `cowpy-${input_file}` par `${prefix}.txt`.**
   C'est ici que Nextflow injectera la valeur de `prefix` déterminée par la ligne ci-dessus.

3. **Dans le bloc `output:`, nous avons remplacé `path("cowpy-${input_file}")` par `path("${prefix}.txt")`.**
   Cela réitère simplement quel sera le chemin du fichier selon ce qui est écrit dans la ligne de commande.

Par conséquent, le nom du fichier de sortie est maintenant construit en utilisant une valeur par défaut sensée (l'identifiant du metamap) combinée avec l'extension de format de fichier appropriée.

#### 1.4.2. Configurer `ext.prefix` dans le fichier `modules.config`

Dans ce cas, la valeur par défaut sensée n'est pas suffisamment expressive à notre goût ; nous voulons utiliser un modèle de dénomination personnalisé qui inclut le nom de l'outil, `cowpy-<id>.txt`, comme nous l'avions auparavant.

Nous ferons cela en configurant `ext.prefix` dans `modules.config`, tout comme nous l'avons fait pour le paramètre `character` avec `ext.args`, sauf que cette fois le bloc `withName: 'COWPY' {}` existe déjà, et nous devons juste ajouter la ligne suivante :

```groovy title="Code à ajouter"
ext.prefix = { "cowpy-${meta.id}" }
```

Cela composera la chaîne que nous voulons.
Notez qu'une fois de plus nous utilisons des accolades, cette fois pour dire à Nextflow d'évaluer la valeur de `meta.id` au moment de l'exécution.

Ajoutons-le.

Ouvrez `conf/modules.config` et ajoutez le code de configuration à l'intérieur du bloc `process {}` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Avant"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Au cas où vous vous poseriez la question, la closure `ext.prefix` a accès à la bonne pièce de métadonnées car la configuration est évaluée dans le contexte de l'exécution du processus, où les métadonnées sont disponibles.

#### 1.4.3. Exécuter le pipeline pour le tester

Testons que le workflow fonctionne toujours comme prévu.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Regardez la sortie dans le répertoire des résultats.
Vous devriez voir le fichier de sortie cowpy avec la même dénomination qu'avant : `cowpy-test.txt`, basée sur le nom de lot par défaut.

??? abstract "Contenu du répertoire"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

N'hésitez pas à modifier la configuration `ext.prefix` dans `conf/modules.config` pour vous assurer que vous pouvez changer le modèle de dénomination sans avoir à apporter de modifications au code du module ou du workflow.

Alternativement, vous pouvez également essayer d'exécuter ceci à nouveau avec un paramètre `--batch` différent spécifié sur la ligne de commande pour vous assurer que cette partie est toujours personnalisable à la volée.

Cela démontre comment `ext.prefix` vous permet de maintenir votre convention de dénomination préférée tout en gardant l'interface du module flexible.

Pour résumer les avantages de cette approche :

- **Dénomination standardisée** : Les fichiers de sortie sont généralement nommés en utilisant les identifiants d'échantillon des métadonnées
- **Configurable** : Les utilisateur·trices peuvent remplacer la dénomination par défaut si nécessaire
- **Cohérent** : Tous les modules nf-core suivent ce modèle
- **Prévisible** : Facile de savoir comment les fichiers de sortie seront appelés

Plutôt bien, non ?
Eh bien, il y a encore une modification importante que nous devons faire pour améliorer notre module pour qu'il corresponde aux directives nf-core.

### 1.5. Centraliser la configuration de publication

Vous avez peut-être remarqué que nous avons publié des sorties dans deux répertoires différents :

- **`results`** — Le répertoire de sortie original que nous utilisons depuis le début pour nos modules locaux, défini individuellement à l'aide de directives `publishDir` par module ;
- **`core-hello-results`** — Le répertoire de sortie défini avec `--outdir` sur la ligne de commande, qui a reçu les journaux nf-core et les résultats publiés par `CAT_CAT`.

C'est désordonné et sous-optimal ; il serait préférable d'avoir un seul emplacement pour tout.
Bien sûr, nous pourrions aller dans chacun de nos modules locaux et mettre à jour la directive `publishDir` manuellement pour utiliser le répertoire `core-hello-results`, mais qu'en est-il la prochaine fois que nous décidons de changer le répertoire de sortie ?

Avoir des modules individuels prendre des décisions de publication n'est clairement pas la voie à suivre, surtout dans un monde où le même module pourrait être utilisé dans beaucoup de pipelines différents, par des personnes qui ont des besoins ou des préférences différents.
Nous voulons pouvoir contrôler où les sorties sont publiées au niveau de la configuration du workflow.

"Hé," pourriez-vous dire, "`CAT_CAT` envoie ses sorties au `--outdir`. Peut-être devrions-nous copier sa directive `publishDir` ?"

Oui, c'est une excellente idée.

Sauf qu'il n'a pas de directive `publishDir`. (Allez-y, regardez le code du module.)

C'est parce que les pipelines nf-core centralisent le contrôle au niveau du workflow en configurant `publishDir` dans `conf/modules.config` plutôt que dans les modules individuels.
Plus précisément, le modèle nf-core déclare une directive `publishDir` par défaut (avec une structure de répertoire prédéfinie) qui s'applique à tous les modules sauf si une directive de remplacement est fournie.

Cela ne semble-t-il pas génial ? Se pourrait-il que pour profiter de cette directive par défaut, tout ce que nous ayons à faire soit de supprimer la directive `publishDir` actuelle de nos modules locaux ?

Essayons cela sur `COWPY` pour voir ce qui se passe, puis nous examinerons le code de la configuration par défaut pour comprendre comment cela fonctionne.

Enfin, nous démontrerons comment remplacer le comportement par défaut si désiré.

#### 1.5.1. Supprimer la directive `publishDir` de `COWPY`

Faisons-le.
Ouvrez le fichier du module `cowpy.nf` (sous `core-hello/modules/local/`) et supprimez la directive `publishDir` comme indiqué ci-dessous.

=== "Après"

    ```groovy title="core-hello/modules/local/cowpy.nf (extrait)" linenums="1"
    #!/usr/bin/env nextflow

    // Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Avant"

    ```groovy title="core-hello/modules/local/cowpy.nf (extrait)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Générer de l'art ASCII avec cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

C'est tout !

#### 1.5.2. Exécuter le pipeline pour le tester

Voyons ce qui se passe si nous exécutons le pipeline maintenant.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Jetez un œil à votre répertoire de travail actuel.
Maintenant, `core-hello-results` contient également les sorties du module `COWPY`.

??? abstract "Contenu du répertoire"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Vous pouvez voir que Nextflow a créé cette hiérarchie de répertoires basée sur les noms du workflow et du module.

Le code responsable se trouve dans le fichier `conf/modules.config`.
Voici la configuration `publishDir` par défaut qui fait partie du modèle nf-core et s'applique à tous les processus :

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Cela peut paraître compliqué, alors regardons chacun des trois composants :

- **`path:`** Détermine le répertoire de sortie en fonction du nom du processus.
  Le nom complet d'un processus contenu dans `task.process` inclut la hiérarchie des importations de workflow et de modules (comme `CORE_HELLO:HELLO:CAT_CAT`).
  Les opérations `tokenize` suppriment cette hiérarchie pour obtenir juste le nom du processus, puis prennent la première partie avant tout underscore (si applicable), et la convertissent en minuscules.
  C'est ce qui détermine que les résultats de `CAT_CAT` sont publiés dans `${params.outdir}/cat/`.
- **`mode:`** Contrôle comment les fichiers sont publiés (copie, lien symbolique, etc.).
  Ceci est configurable via le paramètre `params.publish_dir_mode`.
- **`saveAs:`** Filtre quels fichiers publier.
  Cet exemple exclut les fichiers `versions.yml` en renvoyant `null` pour eux, les empêchant d'être publiés.

Cela fournit une logique cohérente pour organiser les sorties.

La sortie est encore meilleure lorsque tous les modules d'un pipeline adoptent cette convention, alors n'hésitez pas à supprimer les directives `publishDir` des autres modules de votre pipeline.
Cette valeur par défaut sera appliquée même aux modules que nous n'avons pas explicitement modifiés pour suivre les directives nf-core.

Cela dit, vous pouvez décider que vous voulez organiser vos entrées différemment, et la bonne nouvelle est qu'il est facile de le faire.

#### 1.5.3. Remplacer la valeur par défaut

Pour remplacer la directive `publishDir` par défaut, vous pouvez simplement ajouter vos propres directives au fichier `conf/modules.config`.

Par exemple, vous pourriez remplacer la valeur par défaut pour un seul processus en utilisant le sélecteur `withName:`, comme dans cet exemple où nous ajoutons une directive `publishDir` personnalisée pour le processus 'COWPY'.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

Nous n'allons pas réellement faire ce changement, mais n'hésitez pas à jouer avec cela et voir quelle logique vous pouvez implémenter.

Le point est que ce système permet donne le meilleur des deux mondes : cohérence par défaut et flexibilité pour personnaliser la configuration à la demande.

Pour résumer, vous obtenez :

- **Source unique de vérité** : Toute la configuration de publication se trouve dans `modules.config`
- **Valeur par défaut utile** : Les processus fonctionnent prêts à l'emploi sans configuration par module
- **Personnalisation facile** : Remplacez le comportement de publication dans la configuration, pas dans le code du module
- **Modules portables** : Les modules ne codent pas en dur les emplacements de sortie

Cela complète l'ensemble des fonctionnalités de modules nf-core que vous devriez absolument apprendre à utiliser, mais il en existe d'autres sur lesquelles vous pouvez lire dans les [spécifications des modules nf-core](https://nf-co.re/docs/guidelines/components/modules).

### À retenir

Vous savez maintenant comment adapter les modules locaux pour suivre les conventions nf-core :

- Concevez vos modules pour accepter et propager des tuples de métadonnées ;
- Utilisez `ext.args` pour garder les interfaces de modules minimales et portables ;
- Utilisez `ext.prefix` pour une dénomination de fichiers de sortie standardisée et configurable ;
- Adoptez la directive `publishDir` centralisée par défaut pour une structure de répertoire de résultats cohérente.

### Et ensuite ?

Apprenez comment utiliser les outils intégrés de nf-core basés sur des modèles pour créer des modules de manière simple.

---

## 2. Créer un module avec les outils nf-core

Maintenant que vous avez appris les modèles de modules nf-core en les appliquant manuellement, voyons comment vous créeriez des modules en pratique.

### 2.1. Générer un squelette de module à partir d'un modèle

Similaire à ce qui existe pour créer des pipelines, le projet nf-core fournit des outils pour générer des modules correctement structurés basés sur un modèle, avec tous ces modèles intégrés dès le départ.

#### 2.1.1. Exécuter la commande de création de module

La commande `nf-core modules create` génère un modèle de module qui suit déjà toutes les conventions que vous avez apprises.

Créons une nouvelle version du module `COWPY` avec un modèle minimal en exécutant cette commande :

```bash
nf-core modules create --empty-template COWPY
```

Le drapeau `--empty-template` crée un modèle de démarrage propre sans code supplémentaire, facilitant la visualisation de la structure essentielle.

La commande s'exécute de manière interactive, vous guidant à travers la configuration.
Elle recherche automatiquement les informations sur l'outil à partir de dépôts de paquets comme Bioconda et bio.tools pour pré-remplir les métadonnées.

Vous serez invité·e à fournir plusieurs options de configuration :

- **Informations sur l'auteur** : Votre nom d'utilisateur GitHub pour l'attribution
- **Étiquette de ressource** : Un ensemble prédéfini d'exigences de calcul.
  Le projet nf-core fournit des étiquettes standard comme `process_single` pour les outils légers et `process_high` pour les outils exigeants.
  Ces étiquettes aident à gérer l'allocation des ressources dans différents environnements d'exécution.
- **Exigence de métadonnées** : Si le module a besoin d'informations spécifiques aux échantillons via une map `meta` (généralement
