# Partie 2 : Réécrire Hello pour nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette deuxième partie du cours de formation Hello nf-core, nous vous montrons comment créer une version compatible nf-core du pipeline produit par le cours pour débutants [Hello Nextflow](../hello_nextflow/index.md).

Vous aurez remarqué dans la première section de la formation que les pipelines nf-core suivent une structure assez élaborée avec de nombreux fichiers accessoires.
Créer tout cela à partir de zéro serait très fastidieux, donc la communauté nf-core a développé des outils pour le faire à partir d'un modèle à la place, afin d'amorcer le processus.

Nous allons vous montrer comment utiliser ces outils pour créer une structure de pipeline, puis adapter le code du pipeline « régulier » existant sur la structure nf-core.

Si vous n'êtes pas familier avec le pipeline Hello ou si vous avez besoin d'un rappel, consultez [cette page d'information](../info/hello_pipeline.md).

---

## 1. Créer un nouveau projet de pipeline

Tout d'abord, nous créons la structure pour le nouveau pipeline.

!!! note "Note"

    Assurez-vous d'être dans le répertoire `hello-nf-core` dans votre terminal.

### 1.1. Exécuter l'outil de création de pipeline basé sur un modèle

Commençons par créer un nouveau pipeline avec la commande `nf-core pipelines create`.
Cela créera une nouvelle structure de pipeline en utilisant le modèle de base nf-core, personnalisé avec un nom de pipeline, une description et un auteur.

```bash
nf-core pipelines create
```

L'exécution de cette commande ouvrira une Interface Utilisateur Texte (TUI) pour la création du pipeline :

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Cette TUI vous demandera de fournir des informations de base sur votre pipeline et vous offrira un choix de fonctionnalités à inclure ou exclure dans votre structure de pipeline.

- Sur l'écran de bienvenue, cliquez sur **Let's go!**.
- Sur l'écran `Choose pipeline type`, cliquez sur **Custom**.
- Entrez les détails de votre pipeline comme suit (en remplaçant `< VOTRE NOM >` par votre propre nom), puis cliquez sur **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < VOTRE NOM >
```

- Sur l'écran Template features, réglez `Toggle all features` sur **off**, puis **activez** sélectivement les éléments suivants. Vérifiez vos sélections et cliquez sur **Continue**.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- Sur l'écran `Final details`, cliquez sur **Finish**. Attendez que le pipeline soit créé, puis cliquez sur **Continue**.
- Sur l'écran Create GitHub repository, cliquez sur **Finish without creating a repo**. Cela affichera des instructions pour créer un dépôt GitHub ultérieurement. Ignorez-les et cliquez sur **Close**.

Une fois que la TUI se ferme, vous devriez voir la sortie de console suivante.

??? success "Sortie de la commande"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

Il n'y a pas de confirmation explicite dans la sortie de console que la création du pipeline a réussi, mais vous devriez voir un nouveau répertoire appelé `core-hello`.

Visualisez le contenu du nouveau répertoire pour voir combien de travail vous vous êtes épargné·e en utilisant le modèle.

```bash
tree core-hello
```

??? abstract "Contenu du répertoire"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

Cela fait beaucoup de fichiers !

Vous devriez en reconnaître beaucoup comme étant les mêmes que nous avons rencontrés lorsque nous avons exploré la structure du pipeline `nf-core/demo`.
Mais ne vous inquiétez pas si vous vous sentez encore un peu perdu·e ; nous parcourrons ensemble les parties importantes au cours de cette formation.

!!! note "Note"

    Une différence importante par rapport au pipeline `nf-core/demo` que nous avons examiné dans la première partie de cette formation est qu'il n'y a pas de répertoire `modules`.
    C'est parce que nous n'avons pas choisi d'inclure les modules nf-core par défaut.

### 1.2. Tester que la structure est fonctionnelle

Croyez-le ou non, même si vous n'avez pas encore ajouté de modules pour effectuer un vrai travail, la structure du pipeline peut en fait être exécutée en utilisant le profil de test, de la même manière que nous avons exécuté le pipeline `nf-core/demo`.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Cela vous montre que tout le câblage de base est en place.
Alors où sont les sorties ? Y en a-t-il ?

En fait, un nouveau répertoire de résultats appelé `core-hello-results` a été créé contenant les rapports d'exécution standard :

```bash
tree core-hello-results
```

??? abstract "Contenu du répertoire"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Vous pouvez jeter un œil aux rapports pour voir ce qui a été exécuté, et la réponse est : rien du tout !

![rapport de chronologie d'exécution vide](./img/execution_timeline_empty.png)

Regardons ce qu'il y a réellement dans le code.

### 1.3. Examiner le workflow de substitution

Si vous regardez à l'intérieur du fichier `main.nf`, vous verrez qu'il importe un workflow appelé `HELLO` depuis `workflows/hello`.

Ceci est équivalent au workflow `workflows/demo.nf` que nous avons rencontré dans la Partie 1, et sert de workflow de substitution pour notre workflow d'intérêt, avec certaines fonctionnalités nf-core déjà en place.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

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
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Par rapport à un workflow Nextflow basique comme celui développé dans [Hello Nextflow](../hello_nextflow/index.md), vous remarquerez quelques nouveautés ici (lignes surlignées ci-dessus) :

- Le bloc workflow a un nom
- Les entrées du workflow sont déclarées en utilisant le mot-clé `take:` et la construction du canal est déplacée vers le workflow parent
- Le contenu du workflow est placé à l'intérieur d'un bloc `main:`
- Les sorties sont déclarées en utilisant le mot-clé `emit:`

Ce sont des fonctionnalités optionnelles de Nextflow qui rendent le workflow **composable**, ce qui signifie qu'il peut être appelé depuis un autre workflow.

!!! note "Workflows composables en profondeur"

    La [Side Quest Workflows of Workflows](../side_quests/workflows_of_workflows.md) explore la composition de workflows de manière beaucoup plus approfondie, y compris comment composer plusieurs workflows ensemble et gérer des flux de données complexes entre eux. Nous introduisons la composabilité ici car c'est une exigence fondamentale de l'architecture du modèle nf-core, qui utilise des workflows imbriqués pour organiser l'initialisation du pipeline, le workflow d'analyse principal et les tâches de finalisation en composants séparés et réutilisables.

Nous allons devoir intégrer la logique pertinente de notre workflow d'intérêt dans cette structure.
La première étape pour cela est de rendre notre workflow original composable.

### À retenir

Vous savez maintenant comment créer une structure de pipeline en utilisant les outils nf-core.

### Et ensuite ?

Apprenez à rendre un workflow simple composable comme prélude à le rendre compatible nf-core.

---

## 2. Rendre le workflow Hello Nextflow original composable

Maintenant il est temps de se mettre au travail pour intégrer notre workflow dans la structure nf-core.
Pour rappel, nous travaillons avec le workflow présenté dans notre cours de formation [Hello Nextflow](../hello_nextflow/index.md).

!!! tip "Astuce"

    Si vous n'êtes pas familier avec ce pipeline ou si vous avez besoin d'un rappel, consultez [Le pipeline Hello](../info/hello_pipeline.md).

Nous vous fournissons une copie propre et entièrement fonctionnelle du workflow Hello Nextflow terminé dans le répertoire `original-hello` avec ses modules et le fichier CSV par défaut qu'il s'attend à utiliser comme entrée.

```bash
tree original-hello/
```

??? abstract "Contenu du répertoire"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

N'hésitez pas à l'exécuter pour vous assurer qu'il fonctionne :

```bash
nextflow run original-hello/hello.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Ouvrons le fichier de workflow `hello.nf` pour inspecter le code, qui est montré en entier ci-dessous (sans compter les processus, qui sont dans les modules) :

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Paramètres du pipeline
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Inclure les modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // créer un canal pour les entrées depuis un fichier CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // émettre une salutation
  sayHello(greeting_ch)

  // convertir la salutation en majuscules
  convertToUpper(sayHello.out)

  // collecter toutes les salutations dans un seul fichier
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // générer de l'art ASCII des salutations avec cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Comme vous pouvez le voir, ce workflow a été écrit comme un workflow simple sans nom qui peut être exécuté seul.
Pour le rendre exécutable depuis un workflow parent comme l'exige le modèle nf-core, nous devons le rendre **composable**.

Parcourons les modifications nécessaires une par une.

### 2.1. Nommer le workflow

Tout d'abord, donnons un nom au workflow afin de pouvoir y faire référence depuis un workflow parent.

=== "Après"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Avant"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

Les mêmes conventions s'appliquent aux noms de workflow qu'aux noms de modules.

### 2.2. Remplacer la construction de canal par `take`

Maintenant, remplacez la construction de canal par une simple déclaration `take` déclarant les entrées attendues.

=== "Après"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // canal de salutations
        greeting_ch
    ```

=== "Avant"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // créer un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Cela laisse les détails de la façon dont les entrées sont fournies au workflow parent.

Pendant que nous y sommes, nous pouvons également commenter la ligne `params.greeting = 'greetings.csv'`

=== "Après"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Paramètres du pipeline
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Avant"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Paramètres du pipeline
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Note"

    Si vous avez l'extension du serveur de langage Nextflow installée, le vérificateur de syntaxe éclairera votre code avec des lignes ondulées rouges.
    C'est parce que si vous mettez une déclaration `take:`, vous devez aussi avoir un `main:`.

    Nous ajouterons cela à l'étape suivante.

### 2.3. Préfacer les opérations du workflow avec la déclaration `main`

Ensuite, ajoutez une déclaration `main` avant le reste des opérations appelées dans le corps du workflow.

=== "Après"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // émettre une salutation
        sayHello(greeting_ch)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Avant"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // émettre une salutation
        sayHello(greeting_ch)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Cela dit essentiellement « voici ce que ce workflow _fait_ ».

### 2.4. Ajouter la déclaration `emit`

Enfin, ajoutez une déclaration `emit` déclarant quelles sont les sorties finales du workflow.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Il s'agit d'un ajout entièrement nouveau au code par rapport au workflow original.

### 2.5. Récapitulatif des modifications complétées

Si vous avez effectué toutes les modifications comme décrit, votre workflow devrait maintenant ressembler à ceci :

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Paramètres du pipeline
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Inclure les modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // canal de salutations
    greeting_ch

    main:

    // émettre une salutation
    sayHello(greeting_ch)

    // convertir la salutation en majuscules
    convertToUpper(sayHello.out)

    // collecter toutes les salutations dans un seul fichier
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // générer de l'art ASCII des salutations avec cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Cela décrit tout ce dont Nextflow a besoin SAUF ce qu'il faut alimenter dans le canal d'entrée.
Cela va être défini dans le workflow parent, également appelé workflow **d'entrée**.

### 2.6. Créer un workflow d'entrée factice

Avant d'intégrer notre workflow composable dans la structure nf-core complexe, vérifions qu'il fonctionne correctement.
Nous pouvons créer un simple workflow d'entrée factice pour tester le workflow composable de manière isolée.

Créez un fichier vide nommé `main.nf` dans le même répertoire `original-hello`.

```bash
touch original-hello/main.nf
```

Copiez le code suivant dans le fichier `main.nf`.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// importer le code du workflow depuis le fichier hello.nf
include { HELLO } from './hello.nf'

// déclarer le paramètre d'entrée
params.greeting = 'greetings.csv'

workflow {
  // créer un canal pour les entrées depuis un fichier CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // appeler le workflow importé sur le canal de salutations
  HELLO(greeting_ch)

  // visualiser les sorties émises par le workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

Il y a deux observations importantes à faire ici :

- La syntaxe pour appeler le workflow importé est essentiellement la même que la syntaxe pour appeler des modules.
- Tout ce qui est lié à l'acheminement des entrées dans le workflow (paramètre d'entrée et construction de canal) est maintenant déclaré dans ce workflow parent.

!!! note "Note"

    Nommer le fichier de workflow d'entrée `main.nf` est une convention, pas une exigence.

    Si vous suivez cette convention, vous pouvez omettre de spécifier le nom du fichier de workflow dans votre commande `nextflow run`.
    Nextflow recherchera automatiquement un fichier nommé `main.nf` dans le répertoire d'exécution.

    Cependant, vous pouvez nommer le fichier de workflow d'entrée autrement si vous préférez.
    Dans ce cas, assurez-vous de spécifier le nom du fichier de workflow dans votre commande `nextflow run`.

### 2.7. Tester que le workflow s'exécute

Nous avons enfin toutes les pièces dont nous avons besoin pour vérifier que le workflow composable fonctionne.
Exécutons-le !

```bash
nextflow run ./original-hello
```

Ici vous voyez l'avantage d'utiliser la convention de nommage `main.nf`.
Si nous avions nommé le workflow d'entrée `something_else.nf`, nous aurions dû faire `nextflow run original-hello/something_else.nf`.

Si vous avez effectué toutes les modifications correctement, cela devrait s'exécuter jusqu'à la fin.

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Cela signifie que nous avons réussi à mettre à niveau notre workflow HELLO pour le rendre composable.

### À retenir

Vous savez comment rendre un workflow composable en lui donnant un nom et en ajoutant des déclarations `take`, `main` et `emit`, et comment l'appeler depuis un workflow d'entrée.

### Et ensuite ?

Apprenez à greffer un workflow composable basique sur la structure nf-core.

---

## 3. Adapter la logique du workflow mis à jour dans le workflow de substitution

Maintenant que nous avons vérifié que notre workflow composable fonctionne correctement, retournons à la structure du pipeline nf-core que nous avons créée dans la section 1.
Nous voulons intégrer le workflow composable que nous venons de développer dans la structure du modèle nf-core, donc le résultat final devrait ressembler à ceci.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Alors comment faisons-nous cela ? Jetons un œil au contenu actuel du workflow `HELLO` dans `core-hello/workflows/hello.nf` (la structure nf-core).

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

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
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Dans l'ensemble, ce code fait très peu de choses en dehors de quelques tâches administratives liées à la capture de la version de tous les outils logiciels qui sont exécutés dans le pipeline.

Nous devons ajouter le code pertinent de la version composable du workflow original que nous avons développée dans la section 2.

Nous allons aborder cela dans les étapes suivantes :

1. Copier les modules et configurer les imports de modules
2. Laisser la déclaration `take` telle quelle
3. Ajouter la logique du workflow au bloc `main`
4. Mettre à jour le bloc `emit`

!!! note "Note"

    Nous allons ignorer la capture de version pour cette première passe et verrons comment la câbler dans une partie ultérieure de cette formation.

### 3.1. Copier les modules et configurer les imports de modules

Les quatre processus de notre workflow Hello Nextflow sont stockés comme modules dans `original-hello/modules/`.
Nous devons copier ces modules dans la structure du projet nf-core (sous `core-hello/modules/local/`) et ajouter des déclarations d'importation au fichier de workflow nf-core.

Copions d'abord les fichiers de modules de `original-hello/` vers `core-hello/` :

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Vous devriez maintenant voir le répertoire des modules listé sous `core-hello/`.

```bash
tree core-hello/modules
```

??? abstract "Contenu du répertoire"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Maintenant configurons les déclarations d'importation de modules.

Voici les déclarations d'importation dans le workflow `original-hello/hello.nf` :

```groovy title="original-hello/hello.nf" linenums="9"
// Inclure les modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Ouvrez le fichier `core-hello/workflows/hello.nf` et transposez ces déclarations d'importation dedans comme montré ci-dessous.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Deux autres observations intéressantes ici :

- Nous avons adapté le formatage des déclarations d'importation pour suivre la convention de style nf-core.
- Nous avons mis à jour les chemins relatifs vers les modules pour refléter qu'ils sont maintenant stockés à un niveau d'imbrication différent.

### 3.2. Laisser la déclaration `take` telle quelle

Le projet nf-core a beaucoup de fonctionnalités préconçues autour du concept de samplesheet, qui est typiquement un fichier CSV contenant des données en colonnes.
Puisque c'est essentiellement ce qu'est notre fichier `greetings.csv`, nous garderons la déclaration `take` actuelle telle quelle, et mettrons simplement à jour le nom du canal d'entrée à l'étape suivante.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

La gestion des entrées sera effectuée en amont de ce workflow (pas dans ce fichier de code).

### 3.3. Ajouter la logique du workflow au bloc `main`

Maintenant que nos modules sont disponibles pour le workflow, nous pouvons intégrer la logique du workflow dans le bloc `main`.

Pour rappel, voici le code pertinent dans le workflow original, qui n'a pas beaucoup changé lorsque nous l'avons rendu composable (nous avons juste ajouté la ligne `main:`) :

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // émettre une salutation
    sayHello(greeting_ch)

    // convertir la salutation en majuscules
    convertToUpper(sayHello.out)

    // collecter toutes les salutations dans un seul fichier
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // générer de l'art ASCII des salutations avec cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

Nous devons copier le code qui vient après `main:` dans la nouvelle version du workflow.

Il y a déjà du code là-dedans qui a à voir avec la capture des versions des outils qui sont exécutés par le workflow. Nous allons laisser cela tranquille pour l'instant (nous traiterons les versions d'outils plus tard).
Nous garderons l'initialisation `ch_versions = channel.empty()` en haut, puis insérerons notre logique de workflow, en gardant le code de collecte des versions à la fin.
Cet ordonnancement a du sens car dans un vrai pipeline, les processus émettraient des informations de version qui seraient ajoutées au canal `ch_versions` pendant l'exécution du workflow.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // émettre une salutation
        sayHello(greeting_ch)

        // convertir la salutation en majuscules
        convertToUpper(sayHello.out)

        // collecter toutes les salutations dans un seul fichier
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // générer de l'art ASCII des salutations avec cowpy
        cowpy(collectGreetings.out.outfile, params.character)

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
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

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
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

Vous remarquerez que nous avons également ajouté une ligne vide avant `main:` pour rendre le code plus lisible.

Cela semble très bien, mais nous devons encore mettre à jour le nom du canal que nous passons au processus `sayHello()` de `greeting_ch` à `ch_samplesheet` comme montré ci-dessous, pour correspondre à ce qui est écrit sous le mot-clé `take:`.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // émettre une salutation (mis à jour pour utiliser la convention nf-core pour les samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // émettre une salutation
        sayHello(greeting_ch)
    ```

Maintenant la logique du workflow est correctement câblée.

### 3.4. Mettre à jour le bloc `emit`

Enfin, nous devons mettre à jour le bloc `emit` pour inclure la déclaration des sorties finales du workflow.

=== "Après"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Avant"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Ceci conclut les modifications que nous devons apporter au workflow HELLO lui-même.
À ce stade, nous avons atteint la structure globale du code que nous nous étions fixés pour objectif.

### À retenir

Vous savez comment adapter les éléments principaux d'un workflow composable dans un workflow de substitution nf-core.

### Et ensuite ?

Apprenez à adapter la gestion des entrées dans la structure du pipeline nf-core.

---

## 4. Adapter la gestion des entrées

Maintenant que nous avons intégré avec succès notre logique de workflow dans la structure nf-core, nous devons aborder une dernière pièce critique : nous assurer que nos données d'entrée sont traitées correctement.
Le modèle nf-core est livré avec une gestion des entrées sophistiquée conçue pour des ensembles de données génomiques complexes, donc nous devons l'adapter pour qu'il fonctionne avec notre fichier `greetings.csv` plus simple.

### 4.1. Identifier où les entrées sont gérées

La première étape est de déterminer où la gestion des entrées est effectuée.

Vous vous souvenez peut-être que lorsque nous avons réécrit le workflow Hello Nextflow pour le rendre composable, nous avons déplacé la déclaration du paramètre d'entrée d'un niveau vers le haut, dans le workflow d'entrée `main.nf`.
Jetons donc un œil au workflow d'entrée de niveau supérieur `main.nf` qui a été créé dans le cadre de la structure du pipeline :

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Le projet nf-core fait un usage intensif de sous-workflows imbriqués, donc cette partie peut être un peu déroutante à première approche.

Ce qui compte ici, c'est qu'il y a deux workflows définis :

- `CORE_HELLO` est une enveloppe mince pour exécuter le workflow HELLO que nous venons de finir d'adapter dans `core-hello/workflows/hello.nf`.
- Un workflow sans nom qui appelle `CORE_HELLO` ainsi que deux autres sous-workflows, `PIPELINE_INITIALISATION` et `PIPELINE_COMPLETION`.

Voici un diagramme de la façon dont ils sont liés les uns aux autres :

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Fait important, nous ne pouvons trouver aucun code construisant un canal d'entrée à ce niveau, seulement des références à un samplesheet fourni via le paramètre `--input`.

Un peu de fouille révèle que la gestion des entrées est effectuée par le sous-workflow `PIPELINE_INITIALISATION`, comme il se doit, qui est importé depuis `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`.

Si nous ouvrons ce fichier et faisons défiler vers le bas, nous arrivons à ce bloc de code :

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Crée un canal à partir du fichier d'entrée fourni via params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

C'est la fabrique de canaux qui analyse le samplesheet et le transmet sous une forme prête à être consommée par le workflow HELLO.

!!! note "Note"

    La syntaxe ci-dessus est un peu différente de ce que nous avons utilisé précédemment, mais fondamentalement ceci :

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    est équivalent à ceci :

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Ce code implique quelques étapes d'analyse et de validation qui sont très spécifiques à l'exemple de samplesheet inclus avec le modèle de pipeline nf-core, qui au moment de la rédaction est très spécifique au domaine et ne convient pas à notre simple projet de pipeline.

### 4.2. Remplacer le code de canal d'entrée du modèle

La bonne nouvelle est que les besoins de notre pipeline sont beaucoup plus simples, donc nous pouvons remplacer tout cela par le code de construction de canal que nous avons développé dans le workflow Hello Nextflow original.

Pour rappel, voici à quoi ressemblait la construction du canal (comme vu dans le répertoire solutions) :

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // créer un canal pour les entrées depuis un fichier CSV
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Donc nous devons juste intégrer cela dans le workflow d'initialisation, avec des modifications mineures : nous mettons à jour le nom du canal de `greeting_ch` à `ch_samplesheet`, et le nom du paramètre de `params.greeting` à `params.input` (voir ligne surlignée).

=== "Après"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Crée un canal à partir du fichier d'entrée fourni via params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Avant"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Crée un canal à partir du fichier d'entrée fourni via params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Cela termine les modifications que nous devons effectuer pour que le traitement des entrées fonctionne.

Dans sa forme actuelle, cela ne nous permettra pas de profiter des capacités intégrées de nf-core pour la validation de schéma, mais nous pouvons ajouter cela plus tard.
Pour l'instant, nous nous concentrons sur le garder aussi simple que possible pour arriver à quelque chose que nous pouvons exécuter avec succès sur des données de test.

### 4.3. Mettre à jour le profil de test

En parlant de données de test et de paramètres, mettons à jour le profil de test pour ce pipeline afin d'utiliser le mini-samplesheet `greetings.csv` au lieu de l'exemple de samplesheet fourni dans le modèle.

Sous `core-hello/conf`, nous trouvons deux profils de test basés sur le modèle : `test.config` et `test_full.config`, qui sont destinés à tester un petit échantillon de données et un échantillon de taille complète.
Étant donné l'objectif de notre pipeline, il n'y a pas vraiment d'intérêt à configurer un profil de test de taille complète, donc n'hésitez pas à ignorer ou supprimer `test_full.config`.
Nous allons nous concentrer sur la configuration de `test.config` pour s'exécuter sur notre fichier `greetings.csv` avec quelques paramètres par défaut.

#### 4.3.1. Copier le fichier `greetings.csv`

D'abord, nous devons copier le fichier `greetings.csv` dans un endroit approprié de notre projet de pipeline.
Généralement, les petits fichiers de test sont stockés dans le répertoire `assets`, donc copions le fichier depuis notre répertoire de travail.

```bash
cp greetings.csv core-hello/assets/.
```

Maintenant le fichier `greetings.csv` est prêt à être utilisé comme entrée de test.

#### 4.3.2. Mettre à jour le fichier `test.config`

Maintenant nous pouvons mettre à jour le fichier `test.config` comme suit :

=== "Après"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Données d'entrée
        input  = "${projectDir}/assets/greetings.csv"

        // Autres paramètres
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Avant"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Données d'entrée
        // TODO nf-core: Spécifiez les chemins vers vos données de test sur nf-core/test-datasets
        // TODO nf-core: Donnez tous les paramètres requis pour le test afin que les flags de ligne de commande ne soient pas nécessaires
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Points clés :

- **Utilisation de `${projectDir}`** : C'est une variable implicite de Nextflow qui pointe vers le répertoire où se trouve le script de workflow principal (la racine du pipeline). Son utilisation garantit que le chemin fonctionne quel que soit l'endroit d'où le pipeline est exécuté.
- **Chemins absolus** : En utilisant `${projectDir}`, nous créons un chemin absolu, ce qui est important pour les données de test qui sont fournies avec le pipeline.
- **Emplacement des données de test** : les pipelines nf-core stockent généralement les données de test dans le répertoire `assets/` au sein du dépôt du pipeline pour les petits fichiers de test, ou référencent des ensembles de données de test externes pour les fichiers plus volumineux.

Et pendant que nous y sommes, resserrons les limites de ressources par défaut pour garantir que cela s'exécutera sur des machines très basiques (comme les VM minimales dans Github Codespaces) :

=== "Après"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Avant"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Cela termine les modifications de code que nous devons effectuer.

### 4.4. Exécuter le pipeline avec le profil de test

C'était beaucoup, mais nous pouvons enfin essayer d'exécuter le pipeline !
Notez que nous devons ajouter `--validate_params false` à la ligne de commande car nous n'avons pas encore configuré la validation (cela viendra plus tard).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Si vous avez effectué toutes les modifications correctement, cela devrait s'exécuter jusqu'à la fin.

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Comme vous pouvez le voir, cela a produit le résumé typique nf-core au démarrage grâce au sous-workflow d'initialisation, et les lignes pour chaque module montrent maintenant les noms complets PIPELINE:WORKFLOW:module.

### 4.5. Trouver les sorties du pipeline

La question maintenant est : où sont les sorties du pipeline ?
Et la réponse est assez intéressante : il y a maintenant deux endroits différents où chercher les résultats.

Comme vous vous en souvenez peut-être d'avant, notre première exécution du workflow nouvellement créé a produit un répertoire appelé `core-hello-results/` qui contenait divers rapports d'exécution et métadonnées.

```bash
tree core-hello-results
```

??? abstract "Contenu du répertoire"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

Vous voyez que nous avons obtenu un autre ensemble de rapports d'exécution en plus de ceux que nous avons obtenus lors de la première exécution, lorsque le workflow n'était encore qu'un substitut.
Cette fois, vous voyez toutes les tâches qui ont été exécutées comme prévu.

![rapport de chronologie d'exécution pour le pipeline Hello](./img/execution_timeline_hello.png)

!!! note "Note"

    Une fois de plus, les tâches n'ont pas été exécutées en parallèle car nous exécutons sur une machine minimaliste dans Github Codespaces.
    Pour voir celles-ci s'exécuter en parallèle, essayez d'augmenter l'allocation de CPU de votre codespace et les limites de ressources dans la configuration de test.

C'est génial, mais nos résultats de pipeline réels ne sont pas là !

Voici ce qui s'est passé : nous n'avons rien changé aux modules eux-mêmes, donc les sorties gérées par les directives `publishDir` au niveau du module vont toujours dans un répertoire `results` tel que spécifié dans le pipeline original.

```bash
tree results
```

??? abstract "Contenu du répertoire"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, les voilà, mélangés avec les sorties des exécutions précédentes du pipeline Hello original.

Si nous voulons qu'ils soient bien organisés comme les sorties du pipeline demo l'étaient, nous devrons changer la façon dont nous configurons la publication des sorties.
Nous vous montrerons comment faire cela plus tard dans ce cours de formation.

<!-- TODO: Mettre à jour ceci une fois que nous aurons mis à jour Hello Nextflow pour utiliser des sorties au niveau workflow -->

Et voilà ! Cela peut sembler beaucoup de travail pour accomplir le même résultat que le pipeline original, mais vous obtenez tous ces beaux rapports générés automatiquement, et vous avez maintenant une base solide pour profiter de fonctionnalités supplémentaires de nf-core, y compris la validation des entrées et certaines capacités de gestion des métadonnées intéressantes que nous couvrirons dans une section ultérieure.

---

### À retenir

Vous savez comment convertir un pipeline Nextflow ordinaire en un pipeline de style nf-core en utilisant le modèle nf
