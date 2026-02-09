# Partie 2 : Exécuter de vrais pipelines

Dans la Partie 1 de ce cours (Exécuter des opérations de base), nous avons commencé avec un exemple de workflow qui ne comportait que des fonctionnalités minimales afin de maintenir la complexité du code à un niveau faible.
Par exemple, `1-hello.nf` utilisait un paramètre de ligne de commande (`--input`) pour fournir une seule valeur à la fois.

Cependant, la plupart des pipelines du monde réel utilisent des fonctionnalités plus sophistiquées afin de permettre un traitement efficace de grandes quantités de données à grande échelle, et appliquent plusieurs étapes de traitement enchaînées par une logique parfois complexe.

Dans cette partie de la formation, nous démontrons les fonctionnalités clés des pipelines du monde réel en essayant des versions étendues du pipeline Hello World original.

## 1. Traiter des données d'entrée à partir d'un fichier

Dans un pipeline du monde réel, nous voulons généralement traiter plusieurs points de données (ou séries de données) contenus dans un ou plusieurs fichiers d'entrée.
Et dans la mesure du possible, nous voulons exécuter le traitement de données indépendantes en parallèle, afin de réduire le temps d'attente pour l'analyse.

Pour démontrer comment Nextflow fait cela, nous avons préparé un fichier CSV appelé `greetings.csv` qui contient plusieurs messages d'accueil en entrée, imitant le type de données en colonnes que vous pourriez vouloir traiter dans une analyse de données réelle.
Notez que les nombres ne sont pas significatifs, ils sont là uniquement à des fins d'illustration.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Nous avons également écrit une version améliorée du workflow original, maintenant appelée `2a-inputs.nf`, qui lira le fichier CSV, extraira les messages d'accueil et écrira chacun d'eux dans un fichier séparé.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Exécutons d'abord le workflow, puis nous examinerons le code Nextflow pertinent par la suite.

### 1.1. Exécuter le workflow

Exécutez la commande suivante dans votre terminal.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

De manière encourageante, cela semble indiquer que '3 of 3' appels ont été effectués pour le processus, ce qui est encourageant, puisqu'il y avait trois lignes de données dans le CSV que nous avons fourni en entrée.
Cela suggère que le processus `sayHello()` a été appelé trois fois, une fois sur chaque ligne d'entrée.

### 1.2. Trouver les sorties publiées dans le répertoire `results`

Regardons le répertoire 'results' pour voir si notre workflow écrit toujours une copie de nos sorties là-bas.

??? abstract "Contenu du répertoire"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Oui ! Nous voyons un nouveau répertoire appelé `2a-inputs` avec trois fichiers de sortie avec des noms différents, ce qui est assez pratique.

Vous pouvez ouvrir chacun d'eux pour vous assurer qu'ils contiennent la chaîne de message d'accueil appropriée.

??? abstract "Contenu des fichiers"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Cela confirme que chaque message d'accueil dans le fichier d'entrée a été traité de manière appropriée.

### 1.3. Trouver les sorties originales et les journaux

Vous avez peut-être remarqué que la sortie console ci-dessus ne faisait référence qu'à un seul répertoire de tâche.
Cela signifie-t-il que les trois appels à `sayHello()` ont été exécutés dans ce seul répertoire de tâche ?

#### 1.3.1. Examiner le répertoire de tâche donné dans le terminal

Jetons un coup d'œil à l'intérieur de ce répertoire de tâche `8e/0eb066`.

??? abstract "Contenu du répertoire"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Nous ne trouvons que la sortie correspondant à l'un des messages d'accueil (ainsi que les fichiers accessoires si nous activons l'affichage des fichiers cachés).

Alors que se passe-t-il ici ?

Par défaut, le système de journalisation ANSI écrit les informations d'état pour tous les appels au même processus sur la même ligne.
En conséquence, il ne nous a montré qu'un seul des trois chemins de répertoire de tâche (`8e/0eb066`) dans la sortie console.
Il y en a deux autres qui ne sont pas listés là.

#### 1.3.2. Faire afficher plus de détails au terminal

Nous pouvons modifier le comportement de journalisation pour voir la liste complète des appels de processus en ajoutant `-ansi-log false` à la commande comme suit :

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Sortie de la commande"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Cette fois, nous voyons les trois exécutions de processus et leurs sous-répertoires de travail associés listés dans la sortie.
La désactivation de la journalisation ANSI a également empêché Nextflow d'utiliser des couleurs dans la sortie du terminal.

Notez que la façon dont l'état est rapporté est un peu différente entre les deux modes de journalisation.
En mode condensé, Nextflow indique si les appels ont été complétés avec succès ou non.
Dans ce mode étendu, il indique seulement qu'ils ont été soumis.

Cela confirme que le processus `sayHello()` est appelé trois fois, et qu'un répertoire de tâche séparé est créé pour chacun.

Si nous regardons à l'intérieur de chacun des répertoires de tâche listés là, nous pouvons vérifier que chacun correspond à l'un des messages d'accueil.

??? abstract "Contenu du répertoire"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Cela confirme que chaque appel de processus est exécuté de manière isolée de tous les autres.
Cela présente de nombreux avantages, notamment éviter les collisions si le processus produit des fichiers intermédiaires avec des noms non uniques.

!!! tip

    Pour un workflow complexe, ou un grand nombre d'entrées, avoir la liste complète affichée dans le terminal peut devenir un peu écrasant, donc les gens n'utilisent normalement pas `-ansi-log false` dans une utilisation courante.

### 1.4. Examiner le code du workflow

Donc cette version du workflow est capable de lire un fichier CSV d'entrées, de traiter les entrées séparément et de nommer les sorties de manière unique.

Jetons un coup d'œil à ce qui rend cela possible dans le code du workflow.

??? full-code "Fichier de code complet"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Encore une fois, vous n'avez pas besoin de mémoriser la syntaxe du code, mais il est bon d'apprendre à reconnaître les composants clés du workflow qui fournissent des fonctionnalités importantes.

#### 1.4.1. Charger les données d'entrée à partir du CSV

C'est la partie la plus intéressante : comment sommes-nous passés de la prise d'une seule valeur depuis la ligne de commande, à la prise d'un fichier CSV, son analyse et le traitement des messages d'accueil individuels qu'il contient ?

Dans Nextflow, nous faisons cela avec un [**canal**](https://nextflow.io/docs/latest/channel.html) : une construction de file d'attente conçue pour gérer les entrées efficacement et les transporter d'une étape à l'autre dans les workflows multi-étapes, tout en fournissant un parallélisme intégré et de nombreux avantages supplémentaires.

Décomposons cela.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

Ce code crée un canal appelé `greeting_ch` qui lit le fichier CSV, l'analyse et extrait la première colonne de chaque ligne.
Le résultat est un canal contenant `Hello`, `Bonjour` et `Holà`.

??? tip "Comment cela fonctionne-t-il ?"

    Voici ce que cette ligne signifie en français simple :

    - `channel.fromPath` est une **fabrique de canaux** qui crée un canal à partir de chemin(s) de fichier
    - `(params.input)` spécifie que le chemin de fichier est fourni par `--input` sur la ligne de commande

    En d'autres termes, cette ligne dit à Nextflow : prends le chemin de fichier donné avec `--input` et prépare-toi à traiter son contenu comme des données d'entrée.

    Ensuite, les deux lignes suivantes appliquent des **opérateurs** qui effectuent l'analyse réelle du fichier et le chargement des données dans la structure de données appropriée :

    - `.splitCsv()` dit à Nextflow d'analyser le fichier CSV en un tableau représentant les lignes et les colonnes
    - `.map { line -> line[0] }` dit à Nextflow de ne prendre que l'élément de la première colonne de chaque ligne

    Donc en pratique, en partant du fichier CSV suivant :

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Nous avons transformé cela en un tableau qui ressemble à ceci :

    ```txt title="Array contents"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    Et ensuite nous avons pris le premier élément de chacune des trois lignes et les avons chargés dans un canal Nextflow qui contient maintenant : `Hello`, `Bonjour` et `Holà`.

    Si vous voulez comprendre les canaux et les opérateurs en profondeur, y compris comment les écrire vous-même, consultez [Hello Nextflow Partie 2 : Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Appeler le processus sur chaque message d'accueil

Ensuite, dans la dernière ligne du bloc `main:` du workflow, nous fournissons le canal `greeting_ch` chargé comme entrée au processus `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
```

Cela dit à Nextflow d'exécuter le processus individuellement sur chaque élément du canal, _c'est-à-dire_ sur chaque message d'accueil.
Et parce que Nextflow est intelligent comme ça, il exécutera ces appels de processus en parallèle si possible, en fonction de l'infrastructure informatique disponible.

C'est ainsi que vous pouvez obtenir un traitement efficace et évolutif de beaucoup de données (de nombreux échantillons, ou points de données, quelle que soit votre unité de recherche) avec relativement très peu de code.

#### 1.4.3. Comment les sorties sont nommées

Enfin, il vaut la peine de jeter un coup d'œil rapide au code du processus pour voir comment nous obtenons que les fichiers de sortie soient nommés de manière unique.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

Vous voyez que, par rapport à la version de ce processus dans `1-hello.nf`, la déclaration de sortie et la partie pertinente de la commande ont changé pour inclure la valeur du message d'accueil dans le nom du fichier de sortie.

C'est une façon de s'assurer que les noms de fichiers de sortie n'entreront pas en collision lorsqu'ils seront publiés dans le répertoire de résultats commun.

Et c'est le seul changement que nous avons dû faire à l'intérieur de la déclaration du processus !

### À retenir

Vous comprenez à un niveau de base comment les canaux et les opérateurs nous permettent de traiter plusieurs entrées efficacement.

### Et ensuite ?

Découvrez comment les workflows multi-étapes sont construits et comment ils fonctionnent.

---

## 2. Exécuter des workflows multi-étapes

La plupart des workflows du monde réel impliquent plus d'une étape.
Construisons sur ce que nous venons d'apprendre sur les canaux, et regardons comment Nextflow utilise les canaux et les opérateurs pour connecter les processus ensemble dans un workflow multi-étapes.

À cette fin, nous vous fournissons un exemple de workflow qui enchaîne trois étapes séparées et démontre ce qui suit :

1. Faire circuler les données d'un processus au suivant
2. Collecter les sorties de plusieurs appels de processus en un seul appel de processus

Plus précisément, nous avons créé une version étendue du workflow appelée `2b-multistep.nf` qui prend chaque message d'accueil en entrée, le convertit en majuscules, puis collecte tous les messages d'accueil en majuscules dans un seul fichier de sortie.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Comme précédemment, nous allons d'abord exécuter le workflow puis examiner le code pour voir ce qui est nouveau.

### 2.1. Exécuter le workflow

Exécutez la commande suivante dans votre terminal :

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Sortie de la commande"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Vous voyez que comme promis, plusieurs étapes ont été exécutées dans le cadre du workflow ; les deux premières (`sayHello` et `convertToUpper`) ont vraisemblablement été exécutées sur chaque message d'accueil individuel, et la troisième (`collectGreetings`) aura été exécutée une seule fois, sur les sorties des trois appels `convertToUpper`.

### 2.2. Trouver les sorties

Vérifions que c'est bien ce qui s'est passé en regardant dans le répertoire `results`.

??? abstract "Contenu du répertoire"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Comme vous pouvez le voir, nous avons un nouveau répertoire appelé `2b-multistep`, et il contient beaucoup plus de fichiers qu'auparavant.
Certains des fichiers ont été regroupés dans un sous-répertoire appelé `intermediates`, tandis que deux fichiers sont situés au niveau supérieur.

Ces deux-là sont les résultats finaux du workflow multi-étapes.
Prenez une minute pour regarder les noms de fichiers et vérifier leur contenu pour confirmer qu'ils sont ce que vous attendez.

??? abstract "Contenu des fichiers"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

Le premier contient nos trois messages d'accueil, en majuscules et collectés dans un seul fichier comme promis.
Le second est un fichier de rapport qui résume certaines informations sur l'exécution.

### 2.3. Examiner le code

Regardons le code et identifions les modèles clés pour les workflows multi-étapes.

??? full-code "Fichier de code complet"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Il se passe beaucoup de choses là-dedans, mais la différence la plus évidente par rapport à la version précédente du workflow est qu'il y a maintenant plusieurs définitions de processus, et en conséquence, plusieurs appels de processus dans le bloc workflow.

Regardons de plus près et voyons si nous pouvons identifier les éléments les plus intéressants.

#### 2.3.1. Visualiser la structure du workflow

Si vous utilisez VSCode avec l'extension Nextflow, vous pouvez obtenir un diagramme utile de la façon dont les processus sont connectés en cliquant sur le petit lien `DAG preview` affiché juste au-dessus du bloc workflow dans n'importe quel script Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Cela vous donne un bon aperçu de la façon dont les processus sont connectés et de ce qu'ils produisent.

Vous voyez qu'en plus du processus `sayHello` original, nous avons maintenant aussi `convertToUpper` et `collectGreetings`, qui correspondent aux noms des processus que nous avons vus dans la sortie console.
Les deux nouvelles définitions de processus sont structurées de la même manière que le processus `sayHello`, sauf que `collectGreetings` prend un paramètre d'entrée supplémentaire appelé `batch` et produit deux sorties.

Nous n'entrerons pas dans les détails du code pour chacun, mais si vous êtes curieux·se, vous pouvez consulter les détails dans la [Partie 2 de Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Pour l'instant, examinons comment les processus sont connectés les uns aux autres.

#### 2.3.2. Comment les processus sont connectés

La chose vraiment intéressante à regarder ici est comment les appels de processus sont enchaînés dans le bloc `main:` du workflow.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Vous pouvez voir que le premier appel de processus, `sayHello(greeting_ch)`, est inchangé.
Ensuite, l'appel de processus suivant, à `convertToUpper`, fait référence à la sortie de `sayHello` comme `sayHello.out`.

Le modèle est simple : `processName.out` fait référence au canal de sortie d'un processus, qui peut être passé directement au processus suivant.
C'est ainsi que nous transportons les données d'une étape à la suivante dans Nextflow.

#### 2.3.3. Un processus peut prendre plusieurs entrées

Le troisième appel de processus, à `collectGreetings`, est un peu différent.

```groovy title="2b-multistep.nf" linenums="77"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Vous voyez que cet appel reçoit deux entrées, `convertToUpper.out.collect()` et `params.batch`.
En ignorant la partie `.collect()` pour l'instant, nous pouvons généraliser cela comme `collectGreetings(input1, input2)`.

Cela correspond aux deux déclarations d'entrée dans le module de processus :

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Lorsque Nextflow analyse cela, il assignera la première entrée de l'appel à `path input_files`, et la seconde à `val batch_name`.

Donc maintenant vous savez qu'un processus peut prendre plusieurs entrées, et à quoi ressemble l'appel dans le bloc workflow.

Maintenant regardons de plus près cette première entrée, `convertToUpper.out.collect()`.

#### 2.3.4. Ce que fait `collect()` dans l'appel `collectGreetings`

Pour passer la sortie de `sayHello` à `convertToUpper`, nous avons simplement fait référence au canal de sortie de `sayHello` comme `sayHello.out`. Mais pour l'étape suivante, nous voyons une référence à `convertToUpper.out.collect()`.

Qu'est-ce que cette partie `collect()` et que fait-elle ?

C'est un opérateur, bien sûr. Tout comme les opérateurs `splitCsv` et `map` que nous avons rencontrés plus tôt.
Cette fois, l'opérateur s'appelle `collect`, et est appliqué au canal de sortie produit par `convertToUpper`.

L'opérateur `collect` est utilisé pour collecter les sorties de plusieurs appels au même processus et les regrouper en un seul élément de canal.

Dans le contexte de ce workflow, il prend les trois messages d'accueil en majuscules dans le canal `convertToUpper.out` (qui sont trois éléments de canal séparés, et seraient normalement traités dans des appels séparés par le processus suivant) et les regroupe en un seul élément.
C'est ainsi que nous récupérons tous les messages d'accueil dans le même fichier.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

En revanche, si nous n'appliquions pas `collect()` à la sortie de `convertToUpper()` avant de la fournir à `collectGreetings()`, Nextflow exécuterait simplement `collectGreetings()` indépendamment sur chaque message d'accueil, ce qui n'atteindrait pas notre objectif.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Il existe de nombreux autres [opérateurs](https://nextflow.io/docs/latest/reference/operator.html) disponibles pour appliquer des transformations au contenu des canaux entre les appels de processus.

Cela donne aux développeur·ses de pipelines beaucoup de flexibilité pour personnaliser la logique de flux de leur pipeline.
L'inconvénient est que cela peut parfois rendre plus difficile le déchiffrage de ce que fait le pipeline.

#### 2.3.5. Un paramètre d'entrée peut avoir une valeur par défaut

Vous avez peut-être remarqué que `collectGreetings` prend une deuxième entrée, `params.batch` :

```groovy title="2b-multistep.nf" linenums="77"
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Cela passe un paramètre CLI nommé `--batch` au workflow.
Cependant, lorsque nous avons lancé le workflow plus tôt, nous n'avons pas spécifié de paramètre `--batch`.

Que se passe-t-il là ?
Jetez un coup d'œil au bloc `params` :

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Il y a une valeur par défaut configurée dans le workflow, donc nous n'avons pas à la fournir.
Mais si nous en fournissons une sur la ligne de commande, la valeur que nous spécifions sera utilisée à la place de la valeur par défaut.

Essayez :

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Sortie de la commande"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Vous devriez voir de nouvelles sorties finales nommées avec votre nom de lot personnalisé.

??? abstract "Contenu du répertoire"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

C'est un aspect de la configuration des entrées, que nous couvrirons plus en détail dans la Partie 3, mais pour l'instant l'important est de savoir que les paramètres d'entrée peuvent recevoir des valeurs par défaut.

#### 2.3.6. Un processus peut produire plusieurs sorties

Dans la définition du processus `collectGreetings`, nous voyons les déclarations de sortie suivantes :

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Qui sont ensuite référencées par le nom donné avec `emit:` dans le bloc `publish:` :

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Cela facilite ensuite le passage de sorties spécifiques individuellement à d'autres processus dans le workflow, en combinaison avec divers opérateurs.

#### 2.3.7. Les sorties publiées peuvent être organisées

Dans le bloc `output`, nous avons utilisé des chemins personnalisés pour regrouper les résultats intermédiaires afin de faciliter l'identification des sorties finales du workflow.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Il existe des moyens plus sophistiqués d'organiser les sorties publiées ; nous en aborderons quelques-uns dans la partie sur la configuration.

!!! tip "Vous voulez en savoir plus sur la construction de workflows ?"

    Pour une couverture détaillée de la construction de workflows multi-étapes, consultez [Hello Nextflow Partie 3 : Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### À retenir

Vous comprenez à un niveau de base comment les workflows multi-étapes sont construits en utilisant des canaux et des opérateurs et comment ils fonctionnent.
Vous avez également vu que les processus peuvent prendre plusieurs entrées et produire plusieurs sorties, et que celles-ci peuvent être publiées de manière structurée.

### Et ensuite ?

Apprenez comment les pipelines Nextflow peuvent être modularisés pour promouvoir la réutilisation du code et la maintenabilité.

---

## 3. Exécuter des pipelines modularisés

Jusqu'à présent, tous les workflows que nous avons examinés consistaient en un seul fichier de workflow contenant tout le code pertinent.

Cependant, les pipelines du monde réel bénéficient généralement d'être _modularisés_, ce qui signifie que le code est divisé en différents fichiers.
Cela peut rendre leur développement et leur maintenance plus efficaces et durables.

Ici, nous allons démontrer la forme la plus courante de modularité du code dans Nextflow, qui est l'utilisation de **modules**.

Dans Nextflow, un [**module**](https://nextflow.io/docs/latest/module.html) est une définition de processus unique qui est encapsulée par elle-même dans un fichier de code autonome.
Pour utiliser un module dans un workflow, vous ajoutez simplement une instruction d'importation d'une seule ligne à votre fichier de code de workflow ; ensuite vous pouvez intégrer le processus dans le workflow de la même manière que vous le feriez normalement.
Cela permet de réutiliser les définitions de processus dans plusieurs workflows sans produire plusieurs copies du code.

Jusqu'à présent, nous avons exécuté des workflows qui avaient tous leurs processus inclus dans un fichier de code monolithique.
Maintenant nous allons voir à quoi cela ressemble lorsque les processus sont stockés dans des modules individuels.

Nous avons bien sûr encore une fois préparé un workflow approprié à des fins de démonstration, appelé `2c-modules.nf`, ainsi qu'un ensemble de modules situés dans le répertoire `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Contenu du répertoire"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Vous voyez qu'il y a quatre fichiers Nextflow, chacun nommé d'après l'un des processus.
Vous pouvez ignorer le fichier `cowpy.nf` pour l'instant ; nous y reviendrons plus tard.

### 3.1. Examiner le code

Cette fois, nous allons d'abord regarder le code.
Commencez par ouvrir le fichier de workflow `2c-modules.nf`.

??? full-code "Fichier de code complet"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Vous voyez que la logique du workflow est exactement la même que dans la version précédente du workflow.
Cependant, le code du processus a disparu du fichier de workflow, et à la place il y a des instructions `include` pointant vers des fichiers séparés sous `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Ouvrez l'un de ces fichiers et vous trouverez le code pour le processus correspondant.

??? full-code "Fichier de code complet"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

Comme vous pouvez le voir, le code du processus n'a pas changé ; il a juste été copié dans un fichier de module individuel au lieu d'être dans le fichier de workflow principal.
Il en va de même pour les deux autres processus.

Alors voyons à quoi ressemble l'exécution de cette nouvelle version.

### 3.2. Exécuter le workflow

Exécutez cette commande dans votre terminal, avec le flag `-resume` :

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Vous remarquerez que les exécutions de processus ont toutes été mises en cache avec succès, ce qui signifie que Nextflow a reconnu qu'il avait déjà effectué le travail demandé, même si le code a été divisé et que le fichier de workflow principal a été renommé.

Rien de tout cela n'a d'importance pour Nextflow ; ce qui compte, c'est le script de tâche qui est généré une fois que tout le code a été rassemblé et évalué.

!!! tip

    Il est également possible d'encapsuler une section d'un workflow comme un 'sous-workflow' qui peut être importé dans un pipeline plus large, mais cela dépasse le cadre de ce cours.

    Vous pouvez en savoir plus sur le développement de workflows composables dans la Quête secondaire sur les [Workflows de Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### À retenir

Vous savez comment les processus peuvent être stockés dans des modules autonomes pour promouvoir la réutilisation du code et améliorer la maintenabilité.

### Et ensuite ?

Apprenez à utiliser des conteneurs pour gérer les dépendances logicielles.

---

## 4. Utiliser des logiciels conteneurisés

Jusqu'à présent, les workflows que nous avons utilisés comme exemples n'avaient besoin que d'exécuter des opérations de traitement de texte très basiques en utilisant des outils UNIX disponibles dans notre environnement.

Cependant, les pipelines du monde réel nécessitent généralement des outils et des packages spécialisés qui ne sont pas inclus par défaut dans la plupart des environnements.
Habituellement, vous devriez installer ces outils, gérer leurs dépendances et résoudre les conflits éventuels.

Tout cela est très fastidieux et ennuyeux.
Une bien meilleure façon de résoudre ce problème est d'utiliser des **conteneurs**.

Un **conteneur** est une unité de logiciel légère, autonome et exécutable créée à partir d'une **image** de conteneur qui inclut tout ce qui est nécessaire pour exécuter une application, y compris le code, les bibliothèques système et les paramètres.

!!! Tip

    Nous enseignons cela en utilisant la technologie [Docker](https://www.docker.com/get-started/), mais Nextflow prend également en charge plusieurs autres technologies de conteneurs.
    Vous pouvez en savoir plus sur le support de Nextflow pour les conteneurs [ici](https://nextflow.io/docs/latest/container.html).

### 4.1. Utiliser un conteneur directement

Tout d'abord, essayons d'interagir directement avec un conteneur.
Cela aidera à solidifier votre compréhension de ce que sont les conteneurs avant de commencer à les utiliser dans Nextflow.

#### 4.1.1. Télécharger l'image du conteneur

Pour utiliser un conteneur, vous téléchargez généralement ou "tirez" une image de conteneur depuis un registre de conteneurs, puis vous exécutez l'image de conteneur pour créer une instance de conteneur.

La syntaxe générale est la suivante :

```bash title="Syntax"
docker pull '<container>'
```

- `docker pull` est l'instruction au système de conteneurs pour tirer une image de conteneur depuis un dépôt.
- `'<container>'` est l'adresse URI de l'image de conteneur.

À titre d'exemple, téléchargeons une image de conteneur qui contient [cowpy](https://github.com/jeffbuttars/cowpy), une implémentation python d'un outil appelé `cowsay` qui génère de l'art ASCII pour afficher des entrées de texte arbitraires de manière amusante.

Il existe divers dépôts où vous pouvez trouver des conteneurs publiés.
Nous avons utilisé le service [Seqera Containers](https://seqera.io/containers/) pour générer cette image de conteneur Docker à partir du package Conda `cowpy` : `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Exécutez la commande pull complète :

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Sortie de la commande"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
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
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Cela indique au système de télécharger l'image spécifiée.
Une fois le téléchargement terminé, vous avez une copie locale de l'image de conteneur.

#### 4.1.2. Démarrer le conteneur

Les conteneurs peuvent être exécutés comme une commande ponctuelle, mais vous pouvez également les utiliser de manière interactive, ce qui vous donne une invite de shell à l'intérieur du conteneur et vous permet de jouer avec la commande.

La syntaxe générale est la suivante :

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` est l'instruction au système de conteneurs pour démarrer une instance de conteneur à partir d'une image de conteneur et y exécuter une commande.
- `--rm` indique au système d'arrêter l'instance de conteneur après que la commande soit terminée.

Entièrement assemblée, la commande d'exécution du conteneur ressemble à ceci :

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Exécutez cette commande, et vous devriez voir votre invite changer en quelque chose comme `(base) root@b645838b3314:/tmp#`, ce qui indique que vous êtes maintenant à l'intérieur du conteneur.

Vous pouvez vérifier cela en exécutant `ls` pour lister le contenu du répertoire :

```bash
ls /
```

??? success "Sortie de la commande"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Vous voyez que le système de fichiers à l'intérieur du conteneur est différent du système de fichiers sur votre système hôte.

!!! Tip

    Lorsque vous exécutez un conteneur, il est isolé du système hôte par défaut.
    Cela signifie que le conteneur ne peut accéder à aucun fichier sur le système hôte à moins que vous ne l'autorisiez explicitement à le faire en spécifiant que vous voulez monter un volume dans le cadre de la commande `docker run` en utilisant la syntaxe suivante :

    ```bash title="Syntax"
    -v <outside_path>:<inside_path>
    ```

    Cela établit effectivement un tunnel à travers la paroi du conteneur que vous pouvez utiliser pour accéder à cette partie de votre système de fichiers.

    Ceci est couvert plus en détail dans la [Partie 5 de Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Exécuter l'outil `cowpy`

Depuis l'intérieur du conteneur, vous pouvez exécuter la commande `cowpy` directement.

```bash
cowpy "Hello Containers"
```

??? success "Sortie de la commande"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Cela produit de l'art ASCII du personnage de vache par défaut (ou 'cowacter') avec une bulle de dialogue contenant le texte que nous avons spécifié.

Maintenant que vous avez testé l'utilisation de base, vous pouvez essayer de lui donner quelques paramètres.
Par exemple, la documentation de l'outil dit que nous pouvons définir le personnage avec `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Sortie de la commande"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Cette fois, la sortie d'art ASCII montre le pingouin Linux, Tux, parce que nous avons spécifié le paramètre `-c tux`.

Puisque vous êtes à l'intérieur du conteneur, vous pouvez exécuter la commande cowpy autant de fois que vous le souhaitez, en variant les paramètres d'entrée, sans avoir à vous soucier d'installer des bibliothèques sur votre système lui-même.

??? tip "Autres personnages disponibles"

    Utilisez le flag '-c' pour choisir un personnage différent, y compris :

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

N'hésitez pas à jouer avec cela.
Lorsque vous avez terminé, quittez le conteneur en utilisant la commande `exit` :

```bash
exit
```

Vous vous retrouverez dans votre shell normal.

### 4.2. Utiliser un conteneur dans un workflow

Lorsque nous exécutons un pipeline, nous voulons pouvoir dire à Nextflow quel conteneur utiliser à chaque étape, et surtout, nous voulons qu'il gère tout ce travail que nous venons de faire : tirer le conteneur, le démarrer, exécuter la commande et démanteler le conteneur quand c'est terminé.

Bonne nouvelle : c'est exactement ce que Nextflow va faire pour nous.
Nous devons juste spécifier un conteneur pour chaque processus.

Pour démontrer comment cela fonctionne, nous avons créé une autre version de notre workflow qui exécute `cowpy` sur le fichier de messages d'accueil collectés produit dans la troisième étape.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Cela devrait produire un fichier contenant l'art ASCII avec les trois messages d'accueil dans la bulle de dialogue.

#### 4.2.1. Examiner le code

Le workflow est très similaire au précédent, plus l'étape supplémentaire pour exécuter `cowpy`.

??? full-code "Fichier de code complet"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Vous voyez que ce workflow importe un processus `cowpy` depuis un fichier de module, et l'appelle sur la sortie de l'appel `collectGreetings()`, plus un paramètre d'entrée appelé `params.character`.

```groovy title="2d-container.nf" linenums="31"
// generate ASCII art of the greetings with cowpy
cowpy(collectGreetings.out.outfile, params.character)
```

Le processus `cowpy`, qui encapsule la commande cowpy pour générer de l'art ASCII, est défini dans le module `cowpy.nf`.

??? full-code "Fichier de code complet"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

Le processus `cowpy` nécessite deux entrées : le chemin vers un fichier d'entrée contenant le texte à mettre dans la bulle de dialogue (`input_file`), et une valeur pour la variable character.

Surtout, il inclut également la ligne `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, qui pointe vers l'URI du conteneur que nous avons utilisé plus tôt.

#### 4.2.2. Vérifier que Docker est activé dans la configuration

Nous allons légèrement anticiper la Partie 3 de ce cours de formation en introduisant le fichier de configuration `nextflow.config`, qui est l'un des principaux moyens que Nextflow offre pour configurer l'exécution du workflow.
Lorsqu'un fichier nommé `nextflow.config` est présent dans le répertoire actuel, Nextflow le chargera automatiquement et appliquera toute configuration qu'il contient.

À cette fin, nous avons inclus un fichier `nextflow.config` avec une seule ligne de code qui active Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Cette configuration indique à Nextflow d'utiliser Docker pour tout processus qui spécifie un conteneur compatible.

!!! tip

    Il est techniquement possible d'activer l'exécution Docker depuis la ligne de commande, sur une base par exécution, en utilisant le paramètre `-with-docker <container>`.
    Cependant, cela ne nous permet de spécifier qu'un seul conteneur pour l'ensemble du workflow, alors que l'approche que nous venons de vous montrer nous permet de spécifier un conteneur différent par processus.
    Cette dernière est bien meilleure pour la modularité, la maintenance du code et la reproductibilité.

#### 4.2.3. Exécuter le workflow

Juste pour récapituler, voici ce que nous sommes sur le point d'exécuter :

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Pensez-vous que cela va fonctionner ?

Exécutons le workflow avec le flag `-resume`, et spécifions que nous voulons que le personnage soit la dinde.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

Les trois premières étapes ont été mises en cache puisque nous les avons déjà exécutées auparavant, mais le processus `cowpy` est nouveau donc il est réellement exécuté.

Vous pouvez trouver la sortie de l'étape `cowpy` dans le répertoire `results`.

??? abstract "Contenu du fichier"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Vous voyez que le personnage dit tous les messages d'accueil, puisqu'il a été exécuté sur le fichier de messages d'accueil collectés en majuscules.

Plus important encore, nous avons pu exécuter cela dans le cadre de notre pipeline sans avoir à faire une installation appropriée de cowpy et de toutes ses dépendances.
Et nous pouvons maintenant partager le pipeline avec des collaborateur·trices et les faire l'exécuter sur leur infrastructure sans qu'ils aient besoin d'installer quoi que ce soit non plus, à part Docker ou l'une de ses alternatives (comme Singularity/Apptainer) comme mentionné ci-dessus.

#### 4.2.4. Inspecter comment Nextflow a lancé la tâche conteneurisée

Comme coda final à cette section, jetons un coup d'œil au sous-répertoire de travail pour l'un des appels de processus `cowpy` pour obtenir un peu plus d'informations sur la façon dont Nextflow fonctionne avec les conteneurs en coulisses.

Vérifiez la sortie de votre commande `nextflow run` pour trouver le chemin vers le sous-répertoire de travail pour le processus `cowpy`.
En regardant ce que nous avons obtenu pour l'exécution montrée ci-dessus, la ligne de journal de console pour le processus `cowpy` commence par `[7f/caf718]`.
Cela correspond au chemin de répertoire tronqué suivant : `work/7f/caf718`.

Dans ce répertoire, vous trouverez le fichier `.command.run` qui contient toutes les commandes que Nextflow a exécutées en votre nom au cours de l'exécution du pipeline.

??? abstract "Contenu du fichier"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Si vous recherchez `nxf_launch` dans ce fichier, vous devriez voir quelque chose comme ceci :

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Cette commande de lancement montre que Nextflow utilise une commande `docker run` très similaire pour lancer l'appel de processus à celle que nous avons utilisée lorsque nous l'avons exécuté manuellement.
Il monte également le sous-répertoire de travail correspondant dans le conteneur, définit le répertoire de travail à l'intérieur du conteneur en conséquence, et exécute notre script bash modélisé dans le fichier `.command.sh`.

Cela confirme que tout le travail difficile que nous avons dû faire manuellement dans la section précédente est maintenant fait pour nous par Nextflow !

### À retenir

Vous comprenez quel rôle jouent les conteneurs dans la gestion des versions d'outils logiciels et pour assurer la reproductibilité.

Plus généralement, vous avez une compréhension de base de ce que sont les composants principaux des pipelines Nextflow du monde réel et comment ils sont organisés.
Vous connaissez les fondamentaux de la façon dont Nextflow peut traiter plusieurs entrées efficacement, exécuter des workflows composés de plusieurs étapes connectées ensemble, exploiter des composants de code modulaires, et utiliser des conteneurs pour une plus grande reproductibilité et portabilité.

### Et ensuite ?

Prenez une autre pause ! C'était une grosse pile d'informations sur le fonctionnement des pipelines Nextflow.

Dans la dernière section de cette formation, nous allons approfondir le sujet de la configuration.
Vous apprendrez comment configurer l'exécution de votre pipeline pour l'adapter à votre infrastructure ainsi que gérer la configuration des entrées et des paramètres.

---

## Quiz

<quiz>
Pourquoi Nextflow crée-t-il un répertoire de tâche séparé pour chaque appel de processus ?
- [ ] Pour améliorer la vitesse d'exécution
- [ ] Pour réduire l'utilisation de la mémoire
- [x] Pour isoler les exécutions et éviter les collisions entre les sorties
- [ ] Pour activer la compression parallèle des fichiers

En savoir plus : [1.3. Trouver les sorties originales et les journaux](#13-find-the-original-outputs-and-logs)
</quiz>

<quiz>
Que fait l'option `-ansi-log false` lors de l'exécution d'un workflow ?
- [ ] Désactive toute sortie console
- [x] Supprime la couleur de la sortie
- [x] Affiche tous les chemins de répertoire de tâche au lieu de les condenser sur une seule ligne
- [ ] Active le mode de débogage verbeux

En savoir plus : [1.3.2. Faire afficher plus de détails au terminal](#132-make-the-terminal-show-more-details)

Vous pouvez également utiliser l'une des variables d'environnement suivantes si vous préférez ce style :

```bash
export NXF_ANSI_LOG=0
# ou
export NO_COLOR=1
```

</quiz>

<quiz>
Dans le code `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, que fait `#!groovy .map { line -> line[0] }` ?
- [ ] Filtre les lignes vides
- [ ] Trie les lignes par ordre alphabétique
- [x] Extrait la première colonne de chaque ligne CSV
- [ ] Compte le nombre de lignes

En savoir plus : [1.4.1. Charger les données d'entrée à partir du CSV](#141-loading-the-input-data-from-the-csv)
</quiz>

<quiz>
Pourquoi est-il important d'inclure la valeur d'entrée dans les noms de fichiers de sortie (par exemple, `#!groovy "${greeting}-output.txt"`) ?
- [ ] Pour améliorer la vitesse de traitement
- [ ] Pour activer la fonctionnalité de reprise
- [x] Pour empêcher les fichiers de sortie de s'écraser les uns les autres lors du traitement de plusieurs entrées
- [ ] Pour faciliter la compression des fichiers

En savoir plus : [1.4.3. Comment les sorties sont nommées](#143-how-the-outputs-are-named)
</quiz>

<quiz>
Quel est le but de l'instruction `include` dans un workflow modularisé ?
- [ ] Pour copier le code du processus dans le fichier de workflow
- [x] Pour importer une définition de processus depuis un fichier de module externe
- [ ] Pour inclure des paramètres de configuration
- [ ] Pour ajouter des commentaires de documentation

En savoir plus : [3. Exécuter des pipelines modularisés](#3-running-modularized-pipelines)
</quiz>

<quiz>
Lorsque vous modularisez un workflow et l'exécutez avec `-resume`, que se passe-t-il ?
- [ ] La mise en cache est désactivée pour les processus modulaires
- [ ] Toutes les tâches doivent être ré-exécutées
- [x] La mise en cache fonctionne normalement en fonction des scripts de tâche générés
- [ ] Seul le fichier de workflow principal est mis en cache

En savoir plus : [3.2. Exécuter le workflow](#32-run-the-workflow)
</quiz>

<quiz>
Que spécifie la directive `container` dans une définition de processus ?
- [ ] Le répertoire de travail pour le processus
- [ ] L'allocation mémoire maximale
- [x] L'URI de l'image de conteneur à utiliser pour exécuter le processus
- [ ] Le format du fichier de sortie

En savoir plus : [4.2. Utiliser un conteneur dans un workflow](#42-use-a-container-in-a-workflow)
</quiz>

<quiz>
Dans le fichier `.command.run`, que contient la fonction `nxf_launch` ?
- [ ] Les informations de version de Nextflow
- [ ] Les paramètres du workflow
- [x] La commande `docker run` avec les montages de volumes et les paramètres du conteneur
- [ ] Les déclarations d'entrée du processus

En savoir plus : [4.2.4. Inspecter comment Nextflow a lancé la tâche conteneurisée](#424-inspect-how-nextflow-launched-the-containerized-task)
</quiz>

<quiz>
Que gère automatiquement Nextflow lors de l'exécution d'un processus conteneurisé ? (Sélectionnez toutes les réponses qui s'appliquent)
- [x] Tirer l'image du conteneur si nécessaire
- [x] Monter le répertoire de travail dans le conteneur
- [x] Exécuter le script du processus à l'intérieur du conteneur
- [x] Nettoyer l'instance de conteneur après l'exécution

En savoir plus : [4. Utiliser des logiciels conteneurisés](#4-using-containerized-software)
</quiz>
