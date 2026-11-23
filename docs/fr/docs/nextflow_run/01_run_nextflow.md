# Partie 1 : Exécuter Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette partie, nous présentons les concepts fondamentaux de l'exécution de pipelines Nextflow.
Nous commençons par un simple workflow Hello World, puis nous progressons vers un pipeline complet en plusieurs étapes qui traite plusieurs entrées en parallèle à l'aide de conteneurs.

---

## 1. Hello World

Le workflow `1-hello.nf` prend un message de bienvenue via un argument en ligne de commande et l'écrit dans un fichier.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Lancer le workflow

Exécutez la commande suivante dans votre terminal.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

La ligne clé dans la sortie est la ligne d'état du processus :

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Cela nous indique que le processus `sayHello` s'est exécuté avec succès une fois.
Le préfixe `[6d/740edd]` est un chemin tronqué vers le répertoire de travail de la tâche — nous y reviendrons plus bas.
Le bloc `Outputs:` qui suit liste tous les fichiers publiés par le pipeline, étiquetés selon le bloc `output` présenté dans la section [1.4](#14-optional-code-walkthrough) ci-dessous.

### 1.2. Trouver la sortie

Ce workflow est configuré pour publier ses sorties dans un répertoire `results`.
Après l'exécution, vous devriez y trouver la sortie :

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Ouvrez le fichier pour confirmer qu'il contient `Hello World!`.

### 1.3. Explorer le répertoire `work/`

En coulisses, Nextflow crée un répertoire de tâche unique pour chaque appel de processus dans un répertoire nommé `work/`.
Le hash affiché dans la sortie console (`[6d/740edd]`) est le chemin vers ce répertoire.

```bash
ls work/6d/740edd*
```

À l'intérieur, vous trouverez le fichier de sortie ainsi que plusieurs fichiers de log cachés :

- **`.command.sh`** : la commande exacte exécutée par Nextflow
- **`.command.out`** / **`.command.err`** : stdout et stderr du processus
- **`.command.log`** : sortie de log combinée
- **`.exitcode`** : le code de sortie du processus

Le fichier `.command.sh` est particulièrement utile lors du débogage — il montre précisément ce qui a été exécuté.

### 1.4. Facultatif : Exploration du code

Comprendre le code n'est pas indispensable si vous souhaitez simplement exécuter des pipelines, mais si vous êtes curieux·se, cela vaut la peine d'y jeter un œil.

??? optional "Cliquez pour explorer le code associé à cet exercice"

    Ouvrons `1-hello.nf` et examinons ses principaux composants.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Paramètres du pipeline
     */
    params {
        input: String
    }

    workflow {

        main:
        // émet un message de bienvenue
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Nous observons les éléments suivants :

    - une instruction `include` pointant vers un module de `process`
    - un bloc `params` définissant les paramètres du pipeline
    - un bloc `workflow` décrivant le travail à effectuer
    - un bloc `output` décrivant ce qu'il faut faire avec les sorties

    Examinons chacun d'eux.

    ### Le module `process`

    L'instruction `include` indique à Nextflow de charger un élément appelé `sayHello` depuis un fichier de code séparé.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    Dans ce fichier, nous trouvons la définition d'un processus appelé `sayHello` :

    ```groovy title="modules/sayHello.nf" linenums="4"
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

    Un **process** définit une étape unique dans le pipeline.
    Il déclare ses entrées, ses sorties et le script à exécuter.
    Le qualificateur `val` signifie que l'entrée est une valeur simple (string, nombre, etc.).
    Le qualificateur `path` signifie que la sortie est un chemin de fichier.

    Il est possible d'écrire la définition du processus dans le fichier principal du workflow, mais les conserver dans des fichiers de modules séparés les rend réutilisables : le même module peut être importé par plusieurs scripts de workflow.

    ### Le bloc `params`

    Le bloc `params` déclare les paramètres en ligne de commande acceptés par le workflow :

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Tout paramètre déclaré ici devient disponible en ligne de commande avec un double tiret (`--input`).
    Les types supportés incluent `String`, `Integer`, `Float`, `Boolean` et `Path`.

    !!! tip "Astuce"

        Les paramètres du workflow utilisent toujours deux tirets (`--input`) pour les distinguer des options CLI propres à Nextflow, qui n'en utilisent qu'un (par exemple `-resume`).

    ### Le bloc `workflow`

    Le bloc **workflow** définit la logique de flux de données : quels processus exécuter et dans quel ordre.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // émet un message de bienvenue
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Ici, un seul processus est appelé, donc c'est très simple ; nous aborderons des exemples plus réalistes par la suite.

    La section `main:` appelle le processus `sayHello` avec la valeur `--input`.
    La section `publish:` liste les sorties qui doivent être copiées dans le répertoire de résultats.

    ### Le bloc `output`

    Le bloc `output` en bas du fichier spécifie le chemin de destination et le mode de copie.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Chaque entrée nommée correspond à un label `publish:` dans le workflow et le fait correspondre à un sous-répertoire sous `results/`.

### À retenir

Vous savez comment exécuter un pipeline Nextflow et trouver ses sorties, et vous savez que le travail est exécuté dans des répertoires de tâches sous `work/`.

### Et ensuite ?

Découvrez comment Nextflow gère efficacement plusieurs entrées.

---

## 2. Traiter plusieurs entrées

Les pipelines du monde réel traitent généralement de nombreuses données, pas seulement une seule.
Le workflow `2-inputs.nf` lit depuis un fichier CSV et exécute `sayHello` une fois par ligne, en parallèle.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Exécutons d'abord le workflow, puis nous examinerons le mécanisme utilisé par Nextflow pour gérer ces entrées multiples.

### 2.1. Exécuter le workflow

Exécutez la commande suivante dans votre terminal.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

Le `3 of 3` nous indique que le processus `sayHello` a été appelé trois fois, une fois par ligne dans le CSV.

Dans le répertoire `results`, vous devriez maintenant voir trois fichiers de sortie, un par message de bienvenue :

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Ouvrez l'un des fichiers de sortie pour confirmer que chacun contient un message de bienvenue.

La sortie condensée ci-dessus affiche une seule ligne récapitulative pour `sayHello`, mais Nextflow a en réalité lancé trois exécutions de tâches séparées, une par ligne dans le CSV, et les a exécutées en parallèle dès que votre machine disposait des ressources nécessaires.

Tout comme la tâche unique que vous avez explorée dans la section [1.3](#13-explore-the-work-directory), chacune de ces trois exécutions dispose de son propre répertoire de tâche sous `work/`, complètement isolé des autres :

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Chaque `.command.sh` ne contient que la commande pour ce seul message de bienvenue :

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Cet isolement est ce qui rend l'exécution parallèle sûre : trois tâches s'exécutant simultanément ne partagent jamais un répertoire de travail, donc ce qu'une tâche écrit ne peut pas entrer en conflit avec ce qu'une autre tâche écrit, même si elles produisent des fichiers portant le même nom.
C'est aussi pourquoi `-resume` (abordé ensuite) peut mettre en cache et réutiliser des tâches individuelles de manière indépendante : les entrées, sorties et logs de chaque tâche résident entièrement dans son propre répertoire, sans rien de partagé entre les tâches qui pourrait se désynchroniser.

### 2.2. Exécuter le workflow à nouveau avec `-ansi-log false`

Par défaut, Nextflow condense la sortie en une seule ligne récapitulative par processus.
Pour voir chaque appel de processus listé individuellement, ajoutez `-ansi-log false` :

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Cela affiche les trois appels de processus et le sous-répertoire de travail unique créé pour chacun.

### 2.3. Utiliser `-resume` pour ignorer le travail déjà effectué

Passez maintenant au fichier d'entrée étendu, qui ajoute deux messages de bienvenue supplémentaires, et ajoutez `-resume` à la ligne de commande :

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow n'a exécuté que les deux nouvelles entrées.
Les trois messages de bienvenue traités lors de l'exécution précédente ont été mis en cache et réutilisés automatiquement.

Cela fonctionne également pour ignorer l'exécution de processus pour des étapes déjà terminées avec succès dans un pipeline en plusieurs étapes.
Par exemple, si l'exécution d'un pipeline a été interrompue par une erreur système, ou si vous avez ajouté de nouvelles étapes à un pipeline en cours de développement.

La fonctionnalité `-resume` est particulièrement précieuse dans les longs pipelines où la reprise après une défaillance peut économiser un temps et des ressources critiques.

### 2.4. Facultatif : Exploration du code

Comprendre le code n'est pas indispensable si vous souhaitez simplement exécuter des pipelines, mais si vous êtes curieux·se, cela vaut la peine d'y jeter un œil.

??? optional "Cliquez pour explorer le code associé à cet exercice"

    Le changement clé dans `2-inputs.nf` se trouve dans la section `main:` du workflow :

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // crée un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // émet un message de bienvenue
        sayHello(greeting_ch)
    ```

    Ce que vous voyez ici s'appelle un **canal** : une structure de file d'attente qui gère les données d'entrée d'une manière qui facilite la parallélisation des opérations.

    - `channel.fromPath(params.input)` crée un canal à partir du chemin de fichier fourni avec `--input`
    - `.splitCsv()` analyse le CSV en lignes
    - `#!groovy .map { line -> line[0] }` extrait la première colonne de chaque ligne

    Le résultat est un canal contenant `Hello`, `Bonjour` et `Hola`.
    Lorsqu'il est passé à `sayHello(greeting_ch)`, Nextflow appelle automatiquement le processus une fois par élément, en les exécutant en parallèle lorsque les ressources le permettent.

### À retenir

Vous savez comment traiter plusieurs entrées depuis un fichier CSV en parallèle, et comment utiliser `-resume` pour éviter de répéter le travail déjà effectué.

### Et ensuite ?

Apprenez comment un pipeline complet en plusieurs étapes enchaîne les processus à l'aide de canaux, et comment utiliser des conteneurs pour gérer les outils d'analyse et leurs dépendances.

---

## 3. Exécuter un pipeline en plusieurs étapes

Jusqu'à présent, vous avez exécuté un seul processus, puis l'avez exécuté plusieurs fois en parallèle sur un ensemble d'entrées.
Les pipelines réels vont généralement plus loin : ils enchaînent plusieurs processus, en alimentant la sortie de l'un vers le suivant, et s'appuient souvent sur plusieurs logiciels différents.
Le workflow `main.nf` combine ces deux aspects dans un pipeline complet.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Chaque message de bienvenue en entrée passe par les quatre étapes : `sayHello` l'écrit dans un fichier, `convertToUpper` convertit le texte en majuscules, `collectGreetings` fusionne tous les résultats en un seul fichier, et `cowpy` génère de l'art ASCII à partir de la sortie fusionnée à l'aide d'un outil conteneurisé.
Nextflow relie ces étapes entre elles avec des canaux : la sortie d'un processus devient l'entrée du suivant, de sorte que toute la chaîne s'exécute automatiquement au fur et à mesure que les données deviennent disponibles, sans que vous ayez à orchestrer chaque étape manuellement.

Notez que ce workflow utilise des modules : chaque processus est défini dans son propre fichier sous `modules/`, et `main.nf` les importe avec des instructions `include` plutôt que de les définir directement.
Cela rend chaque processus réutilisable dans plusieurs workflows sans dupliquer le code. Pour en savoir plus, consultez la section d'exploration du code ci-dessous.

### 3.1. Exécuter le workflow

Exécutez la commande suivante dans votre terminal.

```bash
nextflow run main.nf --input data/greetings.csv
```

Le paramètre `character` est défini par défaut à `turkey` dans `nextflow.config`, donc l'art ASCII utilise une dinde sauf si vous le remplacez (essayez d'ajouter `--character tux`).

??? success "Sortie de la commande"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Quatre processus se sont exécutés, mais pas le même nombre de fois.
`sayHello` et `convertToUpper` se sont chacun exécutés une fois par entrée (3 of 3) : chaque message de bienvenue doit être écrit et converti en majuscules individuellement.
`collectGreetings` et `cowpy` ne se sont chacun exécutés qu'une seule fois (1 of 1) : la fusion des messages de bienvenue et la génération de l'art ASCII n'ont de sens qu'une fois que tous les résultats individuels sont disponibles.
Cette forme de dispersion puis de convergence — plusieurs tâches parallèles alimentant un nombre réduit de tâches en aval — est courante dans les pipelines réels.

Nextflow n'attend pas qu'une étape entière soit terminée avant de démarrer la suivante.
Dès qu'une sortie de `sayHello` est prête, la tâche `convertToUpper` correspondante peut démarrer, de sorte que les tâches de différents processus s'exécutent de manière concurrente plutôt qu'en lots stricts.
`collectGreetings` et `cowpy` doivent attendre, car chacun dépend de la disponibilité de tous les résultats en amont.

Le répertoire `results` reflète cette convergence, ainsi que ce que l'auteur·trice du pipeline a choisi de publier et où : rappelons le bloc `output` de l'exploration du code en 1.4, qui est ce qui définit cette structure.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

Le répertoire de premier niveau est nommé d'après le paramètre `batch`, dont la valeur par défaut est `batch` ; vous le verrez changer dans les exercices ultérieurs.

Consultez `cowpy-COLLECTED-batch-output.txt` pour le fichier d'art ASCII.

??? abstract "Contenu du fichier"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
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

Tout comme dans la section [2.1](#21-run-the-workflow), chacune de ces 8 exécutions de tâches, réparties sur les quatre processus, dispose de son propre répertoire sous `work/`, complètement isolé des autres.
`collectGreetings` illustre bien pourquoi cela est important : il dépend des sorties des trois tâches `convertToUpper`, qui résident dans trois répertoires de tâches différents. Nextflow crée donc des liens symboliques vers ces fichiers dans le propre répertoire de `collectGreetings`, plutôt que de lui faire lire directement depuis les répertoires de ses tâches en amont :

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Chaque tâche ne voit que les fichiers spécifiques dont elle a besoin, quelle qu'en soit la provenance, et jamais le contenu interne du répertoire d'une autre tâche.
Sur l'ensemble d'un pipeline, ce même isolement que vous avez observé avec un seul processus dans la section [2.1](#21-run-the-workflow) est ce qui permet à Nextflow d'exécuter chaque tâche de chaque processus de manière concurrente, en toute sécurité.

!!! note "Note"

    L'étape `cowpy` s'exécute dans un conteneur Docker plutôt que de s'appuyer sur un logiciel installé localement.
    Un conteneur regroupe une application avec tout ce dont elle a besoin pour s'exécuter, de sorte que vous n'avez pas à installer et gérer les dépendances vous-même, et le pipeline se comporte de la même manière sur n'importe quelle machine capable d'exécuter le conteneur.
    Nextflow prend également en charge Conda comme alternative aux conteneurs ; consultez la [Partie 2](./02_configure_pipeline.md) pour savoir comment basculer entre les deux.

### 3.2. Facultatif : Exploration du code

Comprendre le code n'est pas indispensable si vous souhaitez simplement exécuter des pipelines, mais si vous êtes curieux·se, cela vaut la peine d'y jeter un œil.

??? optional "Cliquez pour explorer le code associé à cet exercice"

    ### Comment les données circulent d'une étape à l'autre

    Chaque processus passe son canal de sortie au suivant :

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // crée un canal pour les entrées depuis un fichier CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    Le motif `processName.out` fait référence au canal de sortie d'un processus.

    L'opérateur `.collect()` rassemble toutes les sorties individuelles de `convertToUpper` en un seul élément de canal avant de les passer à `collectGreetings`.

    ### Utilisation des modules de processus

    `main.nf` ne définit aucun code de processus directement.
    À la place, il importe chaque processus depuis son propre fichier sous `modules/` :

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Chaque fichier de module contient une seule définition de processus, structurée de la même manière que le module `sayHello` dans la section [1.4](#14-optional-code-walkthrough).
    Conserver les processus dans des fichiers séparés les rend réutilisables dans plusieurs workflows sans dupliquer le code.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Utilisation de logiciels conteneurisés

    Le processus `cowpy` s'exécute dans un conteneur Docker spécifié dans son fichier de module :

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
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

    Nextflow extrait automatiquement l'image, exécute le script dans le conteneur et effectue le nettoyage ensuite.
    Docker est activé pour ce projet dans `nextflow.config` :

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Cette seule ligne active Docker pour tout processus du pipeline qui dispose d'un conteneur spécifié.

### À retenir

Vous avez exécuté un pipeline complet en plusieurs étapes qui traite plusieurs entrées en parallèle à l'aide d'un outil conteneurisé.

### Et ensuite ?

Passez à la [Partie 2](./02_configure_pipeline.md), où vous apprendrez à configurer le comportement du pipeline à l'aide de `nextflow.config`.

---

## Résumé

Dans cette partie, vous avez appris à :

- Exécuter un workflow Nextflow et trouver ses sorties
- Explorer le répertoire `work/` et ses fichiers de log
- Traiter plusieurs entrées depuis un fichier CSV en parallèle
- Utiliser `-resume` pour ignorer le travail déjà effectué lors de l'ajout de nouvelles entrées
- Exécuter un pipeline en plusieurs étapes utilisant un outil conteneurisé
