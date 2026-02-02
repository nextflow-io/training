# Partie 1 : Exécuter les opérations de base

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette première partie de la formation Nextflow Run, nous abordons le sujet avec un exemple Hello World très basique et indépendant du domaine, que nous utiliserons pour démontrer les opérations essentielles et pointer les composants de code Nextflow correspondants.

??? info "Qu'est-ce qu'un exemple Hello World ?"

    Un « Hello World! » est un exemple minimaliste destiné à démontrer la syntaxe et la structure de base d'un langage de programmation ou d'un framework logiciel.
    L'exemple consiste généralement à afficher la phrase « Hello, World! » sur le dispositif de sortie, comme la console ou le terminal, ou à l'écrire dans un fichier.

---

## 1. Exécuter un Hello World directement

Démontrons ce concept avec une commande simple que nous exécutons directement dans le terminal, pour montrer ce qu'elle fait avant de l'encapsuler dans Nextflow.

!!! tip "Astuce"

    N'oubliez pas que vous devriez maintenant être dans le répertoire `nextflow-run/` comme décrit sur la page [Démarrage](00_orientation.md).

### 1.1. Faire dire bonjour au terminal

Exécutez la commande suivante dans votre terminal.

```bash
echo 'Hello World!'
```

??? success "Sortie de la commande"

    ```console
    Hello World!
    ```

Cela affiche le texte 'Hello World' directement dans le terminal.

### 1.2. Écrire la sortie dans un fichier

L'exécution de pipelines implique principalement la lecture de données depuis des fichiers et l'écriture de résultats dans d'autres fichiers, alors modifions la commande pour écrire la sortie texte dans un fichier afin de rendre l'exemple un peu plus pertinent.

```bash
echo 'Hello World!' > output.txt
```

??? success "Sortie de la commande"

    ```console

    ```

Cela n'affiche rien dans le terminal.

### 1.3. Trouver la sortie

Le texte 'Hello World' devrait maintenant être dans le fichier de sortie que nous avons spécifié, nommé `output.txt`.
Vous pouvez l'ouvrir dans l'explorateur de fichiers ou depuis la ligne de commande en utilisant l'utilitaire `cat`, par exemple.

??? abstract "Contenu du fichier"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

C'est ce que nous allons essayer de reproduire avec notre tout premier workflow Nextflow.

### Récapitulatif

Vous savez maintenant comment exécuter une commande simple dans le terminal qui affiche du texte, et optionnellement, comment lui faire écrire la sortie dans un fichier.

### Et ensuite ?

Découvrez ce qu'il faut pour exécuter un workflow Nextflow qui atteint le même résultat.

---

## 2. Exécuter le workflow

Nous vous fournissons un script de workflow nommé `1-hello.nf` qui prend une salutation en entrée via un argument de ligne de commande nommé `--input` et produit un fichier texte contenant cette salutation.

Nous n'allons pas regarder le code pour l'instant ; voyons d'abord à quoi ressemble son exécution.

### 2.1. Lancer le workflow et surveiller l'exécution

Dans le terminal, exécutez la commande suivante :

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Sortie de la commande"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Si votre sortie console ressemble à cela, alors félicitations, vous venez d'exécuter votre premier workflow Nextflow !

La sortie la plus importante ici est la dernière ligne, qui est mise en surbrillance dans la sortie ci-dessus :

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Cela nous indique que le process `sayHello` a été exécuté avec succès une fois (`1 of 1 ✔`).

C'est bien, mais vous vous demandez peut-être : où est la sortie ?

### 2.2. Trouver le fichier de sortie dans le répertoire `results`

Ce workflow est configuré pour publier sa sortie dans un répertoire de résultats.
Si vous regardez votre répertoire actuel, vous verrez que lorsque vous avez exécuté le workflow, Nextflow a créé un nouveau répertoire appelé `results`, ainsi qu'un sous-répertoire appelé `1-hello` en dessous, contenant un fichier appelé `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Ouvrez le fichier ; le contenu devrait correspondre à la chaîne que vous avez spécifiée sur la ligne de commande.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

C'est bien, notre workflow a fait ce qu'il était censé faire !

Cependant, soyez conscient que le résultat « publié » est une copie (ou dans certains cas un lien symbolique) de la sortie réelle produite par Nextflow lors de l'exécution du workflow.

Alors maintenant, nous allons jeter un coup d'œil sous le capot pour voir où Nextflow a réellement exécuté le travail.

!!! Warning "Avertissement"

    Tous les workflows ne seront pas configurés pour publier les sorties dans un répertoire de résultats, et/ou les noms de répertoires et la structure peuvent être différents.
    Un peu plus loin dans cette section, nous vous montrerons comment découvrir où ce comportement est spécifié.

### 2.3. Trouver la sortie originale et les logs dans le répertoire `work/`

Lorsque vous exécutez un workflow, Nextflow crée un « répertoire de tâche » distinct pour chaque invocation de chaque process dans le workflow (= chaque étape du pipeline).
Pour chacun, il va préparer les entrées nécessaires, exécuter l'instruction ou les instructions pertinentes et écrire les sorties et les fichiers de log dans ce seul répertoire, qui est nommé automatiquement en utilisant un hash afin de le rendre unique.

Tous ces répertoires de tâches vivront sous un répertoire appelé `work` dans votre répertoire actuel (où vous exécutez la commande).

Cela peut sembler confus, alors voyons à quoi cela ressemble en pratique.

En revenant à la sortie console pour le workflow que nous avons exécuté plus tôt, nous avions cette ligne :

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Vous voyez comment la ligne commence par `[a3/7be2fa]` ?
C'est une forme tronquée du chemin du répertoire de tâche pour cet appel de process, et vous indique où trouver la sortie de l'appel au process `sayHello` dans le chemin du répertoire `work/`.

Vous pouvez trouver le chemin complet en tapant la commande suivante (en remplaçant `a3/7be2fa` par ce que vous voyez dans votre propre terminal) et en appuyant sur la touche tab pour compléter automatiquement le chemin ou en ajoutant un astérisque :

```bash
ls work/a3/7be2fa*
```

Cela devrait donner le chemin complet du répertoire : `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Jetons un coup d'œil à ce qu'il y a dedans.

??? abstract "Contenu du répertoire"

    ```console
    work
    └── a3
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Vous ne voyez pas la même chose ?"

    Les noms exacts des sous-répertoires seront différents sur votre système.

    Si vous parcourez le contenu du sous-répertoire de tâche dans l'explorateur de fichiers VSCode, vous verrez tous les fichiers immédiatement.
    Cependant, les fichiers de log sont définis pour être invisibles dans le terminal, donc si vous voulez utiliser `ls` ou `tree` pour les voir, vous devrez définir l'option appropriée pour afficher les fichiers invisibles.

    ```bash
    tree -a work
    ```

Vous devriez immédiatement reconnaître le fichier `output.txt`, qui est en fait la sortie originale du process `sayHello` qui a été publiée dans le répertoire `results`.
Si vous l'ouvrez, vous retrouverez la salutation `Hello World!`.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt"
Hello World!
```

Alors qu'en est-il de tous ces autres fichiers ?

Ce sont les fichiers auxiliaires et de log que Nextflow a écrits dans le cadre de l'exécution de la tâche :

- **`.command.begin`** : Fichier sentinelle créé dès que la tâche est lancée.
- **`.command.err`** : Messages d'erreur (`stderr`) émis par l'appel du process
- **`.command.log`** : Sortie de log complète émise par l'appel du process
- **`.command.out`** : Sortie régulière (`stdout`) de l'appel du process
- **`.command.run`** : Script complet exécuté par Nextflow pour exécuter l'appel du process
- **`.command.sh`** : La commande qui a été réellement exécutée par l'appel du process
- **`.exitcode`** : Le code de sortie résultant de la commande

Le fichier `.command.sh` est particulièrement utile car il vous montre la commande principale que Nextflow a exécutée, sans inclure toute la comptabilité et la configuration de la tâche/environnement.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Cela confirme donc que le workflow a composé la même commande que nous avons exécutée directement sur la ligne de commande plus tôt.

Lorsque quelque chose ne va pas et que vous devez résoudre ce qui s'est passé, il peut être utile de regarder le script `command.sh` pour vérifier exactement quelle commande Nextflow a composée en fonction des instructions du workflow, de l'interpolation de variables, etc.

### 2.4. Ré-exécuter le workflow avec différentes salutations

Essayez de ré-exécuter le workflow quelques fois avec différentes valeurs pour l'argument `--input`, puis regardez les répertoires de tâches.

??? abstract "Contenu du répertoire"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 67
    │   ├── 134e6317f90726c6c17ad53234a32b
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Vous voyez qu'un nouveau sous-répertoire avec un ensemble complet de fichiers de sortie et de log a été créé pour chaque exécution.

En revanche, si vous regardez le répertoire `results`, il n'y a toujours qu'un seul ensemble de résultats, et le contenu du fichier de sortie correspond à ce que vous avez exécuté en dernier.

??? abstract "Contenu du répertoire"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Cela vous montre que les résultats publiés seront écrasés par les exécutions suivantes, alors que les répertoires de tâches sous `work/` sont préservés.

### Récapitulatif

Vous savez comment exécuter un script Nextflow simple, surveiller son exécution et trouver ses sorties.

### Et ensuite ?

Apprenez à lire un script Nextflow basique et à identifier comment ses composants sont liés à sa fonctionnalité.

---

## 3. Examiner le script de démarrage du workflow Hello World

Ce que nous avons fait là-bas était essentiellement de traiter le script de workflow comme une boîte noire.
Maintenant que nous avons vu ce qu'il fait, ouvrons la boîte et regardons à l'intérieur.

Notre objectif ici n'est pas de mémoriser la syntaxe du code Nextflow, mais de former une intuition de base sur les principaux composants et comment ils sont organisés.

### 3.1. Examiner la structure globale du code

Vous trouverez le script `1-hello.nf` dans votre répertoire actuel, qui devrait être `nextflow-run`. Ouvrez-le dans le panneau de l'éditeur.

??? full-code "Fichier de code complet"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Utilise echo pour imprimer 'Hello World!' dans un fichier
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

    /*
    * Paramètres du pipeline
    */
    params {
        input: String
    }

    workflow {

        main:
        // émettre une salutation
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

Un script de workflow Nextflow inclut généralement une ou plusieurs définitions de **process**, le **workflow** lui-même, et quelques blocs optionnels tels que **params** et **output**.

Chaque **process** décrit quelle(s) opération(s) l'étape correspondante dans le pipeline doit accomplir, tandis que le **workflow** décrit la logique de flux de données qui connecte les différentes étapes.

Examinons de plus près le bloc **process** d'abord, puis nous regarderons le bloc **workflow**.

### 3.2. La définition du `process`

Le premier bloc de code décrit un **process**.
La définition du process commence par le mot-clé `process`, suivi du nom du process et enfin le corps du process délimité par des accolades.
Le corps du process doit contenir un bloc script qui spécifie la commande à exécuter, qui peut être n'importe quoi que vous pourriez exécuter dans un terminal de ligne de commande.

```groovy title="1-hello.nf" linenums="3"
/*
* Utilise echo pour imprimer une salutation dans un fichier
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

Ici, nous avons un **process** appelé `sayHello` qui prend une variable d'**entrée** appelée `greeting` et écrit sa **sortie** dans un fichier nommé `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

C'est une définition de process très minimale qui contient juste une définition d'`input`, une définition d'`output` et le `script` à exécuter.

La définition d'`input` inclut le qualificateur `val`, qui indique à Nextflow d'attendre une valeur de quelque type que ce soit (peut être une chaîne, un nombre, peu importe).

La définition d'`output` inclut le qualificateur `path`, qui indique à Nextflow que cela doit être traité comme un chemin (inclut à la fois les chemins de répertoires et les fichiers).

### 3.3. La définition du `workflow`

Le deuxième bloc de code décrit le **workflow** lui-même.
La définition du workflow commence par le mot-clé `workflow`, suivi d'un nom optionnel, puis le corps du workflow délimité par des accolades.

Ici, nous avons un **workflow** qui consiste en un bloc `main:` et un bloc `publish:`.
Le bloc `main:` est le corps principal du workflow et le bloc `publish:` liste les sorties qui doivent être publiées dans le répertoire `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // émettre une salutation
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

Dans ce cas, le bloc `main:` contient un appel au process `sayHello` et lui donne une entrée appelée `params.input` à utiliser comme salutation.

Comme nous le discuterons plus en détail dans un moment, `params.input` contient la valeur que nous avons donnée au paramètre `--input` dans notre ligne de commande.

Le bloc `publish:` liste la sortie de l'appel au process `sayHello()`, qu'il désigne comme `sayHello.out` et lui donne le nom `first_output` (cela peut être n'importe quoi que l'auteur du workflow souhaite).

C'est une définition de **workflow** très minimale.
Dans un pipeline réel, le workflow contient généralement plusieurs appels à des **processes** connectés par des **channels**, et il peut y avoir des valeurs par défaut configurées pour les entrées variables.

Nous aborderons cela dans la Partie 2 de la formation.
Pour l'instant, examinons de plus près comment notre workflow gère les entrées et les sorties.

### 3.4. Le système `params` de paramètres de ligne de commande

Le `params.input` que nous fournissons à l'appel du process `sayHello()` est un morceau de code Nextflow intéressant qui mérite qu'on s'y attarde une minute supplémentaire.

Comme mentionné ci-dessus, c'est ainsi que nous passons la valeur du paramètre de ligne de commande `--input` à l'appel du process `sayHello()`.
En fait, simplement déclarer `params.someParameterName` est suffisant pour donner au workflow un paramètre nommé `--someParameterName` depuis la ligne de commande.

Ici, nous avons formalisé cette déclaration de paramètre en configurant un bloc `params` qui spécifie le type d'entrée que le workflow attend (Nextflow 25.10.2 et ultérieur).

```groovy title="1-hello.nf" linenums="20"
/*
 * Paramètres du pipeline
 */
params {
    input: String
}
```

Les types pris en charge incluent `String`, `Integer`, `Float`, `Boolean` et `Path`.

!!! tip "Astuce"

    Les paramètres de workflow déclarés en utilisant le système `params` prennent toujours deux tirets sur la ligne de commande (`--`).
    Cela les distingue des paramètres de niveau Nextflow, qui ne prennent qu'un seul tiret (`-`).

### 3.5. La directive `publish`

À l'autre extrémité du workflow, nous avons déjà jeté un coup d'œil au bloc `publish:`.
C'est une moitié du système de gestion des sorties ; l'autre moitié est le bloc `output` situé en dessous.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Cela spécifie que la sortie `first_output` listée dans le bloc `publish:` doit être copiée dans un sous-répertoire appelé `1-hello` sous le répertoire de sortie `results` par défaut.

La ligne `mode 'copy'` remplace le comportement par défaut du système, qui est de faire un lien symbolique (ou symlink) vers le fichier original dans le répertoire `work/` au lieu d'une copie propre.

Il y a plus d'options que celles affichées ici pour contrôler le comportement de publication ; nous en couvrirons quelques-unes plus tard.
Vous verrez également que lorsqu'un workflow génère plusieurs sorties, chacune est listée de cette façon dans le bloc `output`.

??? info "Ancienne syntaxe pour publier les sorties en utilisant `publishDir`"

    Jusqu'à très récemment, la façon établie de publier les sorties était de le faire au niveau de chaque process individuel en utilisant une directive `publishDir`.

    Vous trouverez encore ce modèle de code partout dans les anciens pipelines Nextflow et les modules de process, il est donc important d'en être conscient.

    Au lieu d'avoir un bloc `publish:` dans le workflow et un bloc `output` au niveau supérieur, vous verriez une ligne `publishDir` dans la définition du process `sayHello` :

    ```groovy title="Exemple de syntaxe" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Cependant, nous ne recommandons pas d'utiliser cela dans tout nouveau travail car cela sera éventuellement interdit dans les futures versions du langage Nextflow.

### Récapitulatif

Vous savez maintenant comment un workflow Nextflow simple est structuré, et comment les composants de base sont liés à sa fonctionnalité.

### Et ensuite ?

Apprenez à gérer vos exécutions de workflow de manière pratique.

---

## 4. Gérer les exécutions de workflow

Savoir comment lancer des workflows et récupérer les sorties est bien, mais vous trouverez rapidement qu'il y a quelques autres aspects de la gestion des workflows qui vous faciliteront la vie.

Ici, nous vous montrons comment tirer parti de la fonctionnalité `resume` lorsque vous devez relancer le même workflow, comment inspecter les logs d'exécution avec `nextflow log`, et comment supprimer les anciens répertoires de travail avec `nextflow clean`.

### 4.1. Relancer un workflow avec `-resume`

Parfois, vous allez vouloir ré-exécuter un pipeline que vous avez déjà lancé précédemment sans refaire le travail qui a déjà été complété avec succès.

Nextflow a une option appelée `-resume` qui vous permet de faire cela.
Plus précisément, dans ce mode, tout process qui a déjà été exécuté avec exactement le même code, les mêmes paramètres et les mêmes entrées sera ignoré.
Cela signifie que Nextflow n'exécutera que les processes que vous avez ajoutés ou modifiés depuis la dernière exécution, ou auxquels vous fournissez de nouveaux paramètres ou entrées.

Il y a deux avantages clés à faire cela :

- Si vous êtes en train de développer un pipeline, vous pouvez itérer plus rapidement puisque vous n'avez qu'à exécuter le(s) process(es) sur lesquels vous travaillez activement pour tester vos modifications.
- Si vous exécutez un pipeline en production et que quelque chose ne va pas, dans de nombreux cas, vous pouvez corriger le problème et relancer le pipeline, et il reprendra l'exécution à partir du point d'échec, ce qui peut vous faire gagner beaucoup de temps et de calcul.

Pour l'utiliser, ajoutez simplement `-resume` à votre commande et exécutez-la :

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Sortie de la commande"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

La sortie console devrait sembler familière, mais il y a une chose qui est un peu différente par rapport à avant.

Cherchez la partie `cached:` qui a été ajoutée dans la ligne de statut du process (ligne 5), ce qui signifie que Nextflow a reconnu qu'il a déjà fait ce travail et a simplement réutilisé le résultat de l'exécution réussie précédente.

Vous pouvez également voir que le hash du sous-répertoire de travail est le même que dans l'exécution précédente.
Nextflow vous pointe littéralement vers l'exécution précédente et dit « J'ai déjà fait ça là-bas. »

!!! tip "Astuce"

    Lorsque vous ré-exécutez un pipeline avec `resume`, Nextflow n'écrase pas les fichiers publiés en dehors du répertoire de travail par les exécutions qui ont été exécutées avec succès précédemment.

### 4.2. Inspecter le log des exécutions passées

Chaque fois que vous lancez un workflow Nextflow, une ligne est écrite dans un fichier de log appelé `history`, sous un répertoire caché appelé `.nextflow` dans le répertoire de travail actuel.

??? abstract "Contenu du fichier"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Ce fichier vous donne l'horodatage, le nom d'exécution, le statut, l'ID de révision, l'ID de session et la ligne de commande complète pour chaque exécution Nextflow qui a été lancée depuis le répertoire de travail actuel.

Une façon plus pratique d'accéder à ces informations est d'utiliser la commande `nextflow log`.

```bash
nextflow log
```

??? success "Sortie de la commande"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Cela affichera le contenu du fichier de log dans le terminal, augmenté d'une ligne d'en-tête.

Vous remarquerez que l'ID de session change chaque fois que vous exécutez une nouvelle commande `nextflow run`, SAUF si vous utilisez l'option `-resume`.
Dans ce cas, l'ID de session reste le même.

Nextflow utilise l'ID de session pour regrouper les informations de mise en cache d'exécution sous le répertoire `cache`, également situé sous `.nextflow`.

### 4.3. Supprimer les anciens répertoires de travail

Si vous exécutez beaucoup de pipelines, vous pourriez accumuler de très nombreux fichiers dans de nombreux sous-répertoires.
Puisque les sous-répertoires sont nommés de manière aléatoire, il est difficile de dire d'après leurs noms quelles sont les exécutions plus anciennes par rapport aux plus récentes.

Heureusement, Nextflow inclut une sous-commande `clean` utile qui peut automatiquement supprimer les sous-répertoires de travail des exécutions passées dont vous ne vous souciez plus.

#### 4.3.1. Déterminer les critères de suppression

Il existe plusieurs [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) pour déterminer ce qui doit être supprimé.

Ici, nous vous montrons un exemple qui supprime tous les sous-répertoires des exécutions avant une exécution donnée, spécifiée en utilisant son nom d'exécution.

Recherchez l'exécution réussie la plus récente où vous n'avez pas utilisé `-resume` ; dans notre cas, le nom d'exécution était `backstabbing_swartz`.

Le nom d'exécution est la chaîne en deux parties générée par la machine affichée entre crochets dans la ligne de sortie console `Launching (...)`.
Vous pouvez également utiliser le log Nextflow pour rechercher une exécution basée sur son horodatage et/ou sa ligne de commande.

#### 4.3.2. Faire un test à blanc

D'abord, nous utilisons le drapeau de test à blanc `-n` pour vérifier ce qui sera supprimé avec la commande :

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Sortie de la commande"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Votre sortie aura des noms de répertoires de tâches différents et peut avoir un nombre de lignes différent, mais elle devrait ressembler à l'exemple.

Si vous ne voyez aucune ligne en sortie, vous n'avez soit pas fourni un nom d'exécution valide, soit il n'y a pas d'exécutions passées à supprimer. Assurez-vous de changer `backstabbing_swartz` dans la commande d'exemple par le nom d'exécution le plus récent correspondant dans votre log.

#### 4.3.3. Procéder à la suppression

Si la sortie semble comme attendu et que vous voulez procéder à la suppression, ré-exécutez la commande avec le drapeau `-f` au lieu de `-n` :

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Sortie de la commande"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

La sortie devrait être similaire à avant, mais maintenant disant 'Removed' au lieu de 'Would remove'.
Notez que cela ne supprime pas les sous-répertoires de deux caractères (comme `eb/` ci-dessus) mais vide leur contenu.

!!! Warning "Avertissement"

    Supprimer les sous-répertoires de travail des exécutions passées les supprime du cache de Nextflow et supprime toutes les sorties qui étaient stockées dans ces répertoires.
    Cela signifie que cela casse la capacité de Nextflow à reprendre l'exécution sans ré-exécuter les processes correspondants.

    Vous êtes responsable de sauvegarder toutes les sorties qui vous importent ! C'est la raison principale pour laquelle nous préférons utiliser le mode `copy` plutôt que le mode `symlink` pour la directive `publish`.

### Récapitulatif

Vous savez comment relancer un pipeline sans répéter les étapes qui ont déjà été exécutées de manière identique, inspecter le log d'exécution, et utiliser la commande `nextflow clean` pour nettoyer les anciens répertoires de travail.

### Et ensuite ?

Prenez une petite pause ! Vous venez d'absorber les éléments de base de la syntaxe Nextflow et les instructions d'utilisation de base.

Dans la prochaine section de cette formation, nous allons examiner quatre versions successivement plus réalistes du pipeline Hello World qui démontreront comment Nextflow vous permet de traiter plusieurs entrées efficacement, d'exécuter des workflows composés de plusieurs étapes connectées ensemble, d'exploiter des composants de code modulaires, et d'utiliser des conteneurs pour une plus grande reproductibilité et portabilité.

---

## Quiz

<quiz>
Dans la ligne de sortie console `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, que représente `[a3/7be2fa]` ?
- [ ] Le numéro de version du process
- [ ] Un identifiant d'exécution unique
- [x] Le chemin tronqué vers le répertoire de travail de la tâche
- [ ] La somme de contrôle du fichier de sortie

En savoir plus : [2.3. Trouver la sortie originale et les logs dans le répertoire `work/`](#23-trouver-la-sortie-originale-et-les-logs-dans-le-repertoire-work)
</quiz>

<quiz>
Quel est le but du fichier `.command.sh` dans un répertoire de tâche ?
- [ ] Il stocke les paramètres de configuration de la tâche
- [x] Il montre la commande réelle qui a été exécutée par le process
- [ ] Il contient les messages d'erreur des tâches échouées
- [ ] Il liste les fichiers d'entrée préparés pour la tâche

En savoir plus : [2.3. Trouver la sortie originale et les logs dans le répertoire `work/`](#23-trouver-la-sortie-originale-et-les-logs-dans-le-repertoire-work)
</quiz>

<quiz>
Qu'arrive-t-il aux résultats publiés lorsque vous ré-exécutez un workflow sans `-resume` ?
- [ ] Ils sont préservés dans des répertoires horodatés séparés
- [x] Ils sont écrasés par la nouvelle exécution
- [ ] Nextflow empêche l'écrasement et échoue
- [ ] Ils sont automatiquement sauvegardés

En savoir plus : [2.4. Ré-exécuter le workflow avec différentes salutations](#24-re-executer-le-workflow-avec-differentes-salutations)
</quiz>

<quiz>
Qu'indique cette sortie console ?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] La tâche a échoué et a été ignorée
- [ ] La tâche attend dans une file d'attente
- [x] Nextflow a réutilisé les résultats d'une exécution identique précédente
- [ ] La tâche a été annulée manuellement

En savoir plus : [4.1. Relancer un workflow avec `-resume`](#41-relancer-un-workflow-avec--resume)
</quiz>

<quiz>
Où Nextflow stocke-t-il l'historique d'exécution que la commande `nextflow log` affiche ?
- [ ] Dans le répertoire results
- [ ] Dans le répertoire work
- [x] Dans le fichier `.nextflow/history`
- [ ] Dans `nextflow.config`

En savoir plus : [4.2. Inspecter le log des exécutions passées](#42-inspecter-le-log-des-executions-passees)
</quiz>

<quiz>
Quel est le but du bloc `params` dans un fichier de workflow ?
- [ ] Définir les exigences de ressources du process
- [ ] Configurer l'executor
- [x] Déclarer et typer les paramètres d'entrée du workflow
- [ ] Spécifier les options de publication des sorties

En savoir plus : [3.4. Le système params de paramètres de ligne de commande](#34-le-systeme-params-de-parametres-de-ligne-de-commande)
</quiz>

<quiz>
Dans le bloc `output` du workflow, que fait `mode 'copy'` ?
- [ ] Crée une sauvegarde du répertoire de travail
- [x] Fait une copie complète des fichiers au lieu de liens symboliques
- [ ] Copie le script du workflow vers les résultats
- [ ] Active la copie incrémentielle des fichiers

En savoir plus : [3.5. La directive publish](#35-la-directive-publish)
</quiz>

<quiz>
Quel est le drapeau recommandé à utiliser avec la commande `nextflow clean` avant de supprimer réellement des fichiers ?
- [x] `-n` (test à blanc) pour prévisualiser ce qui serait supprimé
- [ ] `-v` (verbeux) pour voir une sortie détaillée
- [ ] `-a` (tous) pour sélectionner tous les répertoires
- [ ] `-q` (silencieux) pour supprimer les avertissements

En savoir plus : [4.3. Supprimer les anciens répertoires de travail](#43-supprimer-les-anciens-repertoires-de-travail)
</quiz>
