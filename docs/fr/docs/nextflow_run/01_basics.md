# Partie 1 : Exécuter des opérations de base

Dans cette première partie du cours de formation Nextflow Run, nous abordons le sujet en douceur avec un exemple Hello World très basique et indépendant du domaine, que nous utiliserons pour démontrer les opérations essentielles et identifier les composants de code Nextflow correspondants.

??? info "Qu'est-ce qu'un exemple Hello World ?"

    Un « Hello World ! » est un exemple minimaliste destiné à démontrer la syntaxe de base et la structure d'un langage de programmation ou d'un framework logiciel.
    L'exemple consiste généralement à afficher la phrase « Hello, World ! » sur le périphérique de sortie, tel que la console ou le terminal, ou à l'écrire dans un fichier.

---

## 1. Exécuter un Hello World directement

Démontrons ce concept avec une commande simple que nous exécutons directement dans le terminal, pour montrer ce qu'elle fait avant de l'encapsuler dans Nextflow.

!!! tip

    N'oubliez pas que vous devez maintenant être dans le répertoire `nextflow-run/` comme décrit sur la page [Premiers pas](00_orientation.md).

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

L'exécution de pipelines consiste principalement à lire des données à partir de fichiers et à écrire des résultats dans d'autres fichiers, alors modifions la commande pour écrire la sortie texte dans un fichier afin de rendre l'exemple un peu plus pertinent.

```bash
echo 'Hello World!' > output.txt
```

??? success "Sortie de la commande"

    ```console

    ```

Cela n'affiche rien dans le terminal.

### 1.3. Trouver la sortie

Le texte 'Hello World' devrait maintenant se trouver dans le fichier de sortie que nous avons spécifié, nommé `output.txt`.
Vous pouvez l'ouvrir dans l'explorateur de fichiers ou depuis la ligne de commande en utilisant l'utilitaire `cat`, par exemple.

??? abstract "Contenu du fichier"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

C'est ce que nous allons essayer de reproduire avec notre tout premier workflow Nextflow.

### À retenir

Vous savez maintenant comment exécuter une commande simple dans le terminal qui affiche du texte, et éventuellement, comment faire en sorte qu'elle écrive la sortie dans un fichier.

### Et ensuite ?

Découvrez ce qu'il faut pour exécuter un workflow Nextflow qui obtient le même résultat.

---

## 2. Exécuter le workflow

Nous vous fournissons un script de workflow nommé `1-hello.nf` qui prend un message d'accueil en entrée via un argument de ligne de commande nommé `--input` et produit un fichier texte contenant ce message.

Nous n'allons pas encore examiner le code ; voyons d'abord à quoi ressemble son exécution.

### 2.1. Lancer le workflow et surveiller l'exécution

Dans le terminal, exécutez la commande suivante.

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

La sortie la plus importante ici est la dernière ligne, qui est mise en évidence dans la sortie ci-dessus :

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Cela nous indique que le processus `sayHello` a été exécuté avec succès une fois (`1 of 1 ✔`).

C'est excellent, mais vous vous demandez peut-être : où est la sortie ?

### 2.2. Trouver le fichier de sortie dans le répertoire `results`

Ce workflow est configuré pour publier sa sortie dans un répertoire de résultats.
Si vous regardez votre répertoire actuel, vous verrez que lorsque vous avez exécuté le workflow, Nextflow a créé un nouveau répertoire appelé `results`, ainsi qu'un sous-répertoire appelé `1-hello` sous celui-ci, contenant un fichier appelé `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Ouvrez le fichier ; le contenu doit correspondre à la chaîne que vous avez spécifiée sur la ligne de commande.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

C'est parfait, notre workflow a fait ce qu'il était censé faire !

### 2.3. Enregistrer les résultats dans un répertoire différent

Par défaut, Nextflow enregistrera les sorties du pipeline dans un répertoire appelé `results` dans votre chemin actuel.
Pour modifier l'emplacement de publication de vos fichiers, utilisez le flag CLI `-output-dir` (ou `-o` en abrégé)

!!! danger

    Notez que `--input` a deux tirets et `-output-dir` en a un !
    C'est parce que `--input` est un _paramètre_ du pipeline et `-output-dir` est un flag CLI Nextflow de base.
    Nous reviendrons sur ces éléments plus tard.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Vous devriez voir que vos sorties sont maintenant publiées dans un répertoire appelé `hello_results` au lieu de `results` :

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Les fichiers dans ce répertoire sont exactement les mêmes qu'auparavant, c'est juste le répertoire de niveau supérieur qui est différent.
Cependant, sachez que dans les deux cas, le résultat « publié » est une copie (ou dans certains cas un lien symbolique) de la sortie réelle produite par Nextflow lors de l'exécution du workflow.

Nous allons maintenant jeter un œil sous le capot pour voir où Nextflow a réellement exécuté le travail.

!!! Warning

    Tous les workflows ne seront pas configurés pour publier les sorties dans un répertoire de résultats, et/ou les noms et la structure des répertoires peuvent être différents.
    Un peu plus loin dans cette section, nous vous montrerons comment découvrir où ce comportement est spécifié.

### 2.4. Trouver la sortie originale et les logs dans le répertoire `work/`

Lorsque vous exécutez un workflow, Nextflow crée un « répertoire de tâche » distinct pour chaque invocation de chaque processus dans le workflow (= chaque étape du pipeline).
Pour chacun, il préparera les entrées nécessaires, exécutera la ou les instructions pertinentes et écrira les sorties et les fichiers de log dans ce répertoire unique, qui est nommé automatiquement à l'aide d'un hash afin de le rendre unique.

Tous ces répertoires de tâches se trouveront dans un répertoire appelé `work` dans votre répertoire actuel (où vous exécutez la commande).

Cela peut sembler confus, alors voyons à quoi cela ressemble en pratique.

En revenant à la sortie console du workflow que nous avons exécuté précédemment, nous avions cette ligne :

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Vous voyez comment la ligne commence par `[a3/1e1535]` ?
C'est une forme tronquée du chemin du répertoire de tâche pour cet appel de processus, et cela vous indique où trouver la sortie de l'appel du processus `sayHello` dans le chemin du répertoire `work/`.

Vous pouvez trouver le chemin complet en tapant la commande suivante (en remplaçant `a3/1e1535` par ce que vous voyez dans votre propre terminal) et en appuyant sur la touche tab pour compléter automatiquement le chemin ou en ajoutant un astérisque :

```bash
ls work/a3/1e1535*
```

Cela devrait donner le chemin complet du répertoire : `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Jetons un œil à ce qu'il contient.

??? abstract "Contenu du répertoire"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
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
    Cependant, les fichiers de log sont configurés pour être invisibles dans le terminal, donc si vous voulez utiliser `ls` ou `tree` pour les visualiser, vous devrez définir l'option pertinente pour afficher les fichiers invisibles.

    ```bash
    tree -a work
    ```

Il y a deux ensembles de répertoires dans `work/`, provenant des deux exécutions de pipeline différentes que nous avons effectuées.
Chaque exécution de tâche obtient son propre répertoire isolé dans lequel travailler.
Dans ce cas, le pipeline a fait la même chose les deux fois, donc le contenu de chaque répertoire de tâche est identique

Vous devriez immédiatement reconnaître le fichier `output.txt`, qui est en fait la sortie originale du processus `sayHello` qui a été publiée dans le répertoire `results`.
Si vous l'ouvrez, vous trouverez à nouveau le message `Hello World!`.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

Alors qu'en est-il de tous ces autres fichiers ?

Ce sont les fichiers d'aide et de log que Nextflow a écrits dans le cadre de l'exécution de la tâche :

- **`.command.begin`** : Fichier sentinelle créé dès que la tâche est lancée.
- **`.command.err`** : Messages d'erreur (`stderr`) émis par l'appel du processus
- **`.command.log`** : Sortie de log complète émise par l'appel du processus
- **`.command.out`** : Sortie régulière (`stdout`) par l'appel du processus
- **`.command.run`** : Script complet exécuté par Nextflow pour exécuter l'appel du processus
- **`.command.sh`** : La commande qui a été réellement exécutée par l'appel du processus
- **`.exitcode`** : Le code de sortie résultant de la commande

Le fichier `.command.sh` est particulièrement utile car il vous montre la commande principale que Nextflow a exécutée, sans inclure toute la comptabilité et la configuration de la tâche/environnement.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Cela confirme donc que le workflow a composé la même commande que nous avons exécutée directement sur la ligne de commande précédemment.

Lorsque quelque chose ne va pas et que vous devez résoudre le problème, il peut être utile de regarder le script `command.sh` pour vérifier exactement quelle commande Nextflow a composée en fonction des instructions du workflow, de l'interpolation des variables, etc.

### 2.5. Ré-exécuter le workflow avec différents messages

Essayez de ré-exécuter le workflow plusieurs fois avec différentes valeurs pour l'argument `--input`, puis regardez les répertoires de tâches.

??? abstract "Contenu du répertoire"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
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

### À retenir

Vous savez comment exécuter un script Nextflow simple, surveiller son exécution et trouver ses sorties.

### Et ensuite ?

Apprenez à lire un script Nextflow de base et à identifier comment ses composants sont liés à sa fonctionnalité.

---

## 3. Examiner le script de démarrage du workflow Hello World

Ce que nous avons fait là, c'était essentiellement traiter le script de workflow comme une boîte noire.
Maintenant que nous avons vu ce qu'il fait, ouvrons la boîte et regardons à l'intérieur.

Notre objectif ici n'est pas de mémoriser la syntaxe du code Nextflow, mais de former une intuition de base sur les principaux composants et leur organisation.

### 3.1. Examiner la structure globale du code

Vous trouverez le script `1-hello.nf` dans votre répertoire actuel, qui devrait être `nextflow-run`. Ouvrez-le dans le panneau de l'éditeur.

??? full-code "Fichier de code complet"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
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
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // emit a greeting
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

Un script de workflow Nextflow comprend généralement une ou plusieurs définitions de **process**, le **workflow** lui-même, et quelques blocs optionnels tels que **params** et **output**.

Chaque **process** décrit quelle(s) opération(s) l'étape correspondante du pipeline doit accomplir, tandis que le **workflow** décrit la logique de flux de données qui connecte les différentes étapes.

Examinons d'abord de plus près le bloc **process**, puis nous regarderons le bloc **workflow**.

### 3.2. La définition du `process`

Le premier bloc de code décrit un [**process**](https://nextflow.io/docs/latest/process.html).
La définition du processus commence par le mot-clé `process`, suivi du nom du processus et enfin du corps du processus délimité par des accolades.
Le corps du processus doit contenir un bloc script qui spécifie la commande à exécuter, qui peut être n'importe quoi que vous pourriez exécuter dans un terminal en ligne de commande.

```groovy title="1-hello.nf" linenums="3"
/*
* Use echo to print a greeting to a file
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

Ici, nous avons un **process** appelé `sayHello` qui prend une variable **input** appelée `greeting` et écrit sa **output** dans un fichier nommé `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Il s'agit d'une définition de processus très minimale qui contient simplement une définition `input`, une définition `output` et le `script` à exécuter.

La définition `input` inclut le qualificateur `val`, qui indique à Nextflow de s'attendre à une valeur quelconque (peut être une chaîne, un nombre, peu importe).

La définition `output` inclut le qualificateur `path`, qui indique à Nextflow que cela doit être traité comme un chemin (inclut à la fois les chemins de répertoire et les fichiers).

### 3.3. La définition du `workflow`

Le deuxième bloc de code décrit le [**workflow**](https://nextflow.io/docs/latest/workflow.html) lui-même.
La définition du workflow commence par le mot-clé `workflow`, suivi d'un nom optionnel, puis du corps du workflow délimité par des accolades.

Ici, nous avons un **workflow** qui se compose d'un bloc `main:` et d'un bloc `publish:`.
Le bloc `main:` est le corps principal du workflow et le bloc `publish:` liste les sorties qui doivent être publiées dans le répertoire `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

Dans ce cas, le bloc `main:` contient un appel au processus `sayHello` et lui donne une entrée appelée `params.input` à utiliser comme message d'accueil.

Comme nous le verrons plus en détail dans un instant, `params.input` contient la valeur que nous avons donnée au paramètre `--input` dans notre ligne de commande.

Le bloc `publish:` liste la sortie de l'appel du processus `sayHello()`, qu'il désigne comme `sayHello.out` et donne le nom `first_output` (cela peut être n'importe quoi que l'auteur·trice du workflow souhaite).

Il s'agit d'une définition de **workflow** très minimale.
Dans un pipeline réel, le workflow contient généralement plusieurs appels à des **process** connectés par des **channels**, et il peut y avoir des valeurs par défaut configurées pour les entrées variables.

Nous aborderons cela dans la partie 2 du cours.
Pour l'instant, examinons de plus près comment notre workflow gère les entrées et les sorties.

### 3.4. Le système `params` de paramètres de ligne de commande

Le `params.input` que nous fournissons à l'appel du processus `sayHello()` est un élément astucieux du code Nextflow et mérite qu'on y consacre une minute supplémentaire.

Comme mentionné ci-dessus, c'est ainsi que nous transmettons la valeur du paramètre de ligne de commande `--input` à l'appel du processus `sayHello()`.
En fait, déclarer simplement `params.someParameterName` suffit à donner au workflow un paramètre nommé `--someParameterName` depuis la ligne de commande.

Ici, nous avons formalisé cette déclaration de paramètre en configurant un bloc `params` qui spécifie le type d'entrée que le workflow attend (Nextflow 25.10.2 et versions ultérieures).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

Les types pris en charge incluent `String`, `Integer`, `Float`, `Boolean` et `Path`.
Pour en savoir plus, consultez [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) dans la documentation de référence Nextflow.

!!! tip

    N'oubliez pas que les paramètres de _workflow_ déclarés à l'aide du système `params` prennent toujours deux tirets sur la ligne de commande (`--`).
    Cela les distingue des flags CLI de _niveau Nextflow_, qui ne prennent qu'un seul tiret (`-`).

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

La ligne `mode 'copy'` remplace le comportement par défaut du système, qui consiste à créer un lien symbolique (ou symlink) vers le fichier original dans le répertoire `work/` au lieu d'une copie proprement dite.

Il existe plus d'options que celles affichées ici pour contrôler le comportement de publication ; nous en couvrirons quelques-unes plus tard.
Vous verrez également que lorsqu'un workflow génère plusieurs sorties, chacune est listée de cette manière dans le bloc `output`.

Pour en savoir plus, consultez [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) dans la documentation de référence Nextflow.

??? info "Ancienne syntaxe pour publier les sorties en utilisant `publishDir`"

    Jusqu'à très récemment, la méthode établie pour publier les sorties consistait à le faire au niveau de chaque processus individuel en utilisant une directive `publishDir`.

    Vous trouverez encore ce modèle de code partout dans les anciens pipelines Nextflow et les modules de processus, il est donc important d'en être conscient·e.

    Au lieu d'avoir un bloc `publish:` dans le workflow et un bloc `output` au niveau supérieur, vous verriez une ligne `publishDir` dans la définition du processus `sayHello` :

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
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

### À retenir

Vous savez maintenant comment un workflow Nextflow simple est structuré, et comment les composants de base sont liés à sa fonctionnalité.

### Et ensuite ?

Apprenez à gérer vos exécutions de workflow de manière pratique.

---

## 4. Gérer les exécutions de workflow

Savoir comment lancer des workflows et récupérer les sorties est excellent, mais vous découvrirez rapidement qu'il existe quelques autres aspects de la gestion des workflows qui vous faciliteront la vie.

Ici, nous vous montrons comment tirer parti de la fonctionnalité `resume` lorsque vous devez relancer le même workflow, comment inspecter les logs d'exécution avec `nextflow log`, et comment supprimer les anciens répertoires de travail avec `nextflow clean`.

### 4.1. Relancer un workflow avec `-resume`

Parfois, vous voudrez ré-exécuter un pipeline que vous avez déjà lancé précédemment sans refaire le travail qui a déjà été effectué avec succès.

Nextflow dispose d'une option appelée `-resume` qui vous permet de faire cela.
Plus précisément, dans ce mode, tous les processus qui ont déjà été exécutés avec exactement le même code, les mêmes paramètres et les mêmes entrées seront ignorés.
Cela signifie que Nextflow n'exécutera que les processus que vous avez ajoutés ou modifiés depuis la dernière exécution, ou auxquels vous fournissez de nouveaux paramètres ou entrées.

Il y a deux avantages clés à faire cela :

- Si vous êtes en train de développer un pipeline, vous pouvez itérer plus rapidement puisque vous n'avez qu'à exécuter le ou les processus sur lesquels vous travaillez activement afin de tester vos modifications.
- Si vous exécutez un pipeline en production et que quelque chose ne va pas, dans de nombreux cas, vous pouvez corriger le problème et relancer le pipeline, et il reprendra l'exécution à partir du point de défaillance, ce qui peut vous faire économiser beaucoup de temps et de calcul.

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

Recherchez la partie `cached:` qui a été ajoutée dans la ligne d'état du processus (ligne 5), ce qui signifie que Nextflow a reconnu qu'il avait déjà effectué ce travail et a simplement réutilisé le résultat de l'exécution précédente réussie.

Vous pouvez également voir que le hash du sous-répertoire de travail est le même que dans l'exécution précédente.
Nextflow vous indique littéralement l'exécution précédente et dit « J'ai déjà fait ça là-bas. »

!!! tip

    Lorsque vous ré-exécutez un pipeline avec `resume`, Nextflow n'écrase aucun fichier publié en dehors du répertoire de travail par des exécutions qui ont été exécutées avec succès précédemment.

    Pour en savoir plus, consultez [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) dans la documentation de référence Nextflow.

### 4.2. Inspecter le log des exécutions passées

Chaque fois que vous lancez un workflow nextflow, une ligne est écrite dans un fichier de log appelé `history`, dans un répertoire caché appelé `.nextflow` dans le répertoire de travail actuel.

??? abstract "Contenu du fichier"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Ce fichier vous donne l'horodatage, le nom d'exécution, le statut, l'ID de révision, l'ID de session et la ligne de commande complète pour chaque exécution Nextflow qui a été lancée depuis le répertoire de travail actuel.

Une façon plus pratique d'accéder à ces informations est d'utiliser la commande [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

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

Nextflow utilise l'ID de session pour regrouper les informations de mise en cache d'exécution dans le répertoire `cache`, également situé sous `.nextflow`.

### 4.3. Supprimer les anciens répertoires de travail

Si vous exécutez beaucoup de pipelines, vous pouvez finir par accumuler de très nombreux fichiers dans de nombreux sous-répertoires.
Étant donné que les sous-répertoires sont nommés de manière aléatoire, il est difficile de dire d'après leurs noms lesquels sont des exécutions plus anciennes par rapport aux plus récentes.

Heureusement, Nextflow inclut une commande utile appelée [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) qui peut automatiquement supprimer les sous-répertoires de travail des exécutions passées qui ne vous intéressent plus.

#### 4.3.1. Déterminer les critères de suppression

Il existe plusieurs options pour déterminer ce qu'il faut supprimer, que vous pouvez explorer dans la documentation liée ci-dessus.
Ici, nous vous montrons un exemple qui supprime tous les sous-répertoires des exécutions antérieures à une exécution donnée, spécifiée en utilisant son nom d'exécution.

Recherchez l'exécution réussie la plus récente où vous n'avez pas utilisé `-resume` ; dans notre cas, le nom d'exécution était `backstabbing_swartz`.

Le nom d'exécution est la chaîne en deux parties générée par la machine affichée entre crochets dans la ligne de sortie console `Launching (...)`.
Vous pouvez également utiliser le log Nextflow pour rechercher une exécution en fonction de son horodatage et/ou de sa ligne de commande.

#### 4.3.2. Faire un essai à blanc

Nous utilisons d'abord le flag d'essai à blanc `-n` pour vérifier ce qui sera supprimé avec la commande :

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

Si vous ne voyez aucune ligne de sortie, soit vous n'avez pas fourni un nom d'exécution valide, soit il n'y a pas d'exécutions passées à supprimer. Assurez-vous de changer `backstabbing_swartz` dans l'exemple de commande par le nom d'exécution le plus récent correspondant dans votre log.

#### 4.3.3. Procéder à la suppression

Si la sortie semble conforme aux attentes et que vous souhaitez procéder à la suppression, ré-exécutez la commande avec le flag `-f` au lieu de `-n` :

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Sortie de la commande"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

La sortie devrait être similaire à avant, mais disant maintenant 'Removed' au lieu de 'Would remove'.
Notez que cela ne supprime pas les sous-répertoires à deux caractères (comme `eb/` ci-dessus) mais vide leur contenu.

!!! Warning

    La suppression des sous-répertoires de travail des exécutions passées les supprime du cache de Nextflow et supprime toutes les sorties qui étaient stockées dans ces répertoires.
    Cela signifie que cela brise la capacité de Nextflow à reprendre l'exécution sans ré-exécuter les processus correspondants.

    Vous êtes responsable de la sauvegarde de toutes les sorties qui vous intéressent ! C'est la raison principale pour laquelle nous préférons utiliser le mode `copy` plutôt que le mode `symlink` pour la directive `publish`.

### À retenir

Vous savez comment relancer un pipeline sans répéter les étapes qui ont déjà été exécutées de manière identique, inspecter le log d'exécution et utiliser la commande `nextflow clean` pour nettoyer les anciens répertoires de travail.

### Et ensuite ?

Faites une petite pause ! Vous venez d'absorber les éléments de base de la syntaxe Nextflow et les instructions d'utilisation de base.

Dans la section suivante de cette formation, nous allons examiner quatre versions successivement plus réalistes du pipeline Hello World qui démontreront comment Nextflow vous permet de traiter plusieurs entrées efficacement, d'exécuter des workflows composés de plusieurs étapes connectées ensemble, de tirer parti de composants de code modulaires et d'utiliser des conteneurs pour une plus grande reproductibilité et portabilité.

---

## Quiz

<quiz>
Dans la ligne de sortie console `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, que représente `[a3/7be2fa]` ?
- [ ] Le numéro de version du processus
- [ ] Un identifiant d'exécution unique
- [x] Le chemin tronqué vers le répertoire de travail de la tâche
- [ ] La somme de contrôle du fichier de sortie

En savoir plus : [2.3. Trouver la sortie originale et les logs dans le répertoire `work/`](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Quel est le but du fichier `.command.sh` dans un répertoire de tâche ?
- [ ] Il stocke les paramètres de configuration de la tâche
- [x] Il montre la commande réelle qui a été exécutée par le processus
- [ ] Il contient les messages d'erreur des tâches échouées
- [ ] Il liste les fichiers d'entrée préparés pour la tâche

En savoir plus : [2.3. Trouver la sortie originale et les logs dans le répertoire `work/`](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Qu'arrive-t-il aux résultats publiés lorsque vous ré-exécutez un workflow sans `-resume` ?
- [ ] Ils sont préservés dans des répertoires horodatés séparés
- [x] Ils sont écrasés par la nouvelle exécution
- [ ] Nextflow empêche l'écrasement et échoue
- [ ] Ils sont automatiquement sauvegardés

En savoir plus : [2.4. Ré-exécuter le workflow avec différents messages](#24-re-run-the-workflow-with-different-greetings)
</quiz>

<quiz>
Qu'indique cette sortie console ?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] La tâche a échoué et a été ignorée
- [ ] La tâche est en attente dans une file d'attente
- [x] Nextflow a réutilisé les résultats d'une exécution identique précédente
- [ ] La tâche a été annulée manuellement

En savoir plus : [4.1. Relancer un workflow avec `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Où Nextflow stocke-t-il l'historique d'exécution que la commande `nextflow log` affiche ?
- [ ] Dans le répertoire results
- [ ] Dans le répertoire work
- [x] Dans le fichier `.nextflow/history`
- [ ] Dans `nextflow.config`

En savoir plus : [4.2. Inspecter le log des exécutions passées](#42-inspect-the-log-of-past-executions)
</quiz>

<quiz>
Quel est le but du bloc `params` dans un fichier de workflow ?
- [ ] Définir les exigences en ressources des processus
- [ ] Configurer l'exécuteur
- [x] Déclarer et typer les paramètres d'entrée du workflow
- [ ] Spécifier les options de publication de sortie

En savoir plus : [3.4. Le système `params` de paramètres de ligne de commande](#34-the-params-system-of-command-line-parameters)
</quiz>

<quiz>
Dans le bloc `output` du workflow, que fait `mode 'copy'` ?
- [ ] Crée une sauvegarde du répertoire de travail
- [x] Fait une copie complète des fichiers au lieu de liens symboliques
- [ ] Copie le script de workflow dans les résultats
- [ ] Active la copie incrémentale de fichiers

En savoir plus : [3.5. La directive `publish`](#35-the-publish-directive)
</quiz>

<quiz>
Quel est le flag recommandé à utiliser avec la commande `nextflow clean` avant de supprimer réellement des fichiers ?
- [x] `-n` (essai à blanc) pour prévisualiser ce qui serait supprimé
- [ ] `-v` (verbeux) pour voir la sortie détaillée
- [ ] `-a` (tous) pour sélectionner tous les répertoires
- [ ] `-q` (silencieux) pour supprimer les avertissements

En savoir plus : [4.3. Supprimer les anciens répertoires de travail](#43-delete-older-work-directories)
</quiz>
