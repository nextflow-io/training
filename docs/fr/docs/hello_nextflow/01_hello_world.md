# Partie 1 : Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir [la playlist complète](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sur la chaîne YouTube Nextflow.

:green_book: La transcription de la vidéo est disponible [ici](./transcripts/01_hello_world.md).
///
-->

Dans cette première partie du cours de formation Hello Nextflow, nous abordons le sujet avec un exemple Hello World très basique et indépendant du domaine, que nous allons progressivement enrichir pour démontrer l'utilisation de la logique et des composants fondamentaux de Nextflow.

??? info "Qu'est-ce qu'un exemple Hello World ?"

    Un « Hello World ! » est un exemple minimaliste destiné à démontrer la syntaxe de base et la structure d'un langage de programmation ou d'un framework logiciel.
    L'exemple consiste généralement à imprimer la phrase « Hello, World! » sur le périphérique de sortie, comme la console ou le terminal, ou à l'écrire dans un fichier.

---

## 0. Échauffement : Exécuter un exemple Hello World directement

Démontrons ceci avec une commande simple que nous exécutons directement dans le terminal, pour montrer ce qu'elle fait avant de l'encapsuler dans Nextflow.

!!! tip "Astuce"

    N'oubliez pas que vous devriez maintenant être dans le répertoire `hello-nextflow/` comme décrit sur la page [Premiers pas](00_orientation.md).

### 0.1. Faire dire bonjour au terminal

Exécutez la commande suivante dans votre terminal.

```bash
echo 'Hello World!'
```

??? success "Sortie de la commande"

    ```console
    Hello World!
    ```

Cela affiche le texte « Hello World » directement dans le terminal.

### 0.2. Écrire la sortie dans un fichier

L'exécution de pipelines implique principalement la lecture de données à partir de fichiers et l'écriture de résultats dans d'autres fichiers, alors modifions la commande pour écrire la sortie texte dans un fichier afin de rendre l'exemple un peu plus pertinent.

```bash
echo 'Hello World!' > output.txt
```

??? success "Sortie de la commande"

    ```console

    ```

Cela n'affiche rien dans le terminal.

### 0.3. Trouver la sortie

Le texte « Hello World » devrait maintenant être dans le fichier de sortie que nous avons spécifié, nommé `output.txt`.
Vous pouvez l'ouvrir dans l'explorateur de fichiers ou depuis la ligne de commande en utilisant l'utilitaire `cat`, par exemple.

??? abstract "Contenu du fichier"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

C'est ce que nous allons essayer de reproduire avec notre tout premier flux de travail Nextflow.

### À retenir

Vous savez maintenant comment exécuter une commande simple dans le terminal qui affiche du texte, et optionnellement, comment lui faire écrire la sortie dans un fichier.

### Et ensuite ?

Découvrez à quoi cela ressemblerait écrit comme un flux de travail Nextflow.

---

## 1. Examiner le script et l'exécuter

Nous vous fournissons un script de flux de travail pleinement fonctionnel, bien que minimaliste, nommé `hello-world.nf` qui fait la même chose qu'avant (écrire « Hello World! ») mais avec Nextflow.

Pour commencer, ouvrons le script de flux de travail afin que vous puissiez avoir une idée de sa structure.
Ensuite, nous l'exécuterons et chercherons ses sorties.

### 1.1. Examiner le code

Vous trouverez le script `hello-world.nf` dans votre répertoire actuel, qui devrait être `hello-nextflow`. Ouvrez-le dans le panneau de l'éditeur.

??? full-code "Fichier de code complet"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Utilise echo pour imprimer 'Hello World!' dans un fichier
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // émettre une salutation
        sayHello()
    }
    ```

Un script de flux de travail Nextflow comprend généralement une ou plusieurs définitions de **processus** et le **flux de travail** lui-même, plus quelques blocs optionnels (non présents ici) que nous présenterons plus tard.

Chaque **processus** décrit quelle(s) opération(s) l'étape correspondante dans le pipeline devrait accomplir, tandis que le **workflow** décrit la logique de flux de données qui connecte les différentes étapes.

Nous allons d'abord examiner de plus près le bloc **process**, puis nous examinerons le bloc **workflow**.

#### 1.1.1. La définition du `process`

Le premier bloc de code décrit un **processus**.

La définition du processus commence par le mot-clé `process`, suivi du nom du processus et enfin du corps du processus délimité par des accolades.
Le corps du processus doit contenir un bloc script qui spécifie la commande à exécuter, qui peut être n'importe quoi que vous pourriez exécuter dans un terminal de ligne de commande.

```groovy title="hello-world.nf" linenums="3"
/*
* Utilise echo pour imprimer 'Hello World!' dans un fichier
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Ici, nous avons un **processus** appelé `sayHello` qui écrit sa **sortie** dans un fichier nommé `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

C'est une définition de processus très minimale qui contient juste une définition de `output` et le `script` à exécuter.

La définition de `output` inclut le qualificateur `path`, qui indique à Nextflow que cela doit être traité comme un chemin (inclut à la fois les chemins de répertoires et les fichiers).
Un autre qualificateur courant est `val`.

Il est important de noter que la définition de sortie ne _détermine_ pas quelle sortie sera créée.
Elle _déclare_ simplement quelle est la sortie attendue, afin que Nextflow puisse la rechercher une fois l'exécution terminée.
Ceci est nécessaire pour vérifier que la commande a été exécutée avec succès et pour transmettre la sortie aux processus en aval si nécessaire. Les sorties produites qui ne correspondent pas à ce qui est déclaré dans le bloc output ne seront pas transmises aux processus en aval.

!!! warning "Avertissement"

    Cet exemple est fragile car nous avons codé en dur le nom du fichier de sortie à deux endroits séparés (les blocs script et output).
    Si nous changeons l'un mais pas l'autre, le script ne fonctionnera plus.
    Plus tard, vous apprendrez des moyens d'utiliser des variables pour atténuer ce problème.

Dans un pipeline réel, un processus contient généralement des blocs supplémentaires tels que des directives et des entrées, que nous présenterons dans un instant.

#### 1.1.2. La définition du `workflow`

Le deuxième bloc de code décrit le **workflow** lui-même.
La définition du workflow commence par le mot-clé `workflow`, suivi d'un nom optionnel, puis du corps du workflow délimité par des accolades.

Ici, nous avons un **workflow** qui consiste en un bloc `main:` (qui dit « ceci est le corps principal du workflow ») contenant un appel au processus `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // émettre une salutation
    sayHello()
}
```

C'est une définition de **workflow** très minimale.
Dans un pipeline réel, le workflow contient généralement plusieurs appels à des **processus** connectés par des **canaux**, et les processus attendent une ou plusieurs **entrée(s)** variable(s).

Vous apprendrez comment ajouter des entrées variables plus tard dans ce module de formation ; et vous apprendrez comment ajouter plus de processus et les connecter par des canaux dans la Partie 3 de ce cours.

!!! tip "Astuce"

    Techniquement, la ligne `main:` n'est pas requise pour les workflows simples comme celui-ci, vous pouvez donc rencontrer des workflows qui ne l'ont pas.
    Mais nous en aurons besoin pour profiter des sorties au niveau du workflow, donc autant l'inclure dès le départ.

### 1.2. Exécuter le flux de travail

Regarder du code n'est pas aussi amusant que de l'exécuter, alors essayons cela en pratique.

#### 1.2.1. Lancer le flux de travail et surveiller l'exécution

Dans le terminal, exécutez la commande suivante :

```bash
nextflow run hello-world.nf
```

??? success "Sortie de la commande"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Si votre sortie console ressemble à cela, alors félicitations, vous venez d'exécuter votre premier flux de travail Nextflow !

La sortie la plus importante ici est la dernière ligne, qui est mise en évidence dans la sortie ci-dessus :

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Cela nous indique que le processus `sayHello` a été exécuté avec succès une fois (`1 of 1 ✔`).

De manière importante, cette ligne vous indique également où trouver la sortie de l'appel au processus `sayHello`.
Examinons cela maintenant.

#### 1.2.2. Trouver la sortie et les logs dans le répertoire `work`

Lorsque vous exécutez Nextflow pour la première fois dans un répertoire donné, il crée un répertoire appelé `work` où il écrira tous les fichiers (et tous les liens symboliques) générés au cours de l'exécution.

Dans le répertoire `work`, Nextflow organise les sorties et les logs par appel de processus.
Pour chaque appel de processus, Nextflow crée un sous-répertoire imbriqué, nommé avec un hash pour le rendre unique, où il préparera toutes les entrées nécessaires (en utilisant des liens symboliques par défaut), écrira les fichiers d'aide, et écrira les logs et toutes les sorties du processus.

Le chemin vers ce sous-répertoire est affiché sous forme tronquée entre crochets dans la sortie console.
En regardant ce que nous avons obtenu pour l'exécution montrée ci-dessus, la ligne de log console pour le processus sayHello commence par `[65/7be2fa]`. Cela correspond au chemin de répertoire suivant : `work/65/7be2fad5e71e5f49998f795677fd68`

Voyons ce qu'il y a dedans.

??? abstract "Contenu du répertoire"

    ```console
    work
    └── 65
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
    Cependant, les fichiers de log sont configurés pour être invisibles dans le terminal, donc si vous voulez utiliser `ls` ou `tree` pour les voir, vous devrez définir l'option pertinente pour afficher les fichiers invisibles.

    ```bash
    tree -a work
    ```

La première chose que vous voulez regarder est la sortie réelle du flux de travail, c'est-à-dire le fichier `output.txt` produit par le processus `sayHello`.
Ouvrez-le et vous trouverez la salutation `Hello World!`, qui était le but de notre flux de travail minimaliste.

??? abstract "Contenu du fichier"

    ```console title="output.txt"
    Hello World!
    ```

Ça a fonctionné !

Certes, cela peut sembler beaucoup de code d'enveloppe pour un si petit résultat, mais la valeur de tout ce code d'enveloppe deviendra plus évidente une fois que nous commencerons à lire des fichiers d'entrée et à enchaîner plusieurs étapes.

Cela dit, examinons également les autres fichiers dans ce répertoire. Ce sont des fichiers d'aide et de log produits par Nextflow dans le cadre de l'exécution de la tâche.

- **`.command.begin`** : Métadonnées liées au début de l'exécution de l'appel de processus
- **`.command.err`** : Messages d'erreur (`stderr`) émis par l'appel de processus
- **`.command.log`** : Sortie de log complète émise par l'appel de processus
- **`.command.out`** : Sortie régulière (`stdout`) de l'appel de processus
- **`.command.run`** : Script complet exécuté par Nextflow pour exécuter l'appel de processus
- **`.command.sh`** : La commande qui a été réellement exécutée par l'appel de processus
- **`.exitcode`** : Le code de sortie résultant de la commande

Le fichier `.command.sh` est particulièrement utile car il vous indique la commande principale que Nextflow a exécutée, sans inclure toute la comptabilité et la configuration de la tâche/de l'environnement.

??? abstract "Contenu du fichier"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Cela correspond à ce que nous avons exécuté manuellement plus tôt.

Dans ce cas, c'est très simple car la commande du processus était codée en dur, mais plus tard dans le cours, vous verrez des commandes de processus qui impliquent une interpolation de variables.
Cela rend particulièrement précieux de pouvoir voir exactement comment Nextflow a interprété le code et quelle commande a été produite lorsque vous dépannez une exécution échouée.

### 1.3. Exécuter à nouveau le flux de travail

Essayez de ré-exécuter le flux de travail plusieurs fois, puis regardez les répertoires de tâches sous `work/`.

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
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
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
Cela vous montre que l'exécution du même flux de travail plusieurs fois n'écrasera pas les résultats des exécutions précédentes.

### À retenir

Vous savez comment déchiffrer un script Nextflow simple, l'exécuter et trouver la sortie et les fichiers de log pertinents dans le répertoire work.

### Et ensuite ?

Apprendre à publier les sorties du flux de travail dans un emplacement plus pratique.

---

## 2. Publier les sorties

Comme vous venez de l'apprendre, la sortie produite par notre pipeline est enfouie dans un répertoire de travail situé plusieurs niveaux plus bas.
Ceci est fait exprès ; Nextflow contrôle ce répertoire et nous ne sommes pas censés interagir avec lui.
Cependant, cela rend peu pratique la récupération des sorties qui nous importent.

Heureusement, Nextflow fournit un moyen de publier les sorties dans un répertoire désigné en utilisant des [définitions de sortie au niveau du workflow](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Utilisation de base

Cela va impliquer deux nouveaux morceaux de code :

1. Un bloc `publish:` à l'intérieur du corps du `workflow`, déclarant les sorties de processus.
2. Un bloc `output` au script spécifiant les options de sortie telles que le mode et l'emplacement.

#### 2.1.1. Déclarer la sortie du processus `sayHello`

Nous devons ajouter un bloc `publish:` au corps du workflow (même type d'élément de code que le bloc `main:`) et lister la sortie du processus `sayHello()`.

Dans le fichier de script de workflow `hello-world.nf`, ajoutez les lignes de code suivantes :

=== "Après"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // émettre une salutation
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // émettre une salutation
        sayHello()
    }
    ```

Vous voyez que nous pouvons faire référence à la sortie du processus simplement en faisant `sayHello().out`, et lui attribuer un nom arbitraire, `first_output`.

#### 2.1.2. Ajouter un bloc `output:` au script

Maintenant nous avons juste besoin d'ajouter le bloc `output:` où le chemin du répertoire de sortie sera spécifié. Notez que ce nouveau bloc se situe **en dehors** et **en dessous** du bloc `workflow` dans le script.

Dans le fichier de script de workflow `hello-world.nf`, ajoutez les lignes de code suivantes :

=== "Après"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // émettre une salutation
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // émettre une salutation
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Nous pouvons utiliser ceci pour attribuer des chemins spécifiques à toutes les sorties de processus déclarées dans le bloc `workflow`.
Plus tard, vous apprendrez des moyens de générer des structures de répertoires de sortie sophistiquées, mais pour l'instant, nous codons juste en dur un chemin minimal pour plus de simplicité.

#### 2.1.3. Exécuter le flux de travail

Maintenant exécutez le script de workflow modifié :

```bash
nextflow run hello-world.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

La sortie du terminal devrait sembler familière. Extérieurement, rien n'a changé.

Cependant, vérifiez votre explorateur de fichiers : cette fois, Nextflow a créé un nouveau répertoire appelé `results/`.

??? abstract "Contenu du répertoire"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

À l'intérieur du répertoire `results`, nous trouvons un lien symbolique vers le `output.txt` produit dans le répertoire work par la commande que nous venons d'exécuter.

Cela nous permet de récupérer facilement les fichiers de sortie sans avoir à fouiller dans le sous-répertoire work.

### 2.2. Définir un emplacement personnalisé

Avoir un emplacement par défaut est bien, mais vous pourriez vouloir personnaliser où les résultats sont enregistrés et comment ils sont organisés.

Par exemple, vous pouvez vouloir organiser vos sorties dans des sous-répertoires.
La façon la plus simple de le faire est d'attribuer un chemin de sortie spécifique par sortie.

#### 2.2.1. Modifier le chemin de sortie

Encore une fois, modifier le comportement de publication pour une sortie spécifique est vraiment simple.
Pour définir un emplacement personnalisé, modifiez simplement le `path` en conséquence :

=== "Après"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Puisque ceci est défini au niveau de la sortie individuelle, vous pouvez spécifier différents emplacements et sous-répertoires selon vos besoins.

#### 2.2.2. Exécuter à nouveau le flux de travail

Essayons-le.

```bash
nextflow run hello-world.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Cette fois, le résultat est écrit sous le sous-répertoire spécifié.

??? abstract "Contenu du répertoire"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Vous voyez que le résultat de l'exécution précédente est toujours là.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Vous pouvez utiliser autant de niveaux d'imbrication que vous le souhaitez.
Il est également possible d'utiliser le nom du processus ou d'autres variables pour nommer les répertoires utilisés pour organiser les résultats, et il est possible de changer le nom par défaut du répertoire de sortie de niveau supérieur (qui est contrôlé par la variable spéciale `outputDir`).
Nous couvrirons ces options dans des formations ultérieures.

### 2.3. Définir le mode de publication sur copy

Par défaut, les sorties sont publiées comme des liens symboliques depuis le répertoire `work`.
Cela signifie qu'il n'y a qu'un seul fichier sur le système de fichiers.

C'est très bien lorsque vous travaillez avec de très gros fichiers, pour lesquels vous ne voulez pas stocker plusieurs copies.
Cependant, si vous supprimez le répertoire work à un moment donné (nous couvrirons les opérations de nettoyage sous peu), vous perdrez l'accès au fichier.
Vous devez donc avoir un plan pour sauvegarder des copies de tous les fichiers importants dans un endroit sûr.

Une option facile est de basculer le mode de publication sur copy pour les sorties qui vous importent.

#### 2.3.1. Ajouter la directive mode

Cette partie est vraiment simple.
Ajoutez simplement `mode 'copy'` à la définition de sortie au niveau du workflow concernée :

=== "Après"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Cela définit le mode de publication pour cette sortie spécifique.

#### 2.3.2. Exécuter à nouveau le flux de travail

Essayons-le.

```bash
nextflow run hello-world.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Cette fois, si vous regardez les résultats, le fichier est une vraie copie au lieu d'un simple lien symbolique.

??? abstract "Contenu du répertoire"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Puisque ceci aussi est défini au niveau de la sortie individuelle, cela vous permet de définir le mode de publication de manière granulaire.
Cela sera particulièrement pratique plus tard lorsque nous passerons à des pipelines à plusieurs étapes, où vous pourriez vouloir copier uniquement les sorties finales et laisser les sorties intermédiaires comme liens symboliques, par exemple.

Comme noté précédemment, il existe d'autres options plus sophistiquées pour contrôler comment les sorties sont publiées.
Nous vous montrerons comment les utiliser en temps voulu dans votre parcours Nextflow.

### 2.4. Note sur les directives `publishDir` au niveau du processus

Jusqu'à très récemment, la façon établie de publier les sorties était de le faire au niveau de chaque processus individuel en utilisant une directive `publishDir`.

Pour obtenir ce que nous venons de faire pour les sorties du processus `sayHello`, nous aurions plutôt ajouté la ligne suivante à la définition du processus :

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Vous trouverez encore ce modèle de code partout dans les anciens pipelines Nextflow et modules de processus, il est donc important d'en être conscient.
Cependant, nous ne recommandons pas de l'utiliser dans tout nouveau travail car il sera éventuellement interdit dans les futures versions du langage Nextflow.

### À retenir

Vous savez comment publier les sorties du flux de travail dans un emplacement plus pratique.

### Et ensuite ?

Apprendre à fournir une entrée variable via un paramètre de ligne de commande et utiliser efficacement les valeurs par défaut.

---

## 3. Utiliser une entrée variable passée en ligne de commande

Dans son état actuel, notre flux de travail utilise une salutation codée en dur dans la commande du processus.
Nous voulons ajouter de la flexibilité en utilisant une variable d'entrée, afin de pouvoir plus facilement changer la salutation au moment de l'exécution.

Cela nécessite que nous apportions trois ensembles de modifications à notre script :

1. Modifier le processus pour attendre une entrée variable
2. Configurer un paramètre de ligne de commande pour capturer l'entrée utilisateur
3. Passer l'entrée au processus dans le corps du workflow

Faisons ces modifications une à la fois.

### 3.1. Modifier le processus `sayHello` pour attendre une entrée variable

Nous devons modifier la définition du processus pour (1) accepter une variable d'entrée et (2) utiliser cette variable dans la ligne de commande.

#### 3.1.1. Ajouter un bloc input à la définition du processus

D'abord, adaptons la définition du processus pour accepter une entrée appelée `greeting`.

Dans le bloc process, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

La variable `greeting` est préfixée par `val` pour indiquer à Nextflow que c'est une valeur (pas un chemin).

#### 3.1.2. Modifier la commande du processus pour utiliser la variable d'entrée

Maintenant nous échangeons la valeur originale codée en dur pour la valeur de la variable d'entrée que nous attendons de recevoir.

Dans le bloc process, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

Le symbole `$` et les accolades (`{ }`) indiquent à Nextflow que c'est un nom de variable qui doit être remplacé par la valeur d'entrée réelle (=interpolé).

!!! tip "Astuce"

    Les accolades (`{ }`) étaient techniquement optionnelles dans les versions précédentes de Nextflow, vous pourriez donc voir des anciens workflows où ceci est écrit comme `echo '$greeting' > output.txt`.

Maintenant que le processus `sayHello()` est prêt à accepter une entrée variable, nous avons besoin d'un moyen de fournir une valeur d'entrée à l'appel du processus au niveau du workflow.

### 3.2. Configurer un paramètre de ligne de commande pour capturer l'entrée utilisateur

Nous pourrions simplement coder en dur une entrée directement en faisant l'appel de processus `sayHello('Hello World!')`.
Cependant, lorsque nous faisons un vrai travail avec notre flux de travail, nous allons vouloir pouvoir contrôler ses entrées depuis la ligne de commande.

Bonne nouvelle : Nextflow dispose d'un système de paramètres de workflow intégré appelé `params`, qui facilite la déclaration et l'utilisation des paramètres CLI.

La syntaxe générale est de déclarer `params.<nom_paramètre>` pour indiquer à Nextflow d'attendre un paramètre `--<nom_paramètre>` sur la ligne de commande.

Ici, nous voulons créer un paramètre appelé `--input`, donc nous devons déclarer `params.input` quelque part dans le workflow.
En principe, nous pouvons l'écrire n'importe où ; mais puisque nous allons vouloir le donner à l'appel du processus `sayHello()`, nous pouvons le brancher directement là en écrivant `sayHello(params.input)`.

Dans le bloc workflow, faites la modification de code suivante :

=== "Après"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // émettre une salutation
    sayHello(params.input)
    ```

=== "Avant"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // émettre une salutation
    sayHello()
    ```

Cela indique à Nextflow d'exécuter le processus `sayHello` sur la valeur fournie via le paramètre `--input`.

En effet, nous avons accompli les étapes (2) et (3) décrites au début de la section en une seule fois.

### 3.3. Exécuter la commande du flux de travail

Exécutons-le !

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Si vous avez fait toutes ces modifications correctement, vous devriez obtenir une autre exécution réussie.

Assurez-vous d'ouvrir le fichier de sortie pour vérifier que vous avez maintenant la nouvelle version de la salutation.

??? abstract "Contenu du fichier"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà !

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Notez comment la nouvelle exécution a écrasé le fichier de sortie publié dans le répertoire `results`.
Cependant, les résultats des exécutions précédentes sont toujours préservés dans les répertoires de tâches sous `work`.

!!! tip "Astuce"

    Vous pouvez facilement distinguer les paramètres au niveau de Nextflow des paramètres au niveau du pipeline.

    - Les paramètres qui s'appliquent à un pipeline prennent toujours un double tiret (`--`).
    - Les paramètres qui modifient un paramètre Nextflow, _par exemple_ la fonctionnalité `-resume` que nous avons utilisée précédemment, prennent un simple tiret (`-`).

### 3.4. Utiliser des valeurs par défaut pour les paramètres de ligne de commande

D'accord, c'était pratique, mais dans de nombreux cas, il est logique de fournir une valeur par défaut pour un paramètre donné afin de ne pas avoir à le spécifier à chaque exécution.

#### 3.4.1. Définir une valeur par défaut pour le paramètre CLI

Donnons au paramètre `input` une valeur par défaut en le déclarant avant la définition du workflow.

```groovy title="hello-world.nf" linenums="20"
/*
 * Paramètres du pipeline
 */
params {
    input: String = 'Holà mundo!'
}
```

Comme vous pouvez le voir, nous pouvons spécifier le type d'entrée que le workflow attend (Nextflow 25.10.2 et ultérieur).
La syntaxe est `nom: Type = valeur_par_defaut`.
Les types supportés incluent `String`, `Integer`, `Float`, `Boolean` et `Path`.

!!! info "Information"

    Dans les anciens workflows, vous pouvez voir que tout ce bloc `params` est écrit simplement comme `input = 'Holà mundo!'`.

À mesure que vous ajoutez plus de paramètres à votre pipeline, vous devriez tous les ajouter à ce bloc, que vous ayez besoin de leur donner une valeur par défaut ou non.
Cela facilitera la recherche de tous les paramètres configurables en un coup d'œil.

#### 3.4.2. Exécuter à nouveau le flux de travail sans spécifier le paramètre

Maintenant que vous avez une valeur par défaut définie, vous pouvez exécuter à nouveau le flux de travail sans avoir à spécifier une valeur en ligne de commande.

```bash
nextflow run hello-world.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

La sortie sera au même endroit que précédemment, mais le contenu devrait être mis à jour avec le nouveau texte.

??? abstract "Contenu du fichier"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow a utilisé la valeur par défaut du paramètre greeting pour créer la sortie.

#### 3.4.3. Remplacer la valeur par défaut

Si vous fournissez le paramètre en ligne de commande, la valeur CLI remplacera la valeur par défaut.

Essayez-le :

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Encore une fois, vous devriez trouver la sortie mise à jour correspondante dans votre répertoire de résultats.

??? abstract "Contenu du fichier"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Note"

    Dans Nextflow, il y a plusieurs endroits où vous pouvez spécifier des valeurs pour les paramètres.
    Si le même paramètre est défini à des valeurs différentes à plusieurs endroits, Nextflow déterminera quelle valeur utiliser en fonction de l'ordre de priorité qui est décrit [ici](https://www.nextflow.io/docs/latest/config.html).

    Nous couvrirons cela plus en détail dans la Partie 6 (Configuration).

### À retenir

Vous savez comment utiliser une simple entrée variable fournie au moment de l'exécution via un paramètre de ligne de commande, ainsi que comment configurer, utiliser et remplacer les valeurs par défaut.

### Et ensuite ?

Apprendre à gérer les exécutions de flux de travail de manière plus pratique.

---

## 4. Gérer les exécutions de flux de travail

Savoir comment lancer des flux de travail et récupérer les sorties est très bien, mais vous découvrirez rapidement qu'il y a quelques autres aspects de la gestion des flux de travail qui vous faciliteront la vie, surtout si vous développez vos propres flux de travail.

Ici, nous vous montrons comment utiliser la fonctionnalité `resume` lorsque vous devez relancer le même flux de travail, comment inspecter le journal des exécutions passées avec `nextflow log`, et comment supprimer les anciens répertoires de travail avec `nextflow clean`.

<!-- Any other cool options we should include? Added log -->

### 4.1. Relancer un flux de travail avec `-resume`

Parfois, vous allez vouloir relancer un pipeline que vous avez déjà lancé précédemment sans refaire les étapes qui se sont déjà terminées avec succès.

Nextflow a une option appelée `-resume` qui vous permet de faire cela.
Spécifiquement, dans ce mode, tous les processus qui ont déjà été exécutés avec exactement le même code, les mêmes paramètres et les mêmes entrées seront ignorés.
Cela signifie que Nextflow n'exécutera que les processus que vous avez ajoutés ou modifiés depuis la dernière exécution, ou ceux auxquels vous fournissez de nouveaux paramètres ou entrées.

Il y a deux avantages clés à faire cela :

- Si vous êtes en train de développer votre pipeline, vous pouvez itérer plus rapidement puisque vous n'avez qu'à exécuter le(s) processus sur le(s)quel(s) vous travaillez activement pour tester vos modifications.
- Si vous exécutez un pipeline en production et que quelque chose se passe mal, dans de nombreux cas vous pouvez corriger le problème et relancer le pipeline, et il reprendra l'exécution à partir du point d'échec, ce qui peut vous faire gagner beaucoup de temps et de calcul.

Pour l'utiliser, ajoutez simplement `-resume` à votre commande et exécutez-la :

```bash
nextflow run hello-world.nf -resume
```

??? success "Sortie de la commande"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

La sortie console devrait sembler familière, mais il y a une chose qui est un peu différente par rapport à avant.

Cherchez la partie `cached:` qui a été ajoutée dans la ligne d'état du processus (ligne 5), ce qui signifie que Nextflow a reconnu qu'il a déjà fait ce travail et a simplement réutilisé le résultat de l'exécution réussie précédente.

Vous pouvez également voir que le hash du sous-répertoire work est le même que lors de l'exécution précédente.
Nextflow vous indique littéralement l'exécution précédente et dit « J'ai déjà fait cela là-bas. »

!!! tip "Astuce"

    Lorsque vous relancez un pipeline avec `resume`, Nextflow n'écrase pas les fichiers publiés en dehors du répertoire work par les exécutions qui se sont terminées avec succès précédemment.

### 4.2. Inspecter le journal des exécutions passées

Que vous développiez un nouveau pipeline ou que vous exécutiez des pipelines en production, à un moment donné vous aurez probablement besoin de rechercher des informations sur les exécutions passées.
Voici comment faire.

Chaque fois que vous lancez un flux de travail Nextflow, une ligne est écrite dans un fichier journal appelé `history`, sous un répertoire caché appelé `.nextflow` dans le répertoire de travail actuel.

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

Cela affichera le contenu du fichier journal dans le terminal, augmenté d'une ligne d'en-tête.

Vous remarquerez que l'ID de session change chaque fois que vous exécutez une nouvelle commande `nextflow run`, SAUF si vous utilisez l'option `-resume`.
Dans ce cas, l'ID de session reste le même.

Nextflow utilise l'ID de session pour regrouper les informations de mise en cache d'exécution sous le répertoire `cache`, également situé sous `.nextflow`.

### 4.3. Supprimer les anciens répertoires de travail

Pendant le processus de développement, vous exécuterez généralement votre pipeline en cours d'élaboration un grand nombre de fois, ce qui peut conduire à une accumulation de nombreux fichiers dans de nombreux sous-répertoires.

Heureusement, Nextflow inclut une sous-commande `clean` utile qui peut automatiquement supprimer les sous-répertoires de travail des exécutions passées dont vous ne vous souciez plus.

#### 4.3.1. Déterminer les critères de suppression

Il existe plusieurs [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) pour déterminer ce qu'il faut supprimer.

Ici, nous vous montrons un exemple qui supprime tous les sous-répertoires des exécutions avant une exécution donnée, spécifiée en utilisant son nom d'exécution.

Recherchez l'exécution réussie la plus récente où vous n'avez pas utilisé `-resume` ; dans notre cas, le nom d'exécution était `golden_cantor`.

Le nom d'exécution est la chaîne en deux parties générée par la machine affichée entre crochets dans la ligne de sortie console `Launching (...)`.
Vous pouvez également utiliser le journal Nextflow pour rechercher une exécution en fonction de son horodatage et/ou de sa ligne de commande.

#### 4.3.2. Faire un essai à blanc

D'abord, nous utilisons le drapeau d'essai à blanc `-n` pour vérifier ce qui sera supprimé avec la commande :

```bash
nextflow clean -before golden_cantor -n
```

??? success "Sortie de la commande"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Votre sortie aura des noms de répertoires de tâches différents et peut avoir un nombre de lignes différent, mais elle devrait ressembler à l'exemple.

Si vous ne voyez aucune ligne en sortie, soit vous n'avez pas fourni un nom d'exécution valide, soit il n'y a pas d'exécutions passées à supprimer. Assurez-vous de changer `golden_cantor` dans la commande d'exemple par le nom d'exécution le plus récent correspondant dans votre journal.

#### 4.3.3. Procéder à la suppression

Si la sortie semble comme prévu et que vous voulez procéder à la suppression, relancez la commande avec le drapeau `-f` au lieu de `-n` :

```bash
nextflow clean -before golden_cantor -f
```

??? success "Sortie de la commande"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

La sortie devrait être similaire à avant, mais disant maintenant « Removed » au lieu de « Would remove ».
Notez que cela ne supprime pas les sous-répertoires à deux caractères (comme `a3/` ci-dessus) mais vide leur contenu.

!!! warning "Avertissement"

    La suppression des sous-répertoires de travail des exécutions passées les supprime du cache de Nextflow et supprime toutes les sorties qui étaient stockées dans ces répertoires.
    Cela signifie que cela brise la capacité de Nextflow à reprendre l'exécution sans relancer les processus correspondants.

    Vous êtes responsable de sauvegarder toutes les sorties qui vous importent ou sur lesquelles vous prévoyez de compter ! C'est la principale raison pour laquelle nous préférons utiliser le mode `copy` plutôt que le mode `symlink` pour la directive `publish`.

### À retenir

Vous savez comment publier les sorties dans un répertoire spécifique, relancer un pipeline sans répéter les étapes qui ont déjà été exécutées de manière identique, et utiliser la commande `nextflow clean` pour nettoyer les anciens répertoires de travail.

Plus généralement, vous savez comment interpréter un flux de travail Nextflow simple, gérer son exécution et récupérer les sorties.

### Et ensuite ?

Faites une petite pause, vous l'avez bien mérité !

Quand vous êtes prêt·e, passez à la [**Partie 2 : Hello Channels**](./02_hello_channels.md) pour apprendre comment utiliser les canaux pour alimenter les entrées dans votre flux de travail, ce qui vous permettra de profiter du parallélisme de flux de données intégré de Nextflow et d'autres fonctionnalités puissantes.

---

## Quiz

<quiz>
Quels sont les composants minimaux requis d'un processus Nextflow ?
- [ ] Blocs input et output uniquement
- [x] Blocs output et script
- [ ] Blocs input, output et script
- [ ] Uniquement un bloc script

En savoir plus : [1.1.1. La définition du processus](#111-la-definition-du-process)
</quiz>

<quiz>
Quel est le but du bloc output dans un processus ?
- [ ] Imprimer les résultats dans la console
- [ ] Enregistrer les fichiers dans le répertoire work
- [x] Déclarer les sorties attendues du processus
- [ ] Définir les variables d'environnement

En savoir plus : [1.1.1. La définition du processus](#111-la-definition-du-process)
</quiz>

<quiz>
Quelle commande est utilisée pour exécuter un flux de travail Nextflow ?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
En regardant le répertoire work d'une tâche, quel fichier contient la commande réellement exécutée ?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

En savoir plus : [1.2.2. Trouver la sortie et les logs dans le répertoire `work`](#122-trouver-la-sortie-et-les-logs-dans-le-repertoire-work)
</quiz>

<quiz>
Que fait le drapeau `-resume` ?
- [ ] Redémarre le flux de travail depuis le début
- [ ] Met le flux de travail en pause
- [x] Ignore les processus qui se sont déjà terminés avec succès
- [ ] Crée une sauvegarde du flux de travail

En savoir plus : [4.1. Relancer un flux de travail avec `-resume`](#41-relancer-un-flux-de-travail-avec--resume)
</quiz>

<quiz>
Quel est le mode par défaut pour publier les sorties du flux de travail ?
- [ ] Copier les fichiers dans le répertoire de sortie
- [x] Créer des liens symboliques dans le répertoire de sortie
- [ ] Déplacer les fichiers dans le répertoire de sortie
- [ ] Compresser les fichiers dans le répertoire de sortie

En savoir plus : [2.3. Définir le mode de publication sur copy](#23-definir-le-mode-de-publication-sur-copy)
</quiz>

<quiz>
Comment passez-vous une valeur de paramètre à un flux de travail Nextflow depuis la ligne de commande ?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

En savoir plus : [3.2. Configurer un paramètre de ligne de commande pour capturer l'entrée utilisateur](#32-configurer-un-parametre-de-ligne-de-commande-pour-capturer-lentree-utilisateur)
</quiz>

<quiz>
Comment référencez-vous une variable à l'intérieur d'un bloc script Nextflow ?
- [ ] Utiliser la syntaxe `%variable%`
- [x] Utiliser la syntaxe `#!groovy ${variable}`
- [ ] Utiliser la syntaxe `{{variable}}`
- [ ] Utiliser la syntaxe `[variable]`
</quiz>
