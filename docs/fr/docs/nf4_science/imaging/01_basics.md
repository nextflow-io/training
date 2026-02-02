# Partie 1 : Exécuter des opérations de base

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dans cette première partie du cours de formation Nextflow pour la bio-imagerie, nous allons utiliser un exemple Hello World très basique et indépendant du domaine pour démontrer les opérations essentielles et identifier les composants de code Nextflow correspondants.

## 1. Exécuter le workflow

Nous vous fournissons un script de workflow nommé `hello-world.nf` qui prend une entrée via un argument de ligne de commande nommé `--greeting` et produit un fichier texte contenant ce message de salutation.
Nous n'allons pas encore examiner le code ; voyons d'abord à quoi ressemble son exécution.

### 1.1. Lancer le workflow et surveiller l'exécution

Dans le terminal, exécutez la commande suivante :

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

La sortie de votre console devrait ressembler à ceci :

```console title="Sortie" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Félicitations, vous venez d'exécuter votre premier workflow Nextflow !

La sortie la plus importante ici est la dernière ligne (ligne 6) :

```console title="Sortie" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Cela nous indique que le processus `sayHello` a été exécuté avec succès une fois (`1 of 1 ✔`).

C'est parfait, mais vous vous demandez peut-être : où se trouve la sortie ?

### 1.2. Trouver le fichier de sortie dans le répertoire `results`

Ce workflow est configuré pour publier sa sortie dans un répertoire appelé `results`.
Si vous regardez votre répertoire actuel, vous verrez que lorsque vous avez exécuté le workflow, Nextflow a créé un nouveau répertoire appelé `results`, qui contient un fichier appelé `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Ouvrez le fichier ; le contenu devrait correspondre au message de salutation que vous avez spécifié sur la ligne de commande.

<details>
  <summary>Contenu du fichier</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

C'est parfait, notre workflow a fait ce qu'il était censé faire !

Cependant, sachez que le résultat 'publié' est une copie (ou dans certains cas un lien symbolique) de la sortie réelle produite par Nextflow lors de l'exécution du workflow.

Maintenant, nous allons regarder sous le capot pour voir où Nextflow a réellement effectué le travail.

!!! warning "Avertissement"

    Tous les workflows ne seront pas configurés pour publier les sorties dans un répertoire results, et/ou le nom du répertoire peut être différent.
    Un peu plus loin dans cette section, nous vous montrerons comment découvrir où ce comportement est spécifié.

### 1.3. Trouver la sortie originale et les logs dans le répertoire `work/`

Lorsque vous exécutez un workflow, Nextflow crée un 'répertoire de tâche' distinct pour chaque invocation de chaque processus dans le workflow (=chaque étape du pipeline).
Pour chacun, il va préparer les entrées nécessaires, exécuter la ou les instruction(s) pertinente(s) et écrire les sorties et les fichiers de log dans ce répertoire unique, qui est nommé automatiquement à l'aide d'un hash afin de le rendre unique.

Tous ces répertoires de tâches se trouveront dans un répertoire appelé `work` dans votre répertoire actuel (où vous exécutez la commande).

Cela peut sembler déroutant, alors voyons à quoi cela ressemble en pratique.

En revenant à la sortie console du workflow que nous avons exécuté précédemment, nous avions cette ligne :

```console title="Extrait de la sortie de la commande" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Voyez comment la ligne commence par `[a3/7be2fa]` ?
C'est une forme tronquée du chemin du répertoire de tâche pour cet appel de processus particulier, et vous indique où trouver la sortie de l'appel du processus `sayHello` dans le chemin du répertoire `work/`.

Vous pouvez trouver le chemin complet en tapant la commande suivante (en remplaçant `a3/7be2fa` par ce que vous voyez dans votre propre terminal) et en appuyant sur la touche tab pour compléter automatiquement le chemin ou en ajoutant un astérisque :

```bash
tree work/a3/7be2fa*
```

Cela devrait donner le chemin complet du répertoire : `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Voyons ce qu'il y a dedans.

!!! Tip "Astuce"

    Si vous parcourez le contenu du sous-répertoire de tâche dans l'explorateur de fichiers de VSCode, vous verrez tous les fichiers immédiatement.
    Cependant, les fichiers de log sont configurés pour être invisibles dans le terminal, donc si vous voulez utiliser `ls` ou `tree` pour les voir, vous devrez définir l'option pertinente pour afficher les fichiers invisibles.

    ```bash
    tree -a work
    ```

Les noms exacts des sous-répertoires seront différents sur votre système.

<details>
  <summary>Contenu du répertoire</summary>

```console title="work/"
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

</details>

Vous devriez immédiatement reconnaître le fichier `output.txt`, qui est en fait la sortie originale du processus `sayHello` qui a été publiée dans le répertoire `results`.
Si vous l'ouvrez, vous retrouverez le message de salutation `Hello World!`.

<details>
  <summary>Contenu du fichier output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

Alors, qu'en est-il de tous ces autres fichiers ?

Ce sont les fichiers d'aide et de log que Nextflow a écrits dans le cadre de l'exécution de la tâche :

- **`.command.begin`** : Fichier sentinelle créé dès que la tâche est lancée.
- **`.command.err`** : Messages d'erreur (`stderr`) émis par l'appel du processus
- **`.command.log`** : Sortie de log complète émise par l'appel du processus
- **`.command.out`** : Sortie normale (`stdout`) par l'appel du processus
- **`.command.run`** : Script complet exécuté par Nextflow pour exécuter l'appel du processus
- **`.command.sh`** : La commande qui a réellement été exécutée par l'appel du processus
- **`.exitcode`** : Le code de sortie résultant de la commande

Le fichier `.command.sh` est particulièrement utile car il vous montre la commande principale que Nextflow a exécutée, sans inclure toute la gestion administrative et la configuration de la tâche/environnement.

<details>
  <summary>Contenu du fichier</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Astuce"

    Lorsque quelque chose ne va pas et que vous devez résoudre le problème, il peut être utile de regarder le script `command.sh` pour vérifier exactement quelle commande Nextflow a composée en fonction des instructions du workflow, de l'interpolation des variables, etc.

### 1.4. Exercice optionnel : ré-exécuter avec différents messages de salutation

Essayez de ré-exécuter le workflow plusieurs fois avec différentes valeurs pour l'argument `--greeting`, puis regardez à la fois le contenu du répertoire `results/` et les répertoires de tâches.

Observez comment les sorties et les logs des répertoires de tâches isolés sont préservés, alors que le contenu du répertoire `results` est écrasé par la sortie des exécutions suivantes.

### À retenir

Vous savez comment exécuter un script Nextflow simple, surveiller son exécution et trouver ses sorties.

### Et ensuite ?

Apprenez à lire un script Nextflow de base et à identifier comment ses composants sont liés à ses fonctionnalités.

---

## 2. Examiner le script de démarrage du workflow Hello World

Ce que nous avons fait là, c'était essentiellement traiter le script de workflow comme une boîte noire.
Maintenant que nous avons vu ce qu'il fait, ouvrons la boîte et regardons à l'intérieur.

_L'objectif ici n'est pas de mémoriser la syntaxe du code Nextflow, mais de former une intuition de base de quels sont les composants principaux et comment ils sont organisés._

### 2.1. Examiner la structure générale du code

Ouvrons le script `hello-world.nf` dans le volet de l'éditeur.

<details>
  <summary>Code</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Utilise echo pour imprimer une salutation dans un fichier
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // émettre une salutation
    sayHello(params.greeting)
}
```

</details>

Un script Nextflow implique deux types principaux de composants de base : un ou plusieurs **processus**, et le **workflow** lui-même.
Chaque **processus** décrit quelle(s) opération(s) l'étape correspondante du pipeline doit accomplir, tandis que le **workflow** décrit la logique de flux de données qui connecte les différentes étapes.

Examinons d'abord de plus près le bloc **process**, puis nous examinerons le bloc **workflow**.

### 2.2. La définition du `process`

Le premier bloc de code décrit un **processus**.
La définition du processus commence par le mot-clé `process`, suivi du nom du processus et enfin du corps du processus délimité par des accolades.
Le corps du processus doit contenir un bloc script qui spécifie la commande à exécuter, qui peut être n'importe quoi que vous seriez capable d'exécuter dans un terminal de ligne de commande.

Ici, nous avons un **processus** appelé `sayHello` qui prend une variable d'**entrée** appelée `greeting` et écrit sa **sortie** dans un fichier nommé `output.txt`.

<details>
  <summary>Code</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Utilise echo pour imprimer une salutation dans un fichier
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

C'est une définition de processus très minimale qui contient juste une définition d'`input`, une définition d'`output` et le `script` à exécuter.

La définition d'`input` inclut le qualificateur `val`, qui indique à Nextflow de s'attendre à une valeur d'un certain type (peut être une chaîne, un nombre, peu importe).

La définition d'`output` inclut le qualificateur `path`, qui indique à Nextflow que cela doit être traité comme un chemin (inclut à la fois les chemins de répertoires et les fichiers).

!!! Tip "Astuce"

    La définition de sortie ne _détermine_ pas quelle sortie sera créée.
    Elle _déclare_ simplement où trouver le(s) fichier(s) de sortie attendu(s), afin que Nextflow puisse le chercher une fois l'exécution terminée.

    Ceci est nécessaire pour vérifier que la commande a été exécutée avec succès et pour transmettre la sortie aux processus en aval si nécessaire.
    Les sorties produites qui ne correspondent pas à ce qui est déclaré dans le bloc de sortie ne seront pas transmises aux processus en aval.

Dans un pipeline réel, un processus contient généralement des informations supplémentaires telles que des directives de processus, que nous présenterons dans un instant.

### 2.3. La définition du `workflow`

Le deuxième bloc de code décrit le **workflow** lui-même.
La définition du workflow commence par le mot-clé `workflow`, suivi d'un nom optionnel, puis du corps du workflow délimité par des accolades.

Ici, nous avons un **workflow** qui consiste en un appel au processus `sayHello`, qui prend une entrée, `params.greeting`, qui contient la valeur que nous avons donnée au paramètre `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // émettre une salutation
    sayHello(params.greeting)
}
```

C'est une définition de **workflow** très minimale.
Dans un pipeline réel, le workflow contient généralement plusieurs appels à des **processus** connectés par des **canaux**, et il peut y avoir des valeurs par défaut configurées pour les entrées variables.

Nous verrons cela en action lorsque nous exécuterons nf-core/molkart dans la Partie 2 du cours.

### 2.4. Le système `params` de paramètres de ligne de commande

Le `params.greeting` que nous fournissons à l'appel du processus `sayHello()` est un élément intéressant du code Nextflow et mérite qu'on y consacre une minute supplémentaire.

Comme mentionné ci-dessus, c'est ainsi que nous transmettons la valeur du paramètre de ligne de commande `--greeting` à l'appel du processus `sayHello()`.
En fait, le simple fait de déclarer `params.someParameterName` nous permettra de donner au workflow un paramètre nommé `--someParameterName` depuis la ligne de commande.

!!! Tip "Astuce"

    Ces paramètres de workflow déclarés à l'aide du système `params` prennent toujours deux tirets (`--`).
    Cela les distingue des paramètres au niveau de Nextflow, qui ne prennent qu'un seul tiret (`-`).

### À retenir

Vous savez maintenant comment un workflow Nextflow simple est structuré, et comment les composants de base sont liés à ses fonctionnalités.

### Et ensuite ?

Apprenez à gérer vos exécutions de workflow de manière pratique.

---

## 3. Gérer les exécutions de workflow

Savoir comment lancer des workflows et récupérer les sorties est excellent, mais vous découvrirez rapidement qu'il existe quelques autres aspects de la gestion des workflows qui vous faciliteront la vie.

Nous vous montrons ici comment profiter de la fonctionnalité `resume` lorsque vous devez relancer le même workflow, comment inspecter les logs d'exécution avec `nextflow log`, et comment supprimer les anciens répertoires work avec `nextflow clean`.

### 3.1. Relancer un workflow avec `-resume`

Parfois, vous voudrez ré-exécuter un pipeline que vous avez déjà lancé précédemment sans refaire le travail qui a déjà été effectué avec succès.

Nextflow a une option appelée `-resume` qui vous permet de faire cela.
Plus précisément, dans ce mode, tous les processus qui ont déjà été exécutés avec exactement le même code, les mêmes paramètres et les mêmes entrées seront ignorés.
Cela signifie que Nextflow n'exécutera que les processus que vous avez ajoutés ou modifiés depuis la dernière exécution, ou auxquels vous fournissez de nouveaux paramètres ou entrées.

Il y a deux avantages clés à faire cela :

- Si vous êtes en train de développer un pipeline, vous pouvez itérer plus rapidement puisque vous n'avez qu'à exécuter le(s) processus sur le(s)quel(s) vous travaillez activement pour tester vos modifications.
- Si vous exécutez un pipeline en production et que quelque chose ne va pas, dans de nombreux cas vous pouvez corriger le problème et relancer le pipeline, et il reprendra l'exécution à partir du point de défaillance, ce qui peut vous faire gagner beaucoup de temps et de calcul.

Pour l'utiliser, ajoutez simplement `-resume` à votre commande et exécutez-la :

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Recherchez le bit `cached:` qui a été ajouté dans la ligne de statut du processus (ligne 5), ce qui signifie que Nextflow a reconnu qu'il a déjà fait ce travail et a simplement réutilisé le résultat de l'exécution précédente réussie.

Vous pouvez également voir que le hash du sous-répertoire work est le même que lors de l'exécution précédente.
Nextflow vous indique littéralement l'exécution précédente et dit « J'ai déjà fait ça là-bas. »

!!! Tip "Astuce"

    Lorsque vous ré-exécutez un pipeline avec `resume`, Nextflow n'écrase aucun fichier écrit dans un répertoire `publishDir` par un appel de processus qui a été précédemment exécuté avec succès.

### 3.2. Inspecter le log des exécutions passées

Chaque fois que vous lancez un workflow nextflow, une ligne est écrite dans un fichier de log appelé `history`, dans un répertoire caché appelé `.nextflow` dans le répertoire de travail actuel.

Une façon plus pratique d'accéder à cette information est d'utiliser la commande `nextflow log`.

```bash
nextflow log
```

Cela affichera le contenu du fichier de log dans le terminal, vous montrant l'horodatage, le nom de l'exécution, le statut et la ligne de commande complète pour chaque exécution Nextflow qui a été lancée depuis le répertoire de travail actuel.

### 3.3. Supprimer les anciens répertoires work

Pendant le processus de développement, vous exécuterez généralement vos projets de pipelines un grand nombre de fois, ce qui peut conduire à une accumulation de très nombreux fichiers dans de nombreux sous-répertoires.
Comme les sous-répertoires sont nommés aléatoirement, il est difficile de dire d'après leurs noms lesquels sont des exécutions plus anciennes ou plus récentes.

Nextflow inclut une sous-commande `clean` pratique qui peut automatiquement supprimer les sous-répertoires work pour les exécutions passées qui ne vous intéressent plus, avec plusieurs [options](https://www.nextflow.io/docs/latest/reference/cli.html#clean) pour contrôler ce qui sera supprimé.

Vous pouvez utiliser le log Nextflow pour rechercher une exécution en fonction de son horodatage et/ou de sa ligne de commande, puis utiliser `nextflow clean -before <run_name> -f` pour supprimer les répertoires work des exécutions antérieures.

!!! Warning "Avertissement"

    La suppression des sous-répertoires work des exécutions passées les supprime du cache de Nextflow et supprime toutes les sorties qui étaient stockées dans ces répertoires.
    Cela signifie que cela brise la capacité de Nextflow à reprendre l'exécution sans ré-exécuter les processus correspondants.

    Vous êtes responsable de la sauvegarde de toutes les sorties qui vous intéressent ou sur lesquelles vous prévoyez de vous appuyer ! Si vous utilisez la directive `publishDir` à cette fin, assurez-vous d'utiliser le mode `copy`, pas le mode `symlink`.

### À retenir

Vous savez comment relancer un pipeline sans répéter les étapes qui ont déjà été exécutées de manière identique, inspecter le log d'exécution, et utiliser la commande `nextflow clean` pour nettoyer les anciens répertoires work.

### Et ensuite ?

Maintenant que vous comprenez les opérations de base de Nextflow, vous êtes prêt·e à exécuter un véritable pipeline de bio-imagerie avec nf-core/molkart.
