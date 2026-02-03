# Partie 4 : Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Voir [la playlist complète](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) sur la chaîne YouTube de Nextflow.

:green_book: La transcription de la vidéo est disponible [ici](./transcripts/04_hello_modules.md).
///
-->

Cette section couvre comment organiser votre code de workflow pour rendre le développement et la maintenance de votre pipeline plus efficaces et durables.
Plus précisément, nous allons démontrer comment utiliser des **modules**.

Dans Nextflow, un **module** est une définition de processus unique qui est encapsulée par elle-même dans un fichier de code autonome.
Pour utiliser un module dans un workflow, vous ajoutez simplement une instruction d'importation sur une seule ligne à votre fichier de code de workflow ; ensuite vous pouvez intégrer le processus dans le workflow de la même manière que vous le feriez normalement.
Cela permet de réutiliser des définitions de processus dans plusieurs workflows sans produire plusieurs copies du code.

Quand nous avons commencé à développer notre workflow, nous avons tout écrit dans un seul fichier de code.
Maintenant, nous allons déplacer les processus dans des modules individuels.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Cela rendra notre code plus partageable, flexible et maintenable.

??? info "Comment commencer à partir de cette section"

    Cette section du cours suppose que vous avez complété les Parties 1-3 du cours [Hello Nextflow](./index.md), mais si vous êtes à l'aise avec les bases couvertes dans ces sections, vous pouvez commencer ici sans rien faire de spécial.

---

## 0. Échauffement : Exécuter `hello-modules.nf`

Nous allons utiliser le script de workflow `hello-modules.nf` comme point de départ.
Il est équivalent au script produit en travaillant à travers la Partie 3 de ce cours de formation, sauf que nous avons changé les destinations de sortie :

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Juste pour s'assurer que tout fonctionne, exécutez le script une fois avant d'effectuer des modifications :

```bash
nextflow run hello-modules.nf
```

??? success "Sortie de la commande"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Comme précédemment, vous trouverez les fichiers de sortie dans le répertoire spécifié dans le bloc `output` (ici, `results/hello_modules/`).

??? abstract "Contenu du répertoire"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Si cela a fonctionné pour vous, vous êtes prêt·e à apprendre comment modulariser votre code de workflow.

---

## 1. Créer un répertoire pour stocker les modules

Il est recommandé de stocker vos modules dans un répertoire spécifique.
Vous pouvez appeler ce répertoire comme vous le souhaitez, mais la convention est de l'appeler `modules/`.

```bash
mkdir modules
```

!!! tip "Astuce"

    Ici, nous vous montrons comment utiliser des **modules locaux**, c'est-à-dire des modules stockés localement dans le même dépôt que le reste du code du workflow, par opposition aux modules distants, qui sont stockés dans d'autres dépôts (distants).
    Pour plus d'informations sur les **modules distants**, consultez la [documentation](https://www.nextflow.io/docs/latest/module.html).

---

## 2. Créer un module pour `sayHello()`

Dans sa forme la plus simple, transformer un processus existant en module n'est guère plus qu'une opération de copier-coller.
Nous allons créer un fichier stub pour le module, copier le code pertinent puis le supprimer du fichier de workflow principal.

Ensuite, tout ce que nous aurons à faire est d'ajouter une instruction d'importation pour que Nextflow sache récupérer le code pertinent à l'exécution.

### 2.1. Créer un fichier stub pour le nouveau module

Créons un fichier vide pour le module appelé `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

Cela nous donne un endroit pour mettre le code du processus.

### 2.2. Déplacer le code du processus `sayHello` vers le fichier module

Copiez l'ensemble de la définition du processus du fichier de workflow vers le fichier module, en vous assurant de copier également le shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/sayHello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Utilise echo pour imprimer 'Hello World!' dans un fichier
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

Une fois cela fait, supprimez la définition du processus du fichier de workflow, mais assurez-vous de laisser le shebang en place.

### 2.3. Ajouter une déclaration d'importation avant le bloc workflow

La syntaxe pour importer un module local est assez simple :

```groovy title="Syntaxe : Déclaration d'importation"
include { <NOM_DU_MODULE> } from '<chemin_vers_le_module>'
```

Insérons cela au-dessus du bloc `params` et remplissons-le de manière appropriée.

=== "Après"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Paramètres du pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Avant"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Paramètres du pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Vous voyez que nous avons renseigné le nom du module, `sayHello`, et le chemin vers le fichier contenant le code du module, `./modules/sayHello.nf`.

### 2.4. Exécuter le workflow

Nous exécutons le workflow avec essentiellement le même code et les mêmes entrées qu'avant, donc exécutons avec le flag `-resume` et voyons ce qui se passe.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Cela devrait s'exécuter très rapidement car tout est en cache.
N'hésitez pas à vérifier les sorties publiées.

Nextflow a reconnu que c'est toujours le même travail à faire, même si le code est divisé en plusieurs fichiers.

### À retenir

Vous savez comment extraire un processus dans un module local et vous savez que faire cela ne casse pas la reprise du workflow.

### Et ensuite ?

Pratiquez la création de plus de modules.
Une fois que vous en avez fait un, vous pouvez en faire un million de plus...
Mais faisons-en juste deux de plus pour l'instant.

---

## 3. Modulariser le processus `convertToUpper()`

### 3.1. Créer un fichier stub pour le nouveau module

Créez un fichier vide pour le module appelé `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Déplacer le code du processus `convertToUpper` vers le fichier module

Copiez l'ensemble de la définition du processus du fichier de workflow vers le fichier module, en vous assurant de copier également le shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/convertToUpper.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Utilise un outil de remplacement de texte pour convertir la salutation en majuscules
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
```

Une fois cela fait, supprimez la définition du processus du fichier de workflow, mais assurez-vous de laisser le shebang en place.

### 3.3. Ajouter une déclaration d'importation avant le bloc `params`

Insérez la déclaration d'importation au-dessus du bloc `params` et remplissez-la de manière appropriée.

=== "Après"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Paramètres du pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Avant"

    ```groovy title="hello-modules.nf" linenums="23"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Paramètres du pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Cela devrait commencer à sembler très familier.

### 3.4. Exécuter à nouveau le workflow

Exécutez ceci avec le flag `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Cela devrait toujours produire la même sortie que précédemment.

Deux de faits, plus qu'un !

---

## 4. Modulariser le processus `collectGreetings()`

### 4.1. Créer un fichier stub pour le nouveau module

Créez un fichier vide pour le module appelé `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Déplacer le code du processus `collectGreetings` vers le fichier module

Copiez l'ensemble de la définition du processus du fichier de workflow vers le fichier module, en vous assurant de copier également le shebang `#!/usr/bin/env nextflow`.

```groovy title="modules/collectGreetings.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Collecter les salutations en majuscules dans un seul fichier de sortie
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
```

Une fois cela fait, supprimez la définition du processus du fichier de workflow, mais assurez-vous de laisser le shebang en place.

### 4.3. Ajouter une déclaration d'importation avant le bloc `params`

Insérez la déclaration d'importation au-dessus du bloc `params` et remplissez-la de manière appropriée.

=== "Après"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Paramètres du pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Avant"

    ```groovy title="hello-modules.nf" linenums="3"
    // Inclure les modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Paramètres du pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Dernier !

### 4.4. Exécuter le workflow

Exécutez ceci avec le flag `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Cela devrait toujours produire la même sortie que précédemment.

### À retenir

Vous savez comment modulariser plusieurs processus dans un workflow.

Félicitations, vous avez fait tout ce travail et absolument rien n'a changé dans le fonctionnement du pipeline !

Blague à part, maintenant votre code est plus modulaire, et si vous décidez d'écrire un autre pipeline qui fait appel à l'un de ces processus, vous n'avez qu'à taper une courte instruction d'importation pour utiliser le module pertinent.
C'est mieux que de copier-coller le code, parce que si plus tard vous décidez d'améliorer le module, tous vos pipelines hériteront des améliorations.

### Et ensuite ?

Prenez une courte pause si vous le souhaitez.

Quand vous êtes prêt·e, passez à la [**Partie 5 : Hello Containers**](./05_hello_containers.md) pour apprendre comment utiliser les conteneurs pour gérer les dépendances logicielles de manière plus pratique et reproductible.

---

## Quiz

<quiz>
Qu'est-ce qu'un module dans Nextflow ?
- [ ] Un fichier de configuration
- [x] Un fichier autonome contenant une définition de processus unique
- [ ] Une définition de workflow
- [ ] Un opérateur de canal

En savoir plus : [2. Créer un module pour `sayHello()`](#2-creer-un-module-pour-sayhello)
</quiz>

<quiz>
Quelle est la convention de nommage recommandée pour les fichiers de module ?
- [ ] `module_processName.nf`
- [ ] `processName_module.nf`
- [x] `processName.nf`
- [ ] `mod_processName.nf`
</quiz>

<quiz>
Où les fichiers de module doivent-ils être stockés ?
- [ ] Dans le même répertoire que le workflow
- [ ] Dans un répertoire `bin/`
- [x] Dans un répertoire `modules/`
- [ ] Dans un répertoire `lib/`

En savoir plus : [1. Créer un répertoire pour stocker les modules](#1-creer-un-repertoire-pour-stocker-les-modules)
</quiz>

<quiz>
Quelle est la syntaxe correcte pour importer un module ?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

En savoir plus : [2.3. Ajouter une déclaration d'importation](#23-ajouter-une-declaration-dimportation-avant-le-bloc-workflow)
</quiz>

<quiz>
Que se passe-t-il avec la fonctionnalité `-resume` lors de l'utilisation de modules ?
- [ ] Elle ne fonctionne plus
- [ ] Elle nécessite une configuration supplémentaire
- [x] Elle fonctionne comme avant
- [ ] Elle ne fonctionne que pour les modules locaux
</quiz>

<quiz>
Quels sont les avantages de l'utilisation de modules ? (Sélectionnez toutes les réponses applicables)
- [x] Réutilisabilité du code entre les workflows
- [x] Maintenance plus facile
- [x] Meilleure organisation du code du workflow
- [ ] Vitesse d'exécution plus rapide
</quiz>
