# Workflows de Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Lorsque vous développez un pipeline, vous vous retrouvez souvent à créer des séquences similaires de processus pour différents types de données ou étapes d'analyse.
Vous pourriez finir par copier-coller ces séquences de processus, ce qui entraîne du code dupliqué difficile à maintenir ; ou vous pourriez créer un workflow monolithique difficile à comprendre et à modifier.

L'une des fonctionnalités les plus puissantes de Nextflow est sa capacité à composer des pipelines complexes à partir de modules de workflow plus petits et réutilisables.
Cette approche modulaire rend les pipelines plus faciles à développer, tester et maintenir.

### Objectifs d'apprentissage

Dans cette quête secondaire, nous allons explorer comment développer des modules de workflow pouvant être testés et utilisés séparément, composer ces modules en un pipeline plus grand, et gérer le flux de données entre les modules.

À la fin de cette quête secondaire, vous serez en mesure de :

- Décomposer des pipelines complexes en unités logiques et réutilisables
- Tester chaque module de workflow indépendamment
- Combiner des workflows pour créer de nouveaux pipelines
- Partager des modules de workflow communs entre différents pipelines
- Rendre votre code plus maintenable et plus facile à comprendre

Ces compétences vous aideront à construire des pipelines complexes tout en maintenant une structure de code propre et maintenable.

### Prérequis

Avant de vous lancer dans cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../../hello_nextflow/index.md) ou un cours équivalent pour débutant·es.
- Être à l'aise avec les concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs, modules)

---

## 0. Premiers pas

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplacez-vous dans le répertoire où se trouvent les fichiers de ce tutoriel.

```bash
cd side-quests/workflows_of_workflows
```

Vous pouvez configurer VSCode pour qu'il se concentre sur ce répertoire :

```bash
code .
```

L'éditeur s'ouvre avec le répertoire du projet en focus.

#### Examiner les fichiers

Vous trouverez un répertoire `modules` contenant des définitions de processus, un répertoire `workflows` contenant deux scripts de workflow pré-écrits, et un fichier `main.nf` que vous mettrez à jour progressivement :

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

Le répertoire `modules/` contient les définitions de processus individuelles, et le répertoire `workflows/` contient les deux scripts de workflow pré-écrits avec lesquels vous travaillerez dans cette quête secondaire.

#### Examiner l'exercice

Votre défi consiste à assembler ces modules en deux workflows distincts que nous composerons ensuite en un workflow principal :

- Un `GREETING_WORKFLOW` qui valide les noms, crée des salutations et ajoute des horodatages
- Un `TRANSFORM_WORKFLOW` qui convertit le texte en majuscules et l'inverse

Le pipeline final enchaîne les deux workflows en un seul flux de données :

1. **Valider** : vérifier que chaque nom est bien formé
2. **Saluer** : générer une salutation pour chaque nom valide
3. **Horodater** : enregistrer quand chaque salutation a été créée
4. **Majuscules** : convertir la salutation horodatée en majuscules
5. **Inverser** : inverser le texte en majuscules

Les trois premières étapes appartiennent au `GREETING_WORKFLOW`, et les deux dernières au `TRANSFORM_WORKFLOW`.
Les construire comme des workflows séparés et composables vous permet de développer et tester chaque étape indépendamment avant de les assembler.

#### Liste de vérification

Vous pensez être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Ajouter le greeting workflow au pipeline

Le greeting workflow valide les noms et génère des salutations horodatées.

### 1.1. Examiner et exécuter le greeting workflow

Ouvrez `workflows/greeting.nf` et examinez le code :

```groovy title="workflows/greeting.nf" linenums="1"
#!/usr/bin/env nextflow

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chaîne de processus : valider -> créer la salutation -> ajouter l'horodatage
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Il s'agit d'un workflow complet et autonome, avec la même structure que celle que vous avez vue dans le tutoriel 'Hello Nextflow'.
Il code en dur les noms d'entrée, enchaîne trois processus et publie deux sorties.

Exécutez-le pour vérifier que tout fonctionne :

```bash
nextflow run workflows/greeting.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `workflows/greeting.nf` [sharp_kirch] revision: 04c0a33f79
    executor >  local (9)
    [c5/9df3c8] VAL…TE_NAME (validating Alice) | 3 of 3 ✔
    [e0/243874] SAY_HELLO (greeting Alice)     | 3 of 3 ✔
    [c5/e40c8d] TIM…ing timestamp to greeting) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      greetings:
        - Charlie-output.txt
        - Bob-output.txt
        - Alice-output.txt

      timestamped:
        - timestamped_Bob-output.txt
        - timestamped_Charlie-output.txt
        - timestamped_Alice-output.txt
    ```

Pour le rendre composable avec d'autres workflows, quelques modifications sont nécessaires.

### 1.2. Rendre le workflow composable

Pour rendre un workflow composable, trois choses doivent changer :
le workflow reçoit un nom, les entrées sont déplacées dans un bloc `take:`, et les sorties sont déplacées dans un bloc `emit:`
(en remplacement des blocs autonomes `publish:`/`output {}`, qui appartiennent au entry workflow).

Les sections suivantes couvrent ces modifications une par une.

#### 1.2.1. Nommer le workflow

Donnez un nom au workflow pour qu'il puisse être importé depuis un workflow parent.

=== "Après"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Avant"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Avec un nom, le workflow peut être importé dans d'autres scripts.

#### 1.2.2. Déclarer les entrées avec `take:`

Remplacez la déclaration de canal codée en dur par un bloc `take:` qui déclare les entrées attendues par le workflow.
Le bloc `take:` se place avant `main:`, et la ligne `names_ch = channel.of(...)` est supprimée.

=== "Après"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Canal d'entrée avec les noms

        main:
        // Chaîne de processus : valider -> créer la salutation -> ajouter l'horodatage
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Avant"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Chaîne de processus : valider -> créer la salutation -> ajouter l'horodatage
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

Le bloc `take:` déclare le canal par son nom uniquement. Le workflow parent définit ce qui y sera injecté.

#### 1.2.3. Déclarer les sorties avec `emit:`

Remplacez la section `publish:` et supprimez le bloc `output {}`, en les remplaçant par un bloc `emit:` qui nomme les sorties.

=== "Après"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Salutations originales
        timestamped = timestamped_ch // Salutations horodatées
    }
    ```

=== "Avant"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

Le bloc `emit:` expose des sorties nommées auxquelles les workflows parents peuvent accéder via `GREETING_WORKFLOW.out.greetings` et `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Vérifier le résultat et le tester

Après ces trois modifications, le fichier complet devrait ressembler à ceci :

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="7 8 9 11 17 18 19"
#!/usr/bin/env nextflow

include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Canal d'entrée avec les noms

    main:
    // Chaîne de processus : valider -> créer la salutation -> ajouter l'horodatage
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Salutations originales
    timestamped = timestamped_ch // Salutations horodatées
}
```

Essayez maintenant de l'exécuter directement :

```bash
nextflow run workflows/greeting.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `workflows/greeting.nf` [intergalactic_leavitt] revision: 89ac88c6c2
    No entry workflow specified
    ```

Cela introduit un concept clé : le **entry workflow**.
Nextflow utilise un bloc `workflow {}` sans nom comme point d'entrée lorsque vous exécutez un script directement.
`GREETING_WORKFLOW` est nommé, donc Nextflow ne sait pas comment l'exécuter seul.

C'est intentionnel. Les workflows composables sont conçus pour être appelés depuis un entry workflow, et non exécutés directement.
La solution consiste à créer un entry workflow dans `main.nf` qui importe et appelle `GREETING_WORKFLOW`.

### 1.3. Mettre à jour et tester le workflow principal

Mettez maintenant à jour le workflow principal pour appeler le greeting workflow.

#### 1.3.1. Inclure le greeting workflow et l'appeler

Ajoutez l'instruction `include`, mettez à jour le corps du workflow pour appeler `GREETING_WORKFLOW` et remplacez le placeholder `channel.empty()` dans `publish:` :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="3 9 10 13"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Exécuter le greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

Le entry workflow reste sans nom pour que Nextflow l'utilise comme point d'entrée du pipeline.

#### 1.3.2. Mettre à jour le bloc output

Ajoutez une directive `path` pour publier les salutations dans un sous-répertoire `greetings/` :

=== "Après"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Exécuter le workflow

Exécutez le workflow pour vérifier qu'il fonctionne :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `main.nf` [agitated_jones] revision: 6b0e98dd05
    executor >  local (9)
    [00/5ac4ce] GRE…DATE_NAME (validating Bob) | 3 of 3 ✔
    [7d/dde45c] GRE…SAY_HELLO (greeting Alice) | 3 of 3 ✔
    [fd/a044cc] GRE…ing timestamp to greeting) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      greetings:
        - greetings/Bob-output.txt
        - greetings/Charlie-output.txt
        - greetings/Alice-output.txt
    ```

??? abstract "Nouveau contenu ajouté dans `results/`"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

    Les six fichiers de l'exécution autonome de la section 1.1 (`Alice-output.txt`, `timestamped_Alice-output.txt`, etc.) se trouvent toujours dans `results/` aux côtés de ce nouveau sous-répertoire `greetings/`.
    C'est attendu : rien dans cette leçon ne les supprime, et la section 2.1 les lit directement en entrée.

??? abstract "Contenu du fichier"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Les fichiers de salutation sont publiés dans `results/greetings/`.
Le workflow principal appelle `GREETING_WORKFLOW` et connecte directement sa sortie à la section `publish:`.

`GREETING_WORKFLOW.out.timestamped` n'est pas connecté à `publish:` ici.
À partir de la section 2, ce canal devient l'entrée du `TRANSFORM_WORKFLOW` plutôt qu'une sortie publiée en tant que telle.
Comparez cela avec l'exécution autonome de la section 1.1, où `timestamped` n'avait rien d'autre à alimenter, et était donc publié directement.

### À retenir

Dans cette section, vous avez appris plusieurs concepts importants :

- **Workflows nommés** : Créer un workflow nommé (`GREETING_WORKFLOW`) qui peut être importé et réutilisé
- **Interfaces de workflow** : Définir des entrées claires avec `take:` et des sorties avec `emit:` pour créer un workflow composable
- **Points d'entrée** : Comprendre que Nextflow a besoin d'un entry workflow sans nom pour exécuter un script
- **Composition de workflows** : Importer et utiliser un workflow nommé au sein d'un autre workflow
- **Espaces de noms de workflow** : Accéder aux sorties d'un workflow en utilisant l'espace de noms `.out` (`GREETING_WORKFLOW.out.greetings`)

Vous disposez maintenant d'un greeting workflow fonctionnel qui :

- Prend un canal de noms en entrée
- Valide chaque nom
- Crée une salutation pour chaque nom valide
- Ajoute des horodatages aux salutations
- Émet les salutations originales et horodatées, bien que seules les salutations originales soient publiées dans `results/` à ce stade

Cette approche modulaire vous permet de tester le greeting workflow indépendamment ou de l'utiliser comme composant dans des pipelines plus grands.

---

## 2. Ajouter le transform workflow au pipeline

Le transform workflow applique des transformations de texte aux salutations horodatées.

### 2.1. Examiner et exécuter le workflow

Ouvrez `workflows/transform.nf` et examinez le code :

```groovy title="workflows/transform.nf" linenums="1"
#!/usr/bin/env nextflow

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Appliquer les transformations en séquence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Ce workflow autonome lit les fichiers de salutation horodatés depuis le répertoire `results/` produit par `greeting.nf`, les convertit en majuscules, puis inverse le texte.

Exécutez-le pour vérifier qu'il fonctionne avec les résultats du greeting de la section 1.1 :

```bash
nextflow run workflows/transform.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `workflows/transform.nf` [evil_meninsky] revision: 4e35790d32
    executor >  local (6)
    [7d/4ed7de] SAY…estamped_Alice-output.txt) | 3 of 3 ✔
    [70/ddaafa] REV…imestamped_Bob-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      upper:
        - UPPER-timestamped_Bob-output.txt
        - UPPER-timestamped_Charlie-output.txt
        - UPPER-timestamped_Alice-output.txt

      reversed:
        - REVERSED-UPPER-timestamped_Charlie-output.txt
        - REVERSED-UPPER-timestamped_Alice-output.txt
        - REVERSED-UPPER-timestamped_Bob-output.txt
    ```

Pour le rendre composable avec `GREETING_WORKFLOW`, les mêmes trois modifications de la section 1.2 s'appliquent.

### 2.2. Rendre le workflow composable

Appliquez les mêmes trois modifications que dans la section 1.2 : nommez le workflow, remplacez l'entrée codée en dur par `take:`, et remplacez `publish:`/`output {}` par `emit:`.

Le fichier final devrait ressembler à ceci :

```groovy title="workflows/transform.nf" linenums="1" hl_lines="6 7 8 10 15 16 17"
#!/usr/bin/env nextflow

include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Canal d'entrée avec les messages

    main:
    // Appliquer les transformations en séquence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Salutations en majuscules
    reversed = reversed_ch // Salutations en majuscules inversées
}
```

Le transform workflow est maintenant composable et prêt à être importé dans le workflow principal.

### 2.3. Mettre à jour et tester le workflow principal

Mettez maintenant à jour le workflow principal pour appeler le transform workflow.

#### 2.3.1. Inclure le transform workflow et l'appeler

Ajoutez l'instruction include, un appel à `TRANSFORM_WORKFLOW` enchaîné sur les salutations horodatées, et les deux nouvelles entrées `publish:` :

=== "Après"

    ```groovy title="main.nf" linenums="1" hl_lines="4 13 14 18 19"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Exécuter le greeting workflow
        GREETING_WORKFLOW(names)

        // Exécuter le transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Exécuter le greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Cela exécutera le transform workflow sur les salutations horodatées.

#### 2.3.2. Mettre à jour le bloc output

Ajoutez les entrées `upper` et `reversed` au bloc `output {}`, chacune avec une directive `path` pour son sous-répertoire :

=== "Après"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Avant"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Cela publiera les sorties finales dans les répertoires appropriés.

#### 2.3.3. Exécuter le pipeline complet

Exécutez le pipeline pour vérifier que tout fonctionne :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W   ~  version 26.04.4
    Launching `main.nf` [special_magritte] revision: 2a0e54fecd
    executor >  local (15)
    [e7/e051a5] GRE…TE_NAME (validating Alice) | 3 of 3 ✔
    [ca/cf1b18] GRE…SAY_HELLO (greeting Alice) | 3 of 3 ✔
    [e5/aa6be2] GRE…ing timestamp to greeting) | 3 of 3 ✔
    [d2/7afc88] TRA…imestamped_Bob-output.txt) | 3 of 3 ✔
    [5b/70f190] TRA…estamped_Alice-output.txt) | 3 of 3 ✔

    Outputs:

      /workspaces/training/side-quests/workflows_of_workflows/results

      greetings:
        - greetings/Bob-output.txt
        - greetings/Charlie-output.txt
        - greetings/Alice-output.txt

      upper:
        - upper/UPPER-timestamped_Charlie-output.txt
        - upper/UPPER-timestamped_Bob-output.txt
        - upper/UPPER-timestamped_Alice-output.txt

      reversed:
        - reversed/REVERSED-UPPER-timestamped_Charlie-output.txt
        - reversed/REVERSED-UPPER-timestamped_Bob-output.txt
        - reversed/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

??? abstract "Nouveau contenu ajouté dans `results/`"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

    La racine de `results/` contient également les fichiers résiduels des exécutions autonomes des sections 1.1 et 2.1 (`Alice-output.txt`, `timestamped_Alice-output.txt`, `UPPER-timestamped_Alice-output.txt`, etc.).
    C'est attendu : cette leçon ne vous demande jamais de les supprimer.

??? abstract "Contenu du fichier"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]54:84:81 31-90-6202[
    ```

Le pipeline fonctionne de bout en bout : la salutation a été convertie en majuscules et inversée.

### À retenir

Vous devriez maintenant disposer d'un pipeline complet qui :

- Traite les noms à travers le greeting workflow
- Transmet les salutations horodatées au transform workflow
- Produit des versions en majuscules et inversées des salutations

---

## Résumé

Dans cette quête secondaire, nous avons exploré le puissant concept de composition de workflows dans Nextflow, qui nous permet de construire des pipelines complexes à partir de composants plus petits et réutilisables.

Cette approche modulaire offre plusieurs avantages par rapport aux pipelines monolithiques :

- Chaque workflow peut être développé, testé et débogué indépendamment
- Les workflows peuvent être réutilisés dans différents pipelines
- La structure globale du pipeline devient plus lisible et maintenable
- Les modifications apportées à un workflow n'affectent pas nécessairement les autres si les interfaces restent cohérentes
- Les points d'entrée peuvent être configurés pour exécuter différentes parties de votre pipeline selon les besoins

Appeler un workflow est similaire à appeler un processus, mais ce n'est pas tout à fait la même chose.
Vous ne pouvez pas exécuter un workflow N fois en lui passant un canal de taille N : vous passez le canal une seule fois et le workflow itère en interne.

L'application de ces techniques dans votre propre travail vous permettra de construire des pipelines Nextflow plus sophistiqués, capables de gérer des tâches de traitement de données complexes tout en restant maintenables et évolutifs.

### Modèles clés

1.  **Structure du workflow** : Nous avons défini des entrées et des sorties claires pour chaque workflow en utilisant la syntaxe `take:` et `emit:`, créant ainsi des interfaces bien définies entre les composants, et encapsulé la logique du workflow dans le bloc `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
        // Les canaux d'entrée sont déclarés ici
        input_ch

        main:
        // La logique du workflow se trouve ici
        // C'est ici que les processus sont appelés et que les canaux sont manipulés
        result_ch = SOME_PROCESS(input_ch)

        emit:
        // Les canaux de sortie sont déclarés ici
        output_ch = result_ch
    }
    ```

2.  **Imports de workflows** : Nous avons construit deux modules de workflow indépendants et les avons importés dans un pipeline principal à l'aide d'instructions include.

    - Importer un seul workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Importer plusieurs workflows

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Importer avec un alias pour éviter les conflits de noms

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Points d'entrée** : Nextflow nécessite un entry workflow sans nom pour savoir où démarrer l'exécution.
    Ce entry workflow appelle vos workflows nommés.

    - Workflow sans nom (point d'entrée)

    ```groovy
    workflow {
        // C'est le point d'entrée lorsque le script est exécuté
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow nommé (appelé depuis le entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Doit être appelé depuis le entry workflow
    }
    ```

4.  **Gestion du flux de données** : Nous avons appris à accéder aux sorties d'un workflow en utilisant la notation d'espace de noms (`WORKFLOW_NAME.out.channel_name`) et à les transmettre à d'autres workflows.

    ```groovy
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Ressources supplémentaires

- [Documentation Nextflow sur les Workflows](https://www.nextflow.io/docs/latest/workflow.html)
- [Référence des opérateurs de canaux](https://www.nextflow.io/docs/latest/operator.html)
- [Documentation sur les stratégies d'erreur](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](../index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant de la liste.
