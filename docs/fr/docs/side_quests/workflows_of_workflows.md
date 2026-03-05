# Workflows de workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Lorsque vous développez un pipeline, vous vous retrouvez souvent à créer des séquences similaires de processus pour différents types de données ou étapes d'analyse. Vous pourriez finir par copier et coller ces séquences de processus, ce qui conduit à du code dupliqué difficile à maintenir ; ou vous pourriez créer un workflow massif qui est difficile à comprendre et à modifier.

L'une des fonctionnalités les plus puissantes de Nextflow est sa capacité à composer des pipelines complexes à partir de modules de workflow plus petits et réutilisables. Cette approche modulaire rend les pipelines plus faciles à développer, tester et maintenir.

### Objectifs d'apprentissage

Dans cette quête secondaire, nous explorerons comment développer des modules de workflow qui peuvent être testés et utilisés séparément, composer ces modules dans un pipeline plus large, et gérer le flux de données entre les modules.

À la fin de cette quête secondaire, vous serez capable de :

- Décomposer des pipelines complexes en unités logiques et réutilisables
- Tester chaque module de workflow indépendamment
- Mélanger et associer des workflows pour créer de nouveaux pipelines
- Partager des modules de workflow communs entre différents pipelines
- Rendre votre code plus maintenable et plus facile à comprendre

Ces compétences vous aideront à construire des pipelines complexes tout en maintenant une structure de code claire et maintenable.

### Prérequis

Avant d'entreprendre cette quête secondaire, vous devriez :

- Avoir complété le tutoriel [Hello Nextflow](../hello_nextflow/README.md) ou un cours équivalent pour débutants
- Être à l'aise avec l'utilisation des concepts et mécanismes de base de Nextflow (processus, canaux, opérateurs, modules)

---

## 0. Pour commencer

#### Ouvrir le codespace de formation

Si vous ne l'avez pas encore fait, assurez-vous d'ouvrir l'environnement de formation comme décrit dans la [Configuration de l'environnement](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Se déplacer dans le répertoire du projet

Déplaçons-nous dans le répertoire où se trouvent les fichiers pour ce tutoriel.

```bash
cd side-quests/workflows_of_workflows
```

Vous pouvez configurer VSCode pour se concentrer sur ce répertoire :

```bash
code .
```

#### Examiner le matériel

Vous trouverez un répertoire `modules` contenant plusieurs définitions de processus qui s'appuient sur ce que vous avez appris dans 'Hello Nextflow' :

```console title="Contenu du répertoire"
modules/
├── say_hello.nf             # Crée un message de bienvenue (de Hello Nextflow)
├── say_hello_upper.nf       # Convertit en majuscules (de Hello Nextflow)
├── timestamp_greeting.nf    # Ajoute des horodatages aux messages
├── validate_name.nf         # Valide les noms en entrée
└── reverse_text.nf          # Inverse le contenu du texte
```

#### Examiner l'exercice

Votre défi est d'assembler ces modules en deux workflows séparés que nous composerons ensuite dans un workflow principal :

- Un `GREETING_WORKFLOW` qui valide les noms, crée des messages de bienvenue et ajoute des horodatages
- Un `TRANSFORM_WORKFLOW` qui convertit le texte en majuscules et l'inverse

#### Liste de vérification de préparation

Pensez-vous être prêt·e à vous lancer ?

- [ ] Je comprends l'objectif de ce cours et ses prérequis
- [ ] Mon codespace est opérationnel
- [ ] J'ai défini mon répertoire de travail de manière appropriée
- [ ] Je comprends l'exercice

Si vous pouvez cocher toutes les cases, vous êtes prêt·e à commencer.

---

## 1. Créer le workflow de bienvenue

Commençons par créer un workflow qui valide les noms et génère des messages de bienvenue horodatés.

### 1.1. Créer la structure du workflow

```bash title="Créer le répertoire et le fichier du workflow"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Ajouter le code du premier (sous-)workflow

Ajoutez ce code à `workflows/greeting.nf` :

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Enchaîner les processus : valider -> créer message -> ajouter horodatage
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Il s'agit d'un workflow complet, avec une structure similaire à ceux que vous avez vus dans le tutoriel 'Hello Nextflow', que nous pouvons tester indépendamment. Essayons cela maintenant :

```bash
nextflow run workflows/greeting.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Cela fonctionne comme prévu, mais pour le rendre composable, il y a quelques modifications à apporter.

### 1.3. Rendre le workflow composable

Les workflows composables présentent quelques différences par rapport à ceux que vous avez vus dans le tutoriel 'Hello Nextflow' :

- Le bloc workflow doit être nommé
- Les entrées sont déclarées en utilisant le mot-clé `take:`
- Le contenu du workflow est placé à l'intérieur du bloc `main:`
- Les sorties sont déclarées en utilisant le mot-clé `emit:`

Mettons à jour le workflow de bienvenue pour correspondre à cette structure. Modifiez le code comme suit :

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Canal d'entrée avec les noms

    main:
        // Enchaîner les processus : valider -> créer message -> ajouter horodatage
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Messages originaux
        timestamped = timestamped_ch  // Messages horodatés
}
```

Vous pouvez voir que le workflow est maintenant nommé et possède un bloc `take:` et `emit:`, et ce sont les connexions que nous utiliserons pour composer un workflow de niveau supérieur.
Le contenu du workflow est également placé à l'intérieur du bloc `main:`. Notez également que nous avons supprimé la déclaration du canal d'entrée `names_ch`, car il est maintenant passé en argument au workflow.

Testons à nouveau le workflow pour voir s'il fonctionne comme prévu :

```bash
nextflow run workflows/greeting.nf
```

??? failure "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Cela vous informe d'un autre nouveau concept, un 'workflow d'entrée'. Le workflow d'entrée est le workflow qui est appelé lorsque vous exécutez un script Nextflow. Par défaut, Nextflow utilisera un workflow sans nom comme workflow d'entrée, lorsqu'il est présent, et c'est ce que vous avez fait jusqu'à présent, avec des blocs workflow commençant comme ceci :

```groovy title="hello.nf" linenums="1"
workflow {
```

Mais notre workflow de bienvenue n'a pas de workflow sans nom, nous avons plutôt un workflow nommé :

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

C'est pourquoi Nextflow a généré une erreur et n'a pas fait ce que nous voulions.

Nous n'avons pas ajouté la syntaxe `take:`/`emit:` pour pouvoir appeler le workflow directement - nous l'avons fait pour pouvoir le composer avec d'autres workflows. La solution est de créer un script principal avec un workflow d'entrée sans nom qui importe et appelle notre workflow nommé.

### 1.4. Créer et tester le workflow principal

Nous allons maintenant créer un workflow principal qui importe et utilise le workflow `greeting`.

Créez `main.nf` :

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Notez que notre entrée de workflow dans ce fichier est sans nom, et c'est parce que nous allons l'utiliser comme workflow d'entrée.

Exécutez ceci et observez la sortie :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Ça fonctionne ! Nous avons enveloppé le workflow de bienvenue nommé dans un workflow principal avec un bloc `workflow` d'entrée sans nom. Le workflow principal utilise le workflow `GREETING_WORKFLOW` presque (mais pas tout à fait) comme un processus, et passe le canal `names` en argument.

### À retenir

Dans cette section, vous avez appris plusieurs concepts importants :

- **Workflows nommés** : Créer un workflow nommé (`GREETING_WORKFLOW`) qui peut être importé et réutilisé
- **Interfaces de workflow** : Définir des entrées claires avec `take:` et des sorties avec `emit:` pour créer un workflow composable
- **Points d'entrée** : Comprendre que Nextflow a besoin d'un workflow d'entrée sans nom pour exécuter un script
- **Composition de workflows** : Importer et utiliser un workflow nommé dans un autre workflow
- **Espaces de noms de workflow** : Accéder aux sorties de workflow en utilisant l'espace de noms `.out` (`GREETING_WORKFLOW.out.greetings`)

Vous avez maintenant un workflow de bienvenue fonctionnel qui :

- Prend un canal de noms en entrée
- Valide chaque nom
- Crée un message de bienvenue pour chaque nom valide
- Ajoute des horodatages aux messages
- Expose à la fois les messages originaux et horodatés comme sorties

Cette approche modulaire vous permet de tester le workflow de bienvenue indépendamment ou de l'utiliser comme composant dans des pipelines plus grands.

---

## 2. Ajouter le workflow de transformation

Créons maintenant un workflow qui applique des transformations de texte aux messages de bienvenue.

### 2.1. Créer le fichier du workflow

```bash
touch workflows/transform.nf
```

### 2.2. Ajouter le code du workflow

Ajoutez ce code à `workflows/transform.nf` :

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Canal d'entrée avec les messages

    main:
        // Appliquer les transformations en séquence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Messages en majuscules
        reversed = reversed_ch  // Messages en majuscules inversés
}
```

Nous ne répéterons pas l'explication de la syntaxe composable ici, mais notez que le workflow nommé est à nouveau déclaré avec un bloc `take:` et `emit:`, et le contenu du workflow est placé à l'intérieur du bloc `main:`.

### 2.3. Mettre à jour le workflow principal

Mettez à jour `main.nf` pour utiliser les deux workflows :

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Exécuter le workflow de bienvenue
    GREETING_WORKFLOW(names)

    // Exécuter le workflow de transformation
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Afficher les résultats
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Exécutez le pipeline complet :

```bash
nextflow run main.nf
```

??? success "Sortie de la commande"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Si vous examinez l'un de ces fichiers inversés, vous verrez qu'il s'agit de la version en majuscules du message inversé :

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Contenu du fichier inversé"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### À retenir

Vous devriez maintenant avoir un pipeline complet qui :

- Traite les noms à travers le workflow de bienvenue
- Alimente les messages horodatés dans le workflow de transformation
- Produit à la fois les versions en majuscules et inversées des messages

---

## Résumé

Dans cette quête secondaire, nous avons exploré le concept puissant de composition de workflows dans Nextflow, qui nous permet de construire des pipelines complexes à partir de composants plus petits et réutilisables.

Cette approche modulaire offre plusieurs avantages par rapport aux pipelines monolithiques :

- Chaque workflow peut être développé, testé et débogué indépendamment
- Les workflows peuvent être réutilisés dans différents pipelines
- La structure globale du pipeline devient plus lisible et maintenable
- Les modifications apportées à un workflow n'affectent pas nécessairement les autres si les interfaces restent cohérentes
- Les points d'entrée peuvent être configurés pour exécuter différentes parties de votre pipeline selon les besoins

_Il est important de noter cependant que, bien que l'appel de workflows ressemble un peu à l'appel de processus, ce n'est pas réellement la même chose. Vous ne pouvez pas, par exemple, exécuter un workflow N fois en l'appelant avec un canal de taille N - vous devriez passer un canal de taille N au workflow et itérer en interne._

L'application de ces techniques dans votre propre travail vous permettra de construire des pipelines Nextflow plus sophistiqués capables de gérer des tâches bioinformatiques complexes tout en restant maintenables et évolutifs.

### Modèles clés

1.  **Structure de workflow** : Nous avons défini des entrées et sorties claires pour chaque workflow en utilisant la syntaxe `take:` et `emit:`, créant des interfaces bien définies entre les composants, et enveloppé la logique du workflow dans le bloc `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Les canaux d'entrée sont déclarés ici
            input_ch

        main:
            // La logique du workflow va ici
            // C'est ici que les processus sont appelés et les canaux manipulés
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Les canaux de sortie sont déclarés ici
            output_ch = result_ch
    }
    ```

2.  **Imports de workflow :** Nous avons construit deux modules de workflow indépendants et les avons importés dans un pipeline principal avec des instructions include.

    - Inclure un seul workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Inclure plusieurs workflows

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Inclure avec un alias pour éviter les conflits de noms

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Points d'entrée** : Nextflow nécessite un workflow d'entrée sans nom pour savoir où commencer l'exécution. Ce workflow d'entrée appelle vos workflows nommés.

    - Workflow sans nom (point d'entrée)

    ```groovy
    workflow {
        // Ceci est le point d'entrée lorsque le script est exécuté
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow nommé (appelé depuis le workflow d'entrée)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Doit être appelé depuis le workflow d'entrée
    }
    ```

4.  **Gestion du flux de données :** Nous avons appris comment accéder aux sorties de workflow en utilisant la notation d'espace de noms (`WORKFLOW_NAME.out.channel_name`) et les passer à d'autres workflows.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Ressources supplémentaires

- [Documentation Nextflow sur les workflows](https://www.nextflow.io/docs/latest/workflow.html)
- [Référence des opérateurs de canaux](https://www.nextflow.io/docs/latest/operator.html)
- [Documentation sur la stratégie d'erreur](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Et ensuite ?

Retournez au [menu des Quêtes secondaires](./index.md) ou cliquez sur le bouton en bas à droite de la page pour passer au sujet suivant dans la liste.
