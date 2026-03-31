# Partie 5 : Observateurs de trace

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Les observateurs de trace permettent à votre plugin de répondre aux événements du workflow, comme la fin d'une tâche, la publication d'un fichier, ou la fin du pipeline.
Cela permet des cas d'usage tels que des rapports personnalisés, des notifications Slack, la collecte de métriques, ou l'intégration avec des systèmes de surveillance externes.
Dans cette section, vous allez créer un observateur qui compte les tâches terminées et affiche un résumé.

!!! tip "Astuce"

    Si vous commencez à partir de cette partie, copiez la solution de la Partie 4 pour l'utiliser comme point de départ :

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Comprendre l'observateur de trace existant

Le message "Pipeline is starting!" affiché lors de l'exécution du pipeline provient de la classe `GreetingObserver` de votre plugin.

Examinez le code de l'observateur :

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingObserver.groovy
```

```groovy title="Output" hl_lines="30 32-34 37-39"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.trace.TraceObserver

/**
 * Implements an observer that allows implementing custom
 * logic on nextflow execution events.
 */
@Slf4j
@CompileStatic
class GreetingObserver implements TraceObserver {    // (1)!

    @Override
    void onFlowCreate(Session session) {            // (2)!
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {                         // (3)!
        println "Pipeline complete! 👋"
    }
}
```

1. Interface permettant de s'accrocher aux événements du cycle de vie du workflow
2. Appelée au démarrage du workflow ; reçoit la session pour accéder à la configuration
3. Appelée lorsque le workflow se termine avec succès

Deux éléments sont à noter ici :

1. **`class GreetingObserver implements TraceObserver`** : `TraceObserver` est une interface définie par Nextflow. Si votre classe implémente cette interface, Nextflow peut s'y accrocher et appeler vos méthodes lorsque des événements se produisent.
2. **`@Override`** : L'interface `TraceObserver` définit des méthodes comme `onFlowCreate` et `onFlowComplete`. Lorsque vous écrivez des méthodes portant ces noms et que vous ajoutez l'annotation `@Override`, Nextflow les appelle au moment approprié. Toutes les méthodes que vous ne surchargez pas sont ignorées.

L'ensemble complet des événements du cycle de vie auxquels vous pouvez vous accrocher au moment de la rédaction est le suivant :

| Méthode             | Moment d'appel                    |
| ------------------- | --------------------------------- |
| `onFlowCreate`      | Démarrage du workflow             |
| `onFlowComplete`    | Fin du workflow                   |
| `onProcessStart`    | Début d'exécution d'une tâche     |
| `onProcessComplete` | Fin d'une tâche                   |
| `onProcessCached`   | Réutilisation d'une tâche en cache |
| `onFilePublish`     | Publication d'un fichier          |

Pour une liste complète, consultez l'[interface TraceObserver](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) dans le code source de Nextflow.

---

## 2. Ajouter un observateur de comptage de tâches

L'objectif est de créer un observateur qui compte les tâches terminées et affiche un résumé à la fin.
L'ajout d'un nouvel observateur à un plugin nécessite deux choses : écrire la classe de l'observateur, et l'enregistrer dans la factory pour que Nextflow la charge.

### 2.1. Créer un observateur minimal

Créez un nouveau fichier :

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Commencez par l'observateur le plus simple possible, qui affiche un message lorsqu'une tâche se termine :

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observateur qui réagit à la fin des tâches
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Importez les classes requises : `TraceObserver`, `TaskHandler` et `TraceRecord`
2. Créez une classe qui `implements TraceObserver`
3. Surchargez `onProcessComplete` pour exécuter du code lorsqu'une tâche se termine

Voici le minimum nécessaire :

- Importer les classes requises (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Créer une classe qui `implements TraceObserver`
- Surcharger `onProcessComplete` pour effectuer une action lorsqu'une tâche se termine

### 2.2. Enregistrer l'observateur

La `GreetingFactory` crée les observateurs.
Examinez-la :

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingFactory.groovy
```

```groovy title="Output" hl_lines="25 27-29"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }

}
```

Modifiez `GreetingFactory.groovy` pour ajouter le nouvel observateur :

=== "Après"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Avant"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Syntaxe de liste Groovy"

    Nous avons remplacé le style Java `List.<TraceObserver>of(...)` par le littéral de liste Groovy plus simple `[...]`.
    Les deux retournent une `Collection`, mais la syntaxe Groovy est plus lisible lorsqu'on ajoute plusieurs éléments.

### 2.3. Compiler, installer et tester

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Pourquoi `-ansi-log false` ?"

    Par défaut, l'affichage de progression ANSI de Nextflow écrase les lignes précédentes pour afficher une vue propre et actualisée de la progression.
    Cela signifie que vous ne verriez que le *dernier* comptage de tâches, et non les messages intermédiaires.

    L'utilisation de `-ansi-log false` désactive ce comportement et affiche toutes les sorties de manière séquentielle, ce qui est indispensable lors du test d'observateurs qui affichent des messages pendant l'exécution.

Vous devriez voir "✓ Task completed!" affiché cinq fois (une fois par tâche), entremêlé avec la sortie existante du pipeline :

```console title="Output (partial)"
...
[9b/df7630] Submitted process > SAY_HELLO (4)
Decorated: *** Hello ***
✓ Task completed!
✓ Task completed!
Decorated: *** Holà ***
✓ Task completed!
...
Pipeline complete! 👋
```

L'observateur fonctionne.
Chaque fois qu'une tâche se termine, Nextflow appelle `onProcessComplete`, et notre implémentation affiche un message.

??? exercise "Personnaliser le message"

    Essayez de modifier le message dans `onProcessComplete` avec quelque chose de votre choix, recompilez et relancez.
    Cela confirme que le cycle complet modification-compilation-exécution fonctionne pour les observateurs.

### 2.4. Ajouter la logique de comptage

L'observateur minimal prouve que le hook fonctionne, mais il ne suit rien.

Une classe peut contenir des variables (appelées champs ou variables d'instance) qui persistent pendant toute la durée de vie de l'objet.
Cela signifie qu'un observateur peut accumuler un état à travers plusieurs événements pendant l'exécution d'un pipeline.

La version suivante ajoute une variable compteur (`taskCount`) qui commence à zéro.
Chaque fois qu'une tâche se termine, le compteur augmente d'une unité.
Lorsque l'ensemble du workflow se termine, l'observateur affiche le total final.

Mettez à jour `TaskCounterObserver.groovy` avec les modifications mises en évidence :

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observateur qui compte les tâches terminées
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0                // (1)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++                          // (2)!
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFlowComplete() {                  // (3)!
        println "📈 Final task count: ${taskCount}"
    }
}
```

1. `taskCount` est une variable qui appartient à l'objet observateur. Elle conserve sa valeur entre les appels de méthodes, ce qui lui permet d'accumuler un comptage sur l'ensemble de l'exécution du workflow. `private` signifie que seule cette classe peut y accéder.
2. `taskCount++` ajoute un au compteur. Cette ligne s'exécute chaque fois qu'une tâche se termine, de sorte que le comptage augmente au fur et à mesure de la progression du workflow.
3. `onFlowComplete` est un second hook du cycle de vie. Il s'exécute une fois lorsque le workflow se termine, ce qui en fait un bon endroit pour afficher un résumé.

En résumé :

- `taskCount` persiste entre les appels de méthodes, accumulant un comptage sur toute l'exécution
- `onProcessComplete` incrémente le compteur et affiche le total courant chaque fois qu'une tâche se termine
- `onFlowComplete` s'exécute une fois à la fin, affichant le comptage final

Recompilez et testez :

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Sortie"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `greet.nf` [pensive_engelbart] DSL2 - revision: 85fefd90d0
    Pipeline is starting! 🚀
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    [be/bd8e72] Submitted process > SAY_HELLO (2)
    [5b/d24c2b] Submitted process > SAY_HELLO (1)
    [14/1f9dbe] Submitted process > SAY_HELLO (3)
    Decorated: *** Bonjour ***
    Decorated: *** Hello ***
    [85/a6b3ad] Submitted process > SAY_HELLO (4)
    📊 Tasks completed so far: 1
    📊 Tasks completed so far: 2
    Decorated: *** Holà ***
    📊 Tasks completed so far: 3
    Decorated: *** Ciao ***
    [3c/be6686] Submitted process > SAY_HELLO (5)
    📊 Tasks completed so far: 4
    Decorated: *** Hallo ***
    📊 Tasks completed so far: 5
    Pipeline complete! 👋
    📈 Final task count: 5
    ```

    Les messages du compteur sont entremêlés avec les soumissions de tâches car les observateurs s'exécutent au fur et à mesure que les tâches se terminent.

---

## 3. Suivre les fichiers publiés

L'observateur peut également réagir lorsque des fichiers sont publiés.
La méthode `onFilePublish` reçoit les chemins de destination et de source, que vous pouvez utiliser pour journaliser, valider ou traiter les sorties publiées.

### 3.1. Ajouter un répertoire de publication

Tout d'abord, mettez à jour `greet.nf` pour que le processus `SAY_HELLO` publie ses fichiers de sortie :

=== "Après"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Utilise notre fonction de plugin personnalisée pour décorer le message
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Avant"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Utilise notre fonction de plugin personnalisée pour décorer le message
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Ajouter la méthode onFilePublish

Ajoutez une méthode `onFilePublish` et l'import requis à `TaskCounterObserver.groovy` :

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observateur qui compte les tâches terminées
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "📁 Published: ${destination.fileName}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
```

### 3.3. Compiler et tester

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Vous devriez voir des messages "Published:" pour chaque fichier de sortie, aux côtés de la sortie du compteur de tâches :

```console title="Output (partial)"
...
📊 Tasks completed so far: 1
📁 Published: greeting.txt
📊 Tasks completed so far: 2
📁 Published: greeting.txt
...
📈 Final task count: 5
Pipeline complete! 👋
```

La méthode `onFilePublish` se déclenche chaque fois que Nextflow publie un fichier dans le répertoire `results`.
Ce modèle est utile pour créer des journaux d'audit, déclencher des actions en aval, ou valider les sorties au fur et à mesure de leur production.

---

## À retenir

Vous avez appris que :

- Les observateurs de trace s'accrochent aux événements du cycle de vie du workflow comme `onFlowCreate`, `onProcessComplete`, `onFilePublish` et `onFlowComplete`
- On crée des observateurs en implémentant `TraceObserver` et en les enregistrant dans une Factory
- Les observateurs peuvent contenir des variables d'instance pour accumuler un état à travers les événements
- Les observateurs sont utiles pour la journalisation personnalisée, la collecte de métriques, les notifications et les rapports

---

## Et ensuite ?

Le compteur de tâches fonctionne, mais il est toujours actif.
Dans un vrai plugin, les utilisateur·trices devraient pouvoir activer ou désactiver des fonctionnalités, ou ajuster le comportement, depuis `nextflow.config` sans modifier le code source du plugin.
La section suivante montre comment rendre votre observateur configurable et comment partager votre plugin terminé avec d'autres personnes.

[Continuer vers la Partie 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
