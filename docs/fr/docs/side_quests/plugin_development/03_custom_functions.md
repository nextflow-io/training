# Partie 3 : Fonctions personnalisées

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduction assistée par IA - [en savoir plus et suggérer des améliorations](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

À la fin de cette section, vous aurez des fonctions personnalisées dans votre plugin, compilées et installées localement, et fonctionnant dans un vrai workflow.

!!! tip "Astuce"

    Vous commencez à partir de cette partie ? Copiez la solution de la Partie 2 pour l'utiliser comme point de départ :

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Examiner ce que le template a généré

Avant d'écrire vos propres fonctions, examinez la fonction d'exemple créée par le template pour comprendre le modèle utilisé.

Placez-vous dans le répertoire du plugin :

```bash
cd nf-greeting
```

Le template a créé un fichier appelé `GreetingExtension.groovy` dans lequel les fonctions du plugin sont définies.
Ouvrez-le pour voir le point de départ :

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

```groovy title="Output" hl_lines="29 40-43"
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
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Implements a custom function which can be imported by
 * Nextflow scripts.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Say hello to the given target.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. La classe sur laquelle votre extension s'appuie. Nextflow en a besoin pour reconnaître vos fonctions.
2. Appelée au chargement du plugin ; à utiliser pour l'initialisation
3. Rend cette méthode appelable depuis les workflows via `include`

Le template inclut une fonction `sayHello` d'exemple.
L'annotation `@Function` est ce qui rend une méthode appelable depuis les workflows Nextflow.
Sans elle, la méthode n'existe qu'à l'intérieur du code du plugin.

En Groovy (et en Java), les méthodes déclarent le type qu'elles retournent et les types de leurs paramètres.
Par exemple, `String reverseGreeting(String greeting)` déclare une méthode qui prend un paramètre de type `String` et retourne un `String`.
Le mot-clé `void` signifie que la méthode ne retourne rien, comme avec `sayHello` ci-dessus.
C'est différent de Python ou R, où les types n'ont pas besoin d'être déclarés explicitement.

---

## 2. Remplacer sayHello par reverseGreeting

La fonction `sayHello` du template est un exemple de départ.
Remplacez-la par votre propre fonction pour voir le cycle complet d'écriture, de compilation et d'utilisation d'une fonction de plugin.

Modifiez `src/main/groovy/training/plugin/GreetingExtension.groovy` pour remplacer la méthode `sayHello` :

=== "Après"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverse une chaîne de salutation
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Rend la méthode appelable depuis les workflows Nextflow
    2. Prend un String, retourne un String
    3. La méthode d'inversion de chaîne intégrée à Groovy

=== "Avant"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implements a custom function which can be imported by
     * Nextflow scripts.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Say hello to the given target.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Éléments clés de cette fonction :

- **`@Function`** : Rend la méthode appelable depuis les workflows Nextflow
- **`String reverseGreeting(String greeting)`** : Prend un String, retourne un String
- **`greeting.reverse()`** : La méthode d'inversion de chaîne intégrée à Groovy

!!! tip "Méthodes publiques et privées"

    Les méthodes sans `@Function` ne sont pas exposées aux workflows Nextflow.
    Vous pouvez ajouter des méthodes utilitaires à votre classe sans craindre qu'elles ne s'infiltrent dans l'espace de noms du workflow.

---

## 3. Compiler et installer votre plugin

Compilez et installez le plugin :

```bash
make install
```

!!! tip "Si la compilation échoue"

    Lisez attentivement le message d'erreur ; il indique généralement un numéro de ligne et décrit le problème.
    Les causes fréquentes sont les erreurs de syntaxe (accolade ou guillemet manquant), les noms de classes mal orthographiés et les incompatibilités de types.
    Si vous êtes bloqué·e, comparez votre code caractère par caractère avec les exemples.

---

## 4. Utiliser votre fonction dans un workflow

Le plugin est compilé et installé.
L'étape suivante consiste à utiliser `reverseGreeting` dans un workflow pour vérifier qu'il fonctionne de bout en bout.

Retournez dans le répertoire du pipeline :

```bash
cd ..
```

Modifiez `greet.nf` pour importer et utiliser `reverseGreeting` :

=== "Après"

    ```groovy title="greet.nf" hl_lines="4 23-25" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Avant"

    ```groovy title="greet.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Exécutez le pipeline :

```bash
nextflow run greet.nf
```

??? example "Sortie"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Output: Hello
    Output: Bonjour
    Output: Holà
    Output: Ciao
    Output: Hallo
    Pipeline complete! 👋
    ```

Votre première fonction de plugin personnalisée fonctionne dans un vrai workflow.
Le même modèle `include { ... } from 'plugin/...'` que vous avez utilisé avec nf-hello et nf-schema dans la Partie 1 fonctionne avec votre propre plugin.

---

## 5. Ajouter decorateGreeting

Un plugin peut fournir plusieurs fonctions.
Ajoutez-en une deuxième qui encadre une salutation avec des marqueurs décoratifs ; vous la rendrez configurable dans la Partie 6.

Modifiez `GreetingExtension.groovy` pour ajouter `decorateGreeting` après `reverseGreeting`, avant l'accolade fermante de la classe :

=== "Après"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverse une chaîne de salutation
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Décore une salutation avec des marqueurs festifs
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Interpolation de chaîne Groovy : `#!groovy ${...}` insère la valeur de la variable dans la chaîne

=== "Avant"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverse une chaîne de salutation
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Cette fonction utilise l'interpolation de chaîne Groovy (`"*** ${greeting} ***"`) pour intégrer la variable de salutation dans une chaîne.

Compilez, installez et mettez à jour le workflow :

```bash
cd nf-greeting && make install && cd ..
```

Mettez à jour `greet.nf` pour importer et utiliser également `decorateGreeting` :

=== "Après"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Importe les fonctions personnalisées depuis notre plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Utilise notre fonction de plugin personnalisée pour décorer la salutation
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Démontre l'utilisation de la fonction reverseGreeting
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Plusieurs fonctions provenant du même plugin nécessitent des instructions `include` séparées
    2. Les fonctions de plugin fonctionnent également à l'intérieur des blocs `script:` des processus

=== "Avant"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 12 15 28"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

```bash
nextflow run greet.nf
```

??? example "Sortie"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    ```

Les fonctions de plugin fonctionnent aussi bien dans les scripts de processus (comme `decorateGreeting` à l'intérieur de `SAY_HELLO`) que dans les opérations de workflow (comme `reverseGreeting` dans un `map`).

---

## À retenir

Vous avez appris que :

- Les fonctions sont définies avec l'annotation `@Function` dans les sous-classes de `PluginExtensionPoint`
- Les fonctions de plugin importées avec `include` fonctionnent de manière identique qu'elles proviennent de votre propre plugin ou d'un plugin existant
- Les fonctions de plugin fonctionnent aussi bien dans les scripts de processus que dans les opérations de workflow

---

## Et ensuite ?

Vos fonctions fonctionnent, mais jusqu'à présent vous n'avez vérifié cela qu'en exécutant le pipeline complet et en contrôlant la sortie visuellement.
Cette approche ne passe pas à l'échelle : à mesure que vous ajoutez des fonctions, vous avez besoin d'un moyen plus rapide de vérifier que chacune se comporte correctement, notamment après avoir effectué des modifications.
La section suivante présente les tests unitaires, qui vous permettent de vérifier des fonctions individuelles automatiquement sans exécuter de pipeline.

[Continuer vers la Partie 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
