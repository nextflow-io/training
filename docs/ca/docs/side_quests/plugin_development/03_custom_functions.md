# Part 3: Funcions personalitzades

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Al final d'aquesta secció, tindreu funcions personalitzades al vostre plugin, compilades i instal·lades localment, i en execució en un workflow real.

!!! tip "Comenceu des d'aquí?"

    Si us incorporeu en aquesta part, copieu la solució de la Part 2 per utilitzar-la com a punt de partida:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Veieu què ha generat la plantilla

Abans d'escriure les vostres pròpies funcions, observeu la funció d'exemple que ha creat la plantilla per entendre el patró.

Canvieu al directori del plugin:

```bash
cd nf-greeting
```

La plantilla ha creat un fitxer anomenat `GreetingExtension.groovy` on es defineixen les funcions del plugin.
Obriu-lo per veure el punt de partida:

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
 * Implementa una funció personalitzada que pot ser importada per
 * scripts de Nextflow.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Saluda el destinatari indicat.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. La classe sobre la qual es construeix la vostra extensió. Nextflow la requereix per reconèixer les vostres funcions.
2. S'invoca quan es carrega el plugin; utilitzeu-la per a la inicialització.
3. Fa que aquest mètode sigui invocable des dels workflows mitjançant `include`.

La plantilla inclou una funció d'exemple `sayHello`.
L'anotació `@Function` és el que fa que un mètode sigui invocable des dels workflows de Nextflow.
Sense ella, el mètode existeix únicament dins del codi del plugin.

En Groovy (i Java), els mètodes declaren quin tipus retornen i quins tipus tenen els seus paràmetres.
Per exemple, `String reverseGreeting(String greeting)` declara un mètode que pren un paràmetre de tipus `String` i retorna un `String`.
La paraula clau `void` significa que el mètode no retorna res, com en el cas de `sayHello` anterior.
Això és diferent de Python o R, on no cal declarar els tipus explícitament.

---

## 2. Substituïu sayHello per reverseGreeting

La funció `sayHello` de la plantilla és un marcador de posició.
Substituïu-la per la vostra pròpia funció per veure el cicle complet d'escriptura, compilació i ús d'una funció de plugin.

Editeu `src/main/groovy/training/plugin/GreetingExtension.groovy` per substituir el mètode `sayHello`:

=== "Després"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverteix una cadena de salutació
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Fa que el mètode sigui invocable des dels workflows de Nextflow.
    2. Pren un String i retorna un String.
    3. El mètode d'inversió de cadenes integrat de Groovy.

=== "Abans"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implementa una funció personalitzada que pot ser importada per
     * scripts de Nextflow.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Saluda el destinatari indicat.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Parts clau d'aquesta funció:

- **`@Function`**: Fa que el mètode sigui invocable des dels workflows de Nextflow.
- **`String reverseGreeting(String greeting)`**: Pren un String i retorna un String.
- **`greeting.reverse()`**: El mètode d'inversió de cadenes integrat de Groovy.

!!! tip "Mètodes públics i privats"

    Els mètodes sense `@Function` no s'exposen als workflows de Nextflow.
    Podeu afegir mètodes auxiliars a la vostra classe sense preocupar-vos que s'escolin a l'espai de noms del workflow.

---

## 3. Compileu i instal·leu el vostre plugin

Compileu i instal·leu el plugin:

```bash
make install
```

!!! tip "Si la compilació falla"

    Llegiu el missatge d'error amb atenció; normalment inclou un número de línia i descriu el problema.
    Les causes més habituals són errors de sintaxi (parèntesi o cometa que falta), noms de classe mal escrits i incompatibilitats de tipus.
    Si us quedeu encallats, compareu el vostre codi caràcter per caràcter amb els exemples.

---

## 4. Utilitzeu la vostra funció en un workflow

El plugin està compilat i instal·lat.
El pas següent és utilitzar `reverseGreeting` en un workflow per verificar que funciona de cap a cap.

Torneu al directori del pipeline:

```bash
cd ..
```

Editeu `greet.nf` per importar i utilitzar `reverseGreeting`:

=== "Després"

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

=== "Abans"

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

Executeu el pipeline:

```bash
nextflow run greet.nf
```

??? example "Output"

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

La vostra primera funció de plugin personalitzada funciona en un workflow real.
El mateix patró `include { ... } from 'plugin/...'` que heu utilitzat amb nf-hello i nf-schema a la Part 1 funciona amb el vostre propi plugin.

---

## 5. Afegiu decorateGreeting

Un plugin pot proporcionar múltiples funcions.
Afegiu-ne una segona que embolcalli una salutació amb marcadors decoratius; la fareu configurable a la Part 6.

Editeu `GreetingExtension.groovy` per afegir `decorateGreeting` després de `reverseGreeting`, abans de la clau de tancament de la classe:

=== "Després"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverteix una cadena de salutació
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Decora una salutació amb marcadors festius
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Interpolació de cadenes de Groovy: `#!groovy ${...}` insereix el valor de la variable dins la cadena.

=== "Abans"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverteix una cadena de salutació
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Aquesta funció utilitza la interpolació de cadenes de Groovy (`"*** ${greeting} ***"`) per inserir la variable de salutació dins d'una cadena.

Compileu, instal·leu i actualitzeu el workflow:

```bash
cd nf-greeting && make install && cd ..
```

Actualitzeu `greet.nf` per importar i utilitzar també `decorateGreeting`:

=== "Després"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Importa les funcions personalitzades del nostre plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Utilitza la nostra funció de plugin personalitzada per decorar la salutació
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Demostra l'ús de la funció reverseGreeting
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Múltiples funcions del mateix plugin necessiten declaracions `include` separades.
    2. Les funcions del plugin també funcionen dins dels blocs `script:` dels processos.

=== "Abans"

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

??? example "Output"

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

Les funcions del plugin funcionen tant en scripts de processos (com `decorateGreeting` dins de `SAY_HELLO`) com en operacions de workflow (com `reverseGreeting` en un `map`).

---

## Conclusió

Heu après que:

- Les funcions es defineixen amb l'anotació `@Function` en subclasses de `PluginExtensionPoint`.
- Les funcions del plugin importades amb `include` funcionen de manera idèntica tant si provenen del vostre propi plugin com d'un d'existent.
- Les funcions del plugin funcionen tant en scripts de processos com en operacions de workflow.

---

## Què segueix?

Les vostres funcions funcionen, però fins ara només ho heu verificat executant el pipeline complet i comprovant la sortida visualment.
Aquest enfocament no escala: a mesura que afegiu més funcions, necessiteu una manera més ràpida de comprovar que cadascuna es comporta correctament, especialment després de fer canvis.
La secció següent introdueix les proves unitàries, que us permeten verificar funcions individuals de manera automàtica sense executar un pipeline.

[Continueu a la Part 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
