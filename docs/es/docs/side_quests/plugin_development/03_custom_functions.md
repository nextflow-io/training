# Parte 3: Funciones Personalizadas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Al final de esta sección, tendrá funciones personalizadas en su plugin, compiladas e instaladas localmente, ejecutándose en un workflow real.

!!! tip "¿Comenzando desde aquí?"

    Si se une en esta parte, copie la solución de la Parte 2 para usarla como punto de partida:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Ver lo que generó la plantilla

Antes de escribir sus propias funciones, observe la función de ejemplo que creó la plantilla para entender el patrón.

Cambie al directorio del plugin:

```bash
cd nf-greeting
```

La plantilla creó un archivo llamado `GreetingExtension.groovy` donde se definen las funciones del plugin.
Ábralo para ver el punto de partida:

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
 * Implementa una función personalizada que puede ser importada por
 * scripts de Nextflow.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Saluda al destinatario indicado.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. La clase sobre la que se construye su extensión. Nextflow requiere esto para reconocer sus funciones.
2. Se llama cuando el plugin se carga; úselo para la inicialización
3. Hace que este método sea invocable desde workflows mediante `include`

La plantilla incluye una función de ejemplo `sayHello`.
La anotación `@Function` es lo que hace que un método sea invocable desde workflows de Nextflow.
Sin ella, el método existe únicamente dentro del código del plugin.

En Groovy (y Java), los métodos declaran qué tipo retornan y qué tipos tienen sus parámetros.
Por ejemplo, `String reverseGreeting(String greeting)` declara un método que recibe un parámetro `String` y retorna un `String`.
La palabra clave `void` significa que el método no retorna nada, como ocurre con `sayHello` arriba.
Esto es diferente de Python o R, donde los tipos no necesitan declararse explícitamente.

---

## 2. Reemplazar sayHello con reverseGreeting

La función `sayHello` de la plantilla es un marcador de posición.
Reemplácela con su propia función para ver el ciclo completo de escritura, compilación y uso de una función de plugin.

Edite `src/main/groovy/training/plugin/GreetingExtension.groovy` para reemplazar el método `sayHello`:

=== "Después"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Invierte una cadena de saludo
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Hace que el método sea invocable desde workflows de Nextflow
    2. Recibe un String, retorna un String
    3. Método de inversión de cadenas integrado en Groovy

=== "Antes"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implementa una función personalizada que puede ser importada por
     * scripts de Nextflow.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Saluda al destinatario indicado.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Partes clave de esta función:

- **`@Function`**: Hace que el método sea invocable desde workflows de Nextflow
- **`String reverseGreeting(String greeting)`**: Recibe un String, retorna un String
- **`greeting.reverse()`**: Método de inversión de cadenas integrado en Groovy

!!! tip "Métodos públicos y privados"

    Los métodos sin `@Function` no se exponen a los workflows de Nextflow.
    Puede agregar métodos auxiliares a su clase sin preocuparse de que se filtren al espacio de nombres del workflow.

---

## 3. Compilar e instalar su plugin

Compile e instale el plugin:

```bash
make install
```

!!! tip "Si la compilación falla"

    Lea el mensaje de error con atención; generalmente incluye un número de línea y describe el problema.
    Las causas más comunes son errores de sintaxis (corchete o comilla faltante), nombres de clase mal escritos y tipos incompatibles.
    Si está atascado, compare su código carácter por carácter con los ejemplos.

---

## 4. Usar su función en un workflow

El plugin está compilado e instalado.
El siguiente paso es usar `reverseGreeting` en un workflow para verificar que funciona de extremo a extremo.

Regrese al directorio del pipeline:

```bash
cd ..
```

Edite `greet.nf` para importar y usar `reverseGreeting`:

=== "Después"

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

=== "Antes"

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

Ejecute el pipeline:

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

Su primera función de plugin personalizada está funcionando en un workflow real.
El mismo patrón `include { ... } from 'plugin/...'` que usó con nf-hello y nf-schema en la Parte 1 funciona con su propio plugin.

---

## 5. Agregar decorateGreeting

Un plugin puede proveer múltiples funciones.
Agregue una segunda que envuelva un saludo con marcadores decorativos; la hará configurable en la Parte 6.

Edite `GreetingExtension.groovy` para agregar `decorateGreeting` después de `reverseGreeting`, antes de la llave de cierre de la clase:

=== "Después"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Invierte una cadena de saludo
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Decora un saludo con marcadores festivos
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Interpolación de cadenas en Groovy: `#!groovy ${...}` inserta el valor de la variable en la cadena

=== "Antes"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Invierte una cadena de saludo
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Esta función usa la interpolación de cadenas de Groovy (`"*** ${greeting} ***"`) para insertar la variable del saludo dentro de una cadena.

Compile, instale y actualice el workflow:

```bash
cd nf-greeting && make install && cd ..
```

Actualice `greet.nf` para también importar y usar `decorateGreeting`:

=== "Después"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Importar funciones personalizadas de nuestro plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usar nuestra función de plugin personalizada para decorar el saludo
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Demostrar el uso de la función reverseGreeting
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Las múltiples funciones del mismo plugin necesitan declaraciones `include` separadas
    2. Las funciones del plugin también funcionan dentro de bloques `script:` de procesos

=== "Antes"

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

Las funciones del plugin funcionan tanto en scripts de procesos (como `decorateGreeting` dentro de `SAY_HELLO`) como en operaciones de workflow (como `reverseGreeting` en un `map`).

---

## Conclusión

Aprendió que:

- Las funciones se definen con la anotación `@Function` en subclases de `PluginExtensionPoint`
- Las funciones del plugin importadas con `include` funcionan de manera idéntica ya sea que provengan de su propio plugin o de uno existente
- Las funciones del plugin funcionan tanto en scripts de procesos como en operaciones de workflow

---

## ¿Qué sigue?

Sus funciones funcionan, pero hasta ahora solo lo ha verificado ejecutando el pipeline completo y revisando la salida visualmente.
Ese enfoque no escala: a medida que agregue más funciones, necesitará una forma más rápida de verificar que cada una se comporta correctamente, especialmente después de realizar cambios.
La siguiente sección presenta las pruebas unitarias, que le permiten verificar funciones individuales automáticamente sin ejecutar un pipeline.

[Continuar a la Parte 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
