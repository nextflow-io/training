# Parte 5: Observadores de Traza

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Los observadores de traza permiten que su plugin responda a eventos del workflow, como la finalización de una tarea, la publicación de un archivo o la conclusión del pipeline.
Esto habilita casos de uso como reportes personalizados, notificaciones de Slack, recopilación de métricas o integración con sistemas de monitoreo externos.
En esta sección, construirá un observador que cuenta las tareas completadas e imprime un resumen.

!!! tip "¿Empieza desde aquí?"

    Si se une en esta parte, copie la solución de la Parte 4 para usarla como punto de partida:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Comprender el observador de traza existente

El mensaje "Pipeline is starting!" que apareció al ejecutar el pipeline proviene de la clase `GreetingObserver` en su plugin.

Examine el código del observador:

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
 * Implementa un observador que permite ejecutar lógica personalizada
 * en los eventos de ejecución de Nextflow.
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

1. Interfaz para conectarse a los eventos del ciclo de vida del workflow
2. Se llama cuando el workflow inicia; recibe la sesión para acceder a la configuración
3. Se llama cuando el workflow finaliza exitosamente

Hay dos cosas a destacar aquí:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` es una interfaz definida por Nextflow. Si su clase implementa esta interfaz, Nextflow puede conectarse a ella y llamar a sus métodos cuando ocurran eventos.
2. **`@Override`**: La interfaz `TraceObserver` define métodos como `onFlowCreate` y `onFlowComplete`. Cuando escribe métodos con estos nombres y agrega la anotación `@Override`, Nextflow los llama en el momento apropiado. Los métodos que no sobreescriba son ignorados.

El conjunto completo de eventos del ciclo de vida a los que puede conectarse al momento de escribir esto son:

| Método              | Cuándo se llama                 |
| ------------------- | ------------------------------- |
| `onFlowCreate`      | El workflow inicia              |
| `onFlowComplete`    | El workflow finaliza            |
| `onProcessStart`    | Una tarea comienza su ejecución |
| `onProcessComplete` | Una tarea finaliza              |
| `onProcessCached`   | Se reutiliza una tarea en caché |
| `onFilePublish`     | Se publica un archivo           |

Para una lista completa, consulte la [interfaz TraceObserver](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) en el código fuente de Nextflow.

---

## 2. Agregar un observador contador de tareas

El objetivo es construir un observador que cuente las tareas completadas e imprima un resumen al final.
Agregar un nuevo observador a un plugin requiere dos cosas: escribir la clase del observador y registrarla en la fábrica para que Nextflow la cargue.

### 2.1. Crear un observador mínimo

Cree un nuevo archivo:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Comience con el observador más simple posible que imprime un mensaje cuando cualquier tarea se completa:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observador que responde a la finalización de tareas
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Importar las clases requeridas: `TraceObserver`, `TaskHandler` y `TraceRecord`
2. Crear una clase que `implements TraceObserver`
3. Sobreescribir `onProcessComplete` para ejecutar código cuando una tarea finaliza

Esto es lo mínimo necesario:

- Importar las clases requeridas (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Crear una clase que `implements TraceObserver`
- Sobreescribir `onProcessComplete` para hacer algo cuando una tarea finaliza

### 2.2. Registrar el observador

`GreetingFactory` crea los observadores.
Examínelo:

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

Edite `GreetingFactory.groovy` para agregar el nuevo observador:

=== "Después"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Antes"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Sintaxis de listas en Groovy"

    Hemos reemplazado el estilo Java `List.<TraceObserver>of(...)` con el literal de lista más simple de Groovy `[...]`.
    Ambos retornan una `Collection`, pero la sintaxis de Groovy es más legible al agregar múltiples elementos.

### 2.3. Compilar, instalar y probar

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "¿Por qué `-ansi-log false`?"

    Por defecto, la visualización de progreso ANSI de Nextflow sobreescribe las líneas anteriores para mostrar una vista limpia y actualizada del progreso.
    Esto significa que solo vería el recuento *final* de tareas, no los mensajes intermedios.

    Usar `-ansi-log false` deshabilita este comportamiento y muestra toda la salida de forma secuencial, lo cual es esencial al probar observadores que imprimen mensajes durante la ejecución.

Debería ver "✓ Task completed!" impreso cinco veces (una por tarea), intercalado con la salida existente del pipeline:

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

El observador está funcionando.
Cada vez que una tarea finaliza, Nextflow llama a `onProcessComplete`, y nuestra implementación imprime un mensaje.

??? exercise "Personalizar el mensaje"

    Intente cambiar el mensaje en `onProcessComplete` por uno propio, recompile y vuelva a ejecutar.
    Esto confirma que el ciclo completo de edición-compilación-ejecución funciona para los observadores.

### 2.4. Agregar lógica de conteo

El observador mínimo demuestra que el hook funciona, pero no registra nada.

Una clase puede contener variables (llamadas campos o variables de instancia) que persisten durante el tiempo de vida del objeto.
Esto significa que un observador puede acumular estado a través de múltiples eventos durante una ejecución del pipeline.

La siguiente versión agrega una variable contadora (`taskCount`) que comienza en cero.
Cada vez que una tarea se completa, el contador aumenta en uno.
Cuando el workflow completo finaliza, el observador imprime el total final.

Actualice `TaskCounterObserver.groovy` con los cambios resaltados:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observador que cuenta las tareas completadas
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

1. `taskCount` es una variable que pertenece al objeto observador. Mantiene su valor entre llamadas a métodos, por lo que puede acumular un conteo a lo largo de toda la ejecución del workflow. `private` significa que solo esta clase puede acceder a ella.
2. `taskCount++` agrega uno al contador. Esta línea se ejecuta cada vez que una tarea se completa, por lo que el conteo crece a medida que el workflow avanza.
3. `onFlowComplete` es un segundo hook del ciclo de vida. Se ejecuta una vez cuando el workflow finaliza, lo que lo convierte en un buen lugar para imprimir un resumen.

En resumen:

- `taskCount` persiste entre llamadas a métodos, acumulando un conteo durante toda la ejecución
- `onProcessComplete` incrementa el contador e imprime el total acumulado cada vez que una tarea finaliza
- `onFlowComplete` se ejecuta una vez al final, imprimiendo el conteo final

Recompile y pruebe:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Salida"

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

    Los mensajes del contador están intercalados con los envíos de tareas porque los observadores se ejecutan a medida que las tareas se completan.

---

## 3. Rastrear archivos publicados

El observador también puede responder cuando se publican archivos.
El método `onFilePublish` recibe las rutas de destino y origen, que puede usar para registrar, validar o procesar las salidas publicadas.

### 3.1. Agregar un directorio de publicación

Primero, actualice `greet.nf` para que el proceso `SAY_HELLO` publique sus archivos de salida:

=== "Después"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usar nuestra función personalizada del plugin para decorar el saludo
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usar nuestra función personalizada del plugin para decorar el saludo
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Agregar el método onFilePublish

Agregue un método `onFilePublish` y la importación requerida a `TaskCounterObserver.groovy`:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observador que cuenta las tareas completadas
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

### 3.3. Compilar y probar

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Debería ver mensajes "Published:" para cada archivo de salida junto con la salida del contador de tareas:

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

El método `onFilePublish` se activa cada vez que Nextflow publica un archivo en el directorio `results`.
Este patrón es útil para construir registros de auditoría, activar acciones posteriores o validar las salidas a medida que se producen.

---

## Conclusión

Aprendió que:

- Los observadores de traza se conectan a los eventos del ciclo de vida del workflow como `onFlowCreate`, `onProcessComplete`, `onFilePublish` y `onFlowComplete`
- Se crean observadores implementando `TraceObserver` y registrándolos en una fábrica
- Los observadores pueden contener variables de instancia para acumular estado a través de los eventos
- Los observadores son útiles para registro personalizado, recopilación de métricas, notificaciones y reportes

---

## ¿Qué sigue?

El contador de tareas funciona, pero siempre está activo.
En un plugin real, los usuarios deberían poder habilitar o deshabilitar funcionalidades, o ajustar el comportamiento, desde `nextflow.config` sin necesidad de editar el código fuente del plugin.
La siguiente sección muestra cómo hacer que su observador sea configurable y cómo compartir su plugin terminado con otros.

[Continuar a la Parte 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
