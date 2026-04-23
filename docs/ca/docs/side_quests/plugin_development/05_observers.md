# Part 5: Trace Observers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Els trace observers permeten que el vostre plugin respongui a esdeveniments del workflow, com ara la finalització d'una tasca, la publicació d'un fitxer o la finalització del pipeline.
Això permet casos d'ús com ara informes personalitzats, notificacions de Slack, recollida de mètriques o integració amb sistemes de monitoratge externs.
En aquesta secció, construireu un observer que compta les tasques completades i imprimeix un resum.

!!! tip "Comenceu des d'aquí?"

    Si us incorporeu en aquesta part, copieu la solució de la Part 4 per utilitzar-la com a punt de partida:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Comprendre l'observer existent

El missatge "Pipeline is starting!" que apareixia quan executàveu el pipeline provenia de la classe `GreetingObserver` del vostre plugin.

Mireu el codi de l'observer:

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
 * Implementa un observer que permet implementar lògica
 * personalitzada en els esdeveniments d'execució de Nextflow.
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

1. Interfície per connectar-se als esdeveniments del cicle de vida del workflow
2. S'invoca quan el workflow s'inicia; rep la sessió per accedir a la configuració
3. S'invoca quan el workflow finalitza correctament

Hi ha dues coses a destacar:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` és una interfície definida per Nextflow. Si la vostra classe implementa aquesta interfície, Nextflow pot connectar-s'hi i cridar els vostres mètodes quan es produeixen esdeveniments.
2. **`@Override`**: La interfície `TraceObserver` defineix mètodes com `onFlowCreate` i `onFlowComplete`. Quan escriviu mètodes amb aquests noms i afegiu l'anotació `@Override`, Nextflow els crida en el moment adequat. Qualsevol mètode que no sobreescriviu s'ignora.

El conjunt complet d'esdeveniments del cicle de vida als quals podeu connectar-vos en el moment d'escriure aquest document és:

| Mètode              | Quan s'invoca                    |
| ------------------- | -------------------------------- |
| `onFlowCreate`      | El workflow s'inicia             |
| `onFlowComplete`    | El workflow finalitza            |
| `onProcessStart`    | Una tasca comença l'execució     |
| `onProcessComplete` | Una tasca finalitza              |
| `onProcessCached`   | Es reutilitza una tasca en cache |
| `onFilePublish`     | Es publica un fitxer             |

Per a una llista completa, consulteu la [interfície TraceObserver](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) al codi font de Nextflow.

---

## 2. Afegir un observer comptador de tasques

L'objectiu és construir un observer que compti les tasques completades i imprimeixi un resum al final.
Afegir un nou observer a un plugin requereix dues coses: escriure la classe de l'observer i registrar-la a la factory perquè Nextflow la carregui.

### 2.1. Crear un observer mínim

Creeu un fitxer nou:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Comenceu amb l'observer més senzill possible que imprimeixi un missatge quan qualsevol tasca es completi:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer que respon a la finalització de tasques
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Importeu les classes necessàries: `TraceObserver`, `TaskHandler` i `TraceRecord`
2. Creeu una classe que `implements TraceObserver`
3. Sobreescriviu `onProcessComplete` per executar codi quan una tasca finalitza

Això és el mínim necessari:

- Importar les classes necessàries (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Crear una classe que `implements TraceObserver`
- Sobreescriure `onProcessComplete` per fer alguna cosa quan una tasca finalitza

### 2.2. Registrar l'observer

La `GreetingFactory` crea observers.
Mireu-la:

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

Editeu `GreetingFactory.groovy` per afegir el nou observer:

=== "Després"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Abans"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Sintaxi de llistes en Groovy"

    Hem substituït el `List.<TraceObserver>of(...)` d'estil Java per la sintaxi de llista literal `[...]` més senzilla de Groovy.
    Tots dos retornen una `Collection`, però la sintaxi de Groovy és més llegible quan s'afegeixen múltiples elements.

### 2.3. Construir, instal·lar i provar

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Per què `-ansi-log false`?"

    Per defecte, la visualització de progrés ANSI de Nextflow sobreescriu les línies anteriors per mostrar una vista neta i actualitzada del progrés.
    Això significa que només veuríeu el *recompte final* de tasques, no els missatges intermedis.

    Utilitzar `-ansi-log false` desactiva aquest comportament i mostra tota la sortida de manera seqüencial, cosa que és essencial quan es proven observers que imprimeixen missatges durant l'execució.

Hauríeu de veure "✓ Task completed!" imprès cinc vegades (una per tasca), intercalat amb la sortida existent del pipeline:

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

L'observer funciona.
Cada vegada que una tasca finalitza, Nextflow crida `onProcessComplete`, i la nostra implementació imprimeix un missatge.

??? exercise "Personalitzeu el missatge"

    Proveu de canviar el missatge a `onProcessComplete` per alguna cosa pròpia, reconstruïu i torneu a executar.
    Això confirma que el cicle complet d'edició-construcció-execució funciona per als observers.

### 2.4. Afegir lògica de comptatge

L'observer mínim demostra que el hook funciona, però no fa cap seguiment.

Una classe pot contenir variables (anomenades camps o variables d'instància) que persisteixen durant tota la vida de l'objecte.
Això significa que un observer pot acumular estat a través de múltiples esdeveniments durant una execució del pipeline.

La versió següent afegeix una variable comptadora (`taskCount`) que comença a zero.
Cada vegada que una tasca es completa, el comptador augmenta en un.
Quan tot el workflow finalitza, l'observer imprimeix el total final.

Actualitzeu `TaskCounterObserver.groovy` amb els canvis destacats:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer que compta les tasques completades
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

1. `taskCount` és una variable que pertany a l'objecte observer. Manté el seu valor entre crides a mètodes, de manera que pot acumular un recompte durant tota l'execució del workflow. `private` significa que només aquesta classe hi pot accedir.
2. `taskCount++` afegeix un al comptador. Aquesta línia s'executa cada vegada que una tasca es completa, de manera que el recompte creix a mesura que el workflow avança.
3. `onFlowComplete` és un segon hook del cicle de vida. S'executa una vegada quan el workflow finalitza, cosa que el converteix en un bon lloc per imprimir un resum.

En resum:

- `taskCount` persisteix entre crides a mètodes, acumulant un recompte durant tota l'execució
- `onProcessComplete` incrementa el comptador i imprimeix el total acumulat cada vegada que una tasca finalitza
- `onFlowComplete` s'executa una vegada al final, imprimint el recompte final

Reconstruïu i proveu:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Sortida"

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

    Els missatges del comptador s'intercalen amb les submissions de tasques perquè els observers s'executen a mesura que les tasques es completen.

---

## 3. Fer el seguiment dels fitxers publicats

L'observer també pot respondre quan es publiquen fitxers.
El mètode `onFilePublish` rep els camins de destinació i d'origen, que podeu utilitzar per registrar, validar o processar les sortides publicades.

### 3.1. Afegir un directori de publicació

Primer, actualitzeu `greet.nf` perquè el procés `SAY_HELLO` publiqui els seus fitxers de sortida:

=== "Després"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Utilitzeu la nostra funció de plugin personalitzada per decorar la salutació
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Abans"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Utilitzeu la nostra funció de plugin personalitzada per decorar la salutació
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Afegir el mètode onFilePublish

Afegiu un mètode `onFilePublish` i la importació necessària a `TaskCounterObserver.groovy`:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer que compta les tasques completades
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

### 3.3. Construir i provar

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Hauríeu de veure missatges "Published:" per a cada fitxer de sortida juntament amb la sortida del comptador de tasques:

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

El mètode `onFilePublish` s'activa cada vegada que Nextflow publica un fitxer al directori `results`.
Aquest patró és útil per construir registres d'auditoria, activar accions posteriors o validar les sortides a mesura que es produeixen.

---

## Conclusió

Heu après que:

- Els trace observers es connecten als esdeveniments del cicle de vida del workflow com `onFlowCreate`, `onProcessComplete`, `onFilePublish` i `onFlowComplete`
- Es creen observers implementant `TraceObserver` i registrant-los en una Factory
- Els observers poden contenir variables d'instància per acumular estat a través dels esdeveniments
- Els observers són útils per a registres personalitzats, recollida de mètriques, notificacions i informes

---

## Què segueix?

El comptador de tasques funciona, però sempre està actiu.
En un plugin real, els usuaris haurien de poder activar o desactivar funcionalitats, o ajustar el comportament, des de `nextflow.config` sense haver d'editar el codi font del plugin.
La secció següent mostra com fer que el vostre observer sigui configurable i com compartir el vostre plugin acabat amb altres persones.

[Continueu a la Part 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
