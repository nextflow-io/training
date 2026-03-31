# Patrones Esenciales de Scripting en Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow es un lenguaje de programación que se ejecuta en la Máquina Virtual de Java. Aunque Nextflow está construido sobre [Groovy](http://groovy-lang.org/) y comparte gran parte de su sintaxis, Nextflow es más que simplemente "Groovy con extensiones" -- es un lenguaje independiente con una [sintaxis](https://nextflow.io/docs/latest/reference/syntax.html) y una [biblioteca estándar](https://nextflow.io/docs/latest/reference/stdlib.html) completamente especificadas.

Se puede escribir mucho código en Nextflow sin ir más allá de la sintaxis básica para variables, maps y listas. La mayoría de los tutoriales de Nextflow se centran en la orquestación del workflow (canales, procesos y flujo de datos), y se puede avanzar sorprendentemente lejos con solo eso.

Sin embargo, cuando se necesita manipular datos, analizar nombres de archivos complejos, implementar lógica condicional o construir workflows robustos para producción, es útil pensar en dos aspectos distintos del código: **dataflow** (canales, operadores, procesos y workflows) y **scripting** (el código dentro de closures, funciones y scripts de procesos). Aunque esta distinción es algo arbitraria —todo es código Nextflow— proporciona un modelo mental útil para entender cuándo se está orquestando el pipeline versus cuándo se están manipulando datos. Dominar ambos aspectos mejora considerablemente la capacidad de escribir workflows claros y mantenibles.

### Objetivos de aprendizaje

Esta misión secundaria le lleva en un recorrido práctico desde conceptos básicos hasta patrones listos para producción.
Transformaremos un workflow simple de lectura de CSV en un sofisticado pipeline de bioinformática, evolucionándolo paso a paso a través de desafíos realistas:

- **Comprensión de límites:** Distinguir entre operaciones de dataflow y scripting, y entender cómo trabajan juntos
- **Manipulación de datos:** Extraer, transformar y crear subconjuntos de maps y colecciones usando operadores potentes
- **Procesamiento de strings:** Analizar esquemas complejos de nombres de archivos con patrones regex y dominar la interpolación de variables
- **Funciones reutilizables:** Extraer lógica compleja en funciones con nombre para workflows más limpios y mantenibles
- **Lógica dinámica:** Construir procesos que se adapten a diferentes tipos de entrada y usar closures para la asignación dinámica de recursos
- **Enrutamiento condicional:** Enrutar muestras de forma inteligente a través de diferentes procesos según sus características de metadatos
- **Operaciones seguras:** Manejar datos faltantes de forma elegante con operadores null-safe y validar entradas con mensajes de error claros
- **Manejadores basados en configuración:** Usar manejadores de eventos del workflow para registro, notificaciones y gestión del ciclo de vida

### Requisitos previos

Antes de comenzar esta misión secundaria, debe:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirse cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, trabajo con archivos, metadatos)
- Tener familiaridad básica con construcciones de programación comunes (variables, maps, listas)

Este tutorial explicará los conceptos de programación a medida que los encontremos, por lo que no se necesita experiencia extensa en programación.
Comenzaremos con conceptos fundamentales y avanzaremos hasta patrones más complejos.

---

## 0. Primeros pasos

#### Abrir el codespace de capacitación

Si aún no lo ha hecho, asegúrese de abrir el entorno de capacitación como se describe en la [Configuración del entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos de este tutorial.

```bash
cd side-quests/essential_scripting_patterns
```

#### Revisar los materiales

Encontrará un archivo de workflow principal y un directorio `data` que contiene archivos de datos de ejemplo.

```console title="Directory contents"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Nuestro CSV de muestras contiene información sobre muestras biológicas que necesitan diferente procesamiento según sus características:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Usaremos este conjunto de datos realista para explorar técnicas de programación prácticas que encontrará en workflows reales de bioinformática.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo correctamente
<!-- - [ ] I understand the assignment -->

Si puede marcar todas las casillas, está listo para continuar.

---

## 1. Dataflow vs Scripting: Comprendiendo los Límites

### 1.1. Identificando Qué es Qué

Al escribir workflows en Nextflow, es importante distinguir entre **dataflow** (cómo se mueven los datos a través de canales y procesos) y **scripting** (el código que manipula datos y toma decisiones). Construyamos un workflow que demuestre cómo trabajan juntos.

#### 1.1.1. Workflow Básico en Nextflow

Comience con un workflow simple que solo lee el archivo CSV (ya lo hemos hecho por usted en `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

El bloque `workflow` define la estructura de nuestro pipeline, mientras que `channel.fromPath()` crea un canal a partir de una ruta de archivo. El operador `.splitCsv()` procesa el archivo CSV y convierte cada fila en una estructura de datos map.

Ejecute este workflow para ver los datos CSV sin procesar:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Agregando el Operador Map

Ahora vamos a agregar scripting para transformar los datos, usando el operador `.map()` con el que probablemente ya esté familiarizado. Este operador toma un 'closure' donde podemos escribir código para transformar cada elemento.

!!! note "Nota"

    Un **closure** es un bloque de código que puede pasarse y ejecutarse más tarde. Piense en él como una función que define en línea. Los closures se escriben con llaves `{ }` y pueden tomar parámetros. Son fundamentales para el funcionamiento de los operadores de Nextflow y, si lleva un tiempo escribiendo Nextflow, es posible que ya los haya estado usando sin darse cuenta.

Así es como se ve esa operación map:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Este es nuestro primer **closure** -- una función anónima que puede pasarse como argumento (similar a las lambdas en Python o las funciones flecha en JavaScript). Los closures son esenciales para trabajar con operadores de Nextflow.

El closure `{ row -> return row }` toma un parámetro `row` (podría tener cualquier nombre: `item`, `sample`, etc.).

Cuando el operador `.map()` procesa cada elemento del canal, pasa ese elemento a su closure. Aquí, `row` contiene una fila del CSV a la vez.

Aplique este cambio y ejecute el workflow:

```bash
nextflow run main.nf
```

Verá la misma salida que antes, porque simplemente estamos devolviendo la entrada sin cambios. Esto confirma que el operador map está funcionando correctamente. Ahora comencemos a transformar los datos.

#### 1.1.3. Creando una Estructura de Datos Map

Ahora vamos a escribir lógica de **scripting** dentro de nuestro closure para transformar cada fila de datos. Aquí es donde procesamos elementos de datos individuales en lugar de orquestar el flujo de datos.

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting para transformación de datos
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

El map `sample_meta` es una estructura de datos clave-valor (como los diccionarios en Python, los objetos en JavaScript o los hashes en Ruby) que almacena información relacionada: ID de muestra, organismo, tipo de tejido, profundidad de secuenciación y puntuación de calidad.

Usamos métodos de manipulación de strings como `.toLowerCase()` y `.replaceAll()` para limpiar nuestros datos, y métodos de conversión de tipos como `.toInteger()` y `.toDouble()` para convertir los datos de string del CSV en los tipos numéricos apropiados.

Aplique este cambio y ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Agregando Lógica Condicional

Ahora agreguemos más scripting -- esta vez usando un operador ternario para tomar decisiones basadas en valores de datos.

Realice el siguiente cambio:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

El operador ternario es una forma abreviada de una sentencia if/else que sigue el patrón `condición ? valor_si_verdadero : valor_si_falso`. Esta línea significa: "Si la calidad es mayor que 40, usar 'high', de lo contrario usar 'normal'". Su primo, el **operador Elvis** (`?:`), proporciona valores predeterminados cuando algo es null o está vacío -- exploraremos ese patrón más adelante en este tutorial.

El operador de adición de maps `+` crea un **nuevo map** en lugar de modificar el existente. Esta línea crea un nuevo map que contiene todos los pares clave-valor de `sample_meta` más la nueva clave `priority`.

!!! Note "Nota"

    Nunca modifique los maps que se pasan a los closures -- siempre cree nuevos usando `+` (por ejemplo). En Nextflow, los mismos datos a menudo fluyen a través de múltiples operaciones simultáneamente. Modificar un map en el lugar puede causar efectos secundarios impredecibles cuando otras operaciones hacen referencia a ese mismo objeto. Crear nuevos maps garantiza que cada operación tenga su propia copia limpia.

Ejecute el workflow modificado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Hemos agregado con éxito lógica condicional para enriquecer nuestros metadatos con un nivel de prioridad basado en las puntuaciones de calidad.

#### 1.1.5. Creando Subconjuntos de Maps con `.subMap()`

Mientras que el operador `+` agrega claves a un map, a veces necesita hacer lo contrario -- extraer solo claves específicas. El método `.subMap()` es perfecto para esto.

Agreguemos una línea para crear una versión simplificada de nuestros metadatos que solo contenga campos de identificación:

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting para transformación de datos
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting para transformación de datos
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Ejecute el workflow modificado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Esto muestra tanto los metadatos completos mostrados por la operación `view()` como el subconjunto extraído que imprimimos con `println`.

El método `.subMap()` toma una lista de claves y devuelve un nuevo map que contiene solo esas claves. Si una clave no existe en el map original, simplemente no se incluye en el resultado.

Esto es particularmente útil cuando necesita crear diferentes versiones de metadatos para diferentes procesos -- algunos pueden necesitar metadatos completos mientras que otros solo necesitan campos de identificación mínimos.

Ahora elimine esas sentencias println para restaurar su workflow a su estado anterior, ya que no las necesitaremos en adelante.

!!! tip "Consejo"

    **Resumen de Operaciones con Maps**

    - **Agregar claves**: `map1 + [new_key: value]` - Crea un nuevo map con claves adicionales
    - **Extraer claves**: `map1.subMap(['key1', 'key2'])` - Crea un nuevo map con solo las claves especificadas
    - **Ambas operaciones crean nuevos maps** - Los maps originales permanecen sin cambios

#### 1.1.6. Combinando Maps y Devolviendo Resultados

Hasta ahora, solo hemos devuelto lo que la comunidad de Nextflow llama el 'meta map', e hemos ignorado los archivos a los que se refieren esos metadatos. Pero si está escribiendo workflows en Nextflow, probablemente quiera hacer algo con esos archivos.

Generemos una estructura de canal que comprenda una tupla de 2 elementos: el map de metadatos enriquecido y la ruta de archivo correspondiente. Este es un patrón común en Nextflow para pasar datos a los procesos.

=== "Después"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Aplique este cambio y ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Esta estructura de tupla `[meta, file]` es un patrón común en Nextflow para pasar tanto metadatos como archivos asociados a los procesos.

!!! note "Nota"

    **Maps y Metadatos**: Los maps son fundamentales para trabajar con metadatos en Nextflow. Para una explicación más detallada sobre cómo trabajar con maps de metadatos, consulte la misión secundaria [Working with metadata](../metadata/).

Nuestro workflow demuestra el patrón central: las **operaciones de dataflow** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orquestan cómo se mueven los datos a través del pipeline, mientras que el **scripting** (maps `[key: value]`, métodos de string, conversiones de tipos, operadores ternarios) dentro del closure `.map()` maneja la transformación de elementos de datos individuales.

### 1.2. Comprendiendo Diferentes Tipos: Canal vs Lista

Hasta ahora, todo bien, podemos distinguir entre operaciones de dataflow y scripting. ¿Pero qué pasa cuando el mismo nombre de método existe en ambos contextos?

Un ejemplo perfecto es el método `collect`, que existe tanto para tipos de canal como para tipos List en la biblioteca estándar de Nextflow. El método `collect()` en una List transforma cada elemento, mientras que el operador `collect()` en un canal agrupa todas las emisiones del canal en un canal de un solo elemento.

Demostremos esto con algunos datos de muestra, comenzando por refrescar lo que hace el operador `collect()` del canal. Revise `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - agrupa múltiples emisiones del canal en una sola
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Pasos:

- Definir una List de IDs de muestras
- Crear un canal con `fromList()` que emite cada ID de muestra por separado
- Imprimir cada elemento con `view()` a medida que fluye
- Reunir todos los elementos en una sola lista con el operador `collect()` del canal
- Imprimir el resultado recopilado (elemento único que contiene todos los IDs de muestras) con un segundo `view()`

Hemos cambiado la estructura del canal, pero no hemos cambiado los datos en sí.

Ejecute el workflow para confirmarlo:

```bash
nextflow run collect.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` devuelve una salida por cada emisión del canal, por lo que sabemos que esta única salida contiene los 3 elementos originales agrupados en una lista.

Ahora veamos el método `collect` en una List en acción. Modifique `collect.nf` para aplicar el método `collect` de la List a la lista original de IDs de muestras:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones del canal en una sola
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva la estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones del canal en una sola
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

En este nuevo fragmento:

- Definimos una nueva variable `formatted_ids` que usa el método `collect` de la List para transformar cada ID de muestra en la lista original
- Imprimimos el resultado usando `println`

Ejecute el workflow modificado:

```bash
nextflow run collect.nf
```

??? success "Salida del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Esta vez, NO hemos cambiado la estructura de los datos, todavía tenemos 3 elementos en la lista, pero SÍ hemos transformado cada elemento usando el método `collect` de la List para producir una nueva lista con valores modificados. Esto es similar a usar el operador `map` en un canal, pero opera sobre una estructura de datos List en lugar de un canal.

`collect` es un caso extremo que usamos aquí para ilustrar un punto. La lección clave es que al escribir workflows, siempre distinga entre **estructuras de datos** (Lists, Maps, etc.) y **canales** (construcciones de dataflow). Las operaciones pueden compartir nombres pero comportarse de manera completamente diferente dependiendo del tipo sobre el que se llaman.

### 1.3. El Operador Spread (`*.`) - Forma Abreviada para Extracción de Propiedades

Relacionado con el método `collect` de la List está el operador spread (`*.`), que proporciona una forma concisa de extraer propiedades de colecciones. Es esencialmente azúcar sintáctico para un patrón común de `collect`.

Agreguemos una demostración a nuestro archivo `collect.nf`:

=== "Después"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones del canal en una sola
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva la estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Operador spread - acceso conciso a propiedades
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Antes"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones del canal en una sola
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva la estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Ejecute el workflow actualizado:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Salida del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

El operador spread `*.` es una forma abreviada de un patrón común de collect:

```groovy
// Estas son equivalentes:
def ids = samples*.id
def ids = samples.collect { it.id }

// También funciona con llamadas a métodos:
def names = files*.getName()
def names = files.collect { it.getName() }
```

El operador spread es particularmente útil cuando necesita extraer una sola propiedad de una lista de objetos -- es más legible que escribir el closure completo de `collect`.

!!! tip "Consejo"

    **Cuándo Usar Spread vs Collect**

    - **Use spread (`*.`)** para acceso simple a propiedades: `samples*.id`, `files*.name`
    - **Use collect** para transformaciones o lógica compleja: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Conclusión

En esta sección, ha aprendido:

- **Dataflow vs scripting**: Los operadores de canal orquestan cómo fluyen los datos a través de su pipeline, mientras que el scripting transforma elementos de datos individuales
- **Comprensión de tipos**: El mismo nombre de método (como `collect`) puede comportarse de manera diferente dependiendo del tipo sobre el que se llama (Canal vs List)
- **El contexto importa**: Siempre sea consciente de si está trabajando con canales (dataflow) o estructuras de datos (scripting)

Comprender estos límites es esencial para depurar, documentar y escribir workflows mantenibles.

A continuación, profundizaremos en las capacidades de procesamiento de strings, que son esenciales para manejar datos del mundo real.

---

## 2. Procesamiento de Strings y Generación Dinámica de Scripts

Dominar el procesamiento de strings separa los workflows frágiles de los pipelines robustos. Esta sección cubre el análisis de nombres de archivos complejos, la generación dinámica de scripts y la interpolación de variables.

### 2.1. Coincidencia de Patrones y Expresiones Regulares

Los archivos de bioinformática a menudo tienen convenciones de nomenclatura complejas que codifican metadatos. Extraigamos esto automáticamente usando coincidencia de patrones con expresiones regulares.

Vamos a volver a nuestro workflow `main.nf` y agregar lógica de coincidencia de patrones para extraer información adicional de la muestra a partir de los nombres de archivos. Los archivos FASTQ en nuestro conjunto de datos siguen convenciones de nomenclatura estilo Illumina con nombres como `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Pueden parecer crípticos, pero en realidad codifican metadatos útiles como el ID de muestra, el número de carril y la dirección de lectura. Vamos a usar las capacidades de regex para analizar estos nombres.

Realice el siguiente cambio en su workflow `main.nf` existente:

=== "Después"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting para transformación de datos
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Scripting para transformación de datos
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Esto demuestra conceptos clave de **procesamiento de strings**:

1. **Literales de expresión regular** usando la sintaxis `~/patrón/` -- esto crea un patrón regex sin necesidad de escapar barras invertidas
2. **Coincidencia de patrones** con el operador `=~` -- esto intenta hacer coincidir un string con un patrón regex
3. **Objetos Matcher** que capturan grupos con `[0][1]`, `[0][2]`, etc. -- `[0]` se refiere a la coincidencia completa, `[1]`, `[2]`, etc. se refieren a los grupos capturados entre paréntesis

Analicemos el patrón regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Patrón              | Coincide con                           | Captura                                    |
| ------------------- | -------------------------------------- | ------------------------------------------ |
| `^(.+)`             | Nombre de muestra desde el inicio      | Grupo 1: nombre de muestra                 |
| `_S(\d+)`           | Número de muestra `_S1`, `_S2`, etc.   | Grupo 2: número de muestra                 |
| `_L(\d{3})`         | Número de carril `_L001`               | Grupo 3: carril (3 dígitos)                |
| `_(R[12])`          | Dirección de lectura `_R1` o `_R2`     | Grupo 4: dirección de lectura              |
| `_(\d{3})`          | Número de fragmento `_001`             | Grupo 5: fragmento (3 dígitos)             |
| `\.fastq(?:\.gz)?$` | Extensión `.fastq` o `.fastq.gz`       | No capturado (?:  es no capturante)        |

Esto analiza las convenciones de nomenclatura estilo Illumina para extraer metadatos automáticamente.

Ejecute el workflow modificado:

```bash title="Test pattern matching"
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Esto muestra los metadatos enriquecidos a partir de los nombres de archivos.

### 2.2. Generación Dinámica de Scripts en Procesos

Los bloques de script de los procesos son esencialmente strings multilínea que se pasan al shell. Puede usar **lógica condicional** (if/else, operadores ternarios) para generar dinámicamente diferentes strings de script basados en las características de la entrada. Esto es esencial para manejar tipos de entrada diversos -- como lecturas de secuenciación de extremo único vs extremo pareado -- sin duplicar definiciones de procesos.

Agreguemos un proceso a nuestro workflow que demuestre este patrón. Abra `modules/fastp.nf` y échele un vistazo:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

El proceso toma archivos FASTQ como entrada y ejecuta la herramienta `fastp` para recortar adaptadores y filtrar lecturas de baja calidad. Desafortunadamente, la persona que escribió este proceso no contempló las lecturas de extremo único que tenemos en nuestro conjunto de datos de ejemplo. Agreguémoslo a nuestro workflow y veamos qué sucede:

Primero, incluya el módulo en la primera línea de su workflow `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Luego modifique el bloque `workflow` para conectar el canal `ch_samples` al proceso `FASTP`:

=== "Después"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Ejecute este workflow modificado:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Puede ver que el proceso está intentando ejecutar `fastp` con un valor `null` para el segundo archivo de entrada, lo que está causando que falle. Esto se debe a que nuestro conjunto de datos contiene lecturas de extremo único, pero el proceso está codificado para esperar lecturas de extremo pareado (dos archivos de entrada a la vez).

Corrija esto agregando lógica condicional al bloque `script:` del proceso `FASTP`. Una sentencia if/else verifica el número de archivos de lectura y ajusta el comando en consecuencia.

=== "Después"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Detección simple de extremo único vs extremo pareado
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Ahora el workflow puede manejar tanto lecturas de extremo único como de extremo pareado de forma elegante. La lógica condicional verifica el número de archivos de entrada y construye el comando apropiado para `fastp`. Veamos si funciona:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

¡Se ve bien! Si verificamos los comandos reales que se ejecutaron (personalice según el hash de su tarea):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Podemos ver que Nextflow eligió correctamente el comando para lecturas de extremo único:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Otro uso común de la lógica de script dinámica se puede ver en [el módulo de Genómica de Nextflow for Science](../../nf4science/genomics/02_joint_calling). En ese módulo, el proceso GATK que se llama puede tomar múltiples archivos de entrada, pero cada uno debe ir precedido de `-V` para formar una línea de comando correcta. El proceso usa scripting para transformar una colección de archivos de entrada (`all_gvcfs`) en los argumentos de comando correctos:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Estos patrones de uso de scripting en los bloques de script de los procesos son extremadamente poderosos y pueden aplicarse en muchos escenarios -- desde el manejo de tipos de entrada variables hasta la construcción de argumentos de línea de comando complejos a partir de colecciones de archivos, haciendo que sus procesos sean verdaderamente adaptables a los diversos requisitos de los datos del mundo real.

### 2.3. Interpolación de Variables: Variables de Nextflow y del Shell

Los scripts de los procesos mezclan variables de Nextflow, variables del shell y sustituciones de comandos, cada una con diferente sintaxis de interpolación. Usar la sintaxis incorrecta causa errores. Exploremos esto con un proceso que crea un informe de procesamiento.

Eche un vistazo al archivo de módulo `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Este proceso escribe un informe simple con el ID de muestra y el nombre de archivo. Ahora ejecutémoslo para ver qué sucede cuando necesitamos mezclar diferentes tipos de variables.

Incluya el proceso en su `main.nf` y agréguelo al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

Ahora ejecute el workflow y verifique los informes generados en `results/reports/`. Deberían contener información básica sobre cada muestra.

<!-- TODO: add the run command -->

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

¿Pero qué pasa si queremos agregar información sobre cuándo y dónde ocurrió el procesamiento? Modifiquemos el proceso para usar variables del **shell** y algo de sustitución de comandos para incluir el usuario actual, el nombre del host y la fecha en el informe:

=== "Después"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Antes"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Si ejecuta esto, notará un error -- Nextflow intenta interpretar `${USER}` como una variable de Nextflow que no existe.

??? failure "Salida del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Necesitamos escaparla para que Bash pueda manejarla en su lugar.

Corrija esto escapando las variables del shell y las sustituciones de comandos con una barra invertida (`\`):

=== "Después"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Antes"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

¡Ahora funciona! La barra invertida (`\`) le dice a Nextflow "no interpretes esto, pásalo a Bash".

### Conclusión

En esta sección, ha aprendido técnicas de **procesamiento de strings**:

- **Expresiones regulares para análisis de archivos**: Usar el operador `=~` y patrones regex (`~/patrón/`) para extraer metadatos de convenciones complejas de nomenclatura de archivos
- **Generación dinámica de scripts**: Usar lógica condicional (if/else, operadores ternarios) para generar diferentes strings de script basados en las características de la entrada
- **Interpolación de variables**: Entender cuándo Nextflow interpreta los strings versus cuándo lo hace el shell
  - `${var}` - Variables de Nextflow (interpoladas por Nextflow en tiempo de compilación del workflow)
  - `\${var}` - Variables de entorno del shell (escapadas, pasadas a bash en tiempo de ejecución)
  - `\$(cmd)` - Sustitución de comandos del shell (escapada, ejecutada por bash en tiempo de ejecución)

Estos patrones de procesamiento y generación de strings son esenciales para manejar los diversos formatos de archivo y convenciones de nomenclatura que encontrará en workflows reales de bioinformática.

---

## 3. Creando Funciones Reutilizables

La lógica compleja del workflow en línea dentro de los operadores de canal o las definiciones de procesos reduce la legibilidad y el mantenimiento. Las **funciones** le permiten extraer esta lógica en componentes con nombre y reutilizables.

Nuestra operación map se ha vuelto larga y compleja. Extraigámosla en una función reutilizable usando la palabra clave `def`.

Para ilustrar cómo se ve eso con nuestro workflow existente, realice la modificación a continuación, usando `def` para definir una función reutilizable llamada `separateMetadata`:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

Al extraer esta lógica en una función, hemos reducido la lógica real del workflow a algo mucho más limpio:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Esto hace que la lógica del workflow sea mucho más fácil de leer y entender de un vistazo. La función `separateMetadata` encapsula toda la lógica compleja para analizar y enriquecer metadatos, haciéndola reutilizable y comprobable.

Ejecute el workflow para asegurarse de que sigue funcionando:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

La salida debería mostrar ambos procesos completándose con éxito. El workflow ahora es mucho más limpio y fácil de mantener, con toda la lógica compleja de procesamiento de metadatos encapsulada en la función `separateMetadata`.

### Conclusión

En esta sección, ha aprendido sobre **creación de funciones**:

- **Definir funciones con `def`**: La palabra clave para crear funciones con nombre (como `def` en Python o `function` en JavaScript)
- **Alcance de las funciones**: Las funciones definidas a nivel de script son accesibles en todo su workflow de Nextflow
- **Valores de retorno**: Las funciones devuelven automáticamente la última expresión, o use `return` explícito
- **Código más limpio**: Extraer lógica compleja en funciones es una práctica fundamental de ingeniería de software en cualquier lenguaje

A continuación, exploraremos cómo usar closures en las directivas de los procesos para la asignación dinámica de recursos.

---

## 4. Directivas de Recursos Dinámicas con Closures

Hasta ahora hemos usado scripting en el bloque `script` de los procesos. Pero los **closures** (introducidos en la Sección 1.1) también son increíblemente útiles en las directivas de los procesos, especialmente para la asignación dinámica de recursos. Agreguemos directivas de recursos a nuestro proceso FASTP que se adapten según las características de la muestra.

### 4.1. Asignación de recursos específica por muestra

Actualmente, nuestro proceso FASTP usa recursos predeterminados. Hagámoslo más inteligente asignando más CPUs para muestras de alta profundidad. Edite `modules/fastp.nf` para incluir una directiva `cpus` dinámica y una directiva `memory` estática:

=== "Después"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Antes"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

El closure `{ meta.depth > 40000000 ? 2 : 1 }` usa el **operador ternario** (cubierto en la Sección 1.1) y se evalúa para cada tarea, permitiendo la asignación de recursos por muestra. Las muestras de alta profundidad (>40M lecturas) obtienen 2 CPUs, mientras que las demás obtienen 1 CPU.

!!! note "Nota"

    **Acceso a Variables de Entrada en Directivas**

    El closure puede acceder a cualquier variable de entrada (como `meta` aquí) porque Nextflow evalúa estos closures en el contexto de la ejecución de cada tarea.

Ejecute el workflow nuevamente con la opción `-ansi-log false` para facilitar la visualización de los hashes de las tareas.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Puede verificar el comando exacto de `docker` que se ejecutó para ver la asignación de CPU para cualquier tarea dada:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Debería ver algo como:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

En este ejemplo hemos elegido una muestra que solicitó 2 CPUs (`--cpu-shares 2048`), porque era una muestra de alta profundidad, pero debería ver diferentes asignaciones de CPU dependiendo de la profundidad de la muestra. Pruebe esto también para las otras tareas.

### 4.2. Estrategias de reintento

Otro patrón poderoso es usar `task.attempt` para estrategias de reintento. Para mostrar por qué esto es útil, comenzaremos reduciendo la asignación de memoria de FASTP a menos de lo que necesita. Cambie la directiva `memory` en `modules/fastp.nf` a `1.GB`:

=== "Después"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Antes"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... y ejecute el workflow nuevamente:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Esto indica que el proceso fue terminado por exceder los límites de memoria.

Este es un escenario muy común en workflows del mundo real -- a veces simplemente no se sabe cuánta memoria necesitará una tarea hasta que se ejecuta.

Para hacer nuestro workflow más robusto, podemos implementar una estrategia de reintento que aumente la asignación de memoria en cada intento, usando nuevamente un closure de Groovy. Modifique la directiva `memory` para multiplicar la memoria base por `task.attempt`, y agregue las directivas `errorStrategy 'retry'` y `maxRetries 2`:

=== "Después"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Antes"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Ahora si el proceso falla por memoria insuficiente, Nextflow reintentará con más memoria:

- Primer intento: 1 GB (task.attempt = 1)
- Segundo intento: 2 GB (task.attempt = 2)

... y así sucesivamente, hasta el límite de `maxRetries`.

### Conclusión

Las directivas dinámicas con closures le permiten:

- Asignar recursos basados en las características de la entrada
- Implementar estrategias de reintento automático con recursos crecientes
- Combinar múltiples factores (metadatos, número de intento, prioridades)
- Usar lógica condicional para cálculos de recursos complejos

Esto hace que sus workflows sean tanto más eficientes (sin sobreasignación) como más robustos (reintento automático con más recursos).

---

## 5. Lógica Condicional y Control de Procesos

Anteriormente, usamos `.map()` con scripting para transformar datos del canal. Ahora usaremos lógica condicional para controlar qué procesos se ejecutan basándose en los datos -- esencial para workflows flexibles que se adaptan a diferentes tipos de muestras.

Los [operadores de dataflow](https://www.nextflow.io/docs/latest/reference/operator.html) de Nextflow toman closures evaluados en tiempo de ejecución, lo que permite que la lógica condicional impulse las decisiones del workflow basándose en el contenido del canal.

### 5.1. Enrutamiento con `.branch()`

Por ejemplo, supongamos que nuestras muestras de secuenciación necesitan ser recortadas con FASTP solo si son muestras humanas con una cobertura por encima de cierto umbral. Las muestras de ratón o las muestras de baja cobertura deberían ejecutarse con Trimgalore en su lugar (este es un ejemplo artificial, pero ilustra el punto).

Hemos proporcionado un proceso simple de Trimgalore en `modules/trimgalore.nf`, échele un vistazo si lo desea, pero los detalles no son importantes para este ejercicio. El punto clave es que queremos enrutar muestras basándonos en sus metadatos.

Incluya el nuevo módulo de `modules/trimgalore.nf`:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... y luego modifique su workflow `main.nf` para ramificar las muestras según sus metadatos y enrutarlas a través del proceso de recorte apropiado, así:

=== "Después"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Ejecute este workflow modificado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Aquí hemos usado expresiones condicionales pequeñas pero poderosas dentro del operador `.branch{}` para enrutar muestras según sus metadatos. Las muestras humanas con alta cobertura pasan por `FASTP`, mientras que todas las demás muestras pasan por `TRIMGALORE`.

### 5.2. Usando `.filter()` con Veracidad

Otro patrón poderoso para controlar la ejecución del workflow es el operador `.filter()`, que usa un closure para determinar qué elementos deben continuar en el pipeline. Dentro del closure de filter, escribirá **expresiones booleanas** que deciden qué elementos pasan.

Nextflow (como muchos lenguajes dinámicos) tiene un concepto de **"veracidad"** que determina qué valores se evalúan como `true` o `false` en contextos booleanos:

- **Verdadero**: Valores no nulos, strings no vacíos, números distintos de cero, colecciones no vacías
- **Falso**: `null`, strings vacíos `""`, cero `0`, colecciones vacías `[]` o `[:]`, `false`

Esto significa que `meta.id` solo (sin `!= null` explícito) verifica si el ID existe y no está vacío. Usemos esto para filtrar muestras que no cumplen con nuestros requisitos de calidad.

Agregue lo siguiente antes de la operación branch:

=== "Después"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtrar muestras inválidas o de baja calidad
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Ejecute el workflow nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Debido a que hemos elegido un filtro que excluye algunas muestras, se ejecutaron menos tareas.

La expresión de filtro `meta.id && meta.organism && meta.depth >= 25000000` combina veracidad con comparaciones explícitas:

- `meta.id && meta.organism` verifica que ambos campos existan y no estén vacíos (usando veracidad)
- `meta.depth >= 25000000` garantiza una profundidad de secuenciación suficiente con una comparación explícita

!!! note "Nota"

    **Veracidad en la Práctica**

    La expresión `meta.id && meta.organism` es más concisa que escribir:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Esto hace que la lógica de filtrado sea mucho más limpia y fácil de leer.

### Conclusión

En esta sección, ha aprendido a usar lógica condicional para controlar la ejecución del workflow usando las interfaces de closure de los operadores de Nextflow como `.branch{}` y `.filter{}`, aprovechando la veracidad para escribir expresiones condicionales concisas.

Nuestro pipeline ahora enruta inteligentemente las muestras a través de los procesos apropiados, pero los workflows de producción necesitan manejar datos inválidos de forma elegante. Hagamos nuestro workflow robusto contra valores faltantes o nulos.

---

## 6. Operadores de Navegación Segura y Elvis

Nuestra función `separateMetadata` actualmente asume que todos los campos del CSV están presentes y son válidos. ¿Pero qué sucede con datos incompletos? Averigüémoslo.

### 6.1. El Problema: Acceder a Propiedades que No Existen

Supongamos que queremos agregar soporte para información opcional de ejecución de secuenciación. En algunos laboratorios, las muestras pueden tener un campo adicional para el ID de ejecución de secuenciación o el número de lote, pero nuestro CSV actual no tiene esta columna. Intentemos acceder a ella de todas formas.

Modifique la función `separateMetadata` para incluir un campo run_id:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Ahora ejecute el workflow:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Esto falla con una NullPointerException.

El problema es que `row.run_id` devuelve `null` porque la columna `run_id` no existe en nuestro CSV. Cuando intentamos llamar `.toUpperCase()` en `null`, falla. Aquí es donde el operador de navegación segura salva el día.

### 6.2. Operador de Navegación Segura (`?.`)

El operador de navegación segura (`?.`) devuelve `null` en lugar de lanzar una excepción cuando se llama sobre un valor `null`. Si el objeto antes de `?.` es `null`, toda la expresión se evalúa como `null` sin ejecutar el método.

Actualice la función para usar navegación segura:

=== "Después"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Ejecute nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

¡Sin fallos! El workflow ahora maneja el campo faltante de forma elegante. Cuando `row.run_id` es `null`, el operador `?.` previene la llamada a `.toUpperCase()`, y `run_id` se convierte en `null` en lugar de causar una excepción.

### 6.3. Operador Elvis (`?:`) para Valores Predeterminados

El operador Elvis (`?:`) proporciona valores predeterminados cuando el lado izquierdo es "falso" (como se explicó anteriormente). ¡Recibe su nombre de Elvis Presley porque `?:` se parece a su famoso cabello y ojos cuando se ve de lado!

Ahora que estamos usando navegación segura, `run_id` será `null` para las muestras sin ese campo. Usemos el operador Elvis para proporcionar un valor predeterminado y agregarlo a nuestro map `sample_meta`:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

También agregue un operador `view()` en el workflow para ver los resultados:

=== "Después"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

y ejecute el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

¡Perfecto! Ahora todas las muestras tienen un campo `run` con su ID de ejecución real (en mayúsculas) o el valor predeterminado 'UNSPECIFIED'. La combinación de `?.` y `?:` proporciona tanto seguridad (sin fallos) como valores predeterminados sensatos.

Elimine el operador `.view()` ahora que hemos confirmado que funciona.

!!! tip "Consejo"

    **Combinando Navegación Segura y Elvis**

    El patrón `value?.method() ?: 'default'` es común en workflows de producción:

    - `value?.method()` - Llama al método de forma segura, devuelve `null` si `value` es `null`
    - `?: 'default'` - Proporciona un valor alternativo si el resultado es `null`

    Este patrón maneja datos faltantes/incompletos de forma elegante.

Use estos operadores de forma consistente en funciones, closures de operadores (`.map{}`, `.filter{}`), scripts de procesos y archivos de configuración. Previenen fallos al manejar datos del mundo real.

### Conclusión

- **Navegación segura (`?.`)**: Previene fallos en valores nulos -- devuelve null en lugar de lanzar una excepción
- **Operador Elvis (`?:`)**: Proporciona valores predeterminados -- `value ?: 'default'`
- **Combinación**: `value?.method() ?: 'default'` es el patrón común

Estos operadores hacen que los workflows sean resilientes a datos incompletos -- esencial para el trabajo del mundo real.

---

## 7. Validación con `error()` y `log.warn`

A veces necesita detener el workflow inmediatamente si los parámetros de entrada son inválidos. En Nextflow, puede usar funciones integradas como `error()` y `log.warn`, así como construcciones de programación estándar como sentencias `if` y lógica booleana, para implementar lógica de validación. Agreguemos validación a nuestro workflow.

Cree una función de validación antes de su bloque de workflow, llámela desde el workflow y cambie la creación del canal para usar un parámetro para la ruta del archivo CSV. Si el parámetro falta o el archivo no existe, llame a `error()` para detener la ejecución con un mensaje claro.

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Verificar que se proporcionó el parámetro de entrada
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Verificar que el archivo CSV existe
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Ahora intente ejecutar sin el archivo CSV:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

El workflow se detiene inmediatamente con un mensaje de error claro en lugar de fallar misteriosamente más tarde.

Ahora ejecútelo con un archivo inexistente:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Finalmente, ejecútelo con el archivo correcto:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

Esta vez se ejecuta con éxito.

También puede agregar validación dentro de la función `separateMetadata`. Usemos el `log.warn` no fatal para emitir advertencias para muestras con baja profundidad de secuenciación, pero aún permitir que el workflow continúe:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validar que los datos tienen sentido
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Ejecute el workflow nuevamente con el CSV original:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Vemos una advertencia sobre la baja profundidad de secuenciación para una de las muestras.

### Conclusión

- **`error()`**: Detiene el workflow inmediatamente con un mensaje claro
- **`log.warn`**: Emite advertencias sin detener el workflow
- **Validación temprana**: Verificar las entradas antes del procesamiento para fallar rápido con errores útiles
- **Funciones de validación**: Crear lógica de validación reutilizable que puede llamarse al inicio del workflow

La validación adecuada hace que los workflows sean más robustos y fáciles de usar al detectar problemas temprano con mensajes de error claros.

---

## 8. Manejadores de Eventos del Workflow

Hasta ahora, hemos estado escribiendo código en nuestros scripts de workflow y definiciones de procesos. Pero hay una característica más importante que debe conocer: los manejadores de eventos del workflow.

Los manejadores de eventos son closures que se ejecutan en puntos específicos del ciclo de vida de su workflow. Son perfectos para agregar registro, notificaciones u operaciones de limpieza. Estos manejadores deben definirse en su script de workflow junto con su definición de workflow.

### 8.1. El Manejador `onComplete`

El manejador de eventos más utilizado es `onComplete`, que se ejecuta cuando su workflow termina (ya sea que haya tenido éxito o fallado). Agreguemos uno para resumir los resultados de nuestro pipeline.

Agregue el manejador de eventos a su archivo `main.nf`, dentro de su definición de workflow:

=== "Después"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Este closure se ejecuta cuando el workflow se completa. Dentro, tiene acceso al objeto `workflow` que proporciona propiedades útiles sobre la ejecución.

¡Ejecute su workflow y verá este resumen al final!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

Hagámoslo más útil agregando lógica condicional:

=== "Después"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Ahora obtenemos un resumen aún más informativo, incluyendo un mensaje de éxito/fallo y el directorio de salida si se especificó:

<!-- TODO: add run command -->

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

También puede escribir el resumen en un archivo usando operaciones de archivo:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... código de su workflow ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // Escribir en un archivo de registro
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. El Manejador `onError`

Además de `onComplete`, hay otro manejador de eventos que puede usar: `onError`, que se ejecuta solo si el workflow falla:

```groovy title="main.nf - onError handler"
workflow {
    // ... código de su workflow ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Escribir registro de error detallado
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

Puede usar múltiples manejadores juntos en su script de workflow:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... código de su workflow ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### Conclusión

En esta sección, ha aprendido:

- **Closures de manejadores de eventos**: Closures en su script de workflow que se ejecutan en diferentes puntos del ciclo de vida
- **Manejador `onComplete`**: Para resúmenes de ejecución e informes de resultados
- **Manejador `onError`**: Para manejo de errores y registro de fallos
- **Propiedades del objeto workflow**: Acceder a `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Los manejadores de eventos muestran cómo puede usar todo el poder del lenguaje Nextflow dentro de sus scripts de workflow para agregar capacidades sofisticadas de registro y notificación.

---

## Resumen

¡Felicidades, lo logró!

A lo largo de esta misión secundaria, ha construido un pipeline completo de procesamiento de muestras que evolucionó desde el manejo básico de metadatos hasta un workflow sofisticado listo para producción.
Cada sección se construyó sobre la anterior, demostrando cómo las construcciones de programación transforman workflows simples en poderosos sistemas de procesamiento de datos, con los siguientes beneficios:

- **Código más claro**: Comprender dataflow vs scripting le ayuda a escribir workflows más organizados
- **Manejo robusto**: Los operadores de navegación segura y Elvis hacen que los workflows sean resilientes a datos faltantes
- **Procesamiento flexible**: La lógica condicional permite que sus workflows procesen diferentes tipos de muestras apropiadamente
- **Recursos adaptativos**: Las directivas dinámicas optimizan el uso de recursos según las características de la entrada

Esta progresión refleja la evolución real de los pipelines de bioinformática, desde prototipos de investigación que manejan unas pocas muestras hasta sistemas de producción que procesan miles de muestras en laboratorios e instituciones.
Cada desafío que resolvió y patrón que aprendió refleja problemas reales que enfrentan los desarrolladores al escalar workflows de Nextflow.

Aplicar estos patrones en su propio trabajo le permitirá construir workflows robustos y listos para producción.

### Patrones clave

1.  **Dataflow vs Scripting:** Aprendió a distinguir entre operaciones de dataflow (orquestación de canales) y scripting (código que manipula datos), incluyendo las diferencias cruciales entre operaciones en diferentes tipos como `collect` en Canal vs List.

    - Dataflow: orquestación de canales

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: procesamiento de datos en colecciones

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Procesamiento Avanzado de Strings**: Dominó las expresiones regulares para analizar nombres de archivos, la generación dinámica de scripts en procesos y la interpolación de variables (Nextflow vs Bash vs Shell).

    - Coincidencia de patrones

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Función con retorno condicional

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Colección de archivos a argumentos de comando (en bloque script del proceso)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Creación de Funciones Reutilizables**: Aprendió a extraer lógica compleja en funciones con nombre que pueden llamarse desde operadores de canal, haciendo los workflows más legibles y mantenibles.

    - Definir una función con nombre

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Llamar a la función con nombre en un workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Directivas de Recursos Dinámicas con Closures**: Exploró el uso de closures en las directivas de procesos para la asignación adaptativa de recursos basada en las características de la entrada.

    - Closures con nombre y composición

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures con acceso al alcance

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Lógica Condicional y Control de Procesos**: Agregó enrutamiento inteligente usando los operadores `.branch()` y `.filter()`, aprovechando la veracidad para escribir expresiones condicionales concisas.

    - Usar `.branch()` para enrutar datos a través de diferentes ramas del workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Evaluación booleana con Groovy Truth

    ```groovy
    if (sample.files) println "Has files"
    ```

    - Usar `filter()` para crear subconjuntos de datos con 'veracidad'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operadores de Navegación Segura y Elvis**: Hizo el pipeline robusto contra datos faltantes usando `?.` para acceso null-safe a propiedades y `?:` para proporcionar valores predeterminados.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validación con error() y log.warn**: Aprendió a validar entradas temprano y fallar rápido con mensajes de error claros.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Manejadores de Eventos de Configuración**: Aprendió a usar manejadores de eventos del workflow (`onComplete` y `onError`) para registro, notificaciones y gestión del ciclo de vida.

    - Usar `onComplete` para registrar y notificar

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Usar `onError` para tomar acción específicamente en caso de fallo

    ```groovy
    workflow.onError = {
        // Escribir registro de error detallado
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Recursos adicionales

- [Referencia del Lenguaje Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operadores de Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Sintaxis de Script de Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteca Estándar de Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Asegúrese de consultar estos recursos cuando necesite explorar características más avanzadas.

Se beneficiará de practicar y ampliar sus habilidades para:

- Escribir workflows más limpios con una separación adecuada entre dataflow y scripting
- Dominar la interpolación de variables para evitar errores comunes con variables de Nextflow, Bash y shell
- Usar directivas de recursos dinámicas para workflows eficientes y adaptativos
- Transformar colecciones de archivos en argumentos de línea de comando correctamente formateados
- Manejar diferentes convenciones de nomenclatura de archivos y formatos de entrada de forma elegante usando regex y procesamiento de strings
- Construir código reutilizable y mantenible usando patrones avanzados de closure y programación funcional
- Procesar y organizar conjuntos de datos complejos usando operaciones de colección
- Agregar validación, manejo de errores y registro para hacer sus workflows listos para producción
- Implementar gestión del ciclo de vida del workflow con manejadores de eventos

---

## ¿Qué sigue?

Regrese al [menú de Misiones Secundarias](../) o haga clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
