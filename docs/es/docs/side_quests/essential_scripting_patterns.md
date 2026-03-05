# Patrones Esenciales de Scripting en Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow es un lenguaje de programación que se ejecuta en la Máquina Virtual de Java. Aunque Nextflow está construido sobre [Groovy](http://groovy-lang.org/) y comparte gran parte de su sintaxis, Nextflow es más que solo "Groovy con extensiones" -- es un lenguaje independiente con una [sintaxis](https://nextflow.io/docs/latest/reference/syntax.html) y [biblioteca estándar](https://nextflow.io/docs/latest/reference/stdlib.html) completamente especificadas.

Puede escribir mucho código en Nextflow sin aventurarse más allá de la sintaxis básica para variables, mapas y listas. La mayoría de los tutoriales de Nextflow se centran en la orquestación de workflows (canales, procesos y flujo de datos), y puede llegar sorprendentemente lejos con solo eso.

Sin embargo, cuando necesita manipular datos, analizar nombres de archivos complejos, implementar lógica condicional o construir workflows de producción robustos, ayuda pensar en dos aspectos distintos de su código: **flujo de datos** (canales, operadores, procesos y workflows) y **scripting** (el código dentro de closures, funciones y scripts de procesos). Aunque esta distinción es algo arbitraria—todo es código Nextflow—proporciona un modelo mental útil para entender cuándo está orquestando su pipeline versus cuándo está manipulando datos. Dominar ambos mejora dramáticamente su capacidad para escribir workflows claros y mantenibles.

### Objetivos de aprendizaje

Esta misión secundaria lo lleva en un viaje práctico desde conceptos básicos hasta patrones listos para producción.
Transformaremos un workflow simple de lectura de CSV en un pipeline bioinformático sofisticado, evolucionándolo paso a paso a través de desafíos realistas:

- **Entender límites:** Distinguir entre operaciones de flujo de datos y scripting, y comprender cómo trabajan juntos
- **Manipulación de datos:** Extraer, transformar y crear subconjuntos de mapas y colecciones usando operadores poderosos
- **Procesamiento de cadenas:** Analizar esquemas complejos de nomenclatura de archivos con patrones regex y dominar la interpolación de variables
- **Funciones reutilizables:** Extraer lógica compleja en funciones nombradas para workflows más limpios y mantenibles
- **Lógica dinámica:** Construir procesos que se adapten a diferentes tipos de entrada y usar closures para asignación dinámica de recursos
- **Enrutamiento condicional:** Enrutar inteligentemente muestras a través de diferentes procesos según sus características de metadata
- **Operaciones seguras:** Manejar datos faltantes con gracia usando operadores seguros ante nulos y validar entradas con mensajes de error claros
- **Manejadores basados en configuración:** Usar manejadores de eventos de workflow para registro, notificaciones y gestión del ciclo de vida

### Requisitos previos

Antes de emprender esta misión secundaria, debe:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirse cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, trabajo con archivos, metadata)
- Tener familiaridad básica con construcciones de programación comunes (variables, mapas, listas)

Este tutorial explicará conceptos de programación a medida que los encontremos, por lo que no necesita experiencia extensa en programación.
Comenzaremos con conceptos fundamentales y construiremos hasta patrones avanzados.

---

## 0. Primeros pasos

#### Abrir el codespace de capacitación

Si aún no lo ha hecho, asegúrese de abrir el entorno de capacitación como se describe en [Configuración del Entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

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

Nuestro CSV de muestra contiene información sobre muestras biológicas que necesitan diferentes procesamientos según sus características:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Usaremos este conjunto de datos realista para explorar técnicas de programación prácticas que encontrará en workflows bioinformáticos reales.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He configurado mi directorio de trabajo apropiadamente
<!-- - [ ] I understand the assignment -->

Si puede marcar todas las casillas, está listo para comenzar.

---

## 1. Flujo de Datos vs Scripting: Entendiendo los Límites

### 1.1. Identificando Qué es Qué

Al escribir workflows de Nextflow, es importante distinguir entre **flujo de datos** (cómo se mueven los datos a través de canales y procesos) y **scripting** (el código que manipula datos y toma decisiones). Construyamos un workflow que demuestre cómo trabajan juntos.

#### 1.1.1. Workflow Básico de Nextflow

Comience con un workflow simple que solo lee el archivo CSV (ya lo hemos hecho por usted en `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

El bloque `workflow` define la estructura de nuestro pipeline, mientras que `channel.fromPath()` crea un canal desde una ruta de archivo. El operador `.splitCsv()` procesa el archivo CSV y convierte cada fila en una estructura de datos de tipo mapa.

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

Ahora vamos a agregar scripting para transformar los datos, usando el operador `.map()` que probablemente ya conoce. Este operador toma un 'closure' donde podemos escribir código para transformar cada elemento.

!!! note "Nota"

    Un **closure** es un bloque de código que puede pasarse y ejecutarse más tarde. Piense en él como una función que define en línea. Los closures se escriben con llaves `{ }` y pueden tomar parámetros. ¡Son fundamentales para cómo funcionan los operadores de Nextflow y si ha estado escribiendo Nextflow por un tiempo, puede que ya los haya estado usando sin darse cuenta!

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

Este es nuestro primer **closure** - una función anónima que puede pasar como argumento (similar a lambdas en Python o funciones flecha en JavaScript). Los closures son esenciales para trabajar con operadores de Nextflow.

El closure `{ row -> return row }` toma un parámetro `row` (podría ser cualquier nombre: `item`, `sample`, etc.).

Cuando el operador `.map()` procesa cada elemento del canal, pasa ese elemento a su closure. Aquí, `row` contiene una fila CSV a la vez.

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

El mapa `sample_meta` es una estructura de datos clave-valor (como diccionarios en Python, objetos en JavaScript, o hashes en Ruby) que almacena información relacionada: ID de muestra, organismo, tipo de tejido, profundidad de secuenciación y puntuación de calidad.

Usamos métodos de manipulación de cadenas como `.toLowerCase()` y `.replaceAll()` para limpiar nuestros datos, y métodos de conversión de tipos como `.toInteger()` y `.toDouble()` para convertir datos de cadena del CSV en los tipos numéricos apropiados.

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

Ahora agreguemos más scripting - esta vez usando un operador ternario para tomar decisiones basadas en valores de datos.

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

El operador ternario es una forma abreviada de una declaración if/else que sigue el patrón `condición ? valor_si_verdadero : valor_si_falso`. Esta línea significa: "Si la calidad es mayor que 40, use 'high', de lo contrario use 'normal'". Su primo, el **operador Elvis** (`?:`), proporciona valores predeterminados cuando algo es nulo o vacío - exploraremos ese patrón más adelante en este tutorial.

El operador de adición de mapas `+` crea un **nuevo mapa** en lugar de modificar el existente. Esta línea crea un nuevo mapa que contiene todos los pares clave-valor de `sample_meta` más la nueva clave `priority`.

!!! Note "Nota"

    Nunca modifique mapas pasados a closures - siempre cree nuevos usando `+` (por ejemplo). En Nextflow, los mismos datos a menudo fluyen a través de múltiples operaciones simultáneamente. Modificar un mapa in situ puede causar efectos secundarios impredecibles cuando otras operaciones referencian ese mismo objeto. Crear nuevos mapas asegura que cada operación tenga su propia copia limpia.

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

Hemos agregado exitosamente lógica condicional para enriquecer nuestros metadata con un nivel de prioridad basado en puntuaciones de calidad.

#### 1.1.5. Creando Subconjuntos de Mapas con `.subMap()`

Mientras que el operador `+` agrega claves a un mapa, a veces necesita hacer lo contrario - extraer solo claves específicas. El método `.subMap()` es perfecto para esto.

Agreguemos una línea para crear una versión simplificada de nuestros metadata que solo contenga campos de identificación:

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

Esto muestra tanto los metadata completos mostrados por la operación `view()` como el subconjunto extraído que imprimimos con `println`.

El método `.subMap()` toma una lista de claves y devuelve un nuevo mapa que contiene solo esas claves. Si una clave no existe en el mapa original, simplemente no se incluye en el resultado.

Esto es particularmente útil cuando necesita crear diferentes versiones de metadata para diferentes procesos - algunos pueden necesitar metadata completos mientras que otros necesitan solo campos de identificación mínimos.

Ahora elimine esas declaraciones println para restaurar su workflow a su estado anterior, ya que no las necesitamos en adelante.

!!! tip "Resumen de Operaciones con Mapas"

    - **Agregar claves**: `map1 + [new_key: value]` - Crea nuevo mapa con claves adicionales
    - **Extraer claves**: `map1.subMap(['key1', 'key2'])` - Crea nuevo mapa con solo las claves especificadas
    - **Ambas operaciones crean nuevos mapas** - Los mapas originales permanecen sin cambios

#### 1.1.6. Combinando Mapas y Devolviendo Resultados

Hasta ahora, solo hemos estado devolviendo lo que la comunidad de Nextflow llama el 'mapa meta', y hemos estado ignorando los archivos a los que se relacionan esos metadata. Pero si está escribiendo workflows de Nextflow, probablemente quiera hacer algo con esos archivos.

Generemos una estructura de canal que comprenda una tupla de 2 elementos: el mapa de metadata enriquecido y la ruta de archivo correspondiente. Este es un patrón común en Nextflow para pasar datos a procesos.

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

Esta estructura de tupla `[meta, file]` es un patrón común en Nextflow para pasar tanto metadata como archivos asociados a procesos.

!!! note "Nota"

    **Mapas y Metadata**: Los mapas son fundamentales para trabajar con metadata en Nextflow. Para una explicación más detallada sobre cómo trabajar con mapas de metadata, consulte la misión secundaria [Trabajando con metadata](./metadata.md).

Nuestro workflow demuestra el patrón central: **operaciones de flujo de datos** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orquestan cómo se mueven los datos a través del pipeline, mientras que **scripting** (mapas `[key: value]`, métodos de cadenas, conversiones de tipos, operadores ternarios) dentro del closure `.map()` maneja la transformación de elementos de datos individuales.

### 1.2. Entendiendo Diferentes Tipos: Channel vs List

Hasta ahora, bien, podemos distinguir entre operaciones de flujo de datos y scripting. ¿Pero qué pasa cuando el mismo nombre de método existe en ambos contextos?

Un ejemplo perfecto es el método `collect`, que existe tanto para tipos de canal como para tipos List en la biblioteca estándar de Nextflow. El método `collect()` en una List transforma cada elemento, mientras que el operador `collect()` en un canal reúne todas las emisiones del canal en un canal de un solo elemento.

Demostremos esto con algunos datos de muestra, comenzando por refrescar lo que hace el operador `collect()` del canal. Revise `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - agrupa múltiples emisiones de canal en una
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Pasos:

- Definir una List de IDs de muestra
- Crear un canal con `fromList()` que emite cada ID de muestra por separado
- Imprimir cada elemento con `view()` a medida que fluye
- Reunir todos los elementos en una sola lista con el operador `collect()` del canal
- Imprimir el resultado recopilado (elemento único que contiene todos los IDs de muestra) con un segundo `view()`

Hemos cambiado la estructura del canal, pero no hemos cambiado los datos en sí.

Ejecute el workflow para confirmar esto:

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

`view()` devuelve una salida por cada emisión del canal, por lo que sabemos que esta salida única contiene los 3 elementos originales agrupados en una lista.

Ahora veamos el método `collect` en una List en acción. Modifique `collect.nf` para aplicar el método `collect` de List a la lista original de IDs de muestra:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

En este nuevo fragmento:

- Definimos una nueva variable `formatted_ids` que usa el método `collect` de List para transformar cada ID de muestra en la lista original
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

Esta vez, NO hemos cambiado la estructura de los datos, todavía tenemos 3 elementos en la lista, pero SÍ hemos transformado cada elemento usando el método `collect` de List para producir una nueva lista con valores modificados. Esto es similar a usar el operador `map` en un canal, pero está operando en una estructura de datos List en lugar de un canal.

`collect` es un caso extremo que estamos usando aquí para hacer un punto. La lección clave es que cuando está escribiendo workflows, siempre distinga entre **estructuras de datos** (Lists, Maps, etc.) y **canales** (construcciones de flujo de datos). Las operaciones pueden compartir nombres pero comportarse completamente diferente dependiendo del tipo sobre el que se llamen.

### 1.3. El Operador Spread (`*.`) - Abreviatura para Extracción de Propiedades

Relacionado con el método `collect` de List está el operador spread (`*.`), que proporciona una forma concisa de extraer propiedades de colecciones. Es esencialmente azúcar sintáctica para un patrón común de `collect`.

Agreguemos una demostración a nuestro archivo `collect.nf`:

=== "Después"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva estructura
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

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva estructura
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

El operador spread `*.` es una abreviatura para un patrón común de collect:

```groovy
// Estos son equivalentes:
def ids = samples*.id
def ids = samples.collect { it.id }

// También funciona con llamadas a métodos:
def names = files*.getName()
def names = files.collect { it.getName() }
```

El operador spread es particularmente útil cuando necesita extraer una sola propiedad de una lista de objetos - es más legible que escribir el closure `collect` completo.

!!! tip "Cuándo Usar Spread vs Collect"

    - **Use spread (`*.`)** para acceso simple a propiedades: `samples*.id`, `files*.name`
    - **Use collect** para transformaciones o lógica compleja: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Conclusión

En esta sección, ha aprendido:

- **Flujo de datos vs scripting**: Los operadores de canal orquestan cómo fluyen los datos a través de su pipeline, mientras que el scripting transforma elementos de datos individuales
- **Entender tipos**: El mismo nombre de método (como `collect`) puede comportarse diferente dependiendo del tipo sobre el que se llame (Channel vs List)
- **El contexto importa**: Siempre sea consciente de si está trabajando con canales (flujo de datos) o estructuras de datos (scripting)

Entender estos límites es esencial para depuración, documentación y escribir workflows mantenibles.

A continuación profundizaremos en las capacidades de procesamiento de cadenas, que son esenciales para manejar datos del mundo real.

---

## 2. Procesamiento de Cadenas y Generación Dinámica de Scripts

Dominar el procesamiento de cadenas separa workflows frágiles de pipelines robustos. Esta sección cubre el análisis de nombres de archivos complejos, generación dinámica de scripts e interpolación de variables.

### 2.1. Coincidencia de Patrones y Expresiones Regulares

Los archivos bioinformáticos a menudo tienen convenciones de nomenclatura complejas que codifican metadata. Extraigamos esto automáticamente usando coincidencia de patrones con expresiones regulares.

Vamos a volver a nuestro workflow `main.nf` y agregar algo de lógica de coincidencia de patrones para extraer información adicional de muestra de los nombres de archivos. Los archivos FASTQ en nuestro conjunto de datos siguen convenciones de nomenclatura estilo Illumina con nombres como `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Estos pueden parecer crípticos, pero en realidad codifican metadata útiles como ID de muestra, número de carril y dirección de lectura. Vamos a usar capacidades regex para analizar estos nombres.

Realice el siguiente cambio a su workflow `main.nf` existente:

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

Esto demuestra **conceptos clave de procesamiento de cadenas**:

1. **Literales de expresiones regulares** usando sintaxis `~/patrón/` - esto crea un patrón regex sin necesidad de escapar barras invertidas
2. **Coincidencia de patrones** con el operador `=~` - esto intenta hacer coincidir una cadena con un patrón regex
3. **Objetos matcher** que capturan grupos con `[0][1]`, `[0][2]`, etc. - `[0]` se refiere a la coincidencia completa, `[1]`, `[2]`, etc. se refieren a grupos capturados entre paréntesis

Desglosemos el patrón regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Patrón              | Coincide                                    | Captura                            |
| ------------------- | ------------------------------------------- | ---------------------------------- |
| `^(.+)`             | Nombre de muestra desde el inicio           | Grupo 1: nombre de muestra         |
| `_S(\d+)`           | Número de muestra `_S1`, `_S2`, etc.        | Grupo 2: número de muestra         |
| `_L(\d{3})`         | Número de carril `_L001`                    | Grupo 3: carril (3 dígitos)        |
| `_(R[12])`          | Dirección de lectura `_R1` o `_R2`          | Grupo 4: dirección de lectura      |
| `_(\d{3})`          | Número de fragmento `_001`                  | Grupo 5: fragmento (3 dígitos)     |
| `\.fastq(?:\.gz)?$` | Extensión de archivo `.fastq` o `.fastq.gz` | No capturado (?: es no capturante) |

Esto analiza convenciones de nomenclatura estilo Illumina para extraer metadata automáticamente.

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

Esto muestra los metadata enriquecidos desde los nombres de archivos.

### 2.2. Generación Dinámica de Scripts en Procesos

Los bloques script de procesos son esencialmente cadenas multilínea que se pasan al shell. Puede usar **lógica condicional** (if/else, operadores ternarios) para generar dinámicamente diferentes cadenas de script basadas en características de entrada. Esto es esencial para manejar tipos de entrada diversos—como lecturas single-end vs paired-end—sin duplicar definiciones de procesos.

Agreguemos un proceso a nuestro workflow que demuestre este patrón. Abra `modules/fastp.nf` y eche un vistazo:

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

El proceso toma archivos FASTQ como entrada y ejecuta la herramienta `fastp` para recortar adaptadores y filtrar lecturas de baja calidad. Desafortunadamente, la persona que escribió este proceso no permitió las lecturas single-end que tenemos en nuestro conjunto de datos de ejemplo. Agreguémoslo a nuestro workflow y veamos qué sucede:

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

Puede ver que el proceso está intentando ejecutar `fastp` con un valor `null` para el segundo archivo de entrada, lo que está causando que falle. Esto es porque nuestro conjunto de datos contiene lecturas single-end, pero el proceso está codificado para esperar lecturas paired-end (dos archivos de entrada a la vez).

Arregle esto agregando lógica condicional al bloque `script:` del proceso `FASTP`. Una declaración if/else verifica el conteo de archivos de lectura y ajusta el comando en consecuencia.

=== "Después"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Detección simple de single-end vs paired-end
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

Ahora el workflow puede manejar tanto lecturas single-end como paired-end con gracia. La lógica condicional verifica el número de archivos de entrada y construye el comando apropiado para `fastp`. Veamos si funciona:

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

¡Se ve bien! Si verificamos los comandos reales que se ejecutaron (personalice para su hash de tarea):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Podemos ver que Nextflow eligió correctamente el comando correcto para lecturas single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Otro uso común de lógica dinámica de script puede verse en [el módulo de Genómica de Nextflow para Ciencia](../../nf4science/genomics/02_joint_calling). En ese módulo, el proceso GATK que se llama puede tomar múltiples archivos de entrada, pero cada uno debe tener el prefijo `-V` para formar una línea de comando correcta. El proceso usa scripting para transformar una colección de archivos de entrada (`all_gvcfs`) en los argumentos de comando correctos:

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

Estos patrones de usar scripting en bloques script de procesos son extremadamente poderosos y pueden aplicarse en muchos escenarios - desde manejar tipos de entrada variables hasta construir argumentos de línea de comando complejos desde colecciones de archivos, haciendo que sus procesos sean verdaderamente adaptables a los diversos requisitos de datos del mundo real.

### 2.3. Interpolación de Variables: Variables de Nextflow y Shell

Los scripts de procesos mezclan variables de Nextflow, variables de shell y sustituciones de comandos, cada una con diferente sintaxis de interpolación. Usar la sintaxis incorrecta causa errores. Exploremos estos con un proceso que crea un reporte de procesamiento.

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

Este proceso escribe un reporte simple con el ID de muestra y el nombre de archivo. Ahora ejecutémoslo para ver qué sucede cuando necesitamos mezclar diferentes tipos de variables.

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

Ahora ejecute el workflow y verifique los reportes generados en `results/reports/`. Deberían contener información básica sobre cada muestra.

<!-- TODO: add the run command -->

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

¿Pero qué pasa si queremos agregar información sobre cuándo y dónde ocurrió el procesamiento? Modifiquemos el proceso para usar variables **shell** y un poco de sustitución de comandos para incluir el usuario actual, nombre de host y fecha en el reporte:

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

Si ejecuta esto, notará un error - Nextflow intenta interpretar `${USER}` como una variable de Nextflow que no existe.

??? failure "Salida del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Necesitamos escaparlo para que Bash pueda manejarlo en su lugar.

Arregle esto escapando las variables de shell y sustituciones de comandos con una barra invertida (`\`):

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

¡Ahora funciona! La barra invertida (`\`) le dice a Nextflow "no interpretes esto, pásalo a Bash."

### Conclusión

En esta sección, ha aprendido técnicas de **procesamiento de cadenas**:

- **Expresiones regulares para análisis de archivos**: Usando el operador `=~` y patrones regex (`~/patrón/`) para extraer metadata de convenciones complejas de nomenclatura de archivos
- **Generación dinámica de scripts**: Usando lógica condicional (if/else, operadores ternarios) para generar diferentes cadenas de script basadas en características de entrada
- **Interpolación de variables**: Entender cuándo Nextflow interpreta cadenas vs cuándo lo hace el shell
  - `${var}` - Variables de Nextflow (interpoladas por Nextflow en tiempo de compilación del workflow)
  - `\${var}` - Variables de entorno del shell (escapadas, pasadas a bash en tiempo de ejecución)
  - `\$(cmd)` - Sustitución de comandos del shell (escapada, ejecutada por bash en tiempo de ejecución)

Estos patrones de procesamiento y generación de cadenas son esenciales para manejar los diversos formatos de archivo y convenciones de nomenclatura que encontrará en workflows bioinformáticos del mundo real.

---

## 3. Creando Funciones Reutilizables

La lógica compleja de workflow en línea en operadores de canal o definiciones de procesos reduce la legibilidad y mantenibilidad. Las **funciones** le permiten extraer esta lógica en componentes nombrados y reutilizables.

Nuestra operación map ha crecido larga y compleja. Extraigámosla en una función reutilizable usando la palabra clave `def`.

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

Esto hace que la lógica del workflow sea mucho más fácil de leer y entender de un vistazo. La función `separateMetadata` encapsula toda la lógica compleja para analizar y enriquecer metadata, haciéndola reutilizable y comprobable.

Ejecute el workflow para asegurarse de que todavía funciona:

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

La salida debería mostrar ambos procesos completándose exitosamente. El workflow ahora es mucho más limpio y fácil de mantener, con toda la lógica compleja de procesamiento de metadata encapsulada en la función `separateMetadata`.

### Conclusión

En esta sección, ha aprendido **creación de funciones**:

- **Definir funciones con `def`**: La palabra clave para crear funciones nombradas (como `def` en Python o `function` en JavaScript)
- **Alcance de funciones**: Las funciones definidas a nivel de script son accesibles en todo su workflow de Nextflow
- **Valores de retorno**: Las funciones devuelven automáticamente la última expresión, o usan `return` explícito
- **Código más limpio**: Extraer lógica compleja en funciones es una práctica fundamental de ingeniería de software en cualquier lenguaje

A continuación, exploraremos cómo usar closures en directivas de procesos para asignación dinámica de recursos.

---

## 4. Directivas de Recursos Dinámicas con Closures

Hasta ahora hemos usado scripting en el bloque `script` de procesos. Pero los **closures** (introducidos en la Sección 1.1) también son increíblemente útiles en directivas de procesos, especialmente para asignación dinámica de recursos. Agreguemos directivas de recursos a nuestro proceso FASTP que se adapten según las características de la muestra.

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

El closure `{ meta.depth > 40000000 ? 2 : 1 }` usa el **operador ternario** (cubierto en la Sección 1.1) y se evalúa para cada tarea, permitiendo asignación de recursos por muestra. Las muestras de alta profundidad (>40M lecturas) obtienen 2 CPUs, mientras que otras obtienen 1 CPU.

!!! note "Accediendo Variables de Entrada en Directivas"

    El closure puede acceder a cualquier variable de entrada (como `meta` aquí) porque Nextflow evalúa estos closures en el contexto de cada ejecución de tarea.

Ejecute el workflow nuevamente con la opción `-ansi-log false` para facilitar ver los hashes de tareas.

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

En este ejemplo hemos elegido un ejemplo que solicitó 2 CPUs (`--cpu-shares 2048`), porque era una muestra de alta profundidad, pero debería ver diferentes asignaciones de CPU dependiendo de la profundidad de la muestra. Pruebe esto para las otras tareas también.

### 4.2. Estrategias de reintento

Otro patrón poderoso es usar `task.attempt` para estrategias de reintento. Para mostrar por qué esto es útil, vamos a comenzar reduciendo la asignación de memoria a FASTP a menos de lo que necesita. Cambie la directiva `memory` en `modules/fastp.nf` a `1.GB`:

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

Este es un escenario muy común en workflows del mundo real - a veces simplemente no sabe cuánta memoria necesitará una tarea hasta que la ejecute.

Para hacer nuestro workflow más robusto, podemos implementar una estrategia de reintento que aumente la asignación de memoria en cada intento, una vez más usando un closure de Groovy. Modifique la directiva `memory` para multiplicar la memoria base por `task.attempt`, y agregue las directivas `errorStrategy 'retry'` y `maxRetries 2`:

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

Ahora si el proceso falla debido a memoria insuficiente, Nextflow reintentará con más memoria:

- Primer intento: 1 GB (task.attempt = 1)
- Segundo intento: 2.GB (task.attempt = 2)

... y así sucesivamente, hasta el límite de `maxRetries`.

### Conclusión

Las directivas dinámicas con closures le permiten:

- Asignar recursos basados en características de entrada
- Implementar estrategias de reintento automáticas con recursos crecientes
- Combinar múltiples factores (metadata, número de intento, prioridades)
- Usar lógica condicional para cálculos complejos de recursos

Esto hace que sus workflows sean tanto más eficientes (no sobre-asignar) como más robustos (reintento automático con más recursos).

---

## 5. Lógica Condicional y Control de Procesos

Anteriormente, usamos `.map()` con scripting para transformar datos de canal. Ahora usaremos lógica condicional para controlar qué procesos se ejecutan según los datos—esencial para workflows flexibles que se adaptan a diferentes tipos de muestras.

Los [operadores de flujo de datos](https://www.nextflow.io/docs/latest/reference/operator.html) de Nextflow toman closures evaluados en tiempo de ejecución, habilitando lógica condicional para impulsar decisiones de workflow basadas en el contenido del canal.

### 5.1. Enrutamiento con `.branch()`

Por ejemplo, supongamos que nuestras muestras de secuenciación necesitan ser recortadas con FASTP solo si son muestras humanas con una cobertura por encima de cierto umbral. Las muestras de ratón o muestras de baja cobertura deberían ejecutarse con Trimgalore en su lugar (este es un ejemplo artificial, pero ilustra el punto).

Hemos proporcionado un proceso simple de Trimgalore en `modules/trimgalore.nf`, eche un vistazo si lo desea, pero los detalles no son importantes para este ejercicio. El punto clave es que queremos enrutar muestras según sus metadata.

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

... y luego modifique su workflow `main.nf` para ramificar muestras según sus metadata y enrutarlas a través del proceso de recorte apropiado, así:

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

Aquí, hemos usado expresiones condicionales pequeñas pero poderosas dentro del operador `.branch{}` para enrutar muestras según sus metadata. Las muestras humanas con alta cobertura pasan por `FASTP`, mientras que todas las demás muestras pasan por `TRIMGALORE`.

### 5.2. Usando `.filter()` con Truthiness

Otro patrón poderoso para controlar la ejecución del workflow es el operador `.filter()`, que usa un closure para determinar qué elementos deben continuar por el pipeline. Dentro del closure de filtro, escribirá **expresiones booleanas** que deciden qué elementos pasan.

Nextflow (como muchos lenguajes dinámicos) tiene un concepto de **"truthiness"** que determina qué valores se evalúan como `true` o `false` en contextos booleanos:

- **Truthy**: Valores no nulos, cadenas no vacías, números no cero, colecciones no vacías
- **Falsy**: `null`, cadenas vacías `""`, cero `0`, colecciones vacías `[]` o `[:]`, `false`

Esto significa que `meta.id` solo (sin `!= null` explícito) verifica si el ID existe y no está vacío. Usemos esto para filtrar muestras que no cumplan nuestros requisitos de calidad.

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

La expresión de filtro `meta.id && meta.organism && meta.depth >= 25000000` combina truthiness con comparaciones explícitas:

- `meta.id && meta.organism` verifica que ambos campos existan y no estén vacíos (usando truthiness)
- `meta.depth >= 25000000` asegura suficiente profundidad de secuenciación con una comparación explícita

!!! note "Truthiness en la Práctica"

    La expresión `meta.id && meta.organism` es más concisa que escribir:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Esto hace que la lógica de filtrado sea mucho más limpia y fácil de leer.

### Conclusión

En esta sección, ha aprendido a usar lógica condicional para controlar la ejecución del workflow usando las interfaces de closure de operadores de Nextflow como `.branch{}` y `.filter{}`, aprovechando truthiness para escribir expresiones condicionales concisas.

Nuestro pipeline ahora enruta inteligentemente muestras a través de procesos apropiados, pero los workflows de producción necesitan manejar datos inválidos con gracia. Hagamos nuestro workflow robusto contra valores faltantes o nulos.

---

## 6. Operadores de Navegación Segura y Elvis

Nuestra función `separateMetadata` actualmente asume que todos los campos CSV están presentes y son válidos. ¿Pero qué sucede con datos incompletos? Averigüémoslo.

### 6.1. El Problema: Acceder a Propiedades que No Existen

Digamos que queremos agregar soporte para información opcional de ejecución de secuenciación. En algunos laboratorios, las muestras pueden tener un campo adicional para el ID de ejecución de secuenciación o número de lote, pero nuestro CSV actual no tiene esta columna. Intentemos acceder a ella de todos modos.

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

El problema es que `row.run_id` devuelve `null` porque la columna `run_id` no existe en nuestro CSV. Cuando intentamos llamar a `.toUpperCase()` en `null`, falla. Aquí es donde el operador de navegación segura salva el día.

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

¡Sin fallas! El workflow ahora maneja el campo faltante con gracia. Cuando `row.run_id` es `null`, el operador `?.` previene la llamada a `.toUpperCase()`, y `run_id` se convierte en `null` en lugar de causar una excepción.

### 6.3. Operador Elvis (`?:`) para Valores Predeterminados

El operador Elvis (`?:`) proporciona valores predeterminados cuando el lado izquierdo es "falsy" (como se explicó anteriormente). ¡Se llama así por Elvis Presley porque `?:` parece su famoso cabello y ojos cuando se ve de lado!

Ahora que estamos usando navegación segura, `run_id` será `null` para muestras sin ese campo. Usemos el operador Elvis para proporcionar un valor predeterminado y agregarlo a nuestro mapa `sample_meta`:

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

¡Perfecto! Ahora todas las muestras tienen un campo `run` con su ID de ejecución real (en mayúsculas) o el valor predeterminado 'UNSPECIFIED'. La combinación de `?.` y `?:` proporciona tanto seguridad (sin fallas) como valores predeterminados sensatos.

Quite el operador `.view()` ahora que hemos confirmado que funciona.

!!! tip "Combinando Navegación Segura y Elvis"

    El patrón `value?.method() ?: 'default'` es común en workflows de producción:

    - `value?.method()` - Llama al método de forma segura, devuelve `null` si `value` es `null`
    - `?: 'default'` - Proporciona respaldo si el resultado es `null`

    Este patrón maneja datos faltantes/incompletos con gracia.

Use estos operadores consistentemente en funciones, closures de operadores (`.map{}`, `.filter{}`), scripts de procesos y archivos de configuración. Previenen fallas al manejar datos del mundo real.

### Conclusión

- **Navegación segura (`?.`)**: Previene fallas en valores nulos - devuelve null en lugar de lanzar excepción
- **Operador Elvis (`?:`)**: Proporciona valores predeterminados - `value ?: 'default'`
- **Combinando**: `value?.method() ?: 'default'` es el patrón común

Estos operadores hacen que los workflows sean resilientes a datos incompletos - esencial para trabajo del mundo real.

---

## 7. Validación con `error()` y `log.warn`

A veces necesita detener el workflow inmediatamente si los parámetros de entrada son inválidos. En Nextflow, puede usar funciones integradas como `error()` y `log.warn`, así como construcciones de programación estándar como declaraciones `if` y lógica booleana, para implementar lógica de validación. Agreguemos validación a nuestro workflow.

Cree una función de validación antes de su bloque workflow, llámela desde el workflow, y cambie la creación del canal para usar un parámetro para la ruta del archivo CSV. Si el parámetro falta o el archivo no existe, llame a `error()` para detener la ejecución con un mensaje claro.

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Verificar que se proporcione el parámetro de entrada
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Verificar que el archivo CSV exista
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

El workflow se detiene inmediatamente con un mensaje de error claro en lugar de fallar misteriosamente más tarde

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

Esta vez se ejecuta exitosamente.

También puede agregar validación dentro de la función `separateMetadata`. Usemos el no fatal `log.warn` para emitir advertencias para muestras con baja profundidad de secuenciación, pero aún permitir que el workflow continúe:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validar que los datos tengan sentido
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

Vemos una advertencia sobre baja profundidad de secuenciación para una de las muestras.

### Conclusión

- **`error()`**: Detiene el workflow inmediatamente con mensaje claro
- **`log.warn`**: Emite advertencias sin detener el workflow
- **Validación temprana**: Verificar entradas antes de procesar para fallar rápido con errores útiles
- **Funciones de validación**: Crear lógica de validación reutilizable que puede llamarse al inicio del workflow

La validación adecuada hace que los workflows sean más robustos y amigables para el usuario al detectar problemas temprano con mensajes de error claros.

---

## 8. Manejadores de Eventos de Workflow

Hasta ahora, hemos estado escribiendo código en nuestros scripts de workflow y definiciones de procesos. Pero hay una característica más importante que debe conocer: manejadores de eventos de workflow.

Los manejadores de eventos son closures que se ejecutan en puntos específicos del ciclo de vida de su workflow. Son perfectos para agregar registro, notificaciones u operaciones de limpieza. Estos manejadores deben definirse en su script de workflow junto con su definición de workflow.

### 8.1. El Manejador `onComplete`

El manejador de eventos más comúnmente usado es `onComplete`, que se ejecuta cuando su workflow termina (ya sea que haya tenido éxito o fallado). Agreguemos uno para resumir los resultados de nuestro pipeline.

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

¡Ejecute su workflow y verá este resumen aparecer al final!

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

Ahora obtenemos un resumen aún más informativo, incluyendo un mensaje de éxito/falla y el directorio de salida si se especifica:

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
    // ... su código de workflow ...

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
    // ... su código de workflow ...

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
    // ... su código de workflow ...

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
- **Manejador `onComplete`**: Para resúmenes de ejecución y reportes de resultados
- **Manejador `onError`**: Para manejo de errores y registro de fallas
- **Propiedades del objeto workflow**: Accediendo a `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Los manejadores de eventos muestran cómo puede usar el poder completo del lenguaje Nextflow dentro de sus scripts de workflow para agregar capacidades sofisticadas de registro y notificación.

---

## Resumen

¡Felicitaciones, lo logró!

A lo largo de esta misión secundaria, ha construido un pipeline integral de procesamiento de muestras que evolucionó desde el manejo básico de metadata hasta un workflow sofisticado y listo para producción.
Cada sección se construyó sobre la anterior, demostrando cómo las construcciones de programación transforman workflows simples en sistemas poderosos de procesamiento de datos, con los siguientes beneficios:

- **Código más claro**: Entender flujo de datos vs scripting le ayuda a escribir workflows más organizados
- **Manejo robusto**: La navegación segura y los operadores Elvis hacen que los workflows sean resilientes a datos faltantes
- **Procesamiento flexible**: La lógica condicional permite que sus workflows procesen diferentes tipos de muestras apropiadamente
- **Recursos adaptativos**: Las directivas dinámicas optimizan el uso de recursos según las características de entrada

Esta progresión refleja la evolución del mundo real de pipelines bioinformáticos, desde prototipos de investigación que manejan unas pocas muestras hasta sistemas de producción que procesan miles de muestras a través de laboratorios e instituciones.
Cada desafío que resolvió y patrón que aprendió refleja problemas reales que los desarrolladores enfrentan al escalar workflows de Nextflow.

Aplicar estos patrones en su propio trabajo le permitirá construir workflows robustos y listos para producción.

### Patrones clave

1.  **Flujo de Datos vs Scripting:** Aprendió a distinguir entre operaciones de flujo de datos (orquestación de canales) y scripting (código que manipula datos), incluyendo las diferencias cruciales entre operaciones en diferentes tipos como `collect` en Channel vs List.

    - Flujo de datos: orquestación de canales

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: procesamiento de datos en colecciones

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Procesamiento Avanzado de Cadenas**: Dominó expresiones regulares para analizar nombres de archivos, generación dinámica de scripts en procesos e interpolación de variables (Nextflow vs Bash vs Shell).

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

    - Colección de archivos a argumentos de comando (en bloque script de proceso)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Creando Funciones Reutilizables**: Aprendió a extraer lógica compleja en funciones nombradas que pueden llamarse desde operadores de canal, haciendo los workflows más legibles y mantenibles.

    - Definir una función nombrada

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

    - Llamar a la función nombrada en un workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Directivas de Recursos Dinámicas con Closures**: Exploró el uso de closures en directivas de procesos para asignación adaptativa de recursos basada en características de entrada.

    - Closures nombrados y composición

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures con acceso a alcance

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Lógica Condicional y Control de Procesos**: Agregó enrutamiento inteligente usando operadores `.branch()` y `.filter()`, aprovechando truthiness para expresiones condicionales concisas.

    - Usar `.branch()` para enrutar datos a través de diferentes ramas de workflow

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

    - Usar `filter()` para crear subconjuntos de datos con 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operadores de Navegación Segura y Elvis**: Hizo el pipeline robusto contra datos faltantes usando `?.` para acceso seguro a propiedades ante nulos y `?:` para proporcionar valores predeterminados.

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

8.  **Manejadores de Eventos de Configuración**: Aprendió a usar manejadores de eventos de workflow (`onComplete` y `onError`) para registro, notificaciones y gestión del ciclo de vida.

    - Usando `onComplete` para registrar y notificar

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

    - Usando `onError` para tomar acción específicamente en caso de falla

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
- [Sintaxis de Scripts de Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteca Estándar de Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Asegúrese de consultar estos recursos cuando necesite explorar características más avanzadas.

Se beneficiará de practicar y expandir sus habilidades para:

- Escribir workflows más limpios con separación adecuada entre flujo de datos y scripting
- Dominar la interpolación de variables para evitar errores comunes con variables de Nextflow, Bash y shell
- Usar directivas de recursos dinámicas para workflows eficientes y adaptativos
- Transformar colecciones de archivos en argumentos de línea de comando correctamente formateados
- Manejar diferentes convenciones de nomenclatura de archivos y formatos de entrada con gracia usando regex y procesamiento de cadenas
- Construir código reutilizable y mantenible usando patrones avanzados de closures y programación funcional
- Procesar y organizar conjuntos de datos complejos usando operaciones de colecciones
- Agregar validación, manejo de errores y registro para hacer sus workflows listos para producción
- Implementar gestión del ciclo de vida del workflow con manejadores de eventos

---

## ¿Qué sigue?

Regrese al [menú de Misiones Secundarias](./index.md) o haga clic en el botón en la parte inferior derecha de la página para pasar al siguiente tema en la lista.
