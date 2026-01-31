# Patrones Esenciales de Scripting en Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow es un lenguaje de programación que se ejecuta en la Máquina Virtual de Java. Aunque Nextflow está construido sobre [Groovy](http://groovy-lang.org/) y comparte gran parte de su sintaxis, Nextflow es más que solo "Groovy con extensiones" -- es un lenguaje independiente con una [sintaxis](https://nextflow.io/docs/latest/reference/syntax.html) y [biblioteca estándar](https://nextflow.io/docs/latest/reference/stdlib.html) completamente especificadas.

Se puede escribir mucho código en Nextflow sin aventurarse más allá de la sintaxis básica para variables, mapas y listas. La mayoría de los tutoriales de Nextflow se centran en la orquestación del flujo de trabajo (canales, procesos y flujo de datos), y se puede llegar sorprendentemente lejos con solo eso.

Sin embargo, cuando se necesita manipular datos, analizar nombres de archivos complejos, implementar lógica condicional o construir flujos de trabajo de producción robustos, ayuda pensar en dos aspectos distintos de tu código: **flujo de datos** (canales, operadores, procesos y workflows) y **scripting** (el código dentro de closures, funciones y scripts de procesos). Aunque esta distinción es algo arbitraria—todo es código Nextflow—proporciona un modelo mental útil para entender cuándo estás orquestando tu pipeline versus cuándo estás manipulando datos. Dominar ambos aspectos mejora drásticamente tu capacidad de escribir flujos de trabajo claros y mantenibles.

### Objetivos de aprendizaje

Esta misión secundaria te lleva en un viaje práctico desde conceptos básicos hasta patrones listos para producción.
Transformaremos un flujo de trabajo simple de lectura de CSV en un pipeline de bioinformática sofisticado, evolucionándolo paso a paso a través de desafíos realistas:

- **Comprender límites:** Distinguir entre operaciones de flujo de datos y scripting, y entender cómo trabajan juntos
- **Manipulación de datos:** Extraer, transformar y hacer subconjuntos de mapas y colecciones usando operadores poderosos
- **Procesamiento de cadenas:** Analizar esquemas complejos de nombres de archivos con patrones regex y dominar la interpolación de variables
- **Funciones reutilizables:** Extraer lógica compleja en funciones nombradas para flujos de trabajo más limpios y mantenibles
- **Lógica dinámica:** Construir procesos que se adapten a diferentes tipos de entrada y usar closures para asignación dinámica de recursos
- **Enrutamiento condicional:** Enrutar muestras inteligentemente a través de diferentes procesos según sus características de metadatos
- **Operaciones seguras:** Manejar datos faltantes con gracia usando operadores seguros ante null y validar entradas con mensajes de error claros
- **Manejadores basados en configuración:** Usar manejadores de eventos de flujo de trabajo para registro, notificaciones y gestión del ciclo de vida

### Requisitos previos

Antes de emprender esta misión secundaria, debes:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirte cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores, trabajo con archivos, metadatos)
- Tener familiaridad básica con construcciones comunes de programación (variables, mapas, listas)

Este tutorial explicará conceptos de programación a medida que los encontremos, así que no necesitas experiencia extensa en programación.
Comenzaremos con conceptos fundamentales y construiremos hasta patrones avanzados.

---

## 0. Comenzar

#### Abrir el codespace de entrenamiento

Si aún no lo has hecho, asegúrate de abrir el entorno de entrenamiento como se describe en la [Configuración del Entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde están ubicados los archivos para este tutorial.

```bash
cd side-quests/essential_scripting_patterns
```

#### Revisar los materiales

Encontrarás un archivo de flujo de trabajo principal y un directorio `data` que contiene archivos de datos de ejemplo.

```console title="Contenido del directorio"
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

Usaremos este conjunto de datos realista para explorar técnicas prácticas de programación que encontrarás en flujos de trabajo reales de bioinformática.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista de verificación de preparación

¿Crees que estás listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está funcionando
- [ ] He establecido mi directorio de trabajo apropiadamente
<!-- - [ ] I understand the assignment -->

Si puedes marcar todas las casillas, estás listo para continuar.

---

## 1. Flujo de datos vs Scripting: Comprender los límites

### 1.1. Identificar qué es qué

Al escribir flujos de trabajo en Nextflow, es importante distinguir entre **flujo de datos** (cómo los datos se mueven a través de canales y procesos) y **scripting** (el código que manipula datos y toma decisiones). Construyamos un flujo de trabajo que demuestre cómo trabajan juntos.

#### 1.1.1. Flujo de trabajo básico en Nextflow

Comienza con un flujo de trabajo simple que solo lee el archivo CSV (ya lo hemos hecho por ti en `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

El bloque `workflow` define la estructura de nuestro pipeline, mientras que `channel.fromPath()` crea un canal desde una ruta de archivo. El operador `.splitCsv()` procesa el archivo CSV y convierte cada fila en una estructura de datos de tipo mapa.

Ejecuta este flujo de trabajo para ver los datos CSV sin procesar:

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

#### 1.1.2. Agregar el operador map

Ahora vamos a agregar scripting para transformar los datos, usando el operador `.map()` que probablemente ya conoces. Este operador toma un 'closure' donde podemos escribir código para transformar cada elemento.

!!! note

    Un **closure** es un bloque de código que puede pasarse y ejecutarse más tarde. Piénsalo como una función que defines en línea. Los closures se escriben con llaves `{ }` y pueden tomar parámetros. ¡Son fundamentales para cómo funcionan los operadores de Nextflow y si has estado escribiendo Nextflow por un tiempo, puede que ya los hayas estado usando sin darte cuenta!

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

Este es nuestro primer **closure** - una función anónima que puedes pasar como argumento (similar a lambdas en Python o arrow functions en JavaScript). Los closures son esenciales para trabajar con operadores de Nextflow.

El closure `{ row -> return row }` toma un parámetro `row` (podría ser cualquier nombre: `item`, `sample`, etc.).

Cuando el operador `.map()` procesa cada elemento del canal, pasa ese elemento a tu closure. Aquí, `row` contiene una fila del CSV a la vez.

Aplica este cambio y ejecuta el flujo de trabajo:

```bash
nextflow run main.nf
```

Verás la misma salida que antes, porque simplemente estamos devolviendo la entrada sin cambios. Esto confirma que el operador map está funcionando correctamente. Ahora comencemos a transformar los datos.

#### 1.1.3. Crear una estructura de datos de tipo mapa

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

Usamos métodos de manipulación de cadenas como `.toLowerCase()` y `.replaceAll()` para limpiar nuestros datos, y métodos de conversión de tipo como `.toInteger()` y `.toDouble()` para convertir datos de cadena del CSV en los tipos numéricos apropiados.

Aplica este cambio y ejecuta el flujo de trabajo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Agregar lógica condicional

Ahora agreguemos más scripting - esta vez usando un operador ternario para tomar decisiones basadas en valores de datos.

Realiza el siguiente cambio:

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

El operador ternario es una forma abreviada de una declaración if/else que sigue el patrón `condición ? valor_si_verdadero : valor_si_falso`. Esta línea significa: "Si la calidad es mayor que 40, usa 'high', de lo contrario usa 'normal'". Su primo, el **operador Elvis** (`?:`), proporciona valores predeterminados cuando algo es null o está vacío - exploraremos ese patrón más adelante en este tutorial.

El operador de adición de mapas `+` crea un **nuevo mapa** en lugar de modificar el existente. Esta línea crea un nuevo mapa que contiene todos los pares clave-valor de `sample_meta` más la nueva clave `priority`.

!!! Note

    Nunca modifiques mapas pasados a closures - siempre crea nuevos usando `+` (por ejemplo). En Nextflow, los mismos datos a menudo fluyen a través de múltiples operaciones simultáneamente. Modificar un mapa in situ puede causar efectos secundarios impredecibles cuando otras operaciones hacen referencia a ese mismo objeto. Crear nuevos mapas asegura que cada operación tenga su propia copia limpia.

Ejecuta el flujo de trabajo modificado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Hemos agregado con éxito lógica condicional para enriquecer nuestros metadatos con un nivel de prioridad basado en puntuaciones de calidad.

#### 1.1.5. Hacer subconjuntos de mapas con `.subMap()`

Mientras que el operador `+` agrega claves a un mapa, a veces necesitas hacer lo contrario - extraer solo claves específicas. El método `.subMap()` es perfecto para esto.

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
                println "Solo campos de ID: ${id_only}"

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

Ejecuta el flujo de trabajo modificado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    Solo campos de ID: [id:sample_001, organism:human, tissue:liver]
    Solo campos de ID: [id:sample_002, organism:mouse, tissue:brain]
    Solo campos de ID: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Esto muestra tanto los metadatos completos mostrados por la operación `view()` como el subconjunto extraído que imprimimos con `println`.

El método `.subMap()` toma una lista de claves y devuelve un nuevo mapa que contiene solo esas claves. Si una clave no existe en el mapa original, simplemente no se incluye en el resultado.

Esto es particularmente útil cuando necesitas crear diferentes versiones de metadatos para diferentes procesos - algunos podrían necesitar metadatos completos mientras que otros necesitan solo campos mínimos de identificación.

Ahora elimina esas declaraciones println para restaurar tu flujo de trabajo a su estado anterior, ya que no las necesitamos en el futuro.

!!! tip "Resumen de operaciones con mapas"

    - **Agregar claves**: `map1 + [new_key: value]` - Crea nuevo mapa con claves adicionales
    - **Extraer claves**: `map1.subMap(['key1', 'key2'])` - Crea nuevo mapa con solo las claves especificadas
    - **Ambas operaciones crean nuevos mapas** - Los mapas originales permanecen sin cambios

#### 1.1.6. Combinar mapas y devolver resultados

Hasta ahora, solo hemos estado devolviendo lo que la comunidad Nextflow llama el 'meta map', y hemos estado ignorando los archivos a los que esos metadatos se relacionan. Pero si estás escribiendo flujos de trabajo en Nextflow, probablemente quieras hacer algo con esos archivos.

Generemos una estructura de canal que comprenda una tupla de 2 elementos: el mapa de metadatos enriquecido y la ruta de archivo correspondiente. Este es un patrón común en Nextflow para pasar datos a procesos.

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

Aplica este cambio y ejecuta el flujo de trabajo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Esta estructura de tupla `[meta, file]` es un patrón común en Nextflow para pasar tanto metadatos como archivos asociados a procesos.

!!! note

    **Mapas y metadatos**: Los mapas son fundamentales para trabajar con metadatos en Nextflow. Para una explicación más detallada de trabajo con mapas de metadatos, consulta la misión secundaria [Trabajar con metadatos](./metadata.md).

Nuestro flujo de trabajo demuestra el patrón central: **operaciones de flujo de datos** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orquestan cómo los datos se mueven a través del pipeline, mientras que **scripting** (mapas `[key: value]`, métodos de cadenas, conversiones de tipo, operadores ternarios) dentro del closure `.map()` maneja la transformación de elementos de datos individuales.

### 1.2. Comprender diferentes tipos: canal vs lista

Hasta ahora todo bien, podemos distinguir entre operaciones de flujo de datos y scripting. Pero ¿qué pasa cuando el mismo nombre de método existe en ambos contextos?

Un ejemplo perfecto es el método `collect`, que existe tanto para tipos de canal como para tipos List en la biblioteca estándar de Nextflow. El método `collect()` en una List transforma cada elemento, mientras que el operador `collect()` en un canal reúne todas las emisiones del canal en un canal de un solo elemento.

Demostremos esto con algunos datos de muestra, comenzando por refrescar lo que hace el operador `collect()` del canal. Revisa `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - agrupa múltiples emisiones de canal en una
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Elemento de canal individual: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} elementos agrupados en 1)" }
```

Pasos:

- Define una List de IDs de muestra
- Crea un canal con `fromList()` que emite cada ID de muestra por separado
- Imprime cada elemento con `view()` a medida que fluye
- Reúne todos los elementos en una sola lista con el operador `collect()` del canal
- Imprime el resultado recopilado (un solo elemento que contiene todos los IDs de muestra) con un segundo `view()`

Hemos cambiado la estructura del canal, pero no hemos cambiado los datos en sí.

Ejecuta el flujo de trabajo para confirmar esto:

```bash
nextflow run collect.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Elemento de canal individual: sample_001
    Elemento de canal individual: sample_002
    Elemento de canal individual: sample_003
    Resultado de channel.collect(): [sample_001, sample_002, sample_003] (3 elementos agrupados en 1)
    ```

`view()` devuelve una salida por cada emisión del canal, así que sabemos que esta única salida contiene los 3 elementos originales agrupados en una lista.

Ahora veamos el método `collect` en una List en acción. Modifica `collect.nf` para aplicar el método `collect` de List a la lista original de IDs de muestra:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento de canal individual: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} elementos agrupados en 1)" }

    // List.collect() - transforma cada elemento, preserva estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultado de List.collect(): ${formatted_ids} (${sample_ids.size()} elementos transformados en ${formatted_ids.size()})"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento de canal individual: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} elementos agrupados en 1)" }
    ```

En este nuevo fragmento:

- Definimos una nueva variable `formatted_ids` que usa el método `collect` de List para transformar cada ID de muestra en la lista original
- Imprimimos el resultado usando `println`

Ejecuta el flujo de trabajo modificado:

```bash
nextflow run collect.nf
```

??? success "Salida del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    Resultado de List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elementos transformados en 3)
    Elemento de canal individual: sample_001
    Elemento de canal individual: sample_002
    Elemento de canal individual: sample_003
    Resultado de channel.collect(): [sample_001, sample_002, sample_003] (3 elementos agrupados en 1)
    ```

Esta vez, NO hemos cambiado la estructura de los datos, todavía tenemos 3 elementos en la lista, pero SÍ hemos transformado cada elemento usando el método `collect` de List para producir una nueva lista con valores modificados. Esto es similar a usar el operador `map` en un canal, pero está operando en una estructura de datos List en lugar de un canal.

`collect` es un caso extremo que estamos usando aquí para hacer un punto. La lección clave es que cuando estás escribiendo flujos de trabajo, siempre distingue entre **estructuras de datos** (Lists, Maps, etc.) y **canales** (construcciones de flujo de datos). Las operaciones pueden compartir nombres pero comportarse completamente diferente dependiendo del tipo sobre el que se llaman.

### 1.3. El operador de propagación (`*.`) - Atajo para extracción de propiedades

Relacionado con el método `collect` de List está el operador de propagación (`*.`), que proporciona una forma concisa de extraer propiedades de colecciones. Es esencialmente azúcar sintáctica para un patrón común de `collect`.

Agreguemos una demostración a nuestro archivo `collect.nf`:

=== "Después"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento de canal individual: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} elementos agrupados en 1)" }

    // List.collect() - transforma cada elemento, preserva estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultado de List.collect(): ${formatted_ids} (${sample_ids.size()} elementos transformados en ${formatted_ids.size()})"

    // Operador de propagación - acceso conciso a propiedades
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Resultado del operador de propagación: ${all_ids}"
    ```

=== "Antes"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emisiones de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Elemento de canal individual: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} elementos agrupados en 1)" }

    // List.collect() - transforma cada elemento, preserva estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultado de List.collect(): ${formatted_ids} (${sample_ids.size()} elementos transformados en ${formatted_ids.size()})"
    ```

Ejecuta el flujo de trabajo actualizado:

```bash title="Probar operador de propagación"
nextflow run collect.nf
```

??? success "Salida del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    Resultado de List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elementos transformados en 3)
    Resultado del operador de propagación: [s1, s2, s3]
    Elemento de canal individual: sample_001
    Elemento de canal individual: sample_002
    Elemento de canal individual: sample_003
    Resultado de channel.collect(): [sample_001, sample_002, sample_003] (3 elementos agrupados en 1)
    ```

El operador de propagación `*.` es una forma abreviada de un patrón común de collect:

```groovy
// Estos son equivalentes:
def ids = samples*.id
def ids = samples.collect { it.id }

// También funciona con llamadas a métodos:
def names = files*.getName()
def names = files.collect { it.getName() }
```

El operador de propagación es particularmente útil cuando necesitas extraer una sola propiedad de una lista de objetos - es más legible que escribir el closure completo de `collect`.

!!! tip "Cuándo usar propagación vs collect"

    - **Usa propagación (`*.`)** para acceso simple a propiedades: `samples*.id`, `files*.name`
    - **Usa collect** para transformaciones o lógica compleja: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Conclusión

En esta sección, has aprendido:

- **Flujo de datos vs scripting**: Los operadores de canal orquestan cómo los datos fluyen a través de tu pipeline, mientras que el scripting transforma elementos de datos individuales
- **Comprender tipos**: El mismo nombre de método (como `collect`) puede comportarse diferente dependiendo del tipo sobre el que se llama (Channel vs List)
- **El contexto importa**: Siempre sé consciente de si estás trabajando con canales (flujo de datos) o estructuras de datos (scripting)

Comprender estos límites es esencial para depuración, documentación y escribir flujos de trabajo mantenibles.

A continuación profundizaremos en capacidades de procesamiento de cadenas, que son esenciales para manejar datos del mundo real.

---

## 2. Procesamiento de cadenas y generación dinámica de scripts

Dominar el procesamiento de cadenas separa flujos de trabajo frágiles de pipelines robustos. Esta sección cubre análisis de nombres de archivos complejos, generación dinámica de scripts e interpolación de variables.

### 2.1. Coincidencia de patrones y expresiones regulares

Los archivos de bioinformática a menudo tienen convenciones de nomenclatura complejas que codifican metadatos. Extraigamos esto automáticamente usando coincidencia de patrones con expresiones regulares.

Vamos a regresar a nuestro flujo de trabajo `main.nf` y agregar algo de lógica de coincidencia de patrones para extraer información adicional de muestra de los nombres de archivos. Los archivos FASTQ en nuestro conjunto de datos siguen convenciones de nomenclatura estilo Illumina con nombres como `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Estos pueden parecer crípticos, pero en realidad codifican metadatos útiles como ID de muestra, número de carril y dirección de lectura. Vamos a usar capacidades regex para analizar estos nombres.

Realiza el siguiente cambio a tu flujo de trabajo existente `main.nf`:

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

1. **Literales de expresión regular** usando sintaxis `~/patrón/` - esto crea un patrón regex sin necesidad de escapar barras invertidas
2. **Coincidencia de patrones** con el operador `=~` - esto intenta hacer coincidir una cadena con un patrón regex
3. **Objetos matcher** que capturan grupos con `[0][1]`, `[0][2]`, etc. - `[0]` se refiere a la coincidencia completa, `[1]`, `[2]`, etc. se refieren a grupos capturados en paréntesis

Desglosemos el patrón regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Patrón              | Coincide con                                | Captura                            |
| ------------------- | ------------------------------------------- | ---------------------------------- |
| `^(.+)`             | Nombre de muestra desde el inicio           | Grupo 1: nombre de muestra         |
| `_S(\d+)`           | Número de muestra `_S1`, `_S2`, etc.        | Grupo 2: número de muestra         |
| `_L(\d{3})`         | Número de carril `_L001`                    | Grupo 3: carril (3 dígitos)        |
| `_(R[12])`          | Dirección de lectura `_R1` o `_R2`          | Grupo 4: dirección de lectura      |
| `_(\d{3})`          | Número de fragmento `_001`                  | Grupo 5: fragmento (3 dígitos)     |
| `\.fastq(?:\.gz)?$` | Extensión de archivo `.fastq` o `.fastq.gz` | No capturado (?: es no capturante) |

Esto analiza convenciones de nomenclatura estilo Illumina para extraer metadatos automáticamente.

Ejecuta el flujo de trabajo modificado:

```bash title="Probar coincidencia de patrones"
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

Esto muestra los metadatos enriquecidos desde los nombres de archivos.

### 2.2. Generación dinámica de scripts en procesos

Los bloques script de procesos son esencialmente cadenas multilínea que se pasan al shell. Puedes usar **lógica condicional** (if/else, operadores ternarios) para generar dinámicamente diferentes cadenas de script basadas en características de entrada. Esto es esencial para manejar tipos de entrada diversos—como lecturas single-end vs paired-end—sin duplicar definiciones de procesos.

Agreguemos un proceso a nuestro flujo de trabajo que demuestre este patrón. Abre `modules/fastp.nf` y echa un vistazo:

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

El proceso toma archivos FASTQ como entrada y ejecuta la herramienta `fastp` para recortar adaptadores y filtrar lecturas de baja calidad. Desafortunadamente, la persona que escribió este proceso no permitió las lecturas single-end que tenemos en nuestro conjunto de datos de ejemplo. Agreguémoslo a nuestro flujo de trabajo y veamos qué sucede:

Primero, incluye el módulo en la primera línea de tu flujo de trabajo `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Luego modifica el bloque `workflow` para conectar el canal `ch_samples` al proceso `FASTP`:

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

Ejecuta este flujo de trabajo modificado:

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

Puedes ver que el proceso está intentando ejecutar `fastp` con un valor `null` para el segundo archivo de entrada, lo que está causando que falle. Esto es porque nuestro conjunto de datos contiene lecturas single-end, pero el proceso está codificado para esperar lecturas paired-end (dos archivos de entrada a la vez).

Arregla esto agregando lógica condicional al bloque `script:` del proceso `FASTP`. Una declaración if/else verifica el conteo de archivos de lectura y ajusta el comando en consecuencia.

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

Ahora el flujo de trabajo puede manejar tanto lecturas single-end como paired-end con gracia. La lógica condicional verifica el número de archivos de entrada y construye el comando apropiado para `fastp`. Veamos si funciona:

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

¡Se ve bien! Si verificamos los comandos reales que se ejecutaron (personaliza para tu hash de tarea):

```console title="Verificar comandos ejecutados"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Podemos ver que Nextflow correctamente eligió el comando correcto para lecturas single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Otro uso común de lógica dinámica de script puede verse en [el módulo de Genómica de Nextflow para Ciencia](../../nf4science/genomics/02_joint_calling). En ese módulo, el proceso GATK que se está llamando puede tomar múltiples archivos de entrada, pero cada uno debe estar precedido por `-V` para formar una línea de comando correcta. El proceso usa scripting para transformar una colección de archivos de entrada (`all_gvcfs`) en los argumentos de comando correctos:

```groovy title="manipulación de línea de comando para GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Estos patrones de usar scripting en bloques script de procesos son extremadamente poderosos y pueden aplicarse en muchos escenarios - desde manejar tipos de entrada variables hasta construir argumentos complejos de línea de comando desde colecciones de archivos, haciendo tus procesos verdaderamente adaptables a los diversos requisitos de datos del mundo real.

### 2.3. Interpolación de variables: variables de Nextflow y Shell

Los scripts de procesos mezclan variables de Nextflow, variables de shell y sustituciones de comandos, cada una con diferente sintaxis de interpolación. Usar la sintaxis incorrecta causa errores. Exploremos estos con un proceso que crea un reporte de procesamiento.

Echa un vistazo al archivo de módulo `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Procesando ${reads}" > ${meta.id}_report.txt
    echo "Muestra: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Este proceso escribe un reporte simple con el ID de muestra y nombre de archivo. Ahora ejecutémoslo para ver qué sucede cuando necesitamos mezclar diferentes tipos de variables.

Incluye el proceso en tu `main.nf` y agrégalo al flujo de trabajo:

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

Ahora ejecuta el flujo de trabajo y verifica los reportes generados en `results/reports/`. Deberían contener información básica sobre cada muestra.

<!-- TODO: add the run command -->

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

Pero ¿qué pasa si queremos agregar información sobre cuándo y dónde ocurrió el procesamiento? Modifiquemos el proceso para usar variables **shell** y un poco de sustitución de comandos para incluir el usuario actual, nombre del host y fecha en el reporte:

=== "Después"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Procesando ${reads}" > ${meta.id}_report.txt
        echo "Muestra: ${meta.id}" >> ${meta.id}_report.txt
        echo "Procesado por: ${USER}" >> ${meta.id}_report.txt
        echo "Host: $(hostname)" >> ${meta.id}_report.txt
        echo "Fecha: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Antes"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Procesando ${reads}" > ${meta.id}_report.txt
        echo "Muestra: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Si ejecutas esto, notarás un error - Nextflow intenta interpretar `${USER}` como una variable de Nextflow que no existe.

??? failure "Salida del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Procesado por: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Necesitamos escaparlo para que Bash pueda manejarlo en su lugar.

Arregla esto escapando las variables de shell y sustituciones de comandos con una barra invertida (`\`):

=== "Después"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Procesando ${reads}" > ${meta.id}_report.txt
        echo "Muestra: ${meta.id}" >> ${meta.id}_report.txt
        echo "Procesado por: \${USER}" >> ${meta.id}_report.txt
        echo "Host: \$(hostname)" >> ${meta.id}_report.txt
        echo "Fecha: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Antes"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Procesando ${reads}" > ${meta.id}_report.txt
        echo "Muestra: ${meta.id}" >> ${meta.id}_report.txt
        echo "Procesado por: ${USER}" >> ${meta.id}_report.txt
        echo "Host: $(hostname)" >> ${meta.id}_report.txt
        echo "Fecha: $(date)" >> ${meta.id}_report.txt
        """
    ```

¡Ahora funciona! La barra invertida (`\`) le dice a Nextflow "no interpretes esto, pásalo a Bash."

### Conclusión

En esta sección, has aprendido técnicas de **procesamiento de cadenas**:

- **Expresiones regulares para análisis de archivos**: Usando el operador `=~` y patrones regex (`~/patrón/`) para extraer metadatos de convenciones de nomenclatura de archivos complejos
- **Generación dinámica de scripts**: Usando lógica condicional (if/else, operadores ternarios) para generar diferentes cadenas de script basadas en características de entrada
- **Interpolación de variables**: Entender cuándo Nextflow interpreta cadenas vs cuándo lo hace el shell
  - `${var}` - Variables de Nextflow (interpoladas por Nextflow en tiempo de compilación del flujo de trabajo)
  - `\${var}` - Variables de entorno de shell (escapadas, pasadas a bash en tiempo de ejecución)
  - `\$(cmd)` - Sustitución de comandos de shell (escapada, ejecutada por bash
