---
title: Patrones esenciales de scripting
description: Aprende técnicas avanzadas de programación en Nextflow.
weight: 1200
---

# Patrones esenciales de scripting

Esta lección te guiará a través de los patrones de programación esenciales que son fundamentales para crear flujos de trabajo efectivos de Nextflow. Cubriremos patrones para manejar la entrada de datos, transformar valores, controlar la lógica del flujo de trabajo, asignar recursos dinámicamente, y más.

## Objetivos de aprendizaje

- Entender las diferencias entre los paradigmas de flujo de datos y scripting
- Aplicar closures, operadores ternarios, y otras técnicas de Groovy
- Dominar técnicas para manipular metadatos y extraer información de archivos
- Utilizar expresiones regulares y procesamiento de cadenas para analizar nombres de archivos
- Crear funciones reutilizables para lógica compleja
- Implementar asignación dinámica de recursos y estrategias de reintento
- Agregar lógica condicional para controlar la ejecución del flujo de trabajo
- Escribir código robusto utilizando operadores de navegación segura y Elvis
- Validar entradas con mensajes de error claros
- Utilizar event handlers para gestionar la finalización del flujo de trabajo

## Prerrequisitos

- Entendimiento básico de flujos de trabajo de Nextflow
- Familiaridad con la sintaxis de DSL2
- Conocimiento básico de canales y procesos
- Entorno de desarrollo configurado con Nextflow v23.04.0 o posterior
- Acceso a Docker (o Conda) para los contenedores de software

## Cómo empezar

Este tutorial asume que tienes un conocimiento básico de la sintaxis de Nextflow y las operaciones de los canales.

Para comenzar, navega a la carpeta `side-quests/essential_scripting_patterns`:

```bash
cd side-quests/essential_scripting_patterns
```

El repositorio contiene varios archivos:

- `main.nf`: El flujo de trabajo principal
- `modules/fastp.nf`: Módulo para el procesamiento de calidad con FASTP
- `modules/generate_report.nf`: Módulo para generar reportes
- `modules/trimgalore.nf`: Módulo para un programa alternativo de recorte
- `data/samples.csv`: Un CSV con información de muestras
- `data/sequences/*.fastq`: Archivos fastq de ejemplo

Inspeccionemos los archivos:

```bash
cat main.nf
```

Deberías ver un flujo de trabajo mínimo que usa inclusiones de módulos:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
}
```

Los módulos también son simples:

```bash
cat modules/fastp.nf
```

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: fastq
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json

    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
```

```bash
cat modules/generate_report.nf
```

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {
    container 'community.wave.seqera.io/library/bash:5.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_report.txt"), emit: report

    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
}
```

Y finalmente, veamos los datos de muestra:

```bash
cat data/samples.csv
```

```csv title="data/samples.csv"
sample_id,organism,tissue_type,sequencing_depth,quality_score,file_path
SAMPLE_001,human,liver,30000000,38.5,./data/sequences/SAMPLE_001_S1_L001_R1_001.fastq
SAMPLE_002,mouse,brain,25000000,35.2,./data/sequences/SAMPLE_002_S2_L001_R1_001.fastq
SAMPLE_003,human,kidney,45000000,42.1,./data/sequences/SAMPLE_003_S3_L001_R1_001.fastq
```

Ejecutemos el flujo de trabajo base primero para asegurarnos de que funciona:

```bash
nextflow run main.nf
```

Deberías ver el flujo de trabajo ejecutarse y completarse exitosamente.

¡Ahora comencemos a mejorar el flujo de trabajo con patrones de scripting avanzados!

---

## 1. Entendiendo el Flujo de Datos vs. Scripting

Nextflow utiliza **dos paradigmas de programación distintos**:

1. **Flujo de datos**: La orquestación de canales a través de operadores (`.map`, `.filter`, `.branch`)
2. **Scripting**: Código de Groovy ejecutado dentro de closures o bloques de script en procesos

Ambos son cruciales, pero funcionan de manera diferente, por lo que esta primera sección aclarará las diferencias.

### 1.1. Closures y operadores ternarios

Un concepto fundamental en Nextflow es la **closure** - un bloque de código que se puede pasar como un objeto y ejecutar más tarde. Las closures son esenciales para los operadores de los canales (`map`, `filter`, etc.).

Tomemos nuestro operador `.map` actual y mejorémoslo para extraer más metadatos:

```groovy
.map { row ->
    tuple(
        [id: row.sample_id],
        file(row.file_path)
    )
}
```

El bloque `{ row -> ... }` es una closure que recibe un argumento `row`.

Ahora, mejorémosla para extraer más metadatos del CSV:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2-8"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="2-5"
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }
    ```

Ahora, modifiquemos el `main.nf` para utilizar estos metadatos adicionales. Primero, cambiemos la forma en que generamos el informe. Vamos a etiquetar las muestras de alta calidad como de "alta prioridad" usando un **operador ternario**.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="9-10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="9"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

El operador ternario `condition ? value_if_true : value_if_false` es una forma concisa de escribir una declaración if/else. Si `sample_meta.quality > 40`, entonces `priority` será `'high'`; de lo contrario, será `'normal'`.

Modificamos también el proceso `GENERATE_REPORT` para incluir los metadatos adicionales:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-4"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Ahora ejecutemos el flujo de trabajo actualizado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_fourier] DSL2 - revision: d3e76a7fce

    executor >  local (6)
    [a3/6b9e80] process > FASTP (sample_003)           [100%] 3 of 3 ✔
    [29/54f4b6] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Podemos verificar el resultado revisando uno de los archivos de informe:

```console
cat work/29/54f4b6b0eb90fed9e3a673b8e47629/sample_001_report.txt
```

Ahora deberías ver todos los metadatos incluidos:

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.2. La colección 'meta' vs 'reads'

En Nextflow, canales, tuplas y colecciones (como mapas y listas) son estructuras de datos fundamentales. Existe una diferencia importante entre **colecciones en el flujo de datos** (canales/operadores) y **colecciones en bloques de script**.

Para demostrar esto, modifiquemos nuestro flujo de trabajo para extraer metadatos de los nombres de archivo FASTQ utilizando expresiones regulares. Es común que los nombres de los archivos FASTQ sigan una convención como: `SAMPLE_001_S1_L001_R1_001.fastq`.

En este formato:

- `SAMPLE_001`: ID de la muestra
- `S1`: Número de muestra
- `L001`: Número de carril (lane)
- `R1`: Número de lectura (read)
- `001`: Fragmento o trozo (chunk)

Extraigamos estos metadatos del nombre del archivo FASTQ:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="12-20"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            def fastq_path = file(row.file_path)

            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]

            tuple(sample_meta + file_meta + [priority: priority], fastq_path)
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

¡Hay mucho que desempacar aquí!

1. `fastq_path.name` obtiene el nombre del archivo (sin la ruta)
2. El operador `=~` es para coincidencia de patrones regex
3. El patrón `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` captura los componentes
4. `m ? ... : ...` es un operador ternario que maneja el caso donde el nombre del archivo no coincide con el patrón
5. `m[0][2]` accede al segundo grupo de captura (los índices comienzan en 1 para los grupos)
6. La función `toInteger()` convierte la cadena capturada a un número entero
7. `sample_meta + file_meta + [priority: priority]` fusiona los tres mapas en uno solo

Ahora también modifiquemos el proceso `GENERATE_REPORT` para incluir estos nuevos metadatos:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7-10"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Ejecuta el flujo de trabajo de nuevo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_hawking] DSL2 - revision: cd0a5b0d29

    executor >  local (6)
    [b3/1cb89f] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [e6/c2f254] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Verifica uno de los archivos de informe actualizados:

```console
cat work/e6/c2f2542ec23ee6e7f0aa9c66c12e30/sample_001_report.txt
```

Ahora deberías ver los metadatos adicionales extraídos del nombre del archivo:

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Sample number: 1
Lane: 001
Read: R1
Chunk: 001
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.3. Operaciones de colección en cierres vs. manipulación de canales

Un punto confuso común en Nextflow son las operaciones de colecciones dentro de las closures. Los métodos como `collect()` funcionan de manera diferente según el contexto:

1. En los **canales de Nextflow**: `channel.collect()` es un operador que reúne todos los elementos de un canal en una sola lista
2. En **listas y mapas de Groovy**: `list.collect {...}` aplica una función a cada elemento y devuelve una nueva lista

Modifiquemos el proceso FASTP para ilustrar esta diferencia:

```groovy title="modules/fastp.nf" linenums="11" hl_lines="3-5"
script:
"""
# Demostración de operaciones de colección en script
def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
fastp ${options.join(' ')}
"""
```

Este ejemplo usa `options.join(' ')` para unir los elementos de la lista en una sola cadena con espacios entre ellos.

Sin embargo, esto dará error porque estamos tratando de ejecutar código Groovy dentro de un bloque de script bash. Cambiémoslo para mover la lógica de colección fuera del bloque de script:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="1-3"
    def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
    def cmd = "fastp ${options.join(' ')}"

    script:
    """
    $cmd
    """
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="2-3"
    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
    ```

¡Ejecuta el flujo de trabajo de nuevo y debería funcionar! Esto demuestra el uso del scripting de Groovy para la manipulación de colecciones antes de pasar el comando al bloque de script.

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [compassionate_shaw] DSL2 - revision: 3471dc57d9

    executor >  local (6)
    [32/8e94af] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [6e/e3e56d] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

### Conclusión

En esta sección, hemos aprendido sobre las diferencias entre el **flujo de datos** (operaciones de canal) y el **scripting** (código dentro de closures y bloques de script). Hemos usado:

- **Closures** para extraer y transformar metadatos de los datos de entrada
- **Operadores ternarios** (`condition ? true_value : false_value`) para lógica condicional compacta
- **Expresiones regulares** con el operador `=~` para extraer componentes de los nombres de archivos
- **Manipulación de colecciones** como `join()` para crear cadenas de comandos

Estas técnicas son fundamentales para escribir flujos de trabajo de Nextflow que sean limpios, efectivos y fáciles de mantener.

---

## 2. Procesamiento de cadenas para manejar nombres de archivos y metadatos

El procesamiento de cadenas es una tarea común en los flujos de trabajo bioinformáticos, especialmente al extraer metadatos de los nombres de archivos o generar scripts dinámicamente. Nextflow tiene potentes capacidades para manipular cadenas que son cruciales para flujos de trabajo robustos.

### 2.1. Expresiones regulares para analizar nombres de archivos

Ya usamos expresiones regulares para extraer metadatos de los nombres de los archivos FASTQ. Vamos a entender mejor cómo funciona el operador de coincidencia de patrones `=~`.

Cuando escribes `x =~ /pattern/`:

1. Crea un objeto `java.util.regex.Matcher`
2. Si se evalúa en un contexto booleano, comprueba si hay **alguna coincidencia**
3. Si se asigna a una variable, puedes acceder a los **grupos de captura**

Veamos cómo podríamos complicar nuestro regex de fastq para manejar más variaciones:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

Ahora, `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` coincidirá con tanto `.fastq` como `.fastq.gz`.

!!! note "Grupos de captura en regex" - Los paréntesis `(...)` definen "grupos de captura" - `m[0]` es la coincidencia completa - `m[0][1]`, `m[0][2]`, etc. son los grupos de captura (comenzando por 1)

Si tuvieras nombres de archivos con diferentes convenciones, podrías usar el operador OR (`|`) en tu regex:

```groovy title="regex example"
def m = (filename =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/ |
                      /^(.+)_(\d{6})_([ACGT]+)_L(\d{3})_(R[12])\.fastq(?:\.gz)?$/)
```

### 2.2. Generación dinámica de scripts

Otra aplicación poderosa del procesamiento de cadenas es generar dinámicamente scripts bash basados en metadatos o entradas. Esto es especialmente útil para lógica condicional dentro de los procesos.

Modifiquemos el proceso `GENERATE_REPORT` para generar diferentes tipos de informes basados en la prioridad:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="3-11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "HIGH PRIORITY SAMPLE" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Standard Sample" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-13"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Ejecutemos el flujo de trabajo para ver el resultado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [condescending_northcutt] DSL2 - revision: 1a3c16a96f

    executor >  local (6)
    [95/e633a2] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [a8/e4c214] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Revisemos los archivos de informe generados - busquemos específicamente el de alta prioridad:

```console
find work -name "sample_003_report.txt" -exec cat {} \;
```

Deberías ver el encabezado especial:

```
HIGH PRIORITY SAMPLE
===============================================
Sample ID: sample_003
Organism: human
...
```

### 2.3. Interpolación de variables: cuándo Nextflow evalúa vs. cuándo bash evalúa

Un punto sutil pero crucial para entender es cuándo ocurre la interpolación de variables:

1. `${var}` - Interpolado por Nextflow durante la compilación del script
2. `\${var}` - Escapado, pasado literalmente a bash como `${var}` (para variables de entorno bash)

Podemos ver esto en acción actualizando el proceso `GENERATE_REPORT` para usar una variable de entorno:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="15"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "HIGH PRIORITY SAMPLE" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Standard Sample" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "HIGH PRIORITY SAMPLE" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Standard Sample" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Si ejecutamos el flujo de trabajo tal como está, fallará:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

El problema es que Nextflow intenta interpretar `${USER}` como una variable de Nextflow, pero no hay ninguna variable llamada `USER`. Necesitamos escapar el `$` para que se pase a bash, que sí tiene una variable de entorno `USER`:

```groovy title="modules/generate_report.nf" linenums="15" hl_lines="1"
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
```

¡Ahora funciona! La barra invertida (`\`) le dice a Nextflow "no interpretes esto, pásalo a Bash."

### Conclusión

En esta sección, has aprendido técnicas de **procesamiento de cadenas**:

- **Expresiones regulares para análisis de archivos**: Usando el operador `=~` y patrones regex (`~/patrón/`) para extraer metadatos de convenciones de nomenclatura de archivos complejos
- **Generación dinámica de scripts**: Usando lógica condicional (if/else, operadores ternarios) para generar diferentes cadenas de script basadas en características de entrada
- **Interpolación de variables**: Entender cuándo Nextflow interpreta cadenas vs cuándo lo hace el shell
  - `${var}` - Variables de Nextflow (interpoladas por Nextflow en tiempo de compilación del flujo de trabajo)
  - `\${var}` - Variables de entorno de shell (escapadas, pasadas a bash en tiempo de ejecución)
  - `\$(cmd)` - Sustitución de comandos de shell (escapada, ejecutada por bash en tiempo de ejecución)

Estos patrones de procesamiento y generación de cadenas son esenciales para manejar los diversos formatos de archivo y convenciones de nomenclatura que encontrarás en flujos de trabajo bioinformáticos del mundo real.

---

## 3. Creando Funciones Reutilizables

La lógica compleja del flujo de trabajo en línea en los operadores de canales o definiciones de procesos reduce la legibilidad y la mantenibilidad. Las **funciones** te permiten extraer esta lógica en componentes nombrados y reutilizables.

Nuestra operación map ha crecido larga y compleja. Extraigámosla en una función reutilizable usando la palabra clave `def`.

Para ilustrar cómo se ve eso con nuestro flujo de trabajo existente, haz la modificación a continuación, usando `def` para definir una función reutilizable llamada `separateMetadata`:

=== "After"

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

=== "Before"

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

Al extraer esta lógica en una función, hemos reducido la lógica real del flujo de trabajo a algo mucho más limpio:

```groovy title="flujo de trabajo mínimo"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Esto hace que la lógica del flujo de trabajo sea mucho más fácil de leer y entender de un vistazo. La función `separateMetadata` encapsula toda la lógica compleja para analizar y enriquecer metadatos, haciéndola reutilizable y testeable.

Ejecuta el flujo de trabajo para asegurarte de que todavía funciona:

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

La salida debería mostrar que ambos procesos se completan con éxito. El flujo de trabajo es ahora mucho más limpio y fácil de mantener, con toda la lógica compleja de procesamiento de metadatos encapsulada en la función `separateMetadata`.

### Conclusión

En esta sección, has aprendido sobre **creación de funciones**:

- **Definiendo funciones con `def`**: La palabra clave para crear funciones nombradas (como `def` en Python o `function` en JavaScript)
- **Ámbito de las funciones**: Las funciones definidas a nivel del script son accesibles a lo largo de tu flujo de trabajo de Nextflow
- **Valores de retorno**: Las funciones automáticamente devuelven la última expresión, o usa `return` explícitamente
- **Código más limpio**: Extraer lógica compleja en funciones es una práctica fundamental de ingeniería de software en cualquier lenguaje

A continuación, exploraremos cómo usar closures en directivas de procesos para asignación dinámica de recursos.

---

## 4. Directivas de Recursos Dinámicas con Closures

Hasta ahora hemos usado scripting en el bloque `script` de los procesos. Pero las **closures** (introducidas en la Sección 1.1) también son increíblemente útiles en las directivas de los procesos, especialmente para la asignación dinámica de recursos. Agreguemos directivas de recursos a nuestro proceso FASTP que se adapten según las características de la muestra.

### 4.1. Asignación de recursos específicos para cada muestra

Actualmente, nuestro proceso FASTP utiliza recursos predeterminados. Hagámoslo más inteligente asignando más CPUs para muestras de alta profundidad. Edita `modules/fastp.nf` para incluir una directiva `cpus` dinámica y una directiva `memory` estática:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

La closure `{ meta.depth > 40000000 ? 2 : 1 }` usa el **operador ternario** (cubierto en la Sección 1.1) y se evalúa para cada tarea, permitiendo la asignación de recursos por muestra. Las muestras de alta profundidad (>40M lecturas) obtienen 2 CPUs, mientras que otras obtienen 1 CPU.

!!! note "Accediendo a Variables de Entrada en Directivas"

    La closure puede acceder a cualquier variable de entrada (como `meta` aquí) porque Nextflow evalúa estas closures en el contexto de cada ejecución de tarea.

Ejecuta el flujo de trabajo de nuevo con la opción `-ansi-log false` para que sea más fácil ver los hashes de las tareas.

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

Puedes verificar el comando exacto de `docker` que se ejecutó para ver la asignación de CPU para cualquier tarea dada:

```console title="Verifica el comando docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Deberías ver algo como:

```bash title="comando docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

En este ejemplo hemos elegido un ejemplo que solicitó 2 CPUs (`--cpu-shares 2048`), porque era una muestra de alta profundidad, pero deberías ver diferentes asignaciones de CPU dependiendo de la profundidad de la muestra. Prueba esto para las otras tareas también.

### 4.2. Estrategias de reintento

Otro patrón poderoso es usar `task.attempt` para estrategias de reintento. Para mostrar por qué esto es útil, vamos a comenzar reduciendo la asignación de memoria a FASTP a menos de lo que necesita. Cambia la directiva `memory` en `modules/fastp.nf` a `1.GB`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... y ejecuta el flujo de trabajo de nuevo:

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

Este es un escenario muy común en flujos de trabajo del mundo real - a veces simplemente no sabes cuánta memoria necesitará una tarea hasta que la ejecutes.

Para hacer nuestro flujo de trabajo más robusto, podemos implementar una estrategia de reintento que aumente la asignación de memoria en cada intento, una vez más usando una closure de Groovy. Modifica la directiva `memory` para multiplicar la memoria base por `task.attempt`, y agrega directivas `errorStrategy 'retry'` y `maxRetries 2`:

=== "After"

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

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Ahora si el proceso falla debido a memoria insuficiente, Nextflow lo reintentará con más memoria:

- Primer intento: 1 GB (task.attempt = 1)
- Segundo intento: 2.GB (task.attempt = 2)

... y así sucesivamente, hasta el límite de `maxRetries`.

### Conclusión

Las directivas dinámicas con closures te permiten:

- Asignar recursos basados en características de la entrada
- Implementar estrategias automáticas de reintento con recursos crecientes
- Combinar múltiples factores (metadatos, número de intento, prioridades)
- Usar lógica condicional para cálculos complejos de recursos

Esto hace que tus flujos de trabajo sean tanto más eficientes (sin sobreasignar) como más robustos (reintento automático con más recursos).

---

## 5. Lógica Condicional y Control de Procesos

Anteriormente, usamos `.map()` con scripting para transformar datos de canales. Ahora usaremos lógica condicional para controlar qué procesos ejecutar basados en los datos, lo cual es esencial para flujos de trabajo flexibles que se adapten a diferentes tipos de muestras.

Los [operadores de flujo de datos](https://nextflow.io/docs/latest/reference/operator.html) de Nextflow toman closures evaluadas en tiempo de ejecución, permitiendo lógica condicional para dirigir decisiones del flujo de trabajo basadas en el contenido del canal.

### 5.1. Enrutamiento con `.branch()`

Por ejemplo, supongamos que nuestras muestras de secuenciación necesitan ser recortadas con FASTP solo si son muestras humanas con una cobertura por encima de cierto umbral. Las muestras de ratón o de baja cobertura deberían ejecutarse con Trimgalore en su lugar (este es un ejemplo artificioso, pero ilustra el punto).

Hemos proporcionado un proceso simple de Trimgalore en `modules/trimgalore.nf`, échale un vistazo si quieres, pero los detalles no son importantes para este ejercicio. El punto clave es que queremos enrutar las muestras basándonos en sus metadatos.

Incluye el nuevo desde `modules/trimgalore.nf`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... y luego modifica el flujo de trabajo `main.nf` para ramificar las muestras basadas en sus metadatos y enrutarlas a través del proceso de recorte apropiado, así:

=== "After"

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

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Ejecuta este flujo de trabajo modificado:

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

Aquí, hemos usado pequeñas pero poderosas expresiones condicionales dentro del operador `.branch{}` para enrutar las muestras basadas en sus metadatos. Las muestras humanas con alta cobertura pasan por `FASTP`, mientras que todas las demás muestras pasan por `TRIMGALORE`.

### 5.2. Usando `.filter()` con Veracidad

Otro patrón poderoso para controlar la ejecución del flujo de trabajo es el operador `.filter()`, que usa una closure para determinar qué elementos deberían continuar por el pipeline. Dentro de la closure de filtrado, escribirás **expresiones booleanas** que deciden qué elementos pasan.

Nextflow (como muchos lenguajes dinámicos) tiene un concepto de **"veracidad"** que determina qué valores evalúan a `true` o `false` en contextos booleanos:

- **Verdaderos**: Valores no nulos, cadenas no vacías, números no cero, colecciones no vacías
- **Falsos**: `null`, cadenas vacías `""`, cero `0`, colecciones vacías `[]` o `[:]`, `false`

Esto significa que `meta.id` solo (sin un explícito `!= null`) comprueba si el ID existe y no está vacío. Usemos esto para filtrar las muestras que no cumplen con nuestros requisitos de calidad.

Agrega lo siguiente antes de la operación branch:

=== "After"

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

=== "Before"

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

Ejecuta el flujo de trabajo de nuevo:

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

La expresión de filtrado `meta.id && meta.organism && meta.depth >= 25000000` combina veracidad con comparaciones explícitas:

- `meta.id && meta.organism` comprueba que ambos campos existan y no estén vacíos (usando veracidad)
- `meta.depth >= 25000000` asegura una profundidad de secuenciación suficiente con una comparación explícita

!!! note "Veracidad en la Práctica"

    La expresión `meta.id && meta.organism` es más concisa que escribir:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Esto hace que la lógica de filtrado sea mucho más limpia y fácil de leer.

### Conclusión

En esta sección, has aprendido a usar lógica condicional para controlar la ejecución del flujo de trabajo usando las interfaces de closure de operadores Nextflow como `.branch{}` y `.filter{}`, aprovechando la veracidad para escribir expresiones condicionales concisas.

Nuestro pipeline ahora enruta inteligentemente las muestras a través de los procesos apropiados, pero los flujos de trabajo de producción necesitan manejar datos inválidos con elegancia. Hagamos que nuestro flujo de trabajo sea robusto contra valores faltantes o nulos.

---

## 6. Operadores de Navegación Segura y Elvis

Nuestra función `separateMetadata` actualmente asume que todos los campos CSV están presentes y válidos. Pero, ¿qué sucede con datos incompletos? Averigüémoslo.

### 6.1. El Problema: Acceder a Propiedades que No Existen

Digamos que queremos agregar soporte para información opcional del experimento de secuenciación. En algunos laboratorios, las muestras pueden tener un campo adicional para el ID o número de lote del experimento de secuenciación, pero nuestro CSV actual no tiene esta columna. Intentemos acceder a ella de todos modos.

Modifica la función `separateMetadata` para incluir un campo run_id:

=== "After"

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

=== "Before"

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

Ahora ejecuta el flujo de trabajo:

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

El problema es que `row.run_id` devuelve `null` porque la columna `run_id` no existe en nuestro CSV. Cuando tratamos de llamar a `.toUpperCase()` en `null`, falla. Aquí es donde el operador de navegación segura viene al rescate.

### 6.2. Operador de Navegación Segura (`?.`)

El operador de navegación segura (`?.`) devuelve `null` en lugar de lanzar una excepción cuando se llama a un valor `null`. Si el objeto antes de `?.` es `null`, toda la expresión evalúa a `null` sin ejecutar el método.

Actualiza la función para usar navegación segura:

=== "After"

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

=== "Before"

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

Ejecuta de nuevo:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

¡Sin fallo! El flujo de trabajo ahora maneja el campo faltante con elegancia. Cuando `row.run_id` es `null`, el operador `?.` previene la llamada `.toUpperCase()`, y `run_id` se convierte en `null` en vez de causar una excepción.

### 6.3. Operador Elvis (`?:`) para Valores Predeterminados

El operador Elvis (`?:`) proporciona valores predeterminados cuando el lado izquierdo es "falso" (como se explicó anteriormente). Se llama así por Elvis Presley porque `?:` se parece a su famoso cabello y ojos cuando se ve de lado.

Ahora que estamos usando navegación segura, `run_id` será `null` para las muestras sin ese campo. Usemos el operador Elvis para proporcionar un valor predeterminado y agregarlo a nuestro mapa `sample_meta`:

=== "After"

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

=== "Before"

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

También agrega un operador `view()` en el flujo de trabajo para ver los resultados:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

y ejecuta el flujo de trabajo:

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

Quita el operador `.view()` ahora que hemos confirmado que funciona.

!!! tip "Combinando Navegación Segura y Elvis"

    El patrón `value?.method() ?: 'default'` es común en flujos de trabajo de producción:

    - `value?.method()` - Llama al método con seguridad, devuelve `null` si `value` es `null`
    - `?: 'default'` - Proporciona un valor alternativo si el resultado es `null`

    Este patrón maneja datos faltantes/incompletos con elegancia.

Usa estos operadores consistentemente en funciones, closures de operadores (`.map{}`, `.filter{}`), scripts de procesos, y archivos de configuración. Previenen fallos al manejar datos del mundo real.

### Conclusión

- **Navegación segura (`?.`)**: Previene fallos en valores nulos - devuelve nulo en lugar de lanzar excepción
- **Operador Elvis (`?:`)**: Proporciona valores predeterminados - `value ?: 'default'`
- **Combinando**: `value?.method() ?: 'default'` es el patrón común

Estos operadores hacen que los flujos de trabajo sean resistentes a datos incompletos - esencial para el trabajo en el mundo real.

---

## 7. Validación con `error()` y `log.warn`

A veces necesitas detener el flujo de trabajo inmediatamente si los parámetros de entrada son inválidos. En Nextflow, puedes usar funciones incorporadas como `error()` y `log.warn`, así como construcciones de programación estándar como declaraciones `if` y lógica booleana, para implementar lógica de validación. Agreguemos validación a nuestro flujo de trabajo.

Crea una función de validación antes de tu bloque de flujo de trabajo, llámala desde el flujo de trabajo, y cambia la creación del canal para usar un parámetro para la ruta del archivo CSV. Si el parámetro falta o el archivo no existe, llama a `error()` para detener la ejecución con un mensaje claro.

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Verifica que se proporcione el parámetro de entrada
        if (!params.input) {
            error("Ruta del archivo CSV de entrada no proporcionada. Por favor, especifique --input <file.csv>")
        }

        // Verifica que el archivo CSV exista
        if (!file(params.input).exists()) {
            error("Archivo CSV de entrada no encontrado: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Ahora intenta ejecutar sin el archivo CSV:

```bash
nextflow run main.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Ruta del archivo CSV de entrada no proporcionada. Por favor, especifique --input <file.csv>
    ```

El flujo de trabajo se detiene inmediatamente con un mensaje de error claro en lugar de fallar misteriosamente más tarde

Ahora ejecútalo con un archivo inexistente:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Archivo CSV de entrada no encontrado: ./data/nonexistent.csv
    ```

Finalmente, ejecútalo con el archivo correcto:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Salida del comando"

    ```console
    <!-- TODO: output -->
    ```

Esta vez se ejecuta exitosamente.

También puedes agregar validación dentro de la función `separateMetadata`. Usemos el no fatal `log.warn` para emitir advertencias para muestras con baja profundidad de secuenciación, pero aun así permitir que el flujo de trabajo continúe:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validar que los datos tienen sentido
        if (sample_meta.depth < 30000000) {
            log.warn "Baja profundidad de secuenciación para ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Ejecuta el flujo de trabajo de nuevo con el CSV original:

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
    WARN: Baja profundidad de secuenciación para sample_002: 25000000
    ```

Vemos una advertencia sobre la baja profundidad de secuenciación para una de las muestras.

### Conclusión

- **`error()`**: Detiene el flujo de trabajo inmediatamente con un mensaje claro
- **`log.warn`**: Emite advertencias sin detener el flujo de trabajo
- **Validación temprana**: Verifica las entradas antes de procesarlas para fallar rápido con errores útiles
- **Funciones de validación**: Crea lógica de validación reutilizable que puede ser llamada al inicio del flujo de trabajo

La validación adecuada hace que los flujos de trabajo sean más robustos y amigables para el usuario al detectar problemas temprano con mensajes de error claros.

---

## 8. Manejadores de Eventos del Flujo de Trabajo

Hasta ahora, hemos estado escribiendo código en nuestros scripts de flujo de trabajo y definiciones de proceso. Pero hay una característica más importante que deberías conocer: manejadores de eventos de flujo de trabajo.

Los manejadores de eventos son closures que se ejecutan en puntos específicos del ciclo de vida de tu flujo de trabajo. Son perfectos para agregar registro, notificaciones, u operaciones de limpieza. Estos manejadores deben ser definidos en tu script de flujo de trabajo junto a tu definición de flujo de trabajo.

### 8.1. El Manejador `onComplete`

El manejador de eventos más comúnmente usado es `onComplete`, que se ejecuta cuando tu flujo de trabajo termina (ya sea que haya tenido éxito o fallado). Agreguemos uno para resumir los resultados de nuestro pipeline.

Agrega el manejador de eventos a tu archivo `main.nf`, dentro de tu definición de flujo de trabajo:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Resumen de ejecución del pipeline:"
            println "=========================="
            println "Completado en: ${workflow.complete}"
            println "Duración    : ${workflow.duration}"
            println "Éxito       : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "estado salida : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Esta closure se ejecuta cuando el flujo de trabajo se completa. Dentro, tienes acceso al objeto `workflow` que proporciona propiedades útiles sobre la ejecución.

¡Ejecuta tu flujo de trabajo y verás este resumen aparecer al final!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Baja profundidad de secuenciación para sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Resumen de ejecución del pipeline:
    ==========================
    Completado en: 2025-10-10T12:14:24.885384+01:00
    Duración    : 2.9s
    Éxito       : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    estado salida : 0
    ```

Hagámoslo más útil agregando lógica condicional:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Resumen de ejecución del pipeline:"
            println "=========================="
            println "Completado en: ${workflow.complete}"
            println "Duración    : ${workflow.duration}"
            println "Éxito       : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "estado salida : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ ¡Pipeline completado exitosamente!"
            } else {
                println "❌ ¡Pipeline falló!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Resumen de ejecución del pipeline:"
            println "=========================="
            println "Completado en: ${workflow.complete}"
            println "Duración    : ${workflow.duration}"
            println "Éxito       : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "estado salida : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Ahora obtenemos un resumen aún más informativo, incluyendo un mensaje de éxito/fallo y el directorio de salida si se especifica:

<!-- TODO: add run command -->

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Baja profundidad de secuenciación para sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Resumen de ejecución del pipeline:
    ==========================
    Completado en: 2025-10-10T12:16:00.522569+01:00
    Duración    : 3.6s
    Éxito       : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    estado salida : 0

    ✅ ¡Pipeline completado exitosamente!
    ```

También puedes escribir el resumen en un archivo usando operaciones de archivo:

```groovy title="main.nf - Escribiendo resumen a un archivo"
workflow {
    // ... tu código de flujo de trabajo ...

    workflow.onComplete = {
        def summary = """
        Resumen de Ejecución del Pipeline
        ===========================
        Completado: ${workflow.complete}
        Duración : ${workflow.duration}
        Éxito  : ${workflow.success}
        Comando  : ${workflow.commandLine}
        """

        println summary

        // Escribir a un archivo de registro
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. El Manejador `onError`

Además de `onComplete`, hay otro manejador de eventos que puedes usar: `onError`, que se ejecuta solo si el flujo de trabajo falla:

```groovy title="main.nf - manejador onError"
workflow {
    // ... tu código de flujo de trabajo ...

    workflow.onError = {
        println "="* 50
        println "¡Ejecución del pipeline falló!"
        println "Mensaje de error: ${workflow.errorMessage}"
        println "="* 50

        // Escribir registro de error detallado
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Informe de Error del Flujo de Trabajo
        =====================
        Tiempo: ${new Date()}
        Error: ${workflow.errorMessage}
        Informe de error: ${workflow.errorReport ?: 'No hay informe detallado disponible'}
        """

        println "Detalles del error escritos en: ${error_file}"
    }
}
```

Puedes usar múltiples manejadores juntos en tu script de flujo de trabajo:

```groovy title="main.nf - Manejadores combinados"
workflow {
    // ... tu código de flujo de trabajo ...

    workflow.onError = {
        println "Flujo de trabajo falló: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "ÉXITO ✅" : "FALLÓ ❌"

        println """
        Pipeline finalizado: ${status}
        Duración: ${duration_mins} minutos
        """
    }
}
```

### Conclusión

En esta sección, has aprendido:

- **Closures de manejadores de eventos**: Closures en tu script de flujo de trabajo que se ejecutan en diferentes puntos del ciclo de vida
- **Manejador `onComplete`**: Para resúmenes de ejecución y reporte de resultados
- **Manejador `onError`**: Para manejo de errores y registro de fallos
- **Propiedades del objeto workflow**: Accediendo a `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Los manejadores de eventos muestran cómo puedes usar todo el poder del lenguaje Nextflow dentro de tus scripts de flujo de trabajo para agregar capacidades sofisticadas de registro y notificación.

---

## Resumen

¡Felicitaciones, lo lograste!

A lo largo de esta búsqueda lateral, has construido un pipeline completo de procesamiento de muestras que evolucionó desde el manejo básico de metadatos hasta un flujo de trabajo sofisticado y listo para producción.
Cada sección se basó en la anterior, demostrando cómo las construcciones de programación transforman flujos de trabajo simples en sistemas potentes de procesamiento de datos, con los siguientes beneficios:

- **Código más claro**: Entender el flujo de datos vs scripting te ayuda a escribir flujos de trabajo más organizados
- **Manejo robusto**: Los operadores de navegación segura y Elvis hacen que los flujos de trabajo sean resistentes a datos faltantes
- **Procesamiento flexible**: La lógica condicional permite que tus flujos de trabajo procesen diferentes tipos de muestras apropiadamente
- **Recursos adaptativos**: Las directivas dinámicas optimizan el uso de recursos basado en características de entrada

Esta progresión refleja la evolución en el mundo real de los pipelines bioinformáticos, desde prototipos de investigación que manejan unas pocas muestras hasta sistemas de producción que procesan miles de muestras en laboratorios e instituciones.
Cada desafío que resolviste y patrón que aprendiste refleja problemas reales que los desarrolladores enfrentan al escalar flujos de trabajo de Nextflow.

Aplicar estos patrones en tu propio trabajo te permitirá construir flujos de trabajo robustos y listos para producción.

### Patrones clave

1.  **Flujo de Datos vs Scripting:** Aprendiste a distinguir entre operaciones de flujo de datos (orquestación de canales) y scripting (código que manipula datos), incluyendo las diferencias cruciales entre operaciones en diferentes tipos como `collect` en Channel vs List.

    - Flujo de datos: orquestación de canales

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: procesamiento de datos en colecciones

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Procesamiento Avanzado de Cadenas**: Dominaste expresiones regulares para analizar nombres de archivos, generación dinámica de scripts en procesos, e interpolación de variables (Nextflow vs Bash vs Shell).

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

    - Colección de archivos a argumentos de comando (en bloque de script de proceso)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Creando Funciones Reutilizables**: Aprendiste a extraer lógica compleja en funciones nombradas que pueden ser llamadas desde operadores de canales, haciendo los flujos de trabajo más legibles y mantenibles.

    - Define una función nombrada

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* código oculto por brevedad */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* código oculto por brevedad */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Llama a la función nombrada en un flujo de trabajo

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Directivas de Recursos Dinámicas con Closures**: Exploraste el uso de closures en directivas de proceso para asignación adaptativa de recursos basada en características de entrada.

    - Closures nombradas y composición

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures con acceso a ámbito

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Lógica Condicional y Control de Proceso**: Agregaste enrutamiento inteligente usando operadores `.branch()` y `.filter()`, aprovechando la veracidad para expresiones condicionales concisas.

    - Usa `.branch()` para enrutar datos a través de diferentes ramas del flujo de trabajo

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Evaluación booleana con Veracidad de Groovy

    ```groovy
    if (sample.files) println "Tiene archivos"
    ```

    - Usa `filter()` para filtrar datos con 'veracidad'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operadores de Navegación Segura y Elvis**: Hiciste el pipeline robusto contra datos faltantes usando `?.` para acceso seguro a propiedades y `?:` para proporcionar valores predeterminados.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validación con error() y log.warn**: Aprendiste a validar entradas temprano y fallar rápido con mensajes de error claros.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Inválido: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Manejadores de Eventos de Configuración**: Aprendiste a usar manejadores de eventos de flujo de trabajo (`onComplete` y `onError`) para registro, notificaciones, y gestión del ciclo de vida.

    - Usando `onComplete` para registrar y notificar

    ```groovy
    workflow.onComplete = {
        println "Éxito       : ${workflow.success}"
        println "estado salida : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ ¡Pipeline completado exitosamente!"
        } else {
            println "❌ ¡Pipeline falló!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Usando `onError` para tomar acción específicamente en caso de fallo

    ```groovy
    workflow.onError = {
        // Escribir registro de error detallado
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Tiempo: ${new Date()}
        Error: ${workflow.errorMessage}
        Informe de error: ${workflow.errorReport ?: 'No hay informe detallado disponible'}
        """

        println "Detalles del error escritos en: ${error_file}"
    }
    ```

### Recursos adicionales

- [Referencia del Lenguaje Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operadores de Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Sintaxis de Script de Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteca Estándar de Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Asegúrate de consultar estos recursos cuando necesites explorar características más avanzadas.

Te beneficiarás de practicar y expandir tus habilidades para:

- Escribir flujos de trabajo más limpios con separación adecuada entre flujo de datos y scripting
- Dominar la interpolación de variables para evitar escollos comunes con variables de Nextflow, Bash, y shell
- Usar directivas de recursos dinámicas para flujos de trabajo eficientes y adaptativos
- Transformar colecciones de archivos en argumentos de línea de comandos formateados adecuadamente
- Manejar diferentes convenciones de nombres de archivos y formatos de entrada con gracia usando regex y procesamiento de cadenas
- Construir código reutilizable y mantenible usando patrones avanzados de closure y programación funcional
- Procesar y organizar conjuntos de datos complejos usando operaciones de colección
- Agregar validación, manejo de errores, y registro para hacer tus flujos de trabajo listos para producción
- Implementar gestión del ciclo de vida del flujo de trabajo con manejadores de eventos

---

## ¿Qué sigue?

Regresa al [menú de Búsquedas Laterales](./index.md) o haz clic en el botón en la parte inferior derecha de la página para pasar al siguiente tema de la lista.
