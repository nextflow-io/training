# Parte 3: Implementación con múltiples muestras de extremos pareados

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Anteriormente, construyó un pipeline de llamado de variantes por muestra que procesaba los datos de cada muestra de forma independiente.
En esta parte del curso, vamos a llevar nuestro flujo de trabajo simple al siguiente nivel convirtiéndolo en una poderosa herramienta de automatización por lotes para manejar números arbitrarios de muestras.
Y mientras lo hacemos, también vamos a actualizarlo para que espere datos de extremos pareados, que son más comunes en estudios más recientes.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado [Parte 1: Descripción del método](./01_method.md), [Parte 2: Implementación de una sola muestra](./02_single-sample.md) y tiene un pipeline `rnaseq.nf` funcional con archivos de módulos completos.

    Si no completó la Parte 2 o desea comenzar de nuevo para esta parte, puede usar la solución de la Parte 2 como punto de partida.
    Ejecute estos comandos desde dentro del directorio `nf4-science/rnaseq/`:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Esto le proporciona un flujo de trabajo completo de procesamiento de una sola muestra.
    Puede probar que se ejecuta correctamente:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Asignación

En esta parte del curso, vamos a extender el flujo de trabajo para hacer lo siguiente:

1. Leer información de muestras desde una hoja de cálculo CSV
2. Ejecutar control de calidad, recorte y alineamiento por muestra en todas las muestras en paralelo
3. Agregar todos los reportes de control de calidad en un reporte MultiQC completo

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Esto automatiza los pasos de la segunda sección de [Parte 1: Descripción del método](./01_method.md#2-multi-sample-qc-aggregation), donde ejecutó estos comandos manualmente en sus contenedores.

## Plan de lección

Hemos dividido esto en tres etapas:

1. **Hacer que el flujo de trabajo acepte múltiples muestras de entrada.**
   Esto cubre el cambio de una sola ruta de archivo a una hoja de cálculo CSV, analizarla con `splitCsv()`, y ejecutar todos los procesos existentes en múltiples muestras.
2. **Agregar generación de reportes de control de calidad completos.**
   Esto introduce el operador `collect()` para agregar salidas a través de muestras, y agrega un proceso MultiQC para producir un reporte combinado.
3. **Cambiar a datos de RNAseq de extremos pareados.**
   Esto cubre la adaptación de procesos para entradas de extremos pareados (usando tuplas), creación de módulos de extremos pareados, y configuración de un perfil de prueba separado.

Esto implementa el método descrito en [Parte 1: Descripción del método](./01_method.md) (segunda sección que cubre el caso de uso de múltiples muestras) y se construye directamente sobre el flujo de trabajo producido por la Parte 2.

!!! tip "Consejo"

     Asegúrese de estar en el directorio de trabajo correcto:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Hacer que el flujo de trabajo acepte múltiples muestras de entrada

Para ejecutar en múltiples muestras, necesitamos cambiar cómo gestionamos la entrada: en lugar de proporcionar una sola ruta de archivo, leeremos información de muestras desde un archivo CSV.

Proporcionamos un archivo CSV que contiene IDs de muestra y rutas de archivos FASTQ en el directorio `data/`.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Este archivo CSV incluye una línea de encabezado que nombra las columnas.

Observe que estos son todavía datos de lecturas de un solo extremo.

!!! warning "Advertencia"

    Las rutas de archivos en el CSV son rutas absolutas que deben coincidir con su entorno.
    Si no está ejecutando esto en el entorno de capacitación que proporcionamos, necesitará actualizar las rutas para que coincidan con su sistema.

### 1.1. Cambiar la entrada principal a un CSV de rutas de archivos en el perfil de prueba

Primero, necesitamos actualizar el perfil de prueba en `nextflow.config` para proporcionar la ruta del archivo CSV en lugar de la ruta FASTQ única.

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

A continuación, necesitaremos actualizar la creación del canal para leer desde este CSV.

### 1.2. Actualizar la fábrica de canales para analizar entrada CSV

Necesitamos cargar los contenidos del archivo en el canal en lugar de solo la ruta del archivo.

Podemos hacer esto usando el mismo patrón que usamos en [Parte 2 de Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): aplicando el operador [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) para analizar el archivo, luego una operación `map` para extraer la ruta del archivo FASTQ de cada fila.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Crear canal de entrada a partir del contenido de un archivo CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)
    ```

Una cosa que es nueva comparada con lo que encontró en el curso Hello Nextflow es que este CSV tiene una línea de encabezado, así que agregamos `#!groovy header: true` a la llamada `splitCsv()`.
Eso nos permite referenciar columnas por nombre en la operación `map`: `#!groovy row.fastq_path` extrae la ruta del archivo de la columna `fastq_path` de cada fila.

El manejo de entrada está actualizado y el flujo de trabajo está listo para probar.

### 1.3. Ejecutar el flujo de trabajo

El flujo de trabajo ahora lee información de muestras desde un archivo CSV y procesa todas las muestras en paralelo.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Esta vez cada paso se ejecuta 6 veces, una vez para cada muestra en el archivo CSV.

Eso es todo lo que se necesitó para que el flujo de trabajo se ejecute en múltiples archivos.
Nextflow maneja todo el paralelismo por nosotros.

### Conclusión

Usted sabe cómo cambiar de una entrada de archivo único a entrada de múltiples muestras basada en CSV que Nextflow procesa en paralelo.

### ¿Qué sigue?

Agregar un paso de agregación de reportes de control de calidad que combine métricas de todas las muestras.

---

## 2. Agregar métricas de control de calidad de preprocesamiento en un solo reporte MultiQC

Todo esto produce muchos reportes de control de calidad, y no queremos tener que revisar reportes individuales.
Este es el punto perfecto para poner un paso de agregación de reportes MultiQC.

Recuerde el comando `multiqc` de [Parte 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

El comando escanea el directorio actual en busca de archivos de salida de control de calidad reconocidos y los agrega en un solo reporte HTML.
El URI del contenedor era `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Necesitamos configurar un parámetro adicional, preparar las entradas, escribir el proceso, conectarlo, y actualizar el manejo de salida.

### 2.1. Configurar las entradas

El proceso MultiQC necesita un parámetro de nombre de reporte y las salidas de control de calidad recopiladas de todos los pasos anteriores agrupadas juntas.

#### 2.1.1. Agregar un parámetro `report_id`

Agregue un parámetro para nombrar el reporte de salida.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Entrada principal
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path

        // Report ID
        report_id: String
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Entrada principal
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

Agregue el ID de reporte predeterminado al perfil de prueba:

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

A continuación, necesitaremos preparar las entradas para el proceso MultiQC.

#### 2.1.2. Recopilar y combinar salidas de control de calidad de pasos anteriores

Necesitamos darle al proceso `MULTIQC` todas las salidas relacionadas con control de calidad de pasos anteriores agrupadas juntas.

Para eso, usamos el operador `.mix()`, que agrega múltiples canales en uno solo.
Comenzamos desde `channel.empty()` y mezclamos todos los canales de salida que queremos combinar.
Esto es más limpio que encadenar `.mix()` directamente en uno de los canales de salida, porque trata todas las entradas simétricamente.

En nuestro flujo de trabajo, las salidas relacionadas con control de calidad para agregar son:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Las mezclamos en un solo canal, luego usamos `.collect()` para agregar los reportes a través de todas las muestras en una sola lista.

Agregue estas líneas al cuerpo del flujo de trabajo después de la llamada a `HISAT2_ALIGN`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Alineamiento a un genoma de referencia
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Comprehensive QC report generation
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="38"
        // Alineamiento a un genoma de referencia
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

Usar variables intermedias hace que cada paso sea claro: `multiqc_files_ch` contiene todos los archivos de control de calidad individuales mezclados en un canal, y `multiqc_files_list` es el paquete recopilado listo para pasar a MultiQC.

### 2.2. Escribir el proceso de agregación de control de calidad y llamarlo en el flujo de trabajo

Como antes, necesitamos completar la definición del proceso, importar el módulo, y agregar la llamada al proceso.

#### 2.2.1. Completar el módulo para el proceso de agregación de control de calidad

Abra `modules/multiqc.nf` y examine el esquema de la definición del proceso.

Adelante y complete la definición del proceso por usted mismo usando la información proporcionada arriba, luego verifique su trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Este proceso usa `#!groovy path '*'` como el calificador de entrada para los archivos de control de calidad.
El comodín `'*'` le dice a Nextflow que coloque todos los archivos recopilados en el directorio de trabajo sin requerir nombres específicos.
La entrada `val output_name` es una cadena que controla el nombre del archivo del reporte.

El comando `multiqc .` escanea el directorio actual (donde están todos los archivos de control de calidad colocados) y genera el reporte.

Una vez que haya completado esto, el proceso está listo para usar.

#### 2.2.2. Incluir el módulo

Agregue la declaración de importación a `rnaseq.nf`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Declaraciones de inclusión de módulos
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Declaraciones de inclusión de módulos
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Ahora agregue la llamada al proceso al flujo de trabajo.

#### 2.2.3. Agregar la llamada al proceso

Pase los archivos de control de calidad recopilados y el ID de reporte al proceso `MULTIQC`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

El proceso MultiQC ahora está conectado en el flujo de trabajo.

### 2.3. Actualizar el manejo de salida

Necesitamos agregar las salidas de MultiQC a la declaración de publicación y configurar dónde van.

#### 2.3.1. Agregar objetivos de publicación para las salidas de MultiQC

Agregue las salidas de MultiQC a la sección `publish:`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

A continuación, necesitaremos decirle a Nextflow dónde poner estas salidas.

#### 2.3.2. Configurar los nuevos objetivos de salida

Agregue entradas para los objetivos de MultiQC en el bloque `output {}`, publicándolos en un subdirectorio `multiqc/`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

La configuración de salida está completa.

### 2.4. Ejecutar el flujo de trabajo

Usamos `-resume` para que los pasos de procesamiento anteriores estén en caché y solo se ejecute el nuevo paso de MultiQC.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Una sola llamada a MULTIQC se ha agregado después de las llamadas a procesos en caché.

Puede encontrar las salidas de MultiQC en el directorio de resultados.

```bash
tree -L 2 results/multiqc
```

```console title="Salida"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Ese último archivo `all_single-end.html` es el reporte agregado completo, convenientemente empaquetado en un archivo HTML fácil de navegar.

### Conclusión

Usted sabe cómo recopilar salidas de múltiples canales, agruparlas con `.mix()` y `.collect()`, y pasarlas a un proceso de agregación.

### ¿Qué sigue?

Adaptar el flujo de trabajo para manejar datos de RNAseq de extremos pareados.

---

## 3. Habilitar el procesamiento de datos de RNAseq de extremos pareados

Actualmente nuestro flujo de trabajo solo puede manejar datos de RNAseq de un solo extremo.
Es cada vez más común ver datos de RNAseq de extremos pareados, así que queremos poder manejar eso.

Hacer que el flujo de trabajo sea completamente agnóstico del tipo de datos requeriría usar características del lenguaje Nextflow un poco más avanzadas, así que no vamos a hacerlo aquí, pero podemos hacer una versión de procesamiento de extremos pareados para demostrar qué necesita ser adaptado.

### 3.1. Copiar el flujo de trabajo y actualizar las entradas

Comenzamos copiando el archivo del flujo de trabajo de un solo extremo y actualizándolo para datos de extremos pareados.

#### 3.1.1. Copiar el archivo del flujo de trabajo

Cree una copia del archivo del flujo de trabajo para usar como punto de partida para la versión de extremos pareados.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Ahora actualice los parámetros y el manejo de entrada en el nuevo archivo.

#### 3.1.2. Agregar un perfil de prueba de extremos pareados

Proporcionamos un segundo archivo CSV que contiene IDs de muestra y rutas de archivos FASTQ pareados en el directorio `data/`.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Agregue un perfil `test_pe` a `nextflow.config` que apunte a este archivo y use un ID de reporte de extremos pareados.

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

El perfil de prueba para datos de extremos pareados está listo.

#### 3.1.3. Actualizar la fábrica de canales

El operador `.map()` necesita obtener ambas rutas de archivos FASTQ y devolverlas como una lista.

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Crear canal de entrada a partir del contenido de un archivo CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Crear canal de entrada a partir del contenido de un archivo CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

El manejo de entrada está configurado para datos de extremos pareados.

### 3.2. Adaptar el módulo FASTQC para datos de extremos pareados

Copie el módulo para crear una versión de extremos pareados:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

La entrada del proceso FASTQC no necesita cambiar — cuando Nextflow recibe una lista de dos archivos, coloca ambos y `reads` se expande a ambos nombres de archivo.
El único cambio necesario está en el bloque de salida: ya que ahora obtenemos dos reportes FastQC por muestra, cambiamos de patrones basados en `simpleName` a comodines.

=== "Después"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Antes"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Esto generaliza el proceso de una manera que lo hace capaz de manejar datos de un solo extremo o de extremos pareados.

Actualice la importación en `rnaseq_pe.nf` para usar la versión de extremos pareados:

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

El módulo FASTQC y su importación están actualizados para datos de extremos pareados.

### 3.3. Adaptar el módulo TRIM_GALORE para datos de extremos pareados

Copie el módulo para crear una versión de extremos pareados:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Este módulo necesita cambios más sustanciales:

- La entrada cambia de una sola ruta a una tupla de dos rutas
- El comando agrega la bandera `--paired` y toma ambos archivos de lectura
- La salida cambia para reflejar las convenciones de nomenclatura diferentes de Trim Galore para extremos pareados, produciendo reportes FastQC separados para cada archivo de lectura

=== "Después"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Antes"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    ```

Actualice la importación en `rnaseq_pe.nf`:

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

El módulo TRIM_GALORE y su importación están actualizados para datos de extremos pareados.

### 3.4. Adaptar el módulo HISAT2_ALIGN para datos de extremos pareados

Copie el módulo para crear una versión de extremos pareados:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Este módulo necesita cambios similares:

- La entrada cambia de una sola ruta a una tupla de dos rutas
- El comando HISAT2 cambia de `-U` (no pareado) a argumentos de lectura `-1` y `-2` (pareados)
- Todos los usos de `reads.simpleName` cambian a `read1.simpleName` ya que ahora referenciamos un miembro específico del par

=== "Después"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Antes"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    ```

Actualice la importación en `rnaseq_pe.nf`:

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

El módulo HISAT2_ALIGN y su importación están actualizados para datos de extremos pareados.

### 3.5. Actualizar la agregación de MultiQC para salidas de extremos pareados

El proceso `TRIM_GALORE` de extremos pareados ahora produce dos canales de reportes FastQC separados (`fastqc_reports_1` y `fastqc_reports_2`) en lugar de uno.
Actualice el bloque `.mix()` en `rnaseq_pe.nf` para incluir ambos:

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

La agregación de MultiQC ahora incluye ambos conjuntos de reportes FastQC de extremos pareados.

### 3.6. Actualizar el manejo de salida para salidas de extremos pareados

La sección `publish:` y el bloque `output {}` también necesitan reflejar los dos canales de reportes FastQC separados del proceso `TRIM_GALORE` de extremos pareados.

Actualice la sección `publish:` en `rnaseq_pe.nf`:

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Actualice las entradas correspondientes en el bloque `output {}`:

=== "Después"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

El flujo de trabajo de extremos pareados ahora está completamente actualizado y listo para ejecutar.

### 3.7. Ejecutar el flujo de trabajo

No usamos `-resume` ya que esto no usaría la caché, y hay el doble de datos para procesar que antes, pero aún así debería completarse en menos de un minuto.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Ahora tenemos dos versiones ligeramente divergentes de nuestro flujo de trabajo, una para datos de lecturas de un solo extremo y una para datos de extremos pareados.
El siguiente paso lógico sería hacer que el flujo de trabajo acepte cualquier tipo de datos sobre la marcha, lo cual está fuera del alcance de este curso, pero podríamos abordar eso en un seguimiento.

---

### Conclusión

Usted sabe cómo adaptar un flujo de trabajo de una sola muestra para paralelizar el procesamiento de múltiples muestras, generar un reporte de control de calidad completo y adaptar el flujo de trabajo para usar datos de lecturas de extremos pareados.

### ¿Qué sigue?

¡Dese una gran palmada en la espalda! Ha completado el curso de Nextflow para RNAseq.

Diríjase al [resumen del curso](./next_steps.md) final para revisar lo que aprendió y descubrir qué viene después.
