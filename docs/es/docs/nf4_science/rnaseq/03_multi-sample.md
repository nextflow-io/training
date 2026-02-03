# Parte 3: Implementación con múltiples muestras de extremos pareados

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta parte final del curso, vamos a llevar nuestro flujo de trabajo simple al siguiente nivel convirtiéndolo en una poderosa herramienta de automatización por lotes para manejar números arbitrarios de muestras.
Y mientras lo hacemos, también vamos a cambiarlo para que espere datos de extremos pareados, que son más comunes en estudios más recientes.

Lo haremos en tres etapas:

1. Hacer que el flujo de trabajo acepte múltiples muestras de entrada y paralelice la ejecución
2. Agregar generación de reportes de control de calidad completos
3. Cambiar a datos de RNAseq de extremos pareados

---

## 1. Hacer que el flujo de trabajo acepte múltiples muestras de entrada y paralelice la ejecución

Vamos a necesitar cambiar cómo gestionamos la entrada.

### 1.1. Cambiar la entrada principal para que sea un CSV de rutas de archivos en lugar de un solo archivo

Proporcionamos un archivo CSV que contiene IDs de muestra y rutas de archivos FASTQ en el directorio `data/`.
Este archivo CSV incluye una línea de encabezado.
Observe que las rutas de archivos FASTQ son rutas absolutas.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Renombremos el parámetro de entrada principal a `input_csv` y cambiemos el valor predeterminado para que sea la ruta al archivo `single-end.csv`.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Entrada principal
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Actualizar la fábrica de canales de entrada para manejar un CSV como entrada

Vamos a querer cargar los contenidos del archivo en el canal en lugar de solo la ruta del archivo, así que usamos el operador `.splitCsv()` para analizar el formato CSV, luego el operador `.map()` para obtener la información específica que queremos (la ruta del archivo FASTQ).

```groovy title="rnaseq.nf" linenums="16"
    // Crear canal de entrada a partir del contenido de un archivo CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Ejecutar el flujo de trabajo para probar que funciona

```bash
nextflow run rnaseq.nf
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

Esta vez vemos que cada paso se ejecuta 6 veces, en cada uno de los 6 archivos de datos que proporcionamos.

¡Eso es todo lo que se necesitó para que el flujo de trabajo se ejecute en múltiples archivos!
Nextflow maneja todo el paralelismo por nosotros.

---

## 2. Agregar métricas de control de calidad de preprocesamiento en un solo reporte MultiQC

Todo esto produce muchos reportes de control de calidad, y no queremos tener que revisar reportes individuales.
¡Este es el punto perfecto para poner un paso de agregación de reportes MultiQC!

### 2.1. Crear un módulo para el proceso de agregación de control de calidad

Creemos un archivo de módulo llamado `modules/multiqc.nf` para alojar el proceso `MULTIQC`:

```bash
touch modules/multiqc.nf
```

Abra el archivo en el editor de código y copie el siguiente código en él:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. Importar el módulo en el archivo del flujo de trabajo

Agregue la declaración `include { MULTIQC } from './modules/multiqc.nf'` al archivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declaraciones de inclusión de módulos
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Agregar un parámetro `report_id` y darle un valor predeterminado sensato

```groovy title="rnaseq.nf" linenums="9"
params {
    // Entrada principal
    input_csv: Path = "data/single-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 2.4. Llamar al proceso con las salidas de los pasos anteriores

Necesitamos darle al proceso `MULTIQC` todas las salidas relacionadas con control de calidad de los pasos anteriores.

Para eso, vamos a usar el operador `.mix()`, que agrega múltiples canales en uno solo.

Si tuviéramos cuatro procesos llamados A, B, C y D con un canal `.out` simple cada uno, la sintaxis se vería así: `A.out.mix( B.out, C.out, D.out )`. Como puede ver, se aplica al primero de los canales que desea combinar (no importa cuál) y simplemente agrega todos los demás, separados por comas, en el paréntesis que sigue.

En el caso de nuestro flujo de trabajo, tenemos las siguientes salidas para agregar:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Así que el ejemplo de sintaxis se convierte en:

```groovy title="Aplicando .mix() en la llamada a MULTIQC"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

Eso recogerá reportes de control de calidad por muestra.
Pero como queremos agregarlos a través de todas las muestras, necesitamos agregar el operador `collect()` para reunir los reportes de todas las muestras en una sola llamada a `MULTIQC`.
Y también necesitamos darle el parámetro `report_id`.

Esto nos da lo siguiente:

```groovy title="La llamada completa a MULTIQC" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

En el contexto del bloque de flujo de trabajo completo, termina viéndose así:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Crear canal de entrada a partir del contenido de un archivo CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    /// Initial quality control
    FASTQC(read_ch)

    // Adapter trimming and post-trimming QC
    TRIM_GALORE(read_ch)

    // Alineamiento a un genoma de referencia
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Ejecutar el flujo de trabajo para probar que funciona

```bash
nextflow run rnaseq.nf -resume
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

Esta vez vemos una sola llamada a MULTIQC agregada después de las llamadas a procesos en caché:

Puede encontrar las salidas en `results/trimming` según lo especificado en el proceso `TRIM_GALORE` por la directiva `publishDir`.

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

---

## 3. Habilitar el procesamiento de datos de RNAseq de extremos pareados

Actualmente nuestro flujo de trabajo solo puede manejar datos de RNAseq de un solo extremo.
Es cada vez más común ver datos de RNAseq de extremos pareados, así que queremos poder manejar eso.

Hacer que el flujo de trabajo sea completamente agnóstico del tipo de datos requeriría usar características del lenguaje Nextflow un poco más avanzadas, así que no vamos a hacerlo aquí, pero podemos hacer una versión de procesamiento de extremos pareados para demostrar qué necesita ser adaptado.

### 3.1. Hacer una copia del flujo de trabajo llamada `rnaseq_pe.nf`

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Modificar el `input_csv` predeterminado para que apunte a los datos de extremos pareados

Proporcionamos un segundo archivo CSV que contiene IDs de muestra y rutas de archivos FASTQ pareados en el directorio `data/`

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Cambiemos el valor predeterminado de `input_csv` para que sea la ruta al archivo `paired-end.csv`.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Entrada principal
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_single-end"
}
```

### 3.3. Actualizar la fábrica de canales

Necesitamos decirle al operador `.map()` que obtenga ambas rutas de archivos FASTQ ahora.

Así que `row -> file(row.fastq_path)` se convierte en `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Crear canal de entrada a partir del contenido de un archivo CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Hacer una versión de extremos pareados del proceso FASTQC

Hagamos una copia del módulo para que podamos tener ambas versiones a mano.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Abra el nuevo archivo de módulo `fastqc_pe.nf` en el editor de código y haga los siguientes cambios de código:

- Cambie `fastqc $reads` a `fastqc ${reads}` en el bloque `script` (línea 17) para que la entrada `reads` se desempaquete, ya que ahora es una tupla de dos rutas en lugar de una sola ruta.
- Reemplace `${reads.simpleName}` con un comodín (`*`) para evitar tener que manejar los archivos de salida individualmente.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Técnicamente esto generaliza el proceso `FASTQC` de una manera que lo hace capaz de manejar datos de RNAseq de un solo extremo o de extremos pareados.

Finalmente, actualice la declaración de importación del módulo para usar la versión de extremos pareados del módulo.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Hacer una versión de extremos pareados del proceso TRIM_GALORE

Haga una copia del módulo para que podamos tener ambas versiones a mano.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Abra el nuevo archivo de módulo `trim_galore_pe.nf` en el editor de código y haga los siguientes cambios de código:

- Cambie la declaración de entrada de `path reads` a `tuple path(read1), path(read2)`
- Actualice el comando en el bloque `script`, reemplazando `$reads` con `--paired ${read1} ${read2}`
- Actualice las declaraciones de salida para reflejar los archivos agregados y las convenciones de nomenclatura diferentes, usando comodines para evitar tener que listar todo.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Finalmente, actualice la declaración de importación del módulo para usar la versión de extremos pareados del módulo.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Actualizar la llamada al proceso MULTIQC para esperar dos reportes de TRIM_GALORE

El proceso `TRIM_GALORE` ahora produce un canal de salida adicional, así que necesitamos alimentar eso a MultiQC.

Reemplace `TRIM_GALORE.out.fastqc_reports,` con `TRIM_GALORE.out.fastqc_reports_1,` más `TRIM_GALORE.out.fastqc_reports_2,`:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Comprehensive QC report generation
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Mientras estamos en MultiQC, actualicemos también el valor predeterminado del parámetro `report_id` de `"all_single-end"` a `"all_paired-end"`.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Entrada principal
    input_csv: Path = "data/paired-end.csv"

    // Reference genome archive
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // Report ID
    report_id: String = "all_paired-end"
}
```

### 3.7. Hacer una versión de extremos pareados del proceso HISAT2_ALIGN

Haga una copia del módulo para que podamos tener ambas versiones a mano.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Abra el nuevo archivo de módulo `hisat2_align_pe.nf` en el editor de código y haga los siguientes cambios de código:

- Cambie la declaración de entrada de `path reads` a `tuple path(read1), path(read2)`
- Actualice el comando en el bloque `script`, reemplazando `-U $reads` con `-1 ${read1} -2 ${read2}`
- Reemplace todas las instancias de `${reads.simpleName}` con `${read1.simpleName}` en el comando en el bloque `script` así como en las declaraciones de salida.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Finalmente, actualice la declaración de importación del módulo para usar la versión de extremos pareados del módulo.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Ejecutar el flujo de trabajo para probar que funciona

No usamos `-resume` ya que esto no usaría la caché, y hay el doble de datos para procesar que antes, pero aún así debería completarse en menos de un minuto.

```bash
nextflow run rnaseq_pe.nf
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

¡Y eso es todo! Ahora tenemos dos versiones ligeramente divergentes de nuestro flujo de trabajo, una para datos de lecturas de un solo extremo y una para datos de extremos pareados.
El siguiente paso lógico sería hacer que el flujo de trabajo acepte cualquier tipo de datos sobre la marcha, lo cual está fuera del alcance de este curso, pero podríamos abordar eso en un seguimiento.

---

### Conclusión

Usted sabe cómo adaptar un flujo de trabajo de una sola muestra para paralelizar el procesamiento de múltiples muestras, generar un reporte de control de calidad completo y adaptar el flujo de trabajo para usar datos de lecturas de extremos pareados si es necesario.

### ¿Qué sigue?

¡Felicitaciones, ha completado el mini-curso de Nextflow para RNAseq! ¡Celebre su éxito y tome un merecido descanso!

A continuación, le pedimos que complete una encuesta muy breve sobre su experiencia con este curso de entrenamiento, luego lo llevaremos a una página con enlaces a recursos de entrenamiento adicionales y enlaces útiles.
