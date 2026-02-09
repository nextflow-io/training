# Parte 2: Implementación de muestra única

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta parte del curso, vamos a escribir el workflow más simple posible que envuelva todos los comandos que ejecutamos en la Parte 1 para automatizar su ejecución, y nos enfocaremos en procesar una muestra a la vez.

Haremos esto en tres etapas:

1. Escribir un workflow de una sola etapa que ejecute el paso inicial de control de calidad
2. Agregar el recorte de adaptadores y el control de calidad posterior al recorte
3. Agregar el alineamiento al genoma de referencia

!!! warning "Requisito previo"

    Debes completar la Parte 1 del curso antes de comenzar esta lección.
    Específicamente, trabajar en las secciones 2.1-3 crea el archivo de índice del genoma (`data/genome_index.tar.gz`) requerido para el paso de alineamiento en esta lección.

---

## 1. Escribir un workflow de una sola etapa que ejecute el control de calidad inicial

Comencemos escribiendo un workflow simple que ejecute la herramienta FastQC en un archivo FASTQ que contiene lecturas de RNAseq de extremo único.

Te proporcionamos un archivo de workflow, `rnaseq.nf`, que describe las partes principales del workflow.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Declaraciones INCLUDE de módulos

/*
 * Parámetros del pipeline
 */

// Entrada principal

workflow {

    // Crear canal de entrada

    // Llamar procesos

}
```

Ten en cuenta que este código de workflow es correcto pero no es funcional; su propósito es solo servir como esqueleto que usarás para escribir el workflow real.

### 1.1. Crear un directorio para almacenar módulos

Crearemos módulos independientes para cada proceso para facilitar su gestión y reutilización, así que creemos un directorio para almacenarlos.

```bash
mkdir modules
```

### 1.2. Crear un módulo para el proceso de recolección de métricas de control de calidad

Creemos un archivo de módulo llamado `modules/fastqc.nf` para alojar el proceso `FASTQC`:

```bash
touch modules/fastqc.nf
```

Abre el archivo en el editor de código y copia el siguiente código en él:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

Deberías reconocer todas las partes de lo que aprendiste en la Parte 1 y Parte 2 de esta serie de capacitación; el único cambio notable es que esta vez estamos usando `mode: symlink` para la directiva `publishDir`, y estamos usando un parámetro para definir el `publishDir`.

!!! note

    Aunque los archivos de datos que estamos usando aquí son muy pequeños, en genómica pueden llegar a ser muy grandes. Para fines de demostración en el entorno de enseñanza, estamos usando el modo de publicación 'symlink' para evitar copias de archivos innecesarias. No deberías hacer esto en tus workflows finales, ya que perderás resultados cuando limpies tu directorio `work`.

### 1.3. Importar el módulo en el archivo de workflow

Agrega la declaración `include { FASTQC } from './modules/fastqc.nf'` al archivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declaraciones INCLUDE de módulos
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Agregar una declaración de entrada

Declara un parámetro de entrada con un valor predeterminado:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Entrada principal
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Crear un canal de entrada en el bloque workflow

Usa una factory de canal básica `.fromPath()` para crear el canal de entrada:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Crear canal de entrada desde una ruta de archivo
    read_ch = channel.fromPath(params.reads)

    // Llamar procesos

}
```

### 1.6. Llamar al proceso `FASTQC` en el canal de entrada

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Crear canal de entrada desde una ruta de archivo
    read_ch = channel.fromPath(params.reads)

    // Control de calidad inicial
    FASTQC(read_ch)

}
```

### 1.7. Ejecutar el workflow para probar que funciona

Podríamos usar el parámetro `--reads` para especificar una entrada desde la línea de comandos, pero durante el desarrollo podemos ser perezosos y simplemente usar el valor predeterminado de prueba que configuramos.

```bash
nextflow run rnaseq.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Esto debería ejecutarse muy rápidamente si trabajaste en la Parte 1 y ya has descargado el contenedor.
Si te la saltaste, Nextflow descargará el contenedor por ti; no tienes que hacer nada para que suceda, pero es posible que debas esperar hasta un minuto.

Puedes encontrar las salidas en `results/fastqc` como se especifica en el proceso `FASTQC` por la directiva `publishDir`.

```bash
ls results/fastqc
```

```console title="Salida"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Agregar el recorte de adaptadores y el control de calidad posterior al recorte

Vamos a usar el wrapper Trim_Galore, que incluye Cutadapt para el recorte en sí y FastQC para el control de calidad posterior al recorte.

### 2.1. Crear un módulo para el proceso de recorte y control de calidad

Creemos un archivo de módulo llamado `modules/trim_galore.nf` para alojar el proceso `TRIM_GALORE`:

```bash
touch modules/trim_galore.nf
```

Abre el archivo en el editor de código y copia el siguiente código en él:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Importar el módulo en el archivo de workflow

Agrega la declaración `include { TRIM_GALORE } from './modules/trim_galore.nf'` al archivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declaraciones INCLUDE de módulos
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Llamar al proceso en el canal de entrada

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Crear canal de entrada desde una ruta de archivo
    read_ch = channel.fromPath(params.reads)

    // Control de calidad inicial
    FASTQC(read_ch)

    // Recorte de adaptadores y control de calidad posterior al recorte
    TRIM_GALORE(read_ch)
}
```

### 2.4. Ejecutar el workflow para probar que funciona

```bash
nextflow run rnaseq.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

Esto también debería ejecutarse muy rápidamente, ya que estamos ejecutando en un archivo de entrada tan pequeño.

Puedes encontrar las salidas en `results/trimming` como se especifica en el proceso `TRIM_GALORE` por la directiva `publishDir`.

```bash
ls results/trimming
```

```console title="Salida"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Alinear las lecturas al genoma de referencia

Finalmente podemos ejecutar el paso de alineamiento del genoma usando Hisat2, que también emitirá métricas de control de calidad estilo FastQC.

### 3.1. Crear un módulo para el proceso HiSat2

Creemos un archivo de módulo llamado `modules/hisat2_align.nf` para alojar el proceso `HISAT2_ALIGN`:

```bash
touch modules/hisat2_align.nf
```

Abre el archivo en el editor de código y copia el siguiente código en él:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Importar el módulo en el archivo de workflow

Agrega la declaración `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` al archivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declaraciones INCLUDE de módulos
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Agregar una declaración de parámetro para proporcionar el índice del genoma

Declara un parámetro de entrada con un valor predeterminado:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Entrada principal
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Archivo del genoma de referencia
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Llamar al proceso `HISAT2_ALIGN` en las lecturas recortadas generadas por `TRIM_GALORE`

Las lecturas recortadas están en el canal `TRIM_GALORE.out.trimmed_reads` generado por el paso anterior.

Además, usamos `file (params.hisat2_index_zip)` para proporcionar a la herramienta Hisat2 el tarball comprimido del índice del genoma.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Crear canal de entrada desde una ruta de archivo
    read_ch = channel.fromPath(params.reads)

    // Control de calidad inicial
    FASTQC(read_ch)

    // Recorte de adaptadores y control de calidad posterior al recorte
    TRIM_GALORE(read_ch)

    // Alineamiento a un genoma de referencia
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Ejecutar el workflow para probar que funciona

```bash
nextflow run rnaseq.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Puedes encontrar las salidas en `results/align` como se especifica en el proceso `HISAT2_ALIGN` por la directiva `publishDir`.

```bash
ls results/align
```

```console title="Salida"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Esto completa el procesamiento básico que necesitamos aplicar a cada muestra.

_Agregaremos la agregación de reportes MultiQC en la Parte 2, después de que hayamos modificado el workflow para aceptar múltiples muestras a la vez._

---

### Conclusión

Sabes cómo envolver todos los pasos principales para procesar muestras de RNAseq de extremo único individualmente.

### ¿Qué sigue?

Aprende cómo modificar el workflow para procesar múltiples muestras en paralelo, agregar reportes de control de calidad en todos los pasos para todas las muestras, y habilitar la ejecución del workflow en datos de RNAseq de extremo pareado.
