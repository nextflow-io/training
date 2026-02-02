# Parte 1: Descripción general del método y pruebas manuales

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Existen múltiples métodos válidos para procesar y analizar datos de RNAseq en bulk.
Para este curso, seguimos el método descrito [aquí](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) por los Drs. Simon Andrews y Laura Biggins en el [Babraham Institute](https://www.babraham.ac.uk/).

Nuestro objetivo es desarrollar un flujo de trabajo que implemente los siguientes pasos de procesamiento: ejecutar control de calidad inicial en las lecturas de una muestra de RNAseq en bulk, recortar secuencias de adaptadores de las lecturas, alinear las lecturas a un genoma de referencia y producir un informe completo de control de calidad (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/preprocess.svg"
</figure>

- **FASTQC:** Realizar QC en los datos de lectura antes del recorte usando FastQC
- **TRIM_GALORE:** Recortar secuencias de adaptadores y realizar QC después del recorte usando Trim Galore (agrupa Cutadapt y FastQC)
- **HISAT2_ALIGN:** Alinear lecturas al genoma de referencia usando Hisat2
- **MULTIQC:** Generar un informe QC completo usando MultiQC

Sin embargo, antes de comenzar a escribir cualquier código de flujo de trabajo, vamos a probar los comandos manualmente con algunos datos de prueba.
Las herramientas que necesitamos no están instaladas en el entorno de GitHub Codespaces, por lo que las usaremos a través de contenedores (ver [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Asegúrese de estar en el directorio `nf4-science/rnaseq`. La última parte de la ruta que se muestra cuando escribe `pwd` debe ser `rnaseq`.

---

## 1. QC inicial y recorte de adaptadores

Vamos a descargar una imagen de contenedor que tiene instalados tanto `fastqc` como `trim_galore`, iniciarlo de forma interactiva y ejecutar los comandos de recorte y QC en uno de los archivos de datos de ejemplo.

### 1.1. Descargar el contenedor

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Esto le proporciona la siguiente salida de consola mientras el sistema descarga la imagen:

??? success "Salida del comando"

    ```console
    0.6.10--1bf8ca4e1967cd18: Pulling from library/trim-galore
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    32ec762be2d0: Pull complete
    d2cb90387285: Pull complete
    Digest: sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907
    Status: Downloaded newer image for community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
    ```

### 1.2. Iniciar el contenedor de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

Su prompt cambiará a algo como `(base) root@b645838b3314:/tmp#`, lo que indica que ahora está dentro del contenedor.

La parte `-v ./data:/data` del comando nos permitirá acceder al contenido del directorio `data/` desde dentro del contenedor.

```bash
ls /data/reads
```

??? success "Salida del comando"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gzO
    ```

### 1.3. Ejecutar el primer comando `fastqc`

Ejecutemos `fastqc` para recopilar métricas de control de calidad en los datos de lectura.

```bash
fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Salida del comando"

    ```console
    application/gzip
    Started analysis of ENCSR000COQ1_1.fastq.gz
    Approx 5% complete for ENCSR000COQ1_1.fastq.gz
    Approx 10% complete for ENCSR000COQ1_1.fastq.gz
    Approx 15% complete for ENCSR000COQ1_1.fastq.gz
    Approx 20% complete for ENCSR000COQ1_1.fastq.gz
    Approx 25% complete for ENCSR000COQ1_1.fastq.gz
    Approx 30% complete for ENCSR000COQ1_1.fastq.gz
    Approx 35% complete for ENCSR000COQ1_1.fastq.gz
    Approx 40% complete for ENCSR000COQ1_1.fastq.gz
    Approx 45% complete for ENCSR000COQ1_1.fastq.gz
    Approx 50% complete for ENCSR000COQ1_1.fastq.gz
    Approx 55% complete for ENCSR000COQ1_1.fastq.gz
    Approx 60% complete for ENCSR000COQ1_1.fastq.gz
    Approx 65% complete for ENCSR000COQ1_1.fastq.gz
    Approx 70% complete for ENCSR000COQ1_1.fastq.gz
    Approx 75% complete for ENCSR000COQ1_1.fastq.gz
    Approx 80% complete for ENCSR000COQ1_1.fastq.gz
    Approx 85% complete for ENCSR000COQ1_1.fastq.gz
    Approx 90% complete for ENCSR000COQ1_1.fastq.gz
    Approx 95% complete for ENCSR000COQ1_1.fastq.gz
    Analysis complete for ENCSR000COQ1_1.fastq.gz
    ```

Esto debería ejecutarse muy rápidamente.
Puede encontrar los archivos de salida en el mismo directorio que los datos originales:

```bash
ls /data/reads/ENCSR000COQ1_1_fastqc*
```

<!-- switch to tree -->

```console title="Salida"
/data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
```

### 1.4. Recortar secuencias de adaptadores con `trim_galore`

Ahora ejecutemos `trim_galore`, que agrupa Cutadapt y FastQC, para recortar las secuencias de adaptadores y recopilar métricas QC posteriores al recorte.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

La opción `--fastqc` hace que el comando ejecute automáticamente un paso de recopilación de QC después de que se complete el recorte.

_La salida es muy detallada, por lo que lo siguiente está abreviado._

??? success "Salida del comando"

    ```console
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    <...>

    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

Puede encontrar los archivos de salida en el directorio de trabajo:

```bash
ls ENCSR000COQ1_1*
```

```console title="Salida"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
```

### 1.5. Mover los archivos de salida al sistema de archivos fuera del contenedor

Cualquier cosa que permanezca dentro del contenedor será inaccesible para trabajo futuro, así que movamos estos archivos a un nuevo directorio.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

### 1.6. Salir del contenedor

```bash
exit
```

---

## 2. Alinear las lecturas al genoma de referencia

Vamos a descargar una imagen de contenedor que tiene instalado `hisat2`, iniciarlo de forma interactiva y ejecutar el comando de alineamiento para alinear los datos de RNAseq a un genoma de referencia.

### 2.1. Descargar el contenedor de `hisat2`

```bash
docker pull community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

??? success "Salida del comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e' locally
    5e49f68a37dc010e: Pulling from library/hisat2_samtools
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    e74ed5dd390b: Pull complete
    abfcf0185e51: Pull complete
    Digest: sha256:29d8e1a3172a2bdde7be813f7ebec22d331388194a7c0de872b4ccca4bed8f45
    Status: Downloaded newer image for community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
    ```

### 2.2. Iniciar el contenedor de `hisat2` de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

El comando es el mismo que antes, con el URI del contenedor correspondiente intercambiado.

### 2.3. Crear los archivos de índice del genoma para Hisat2

Hisat2 requiere que la referencia del genoma se proporcione en un formato muy específico, y no puede simplemente consumir el archivo FASTA `genome.fa` que proporcionamos, por lo que vamos a aprovechar esta oportunidad para crear los recursos relevantes.

```bash
hisat2-build /data/genome.fa genome_index
```

La salida es muy detallada, por lo que lo siguiente está abreviado:

<!-- TODO: switch to full output -->

??? success "Salida del comando"

    ```console
    Settings:
      Output files: "genome_index.*.ht2"
    <...>
    Total time for call to driver() for forward index: 00:00:16
    ```

Esto crea múltiples archivos de índice del genoma, que puede encontrar en el directorio de trabajo.

```bash
ls genome_index.*
```

```console title="Salida"
genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
```

Los usaremos en un momento, pero primero generemos un archivo tar comprimido con estos archivos de índice del genoma; los necesitaremos más adelante y generar estos no es típicamente algo que queramos hacer como parte de un flujo de trabajo.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

Esto almacena un archivo tar `genome_index.tar.gz` que contiene los archivos de índice del genoma en el directorio `data/` en nuestro sistema de archivos, lo cual será útil en la Parte 2 de este curso.

### 2.4. Ejecutar el comando `hisat2`

Ahora podemos ejecutar el comando de alineamiento, que realiza el paso de alineamiento con `hisat2` y luego dirige la salida a `samtools` para escribir la salida como un archivo BAM.

La entrada de datos de lectura es el archivo `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que generamos con `trim_galore` en el paso anterior.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Salida del comando"

    ```console
    HISAT2 summary stats:
            Total reads: 27816
                    Aligned 0 time: 1550 (5.57%)
                    Aligned 1 time: 25410 (91.35%)
                    Aligned >1 times: 856 (3.08%)
            Overall alignment rate: 94.43%
    ```

Esto se ejecuta casi instantáneamente porque es un archivo de prueba muy pequeño.
A escala real, esto podría tardar mucho más.

Una vez más, puede encontrar los archivos de salida en el directorio de trabajo:

```bash
ls ENCSR000COQ1_1*
```

```console title="Salida"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

### 2.5. Mover los archivos de salida al sistema de archivos fuera del contenedor

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

### 2.6. Salir del contenedor

```bash
exit
```

---

## 3. Generar un informe QC completo

Vamos a descargar una imagen de contenedor que tiene instalado `multiqc`, iniciarlo de forma interactiva y ejecutar un comando de generación de informes en los archivos de informe FastQC antes/después.

### 3.1. Descargar el contenedor de `multiqc`

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Salida del comando"

    ```console
    ad8f247edb55897c: Pulling from library/pip_multiqc
    dafa2b0c44d2: Already exists
    dec6b097362e: Already exists
    f88da01cff0b: Already exists
    4f4fb700ef54: Already exists
    92dc97a3ef36: Already exists
    403f74b0f85e: Already exists
    10b8c00c10a5: Already exists
    17dc7ea432cc: Already exists
    bb36d6c3110d: Already exists
    0ea1a16bbe82: Already exists
    030a47592a0a: Already exists
    3f229294c69a: Pull complete
    5a5ad47fd84c: Pull complete
    Digest: sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

### 3.2. Iniciar el contenedor de `multiqc` de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

### 3.3. Ejecutar el comando `multiqc`

```bash
multiqc /data/reads /data/trimmed /data/aligned -n ENCSR000COQ1_1_QC
```

??? success "Salida del comando"

    ```console

    /// MultiQC 🔍 v1.27.1

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 20/20
                hisat2 | Found 1 reports
              cutadapt | Found 1 reports
                fastqc | Found 1 reports
        write_results | Data        : ENCSR000COQ1_1_QC_data
        write_results | Report      : ENCSR000COQ1_1_QC.html
              multiqc | MultiQC complete
    ```

MultiQC puede buscar en directorios informes QC compatibles y agregará todo lo que encuentre.

Aquí vemos que la herramienta encontró los tres informes QC que generamos: el QC inicial que hicimos con `fastqc`, el informe posterior al recorte de `cutadapt` (hecho a través de `trim_galore`) y el QC posterior al alineamiento producido por `hisat2`.

Los archivos de salida están una vez más en el directorio de trabajo:

```bash
ls ENCSR000COQ1_1_QC*
```

```console title="Salida"
ENCSR000COQ1_1_QC.html

ENCSR000COQ1_1_QC_data:
cutadapt_filtered_reads_plot.txt                     fastqc_top_overrepresented_sequences_table.txt
cutadapt_trimmed_sequences_plot_3_Counts.txt         hisat2_se_plot.txt
cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc.log
fastqc-status-check-heatmap.txt                      multiqc_citations.txt
fastqc_adapter_content_plot.txt                      multiqc_cutadapt.txt
fastqc_per_base_n_content_plot.txt                   multiqc_data.json
fastqc_per_base_sequence_quality_plot.txt            multiqc_fastqc.txt
fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_general_stats.txt
fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_hisat2.txt
fastqc_per_sequence_quality_scores_plot.txt          multiqc_software_versions.txt
fastqc_sequence_counts_plot.txt                      multiqc_sources.txt
fastqc_sequence_duplication_levels_plot.txt
```

### 3.4. Mover los archivos de salida al sistema de archivos fuera del contenedor

```bash
mkdir /data/final_qc
mv ENCSR000COQ1_1_QC** /data/final_qc
```

### 3.5. Salir del contenedor

```bash
exit
```

---

### Conclusión

Ha probado todos los comandos individuales de forma interactiva en los contenedores correspondientes.

### ¿Qué sigue?

Aprenda a envolver esos mismos comandos en un flujo de trabajo de múltiples pasos que usa contenedores para ejecutar el trabajo.
