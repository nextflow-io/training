# Parte 1: Descripción general del método

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Existen múltiples métodos válidos para procesar y analizar datos de RNAseq en bulk.
Para este curso, seguimos el método descrito [aquí](https://www.bioinformatics.babraham.ac.uk/training/RNASeq_Course/Analysing%20RNA-Seq%20data%20Exercise.pdf) por los Drs. Simon Andrews y Laura Biggins en el [Babraham Institute](https://www.babraham.ac.uk/).

Nuestro objetivo es desarrollar un flujo de trabajo que implemente los siguientes pasos de procesamiento: ejecutar control de calidad inicial en las lecturas de una muestra de RNAseq en bulk, recortar secuencias de adaptadores de las lecturas, alinear las lecturas a un genoma de referencia y producir un informe completo de control de calidad (QC).

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

- **FASTQC:** Realizar QC en los datos de lectura antes del recorte usando FastQC
- **TRIM_GALORE:** Recortar secuencias de adaptadores y realizar QC después del recorte usando Trim Galore (agrupa Cutadapt y FastQC)
- **HISAT2_ALIGN:** Alinear lecturas al genoma de referencia usando Hisat2
- **MULTIQC:** Generar un informe QC completo usando MultiQC

### Métodos

Vamos a mostrarle cómo aplicar estos pasos de procesamiento en dos fases.
Primero comenzaremos con **procesamiento de muestra única** que ejecuta las herramientas de QC, recorte y alineamiento en una muestra.
Luego extenderemos a **procesamiento de múltiples muestras** que ejecuta las mismas herramientas en múltiples muestras y genera un informe de control de calidad agregado.

Antes de comenzar a escribir cualquier código de flujo de trabajo para cualquiera de los enfoques, vamos a probar los comandos manualmente con algunos datos de prueba.

### Conjunto de datos

Proporcionamos los siguientes datos y recursos relacionados:

- **Datos de RNAseq** (`reads/`): archivos FASTQ de seis muestras, reducidos a una región pequeña para mantener los tamaños de archivo bajos. Cada muestra tiene lecturas de extremos emparejados (dos archivos por muestra), aunque comenzamos trabajando solo con lecturas de extremo único.
- **Un genoma de referencia** (`genome.fa`): una región pequeña del cromosoma humano 20 (de hg19/b37).
- **Hojas de cálculo CSV** (`single-end.csv` y `paired-end.csv`): archivos que listan los IDs y rutas de los archivos de datos de ejemplo.

### Software

Las cuatro herramientas principales involucradas son [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) para recopilación de métricas de control de calidad, [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) para recorte de adaptadores (agrupa Cutadapt y FastQC para QC posterior al recorte), [HISAT2](http://daehwankimlab.github.io/hisat2/) para alineamiento con empalmes a un genoma de referencia, y [MultiQC](https://multiqc.info/) para generación de informes QC agregados.

Estas herramientas no están instaladas en el entorno de GitHub Codespaces, por lo que las usaremos a través de contenedores obtenidos mediante el servicio Seqera Containers (ver [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! tip "Consejo"

     Asegúrese de estar en el directorio `nf4-science/rnaseq`. La última parte de la ruta que se muestra cuando escribe `pwd` debe ser `rnaseq`.

---

## 1. Procesamiento de muestra única

En esta sección probamos los comandos que procesan una única muestra de RNAseq: control de calidad, recorte de adaptadores y alineamiento a un genoma de referencia.
Estos son los comandos que envolveremos en un flujo de trabajo de Nextflow en la Parte 2 de este curso.

1. Ejecutar QC inicial en un archivo FASTQ usando FastQC
2. Recortar secuencias de adaptadores y ejecutar QC posterior al recorte usando Trim Galore
3. Alinear las lecturas recortadas al genoma de referencia usando HISAT2

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-02.svg"
</figure>

Comenzamos probando estos comandos en solo una muestra.

### 1.1. QC y recorte de adaptadores

Primero, queremos ejecutar los comandos de QC y recorte en uno de los archivos de datos de ejemplo.

#### 1.1.1. Descargar el contenedor

Descarguemos una imagen de contenedor que tiene instalados tanto `fastqc` como `trim_galore`:

```bash
docker pull community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

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

Si no ha descargado esta imagen antes, puede tardar un minuto en completarse.
Una vez que termine, tendrá una copia local de la imagen del contenedor.

#### 1.1.2. Iniciar el contenedor de forma interactiva

Para ejecutar el contenedor de forma interactiva, use `docker run` con las opciones `-it`.
La opción `-v ./data:/data` monta nuestro directorio local `data/` para que podamos acceder a los archivos de entrada desde dentro del contenedor.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

??? success "Salida del comando"

    ```console
    (base) root@b645838b3314:/tmp#
    ```

Su prompt cambiará a algo como `(base) root@b645838b3314:/tmp#`, lo que indica que ahora está dentro del contenedor.

Verifique que puede ver los archivos de datos de secuencia bajo `/data/reads`:

```bash
ls /data/reads
```

??? abstract "Contenido del directorio"

    ```console
    ENCSR000COQ1_1.fastq.gz  ENCSR000COQ2_2.fastq.gz  ENCSR000COR2_1.fastq.gz  ENCSR000CPO1_2.fastq.gz
    ENCSR000COQ1_2.fastq.gz  ENCSR000COR1_1.fastq.gz  ENCSR000COR2_2.fastq.gz  ENCSR000CPO2_1.fastq.gz
    ENCSR000COQ2_1.fastq.gz  ENCSR000COR1_2.fastq.gz  ENCSR000CPO1_1.fastq.gz  ENCSR000CPO2_2.fastq.gz
    ```

Con eso, está listo para probar su primer comando.

#### 1.1.3. Ejecutar el comando FastQC

El método referenciado anteriormente nos proporciona la línea de comando para ejecutar QC en un solo archivo.
Solo necesitamos proporcionar el archivo de entrada; la herramienta generará automáticamente archivos de salida en el mismo directorio que los datos originales.

Ejecute el comando `fastqc` en un archivo de datos:

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

??? abstract "Contenido del directorio"

    ```console
    /data/reads/ENCSR000COQ1_1_fastqc.html  /data/reads/ENCSR000COQ1_1_fastqc.zip
    ```

Debería ver un informe HTML y un archivo ZIP que contiene las métricas QC.
Eso completa la prueba del primer paso.

#### 1.1.4. Recortar secuencias de adaptadores con Trim Galore

Ahora ejecutemos `trim_galore`, que agrupa Cutadapt y FastQC, para recortar las secuencias de adaptadores y recopilar métricas QC posteriores al recorte.
Como se señaló anteriormente, el software está incluido en el mismo contenedor, por lo que no se necesita ningún cambio allí.

El comando es sencillo; simplemente necesitamos agregar la opción `--fastqc` para ejecutar automáticamente un paso de recopilación de QC después de que se complete el recorte.

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ1_1.fastq.gz
```

??? success "Salida del comando"

    ```console hl_lines="54 55 56 58 59 60"
    Multicore support not enabled. Proceeding with single-core trimming.
    Path to Cutadapt set as: 'cutadapt' (default)
    Cutadapt seems to be working fine (tested command 'cutadapt --version')
    Cutadapt version: 4.9
    single-core operation.
    igzip command line interface 2.31.0
    igzip detected. Using igzip for decompressing

    No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)



    AUTO-DETECTING ADAPTER TYPE
    ===========================
    Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /data/reads/ENCSR000COQ1_1.fastq.gz <<)

    Found perfect matches for the following adapter sequences:
    Adapter type	Count	Sequence	Sequences analysed	Percentage
    Illumina	9	AGATCGGAAGAGC	27816	0.03
    smallRNA	0	TGGAATTCTCGG	27816	0.00
    Nextera	0	CTGTCTCTTATA	27816	0.00
    Using Illumina adapter for trimming (count: 9). Second best hit was smallRNA (count: 0)

    Writing report to 'ENCSR000COQ1_1.fastq.gz_trimming_report.txt'

    SUMMARISING RUN PARAMETERS
    ==========================
    Input filename: /data/reads/ENCSR000COQ1_1.fastq.gz
    Trimming mode: single-end
    Trim Galore version: 0.6.10
    Cutadapt version: 4.9
    Number of cores used for trimming: 1
    Quality Phred score cutoff: 20
    Quality encoding type selected: ASCII+33
    Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
    Maximum trimming error rate: 0.1 (default)
    Minimum required adapter overlap (stringency): 1 bp
    Minimum required sequence length before a sequence gets removed: 20 bp
    Running FastQC on the data once trimming has completed
    Output file(s) will be GZIP compressed

    Cutadapt seems to be fairly up-to-date (version 4.9). Setting -j 1
    Writing final adapter and quality trimmed output to ENCSR000COQ1_1_trimmed.fq.gz


      >>> Now performing quality (cutoff '-q 20') and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /data/reads/ENCSR000COQ1_1.fastq.gz <<<
    This is cutadapt 4.9 with Python 3.12.7
    Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /data/reads/ENCSR000COQ1_1.fastq.gz
    Processing single-end reads on 1 core ...
    Finished in 0.373 s (13.399 µs/read; 4.48 M reads/minute).

    === Summary ===

    Total reads processed:                  27,816
    Reads with adapters:                     9,173 (33.0%)
    Reads written (passing filters):        27,816 (100.0%)

    Total basepairs processed:     2,114,016 bp
    Quality-trimmed:                       0 bp (0.0%)
    Total written (filtered):      2,100,697 bp (99.4%)

    === Adapter 1 ===

    Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 9173 times

    Minimum overlap: 1
    No. of allowed errors:
    1-9 bp: 0; 10-13 bp: 1

    Bases preceding removed adapters:
      A: 27.4%
      C: 37.4%
      G: 20.9%
      T: 14.3%
      none/other: 0.0%

    Overview of removed sequences
    length	count	expect	max.err	error counts
    1	6229	6954.0	0	6229
    2	2221	1738.5	0	2221
    3	581	434.6	0	581
    4	88	108.7	0	88
    5	33	27.2	0	33
    6	2	6.8	0	2
    7	1	1.7	0	1
    9	1	0.1	0	1
    10	2	0.0	1	2
    12	1	0.0	1	0 1
    14	4	0.0	1	3 1
    16	1	0.0	1	1
    19	1	0.0	1	1
    22	1	0.0	1	1
    29	4	0.0	1	0 4
    33	3	0.0	1	3

    RUN STATISTICS FOR INPUT FILE: /data/reads/ENCSR000COQ1_1.fastq.gz
    =============================================
    27816 sequences processed in total
    Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)


      >>> Now running FastQC on the data <<<

    application/gzip
    Started analysis of ENCSR000COQ1_1_trimmed.fq.gz
    Approx 5% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 10% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 15% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 20% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 25% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 30% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 35% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 40% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 45% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 50% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 55% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 60% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 65% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 70% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 75% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 80% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 85% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 90% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Approx 95% complete for ENCSR000COQ1_1_trimmed.fq.gz
    Analysis complete for ENCSR000COQ1_1_trimmed.fq.gz
    ```

La salida es muy detallada, por lo que hemos resaltado las líneas más relevantes en el ejemplo anterior.
Puede encontrar los archivos de salida en el directorio de trabajo:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contenido del directorio"

    ```console
    ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.html
    ENCSR000COQ1_1_trimmed.fq.gz                 ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Esto incluye las lecturas recortadas, el informe de recorte y los archivos QC posteriores al recorte.

#### 1.1.5. Mover los archivos de salida

Cualquier cosa que permanezca dentro del contenedor será inaccesible para trabajo futuro, por lo que necesitamos mover estos archivos a un directorio en el sistema de archivos montado.

```bash
mkdir /data/trimmed
mv ENCSR000COQ1_1* /data/trimmed
```

??? abstract "Contenido del directorio"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    └── ENCSR000COQ1_1_trimmed_fastqc.zip
    ```

Los archivos ahora son accesibles en su sistema de archivos normal.

#### 1.1.6. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad; eso completa la prueba de los dos primeros pasos.

### 1.2. Alinear las lecturas al genoma de referencia

A continuación, queremos ejecutar el comando de alineamiento para alinear las lecturas de RNAseq recortadas a un genoma de referencia.

#### 1.2.1. Descargar el contenedor

Descarguemos una imagen de contenedor que tiene instalados `hisat2` y `samtools`:

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

Notará que algunas capas muestran `Already exists` porque se comparten con la imagen del contenedor Trim Galore que descargamos anteriormente.
Como resultado, esta descarga debería ser más rápida que la primera.

#### 1.2.2. Iniciar el contenedor de forma interactiva

Inicie el contenedor de forma interactiva, usando el mismo enfoque que antes con el URI del contenedor relevante intercambiado.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Su prompt cambiará nuevamente para indicar que está dentro del contenedor.

#### 1.2.3. Crear los archivos de índice del genoma

HISAT2 requiere que la referencia del genoma se proporcione en un formato muy específico, y no puede simplemente consumir el archivo FASTA `genome.fa` que proporcionamos, por lo que vamos a aprovechar esta oportunidad para crear los recursos relevantes.

```bash
hisat2-build /data/genome.fa genome_index
```

??? success "Salida del comando"

    ```console hl_lines="1 2 218"
    Settings:
      Output files: "genome_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /data/genome.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 6542727 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6542727 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:01
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 7; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 4.98493e+06 (target: 6542726)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Getting block 1 of 7
      Reserving size (6542727) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 3540952 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 3540953 for bucket 1
    Getting block 2 of 7
      Reserving size (6542727) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 6195795 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6195796 for bucket 2
    Getting block 3 of 7
      Reserving size (6542727) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 6199288 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 6199289 for bucket 3
    Getting block 4 of 7
      Reserving size (6542727) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 6454986 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 6454987 for bucket 4
    Getting block 5 of 7
      Reserving size (6542727) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 3493181 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3493182 for bucket 5
    Getting block 6 of 7
      Reserving size (6542727) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 5875908 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 5875909 for bucket 6
    Getting block 7 of 7
      Reserving size (6542727) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 3134429 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 3134430 for bucket 7
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 9094775
    fchr[G]: 17470759
    fchr[T]: 25839994
    fchr[$]: 34894545
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 15826295 bytes to primary GFM file: genome_index.1.ht2
    Wrote 8723644 bytes to secondary GFM file: genome_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 15353415 bytes to primary GFM file: genome_index.5.ht2
    Wrote 8883598 bytes to secondary GFM file: genome_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 34894545
        gbwtLen: 34894546
        nodes: 34894546
        sz: 8723637
        gbwtSz: 8723637
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2180910
        offsSz: 8723640
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 181743
        numLines: 181743
        gbwtTotLen: 11631552
        gbwtTotSz: 11631552
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:12
    ```

La salida es muy detallada, por lo que hemos resaltado algunas líneas relevantes en el ejemplo anterior.

Esto crea múltiples archivos de índice del genoma, que puede encontrar en el directorio de trabajo.

```bash
ls genome_index.*
```

??? abstract "Contenido del directorio"

    ```console
    genome_index.1.ht2  genome_index.3.ht2  genome_index.5.ht2  genome_index.7.ht2
    genome_index.2.ht2  genome_index.4.ht2  genome_index.6.ht2  genome_index.8.ht2
    ```

Necesitaremos estos archivos más adelante, y generar estos no es típicamente algo que queramos hacer como parte de un flujo de trabajo, por lo que vamos a generar un archivo tar comprimido que contenga los archivos de índice del genoma que podamos pasar fácilmente según sea necesario.

```bash
tar -czvf /data/genome_index.tar.gz genome_index.*
```

??? success "Salida del comando"

    ```console
    genome_index.1.ht2
    genome_index.2.ht2
    genome_index.3.ht2
    genome_index.4.ht2
    genome_index.5.ht2
    genome_index.6.ht2
    genome_index.7.ht2
    genome_index.8.ht2
    ```

Moveremos el archivo tar resultante `genome_index.tar.gz` que contiene los archivos de índice del genoma al directorio `data/` en nuestro sistema de archivos en unos minutos.
Eso será útil en la Parte 2 de este curso.

#### 1.2.4. Ejecutar el comando de alineamiento

Ahora podemos ejecutar el comando de alineamiento, que realiza el paso de alineamiento con `hisat2` y luego dirige la salida a `samtools` para escribir la salida como un archivo BAM.

La entrada de datos de lectura es el archivo `/data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz` que generamos con `trim_galore` en el paso anterior.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ1_1_trimmed.bam
```

??? success "Salida del comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 27816
    		Aligned 0 time: 1550 (5.57%)
    		Aligned 1 time: 25410 (91.35%)
    		Aligned >1 times: 856 (3.08%)
    	Overall alignment rate: 94.43%
    ```

Esto se ejecuta casi instantáneamente porque es un archivo de prueba muy pequeño.
A escala real, esto podría tardar mucho más.

Una vez más puede encontrar los archivos de salida en el directorio de trabajo:

```bash
ls ENCSR000COQ1_1*
```

??? abstract "Contenido del directorio"

    ```console title="Salida"
    ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
    ```

El alineamiento produjo un archivo BAM y un archivo de registro con estadísticas de alineamiento.

#### 1.2.5. Mover los archivos de salida

Como antes, mueva los archivos de salida a un directorio en el sistema de archivos montado para que permanezcan accesibles después de salir del contenedor.

```bash
mkdir /data/aligned
mv ENCSR000COQ1_1* /data/aligned
```

Con eso hecho, tenemos todo lo que necesitamos.

#### 1.2.6. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad.
Eso concluye la ejecución de prueba del procesamiento de muestra única.

!!! example "¡Escríbalo como un flujo de trabajo!"

    Siéntase libre de pasar a la [Parte 2](./02_single-sample.md) de inmediato si desea comenzar a implementar este análisis como un flujo de trabajo de Nextflow.
    Solo necesitará regresar para completar la segunda ronda de pruebas antes de pasar a la Parte 3.

---

## 2. Agregación de QC de múltiples muestras

Los comandos que acabamos de probar procesan una muestra a la vez.
En la práctica, típicamente necesitamos procesar muchas muestras y luego agregar los resultados de QC en todas ellas para evaluar la calidad del conjunto de datos general.

[MultiQC](https://multiqc.info/) es una herramienta que busca en directorios informes QC de muchas herramientas bioinformáticas comunes y los agrega en un único informe HTML completo.
Puede reconocer salida de FastQC, Cutadapt (a través de Trim Galore) y HISAT2, entre muchas otras.

Aquí procesamos dos muestras adicionales a través de las mismas herramientas por muestra, luego usamos MultiQC para agregar informes QC en las tres muestras.
Estos son los comandos que envolveremos en un flujo de trabajo de Nextflow en la Parte 3 de este curso.

1. Ejecutar QC y recorte en muestras adicionales usando Trim Galore
2. Ejecutar alineamiento en muestras adicionales usando HISAT2
3. Agregar todos los informes QC en un informe completo usando MultiQC

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-method-03.svg"
</figure>

### 2.1. QC y recorte de muestras adicionales

Los comandos de QC y recorte por muestra son idénticos a lo que ejecutamos en la sección 1.1.
Ya descargamos la imagen del contenedor, por lo que podemos iniciarlo directamente.

#### 2.1.1. Iniciar el contenedor

Ya descargamos esta imagen de contenedor en la sección 1.1, por lo que podemos iniciarla directamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18
```

Su prompt cambia para indicar que está dentro del contenedor.

#### 2.1.2. Ejecutar QC y recorte en muestras adicionales

Ejecute FastQC y Trim Galore en dos muestras más, una tras otra.

```bash
fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

```bash
trim_galore --fastqc /data/reads/ENCSR000COQ2_1.fastq.gz
trim_galore --fastqc /data/reads/ENCSR000COR1_1.fastq.gz
```

Una vez que esto se complete, debería tener archivos de salida de Trim Galore para ambas muestras en el directorio de trabajo.

#### 2.1.3. Mover los archivos de salida

Mueva los archivos de salida de Trim Galore al mismo directorio que usamos en la sección 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/trimmed
```

??? abstract "Contenido del directorio"

    ```console
    /data/trimmed
    ├── ENCSR000COQ1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ1_1_trimmed.fq.gz
    ├── ENCSR000COQ1_1_trimmed_fastqc.html
    ├── ENCSR000COQ1_1_trimmed_fastqc.zip
    ├── ENCSR000COQ2_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COQ2_1_trimmed.fq.gz
    ├── ENCSR000COQ2_1_trimmed_fastqc.html
    ├── ENCSR000COQ2_1_trimmed_fastqc.zip
    ├── ENCSR000COR1_1.fastq.gz_trimming_report.txt
    ├── ENCSR000COR1_1_trimmed.fq.gz
    ├── ENCSR000COR1_1_trimmed_fastqc.html
    └── ENCSR000COR1_1_trimmed_fastqc.zip
    ```

Los archivos ahora son accesibles en su sistema de archivos normal.

#### 2.1.4. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad.

### 2.2. Alinear muestras adicionales

Los comandos de alineamiento son idénticos a lo que ejecutamos en la sección 1.2.
Necesitamos extraer el índice del genoma del archivo tar que guardamos anteriormente, ya que los archivos de índice originales se crearon dentro de un contenedor que ya no existe.

#### 2.2.1. Iniciar el contenedor

Ya descargamos esta imagen de contenedor en la sección 1.2, por lo que podemos iniciarla directamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e
```

Su prompt cambia para indicar que está dentro del contenedor.

#### 2.2.2. Extraer el índice del genoma

Extraiga los archivos de índice del genoma del archivo tar que guardamos en el sistema de archivos montado:

```bash
tar -xzf /data/genome_index.tar.gz
```

Esto restaura los archivos `genome_index.*` en el directorio de trabajo.

#### 2.2.3. Ejecutar alineamiento en muestras adicionales

Ejecute el alineamiento HISAT2 en las dos muestras recién recortadas, una tras otra.

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COQ2_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COQ2_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COQ2_1_trimmed.bam
```

??? success "Salida del comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 18736
    		Aligned 0 time: 1531 (8.17%)
    		Aligned 1 time: 16726 (89.27%)
    		Aligned >1 times: 479 (2.56%)
    	Overall alignment rate: 91.83%
    ```

```bash
hisat2 -x genome_index -U /data/trimmed/ENCSR000COR1_1_trimmed.fq.gz \
    --new-summary --summary-file ENCSR000COR1_1_trimmed.hisat2.log | \
    samtools view -bS -o ENCSR000COR1_1_trimmed.bam
```

??? success "Salida del comando"

    ```console hl_lines="6"
    HISAT2 summary stats:
    	Total reads: 38056
    		Aligned 0 time: 2311 (6.07%)
    		Aligned 1 time: 33289 (87.47%)
    		Aligned >1 times: 2456 (6.45%)
    	Overall alignment rate: 93.93%
    ```

Una vez que esto se complete, debería tener archivos BAM y de registro para ambas muestras en el directorio de trabajo.

#### 2.2.4. Mover los archivos de salida

Mueva los archivos de salida del alineamiento al mismo directorio que usamos en la sección 1.

```bash
mv ENCSR000COQ2_1* ENCSR000COR1_1* /data/aligned
```

??? abstract "Contenido del directorio"

    ```console
    /data/aligned
    ├── ENCSR000COQ1_1_trimmed.bam
    ├── ENCSR000COQ1_1_trimmed.hisat2.log
    ├── ENCSR000COQ2_1_trimmed.bam
    ├── ENCSR000COQ2_1_trimmed.hisat2.log
    ├── ENCSR000COR1_1_trimmed.bam
    └── ENCSR000COR1_1_trimmed.hisat2.log
    ```

Los archivos ahora son accesibles en su sistema de archivos normal.

#### 2.2.5. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad.

### 2.3. Generar un informe QC completo

Ahora que tenemos salida de QC, recorte y alineamiento para tres muestras, podemos usar MultiQC para agregarlas en un único informe.
MultiQC busca en directorios informes QC compatibles y agrega todo lo que encuentra.

#### 2.3.1. Descargar el contenedor

Descarguemos una imagen de contenedor que tiene instalado `multiqc`:

```bash
docker pull community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

??? success "Salida del comando"

    ```console
    a3c26f6199d64b7c: Pulling from library/pip_multiqc
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
    2ed162b168e8: Pull complete
    ca06fe148f21: Pull complete
    Digest: sha256:af0e9de56896805aa2a065f7650362956f4213d99e95314f6fec472c6a3bf091
    Status: Downloaded newer image for community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
    ```

Notará que algunas capas muestran `Already exists` porque se comparten con las imágenes de contenedor que descargamos anteriormente.
Como resultado, esta descarga debería ser más rápida que las anteriores.

#### 2.3.2. Iniciar el contenedor de forma interactiva

Inicie el contenedor de forma interactiva con el directorio de datos montado, como antes.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c
```

Su prompt cambiará para indicar que está dentro del contenedor.

#### 2.3.3. Ejecutar el comando MultiQC

Ejecute `multiqc`, apuntándolo a los directorios donde almacenamos archivos de salida relacionados con QC para las tres muestras.
La opción `-n` establece el nombre del informe de salida.

```bash
multiqc /data/reads /data/trimmed /data/aligned -n all_samples_QC
```

??? success "Salida del comando"

    ```console hl_lines="8 9 10 11 12"

    /// MultiQC 🔍 v1.32

          file_search | Search path: /data/reads
          file_search | Search path: /data/trimmed
          file_search | Search path: /data/aligned
            searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 36/36
               hisat2 | Found 3 reports
             cutadapt | Found 3 reports
               fastqc | Found 3 reports
        write_results | Data        : all_samples_QC_data
        write_results | Report      : all_samples_QC.html
              multiqc | MultiQC complete
    ```

Aquí vemos que la herramienta encontró informes QC para las tres muestras: el QC inicial de `fastqc`, los informes posteriores al recorte de `cutadapt` (a través de `trim_galore`) y los resúmenes de alineamiento producidos por `hisat2`.

Los archivos de salida están en el directorio de trabajo:

```bash
ls all_samples_QC*
```

??? abstract "Contenido del directorio"

    ```console
    all_samples_QC.html

    all_samples_QC_data:
    cutadapt_filtered_reads_plot.txt                     multiqc.log
    cutadapt_trimmed_sequences_plot_3_Counts.txt         multiqc.parquet
    cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt        multiqc_citations.txt
    fastqc-status-check-heatmap.txt                      multiqc_cutadapt.txt
    fastqc_adapter_content_plot.txt                      multiqc_data.json
    fastqc_overrepresented_sequences_plot.txt            multiqc_fastqc.txt
    fastqc_per_base_n_content_plot.txt                   multiqc_general_stats.txt
    fastqc_per_base_sequence_quality_plot.txt            multiqc_hisat2.txt
    fastqc_per_sequence_gc_content_plot_Counts.txt       multiqc_software_versions.txt
    fastqc_per_sequence_gc_content_plot_Percentages.txt  multiqc_sources.txt
    fastqc_per_sequence_quality_scores_plot.txt
    fastqc_sequence_counts_plot.txt
    fastqc_sequence_duplication_levels_plot.txt
    fastqc_top_overrepresented_sequences_table.txt
    hisat2_se_plot.txt
    llms-full.txt
    ```

La salida principal es el informe `all_samples_QC.html`, acompañado de un directorio de datos que contiene las métricas subyacentes.

#### 2.3.4. Mover los archivos de salida

Mueva el informe y su directorio de datos al sistema de archivos montado.

```bash
mkdir /data/multiqc
mv all_samples_QC* /data/multiqc
```

Los archivos ahora son accesibles en su sistema de archivos normal.

#### 2.3.5. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad.
Eso concluye la prueba de todos los comandos de procesamiento de RNAseq.

---

### Conclusión

Sabe cómo ejecutar los comandos FastQC, Trim Galore, HISAT2 y MultiQC en sus respectivos contenedores, incluyendo cómo procesar múltiples muestras y agregar informes QC.

### ¿Qué sigue?

Tome un descanso, luego diríjase a la [Parte 2](./02_single-sample.md) para aprender cómo envolver esos mismos comandos en flujos de trabajo que usan contenedores para ejecutar el trabajo.
