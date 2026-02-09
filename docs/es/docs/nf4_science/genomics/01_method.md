# Parte 1: Descripción del método y pruebas manuales

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

El llamado de variantes es un método de análisis genómico que tiene como objetivo identificar variaciones en una secuencia genómica en relación con un genoma de referencia.
Aquí vamos a usar herramientas y métodos diseñados para llamar variantes germinales cortas, _es decir_ SNPs e indels, en datos de secuenciación de genoma completo.

![Pipeline GATK](img/gatk-pipeline.png)

Un pipeline completo de llamado de variantes típicamente involucra muchos pasos, incluyendo el mapeo a la referencia (a veces referido como alineamiento del genoma) y el filtrado y priorización de variantes.
Por simplicidad, en este curso nos vamos a enfocar solo en la parte de llamado de variantes.

### Métodos

Vamos a mostrarte dos formas de aplicar el llamado de variantes a muestras de secuenciación de genoma completo para identificar SNPs e indels germinales.
Primero comenzaremos con un **enfoque simple por muestra** que llama variantes independientemente de cada muestra.
Luego te mostraremos un **enfoque de llamado conjunto** más sofisticado que analiza múltiples muestras juntas, produciendo resultados más precisos e informativos.

Antes de escribir cualquier código de workflow para cualquiera de los dos enfoques, vamos a probar los comandos manualmente en algunos datos de prueba.

### Conjunto de datos

Proporcionamos los siguientes datos y recursos relacionados:

- **Un genoma de referencia** que consiste en una pequeña región del cromosoma humano 20 (de hg19/b37) y sus archivos accesorios (índice y diccionario de secuencia).
- **Tres muestras de secuenciación de genoma completo** correspondientes a un trío familiar (madre, padre e hijo), que han sido reducidas a una pequeña porción de datos en el cromosoma 20 para mantener los tamaños de archivo pequeños.
  Estos son datos de secuenciación Illumina de lecturas cortas que ya han sido mapeados al genoma de referencia, proporcionados en formato [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, una versión comprimida de SAM, Sequence Alignment Map).
- **Una lista de intervalos genómicos**, es decir, coordenadas en el genoma donde nuestras muestras tienen datos adecuados para llamar variantes, proporcionada en formato BED.

### Software

Las dos herramientas principales involucradas son [Samtools](https://www.htslib.org/), un conjunto de herramientas ampliamente utilizado para manipular archivos de alineamiento de secuencias, y [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit), un conjunto de herramientas para el descubrimiento de variantes desarrollado en el Broad Institute.

Estas herramientas no están instaladas en el entorno de GitHub Codespaces, así que las usaremos a través de contenedores (ver [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

     Asegúrate de estar en el directorio `nf4-science/genomics` para que la última parte de la ruta mostrada cuando escribes `pwd` sea `genomics`.

---

## 1. Llamado de variantes por muestra

El llamado de variantes por muestra procesa cada muestra independientemente: el llamador de variantes examina los datos de secuenciación de una muestra a la vez e identifica posiciones donde la muestra difiere de la referencia.

En esta sección probamos los dos comandos que conforman el enfoque de llamado de variantes por muestra: indexar un archivo BAM con Samtools y llamar variantes con GATK HaplotypeCaller.
Estos son los comandos que envolveremos en un workflow de Nextflow en la Parte 2 de este curso.

1. Generar un archivo índice para un archivo de entrada BAM usando [Samtools](https://www.htslib.org/)
2. Ejecutar GATK HaplotypeCaller en el archivo BAM indexado para generar llamados de variantes por muestra en VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Comenzamos probando los dos comandos en solo una muestra.

### 1.1. Indexar un archivo de entrada BAM con Samtools

Los archivos índice son una característica común de los formatos de archivo bioinformáticos; contienen información sobre la estructura del archivo principal que permite a herramientas como GATK acceder a un subconjunto de los datos sin tener que leer todo el archivo.
Esto es importante debido a lo grandes que pueden llegar a ser estos archivos.

Los archivos BAM a menudo se proporcionan sin un índice, por lo que el primer paso en muchos workflows de análisis es generar uno usando `samtools index`.

Vamos a descargar un contenedor de Samtools, iniciarlo de forma interactiva y ejecutar el comando `samtools index` en uno de los archivos BAM.

#### 1.1.1. Descargar el contenedor de Samtools

Ejecuta el comando `docker pull` para descargar la imagen del contenedor de Samtools:

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

??? success "Salida del comando"

    ```console
    1.20--b5dfbd93de237464: Pulling from library/samtools
    6360b3717211: Pull complete
    2ec3f7ad9b3c: Pull complete
    7716ca300600: Pull complete
    4f4fb700ef54: Pull complete
    8c61d418774c: Pull complete
    03dae77ff45c: Pull complete
    aab7f787139d: Pull complete
    4f4fb700ef54: Pull complete
    837d55536720: Pull complete
    897362c12ca7: Pull complete
    3893cbe24e91: Pull complete
    d1b61e94977b: Pull complete
    c72ff66fb90f: Pull complete
    0e0388f29b6d: Pull complete
    Digest: sha256:bbfc45b4f228975bde86cba95e303dd94ecf2fdacea5bfb2e2f34b0d7b141e41
    Status: Downloaded newer image for community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
    ```

Si no has descargado esta imagen antes, puede tardar un minuto en completarse.
Una vez que termine, tienes una copia local de la imagen del contenedor.

#### 1.1.2. Iniciar el contenedor de Samtools de forma interactiva

Para ejecutar el contenedor de forma interactiva, usa `docker run` con las opciones `-it`.
La opción `-v ./data:/data` monta el directorio local `data` dentro del contenedor para que las herramientas puedan acceder a los archivos de entrada.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Tu prompt cambia a algo como `(base) root@a1b2c3d4e5f6:/tmp#`, indicando que ahora estás dentro del contenedor.
Los archivos de datos son accesibles bajo `/data`.

#### 1.1.3. Ejecutar el comando de indexación

La [documentación de Samtools](https://www.htslib.org/doc/samtools-index.html) nos da la línea de comando para ejecutar para indexar un archivo BAM.

Solo necesitamos proporcionar el archivo de entrada; la herramienta generará automáticamente un nombre para la salida agregando `.bai` al nombre del archivo de entrada.

```bash
samtools index /data/bam/reads_mother.bam
```

??? abstract "Contenidos del directorio"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

Ahora deberías ver un archivo llamado `reads_mother.bam.bai` en el mismo directorio que el archivo de entrada BAM original.

#### 1.1.4. Salir del contenedor de Samtools

Para salir del contenedor, escribe `exit`.

```bash
exit
```

Tu prompt ahora debería volver a lo que era antes de iniciar el contenedor.

### 1.2. Llamar variantes con GATK HaplotypeCaller

Vamos a descargar un contenedor de GATK, iniciarlo de forma interactiva y ejecutar el comando `gatk HaplotypeCaller` en el archivo BAM que acabamos de indexar.

#### 1.2.1. Descargar el contenedor de GATK

Ejecuta el comando `docker pull` para descargar la imagen del contenedor de GATK:

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

??? success "Salida del comando"

    Algunas capas muestran `Already exists` porque son compartidas con la imagen del contenedor de Samtools que descargamos anteriormente.

    ```console
    4.5.0.0--730ee8817e436867: Pulling from library/gatk4
    6360b3717211: Already exists
    2ec3f7ad9b3c: Already exists
    7716ca300600: Already exists
    4f4fb700ef54: Already exists
    8c61d418774c: Already exists
    03dae77ff45c: Already exists
    aab7f787139d: Already exists
    4f4fb700ef54: Already exists
    837d55536720: Already exists
    897362c12ca7: Already exists
    3893cbe24e91: Already exists
    d1b61e94977b: Already exists
    e5c558f54708: Pull complete
    087cce32d294: Pull complete
    Digest: sha256:e33413b9100f834fcc62fd5bc9edc1e881e820aafa606e09301eac2303d8724b
    Status: Downloaded newer image for community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
    ```

Esto debería ser más rápido que la primera descarga porque las dos imágenes de contenedor comparten la mayoría de sus capas.

#### 1.2.2. Iniciar el contenedor de GATK de forma interactiva

Inicia el contenedor de GATK de forma interactiva con el directorio de datos montado, tal como lo hicimos para Samtools.

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Tu prompt cambia para indicar que ahora estás dentro del contenedor de GATK.

#### 1.2.3. Ejecutar el comando de llamado de variantes

La [documentación de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) nos da la línea de comando para ejecutar para realizar el llamado de variantes en un archivo BAM.

Necesitamos proporcionar el archivo de entrada BAM (`-I`) así como el genoma de referencia (`-R`), un nombre para el archivo de salida (`-O`) y una lista de intervalos genómicos a analizar (`-L`).

Sin embargo, no necesitamos especificar la ruta al archivo índice; la herramienta lo buscará automáticamente en el mismo directorio, basándose en la convención establecida de nomenclatura y ubicación.
Lo mismo aplica para los archivos accesorios del genoma de referencia (archivos de índice y diccionario de secuencia, `*.fai` y `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

??? success "Salida del comando"

    La herramienta produce una salida de registro detallada. Las líneas resaltadas confirman la finalización exitosa.

    ```console hl_lines="37 51 56 57"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.vcf -L /data/ref/intervals.bed
    00:27:50.687 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:27:50.854 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.858 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:27:50.858 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:27:50.858 INFO  HaplotypeCaller - Executing as root@a1fe8ff42d07 on Linux v6.10.14-linuxkit amd64
    00:27:50.858 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:27:50.859 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:27:50 AM GMT
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.859 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:27:50.861 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:27:50.861 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:27:50.861 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:27:50.862 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:27:50.863 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:27:50.864 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:27:50.864 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:27:50.864 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:27:50.864 INFO  HaplotypeCaller - Requester pays: disabled
    00:27:50.865 INFO  HaplotypeCaller - Initializing engine
    00:27:50.991 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:27:51.016 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:27:51.029 INFO  HaplotypeCaller - Done initializing engine
    00:27:51.040 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:27:51.042 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:27:51.042 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:27:51.046 INFO  HaplotypeCallerEngine - Disabling physical phasing, which is supported only for reference-model confidence output
    00:27:51.063 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:27:51.085 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:27:51.086 INFO  IntelPairHmm - Available threads: 10
    00:27:51.086 INFO  IntelPairHmm - Requested threads: 4
    00:27:51.086 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:27:51.128 INFO  ProgressMeter - Starting traversal
    00:27:51.136 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:27:51.882 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:27:52.969 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:27:52.971 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35           1145.7
    00:27:52.971 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:27:52.976 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003346916
    00:27:52.976 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.045731709
    00:27:52.977 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:27:52.981 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:27:52 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.04 minutes.
    Runtime.totalMemory()=203423744
    ```

El archivo de salida `reads_mother.vcf` se crea dentro de tu directorio de trabajo en el contenedor, por lo que no lo verás en el explorador de archivos de VS Code a menos que cambies la ruta del archivo de salida.
Sin embargo, es un archivo de prueba pequeño, por lo que puedes usar `cat` para abrirlo y ver el contenido.
Si te desplazas hasta el inicio del archivo, encontrarás un encabezado compuesto de muchas líneas de metadatos, seguido de una lista de llamados de variantes, uno por línea.

??? abstract "Contenidos del archivo"

    ```console title="reads_mother.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Cada línea describe una posible variante identificada en los datos de secuenciación de la muestra. Para orientación sobre cómo interpretar el formato VCF, consulta [este artículo útil](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

El archivo de salida VCF está acompañado de un archivo índice llamado `reads_mother.vcf.idx` que fue creado automáticamente por GATK.
Tiene la misma función que el archivo índice BAM, permitir que las herramientas busquen y recuperen subconjuntos de datos sin cargar todo el archivo.

#### 1.2.4. Salir del contenedor de GATK

Para salir del contenedor, escribe `exit`.

```bash
exit
```

Tu prompt debería volver a la normalidad.
Esto concluye la prueba de llamado de variantes por muestra.

---

## 2. Llamado conjunto en una cohorte

El enfoque de llamado de variantes que acabamos de usar genera llamados de variantes por muestra.
Eso está bien para observar variantes de cada muestra de forma aislada, pero produce información limitada.
A menudo es más interesante observar cómo difieren los llamados de variantes entre múltiples muestras.
GATK ofrece un método alternativo llamado llamado conjunto de variantes para este propósito.

El llamado conjunto de variantes implica generar un tipo especial de salida de variantes llamada GVCF (Genomic VCF) para cada muestra, luego combinar los datos GVCF de todas las muestras y ejecutar un análisis estadístico de 'genotipado conjunto'.

![Análisis conjunto](img/joint-calling.png)

Lo que es especial sobre el GVCF de una muestra es que contiene registros que resumen estadísticas de datos de secuencia sobre todas las posiciones en el área objetivo del genoma, no solo las posiciones donde el programa encontró evidencia de variación.
Esto es crítico para el cálculo de genotipado conjunto ([lectura adicional](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

El GVCF es producido por GATK HaplotypeCaller, la misma herramienta que acabamos de probar, con un parámetro adicional (`-ERC GVCF`).
La combinación de los GVCFs se realiza con GATK GenomicsDBImport, que combina los llamados por muestra en un almacén de datos (análogo a una base de datos).
El análisis de 'genotipado conjunto' propiamente dicho se realiza entonces con GATK GenotypeGVCFs.

Aquí probamos los comandos necesarios para generar GVCFs y ejecutar el genotipado conjunto.
Estos son los comandos que envolveremos en un workflow de Nextflow en la Parte 3 de este curso.

1. Generar un archivo índice para cada archivo de entrada BAM usando Samtools
2. Ejecutar GATK HaplotypeCaller en cada archivo de entrada BAM para generar un GVCF de llamados de variantes genómicas por muestra
3. Recolectar todos los GVCFs y combinarlos en un almacén de datos GenomicsDB
4. Ejecutar el genotipado conjunto en el almacén de datos GVCF combinado para producir un VCF a nivel de cohorte

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Ahora necesitamos probar todos estos comandos, comenzando con la indexación de los tres archivos BAM.

### 2.1. Indexar archivos BAM para las tres muestras

En la primera sección anterior, solo indexamos un archivo BAM.
Ahora necesitamos indexar las tres muestras para que GATK HaplotypeCaller pueda procesarlas.

#### 2.1.1. Iniciar el contenedor de Samtools de forma interactiva

Ya descargamos la imagen del contenedor de Samtools, así que podemos iniciarlo directamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

Tu prompt cambia para indicar que estás dentro del contenedor, con el directorio de datos montado como antes.

#### 2.1.2. Ejecutar el comando de indexación en las tres muestras

Ejecuta el comando de indexación en cada uno de los tres archivos BAM:

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

??? abstract "Contenidos del directorio"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Esto debería producir los archivos índice en el mismo directorio que los archivos BAM correspondientes.

#### 2.1.3. Salir del contenedor de Samtools

Para salir del contenedor, escribe `exit`.

```bash
exit
```

Tu prompt debería volver a la normalidad.

### 2.2. Generar GVCFs para las tres muestras

Para ejecutar el paso de genotipado conjunto, necesitamos GVCFs para las tres muestras.

#### 2.2.1. Iniciar el contenedor de GATK de forma interactiva

Ya descargamos la imagen del contenedor de GATK anteriormente, así que podemos iniciarlo directamente:

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

Tu prompt cambia para indicar que estás dentro del contenedor de GATK.

#### 2.2.2. Ejecutar el comando de llamado de variantes con la opción GVCF

Para producir un VCF genómico (GVCF), agregamos la opción `-ERC GVCF` al comando base, que activa el modo GVCF de HaplotypeCaller.

También cambiamos la extensión del archivo de salida de `.vcf` a `.g.vcf`.
Esto técnicamente no es un requisito, pero es una convención fuertemente recomendada.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

??? success "Salida del comando"

    ```console hl_lines="39 53 58 59"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R /data/ref/ref.fasta -I /data/bam/reads_mother.bam -O reads_mother.g.vcf -L /data/ref/intervals.bed -ERC GVCF
    00:28:03.593 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:03.765 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.768 INFO  HaplotypeCaller - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:03.768 INFO  HaplotypeCaller - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:03.768 INFO  HaplotypeCaller - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:03.768 INFO  HaplotypeCaller - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:03.769 INFO  HaplotypeCaller - Start Date/Time: February 8, 2026 at 12:28:03 AM GMT
    00:28:03.769 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.770 INFO  HaplotypeCaller - ------------------------------------------------------------
    00:28:03.772 INFO  HaplotypeCaller - HTSJDK Version: 4.1.0
    00:28:03.773 INFO  HaplotypeCaller - Picard Version: 3.1.1
    00:28:03.773 INFO  HaplotypeCaller - Built for Spark Version: 3.5.0
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:03.773 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:03.774 INFO  HaplotypeCaller - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:03.774 INFO  HaplotypeCaller - Deflater: IntelDeflater
    00:28:03.774 INFO  HaplotypeCaller - Inflater: IntelInflater
    00:28:03.775 INFO  HaplotypeCaller - GCS max retries/reopens: 20
    00:28:03.775 INFO  HaplotypeCaller - Requester pays: disabled
    00:28:03.776 INFO  HaplotypeCaller - Initializing engine
    00:28:03.896 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:03.919 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:03.934 INFO  HaplotypeCaller - Done initializing engine
    00:28:03.935 INFO  HaplotypeCallerEngine - Tool is in reference confidence mode and the annotation, the following changes will be made to any specified annotations: 'StrandBiasBySample' will be enabled. 'ChromosomeCounts', 'FisherStrand', 'StrandOddsRatio' and 'QualByDepth' annotations have been disabled
    00:28:03.943 INFO  NativeLibraryLoader - Loading libgkl_utils.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_utils.so
    00:28:03.945 INFO  NativeLibraryLoader - Loading libgkl_smithwaterman.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_smithwaterman.so
    00:28:03.946 INFO  SmithWatermanAligner - Using AVX accelerated SmithWaterman implementation
    00:28:03.955 INFO  HaplotypeCallerEngine - Standard Emitting and Calling confidence set to -0.0 for reference-model confidence output
    00:28:03.956 INFO  HaplotypeCallerEngine - All sites annotated with PLs forced to true for reference-model confidence output
    00:28:03.972 INFO  NativeLibraryLoader - Loading libgkl_pairhmm_omp.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_pairhmm_omp.so
    00:28:03.993 INFO  IntelPairHmm - Flush-to-zero (FTZ) is enabled when running PairHMM
    00:28:03.994 INFO  IntelPairHmm - Available threads: 10
    00:28:03.994 INFO  IntelPairHmm - Requested threads: 4
    00:28:03.994 INFO  PairHMM - Using the OpenMP multi-threaded AVX-accelerated native PairHMM implementation
    00:28:04.044 INFO  ProgressMeter - Starting traversal
    00:28:04.070 INFO  ProgressMeter -        Current Locus  Elapsed Minutes     Regions Processed   Regions/Minute
    00:28:04.874 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    00:28:06.535 INFO  HaplotypeCaller - 7 read(s) filtered by: MappingQualityReadFilter
    0 read(s) filtered by: MappingQualityAvailableReadFilter
    0 read(s) filtered by: MappedReadFilter
    0 read(s) filtered by: NotSecondaryAlignmentReadFilter
    0 read(s) filtered by: NotDuplicateReadFilter
    0 read(s) filtered by: PassesVendorQualityCheckReadFilter
    0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
    0 read(s) filtered by: GoodCigarReadFilter
    0 read(s) filtered by: WellformedReadFilter
    7 total reads filtered out of 1867 reads processed
    00:28:06.537 INFO  ProgressMeter - 20_10037292_10066351:13499              0.0                    35            851.6
    00:28:06.538 INFO  ProgressMeter - Traversal complete. Processed 35 total regions in 0.0 minutes.
    00:28:06.543 INFO  VectorLoglessPairHMM - Time spent in setup for JNI call : 0.003648749
    00:28:06.544 INFO  PairHMM - Total compute time in PairHMM computeLogLikelihoods() : 0.031498916
    00:28:06.544 INFO  SmithWatermanAligner - Total compute time in native Smith-Waterman : 0.02 sec
    00:28:06.547 INFO  HaplotypeCaller - Shutting down engine
    [February 8, 2026 at 12:28:06 AM GMT] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 0.05 minutes.
    Runtime.totalMemory()=281018368
    ```

Esto crea el archivo de salida GVCF `reads_mother.g.vcf` en el directorio de trabajo actual en el contenedor.

Si usas `cat` para ver el contenido, verás que es mucho más largo que el VCF equivalente que generamos en la sección 1. Ni siquiera puedes desplazarte hasta el inicio del archivo, y la mayoría de las líneas se ven bastante diferentes de lo que vimos en el VCF.

??? abstract "Contenidos del archivo"

    ```console title="reads_mother.g.vcf" linenums="1674"
    20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
    20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
    20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
    ```

Estas representan regiones no variantes donde el llamador de variantes no encontró evidencia de variación, por lo que capturó algunas estadísticas que describen su nivel de confianza en la ausencia de variación.
Esto hace posible distinguir entre dos cifras de casos muy diferentes: (1) hay datos de buena calidad que muestran que la muestra es homocigota-referencia, y (2) no hay suficientes datos buenos disponibles para hacer una determinación de cualquier manera.

En un GVCF, típicamente hay muchas de estas líneas no variantes, con un número menor de registros de variantes dispersos entre ellas.
Intenta ejecutar `head -176` en el GVCF para cargar solo las primeras 176 líneas del archivo para encontrar un llamado de variante real.

??? abstract "Contenidos del archivo"

    ```console title="reads_mother.g.vcf" linenums="174"
    20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
    20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
    20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
    ```

La segunda línea muestra el primer registro de variante en el archivo, que corresponde a la primera variante en el archivo VCF que observamos anteriormente.

Al igual que el VCF original, el archivo de salida GVCF también está acompañado de un archivo índice, llamado `reads_mother.g.vcf.idx`.

#### 2.2.3. Repetir el proceso en las otras dos muestras

Genera GVCFs para las dos muestras restantes ejecutando los comandos a continuación, uno tras otro.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Una vez que esto se complete, deberías tener tres archivos que terminan en `.g.vcf` en tu directorio actual (uno por muestra) y sus respectivos archivos índice que terminan en `.g.vcf.idx`.

¡Pero no salgas del contenedor!
Vamos a usar el mismo contenedor en el siguiente paso.

### 2.3. Ejecutar el genotipado conjunto

Ahora que tenemos todos los GVCFs, podemos probar el enfoque de genotipado conjunto para generar llamados de variantes para una cohorte de muestras.
Es un método de dos pasos que consiste en combinar los datos de todos los GVCFs en un almacén de datos, luego ejecutar el análisis de genotipado conjunto propiamente dicho para generar el VCF final de variantes llamadas conjuntamente.

#### 2.3.1. Combinar todos los GVCFs por muestra

Este primer paso usa otra herramienta de GATK, llamada GenomicsDBImport, para combinar los datos de todos los GVCFs en un almacén de datos GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

??? success "Salida del comando"

    ```console hl_lines="33 36 37 39 40"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenomicsDBImport -V reads_mother.g.vcf -V reads_father.g.vcf -V reads_son.g.vcf -L /data/ref/intervals.bed --genomicsdb-workspace-path family_trio_gdb
    00:28:20.772 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:20.914 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.917 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:20.917 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:20.917 INFO  GenomicsDBImport - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:20.917 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:20.918 INFO  GenomicsDBImport - Start Date/Time: February 8, 2026 at 12:28:20 AM GMT
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.918 INFO  GenomicsDBImport - ------------------------------------------------------------
    00:28:20.920 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    00:28:20.921 INFO  GenomicsDBImport - Picard Version: 3.1.1
    00:28:20.921 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    00:28:20.922 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:20.923 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:20.923 INFO  GenomicsDBImport - Deflater: IntelDeflater
    00:28:20.924 INFO  GenomicsDBImport - Inflater: IntelInflater
    00:28:20.924 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    00:28:20.924 INFO  GenomicsDBImport - Requester pays: disabled
    00:28:20.925 INFO  GenomicsDBImport - Initializing engine
    00:28:21.144 INFO  FeatureManager - Using codec BEDCodec to read file file:///data/ref/intervals.bed
    00:28:21.152 INFO  IntervalArgumentCollection - Processing 6369 bp from intervals
    00:28:21.157 INFO  GenomicsDBImport - Done initializing engine
    00:28:21.287 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:21.290 INFO  GenomicsDBImport - Vid Map JSON file will be written to /tmp/family_trio_gdb/vidmap.json
    00:28:21.290 INFO  GenomicsDBImport - Callset Map JSON file will be written to /tmp/family_trio_gdb/callset.json
    00:28:21.291 INFO  GenomicsDBImport - Complete VCF Header will be written to /tmp/family_trio_gdb/vcfheader.vcf
    00:28:21.291 INFO  GenomicsDBImport - Importing to workspace - /tmp/family_trio_gdb
    00:28:21.453 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.757 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.859 INFO  GenomicsDBImport - Importing batch 1 with 3 samples
    00:28:21.979 INFO  GenomicsDBImport - Done importing batch 1/1
    00:28:21.988 INFO  GenomicsDBImport - Import completed!
    00:28:21.988 INFO  GenomicsDBImport - Shutting down engine
    [February 8, 2026 at 12:28:21 AM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.02 minutes.
    Runtime.totalMemory()=305135616
    ```

La salida de este paso es efectivamente un directorio que contiene un conjunto de directorios anidados adicionales que contienen los datos de variantes combinados en forma de múltiples archivos diferentes.
Puedes explorarlo pero rápidamente verás que este formato de almacén de datos no está diseñado para ser leído directamente por humanos.

!!! note "Nota"

    GATK incluye herramientas que hacen posible inspeccionar y extraer datos de llamados de variantes del almacén de datos según sea necesario.

#### 2.3.2. Ejecutar el análisis de genotipado conjunto propiamente dicho

Este segundo paso usa otra herramienta de GATK, llamada GenotypeGVCFs, para recalcular las estadísticas de variantes y los genotipos individuales a la luz de los datos disponibles en todas las muestras de la cohorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

??? success "Salida del comando"

    ```console hl_lines="30 35 37 38"
    Using GATK jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R /data/ref/ref.fasta -V gendb://family_trio_gdb -O family_trio.vcf
    00:28:24.625 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/conda/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    00:28:24.798 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.801 INFO  GenotypeGVCFs - The Genome Analysis Toolkit (GATK) v4.5.0.0
    00:28:24.801 INFO  GenotypeGVCFs - For support and documentation go to https://software.broadinstitute.org/gatk/
    00:28:24.801 INFO  GenotypeGVCFs - Executing as root@8515e5a0598e on Linux v6.10.14-linuxkit amd64
    00:28:24.801 INFO  GenotypeGVCFs - Java runtime: OpenJDK 64-Bit Server VM v17.0.11-internal+0-adhoc..src
    00:28:24.802 INFO  GenotypeGVCFs - Start Date/Time: February 8, 2026 at 12:28:24 AM GMT
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.802 INFO  GenotypeGVCFs - ------------------------------------------------------------
    00:28:24.804 INFO  GenotypeGVCFs - HTSJDK Version: 4.1.0
    00:28:24.804 INFO  GenotypeGVCFs - Picard Version: 3.1.1
    00:28:24.804 INFO  GenotypeGVCFs - Built for Spark Version: 3.5.0
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    00:28:24.805 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    00:28:24.806 INFO  GenotypeGVCFs - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    00:28:24.806 INFO  GenotypeGVCFs - Deflater: IntelDeflater
    00:28:24.806 INFO  GenotypeGVCFs - Inflater: IntelInflater
    00:28:24.807 INFO  GenotypeGVCFs - GCS max retries/reopens: 20
    00:28:24.807 INFO  GenotypeGVCFs - Requester pays: disabled
    00:28:24.808 INFO  GenotypeGVCFs - Initializing engine
    00:28:25.023 INFO  GenomicsDBLibLoader - GenomicsDB native library version : 1.5.1-84e800e
    00:28:25.081 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field InbreedingCoeff  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAC  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.082 INFO  NativeGenomicsDB - pid=162 tid=163 No valid combination operation found for INFO field MLEAF  - the field will NOT be part of INFO fields in the generated VCF records
    00:28:25.109 INFO  GenotypeGVCFs - Done initializing engine
    00:28:25.184 INFO  ProgressMeter - Starting traversal
    00:28:25.187 INFO  ProgressMeter -        Current Locus  Elapsed Minutes    Variants Processed  Variants/Minute
    00:28:25.446 WARN  InbreedingCoeff - InbreedingCoeff will not be calculated at position 20_10037292_10066351:3480 and possibly subsequent; at least 10 samples must have called genotypes
    GENOMICSDB_TIMER,GenomicsDB iterator next() timer,Wall-clock time(s),0.15034835899999904,Cpu time(s),0.1355218420000006
    00:28:26.189 INFO  ProgressMeter - 20_10037292_10066351:13953              0.0                  3390         202994.0
    00:28:26.190 INFO  ProgressMeter - Traversal complete. Processed 3390 total variants in 0.0 minutes.
    00:28:26.194 INFO  GenotypeGVCFs - Shutting down engine
    [February 8, 2026 at 12:28:26 AM GMT] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.03 minutes.
    Runtime.totalMemory()=296747008
    ```

Esto crea el archivo de salida VCF `family_trio.vcf` en el directorio de trabajo actual en el contenedor.
Es otro archivo razonablemente pequeño, así que puedes usar `cat` en este archivo para ver su contenido, y desplazarte hacia arriba para encontrar las primeras líneas de variantes.

??? abstract "Contenidos del archivo"

    ```console title="family_trio.vcf" linenums="40"
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
    20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
    20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
    20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
    ```

Esto se ve similar al VCF que generamos anteriormente, excepto que esta vez tenemos información a nivel de genotipo para las tres muestras.
Las últimas tres columnas en el archivo son los bloques de genotipo para las muestras, listados en orden alfabético.

Si observamos los genotipos llamados para nuestro trío familiar de prueba para la primera variante, vemos que el padre es heterocigoto-variante (`0/1`), y la madre y el hijo son ambos homocigotos-variante (`1/1`).

¡Esa es en última instancia la información que estamos buscando extraer del conjunto de datos!

#### 2.3.3. Salir del contenedor de GATK

Para salir del contenedor, escribe `exit`.

```bash
exit
```

Tu prompt debería volver a la normalidad.
Esto concluye las pruebas manuales de los comandos de llamado de variantes.

---

### Conclusión

Sabes cómo probar los comandos de indexación de Samtools y llamado de variantes de GATK en sus respectivos contenedores, incluyendo cómo generar GVCFs y ejecutar el genotipado conjunto en múltiples muestras.

### ¿Qué sigue?

Aprende cómo envolver esos mismos comandos en workflows que usan contenedores para ejecutar el trabajo.
