# Parte 1: Llamado de variantes por muestra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la primera parte de este curso, le mostramos cómo construir un pipeline simple de llamado de variantes que aplica el llamado de variantes de GATK a muestras de secuenciación individuales.

### Descripción general del método

El llamado de variantes es un método de análisis genómico que tiene como objetivo identificar variaciones en una secuencia genómica en relación con un genoma de referencia.
Aquí vamos a utilizar herramientas y métodos diseñados para llamar variantes cortas, _es decir_, SNPs e indels.

![Pipeline GATK](img/gatk-pipeline.png)

Un pipeline completo de llamado de variantes típicamente involucra muchos pasos, incluyendo el mapeo a la referencia (a veces denominado alineamiento del genoma) y el filtrado y priorización de variantes.
Para simplificar, en esta parte del curso nos vamos a enfocar solo en la parte del llamado de variantes.

### Conjunto de datos

Proporcionamos los siguientes datos y recursos relacionados:

- **Un genoma de referencia** que consiste en una pequeña región del cromosoma 20 humano (de hg19/b37) y sus archivos accesorios (índice y diccionario de secuencias).
- **Tres muestras de secuenciación de genoma completo** correspondientes a un trío familiar (madre, padre e hijo), que han sido reducidas a un pequeño fragmento de datos en el cromosoma 20 para mantener los tamaños de archivo pequeños.
  Estos son datos de secuenciación Illumina de lecturas cortas que ya han sido mapeados al genoma de referencia, proporcionados en formato [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, una versión comprimida de SAM, Sequence Alignment Map).
- **Una lista de intervalos genómicos**, es decir, coordenadas en el genoma donde nuestras muestras tienen datos adecuados para llamar variantes, proporcionados en formato BED.

### Workflow

En esta parte del curso, vamos a desarrollar un workflow que hace lo siguiente:

1. Generar un archivo de índice para cada archivo BAM de entrada usando [Samtools](https://www.htslib.org/)
2. Ejecutar GATK HaplotypeCaller en cada archivo BAM de entrada para generar llamados de variantes por muestra en VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note

    Los archivos de índice son una característica común de los formatos de archivo bioinformáticos; contienen información sobre la estructura del archivo principal que permite a herramientas como GATK acceder a un subconjunto de los datos sin tener que leer todo el archivo.
    Esto es importante debido a lo grandes que pueden llegar a ser estos archivos.

---

## 0. Calentamiento: Probar los comandos de Samtools y GATK de forma interactiva

Primero queremos probar los comandos manualmente antes de intentar envolverlos en un workflow.
Las herramientas que necesitamos (Samtools y GATK) no están instaladas en el entorno de GitHub Codespaces, así que las usaremos a través de contenedores (ver [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

     Asegúrese de estar en el directorio `nf4-science/genomics` para que la última parte de la ruta mostrada cuando escribe `pwd` sea `genomics`.

### 0.1. Indexar un archivo BAM de entrada con Samtools

Vamos a descargar un contenedor de Samtools, iniciarlo de forma interactiva y ejecutar el comando `samtools index` en uno de los archivos BAM.

#### 0.1.1. Descargar el contenedor de Samtools

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

#### 0.1.2. Iniciar el contenedor de Samtools de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

#### 0.1.3. Ejecutar el comando de indexación

La [documentación de Samtools](https://www.htslib.org/doc/samtools-index.html) nos da la línea de comando para ejecutar para indexar un archivo BAM.

Solo necesitamos proporcionar el archivo de entrada; la herramienta generará automáticamente un nombre para la salida agregando `.bai` al nombre del archivo de entrada.

```bash
samtools index /data/bam/reads_mother.bam
```

Esto debería completarse de inmediato, y ahora debería ver un archivo llamado `reads_mother.bam.bai` en el mismo directorio que el archivo BAM de entrada original.

??? abstract "Contenido del directorio"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Salir del contenedor de Samtools

```bash
exit
```

### 0.2. Llamar variantes con GATK HaplotypeCaller

Vamos a descargar un contenedor de GATK, iniciarlo de forma interactiva y ejecutar el comando `gatk HaplotypeCaller` en el archivo BAM que acabamos de indexar.

#### 0.2.1. Descargar el contenedor de GATK

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

#### 0.2.2. Iniciar el contenedor de GATK de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

#### 0.2.3. Ejecutar el comando de llamado de variantes

La [documentación de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) nos da la línea de comando para ejecutar para realizar el llamado de variantes en un archivo BAM.

Necesitamos proporcionar el archivo BAM de entrada (`-I`) así como el genoma de referencia (`-R`), un nombre para el archivo de salida (`-O`) y una lista de intervalos genómicos a analizar (`-L`).

Sin embargo, no necesitamos especificar la ruta al archivo de índice; la herramienta lo buscará automáticamente en el mismo directorio, basándose en la convención establecida de nomenclatura y co-ubicación.
Lo mismo se aplica a los archivos accesorios del genoma de referencia (archivos de índice y diccionario de secuencias, `*.fai` y `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

El archivo de salida `reads_mother.vcf` se crea dentro de su directorio de trabajo en el contenedor, por lo que no lo verá en el explorador de archivos de VS Code a menos que cambie la ruta del archivo de salida.
Sin embargo, es un archivo de prueba pequeño, así que puede usar `cat` para abrirlo y ver el contenido.
Si se desplaza hasta el inicio del archivo, encontrará un encabezado compuesto de muchas líneas de metadatos, seguido de una lista de llamados de variantes, uno por línea.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Cada línea describe una posible variante identificada en los datos de secuenciación de la muestra. Para orientación sobre cómo interpretar el formato VCF, vea [este artículo útil](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

El archivo VCF de salida está acompañado de un archivo de índice llamado `reads_mother.vcf.idx` que fue creado automáticamente por GATK.
Tiene la misma función que el archivo de índice BAM, permitir a las herramientas buscar y recuperar subconjuntos de datos sin cargar el archivo completo.

#### 0.2.4. Salir del contenedor de GATK

```bash
exit
```

### Conclusión

Sabe cómo probar los comandos de indexación de Samtools y llamado de variantes de GATK en sus respectivos contenedores.

### ¿Qué sigue?

Aprenda cómo envolver esos mismos comandos en un workflow de dos pasos que usa contenedores para ejecutar el trabajo.

---

## 1. Escribir un workflow de una etapa que ejecuta Samtools index en un archivo BAM

Le proporcionamos un archivo de workflow, `genomics-1.nf`, que describe las partes principales del workflow.
No es funcional; su propósito es solo servir como un esqueleto que usará para escribir el workflow real.

### 1.1. Definir el proceso de indexación

Comencemos escribiendo un proceso, que llamaremos `SAMTOOLS_INDEX`, describiendo la operación de indexación.

```groovy title="genomics-1.nf" linenums="9"
/*
 * Generar archivo de índice BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    path "${input_bam}.bai"

    script:
    """
    samtools index '$input_bam'
    """
}
```

Debería reconocer todas las partes de lo que aprendió en la Parte 1 y Parte 2 de esta serie de entrenamiento.

Este proceso va a requerir que pasemos una ruta de archivo a través de la entrada `input_bam`, así que configuremos eso a continuación.

### 1.2. Agregar una declaración de parámetro de entrada

En la parte superior del archivo, bajo la sección `Pipeline parameters`, declaramos un parámetro CLI llamado `reads_bam` y le damos un valor predeterminado.
De esa manera, podemos ser perezosos y no especificar la entrada cuando escribimos el comando para lanzar el pipeline (con fines de desarrollo).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Parámetros del pipeline
 */
params {
    // Entrada principal
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Ahora tenemos un proceso listo, así como un parámetro para darle una entrada para ejecutar, así que conectemos esas cosas juntas.

!!! note

    `${projectDir}` es una variable integrada de Nextflow que apunta al directorio donde se encuentra el script de workflow actual de Nextflow (`genomics-1.nf`).

    Esto facilita hacer referencia a archivos, directorios de datos y otros recursos incluidos en el repositorio del workflow sin codificar rutas absolutas.

### 1.3. Agregar bloque workflow para ejecutar SAMTOOLS_INDEX

En el bloque `workflow`, necesitamos configurar un **canal** para alimentar la entrada al proceso `SAMTOOLS_INDEX`; luego podemos llamar al proceso mismo para que se ejecute en el contenido de ese canal.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Crear canal de entrada (archivo único vía parámetro CLI)
    reads_ch = channel.fromPath(params.reads_bam)

    // Crear archivo de índice para el archivo BAM de entrada
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

El bloque workflow tiene dos secciones:

- `main:` contiene las operaciones de canal y llamadas a procesos
- `publish:` declara qué salidas deben publicarse, asignándolas a destinos con nombre

Notará que estamos usando la misma channel factory `.fromPath` que usamos en [Hello Channels](../../hello_nextflow/02_hello_channels.md).
De hecho, estamos haciendo algo muy similar.
La diferencia es que le estamos diciendo a Nextflow que solo cargue la ruta del archivo en sí en el canal como un elemento de entrada, en lugar de leer su contenido.

### 1.4. Agregar un bloque output para definir dónde se publican los resultados

Después del bloque workflow, agregamos un bloque `output` que especifica dónde publicar las salidas del workflow.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

Cada destino nombrado de la sección `publish:` (como `bam_index`) obtiene su propio bloque donde puede configurar la ruta de salida relativa al directorio de salida base.

!!! note

    Aunque los archivos de datos que estamos usando aquí son muy pequeños, en genómica pueden llegar a ser muy grandes.
    Por defecto, Nextflow crea enlaces simbólicos a los archivos de salida en el directorio de publicación, lo que evita copias de archivos innecesarias.
    Puede cambiar este comportamiento usando la opción `mode` (por ejemplo, `mode 'copy'`) para crear copias reales en su lugar.
    Tenga en cuenta que los enlaces simbólicos se romperán cuando limpie su directorio `work`, por lo que para workflows de producción es posible que desee usar `mode 'copy'`.

### 1.5. Configurar el directorio de salida

El directorio de salida base se establece mediante la opción de configuración `outputDir`. Agréguela a `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Ejecutar el workflow para verificar que el paso de indexación funciona

¡Ejecutemos el workflow! Como recordatorio, no necesitamos especificar una entrada en la línea de comando porque establecimos un valor predeterminado para la entrada cuando declaramos el parámetro de entrada.

```bash
nextflow run genomics-1.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Puede verificar que el archivo de índice se ha generado correctamente mirando en el directorio de trabajo o en el directorio de resultados.

??? abstract "Contenido del directorio de trabajo"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Contenido del directorio de resultados"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

¡Ahí está!

### Conclusión

Sabe cómo envolver una herramienta genómica en un workflow Nextflow de un solo paso y hacer que se ejecute usando un contenedor.

### ¿Qué sigue?

Agregar un segundo paso que consuma la salida del primero.

---

## 2. Agregar un segundo proceso para ejecutar GATK HaplotypeCaller en el archivo BAM indexado

Ahora que tenemos un índice para nuestro archivo de entrada, podemos pasar a configurar el paso de llamado de variantes, que es la parte interesante del workflow.

### 2.1. Definir el proceso de llamado de variantes

Escribamos un proceso, que llamaremos `GATK_HAPLOTYPECALLER`, describiendo la operación de llamado de variantes.

```groovy title="genomics-1.nf" linenums="44"
/*
 * Llamar variantes con GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path input_bam
    path input_bam_index
    path ref_fasta
    path ref_index
    path ref_dict
    path interval_list

    output:
    path "${input_bam}.vcf"     , emit: vcf
    path "${input_bam}.vcf.idx" , emit: idx

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
}
```

Notará que hemos introducido una nueva sintaxis aquí (`emit:`) para nombrar de forma única cada uno de nuestros canales de salida, y las razones para esto quedarán claras pronto.

Este comando toma bastantes más entradas, porque GATK necesita más información para realizar el análisis en comparación con un simple trabajo de indexación.
Pero notará que hay incluso más entradas definidas en el bloque de entradas de las que se enumeran en el comando de GATK. ¿Por qué es eso?

!!! note

    GATK sabe buscar el archivo de índice BAM y los archivos accesorios del genoma de referencia porque está al tanto de las convenciones que rodean esos archivos.
    Sin embargo, Nextflow está diseñado para ser agnóstico del dominio y no sabe nada sobre los requisitos de formato de archivo bioinformático.

Necesitamos decirle a Nextflow explícitamente que tiene que preparar esos archivos en el directorio de trabajo en tiempo de ejecución; de lo contrario no lo hará, y GATK (correctamente) lanzará un error sobre los archivos de índice que faltan.

De manera similar, tenemos que listar el archivo de índice del VCF de salida (el archivo `"${input_bam}.vcf.idx"`) explícitamente para que Nextflow sepa hacer un seguimiento de ese archivo en caso de que sea necesario en pasos posteriores.

### 2.2. Agregar definiciones para entradas accesorias

Dado que nuestro nuevo proceso espera que se proporcionen un puñado de archivos adicionales, configuramos algunos parámetros CLI para ellos bajo la sección `Pipeline parameters`, junto con algunos valores predeterminados (por las mismas razones que antes).

```groovy title="genomics-1.nf" linenums="8"
    // Archivos accesorios
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Crear variables para contener las rutas de archivos accesorios

Mientras que las entradas de datos principales se transmiten dinámicamente a través de canales, hay dos enfoques para manejar archivos accesorios. El enfoque recomendado es crear canales explícitos, lo que hace que el flujo de datos sea más claro y consistente. Alternativamente, se puede usar la función file() para crear variables en casos más simples, particularmente cuando necesita hacer referencia al mismo archivo en múltiples procesos, aunque tenga en cuenta que esto aún crea canales implícitamente. <!-- TODO: Aclarar: ¿esto sigue siendo necesario con entradas tipadas? -->

Agregue esto al bloque workflow (después de la creación de `reads_ch`, dentro de la sección `main:`):

```groovy title="genomics-1.nf" linenums="79"
    // Cargar las rutas de archivo para los archivos accesorios (referencia e intervalos)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

Esto hará que las rutas de archivos accesorios estén disponibles para proporcionar como entrada a cualquier proceso que las necesite.

### 2.4. Agregar una llamada al bloque workflow para ejecutar GATK_HAPLOTYPECALLER

Ahora que tenemos nuestro segundo proceso configurado y todas las entradas y archivos accesorios están listos y disponibles, podemos agregar una llamada al proceso `GATK_HAPLOTYPECALLER` en el cuerpo del workflow.

```groovy title="genomics-1.nf" linenums="88"
    // Llamar variantes del archivo BAM indexado
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

Debería reconocer la sintaxis `*.out` de la Parte 1 de esta serie de entrenamiento; le estamos diciendo a Nextflow que tome la salida de canal de `SAMTOOLS_INDEX` y la conecte a la llamada del proceso `GATK_HAPLOTYPECALLER`.

!!! note

    Notará que las entradas se proporcionan en el mismo orden exacto en la llamada al proceso que se enumeran en el bloque de entrada del proceso.
    En Nextflow, las entradas son posicionales, lo que significa que _debe_ seguir el mismo orden; y por supuesto tiene que haber el mismo número de elementos.

### 2.5. Actualizar la sección publish y el bloque output

Necesitamos actualizar la sección `publish:` para incluir las salidas VCF, y agregar los destinos correspondientes en el bloque `output`.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Ejecutar el workflow para verificar que el paso de llamado de variantes funciona

Ejecutemos el workflow expandido con `-resume` para que no tengamos que ejecutar el paso de indexación nuevamente.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Ahora si miramos la salida de la consola, vemos los dos procesos listados.

El primer proceso fue omitido gracias al almacenamiento en caché, como se esperaba, mientras que el segundo proceso se ejecutó ya que es completamente nuevo.

Encontrará los archivos de salida en el directorio de resultados (como enlaces simbólicos al directorio de trabajo).

??? abstract "Contenido del directorio"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

Si abre el archivo VCF, debería ver el mismo contenido que en el archivo que generó al ejecutar el comando GATK directamente en el contenedor.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Esta es la salida que nos interesa generar para cada muestra en nuestro estudio.

### Conclusión

Sabe cómo hacer un workflow muy básico de dos pasos que hace un trabajo de análisis real y es capaz de lidiar con las idiosincrasias del formato de archivo genómico como los archivos accesorios.

### ¿Qué sigue?

Hacer que el workflow maneje múltiples muestras en lote.

---

## 3. Adaptar el workflow para ejecutarse en un lote de muestras

Está muy bien tener un workflow que pueda automatizar el procesamiento en una sola muestra, pero ¿qué pasa si tiene 1000 muestras?
¿Necesita escribir un script bash que haga un bucle a través de todas sus muestras?

¡No, gracias a Dios! Solo haga un pequeño ajuste al código y Nextflow también manejará eso por usted.

### 3.1. Convertir la declaración de parámetro de entrada en un array que liste las tres muestras

Convirtamos esa ruta de archivo predeterminada en la declaración del archivo BAM de entrada en un array que liste las rutas de archivo para nuestras tres muestras de prueba, arriba bajo la sección `Pipeline parameters`.

=== "Después"

    ```groovy title="genomics-1.nf" linenums="7"
    // Entrada principal (array de tres muestras)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="7"
        // Entrada principal
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note

    Al usar declaraciones de parámetros tipadas (como `reads_bam: Path`), no puede asignar un valor de array.
    Para arrays, omita la anotación de tipo.

Y eso es realmente todo lo que necesitamos hacer, porque la channel factory que usamos en el cuerpo del workflow (`.fromPath`) está tan feliz de aceptar múltiples rutas de archivo para cargar en el canal de entrada como lo estaba para cargar una sola.

!!! note

    Normalmente, no querría codificar la lista de muestras en su archivo de workflow, pero lo estamos haciendo aquí para mantener las cosas simples.
    Presentaremos formas más elegantes de manejar entradas más adelante en esta serie de entrenamiento.

### 3.2. Ejecutar el workflow para verificar que se ejecuta en las tres muestras

Intentemos ejecutar el workflow ahora que la tubería está configurada para ejecutarse en las tres muestras de prueba.

```bash
nextflow run genomics-1.nf -resume
```

Cosa graciosa: esto _podría funcionar_, O _podría fallar_. Por ejemplo, aquí hay una ejecución que tuvo éxito:

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Si su ejecución de workflow tuvo éxito, ejecútelo nuevamente hasta que obtenga un error como este:

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Si mira la salida de error del comando GATK, habrá una línea como esta:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Bueno, eso es extraño, considerando que indexamos explícitamente los archivos BAM en el primer paso del workflow. ¿Podría haber algo mal con la tubería?

#### 3.2.1. Verificar los directorios de trabajo para las llamadas relevantes

Echemos un vistazo dentro del directorio de trabajo para la llamada de proceso `GATK_HAPLOTYPECALLER` fallida listada en la salida de la consola.

??? abstract "Contenido del directorio"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Preste especial atención a los nombres del archivo BAM y el índice BAM que se enumeran en este directorio: `reads_son.bam` y `reads_father.bam.bai`.

¿Qué demonios? Nextflow ha preparado un archivo de índice en el directorio de trabajo de esta llamada de proceso, pero es el incorrecto. ¿Cómo pudo haber sucedido esto?

#### 3.2.2. Usar el [operador view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) para inspeccionar el contenido del canal

Agregue estas dos líneas en el cuerpo del workflow antes de la llamada al proceso `GATK_HAPLOTYPER`:

```groovy title="genomics-1.nf" linenums="84"
    // diagnósticos temporales
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Luego ejecute el comando del workflow nuevamente.

```bash
nextflow run genomics-1.nf
```

Una vez más, esto puede tener éxito o fallar. Aquí hay una ejecución exitosa:

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

Y aquí hay una fallida:

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Es posible que necesite ejecutarlo varias veces para que falle nuevamente.
Este error no se reproducirá de manera consistente porque depende de cierta variabilidad en los tiempos de ejecución de las llamadas de proceso individuales.

Esto es lo que se ve la salida de las dos llamadas `.view()` que agregamos para una ejecución fallida:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Las primeras tres líneas corresponden al canal de entrada y las segundas, al canal de salida.
¡Puede ver que los archivos BAM y los archivos de índice para las tres muestras no están listados en el mismo orden!

!!! note

    Cuando llama a un proceso Nextflow en un canal que contiene múltiples elementos, Nextflow intentará paralelizar la ejecución tanto como sea posible, y recopilará salidas en el orden en que estén disponibles.
    La consecuencia es que las salidas correspondientes pueden recopilarse en un orden diferente al que se alimentaron las entradas originales.

Tal como está escrito actualmente, nuestro script de workflow asume que los archivos de índice saldrán del paso de indexación listados en el mismo orden madre/padre/hijo que se dieron las entradas.
Pero eso no está garantizado, razón por la cual a veces (aunque no siempre) los archivos incorrectos se emparejan en el segundo paso.

Para solucionar esto, necesitamos asegurarnos de que los archivos BAM y sus archivos de índice viajen juntos a través de los canales.

!!! tip

    Las declaraciones `view()` en el código del workflow no hacen nada, por lo que no es un problema dejarlas.
    Sin embargo, desordenarán su salida de consola, por lo que recomendamos eliminarlas cuando termine de solucionar el problema.

### 3.3. Cambiar la salida del proceso SAMTOOLS_INDEX en una tupla que mantenga el archivo de entrada y su índice juntos

La forma más simple de asegurar que un archivo BAM y su índice permanezcan estrechamente asociados es empaquetarlos juntos en una tupla que salga de la tarea de índice.

!!! note

    Una **tupla** es una lista ordenada y finita de elementos que se usa comúnmente para devolver múltiples valores de una función. Las tuplas son particularmente útiles para pasar múltiples entradas o salidas entre procesos mientras se preserva su asociación y orden.

Primero, cambiemos la salida del proceso `SAMTOOLS_INDEX` para incluir el archivo BAM en su declaración de salida.

=== "Después"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

De esta manera, cada archivo de índice estará estrechamente acoplado con su archivo BAM original, y la salida general del paso de indexación será un solo canal que contiene pares de archivos.

### 3.4. Cambiar la entrada al proceso GATK_HAPLOTYPECALLER para que sea una tupla

Dado que hemos cambiado la 'forma' de la salida del primer proceso en el workflow, necesitamos actualizar la definición de entrada del segundo proceso para que coincida.

Específicamente, donde previamente declaramos dos rutas de entrada separadas en el bloque de entrada del proceso `GATK_HAPLOTYPECALLER`, ahora declaramos una sola entrada que coincide con la estructura de la tupla emitida por `SAMTOOLS_INDEX`.

=== "Después"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

Por supuesto, dado que ahora hemos cambiado la forma de las entradas que `GATK_HAPLOTYPECALLER` espera, necesitamos actualizar la llamada del proceso en consecuencia en el cuerpo del workflow.

### 3.5. Actualizar la llamada a GATK_HAPLOTYPECALLER en el bloque workflow

Ya no necesitamos proporcionar el `reads_ch` original al proceso `GATK_HAPLOTYPECALLER`, ya que el archivo BAM ahora está empaquetado en la salida del canal por `SAMTOOLS_INDEX`.

Como resultado, simplemente podemos eliminar esa línea.

=== "Después"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

Esa es toda la reconexión que es necesaria para resolver el problema de desajuste de índice.

### 3.6. Actualizar la sección publish y el bloque output para la tupla

Dado que `SAMTOOLS_INDEX.out` ahora es una tupla que contiene tanto el BAM como su índice, ambos archivos se publicarán juntos.
Renombramos el destino de `bam_index` a `indexed_bam` para reflejar que ahora contiene ambos archivos.

=== "Después"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Antes"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

También necesitamos actualizar el bloque output para usar el nuevo nombre de destino:

=== "Después"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Antes"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Ejecutar el workflow para verificar que funciona correctamente en las tres muestras cada vez

Por supuesto, la prueba está en el resultado, así que ejecutemos el workflow nuevamente algunas veces para asegurarnos de que esto funcione de manera confiable en el futuro.

```bash
nextflow run genomics-1.nf
```

Esta vez (y cada vez) todo debería ejecutarse correctamente:

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

El directorio de resultados ahora contiene archivos BAM y BAI para cada muestra (de la tupla), junto con las salidas VCF:

??? abstract "Contenido del directorio de resultados"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

Si lo desea, puede usar `.view()` nuevamente para ver cómo se ve el contenido del canal de salida de `SAMTOOLS_INDEX`:

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Verá que el canal contiene las tres tuplas esperadas (rutas de archivo truncadas para legibilidad).

```console title="Salida"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

Eso será mucho más seguro, en adelante.

### Conclusión

Sabe cómo hacer que su workflow se ejecute en múltiples muestras (independientemente).

### ¿Qué sigue?

Facilitar el manejo de muestras en lote.

---

## 4. Hacer que el workflow acepte un archivo de texto que contenga un lote de archivos de entrada

Una forma muy común de proporcionar múltiples archivos de datos de entrada a un workflow es hacerlo con un archivo de texto que contenga las rutas de los archivos.
Puede ser tan simple como un archivo de texto que liste una ruta de archivo por línea y nada más, o el archivo puede contener metadatos adicionales, en cuyo caso a menudo se llama samplesheet.

Aquí vamos a mostrarle cómo hacer el caso simple.

### 4.1. Examinar el archivo de texto proporcionado que lista las rutas de los archivos de entrada

Ya hicimos un archivo de texto que lista las rutas de los archivos de entrada, llamado `sample_bams.txt`, que puede encontrar en el directorio `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Como puede ver, listamos una ruta de archivo por línea, y son rutas absolutas.

!!! note

    Los archivos que estamos usando aquí están simplemente en el sistema de archivos local de su GitHub Codespaces, pero también podríamos apuntar a archivos en almacenamiento en la nube.

### 4.2. Actualizar el valor predeterminado del parámetro

Cambiemos el valor predeterminado para nuestro parámetro de entrada `reads_bam` para que apunte al archivo `sample_bams.txt`.

=== "Después"

    ```groovy title="genomics-1.nf" linenums="7"
        // Entrada principal (archivo de archivos de entrada, uno por línea)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="7"
    // Entrada principal (array de tres muestras)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

De esta manera podemos continuar siendo perezosos, pero la lista de archivos ya no vive en el código del workflow en sí, lo cual es un gran paso en la dirección correcta.

### 4.3. Actualizar la channel factory para leer líneas de un archivo

Actualmente, nuestra channel factory de entrada trata cualquier archivo que le demos como las entradas de datos que queremos alimentar al proceso de indexación.
Dado que ahora le estamos dando un archivo que lista las rutas de los archivos de entrada, necesitamos cambiar su comportamiento para analizar el archivo y tratar las rutas de archivo que contiene como las entradas de datos.

Afortunadamente, podemos hacer eso muy simplemente, solo agregando el [operador `.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) al paso de construcción del canal.

=== "Después"

    ```groovy title="genomics-1.nf" linenums="68"
        // Crear canal de entrada desde un archivo de texto que lista rutas de archivos de entrada
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Antes"

    ```groovy title="genomics-1.nf" linenums="68"
        // Crear canal de entrada (archivo único vía parámetro CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip

    Esta es otra gran oportunidad para usar el operador `.view()` para ver cómo se ve el contenido del canal antes y después de aplicar un operador.

### 4.4. Ejecutar el workflow para verificar que funciona correctamente

Ejecutemos el workflow una vez más. Esto debería producir el mismo resultado que antes, ¿verdad?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

¡Sí! De hecho, Nextflow detecta correctamente que las llamadas de proceso son exactamente las mismas, y ni siquiera se molesta en volver a ejecutar todo, ya que estábamos ejecutando con `-resume`.

¡Y eso es todo! Nuestro workflow simple de llamado de variantes tiene todas las características básicas que queríamos.

### Conclusión

Sabe cómo hacer un workflow lineal de múltiples pasos para indexar un archivo BAM y aplicar el llamado de variantes por muestra usando GATK.

De manera más general, ha aprendido cómo usar componentes y lógica esenciales de Nextflow para construir un pipeline genómico simple que hace trabajo real, teniendo en cuenta las idiosincrasias de los formatos de archivo genómicos y los requisitos de las herramientas.

### ¿Qué sigue?

¡Celebre su éxito y tome un descanso extra largo!

En la próxima parte de este curso, aprenderá cómo usar algunas características adicionales de Nextflow (incluyendo más operadores de canal) para aplicar el llamado de variantes conjunto a los datos.
