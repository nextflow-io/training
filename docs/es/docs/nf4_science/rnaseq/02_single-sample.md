# Parte 2: Implementación de muestra única

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta parte del curso, vamos a desarrollar un flujo de trabajo que automatice los comandos que ejecutamos en la Parte 1, enfocándonos en procesar una muestra a la vez.

!!! warning "Requisito previo"

    Debe completar la [Parte 1: Descripción del método](./01_method.md) antes de comenzar esta lección.
    Específicamente, trabajar en la sección 1.2.3 crea el archivo de índice del genoma (`data/genome_index.tar.gz`) requerido para el paso de alineamiento en esta lección.

## Asignación

En esta parte del curso, vamos a desarrollar un flujo de trabajo que haga lo siguiente:

1. Ejecutar control de calidad (FastQC) en las lecturas de entrada
2. Recortar adaptadores y ejecutar control de calidad posterior al recorte (Trim Galore)
3. Alinear las lecturas recortadas a un genoma de referencia (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

Esto automatiza los pasos de la primera sección de la [Parte 1: Descripción del método](./01_method.md#1-single-sample-processing), donde ejecutó estos comandos manualmente en sus contenedores.

Como punto de partida, le proporcionamos un archivo de flujo de trabajo, `rnaseq.nf`, que describe las partes principales del flujo de trabajo, así como cuatro archivos de módulo en el directorio `modules/` (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf` y `multiqc.nf`) que describen la estructura de cada proceso.

??? full-code "Archivos de esqueleto"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Module INCLUDE statements

    /*
     * Pipeline parameters
     */

    // Primary input

    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }

    output {
        // Configure publish targets
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
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

Estos archivos no son funcionales; su propósito es solo servir como esqueletos para que usted complete con las partes interesantes del código.

## Plan de la lección

Para hacer el proceso de desarrollo más educativo, lo hemos dividido en tres etapas:

1. **Escribir un flujo de trabajo de una sola etapa que ejecute el paso de control de calidad inicial.**
   Esto cubre la configuración de un parámetro CLI, la creación de un canal de entrada, la escritura de un módulo de proceso y la configuración de publicación de salidas.
2. **Agregar recorte de adaptadores y control de calidad posterior al recorte.**
   Esto introduce el encadenamiento de procesos conectando la salida de un proceso con la entrada de otro.
3. **Agregar alineamiento al genoma de referencia.**
   Esto cubre el manejo de entradas de referencia adicionales y el trabajo con archivos comprimidos.

Cada paso se enfoca en un aspecto específico del desarrollo de flujos de trabajo.

!!! tip "Consejo"

     Asegúrese de estar en el directorio de trabajo correcto:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Escribir un flujo de trabajo de una sola etapa que ejecute el control de calidad inicial

Este primer paso se enfoca en lo básico: cargar un archivo FASTQ y ejecutar control de calidad en él.

Recuerde el comando `fastqc` de la [Parte 1](01_method.md):

```bash
fastqc <reads>
```

El comando toma un archivo FASTQ como entrada y produce un informe de control de calidad como un archivo `.zip` y un resumen `.html`.
El URI del contenedor era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Vamos a tomar esta información y envolverla en Nextflow en tres etapas:

1. Configurar la entrada
2. Escribir el proceso de control de calidad y llamarlo en el flujo de trabajo
3. Configurar el manejo de salidas

### 1.1. Configurar la entrada

Necesitamos declarar un parámetro de entrada, crear un perfil de prueba para proporcionar un valor predeterminado conveniente y crear un canal de entrada.

#### 1.1.1. Agregar una declaración de parámetro de entrada

En `rnaseq.nf`, bajo la sección `Pipeline parameters`, declare un parámetro llamado `reads` con el tipo `Path`.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        input: Path
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Eso configura el parámetro CLI, pero no queremos escribir la ruta del archivo cada vez que ejecutamos el flujo de trabajo durante el desarrollo.
Hay múltiples opciones para proporcionar un valor predeterminado; aquí usamos un perfil de prueba.

#### 1.1.2. Crear un perfil de prueba con un valor predeterminado en `nextflow.config`

Un perfil de prueba proporciona valores predeterminados convenientes para probar un flujo de trabajo sin especificar entradas en la línea de comandos.
Esta es una convención común en el ecosistema Nextflow (consulte [Hello Config](../../hello_nextflow/06_hello_config.md) para más detalles).

Agregue un bloque `profiles` a `nextflow.config` con un perfil `test` que establezca el parámetro `reads` a uno de los archivos FASTQ de prueba.

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aquí, estamos usando `#!groovy ${projectDir}`, una variable integrada de Nextflow que apunta al directorio donde se encuentra el script del flujo de trabajo.
Esto facilita referenciar archivos de datos y otros recursos sin codificar rutas absolutas.

El parámetro ahora tiene un valor predeterminado conveniente. A continuación, necesitamos crear un canal a partir de él.

#### 1.1.3. Configurar el canal de entrada

En el bloque workflow, cree un canal de entrada a partir del valor del parámetro usando la fábrica de canal `.fromPath` (como se usa en [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Después"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

A continuación, necesitaremos crear el proceso para ejecutar control de calidad en esta entrada.

### 1.2. Escribir el proceso de control de calidad y llamarlo en el flujo de trabajo

Necesitamos completar la definición del proceso en el archivo de módulo, importarlo al flujo de trabajo usando una declaración include y llamarlo en la entrada.

#### 1.2.1. Completar el módulo para el proceso de control de calidad

Abra `modules/fastqc.nf` y examine el esquema de la definición del proceso.
Debería reconocer los elementos estructurales principales; si no, considere leer [Hello Nextflow](../../hello_nextflow/01_hello_world.md) para refrescar la memoria.

Adelante, complete la definición del proceso por su cuenta usando la información proporcionada arriba, luego verifique su trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

El accesor `simpleName` elimina todas las extensiones del nombre del archivo, por lo que `ENCSR000COQ1_1.fastq.gz` se convierte en `ENCSR000COQ1_1`.
Usamos la sintaxis `emit:` para asignar nombres a cada canal de salida, lo que será útil para conectar las salidas al bloque publish.

Una vez que haya completado esto, el proceso está completo.
Para usarlo en el flujo de trabajo, necesitará importar el módulo y agregar una llamada al proceso.

#### 1.2.2. Incluir el módulo

En `rnaseq.nf`, agregue una declaración `include` para hacer que el proceso esté disponible para el flujo de trabajo:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    ```

El proceso ahora está disponible en el ámbito del flujo de trabajo.

#### 1.2.3. Llamar al proceso de control de calidad en la entrada

Agregue una llamada a `FASTQC` en el bloque workflow, pasando el canal de entrada como argumento.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Control de calidad inicial
        FASTQC(read_ch)

        publish:
        // Declare outputs to publish
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

El flujo de trabajo ahora carga la entrada y ejecuta el proceso de control de calidad en ella.
A continuación, necesitamos configurar cómo se publica la salida.

### 1.3. Configurar el manejo de salidas

Necesitamos declarar qué salidas de proceso publicar y especificar dónde deben ir.

#### 1.3.1. Declarar salidas en la sección `publish:`

La sección `publish:` dentro del bloque workflow declara qué salidas de proceso deben publicarse.
Asigne las salidas de `FASTQC` a destinos nombrados.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Declare outputs to publish
    }
    ```

A continuación, necesitaremos decirle a Nextflow dónde colocar las salidas publicadas.

#### 1.3.2. Configurar los destinos de salida en el bloque `output {}`

El bloque `output {}` se encuentra fuera del flujo de trabajo y especifica dónde se publica cada destino nombrado.
Configure ambos destinos para publicar en un subdirectorio `fastqc/`.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configure publish targets
    }
    ```

!!! note "Nota"

    Por defecto, Nextflow publica archivos de salida como enlaces simbólicos, lo que evita duplicación innecesaria.
    Aunque los archivos de datos que estamos usando aquí son muy pequeños, en genómica pueden volverse muy grandes.
    Los enlaces simbólicos se romperán cuando limpie su directorio `work`, por lo que para flujos de trabajo de producción es posible que desee anular el modo de publicación predeterminado a `'copy'`.

### 1.4. Ejecutar el flujo de trabajo

En este punto, tenemos un flujo de trabajo de control de calidad de un paso que debería ser completamente funcional.

Ejecutamos con `-profile test` para usar el valor predeterminado configurado en el perfil de prueba, evitando la necesidad de escribir la ruta en la línea de comandos.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

Esto debería ejecutarse muy rápidamente si trabajó en la Parte 1 y ya ha descargado el contenedor.
Si la omitió, Nextflow descargará el contenedor por usted; no tiene que hacer nada para que suceda, pero es posible que deba esperar hasta un minuto.

Puede verificar las salidas en el directorio de resultados.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

Los informes de control de calidad para la muestra ahora están publicados en el subdirectorio `fastqc/`.

### Conclusión

Sabe cómo crear un módulo que contenga un proceso, importarlo a un flujo de trabajo, llamarlo con un canal de entrada y publicar los resultados usando el bloque de salida a nivel de flujo de trabajo.

### ¿Qué sigue?

Agregue recorte de adaptadores con control de calidad posterior al recorte como un segundo paso en el flujo de trabajo.

---

## 2. Agregar recorte de adaptadores y control de calidad posterior al recorte

Ahora que tenemos el control de calidad inicial en su lugar, podemos agregar el paso de recorte de adaptadores con su control de calidad posterior al recorte integrado.

Recuerde el comando `trim_galore` de la [Parte 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

El comando recorta adaptadores de un archivo FASTQ y ejecuta FastQC en la salida recortada.
Produce lecturas recortadas, un informe de recorte e informes FastQC para las lecturas recortadas.
El URI del contenedor era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Solo necesitamos escribir la definición del proceso, importarlo, llamarlo en el flujo de trabajo y actualizar el manejo de salidas.

### 2.1. Escribir el proceso de recorte y llamarlo en el flujo de trabajo

Como antes, necesitamos completar la definición del proceso, importar el módulo y agregar la llamada al proceso.

#### 2.1.1. Completar el módulo para el proceso de recorte

Abra `modules/trim_galore.nf` y examine el esquema de la definición del proceso.

Adelante, complete la definición del proceso por su cuenta usando la información proporcionada arriba, luego verifique su trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

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
    }
    ```

Este proceso tiene tres salidas nombradas: las lecturas recortadas que alimentan el paso de alineamiento, el informe de recorte y los informes FastQC posteriores al recorte.
La bandera `--fastqc` le dice a Trim Galore que ejecute automáticamente FastQC en la salida recortada.

#### 2.1.2. Incluir el módulo

Actualice `rnaseq.nf` para importar el nuevo módulo:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

A continuación, agregaremos la llamada al proceso al flujo de trabajo.

#### 2.1.3. Llamar al proceso de recorte en la entrada

Agregue la llamada al proceso en el bloque workflow:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Control de calidad inicial
        FASTQC(read_ch)

        // Recorte de adaptadores y control de calidad posterior al recorte
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Control de calidad inicial
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

El proceso de recorte ahora está conectado al flujo de trabajo.

### 2.2. Actualizar el manejo de salidas

Necesitamos agregar las salidas de recorte a la declaración de publicación y configurar dónde van.

#### 2.2.1. Agregar destinos de publicación para las salidas de recorte

Agregue las salidas de recorte a la sección `publish:`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

A continuación, necesitaremos decirle a Nextflow dónde colocar estas salidas.

#### 2.2.2. Configurar los nuevos destinos de salida

Agregue entradas para los destinos de recorte en el bloque `output {}`, publicándolos en un subdirectorio `trimming/`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

La configuración de salida está completa.

### 2.3. Ejecutar el flujo de trabajo

El flujo de trabajo ahora incluye tanto el control de calidad inicial como el recorte de adaptadores.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

Esto también debería ejecutarse muy rápidamente, ya que estamos ejecutando en un archivo de entrada tan pequeño.

Puede encontrar las salidas de recorte en el directorio de resultados.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

Las salidas de recorte y los informes de control de calidad posterior al recorte ahora están en el subdirectorio `trimming/`.

### Conclusión

Sabe cómo agregar un segundo paso de procesamiento que se ejecuta independientemente en la misma entrada, produciendo múltiples salidas nombradas.

### ¿Qué sigue?

Agregue el paso de alineamiento que se encadena a partir de la salida de lecturas recortadas.

---

## 3. Agregar alineamiento al genoma de referencia

Finalmente podemos agregar el paso de alineamiento del genoma usando HISAT2.

Recuerde el comando de alineamiento de la [Parte 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

El comando alinea lecturas a un genoma de referencia y convierte la salida a formato BAM.
Requiere un archivo de índice del genoma preconstruido y produce un archivo BAM y un registro de resumen de alineamiento.
El URI del contenedor era `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

Este proceso requiere una entrada adicional (el archivo de índice del genoma), por lo que necesitamos configurar eso primero, luego escribir y conectar el proceso.

### 3.1. Configurar las entradas

Necesitamos declarar un parámetro para el archivo de índice del genoma.

#### 3.1.1. Agregar un parámetro para el índice del genoma

Agregue una declaración de parámetro para el archivo de índice del genoma en `rnaseq.nf`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Primary input
        input: Path

        // Reference genome archive
        hisat2_index_zip: Path
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primary input
        input: Path
    }
    ```

#### 3.1.2. Agregar el índice del genoma predeterminado al perfil de prueba

Tal como hicimos para `reads` en la sección 1.1.2, agregue un valor predeterminado para el índice del genoma al perfil de prueba en `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
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
        }
    }
    ```

El parámetro está listo; ahora podemos crear el proceso de alineamiento.

### 3.2. Escribir el proceso de alineamiento y llamarlo en el flujo de trabajo

Como antes, necesitamos completar la definición del proceso, importar el módulo y agregar la llamada al proceso.

#### 3.2.1. Completar el módulo para el proceso de alineamiento

Abra `modules/hisat2_align.nf` y examine el esquema de la definición del proceso.

Adelante, complete la definición del proceso por su cuenta usando la información proporcionada arriba, luego verifique su trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

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
    }
    ```

Este proceso toma dos entradas: las lecturas y el archivo de índice del genoma.
El bloque script primero extrae el índice del archivo, luego ejecuta el alineamiento HISAT2 canalizado a `samtools view` para convertir la salida a formato BAM.
El accesor `simpleName` en `index_zip` extrae el nombre base del archivo (`genome_index`) para usar como prefijo del índice.

#### 3.2.2. Incluir el módulo

Actualice `rnaseq.nf` para importar el nuevo módulo:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

A continuación, agregaremos la llamada al proceso al flujo de trabajo.

#### 3.2.3. Llamar al proceso de alineamiento

Las lecturas recortadas están en el canal `TRIM_GALORE.out.trimmed_reads` de salida del paso anterior.
Usamos `#!groovy file(params.hisat2_index_zip)` para proporcionar el archivo de índice del genoma.

=== "Después"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Control de calidad inicial
        FASTQC(read_ch)

        // Recorte de adaptadores y control de calidad posterior al recorte
        TRIM_GALORE(read_ch)

        // Alineamiento a un genoma de referencia
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Crear canal de entrada a partir de una ruta de archivo
        read_ch = channel.fromPath(params.input)

        // Control de calidad inicial
        FASTQC(read_ch)

        // Recorte de adaptadores y control de calidad posterior al recorte
        TRIM_GALORE(read_ch)
    ```

El proceso de alineamiento ahora está conectado al flujo de trabajo.

### 3.3. Actualizar el manejo de salidas

Necesitamos agregar las salidas de alineamiento a la declaración de publicación y configurar dónde van.

#### 3.3.1. Agregar destinos de publicación para las salidas de alineamiento

Agregue las salidas de alineamiento a la sección `publish:`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
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

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

A continuación, necesitaremos decirle a Nextflow dónde colocar estas salidas.

#### 3.3.2. Configurar los nuevos destinos de salida

Agregue entradas para los destinos de alineamiento en el bloque `output {}`, publicándolos en un subdirectorio `align/`:

=== "Después"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

La configuración de salida está completa.

### 3.4. Ejecutar el flujo de trabajo

El flujo de trabajo ahora incluye los tres pasos de procesamiento: control de calidad, recorte y alineamiento.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

Puede encontrar las salidas de alineamiento en el directorio de resultados.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Esto completa el procesamiento básico que necesitamos aplicar a cada muestra.

_Agregaremos la agregación de informes MultiQC en la Parte 3, después de que hayamos modificado el flujo de trabajo para aceptar múltiples muestras a la vez._

---

### Conclusión

Sabe cómo envolver todos los pasos principales para procesar muestras de RNAseq de extremo simple individualmente.

### ¿Qué sigue?

¡Tome un descanso! Eso fue mucho.

Cuando se sienta renovado, diríjase a la [Parte 3](./03_multi-sample.md), donde aprenderá cómo modificar el flujo de trabajo para procesar múltiples muestras en paralelo, agregar informes de control de calidad en todos los pasos para todas las muestras y habilitar la ejecución del flujo de trabajo en datos de RNAseq de extremo pareado.
