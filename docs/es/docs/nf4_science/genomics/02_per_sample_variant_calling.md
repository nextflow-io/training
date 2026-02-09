# Parte 2: Llamado de variantes por muestra

En la Parte 1, probaste los comandos de Samtools y GATK manualmente en sus respectivos contenedores.
Ahora vamos a envolver esos mismos comandos en un workflow de Nextflow.

## Asignación

En esta parte del curso, vamos a desarrollar un workflow que hace lo siguiente:

1. Generar un archivo de índice para cada archivo BAM de entrada usando [Samtools](https://www.htslib.org/)
2. Ejecutar GATK HaplotypeCaller en cada archivo BAM de entrada para generar llamados de variantes por muestra en formato VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Esto replica los pasos de la Parte 1, donde ejecutaste estos comandos manualmente en sus contenedores.

Como punto de partida, te proporcionamos un archivo de workflow, `genomics.nf`, que describe las partes principales del workflow, así como dos archivos de módulo, samtools_index.nf y gatk_haplotypecaller.nf, que describen la estructura de los módulos.
Estos archivos no son funcionales; su propósito es simplemente servir como esqueletos para que los completes con las partes interesantes del código.

## Plan de la lección

Para hacer el proceso de desarrollo más educativo, lo hemos dividido en cuatro pasos:

1. **Escribir un workflow de una sola etapa que ejecute Samtools index en un archivo BAM.**
   Esto cubre crear un módulo, importarlo y llamarlo en un workflow.
2. **Agregar un segundo proceso para ejecutar GATK HaplotypeCaller en el archivo BAM indexado.**
   Esto introduce el encadenamiento de salidas de procesos a entradas y el manejo de archivos accesorios.
3. **Adaptar el workflow para ejecutarse en un lote de muestras.**
   Esto cubre la ejecución paralela e introduce tuplas para mantener archivos asociados juntos.
4. **Hacer que el workflow acepte un archivo de texto que contenga un lote de archivos de entrada.**
   Esto demuestra un patrón común para proporcionar entradas en lote.

Cada paso se enfoca en un aspecto específico del desarrollo de workflows.

---

## 1. Escribir un workflow de una sola etapa que ejecute Samtools index en un archivo BAM

Este primer paso se enfoca en lo básico: cargar un archivo BAM y generar un índice para él.

Recuerda el comando `samtools index` de la [Parte 1](01_method.md):

```bash
samtools index '<input_bam>'
```

El comando toma un archivo BAM como entrada y produce un archivo de índice `.bai` junto a él.
El URI del contenedor era `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Vamos a tomar esta información y envolverla en Nextflow en tres etapas:

1. Configurar la entrada
2. Escribir el proceso de indexación y llamarlo en el workflow
3. Configurar el manejo de salida

### 1.1. Configurar la entrada

Necesitamos declarar un parámetro de entrada, crear un perfil de prueba para proporcionar un valor predeterminado conveniente, y crear un canal de entrada.

#### 1.1.1. Agregar una declaración de parámetro de entrada

En el archivo principal de workflow `genomics.nf`, bajo la sección `Pipeline parameters`, declara un parámetro CLI llamado `reads_bam`.

=== "Después"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        reads_bam: Path
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Eso configura el parámetro CLI, pero no queremos escribir la ruta del archivo cada vez que ejecutemos el workflow durante el desarrollo.
Hay múltiples opciones para proporcionar un valor predeterminado; aquí usamos un perfil de prueba.

#### 1.1.2. Crear un perfil de prueba con un valor predeterminado en `nextflow.config`

Un perfil de prueba proporciona valores predeterminados convenientes para probar un workflow sin especificar entradas en la línea de comandos.
Esta es una convención común en el ecosistema de Nextflow (consulta [Hello Config](../../hello_nextflow/06_hello_config.md) para más detalles).

Agrega un bloque `profiles` a `nextflow.config` con un perfil `test` que establezca el parámetro `reads_bam` a uno de los archivos BAM de prueba.

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aquí estamos usando `${projectDir}`, una variable integrada de Nextflow que apunta al directorio donde se encuentra el script del workflow.
Esto facilita referenciar archivos de datos y otros recursos sin codificar rutas absolutas.

#### 1.1.3. Configurar el canal de entrada

En el bloque workflow, crea un canal de entrada desde el valor del parámetro usando el channel factory `.fromPath` (como se usa en [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Después"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Create input channel
    ```

Ahora necesitamos crear el proceso para ejecutar la indexación en esta entrada.

### 1.2. Escribir el proceso de indexación y llamarlo en el workflow

Necesitamos escribir la definición del proceso en el archivo del módulo, importarlo al workflow usando una declaración include, y llamarlo en la entrada.

#### 1.2.1. Completar el módulo para el proceso de indexación

Abre `modules/samtools_index.nf` y examina el esquema de la definición del proceso.
Deberías reconocer los elementos estructurales principales; si no, considera leer [Hello Nextflow](../../hello_nextflow/01_hello_world.md) para refrescar la memoria.

Adelante, completa la definición del proceso por ti mismo usando la información proporcionada arriba, luego verifica tu trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Generate BAM index file
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

Una vez que hayas completado esto, el proceso está completo.
Para usarlo en el workflow, necesitarás importar el módulo y agregar una llamada al proceso.

#### 1.2.2. Incluir el módulo

En `genomics.nf`, agrega una declaración `include` para hacer el proceso disponible al workflow:

=== "Después"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    ```

El proceso ahora está disponible en el ámbito del workflow.

#### 1.2.3. Llamar el proceso de indexación en la entrada

Ahora, agreguemos una llamada a `SAMTOOLS_INDEX` en el bloque workflow, pasando el canal de entrada como argumento.

=== "Después"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Call processes
    ```

El workflow ahora carga la entrada y ejecuta el proceso de indexación en ella.
A continuación, necesitamos configurar cómo se publica la salida.

### 1.3. Configurar el manejo de salida

Necesitamos declarar qué salidas de proceso publicar y especificar dónde deben ir.

#### 1.3.1. Declarar una salida en la sección `publish:`

La sección `publish:` dentro del bloque workflow declara qué salidas de proceso deben publicarse.
Asigna la salida de `SAMTOOLS_INDEX` a un objetivo nombrado llamado `bam_index`.

=== "Después"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declare outputs to publish
    }
    ```

Ahora necesitamos decirle a Nextflow dónde colocar la salida publicada.

#### 1.3.2. Configurar el objetivo de salida en el bloque `output {}`

El bloque `output {}` se encuentra fuera del workflow y especifica dónde se publica cada objetivo nombrado.
Agreguemos un objetivo para `bam_index` que publique en un subdirectorio `bam/`.

=== "Después"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configure publish targets
    }
    ```

!!! note "Nota"

    Por defecto, Nextflow publica archivos de salida como enlaces simbólicos, lo que evita duplicación innecesaria.
    Aunque los archivos de datos que estamos usando aquí son muy pequeños, en genómica pueden ser muy grandes.
    Los enlaces simbólicos se romperán cuando limpies tu directorio `work`, así que para workflows de producción puedes querer sobrescribir el modo de publicación predeterminado a `'copy'`.

### 1.4. Ejecutar el workflow

En este punto, tenemos un workflow de indexación de un paso que debería ser completamente funcional. ¡Probemos que funciona!

Podemos ejecutarlo con `-profile test` para usar el valor predeterminado configurado en el perfil de prueba y evitar tener que escribir la ruta en la línea de comandos.

```bash
nextflow run genomics.nf -profile test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Puedes verificar que el archivo de índice se haya generado correctamente mirando en el directorio work o en el directorio de resultados.

??? abstract "Contenido del directorio work"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Contenido del directorio de resultados"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

¡Ahí está!

### Conclusión

Sabes cómo crear un módulo que contiene un proceso, importarlo a un workflow, llamarlo con un canal de entrada y publicar los resultados.

### ¿Qué sigue?

Agregar un segundo paso que tome la salida del proceso de indexación y lo use para ejecutar el llamado de variantes.

---

## 2. Agregar un segundo proceso para ejecutar GATK HaplotypeCaller en el archivo BAM indexado

Ahora que tenemos un índice para nuestro archivo de entrada, podemos pasar a configurar el paso de llamado de variantes.

Recuerda el comando `gatk HaplotypeCaller` de la [Parte 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

El comando toma un archivo BAM (`-I`), un genoma de referencia (`-R`), y un archivo de intervalos (`-L`), y produce un archivo VCF (`-O`) junto con su índice.
La herramienta también espera que el índice BAM, el índice de referencia y el diccionario de referencia estén ubicados junto a sus respectivos archivos.
El URI del contenedor era `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Seguimos las mismas tres etapas que antes:

1. Configurar las entradas
2. Escribir el proceso de llamado de variantes y llamarlo en el workflow
3. Configurar el manejo de salida

### 2.1. Configurar las entradas

El paso de llamado de variantes requiere varios archivos de entrada adicionales.
Necesitamos declarar parámetros para ellos, agregar valores predeterminados al perfil de prueba, y crear variables para cargarlos.

#### 2.1.1. Agregar declaraciones de parámetros para entradas accesorias

Dado que nuestro nuevo proceso espera un puñado de archivos adicionales, agrega declaraciones de parámetros para ellos en `genomics.nf` bajo la sección `Pipeline parameters`:

=== "Después"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Primary input
        reads_bam: Path

        // Accessory files
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Primary input
        reads_bam: Path
    }
    ```

Como antes, proporcionamos valores predeterminados a través del perfil de prueba en lugar de en línea.

#### 2.1.2. Agregar valores predeterminados de archivos accesorios al perfil de prueba

Así como hicimos para `reads_bam` en la sección 1.1.2, agrega valores predeterminados para los archivos accesorios al perfil de prueba en `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

Ahora necesitamos crear variables que carguen estas rutas de archivos para usar en el workflow.

#### 2.1.3. Crear variables para los archivos accesorios

Agrega variables para las rutas de archivos accesorios dentro del bloque workflow:

=== "Después"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Load the file paths for the accessory files (reference and intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)

        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

La sintaxis `file()` le dice a Nextflow explícitamente que maneje estas entradas como rutas de archivos.
Puedes aprender más sobre esto en la Misión Secundaria [Working with files](../../side_quests/working_with_files.md).

### 2.2. Escribir el proceso de llamado de variantes y llamarlo en el workflow

Necesitamos escribir la definición del proceso en el archivo del módulo, importarlo al workflow usando una declaración include, y llamarlo en las lecturas de entrada más la salida del paso de indexación y los archivos accesorios.

#### 2.2.1. Completar el módulo para el proceso de llamado de variantes

Abre `modules/gatk_haplotypecaller.nf` y examina el esquema de la definición del proceso.

Adelante, completa la definición del proceso por ti mismo usando la información proporcionada arriba, luego verifica tu trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Call variants with GATK HaplotypeCaller
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

Notarás que este proceso tiene más entradas de las que el comando GATK en sí requiere.
El GATK sabe buscar el archivo de índice BAM y los archivos accesorios del genoma de referencia basándose en convenciones de nomenclatura, pero Nextflow es independiente del dominio y no conoce estas convenciones.
Necesitamos listarlos explícitamente para que Nextflow los prepare en el directorio de trabajo en tiempo de ejecución; de lo contrario, GATK arrojará un error sobre archivos faltantes.

De manera similar, listamos el archivo de índice del VCF de salida (`"${input_bam}.vcf.idx"`) explícitamente para que Nextflow lleve un seguimiento de él para pasos subsecuentes.
Usamos la sintaxis `emit:` para asignar un nombre a cada canal de salida, lo que será útil cuando conectemos las salidas al bloque publish.

Una vez que hayas completado esto, el proceso está completo.
Para usarlo en el workflow, necesitarás importar el módulo y agregar una llamada al proceso.

#### 2.2.2. Importar el nuevo módulo

Actualiza `genomics.nf` para importar el nuevo módulo:

=== "Después"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="3"
    // Module INCLUDE statements
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

El proceso ahora está disponible en el ámbito del workflow.

#### 2.2.3. Agregar la llamada al proceso

Agrega la llamada al proceso en el cuerpo del workflow, bajo `main:`:

=== "Después"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="33"
        // Create index file for input BAM file
        SAMTOOLS_INDEX(reads_ch)
    ```

Deberías reconocer la sintaxis `*.out` de la serie de capacitación Hello Nextflow; le estamos diciendo a Nextflow que tome el canal de salida de `SAMTOOLS_INDEX` y lo conecte a la llamada del proceso `GATK_HAPLOTYPECALLER`.

!!! note "Nota"

    Observa que las entradas se proporcionan en exactamente el mismo orden en la llamada al proceso que están listadas en el bloque input del proceso.
    En Nextflow, las entradas son posicionales, lo que significa que _debes_ seguir el mismo orden; y por supuesto debe haber el mismo número de elementos.

### 2.3. Configurar el manejo de salida

Necesitamos agregar las nuevas salidas a la declaración publish y configurar dónde van.

#### 2.3.1. Agregar objetivos de publicación para las salidas de llamado de variantes

Agrega las salidas VCF e índice a la sección `publish:`:

=== "Después"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

Ahora necesitamos decirle a Nextflow dónde colocar las nuevas salidas.

#### 2.3.2. Configurar los nuevos objetivos de salida

Agrega entradas para los objetivos `vcf` y `vcf_idx` en el bloque `output {}`, publicando ambos en un subdirectorio `vcf/`:

=== "Después"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

El VCF y su índice se publican como objetivos separados que ambos van al subdirectorio `vcf/`.

### 2.4. Ejecutar el workflow

Ejecuta el workflow expandido, agregando `-resume` esta vez para que no tengamos que ejecutar el paso de indexación nuevamente.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Ahora si miramos la salida de la consola, vemos los dos procesos listados.

El primer proceso fue omitido gracias al caché, como se esperaba, mientras que el segundo proceso se ejecutó ya que es completamente nuevo.

Encontrarás los archivos de salida en el directorio de resultados (como enlaces simbólicos al directorio work).

??? abstract "Contenido del directorio"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Si abres el archivo VCF, deberías ver el mismo contenido que en el archivo que generaste ejecutando el comando GATK directamente en el contenedor.

??? abstract "Contenido del archivo"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Esta es la salida que nos importa generar para cada muestra en nuestro estudio.

### Conclusión

Sabes cómo hacer un workflow modular de dos pasos que hace trabajo de análisis real y es capaz de lidiar con las idiosincrasias de los formatos de archivo de genómica como los archivos accesorios.

### ¿Qué sigue?

Hacer que el workflow maneje múltiples muestras en lote.

---

## 3. Adaptar el workflow para ejecutarse en un lote de muestras

Está muy bien tener un workflow que pueda automatizar el procesamiento en una sola muestra, pero ¿qué pasa si tienes 1000 muestras?
¿Necesitas escribir un script bash que haga un ciclo a través de todas tus muestras?

¡No, gracias a Dios! Solo haz un pequeño ajuste al código y Nextflow manejará eso por ti también.

### 3.1. Actualizar la entrada para listar tres muestras

Para ejecutar en múltiples muestras, actualiza el perfil de prueba para proporcionar un arreglo de rutas de archivos en lugar de una sola.
Esta es una forma rápida de probar la ejecución multi-muestra; en el siguiente paso cambiaremos a un enfoque más escalable usando un archivo de entradas.

Primero, comenta la anotación de tipo en la declaración del parámetro, ya que los arreglos no pueden usar declaraciones tipadas:

=== "Después"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (array of three samples)
        reads_bam //: Path
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input
        reads_bam: Path
    ```

Luego actualiza el perfil de prueba para listar las tres muestras:

=== "Después"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

El channel factory en el cuerpo del workflow (`.fromPath`) acepta múltiples rutas de archivos tan bien como una sola, así que no se necesitan otros cambios.

### 3.2. Ejecutar el workflow

Intenta ejecutar el workflow ahora que la plomería está configurada para ejecutarse en las tres muestras de prueba.

```bash
nextflow run genomics.nf -profile test -resume
```

Cosa curiosa: esto _podría funcionar_, O _podría fallar_. Por ejemplo, aquí hay una ejecución que tuvo éxito:

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Si tu ejecución del workflow tuvo éxito, ejecútalo nuevamente hasta que obtengas un error como este:

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

Si miras la salida de error del comando GATK, habrá una línea como esta:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Bueno, eso es raro, considerando que indexamos explícitamente los archivos BAM en el primer paso del workflow. ¿Podría haber algo mal con la plomería?

### 3.3. Solucionar el problema

Inspeccionaremos los directorios work y usaremos el operador `view()` para averiguar qué salió mal.

#### 3.3.1. Verificar los directorios work para las llamadas relevantes

Echa un vistazo dentro del directorio work para la llamada al proceso `GATK_HAPLOTYPECALLER` fallida listada en la salida de la consola.

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

Presta particular atención a los nombres del archivo BAM y el índice BAM que están listados en este directorio: `reads_son.bam` y `reads_father.bam.bai`.

¿Qué demonios? Nextflow ha preparado un archivo de índice en el directorio work de esta llamada al proceso, pero es el incorrecto. ¿Cómo pudo haber sucedido esto?

#### 3.3.2. Usar el [operador view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) para inspeccionar contenidos de canales

Agrega estas dos líneas en el cuerpo del workflow antes de la llamada al proceso `GATK_HAPLOTYPECALLER` para ver el contenido del canal:

=== "Después"

    ```groovy title="genomics.nf" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // temporary diagnostics
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

=== "Antes"

    ```groovy title="genomics.nf"
        SAMTOOLS_INDEX(reads_ch)

        // Call variants from the indexed BAM file
        GATK_HAPLOTYPECALLER(
    ```

Luego ejecuta el comando del workflow nuevamente.

```bash
nextflow run genomics.nf -profile test
```

Una vez más, esto puede tener éxito o fallar. Aquí está lo que muestra la salida de las dos llamadas `.view()` para una ejecución fallida:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Las primeras tres líneas corresponden al canal de entrada y las segundas, al canal de salida.
Puedes ver que los archivos BAM y los archivos de índice para las tres muestras no están listados en el mismo orden.

!!! note "Nota"

    Cuando llamas un proceso de Nextflow en un canal que contiene múltiples elementos, Nextflow intentará paralelizar la ejecución tanto como sea posible, y recolectará las salidas en cualquier orden en que estén disponibles.
    La consecuencia es que las salidas correspondientes pueden recolectarse en un orden diferente al que se proporcionaron las entradas originales.

Como está escrito actualmente, nuestro script de workflow asume que los archivos de índice saldrán del paso de indexación listados en el mismo orden madre/padre/hijo que se dieron las entradas.
Pero eso no está garantizado que sea el caso, razón por la cual a veces (aunque no siempre) los archivos incorrectos se emparejan en el segundo paso.

Para arreglar esto, necesitamos asegurarnos de que los archivos BAM y sus archivos de índice viajen juntos a través de los canales.

!!! tip "Consejo"

    Las declaraciones `view()` en el código del workflow no hacen nada, así que no es un problema dejarlas.
    Sin embargo, desordenarán tu salida de consola, así que recomendamos eliminarlas cuando hayas terminado de solucionar el problema.

### 3.4. Actualizar el workflow para manejar los archivos de índice correctamente

La solución es agrupar cada archivo BAM con su índice en una tupla, luego actualizar el proceso downstream y la plomería del workflow para que coincidan.

#### 3.4.1. Cambiar la salida del módulo SAMTOOLS_INDEX a una tupla

La forma más simple de asegurar que un archivo BAM y su índice permanezcan estrechamente asociados es empaquetarlos juntos en una tupla que sale de la tarea de índice.

!!! note "Nota"

    Una **tupla** es una lista finita y ordenada de elementos que se usa comúnmente para devolver múltiples valores desde una función. Las tuplas son particularmente útiles para pasar múltiples entradas o salidas entre procesos mientras se preserva su asociación y orden.

Actualiza la salida en `modules/samtools_index.nf` para incluir el archivo BAM:

=== "Después"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Antes"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

De esta manera, cada archivo de índice estará estrechamente acoplado con su archivo BAM original, y la salida general del paso de indexación será un solo canal que contiene pares de archivos.

#### 3.4.2. Cambiar la entrada del módulo GATK_HAPLOTYPECALLER para aceptar una tupla

Ya que hemos cambiado la 'forma' de la salida del primer proceso, necesitamos actualizar la definición de entrada del segundo proceso para que coincida.

Actualiza `modules/gatk_haplotypecaller.nf`:

=== "Después"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

Ahora necesitamos actualizar el workflow para reflejar la nueva estructura de tupla en la llamada al proceso y los objetivos de publicación.

#### 3.4.3. Actualizar la llamada a GATK_HAPLOTYPECALLER en el workflow

Ya no necesitamos proporcionar el `reads_ch` original al proceso `GATK_HAPLOTYPECALLER`, ya que el archivo BAM ahora está empaquetado en la salida del canal por `SAMTOOLS_INDEX`.

Actualiza la llamada en `genomics.nf`:

=== "Después"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Finalmente, necesitamos actualizar los objetivos de publicación para reflejar la nueva estructura de salida.

#### 3.4.4. Actualizar el objetivo de publicación para la salida del BAM indexado

Ya que la salida de SAMTOOLS_INDEX ahora es una tupla que contiene tanto el archivo BAM como su índice, renombra el objetivo de publicación de `bam_index` a `indexed_bam` para reflejar mejor su contenido:

=== "Después"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Con estos cambios, el BAM y su índice están garantizados de viajar juntos, así que el emparejamiento siempre será correcto.

### 3.5. Ejecutar el workflow corregido

Ejecuta el workflow nuevamente para asegurarte de que esto funcionará de manera confiable en el futuro.

```bash
nextflow run genomics.nf -profile test
```

Esta vez (y cada vez) todo debería ejecutarse correctamente:

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

El directorio de resultados ahora contiene tanto archivos BAM como BAI para cada muestra (de la tupla), junto con las salidas VCF:

??? abstract "Contenido del directorio de resultados"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

Al agrupar archivos asociados en tuplas, aseguramos que los archivos correctos siempre viajen juntos a través del workflow.
El workflow ahora procesa cualquier número de muestras de manera confiable, pero listarlas individualmente en el config no es muy escalable.
En el siguiente paso, cambiaremos a leer entradas desde un archivo.

### Conclusión

Sabes cómo hacer que tu workflow se ejecute en múltiples muestras (independientemente).

### ¿Qué sigue?

Hacer más fácil manejar muestras en lote.

---

## 4. Hacer que el workflow acepte un archivo de texto que contenga un lote de archivos de entrada

Una forma muy común de proporcionar múltiples archivos de datos de entrada a un workflow es hacerlo con un archivo de texto que contiene las rutas de archivos.
Puede ser tan simple como un archivo de texto listando una ruta de archivo por línea y nada más, o el archivo puede contener metadata adicional, en cuyo caso a menudo se llama una hoja de muestras (samplesheet).

Aquí vamos a mostrarte cómo hacer el caso simple.

### 4.1. Examinar el archivo de texto proporcionado que lista las rutas de archivos de entrada

Ya hicimos un archivo de texto listando las rutas de archivos de entrada, llamado `sample_bams.txt`, que puedes encontrar en el directorio `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Como puedes ver, listamos una ruta de archivo por línea, y son rutas absolutas.

!!! note "Nota"

    Los archivos que estamos usando aquí están simplemente en el sistema de archivos local de tu GitHub Codespaces, pero también podríamos apuntar a archivos en almacenamiento en la nube.
    Si no estás usando el entorno de Codespaces proporcionado, es posible que necesites adaptar las rutas de archivos para que coincidan con tu configuración local.

### 4.2. Actualizar el parámetro y el perfil de prueba

Cambia el parámetro `reads_bam` para que apunte al archivo `sample_bams.txt` en lugar de listar muestras individuales.

Restaura la anotación de tipo en el bloque params (ya que es una sola ruta nuevamente):

=== "Después"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Primary input (file of input files, one per line)
        reads_bam: Path
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="10"
        // Primary input (array of three samples)
        reads_bam
    ```

Luego actualiza el perfil de prueba para que apunte al archivo de texto:

=== "Después"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.reads_bam = "${projectDir}/data/sample_bams.txt"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

La lista de archivos ya no vive en el código en absoluto, lo cual es un gran paso en la dirección correcta.

### 4.3. Actualizar el channel factory para leer líneas de un archivo

Actualmente, nuestro channel factory de entrada trata cualquier archivo que le demos como las entradas de datos que queremos alimentar al proceso de indexación.
Ya que ahora le estamos dando un archivo que lista rutas de archivos de entrada, necesitamos cambiar su comportamiento para analizar el archivo y tratar las rutas de archivos que contiene como las entradas de datos.

Podemos hacer esto usando el mismo patrón que usamos en la [Parte 2 de Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): aplicando el operador [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) para analizar el archivo, luego una operación `map` para seleccionar el primer campo de cada línea.

=== "Después"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Create input channel from a CSV file listing input file paths
        reads_ch = Channel.fromPath(params.reads_bam)
                .splitCsv()
                .map { line -> file(line[0]) }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="24"
        // Create input channel (single file via CLI parameter)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

Técnicamente podríamos hacer esto de manera más simple usando el operador [`.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext), ya que nuestro archivo de entrada actualmente solo contiene rutas de archivos.
Sin embargo, al usar el operador más versátil `splitCsv` (complementado con `map`), podemos hacer nuestro workflow a prueba de futuro en caso de que decidamos agregar metadata al archivo que contiene rutas de archivos.

!!! tip "Consejo"

    Si no estás seguro de entender lo que están haciendo los operadores aquí, esta es otra gran oportunidad para usar el operador `.view()` para ver cómo se ven los contenidos del canal antes y después de aplicarlos.

### 4.4. Ejecutar el workflow

Ejecuta el workflow una vez más. Esto debería producir el mismo resultado que antes, ¿verdad?

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

¡Sí! De hecho, Nextflow detecta correctamente que las llamadas a procesos son exactamente las mismas, y ni siquiera se molesta en volver a ejecutar todo, ya que estábamos ejecutando con `-resume`.

¡Y eso es todo! Nuestro simple workflow de llamado de variantes tiene todas las características básicas que queríamos.

### Conclusión

Sabes cómo hacer un workflow modular de múltiples pasos para indexar un archivo BAM y aplicar llamado de variantes por muestra usando GATK.

Más generalmente, has aprendido cómo usar componentes esenciales de Nextflow y lógica para construir un pipeline de genómica simple que hace trabajo real, teniendo en cuenta las idiosincrasias de los formatos de archivo de genómica y los requisitos de herramientas.

### ¿Qué sigue?

¡Celebra tu éxito y toma un descanso extra largo!

En la siguiente parte de este curso, aprenderás cómo transformar este simple workflow de llamado de variantes por muestra para aplicar llamado de variantes conjunto a los datos.
