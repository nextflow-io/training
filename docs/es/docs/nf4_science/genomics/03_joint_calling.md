# Parte 3: Llamado conjunto de variantes en una cohorte

En la Parte 2, construiste un pipeline de llamado de variantes por muestra que procesaba los datos de cada muestra de forma independiente.
Ahora vamos a extenderlo para implementar el llamado conjunto de variantes, como se cubrió en la [Parte 1](01_method.md).

## Asignación

En esta parte del curso, vamos a extender el workflow para hacer lo siguiente:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generar un archivo índice para cada archivo BAM de entrada usando Samtools
2. Ejecutar GATK HaplotypeCaller en cada archivo BAM de entrada para generar un GVCF de llamados de variantes genómicas por muestra
3. Recolectar todos los GVCFs y combinarlos en un almacén de datos GenomicsDB
4. Ejecutar genotipado conjunto en el almacén de datos GVCF combinado para producir un VCF a nivel de cohorte

Esta parte se construye directamente sobre el workflow producido por la Parte 2.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que has completado la [Parte 2: Llamado de variantes por muestra](./02_per_sample_variant_calling.md) y tienes un pipeline `genomics.nf` funcional.

    Si no completaste la Parte 2 o quieres comenzar de nuevo para esta parte, puedes usar la solución de la Parte 2 como punto de partida.
    Ejecuta estos comandos desde dentro del directorio `nf4-science/genomics/`:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Esto te proporciona un workflow completo de llamado de variantes por muestra.
    Puedes verificar que se ejecuta correctamente ejecutando el siguiente comando:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Plan de la lección

Hemos dividido esto en dos pasos:

1. **Modificar el paso de llamado de variantes por muestra para producir un GVCF.**
   Esto cubre la actualización de comandos y salidas del proceso.
2. **Agregar un paso de genotipado conjunto que combine y genotipe los GVCFs por muestra.**
   Esto introduce el operador `collect()`, closures de Groovy para construcción de líneas de comando y procesos con múltiples comandos.

!!! note

     Asegúrate de estar en el directorio de trabajo correcto:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modificar el paso de llamado de variantes por muestra para producir un GVCF

El pipeline de la Parte 2 produce archivos VCF, pero el llamado conjunto requiere archivos GVCF.
Necesitamos activar el modo de llamado de variantes GVCF y actualizar la extensión del archivo de salida.

Recuerda el comando de llamado de variantes GVCF de la [Parte 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Comparado con el comando base de HaplotypeCaller que envolvimos en la Parte 2, las diferencias son el parámetro `-ERC GVCF` y la extensión de salida `.g.vcf`.

### 1.1. Indicar a HaplotypeCaller que emita un GVCF y actualizar la extensión de salida

Abre el archivo de módulo `modules/gatk_haplotypecaller.nf` para hacer dos cambios:

- Agregar el parámetro `-ERC GVCF` al comando GATK HaplotypeCaller;
- Actualizar la ruta del archivo de salida para usar la extensión correspondiente `.g.vcf`, según la convención de GATK.

Asegúrate de agregar una barra invertida (`\`) al final de la línea anterior cuando agregues `-ERC GVCF`.

=== "Después"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

También necesitamos actualizar el bloque output para que coincida con la nueva extensión de archivo.
Dado que cambiamos la salida del comando de `.vcf` a `.g.vcf`, el bloque `output:` del proceso debe reflejar el mismo cambio.

### 1.2. Actualizar la extensión del archivo de salida en el bloque de salidas del proceso

=== "Después"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Antes"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

También necesitamos actualizar la configuración de publicación y salida del workflow para reflejar las nuevas salidas GVCF.

### 1.3. Actualizar los destinos de publicación para las nuevas salidas GVCF

Dado que ahora estamos produciendo GVCFs en lugar de VCFs, deberíamos actualizar la sección `publish:` del workflow para usar nombres más descriptivos.
También organizaremos los archivos GVCF en su propio subdirectorio para mayor claridad.

=== "Después"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Ahora actualiza el bloque output para que coincida.

### 1.4. Actualizar el bloque output para la nueva estructura de directorios

También necesitamos actualizar el bloque `output` para colocar los archivos GVCF en un subdirectorio `gvcf`.

=== "Después"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="53"
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

Con el módulo, los destinos de publicación y el bloque output actualizados, podemos probar los cambios.

### 1.5. Ejecutar el pipeline

Ejecuta el workflow para verificar que los cambios funcionen.

```bash
nextflow run genomics.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

La salida de Nextflow se ve igual que antes, pero los archivos `.g.vcf` y sus archivos índice ahora están organizados en subdirectorios.

??? abstract "Contenido del directorio (enlaces simbólicos acortados)"

    ```console
    results/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Si abres uno de los archivos GVCF y lo revisas, puedes verificar que GATK HaplotypeCaller produjo archivos GVCF según lo solicitado.

### Conclusión

Cuando cambias el nombre del archivo de salida de un comando de herramienta, el bloque `output:` del proceso y la configuración de publicación/salida deben actualizarse para coincidir.

### ¿Qué sigue?

Aprende a recolectar el contenido de un canal y pasarlo al siguiente proceso como una única entrada.

---

## 2. Agregar un paso de genotipado conjunto

Ahora necesitamos recolectar los GVCFs por muestra, combinarlos en un almacén de datos GenomicsDB y ejecutar el genotipado conjunto para producir un VCF a nivel de cohorte.
Como se cubrió en la [Parte 1](01_method.md), esta es una operación de dos herramientas: GenomicsDBImport combina los GVCFs, luego GenotypeGVCFs produce los llamados de variantes finales.
Envolveremos ambas herramientas en un único proceso llamado `GATK_JOINTGENOTYPING`.

Recuerda los dos comandos de la [Parte 1](01_method.md):

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

El primer comando toma los GVCFs por muestra y un archivo de intervalos, y produce un almacén de datos GenomicsDB.
El segundo toma ese almacén de datos, un genoma de referencia, y produce el VCF final a nivel de cohorte.
La URI del contenedor es la misma que para HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Configurar las entradas

El proceso de genotipado conjunto necesita dos tipos de entradas que aún no tenemos: un nombre arbitrario de cohorte y las salidas GVCF recolectadas de todas las muestras agrupadas juntas.

#### 2.1.1. Agregar un parámetro `cohort_name`

Necesitamos proporcionar un nombre arbitrario para la cohorte.
Más adelante en la serie de capacitación aprenderás cómo usar metadatos de muestras para este tipo de cosas, pero por ahora simplemente declaramos un parámetro CLI usando `params` y le damos un valor predeterminado por conveniencia.

=== "Después"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Reunir las salidas de HaplotypeCaller entre muestras

Si conectáramos el canal de salida de `GATK_HAPLOTYPECALLER` directamente al nuevo proceso, Nextflow llamaría al proceso en cada GVCF de muestra por separado.
Queremos agrupar los tres GVCFs (y sus archivos índice) para que Nextflow los entregue todos juntos a una única llamada de proceso.

Podemos hacer eso usando el operador de canal `collect()`.
Agrega las siguientes líneas al cuerpo del `workflow`, justo después de la llamada a GATK_HAPLOTYPECALLER:

=== "Después"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Antes"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Desglosando esto:

1. Tomamos el canal de salida de `GATK_HAPLOTYPECALLER` usando la propiedad `.out`.
2. Debido a que nombramos las salidas usando `emit:` en la sección 1, podemos seleccionar los GVCFs con `.vcf` y los archivos índice con `.idx`. Sin salidas nombradas, tendríamos que usar `.out[0]` y `.out[1]`.
3. El operador `collect()` agrupa todos los archivos en un único elemento, por lo que `all_gvcfs_ch` contiene los tres GVCFs juntos, y `all_idxs_ch` contiene los tres archivos índice juntos.

Podemos recolectar los GVCFs y sus archivos índice por separado (en lugar de mantenerlos juntos en tuplas) porque Nextflow preparará todos los archivos de entrada juntos para la ejecución, por lo que los archivos índice estarán presentes junto a los GVCFs.

!!! tip

    Puedes usar el operador `view()` para inspeccionar el contenido de los canales antes y después de aplicar operadores de canal.

### 2.2. Escribir el proceso de genotipado conjunto y llamarlo en el workflow

Siguiendo el mismo patrón que usamos en la Parte 2, escribiremos la definición del proceso en un archivo de módulo, lo importaremos en el workflow y lo llamaremos sobre las entradas que acabamos de preparar.

#### 2.2.1. Construir una cadena para dar a cada GVCF un argumento `-V`

Antes de comenzar a completar la definición del proceso, hay una cosa que resolver.
El comando GenomicsDBImport espera un argumento `-V` separado para cada archivo GVCF, así:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Si escribiéramos `-V ${all_gvcfs_ch}`, Nextflow simplemente concatenaría los nombres de archivo y esa parte del comando se vería así:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Pero necesitamos que la cadena se vea así:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Importante, necesitamos construir esta cadena dinámicamente a partir de los archivos que estén en el canal recolectado.
Nextflow (vía Groovy) proporciona una forma concisa de hacer esto:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Desglosando esto:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` itera sobre cada ruta de archivo y antepone `-V `, produciendo `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` los concatena con espacios: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. El resultado se asigna a una variable local `gvcfs_line` (definida con `def`), que podemos interpolar en la plantilla de comando.

Esta línea va dentro del bloque `script:` del proceso, antes de la plantilla de comando.
Puedes colocar código Groovy arbitrario entre `script:` y el `"""` de apertura de la plantilla de comando.

Entonces podrás referirte a toda esa cadena como `gvcfs_line` en el bloque `script:` del proceso.

#### 2.2.2. Completar el módulo para el proceso de genotipado conjunto

Ahora podemos abordar la escritura del proceso completo.

Abre `modules/gatk_jointgenotyping.nf` y examina el esquema de la definición del proceso.

Adelante, completa la definición del proceso usando la información proporcionada arriba, luego verifica tu trabajo contra la solución en la pestaña "Después" a continuación.

=== "Antes"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Combinar GVCFs en almacén de datos GenomicsDB y ejecutar genotipado conjunto para producir llamados a nivel de cohorte
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Después"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Combinar GVCFs en almacén de datos GenomicsDB y ejecutar genotipado conjunto para producir llamados a nivel de cohorte
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    }
    ```

Hay varias cosas que vale la pena destacar aquí.

Como antes, varias entradas están listadas aunque los comandos no las referencien directamente: `all_idxs`, `ref_index` y `ref_dict`.
Listarlas asegura que Nextflow prepare estos archivos en el directorio de trabajo junto con los archivos que sí aparecen en los comandos, que GATK espera encontrar según convenciones de nomenclatura.

La variable `gvcfs_line` usa el closure de Groovy descrito arriba para construir los argumentos `-V` para GenomicsDBImport.

Este proceso ejecuta dos comandos en serie, tal como lo harías en la terminal.
GenomicsDBImport combina los GVCFs por muestra en un almacén de datos, luego GenotypeGVCFs lee ese almacén de datos y produce el VCF final a nivel de cohorte.
El almacén de datos GenomicsDB (`${cohort_name}_gdb`) es un artefacto intermedio usado solo dentro del proceso; no aparece en el bloque output.

Una vez que hayas completado esto, el proceso está listo para usar.
Para usarlo en el workflow, necesitarás importar el módulo y agregar una llamada al proceso.

#### 2.2.3. Importar el módulo

Agrega la declaración de importación a `genomics.nf`, debajo de las declaraciones de importación existentes:

=== "Después"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Antes"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

El proceso ahora está disponible en el ámbito del workflow.

#### 2.2.4. Agregar la llamada al proceso

Agrega la llamada a `GATK_JOINTGENOTYPING` en el cuerpo del workflow, después de las líneas `collect()`:

=== "Después"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
        GATK_JOINTGENOTYPING(
            all_gvcfs_ch,
            all_idxs_ch,
            intervals_file,
            params.cohort_name,
            ref_file,
            ref_index_file,
            ref_dict_file
        )
    ```

=== "Antes"

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

El proceso ahora está completamente conectado.
A continuación, configuramos cómo se publican las salidas.

### 2.3. Configurar el manejo de salidas

Necesitamos publicar las salidas del VCF conjunto.
Agrega destinos de publicación y entradas de bloque output para los resultados de genotipado conjunto.

#### 2.3.1. Agregar destinos de publicación para el VCF conjunto

Agrega el VCF conjunto y su índice a la sección `publish:` del workflow:

=== "Después"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Antes"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Ahora actualiza el bloque output para que coincida.

#### 2.3.2. Agregar entradas de bloque output para el VCF conjunto

Agrega entradas para los archivos VCF conjunto.
Los colocaremos en la raíz del directorio de resultados ya que esta es la salida final.

=== "Después"

    ```groovy title="genomics.nf" hl_lines="11-16"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Antes"

    ```groovy title="genomics.nf"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

Con el proceso, los destinos de publicación y el bloque output todos en su lugar, podemos probar el workflow completo.

### 2.4. Ejecutar el workflow

Ejecuta el workflow para verificar que todo funcione.

```bash
nextflow run genomics.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Los primeros dos pasos están almacenados en caché de la ejecución anterior, y el nuevo paso `GATK_JOINTGENOTYPING` se ejecuta una vez sobre las entradas recolectadas de las tres muestras.
El archivo de salida final, `family_trio.joint.vcf` (y su índice), están en el directorio de resultados.

??? abstract "Contenido del directorio (enlaces simbólicos acortados)"

    ```console
    results/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Si abres el archivo VCF conjunto, puedes verificar que el workflow produjo los llamados de variantes esperados.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

¡Ahora tienes un workflow de llamado conjunto de variantes automatizado y completamente reproducible!

!!! note

    Ten en cuenta que los archivos de datos que te proporcionamos cubren solo una pequeña porción del cromosoma 20.
    El tamaño real de un conjunto de llamados de variantes se contaría en millones de variantes.
    ¡Por eso usamos solo pequeños subconjuntos de datos para fines de capacitación!

### Conclusión

Sabes cómo recolectar salidas de un canal y agruparlas como una única entrada para otro proceso.
También sabes cómo construir una línea de comando usando closures de Groovy, y cómo ejecutar múltiples comandos en un único proceso.

### ¿Qué sigue?

Celebra tu éxito y tómate un merecido descanso.

En la siguiente parte de este curso, aprenderás cómo ejecutar un pipeline de llamado de variantes listo para producción desde nf-core y compararlo con el pipeline que construiste manualmente.
