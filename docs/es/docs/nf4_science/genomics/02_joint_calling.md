# Parte 2: Llamado conjunto en una cohorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la primera parte de este curso, construyó un pipeline de llamado de variantes que era completamente lineal y procesaba los datos de cada muestra independientemente de las demás.
Sin embargo, en un caso de uso genómico real, típicamente necesitará examinar los llamados de variantes de múltiples muestras juntas.

En esta segunda parte, le mostramos cómo usar canales y operadores de canal para implementar el llamado conjunto de variantes con GATK, basándose en el pipeline de la Parte 1.

### Descripción general del método

El método de llamado de variantes de GATK que usamos en la primera parte de este curso simplemente generaba llamados de variantes por muestra.
Eso está bien si solo desea examinar las variantes de cada muestra en aislamiento, pero eso produce información limitada.
A menudo es más interesante observar cómo difieren los llamados de variantes entre múltiples muestras, y para hacerlo, GATK ofrece un método alternativo llamado llamado conjunto de variantes, que demostramos aquí.

El llamado conjunto de variantes implica generar un tipo especial de salida de variantes llamado GVCF (por Genomic VCF) para cada muestra, luego combinar los datos GVCF de todas las muestras y finalmente, ejecutar un análisis estadístico de 'genotipado conjunto'.

![Análisis conjunto](img/joint-calling.png)

Lo especial del GVCF de una muestra es que contiene registros que resumen las estadísticas de datos de secuencia sobre todas las posiciones en el área objetivo del genoma, no solo las posiciones donde el programa encontró evidencia de variación.
Esto es crítico para el cálculo del genotipado conjunto ([lectura adicional](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

El GVCF es producido por GATK HaplotypeCaller, la misma herramienta que usamos en la Parte 1, con un parámetro adicional (`-ERC GVCF`).
La combinación de los GVCFs se realiza con GATK GenomicsDBImport, que combina los llamados por muestra en un almacén de datos (análogo a una base de datos), luego el análisis de 'genotipado conjunto' propiamente dicho se realiza con GATK GenotypeGVCFs.

### Flujo de trabajo

Entonces, para recapitular, en esta parte del curso, vamos a desarrollar un flujo de trabajo que hace lo siguiente:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generar un archivo de índice para cada archivo BAM de entrada usando Samtools
2. Ejecutar GATK HaplotypeCaller en cada archivo BAM de entrada para generar un GVCF de llamados de variantes genómicas por muestra
3. Recopilar todos los GVCFs y combinarlos en un almacén de datos GenomicsDB
4. Ejecutar genotipado conjunto en el almacén de datos GVCF combinado para producir un VCF a nivel de cohorte

Aplicaremos esto al mismo conjunto de datos que en la Parte 1.

---

## 0. Calentamiento: Ejecutar Samtools y GATK directamente

Al igual que anteriormente, queremos probar los comandos manualmente antes de intentar envolverlos en un flujo de trabajo.

!!! note

     Asegúrese de estar en el directorio de trabajo correcto:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Indexar un archivo BAM de entrada con Samtools

Este primer paso es el mismo que en la Parte 1, por lo que debería sentirse muy familiar, pero esta vez necesitamos hacerlo para las tres muestras.

!!! note

    Técnicamente ya hemos generado archivos de índice para las tres muestras a través de nuestro pipeline, por lo que podríamos ir a buscarlos en el directorio de resultados. Sin embargo, es más limpio simplemente rehacerlo manualmente, y solo tomará un minuto.

#### 0.1.1. Iniciar el contenedor de Samtools de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

#### 0.1.2. Ejecutar el comando de indexación para las tres muestras

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Al igual que anteriormente, esto debería producir los archivos de índice en el mismo directorio que los archivos BAM correspondientes.

??? abstract "Contenido del directorio"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Ahora que tenemos archivos de índice para las tres muestras, podemos proceder a generar los GVCFs para cada una de ellas.

#### 0.1.3. Salir del contenedor de Samtools

```bash
exit
```

### 0.2. Llamar variantes con GATK HaplotypeCaller en modo GVCF

Este segundo paso es muy similar a lo que hicimos en la Parte 1: Hello Genomics, pero ahora vamos a ejecutar GATK en 'modo GVCF'.

#### 0.2.1. Iniciar el contenedor de GATK de forma interactiva

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

#### 0.2.2. Ejecutar el comando de llamado de variantes con la opción GVCF

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

<!--
??? success "Salida del comando"

    ```console

    ```
-->

Esto crea el archivo de salida GVCF `reads_mother.g.vcf` en el directorio de trabajo actual en el contenedor.

Si ejecuta `cat` para ver su contenido, verá que es mucho más largo que el VCF equivalente que generamos en la Parte 1. Ni siquiera puede desplazarse hacia arriba hasta el inicio del archivo, y la mayoría de las líneas se ven bastante diferentes de lo que vimos en el VCF en la Parte 1.

```console title="Salida" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Estas representan regiones no variantes donde el llamador de variantes no encontró evidencia de variación, por lo que capturó algunas estadísticas que describen su nivel de confianza en la ausencia de variación. Esto hace posible distinguir entre dos casos muy diferentes: (1) hay datos de buena calidad que muestran que la muestra es homocigota-referencia, y (2) no hay suficientes datos buenos disponibles para hacer una determinación de cualquier manera.

En un GVCF, típicamente hay muchas de estas líneas no variantes, con un número menor de registros de variantes intercalados entre ellas. Intente ejecutar `head -176` en el GVCF para cargar solo las primeras 176 líneas del archivo para encontrar un llamado de variante real.

```console title="Salida" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

La segunda línea muestra el primer registro de variante en el archivo, que corresponde a la primera variante en el archivo VCF que examinamos en la Parte 1.

Al igual que el VCF original, el archivo GVCF de salida también está acompañado por un archivo de índice, llamado `reads_mother.g.vcf.idx`.

#### 0.2.3. Repetir el proceso en las otras dos muestras

Para probar el paso de genotipado conjunto, necesitamos GVCFs para las tres muestras, así que generémoslos manualmente ahora.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

Una vez que esto se complete, debería tener tres archivos que terminan en `.g.vcf` en su directorio actual (uno por muestra) y sus respectivos archivos de índice que terminan en `.g.vcf.idx`.

### 0.3. Ejecutar genotipado conjunto

Ahora que tenemos todos los GVCFs, finalmente podemos probar el enfoque de genotipado conjunto para generar llamados de variantes para una cohorte de muestras.
Como recordatorio, es un método de dos pasos que consiste en combinar los datos de todos los GVCFs en un almacén de datos, luego ejecutar el análisis de genotipado conjunto propiamente dicho para generar el VCF final de variantes llamadas conjuntamente.

#### 0.3.1. Combinar todos los GVCFs por muestra

Este primer paso usa otra herramienta de GATK, llamada GenomicsDBImport, para combinar los datos de todos los GVCFs en un almacén de datos GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

La salida de este paso es efectivamente un directorio que contiene un conjunto de directorios más anidados que contienen los datos de variantes combinados en forma de múltiples archivos diferentes.
Puede explorar dentro de él pero rápidamente verá que este formato de almacén de datos no está destinado a ser leído directamente por humanos.

!!! note

    GATK incluye herramientas que hacen posible inspeccionar y extraer datos de llamados de variantes del almacén de datos según sea necesario.

#### 0.3.2. Ejecutar el análisis de genotipado conjunto propiamente dicho

Este segundo paso usa otra herramienta de GATK, llamada GenotypeGVCFs, para recalcular las estadísticas de variantes y los genotipos individuales a la luz de los datos disponibles en todas las muestras de la cohorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Salida del comando"

    ```console

    ```
-->

Esto crea el archivo de salida VCF `family_trio.vcf` en el directorio de trabajo actual en el contenedor.
Es otro archivo razonablemente pequeño, por lo que puede ejecutar `cat` en este archivo para ver su contenido, y desplazarse hacia arriba para encontrar las primeras líneas de variantes.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

Esto se parece más al VCF original que generamos en la Parte 1, excepto que esta vez tenemos información a nivel de genotipo para las tres muestras.
Las últimas tres columnas en el archivo son los bloques de genotipo para las muestras, listadas en orden alfabético.

Si observamos los genotipos llamados para nuestro trío familiar de prueba para la primera variante, vemos que el padre es heterocigoto-variante (`0/1`), y la madre y el hijo son ambos homocigotos-variante (`1/1`).

¡Esa es en última instancia la información que buscamos extraer del conjunto de datos! Así que envolvamos todo esto en un flujo de trabajo de Nextflow para poder hacerlo a escala.

#### 0.3.3. Salir del contenedor de GATK

```bash
exit
```

### Conclusión

Sabe cómo ejecutar los comandos individuales involucrados en el llamado conjunto de variantes en la terminal para verificar que producirán la información que desea.

### ¿Qué sigue?

Envolver estos comandos en un pipeline real.

---

## 1. Modificar el paso de llamado de variantes por muestra para producir un GVCF

La buena noticia es que no necesitamos comenzar desde cero, ya que escribimos un flujo de trabajo que hace parte de este trabajo en la Parte 1.
Sin embargo, ese pipeline produce archivos VCF, mientras que ahora queremos archivos GVCF para hacer el genotipado conjunto.
Entonces necesitamos comenzar activando el modo de llamado de variantes GVCF y actualizando la extensión del archivo de salida.

!!! note

    Por conveniencia, vamos a trabajar con una copia nueva del flujo de trabajo de GATK tal como está al final de la Parte 1, pero con un nombre diferente: `genomics-2.nf`.

### 1.1. Indicar a HaplotypeCaller que emita un GVCF y actualizar la extensión de salida

Abramos el archivo `genomics-2.nf` en el editor de código.
Debería verse muy familiar, pero siéntase libre de ejecutarlo si desea asegurarse de que funciona como se espera.

Vamos a comenzar haciendo dos cambios:

- Agregar el parámetro `-ERC GVCF` al comando GATK HaplotypeCaller;
- Actualizar la ruta del archivo de salida para usar la extensión `.g.vcf` correspondiente, según la convención de GATK.

Asegúrese de agregar una barra invertida (`\`) al final de la línea anterior cuando agregue `-ERC GVCF`.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
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

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Y eso es todo lo que se necesita para cambiar HaplotypeCaller a generar GVCFs en lugar de VCFs, ¿verdad?

### 1.2. Ejecutar el pipeline para verificar que puede generar GVCFs

El comando de ejecución de Nextflow es el mismo que antes, salvo por el nombre del archivo de flujo de trabajo en sí.
Asegúrese de actualizarlo apropiadamente.

```bash
nextflow run genomics-2.nf
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

Y la salida es... ¡todo rojo! Oh no.

El comando que se ejecutó es correcto, así que teníamos razón en que eso era suficiente para cambiar el comportamiento de la herramienta GATK.
Pero observe esa línea sobre el archivo de salida faltante. ¿Nota algo?

Así es, olvidamos decirle a Nextflow que espere un nuevo nombre de archivo. Ups.

### 1.3. Actualizar la extensión del archivo de salida en el bloque de salidas del proceso también

Porque no es suficiente solo cambiar la extensión del archivo en el comando de la herramienta en sí, también tiene que decirle a Nextflow que el nombre del archivo de salida esperado ha cambiado.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Actualizar los destinos de publicación para las nuevas salidas GVCF

Ya que ahora estamos produciendo GVCFs en lugar de VCFs, deberíamos actualizar la sección `publish:` del flujo de trabajo para usar nombres más descriptivos.
También organizaremos los archivos GVCF en su propio subdirectorio para mayor claridad.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Actualizar el bloque de salida para la nueva estructura de directorios

También necesitamos actualizar el bloque `output` para colocar los archivos GVCF en un subdirectorio `gvcf`.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
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

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
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

### 1.6. Ejecutar el pipeline nuevamente

Ejecutémoslo con `-resume` esta vez.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Esta vez funciona.

La salida de Nextflow en sí no se ve diferente (comparada con una ejecución exitosa en modo VCF normal), pero ahora podemos encontrar los archivos `.g.vcf` y sus respectivos archivos de índice, para las tres muestras, organizados en subdirectorios.

??? abstract "Contenido del directorio (symlinks acortados)"

    ```console
    results_genomics/
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

Si abre uno de los archivos GVCF y se desplaza a través de él, puede verificar que GATK HaplotypeCaller produjo archivos GVCF como se solicitó.

### Conclusión

Bien, este fue mínimo en términos de aprendizaje de Nextflow...
¡Pero fue una buena oportunidad para reiterar la importancia del bloque de salida del proceso!

### ¿Qué sigue?

Aprender a recopilar el contenido de un canal y pasarlos al siguiente proceso como una única entrada.

---

## 2. Recopilar y combinar los datos GVCF de todas las muestras

Ahora necesitamos combinar los datos de todos los GVCFs por muestra en una forma que soporte el análisis de genotipado conjunto que queremos hacer.

### 2.1. Definir el proceso que combinará los GVCFs

Como recordatorio de lo que hicimos anteriormente en la sección de calentamiento, combinar los GVCFs es un trabajo para la herramienta GATK GenomicsDBImport, que producirá un almacén de datos en el llamado formato GenomicsDB.

Escribamos un nuevo proceso para definir cómo va a funcionar eso, basándonos en el comando que usamos anteriormente en la sección de calentamiento.

```groovy title="genomics-2.nf" linenums="66"
/*
 * Combinar GVCFs en almacén de datos GenomicsDB
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

¿Qué piensa, se ve razonable?

Conectémoslo y veamos qué sucede.

### 2.2. Agregar un parámetro `cohort_name` con un valor predeterminado

Necesitamos proporcionar un nombre arbitrario para la cohorte.
Más adelante en la serie de entrenamiento aprenderá cómo usar metadatos de muestras para este tipo de cosas, pero por ahora solo declaramos un parámetro CLI usando `params` y le damos un valor predeterminado por conveniencia.

```groovy title="genomics-2.nf" linenums="16"
    // Nombre base para el archivo de salida final
    cohort_name: String = "family_trio"
```

### 2.3. Reunir las salidas de GATK_HAPLOTYPECALLER entre muestras

Si simplemente conectáramos el canal de salida del proceso `GATK_HAPLOTYPECALLER` tal como está, Nextflow llamaría al proceso en cada GVCF de muestra por separado.
Sin embargo, queremos agrupar los tres GVCFs (y sus archivos de índice) de tal manera que Nextflow los entregue todos juntos a una sola llamada de proceso.

Buenas noticias: podemos hacer eso usando el operador de canal `collect()`. Agreguemos las siguientes líneas al cuerpo del `workflow`, justo después de la llamada a GATK_HAPLOTYPECALLER:

```groovy title="genomics-2.nf" linenums="118"
// Recopilar salidas de llamado de variantes entre muestras
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

¿Parece un poco complicado? Desglosemos esto y traduzcámoslo a lenguaje sencillo.

1. Estamos tomando el canal de salida del proceso `GATK_HAPLOTYPECALLER`, referido usando la propiedad `.out`.
2. Cada 'elemento' que sale del canal es un par de archivos: el GVCF y su archivo de índice, en ese orden porque ese es el orden en que están listados en el bloque de salida del proceso. Convenientemente, debido a que en la última sesión nombramos las salidas de este proceso (usando `emit:`), podemos seleccionar los GVCFs por un lado agregando `.vcf` y los archivos de índice por otro agregando `.idx` después de la propiedad `.out`. Si no hubiéramos nombrado esas salidas, habríamos tenido que referirnos a ellas como `.out[0]` y `.out[1]`, respectivamente.
3. Agregamos el operador de canal `collect()` para agrupar todos los archivos GVCF juntos en un solo elemento en un nuevo canal llamado `all_gvcfs_ch`, y hacemos lo mismo con los archivos de índice para formar el nuevo canal llamado `all_idxs_ch`.

!!! tip

    Si tiene dificultades para visualizar exactamente qué está sucediendo aquí, recuerde que puede usar el operador `view()` para inspeccionar el contenido de los canales antes y después de aplicar operadores de canal.

Los canales `all_gvcfs_ch` y `all_idxs_ch` resultantes son lo que vamos a conectar en el proceso `GATK_GENOMICSDB` que acabamos de escribir.

!!! note

    En caso de que se lo estuviera preguntando, recopilamos los GVCFs y sus archivos de índice por separado porque el comando GATK GenomicsDBImport solo quiere ver las rutas de los archivos GVCF. Afortunadamente, dado que Nextflow preparará todos los archivos juntos para la ejecución, no tenemos que preocuparnos por el orden de los archivos como lo hicimos para los BAMs y su índice en la Parte 1.

### 2.4. Agregar una llamada al bloque de flujo de trabajo para ejecutar GATK_GENOMICSDB

Tenemos un proceso y tenemos canales de entrada. Solo necesitamos agregar la llamada al proceso.

```groovy title="genomics-2.nf" linenums="122"
    // Combinar GVCFs en un almacén de datos GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, todo está conectado.

### 2.5. Ejecutar el flujo de trabajo

Veamos si esto funciona.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

Se ejecuta bastante rápido, ya que estamos ejecutando con `-resume`, ¡pero falla!

Ah. Por el lado positivo, vemos que Nextflow ha recogido el proceso `GATK_GENOMICSDB`, y específicamente lo llamó solo una vez.
Eso sugiere que el enfoque de `collect()` funcionó, hasta cierto punto.
Pero, y es grande, la llamada al proceso falló.

Cuando profundizamos en la salida de la consola anterior, podemos ver que el comando ejecutado no es correcto.

¿Puede detectar el error?
Observe este fragmento: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

Le dimos a `gatk GenomicsDBImport` múltiples archivos GVCF para un solo argumento `-V`, pero la herramienta espera un argumento `-V` separado para cada archivo GVCF.

Como recordatorio, este fue el comando que ejecutamos en el contenedor:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Entonces eso significa que necesitamos de alguna manera transformar nuestro paquete de archivos GVCF en una cadena de comando formateada correctamente.

### 2.6. Construir una línea de comando con un argumento `-V` separado para cada GVCF de entrada

Aquí es donde Nextflow estar basado en Groovy resulta útil, porque nos permitirá usar algunas manipulaciones de cadenas bastante directas para construir la cadena de comando necesaria.

Específicamente, usando esta sintaxis: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Una vez más, desglosémoslo en sus componentes.

1. Primero, tomamos el contenido del canal de entrada `all_gvcfs` y aplicamos `.collect()` en él (al igual que antes).
2. Eso nos permite pasar cada ruta de archivo GVCF individual en el paquete a la **closure**, `{ gvcf -> "-V ${gvcf}" }`, donde `gvcf` se refiere a esa ruta de archivo GVCF.
   La closure es una mini-función que usamos para anteponer `-V ` a la ruta del archivo, en la forma de `"-V ${gvcf}"`.
3. Luego usamos `.join(' ')` para concatenar las tres cadenas con un solo espacio como separador.

Con un ejemplo concreto, se ve así:

1. Tenemos tres archivos:

   `[A.ext, B.ext, C.ext]`

2. La closure modifica cada uno para crear las cadenas:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. La operación `.join(' ')` genera la cadena final:

   `"-V A.ext -V B.ext -V C.ext"`

Una vez que tengamos esa cadena, podemos asignarla a una variable local, `gvcfs_line`, definida con la palabra clave `def`:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, entonces tenemos nuestra cosa de manipulación de cadenas. ¿Dónde la ponemos?

Queremos que esto vaya dentro de la definición del proceso en alguna parte, porque queremos hacerlo _después_ de haber canalizado las rutas de archivos GVCF en el proceso.
Eso es porque Nextflow debe verlas como rutas de archivos para preparar los archivos en sí correctamente para la ejecución.

¿Pero _dónde_ en el proceso podemos agregar esto?

Dato curioso: ¡puede agregar código arbitrario después de `script:` y antes de las `"""`!

Genial, agreguemos nuestra línea de manipulación de cadenas allí entonces, y actualicemos el comando `gatk GenomicsDBImport` para usar la cadena concatenada que produce.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Eso debería ser todo lo necesario para proporcionar las entradas a `gatk GenomicsDBImport` correctamente.

!!! tip

    Cuando actualice el comando `gatk GenomicsDBImport`, asegúrese de eliminar el prefijo `-V ` cuando intercambie la variable `${gvcfs_line}`.

### 2.7. Ejecutar el flujo de trabajo para verificar que genera la salida GenomicsDB como se esperaba

Muy bien, veamos si eso abordó el problema.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

¡Ajá! Parece estar funcionando ahora.

Los primeros dos pasos se omitieron exitosamente, y el tercer paso funcionó como un encanto esta vez.
El almacén de datos GenomicsDB se crea en el directorio de trabajo pero no se publica en los resultados, ya que es solo un formato intermedio que usaremos para el genotipado conjunto.

Por cierto, no tuvimos que hacer nada especial para manejar que la salida sea un directorio en lugar de un solo archivo.

### Conclusión

Ahora sabe cómo recopilar salidas de un canal y agruparlas como una única entrada a otro proceso.
También sabe cómo construir una línea de comando para proporcionar entradas a una herramienta dada con la sintaxis apropiada.

### ¿Qué sigue?

Aprender cómo agregar un segundo comando al mismo proceso.

---

## 3. Ejecutar el paso de genotipado conjunto como parte del mismo proceso

Ahora que tenemos los llamados de variantes genómicas combinados, podemos ejecutar la herramienta de genotipado conjunto, que producirá la salida final que realmente nos interesa: el VCF de llamados de variantes a nivel de cohorte.

Por razones logísticas, decidimos incluir el genotipado conjunto dentro del mismo proceso.

### 3.1. Renombrar el proceso de GATK_GENOMICSDB a GATK_JOINTGENOTYPING

Dado que el proceso ejecutará más de una herramienta, cambiamos su nombre para referirse a la operación general en lugar de un solo nombre de herramienta.

=== "Después"

    ```groovy title="genomics-2.nf"
    /*
     * Combinar GVCFs en almacén de datos GenomicsDB y ejecutar genotipado conjunto para producir llamados a nivel de cohorte
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Antes"

    ```groovy title="genomics-2.nf"
    /*
     * Combinar GVCFs en almacén de datos GenomicsDB
     */
    process GATK_GENOMICSDB {
    ```

¡Recuerde mantener sus nombres de procesos lo más descriptivos posible, para maximizar la legibilidad para sus colegas —y su yo futuro!

### 3.2. Agregar el comando de genotipado conjunto al proceso GATK_JOINTGENOTYPING

Simplemente agregue el segundo comando después del primero dentro de la sección de script.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
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
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Los dos comandos se ejecutarán en serie, de la misma manera que lo harían si los ejecutáramos manualmente en la terminal.

### 3.3. Agregar los archivos del genoma de referencia a las definiciones de entrada del proceso GATK_JOINTGENOTYPING

El segundo comando requiere los archivos del genoma de referencia, así que necesitamos agregarlos a las entradas del proceso.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Puede parecer molesto escribirlos, pero recuerde, solo los escribe una vez, y luego puede ejecutar el flujo de trabajo un millón de veces. ¿Vale la pena?

### 3.4. Actualizar la definición de salida del proceso para emitir el VCF de llamados de variantes a nivel de cohorte

Realmente no nos importa guardar el almacén de datos GenomicsDB, que es solo un formato intermedio que solo existe por razones logísticas, así que podemos simplemente eliminarlo del bloque de salida si queremos.

La salida que realmente nos interesa es el VCF producido por el comando de genotipado conjunto.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Antes"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

¡Casi terminamos!

### 3.5. Actualizar la llamada al proceso de GATK_GENOMICSDB a GATK_JOINTGENOTYPING

No olvidemos renombrar la llamada al proceso en el cuerpo del flujo de trabajo de GATK_GENOMICSDB a GATK_JOINTGENOTYPING. Y mientras estamos en eso, también deberíamos agregar los archivos del genoma de referencia como entradas, ya que necesitamos proporcionarlos a la herramienta de genotipado conjunto.

=== "Después"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combinar GVCFs en un almacén de datos GenomicsDB y aplicar genotipado conjunto
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

    ```groovy title="genomics-2.nf" linenums="126"
    // Combinar GVCFs en un almacén de datos GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Ahora el proceso está completamente conectado.

### 3.6. Agregar el VCF conjunto a la sección de publicación

Necesitamos publicar las salidas del VCF conjunto del nuevo proceso.
Agregue estas líneas a la sección `publish:` del flujo de trabajo:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Agregar los destinos del VCF conjunto al bloque de salida

Finalmente, agregue destinos de salida para los archivos VCF conjuntos.
Los colocaremos en la raíz del directorio de resultados ya que esta es la salida final.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Ahora todo debería estar completamente conectado.

### 3.8. Ejecutar el flujo de trabajo

Finalmente, podemos ejecutar el flujo de trabajo modificado...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

¡Y funciona!

Encontrará el archivo de salida final, `family_trio.joint.vcf` (y su índice de archivo), en el directorio de resultados.

??? abstract "Contenido del directorio (symlinks acortados)"

    ```console
    results_genomics/
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

Si es del tipo escéptico, puede hacer clic en el archivo VCF conjunto para abrirlo y verificar que el flujo de trabajo haya generado los mismos llamados de variantes que obtuvo al ejecutar las herramientas manualmente al comienzo de esta sección.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

¡Ahora tiene un flujo de trabajo de llamado conjunto de variantes automatizado y completamente reproducible!

!!! note

    Tenga en cuenta que los archivos de datos que le dimos cubren solo una pequeña porción del cromosoma 20.
    El tamaño real de un conjunto de llamados de variantes se contaría en millones de variantes.
    ¡Por eso usamos solo subconjuntos pequeños de datos para propósitos de entrenamiento!

### Conclusión

Sabe cómo usar algunos operadores comunes así como closures de Groovy para controlar el flujo de datos en su flujo de trabajo.

### ¿Qué sigue?

Celebre su éxito y tome un merecido descanso.

En la próxima parte de este curso, aprenderá cómo modularizar su flujo de trabajo extrayendo definiciones de procesos en módulos reutilizables.
