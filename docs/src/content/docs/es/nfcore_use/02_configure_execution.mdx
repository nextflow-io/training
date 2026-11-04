# Parte 2: Configurar la ejecución del pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la [Parte 1](./01_run_demo.md), encontró y ejecutó el pipeline nf-core/demo usando su perfil de prueba.
Ahora veremos cómo configurar la ejecución del pipeline: establecer parámetros, comprender la validación y personalizar la asignación de recursos y los argumentos de las herramientas.

Como se explica en [Hello Config](../hello_nextflow/06_hello_config.md), queremos poder cambiar con qué datos se ejecutará nuestro pipeline y cómo se ejecutará sin modificar el código del pipeline en sí.
Para ello, Nextflow admite múltiples formas de controlar la configuración del pipeline, lo que puede resultar un poco abrumador.

El proyecto nf-core especifica convenciones para organizar los elementos de configuración, distinguiendo dos tipos de configuración en el nivel superior: **parámetros del pipeline** y **configuración** en sentido estricto.

- **Parámetros del pipeline** (establecidos a través del sistema `params`) generalmente incluyen cosas como archivos de entrada, indicadores de comportamiento de herramientas y parámetros de análisis.
- **Configuración** en sentido estricto se refiere a la logística de cómo se ejecuta el pipeline, es decir, el executor, la asignación de recursos de cómputo, etc.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Comencemos abordando los parámetros del pipeline y luego veremos la configuración en sentido estricto.

---

## 1. Parámetros del pipeline

Para todos los pipelines de nf-core, puede obtener una lista completa de los parámetros del pipeline directamente desde la línea de comandos usando el indicador `--help`, que en sí mismo es un parámetro del pipeline.

### 1.1. Obtener la lista de parámetros con `--help`

Ejecute el comando de ayuda para el pipeline demo:

```bash
nextflow run nf-core/demo --help
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Como puede ver, la salida agrupa los parámetros en categorías (opciones de entrada/salida, opciones de genoma de referencia, etc.) con tipos y descripciones para cada uno.

Esta categorización está determinada por un archivo de esquema, que se explica más adelante.
En pipelines de Nextflow simples, `--help` solo funciona si el desarrollador lo implementó manualmente.

!!! tip "Consejo"

    Use `--help --show_hidden` para ver parámetros adicionales que están ocultos por defecto, como `--publish_dir_mode` o `--monochrome_logs`.

### 1.2. Establecer valores de parámetros

Como se explica en [Hello Config](../hello_nextflow/06_hello_config.md), puede establecer valores de parámetros en la línea de comandos con `--nombre_param` o recopilar un conjunto de parámetros en un archivo YAML y pasarlo con `-params-file`.
Ambos enfoques funcionan de la misma manera con los pipelines de nf-core.

Por ejemplo, para omitir el paso de recorte, queremos establecer el parámetro booleano `skip_trim` en `true`.
En su directorio de trabajo se proporciona un archivo de parámetros llamado `my_params.yml` con ese valor ya establecido:

```yaml title="my_params.yml"
skip_trim: true
```

Páselo con `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El proceso `SEQTK_TRIM` ya no aparece en la salida.

!!! warning "Advertencia: Limitaciones importantes sobre las entradas de parámetros"

    **Establecer parámetros booleanos en la línea de comandos**

    A partir de la versión 26.04 de Nextflow, todos los valores proporcionados en la línea de comandos se tipifican como strings.
    Para un parámetro booleano como `skip_trim`, pasarlo como un indicador simple (`--skip_trim`) o como `--skip_trim true` se evalúa como el **string** `"true"`, lo que falla en la validación del esquema:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Para establecer un parámetro booleano con un valor genuino `true`/`false`, use un `-params-file` como se muestra arriba, o establézcalo en un archivo de configuración.
    Los parámetros de tipo string, integer y ruta de archivo no se ven afectados y aún pueden establecerse directamente en la línea de comandos.
    Este curso usa este patrón en todo momento para los parámetros booleanos.

    **Usar archivos de configuración personalizados**

    Aunque técnicamente es posible establecer parámetros del pipeline en un archivo de configuración personalizado pasado con `-c`, esto puede no anular los valores predeterminados ya establecidos en el propio `nextflow.config` del pipeline, dependiendo de las reglas de precedencia de configuración de Nextflow.
    Usar `--nombre_param` en la línea de comandos o `-params-file` es más confiable, ya que estos siempre tienen precedencia.

    Como regla general: si aparece en la salida de `--help`, establézcalo mediante la línea de comandos o un archivo de parámetros en lugar de un archivo de configuración.

### 1.3. Validación de parámetros

Dato curioso: el comando `--help` funciona para todos los pipelines de nf-core porque el proyecto nf-core requiere que los desarrolladores definan formalmente todos los parámetros del pipeline en un archivo de esquema JSON (`nextflow_schema.json`).
Este esquema registra el tipo, la descripción, el valor predeterminado y la agrupación de cada parámetro.

Además de potenciar la salida de `--help`, el archivo de esquema también permite la validación automatizada en el momento del lanzamiento.
Esto significa que Nextflow puede verificar que cada parámetro que pase exista y haya recibido un valor apropiado (del tipo apropiado, dentro del rango de valores permitidos, etc.).

Cubrimos esto con más detalle en la [sección de validación de entradas](../nfcore_build/04_input_validation.md), pero ya puede verlo en acción proporcionando al pipeline demo alguna entrada de parámetros no válida.

#### 1.3.1. Parámetros no reconocidos

Intente pasar un parámetro que no existe:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

La salida de la consola incluye una advertencia:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

El pipeline sigue ejecutándose, pero la advertencia le alerta de inmediato que `--foobar` no es un parámetro reconocido.
Esto está diseñado para llamar su atención sobre errores tipográficos no críticos, como usar `--outDir` en lugar de `--outdir`, lo que puede ayudarle a evitar desperdiciar tiempo y recursos de cómputo.

#### 1.3.2. Valores de parámetros no válidos

La validación también verifica los **valores** de los parámetros.
El parámetro `--skip_trim` es un indicador booleano, por lo que pasar un valor de tipo string hace que el pipeline falle inmediatamente:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

El pipeline se detiene antes de que se ejecute cualquier proceso, evitando una ejecución fallida o incorrecta.
Como se señala en [1.2](#12-set-parameter-values), los parámetros booleanos deben establecerse con un valor genuino `true`/`false` en un archivo de parámetros en lugar de pasarse en la línea de comandos, ya que los valores de la línea de comandos se tipifican como strings.

### 1.4. Validación de entradas

La misma lógica de validación también puede usarse para verificar la validez de los archivos de entrada.
Por ejemplo, si un pipeline espera una hoja de muestras como su entrada de datos principal (que es el caso de muchos, si no la mayoría, de los pipelines de nf-core), el desarrollador puede proporcionar un esquema de entrada (distinto del esquema de parámetros) que describa cómo debe estar estructurado el archivo de entrada.

Luego, en tiempo de ejecución, Nextflow puede verificar que el archivo de entrada proporcionado sea válido.

También cubrimos esto con más detalle en la [sección de validación de entradas](../nfcore_build/04_input_validation.md), pero ya puede verlo en acción proporcionando al pipeline demo una hoja de muestras de entrada no válida.

El pipeline `nf-core/demo` espera un archivo CSV con las columnas `sample`, `fastq_1` y `fastq_2`.
Esto está definido en un archivo de esquema (`assets/schema_input.json`) que especifica la estructura esperada, los tipos de columnas y las restricciones.

??? abstract "Archivo de esquema para entradas"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

El esquema especifica que `sample` y `fastq_1` son obligatorios, mientras que `fastq_2` es opcional (admitiendo datos tanto de extremos pareados como de extremo único).
Las rutas de archivo se validan en cuanto a su existencia y patrón de extensión.

Para demostrar esto, proporcionamos una hoja de muestras malformada llamada `malformed_samplesheet.csv` en su directorio de trabajo:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Esta hoja de muestras no tiene la columna obligatoria `fastq_1` y tiene una ruta de archivo inexistente en `fastq_2`.

Ejecute el pipeline demo usando `malformed_samplesheet.csv` como entrada:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Como puede ver, el pipeline falla inmediatamente e informa **todos** los errores de validación a la vez.
nf-schema no se detiene en el primer error — recopila cada problema y los lista juntos, para que pueda corregir todo de una vez en lugar de descubrir los problemas uno por uno.

Cada error identifica la entrada y el campo exactos que causaron el problema, para que pueda corregir su hoja de muestras y luego volver a lanzar el pipeline con la confianza de que no va a fallar en algún punto posterior cuando Nextflow intente acceder a la ruta del archivo.

Para los desarrolladores, todo esto se cubre con más detalle en la [Parte 4 de Build with nf-core](../nfcore_build/04_input_validation.md).

### Conclusión

Sabe cómo obtener una lista completa de los parámetros de un pipeline con `--help`, establecerlos mediante la línea de comandos o un archivo de parámetros, y cómo Nextflow valida tanto los valores de los parámetros como los archivos de entrada contra los esquemas del pipeline.

### ¿Qué sigue?

Aprenda sobre el otro tipo de configuración: cómo se ejecuta el pipeline, cubriendo la asignación de recursos y los argumentos de las herramientas.

---

## 2. Configuración

La configuración en sentido estricto controla **cómo** se ejecuta el pipeline: asignación de recursos, argumentos específicos de las herramientas, dónde se ejecutan las tareas y qué sistema de empaquetado de software se usa.

Los pipelines de nf-core incluyen configuración predeterminada en `nextflow.config` y el directorio `conf/`.
Antes de anular cualquier configuración, es útil saber dónde se encuentran los valores predeterminados.

### 2.1. Explorar los archivos de configuración

Ya vio en la [Parte 1](./01_run_demo.md) que el código fuente del pipeline se encuentra en `$NXF_HOME/assets`.
Usando el enlace simbólico `pipelines` que creó en la [Parte 1](./01_run_demo.md), liste los archivos de configuración para ver qué está disponible:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

Los archivos de configuración más importantes son:

- **`conf/base.config`**: Define etiquetas de recursos (`process_low`, `process_medium`, `process_high`) que asignan CPUs, memoria y tiempo a los procesos. Cuando vea que un proceso usa más recursos de los esperados, aquí es donde se originan esos valores predeterminados.
- **`conf/modules.config`**: Establece los argumentos de las herramientas por proceso (`ext.args`) y la configuración de publicación de salidas (`publishDir`). Abra este archivo para ver qué argumentos recibe cada herramienta por defecto.
- **`conf/test.config`**: El perfil de prueba que usó en la [Parte 1](./01_run_demo.md), que limita los recursos mediante `resourceLimits` y establece una hoja de muestras de prueba. Se activa con `-profile test`.
  También existe un `conf/test_full.config` para ejecutar con un conjunto de datos de prueba de tamaño completo, útil para benchmarking.

El `nextflow.config` central carga todo lo anterior y establece los valores predeterminados apropiados para todo.

Si desea modificar alguna de las configuraciones especificadas en estos archivos, no modifique ninguno de ellos directamente.
En su lugar, cree su propio archivo de configuración y páselo con `-c`.
Los valores que especifique anularán los valores predeterminados establecidos en esos otros archivos.

Probemos esto en la práctica.

### 2.2. Personalizar los recursos de los procesos y los argumentos de las herramientas

Los módulos de nf-core admiten dos tipos comunes de anulación de configuración: **asignación de recursos** (CPUs, memoria, tiempo) y **argumentos de herramientas** mediante `ext.args`.

Muchas herramientas de línea de comandos tienen argumentos que no se usan con suficiente frecuencia como para exponerlos como parámetros del pipeline.
La convención `ext.args` le permite pasar estos argumentos a la herramienta subyacente a través de un archivo de configuración.

El archivo `custom.config` proporcionado en su directorio de trabajo demuestra ambas anulaciones:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

El primer bloque anula la asignación de recursos de `FASTQC`.
Por defecto, `FASTQC` usa la etiqueta `process_medium` de `base.config`, que asigna 6 CPUs y 36 GB de memoria; aquí lo limitamos a 2 CPUs y 4 GB.

El segundo bloque pasa un argumento adicional a `SEQTK_TRIM` mediante `ext.args`.
El indicador `-b 5` le dice a `seqtk trimfq` que recorte 5 bases del inicio de cada lectura además del recorte por calidad.

Ejecute el pipeline con esta configuración:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Salida del comando"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El indicador `-c` agrega su configuración sobre la configuración integrada del pipeline.

Para verificar que la anulación de `ext.args` tuvo efecto, encuentre el hash del directorio de trabajo de `SEQTK_TRIM` en la salida de la ejecución (p. ej., `work/17/428668...`) y verifique el archivo `.command.sh` dentro de él:

```bash
cat work/17/428668/.command.sh
```

??? success "Salida del comando"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Debería ver `-b 5` en el comando `seqtk trimfq`.

Una cosa importante que debe saber sobre `ext.args`: si un módulo ya tiene un valor predeterminado establecido, su valor lo **reemplazará completamente** en lugar de agregarse a él.
Por ejemplo, `FASTQC` tiene `ext.args = '--quiet'` establecido por defecto en `conf/modules.config`:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Si establece `ext.args = '--kmers 8'` para `FASTQC`, el indicador `--quiet` ya no se aplicará.
Para mantener ambos, establezca `ext.args = '--quiet --kmers 8'`.

Siempre debe verificar la configuración predeterminada de un módulo antes de anular `ext.args`.

### Conclusión

Sabe dónde se encuentran los valores predeterminados de configuración de los pipelines de nf-core y cómo anular las asignaciones de recursos y los argumentos de las herramientas con un archivo de configuración personalizado.

### ¿Qué sigue?

Continúe con la [Parte 3](./03_run_production_pipeline.md), donde aplicará lo que ha aprendido a un pipeline de producción real.

---

## Resumen

En esta parte aprendió a:

- Obtener ayuda, establecer parámetros y comprender la validación de parámetros y entradas
- Personalizar la asignación de recursos y los argumentos de las herramientas mediante archivos de configuración
