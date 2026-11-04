# Parte 1: Ejecutar un pipeline de demostración

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta primera parte del curso de capacitación "Use nf-core", le mostramos cómo encontrar un pipeline de nf-core y probarlo usando su perfil de prueba integrado.

Vamos a usar un pipeline llamado nf-core/demo que es mantenido por el proyecto nf-core como parte de su inventario de pipelines para fines de demostración y capacitación.

Asegúrese de que su directorio de trabajo esté configurado en `nfcore-use/` como se indica en la página de [Primeros pasos](./00_orientation.md).

---

## 1. Encontrar y recuperar el pipeline nf-core/demo

Comencemos localizando el pipeline nf-core/demo en el sitio web del proyecto en [nf-co.re](https://nf-co.re), que centraliza toda la información, como: documentación general y artículos de ayuda, documentación para cada uno de los pipelines, publicaciones de blog, anuncios de eventos, entre otros.

### 1.1. Encontrar el pipeline en el sitio web

En su navegador web, vaya a [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) y escriba `demo` en la barra de búsqueda.

![resultados de búsqueda](./img/search-results.png)

Haga clic en el nombre del pipeline, `demo`, para acceder a la página de documentación del pipeline.

Cada pipeline publicado tiene una página dedicada que incluye las siguientes secciones de documentación:

- **Introduction:** Una introducción y descripción general del pipeline
- **Usage:** Descripciones de cómo ejecutar el pipeline
- **Parameters:** Parámetros del pipeline agrupados con descripciones
- **Output:** Descripciones y ejemplos de los archivos de salida esperados
- **Results:** Archivos de salida de ejemplo generados a partir del conjunto de datos de prueba completo
- **Releases & Statistics:** Historial de versiones del pipeline y estadísticas

Siempre que esté considerando adoptar un nuevo pipeline, debe leer la documentación del pipeline cuidadosamente primero para entender qué hace y cómo debe configurarse antes de intentar ejecutarlo.

Revíselo ahora y vea si puede descubrir:

- Qué herramientas ejecutará el pipeline (Consulte la pestaña: `Introduction`)
- Qué entradas y parámetros acepta o requiere el pipeline (Consulte la pestaña: `Parameters`)
- Cuáles son las salidas producidas por el pipeline (Consulte la pestaña: `Output`)

#### 1.1.1. Descripción general del pipeline

La pestaña `Introduction` proporciona una descripción general del pipeline, incluyendo una representación visual (llamada mapa de metro) y una lista de herramientas que se ejecutan como parte del pipeline.

![mapa de metro del pipeline](./img/nf-core-demo-subway-cropped.png)

1. Control de calidad de lecturas ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Recorte de adaptadores y calidad ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Presentación del control de calidad para lecturas sin procesar ([MULTIQC](http://multiqc.info/))
4. Generación de un mensaje de texto divertido desde una vaca ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Ejemplo de línea de comandos

La documentación también proporciona un archivo de entrada de ejemplo (que se analiza más adelante) y un ejemplo de línea de comandos.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Notará que el comando de ejemplo NO especifica un archivo de workflow, solo la referencia al repositorio del pipeline, `nf-core/demo`.

Cuando se invoca de esta manera, Nextflow asumirá que el código está organizado de cierta forma.
Recuperemos el código para poder examinar esta estructura.

### 1.2. Recuperar el código del pipeline

Una vez que hemos determinado que el pipeline parece ser adecuado para nuestros propósitos, vamos a probarlo.
Afortunadamente, Nextflow facilita la recuperación de pipelines desde repositorios correctamente formateados sin tener que descargar nada manualmente.

#### 1.2.1. Usar `nextflow pull`

Volvamos a la terminal y ejecutemos lo siguiente:

```bash
nextflow pull nf-core/demo
```

??? success "Salida del comando"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow realiza un `pull` del código del pipeline, lo que significa que descarga el repositorio completo en su unidad local.

Para ser claros, puede hacer esto con cualquier pipeline de Nextflow que esté configurado apropiadamente en GitHub, no solo con pipelines de nf-core.
Sin embargo, nf-core es la colección de código abierto más grande de pipelines de Nextflow.

#### 1.2.2. Usar `nextflow list`

Puede pedirle a Nextflow que le muestre una lista de los pipelines que ha recuperado de esta manera:

```bash
nextflow list
```

??? success "Salida del comando"

    ```console
    nf-core/demo
    ```

Puede intentar descargar algunos otros pipelines para ver cómo se listan cuando tiene más de uno.

#### 1.2.3. Encontrar dónde se descargó el pipeline

Notará que los archivos no están en su directorio de trabajo actual.
Por defecto, Nextflow guarda los pipelines descargados en `$NXF_HOME/assets`.

Para encontrar dónde vive un pipeline específico, pregúntele directamente a Nextflow:

```bash
nextflow info nf-core/demo
```

??? success "Salida del comando"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Info"

    La ruta completa puede diferir en su sistema si no está usando nuestro entorno de capacitación.

Nextflow mantiene el código fuente descargado intencionalmente "fuera del camino" bajo el principio de que estos pipelines deben usarse más como bibliotecas que como código con el que interactuaría directamente.

Internamente, Nextflow almacena cada pipeline descargado como un repositorio git en `$NXF_HOME/assets/.repos/`, y extrae el código para cada revisión en un subdirectorio `clones/<commit>/`.
Dado que `.repos` es un directorio oculto, un simple `tree -L 2 $NXF_HOME/assets/` se verá vacío.

#### 1.2.4. Crear un enlace simbólico para acceder fácilmente al código fuente

No vamos a examinar el código en detalle, pero echemos un vistazo rápido solo para tener una idea de cómo se ve la organización general.

Para facilitar la exploración del código fuente del pipeline, cree un enlace simbólico que apunte a la copia extraída del pipeline:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Esto crea un acceso directo para que pueda explorar el código con `tree -L 2 pipelines/nf-core/demo` o abrir archivos directamente.

#### 1.2.5. Descripción general de la organización del código

Puede usar `tree` o el explorador de archivos para encontrar y abrir el directorio `nf-core/demo`.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Contenido del directorio"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows

    7 directories, 12 files
    ```

Como puede ver, hay mucho contenido allí, la mayor parte de lo cual no necesita preocuparse.

Brevemente, notemos que en el nivel superior puede encontrar un archivo README con información resumida, así como archivos accesorios que resumen información del proyecto, como licencias, pautas de contribución, citas y código de conducta.
La documentación detallada del pipeline se encuentra en el directorio `docs`.
Todo este contenido se usa para generar las páginas web del sitio de nf-core de forma programática, por lo que siempre están actualizadas con el código.

Para el resto, podemos distinguir tres grupos funcionales de archivos de código:

1. Componentes del código del pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuración del pipeline
3. Parámetros / entradas del pipeline y validación

No revisaremos los componentes del código del pipeline en esta parte del curso, pero sí abordaremos elementos de configuración y validación que probablemente sean relevantes para usted como usuario final de los pipelines de nf-core.

!!! tip "Consejo"

    También puede explorar el código fuente de cualquier pipeline de nf-core en GitHub, por ejemplo, [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Cada pipeline de nf-core sigue el mismo diseño de directorios, por lo que una vez que conozca la estructura, puede encontrar archivos de configuración, módulos y workflows para cualquier pipeline de la misma manera.

¡Por ahora, a ejecutar el pipeline!

### Conclusión

Ahora sabe cómo encontrar un pipeline a través del sitio web de nf-core y recuperar una copia local del código fuente.

### ¿Qué sigue?

Aprenda cómo probar un pipeline de nf-core con el mínimo esfuerzo.

---

## 2. Probar el pipeline con su perfil de prueba

Convenientemente, cada pipeline de nf-core viene con un perfil de prueba.
Este es un conjunto mínimo de configuraciones para que el pipeline se ejecute usando un pequeño conjunto de datos de prueba alojado en el repositorio [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
Es una excelente manera de probar rápidamente un pipeline a pequeña escala.

!!! tip "Consejo"

    El sistema de perfiles de configuración de Nextflow le permite cambiar fácilmente entre diferentes motores de contenedores o entornos de ejecución.
    Para más detalles, consulte [Hello Nextflow Parte 6: Configuración](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar el perfil de prueba

Es una buena práctica verificar qué especifica el perfil de prueba de un pipeline antes de ejecutarlo.
El perfil `test` para `nf-core/demo` se encuentra en el archivo de configuración `conf/test.config`.
Puede encontrarlo localmente dentro del código fuente del pipeline que descargó `nextflow pull`, a través del enlace simbólico `pipelines` creado en la sección 1.2.4:

```bash
code pipelines/nf-core/demo/conf/test.config
```

Aquí está el contenido de ese archivo:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Datos de entrada
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Notará de inmediato que el bloque de comentarios en la parte superior incluye un ejemplo de uso que muestra cómo ejecutar el pipeline con este perfil de prueba.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Las únicas cosas que necesitamos proporcionar son las que se muestran entre corchetes angulares en el comando de ejemplo: `<docker/singularity>` y `<OUTDIR>`.

Como recordatorio, `<docker/singularity>` se refiere a la elección del sistema de contenedores. Todos los pipelines de nf-core están diseñados para ser utilizables con contenedores (Docker, Singularity, etc.) para garantizar la reproducibilidad y eliminar problemas de instalación de software.
Por lo tanto, necesitaremos especificar si queremos usar Docker o Singularity para probar el pipeline.

La parte `--outdir <OUTDIR>` se refiere al directorio donde Nextflow escribirá las salidas del pipeline.
Necesitamos proporcionar un nombre para él, que podemos simplemente inventar.
Si no existe ya, Nextflow lo creará por nosotros en tiempo de ejecución.

Continuando con la sección después del bloque de comentarios, el perfil de prueba nos muestra qué se ha preconfigurado para las pruebas: en particular, el parámetro `input` ya está configurado para apuntar a un conjunto de datos de prueba, por lo que no necesitamos proporcionar nuestros propios datos.
Si sigue el enlace a la entrada preconfigurada, verá que es un archivo CSV que contiene identificadores de muestras y rutas de archivos para varias muestras experimentales.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Esto se llama una hoja de muestras (samplesheet), y es la forma más común de entrada para los pipelines de nf-core.
No se preocupe si no está familiarizado con los formatos y tipos de datos, no es importante para lo que sigue.

Ahora tenemos todo lo que necesitamos para probar el pipeline.

### 2.2. Ejecutar el pipeline

Como se indicó anteriormente, podemos usar el comando de prueba de ejemplo casi tal como está; solo necesitamos especificar qué empaquetado de software usar y cómo nombrar el directorio de salida.
Aquí usaremos Docker para el sistema de contenedores y `demo-results`, respectivamente.

Con eso, podemos ejecutar el comando de prueba:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


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
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
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

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Si su salida coincide con esa, ¡felicidades! Acaba de ejecutar su primer pipeline de nf-core.

Notará que hay mucha más salida en la consola que cuando ejecuta un pipeline básico de Nextflow.
Hay un encabezado que incluye un resumen de la versión del pipeline, entradas y salidas, y algunos elementos de configuración.

!!! info "Info"

    Su salida mostrará diferentes marcas de tiempo, nombres de ejecución y rutas de archivos, pero la estructura general y la ejecución de los procesos deberían ser similares.

Observe la línea cerca de la parte superior de la salida:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Esto le indica qué revisión del pipeline se utilizó.
Dado que no especificamos una versión, Nextflow usó el último commit en `master`.
Para ejecuciones reproducibles, debe fijar una versión específica usando la bandera `-r`:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

Esto garantiza que se use el mismo código del pipeline cada vez, independientemente de nuevos commits o versiones.
Para esta capacitación omitimos `-r` por simplicidad, pero en producción siempre debe especificarlo.

Continuando con la salida de ejecución, veamos las líneas que nos indican qué procesos se ejecutaron:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Esto nos indica que se ejecutaron cuatro procesos, correspondientes a las cuatro herramientas mostradas en la página de documentación del pipeline en el sitio web de nf-core: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` y `COWPY`.

Los nombres completos de los procesos como se muestran aquí, como `NFCORE_DEMO:DEMO:MULTIQC`, son más largos que lo que puede haber visto en el material introductorio de Hello Nextflow.
Estos incluyen los nombres de sus workflows padre y reflejan la modularidad del código del pipeline.
Si desea aprender a desarrollar pipelines al estilo nf-core usted mismo, consulte el curso [Build with nf-core](../nfcore_build/index.md).

### 2.3. Examinar las salidas del pipeline

Finalmente, echemos un vistazo al directorio `demo-results` producido por el pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contenido del directorio"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Puede parecer mucho.
Para obtener más información sobre las salidas del pipeline `nf-core/demo`, consulte su [página de documentación](https://nf-co.re/demo/1.2.0/docs/output/).

En esta etapa, lo importante a observar es que los resultados están organizados por módulo, y hay además un directorio llamado `pipeline_info` que contiene varios informes con marcas de tiempo sobre la ejecución del pipeline.

Por ejemplo, el archivo `execution_timeline_*` le muestra qué procesos se ejecutaron, en qué orden y cuánto tiempo tardaron en ejecutarse:

![informe de línea de tiempo de ejecución](./img/execution_timeline.png)

!!! info "Info"

    Aquí las tareas no se ejecutaron en paralelo porque estamos ejecutando en una máquina minimalista en Github Codespaces.
    Para ver estas ejecuciones en paralelo, intente aumentar la asignación de CPU de su codespace y los límites de recursos en la configuración de prueba.

Estos informes se generan automáticamente para todos los pipelines de nf-core.

### Conclusión

Sabe cómo ejecutar un pipeline de nf-core usando su perfil de prueba integrado y dónde encontrar sus salidas.

### ¿Qué sigue?

Continúe con la [Parte 2](./02_configure_execution.md), donde aprenderá cómo configurar la ejecución del pipeline.

---

## Resumen

En esta parte aprendió a:

- Encontrar y recuperar un pipeline de nf-core y examinar su estructura de código
- Ejecutar un pipeline usando su perfil de prueba integrado
