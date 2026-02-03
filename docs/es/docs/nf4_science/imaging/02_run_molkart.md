# Parte 2: Ejecutar nf-core/molkart

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la Parte 1, ejecutamos un flujo de trabajo simple de Hello World para comprender los conceptos básicos de la ejecución de Nextflow.
Ahora vamos a ejecutar un pipeline de bioimagen del mundo real: **nf-core/molkart**.

Este pipeline procesa datos de transcriptómica espacial de Molecular Cartography de Resolve Bioscience.
Sin embargo, los patrones de Nextflow que aprenderá aquí se aplican a cualquier pipeline de nf-core o flujo de trabajo de producción.

## 1. Comprender los pipelines de nf-core

Antes de ejecutar el pipeline, entendamos qué es nf-core y por qué es importante para ejecutar flujos de trabajo.

### 1.1. ¿Qué es nf-core?

[nf-core](https://nf-co.re/) es una colección impulsada por la comunidad de pipelines de Nextflow de alta calidad.
Todos los pipelines de nf-core siguen la misma estructura y convenciones, lo que significa que una vez que aprende a ejecutar uno, puede ejecutar cualquiera de ellos.

Características clave de los pipelines de nf-core:

- **Estructura estandarizada**: Todos los pipelines tienen nombres de parámetros y patrones de uso consistentes
- **Datos de prueba integrados**: Cada pipeline incluye perfiles de prueba para validación rápida
- **Documentación completa**: Instrucciones de uso detalladas y descripciones de parámetros
- **Control de calidad**: Informes de QC automatizados usando MultiQC
- **Soporte de contenedores**: Contenedores pre-construidos para reproducibilidad

!!! tip "¿Quiere aprender más sobre nf-core?"

    Para una introducción detallada al desarrollo de pipelines de nf-core, consulte el curso de entrenamiento [Hello nf-core](../../hello_nf-core/index.md).
    Cubre cómo crear y personalizar pipelines de nf-core desde cero.

### 1.2. El pipeline molkart

![Pipeline nf-core/molkart](img/molkart.png)

El pipeline [nf-core/molkart](https://nf-co.re/molkart) procesa datos de imágenes de transcriptómica espacial a través de varias etapas:

1. **Preprocesamiento de imágenes**: Relleno de patrón de cuadrícula y mejora de contraste opcional
2. **Segmentación celular**: Múltiples opciones de algoritmos (Cellpose, Mesmer, ilastik, Stardist)
3. **Asignación de puntos**: Asignar puntos de transcripción a células segmentadas
4. **Control de calidad**: Generar informes de QC completos

Las salidas clave son:

- Tablas de conteo de células por transcripción
- Máscaras de segmentación
- Informe de control de calidad MultiQC

---

## 2. Ejecutar molkart con datos de prueba

Antes de comenzar, clonemos el repositorio de molkart localmente para que podamos inspeccionar su código:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Esto crea un directorio `molkart/` que contiene el código fuente completo del pipeline.

!!! note "¿Por qué estamos clonando localmente?"

    Típicamente, ejecutaría pipelines de nf-core directamente desde GitHub usando `nextflow run nf-core/molkart -r 1.2.0`.
    Nextflow descarga automáticamente la versión del pipeline solicitada para usted en `$HOME/.nextflow/assets/nf-core/molkart` y lo ejecuta desde allí.
    Sin embargo, para este entrenamiento, estamos clonando el pipeline a un directorio local diferente para que podamos inspeccionar el código más fácilmente.

### 2.1. Comprender los requisitos de contenedores

Antes de ejecutar el pipeline completo, aprendamos por qué los contenedores son esenciales para los pipelines de nf-core.

Intentemos ejecutar el pipeline usando el conjunto de datos de prueba y los parámetros de la configuración de prueba de molkart:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "mesmer,cellpose,stardist" \
  --outdir results
```

Desglosemos estos parámetros:

- `--input`: Ruta a la hoja de muestras que contiene metadatos de muestras
- `--mindagap_tilesize`, `--mindagap_boxsize`, `--mindagap_loopnum`: Parámetros para el relleno de patrón de cuadrícula
- `--clahe_pyramid_tile`: Tamaño del kernel para mejora de contraste
- `--segmentation_method`: Qué algoritmo(s) usar para la segmentación celular
- `--outdir`: Dónde guardar los resultados

!!! Warning "¡Este comando fallará - eso es intencional!"

    Estamos ejecutando esto deliberadamente sin contenedores para demostrar por qué son necesarios.

Después de unos momentos, verá un error como este:

??? failure "Salida del comando"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**¿Qué está sucediendo aquí?**

El error `command not found` (estado de salida 127) significa que Nextflow intentó ejecutar `duplicate_finder.py` pero no pudo encontrarlo en su sistema.
Esto se debe a que:

1. El pipeline espera que el software de bioinformática especializado esté instalado
2. Estas herramientas (como `duplicate_finder.py`, `apply_clahe.dask.py`, etc.) no son parte de las distribuciones estándar de Linux
3. Sin contenedores, Nextflow intenta ejecutar comandos directamente en su máquina local

**¿De dónde se supone que vienen estas herramientas?**

Inspeccionemos uno de los módulos de proceso para ver cómo declara sus requisitos de software.

Abra el módulo de preprocesamiento CLAHE:

```bash
code molkart/modules/local/clahe/main.nf
```

Mire la línea 5 - verá:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Esta línea le dice a Nextflow: "Para ejecutar este proceso, use la imagen de Docker `ghcr.io/schapirolabor/molkart-local:v0.0.4`, que contiene todo el software requerido."

Cada proceso declara qué imagen de contenedor proporciona sus herramientas requeridas.
¡Sin embargo, Nextflow solo usa estos contenedores si usted se lo indica!

**La solución: Habilitar Docker en la configuración**

### 2.2. Configurar Docker e iniciar el pipeline

Para habilitar Docker, necesitamos cambiar `docker.enabled` de `false` a `true` en el archivo `nextflow.config`.

Abra el archivo de configuración:

```bash
code nextflow.config
```

Cambie `docker.enabled = false` a `docker.enabled = true`:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Ahora ejecute el pipeline nuevamente con el mismo comando:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose,mesmer,stardist" \
  --outdir results
```

Esta vez, Nextflow:

1. Leerá la configuración `docker.enabled = true` del config
2. Descargará las imágenes de Docker requeridas (solo la primera vez)
3. Ejecutará cada proceso dentro de su contenedor especificado
4. Se ejecutará exitosamente porque todas las herramientas están disponibles dentro de los contenedores

!!! Tip "Por qué importan los contenedores"

    La mayoría de los pipelines de nf-core **requieren** contenedorización (Docker, Singularity, Podman, etc.) porque:

    - Utilizan software de bioinformática especializado no disponible en entornos estándar
    - Los contenedores aseguran reproducibilidad - las mismas versiones exactas de software se ejecutan en todas partes
    - No necesita instalar manualmente docenas de herramientas y sus dependencias

    Para más detalles sobre contenedores en Nextflow, consulte [Hello Containers](../../hello_nextflow/05_hello_containers.md) del entrenamiento Hello Nextflow.

### 2.3. Monitorear la ejecución

Mientras se ejecuta el pipeline, verá una salida similar a esta:

??? success "Salida del comando"

    ```console
    Nextflow 25.04.8 is available - Please consider updating your version to it

    N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/molkart` [soggy_kalam] DSL2 - revision: 5e54b29cb3 [dev]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0dev
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method       : mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize          : 7
      mindagap_loopnum          : 100
      clahe_kernel              : 25
      mindagap_tilesize         : 90
      clahe_pyramid_tile        : 368

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv
      outdir                    : results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-10-18_22-22-21

    Core Nextflow options
      revision                  : dev
      runName                   : soggy_kalam
      containerEngine           : docker
      launchDir                 : /workspaces/training/nf4-science/imaging
      workDir                   : /workspaces/training/nf4-science/imaging/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/molkart
      userName                  : root
      profile                   : docker,test
      configFiles               :

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [c1/da5009] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2 ✔
    [73/8f5e8a] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2 ✔
    [ec/8f84d5] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1 ✔
    [a2/99349b] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1 ✔
    [95/c9b4b1] NFCORE_MOLKART:MOLKART:DEEPCELL_MESMER (mem_only)          [100%] 1 of 1 ✔
    [d4/1ebd1e] NFCORE_MOLKART:MOLKART:STARDIST (mem_only)                 [100%] 1 of 1 ✔
    [3e/3c0736] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1 ✔
    [a0/415c6a] NFCORE_MOLKART:MOLKART:MASKFILTER (mem_only)               [100%] 3 of 3 ✔
    [14/a830c9] NFCORE_MOLKART:MOLKART:SPOT2CELL (mem_only)                [100%] 3 of 3 ✔
    [b5/391836] NFCORE_MOLKART:MOLKART:CREATE_ANNDATA (mem_only)           [100%] 3 of 3 ✔
    [77/aed558] NFCORE_MOLKART:MOLKART:MOLKARTQC (mem_only)                [100%] 3 of 3 ✔
    [e6/b81475] NFCORE_MOLKART:MOLKART:MULTIQC                             [100%] 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 19-Oct-2025 22:23:01
    Duration    : 2m 52s
    CPU hours   : 0.1
    Succeeded   : 22
    ```

Note cómo esta salida es más detallada que nuestro ejemplo de Hello World debido a las convenciones de nf-core que sigue el pipeline:

- El pipeline muestra su versión y logo
- Se muestran los parámetros de configuración
- Múltiples procesos se ejecutan en paralelo (indicado por múltiples líneas de proceso)
- Los nombres de proceso incluyen la ruta completa del módulo (por ejemplo, `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Comprender la ejecución de procesos

La línea del executor `executor > local (22)` le dice:

- **executor**: Qué entorno de cómputo se está usando (`local` = su máquina)
- **(22)**: Número total de tareas lanzadas

Cada línea de proceso muestra:

- **Hash** (`[1a/2b3c4d]`): Identificador del directorio de trabajo (como antes)
- **Nombre del proceso**: Ruta completa del módulo y nombre del proceso
- **Identificador de entrada**: Nombre de muestra entre paréntesis
- **Progreso**: Porcentaje completo y conteo (por ejemplo, `1 of 1 ✔`)

### Conclusión

Sabe cómo lanzar un pipeline de nf-core con datos de prueba e interpretar su salida de ejecución.

### ¿Qué sigue?

Aprenda dónde encontrar los resultados y cómo interpretarlos.

---

## 3. Encontrar y examinar las salidas

Cuando el pipeline se complete exitosamente, verá un mensaje de finalización y un resumen de ejecución.

### 3.1. Localizar el directorio de resultados

Por defecto, los pipelines de nf-core escriben las salidas en un directorio especificado por el parámetro `outdir`, que establecimos en `results/`.

Liste el contenido:

```bash
tree results/
```

Debería ver varios subdirectorios:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Cada subdirectorio contiene salidas de una etapa específica del pipeline:

- **mindagap/**: Imágenes con cuadrícula rellenada del paso de preprocesamiento MindaGap
- **clahe/**: Imágenes con contraste mejorado del preprocesamiento CLAHE
- **stack/**: Pilas de imágenes multicanal creadas para segmentación
- **segmentation/**: Resultados de segmentación de diferentes algoritmos (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Tablas de conteo de células por transcripción
- **anndata/**: Objetos AnnData que contienen matrices de células por transcripción y coordenadas espaciales
- **molkartqc/**: Métricas de control de calidad para asignación de puntos
- **multiqc/**: Informe de control de calidad completo
- **pipeline_info/**: Informes de ejecución y logs

### 3.2. Examinar el informe MultiQC

El informe MultiQC es un archivo HTML completo que agrega métricas de calidad de todos los pasos del pipeline.

Abra el informe en el explorador de archivos y luego haga clic en el botón "Show Preview" para verlo renderizado directamente en VS Code.

El informe incluye:

- Estadísticas generales para todas las muestras
- Métricas de preprocesamiento
- Métricas de calidad de segmentación
- Número de células y puntos detectados

!!! Tip

    Los informes MultiQC se incluyen típicamente en todos los pipelines de nf-core.
    Siempre proporcionan una visión general de alto nivel de la ejecución del pipeline y la calidad de los datos.

### 3.3. Examinar las tablas de células por transcripción

La salida científica más importante es la tabla de conteo de células por transcripción.
Esto le dice cuántos de cada transcripción se detectaron en cada célula.

Navegue al directorio spot2cell:

```bash
ls results/spot2cell/
```

Encontrará archivos como:

- `cellxgene_mem_only_cellpose.csv`: Tabla de células por transcripción usando segmentación Cellpose
- `cellxgene_mem_only_mesmer.csv`: Tabla de células por transcripción usando segmentación Mesmer
- `cellxgene_mem_only_stardist.csv`: Tabla de células por transcripción usando segmentación Stardist

Solo ejecutamos 1 muestra en este conjunto de datos de prueba, pero en un experimento real tendríamos estas tablas para cada muestra.
Note cómo Nextflow puede procesar múltiples métodos de segmentación en paralelo, haciendo fácil comparar resultados.

### 3.4. Ver informes de ejecución

Nextflow genera varios informes de ejecución automáticamente.

Verifique el directorio pipeline_info:

```bash
ls results/pipeline_info/
```

Archivos clave:

- **execution_report.html**: Visualización de línea de tiempo y uso de recursos
- **execution_timeline.html**: Gráfico de Gantt de ejecución de procesos
- **execution_trace.txt**: Métricas detalladas de ejecución de tareas
- **pipeline_dag.html**: Grafo acíclico dirigido que muestra la estructura del flujo de trabajo

Abra el informe de ejecución para ver el uso de recursos:

```bash
code results/pipeline_info/execution_report.html
```

Esto muestra:

- Cuánto tiempo tomó cada proceso
- Uso de CPU y memoria
- Qué tareas fueron cacheadas vs. ejecutadas

!!! Tip

    Estos informes son increíblemente útiles para optimizar la asignación de recursos y solucionar problemas de rendimiento.

### Conclusión

Sabe cómo localizar las salidas del pipeline, examinar informes de control de calidad y acceder a métricas de ejecución.

### ¿Qué sigue?

Aprenda sobre el directorio de trabajo y cómo Nextflow gestiona los archivos intermedios.

---

## 4. Explorar el directorio de trabajo

Al igual que con nuestro ejemplo de Hello World, todo el trabajo real ocurre en el directorio `work/`.

### 4.1. Comprender la estructura del directorio de trabajo

El directorio de trabajo contiene un subdirectorio para cada tarea que fue ejecutada.
Para este pipeline con 12 tareas, habrá 12 subdirectorios de trabajo.

Liste el directorio de trabajo:

```bash
ls -d work/*/*/ | head -5
```

Esto muestra los primeros 5 directorios de tareas.

### 4.2. Inspeccionar un directorio de tareas

Elija uno de los hashes de proceso de segmentación de la salida de consola (por ejemplo, `[3m/4n5o6p]`) y mire dentro:

```bash
ls -la work/3m/4n5o6p*/
```

Verá:

- **Archivos .command.\***: Scripts de ejecución de Nextflow y logs (como antes)
- **Archivos de entrada preparados**: Enlaces simbólicos a los archivos de entrada reales
- **Archivos de salida**: Máscaras de segmentación, resultados intermedios, etc.

La diferencia clave con Hello World:

- Los pipelines reales preparan archivos de entrada grandes (imágenes, datos de referencia)
- Los archivos de salida pueden ser bastante grandes (máscaras de segmentación, imágenes procesadas)
- Múltiples archivos de entrada y salida por tarea

!!! Tip

    Si un proceso falla, puede navegar a su directorio de trabajo, examinar `.command.err` para mensajes de error e incluso re-ejecutar `.command.sh` manualmente para depurar el problema.

### 4.3. Limpieza del directorio de trabajo

El directorio de trabajo puede volverse bastante grande con múltiples ejecuciones de pipeline.
Como aprendimos en la Parte 1, puede usar `nextflow clean` para eliminar directorios de trabajo de ejecuciones antiguas.

Sin embargo, para pipelines de nf-core con archivos intermedios grandes, es especialmente importante limpiar regularmente.

### Conclusión

Comprende cómo los pipelines de nf-core organizan sus directorios de trabajo y cómo inspeccionar tareas individuales para depuración.

### ¿Qué sigue?

Aprenda sobre el caché de Nextflow y cómo reanudar ejecuciones de pipeline fallidas.

---

## 5. Reanudar una ejecución de pipeline

Una de las características más poderosas de Nextflow es la capacidad de reanudar un pipeline desde el punto de falla.

### 5.1. El mecanismo de caché

Cuando ejecuta un pipeline con `-resume`, Nextflow:

1. Verifica el caché para cada tarea
2. Si las entradas, código y parámetros son idénticos, reutiliza el resultado cacheado
3. Solo re-ejecuta tareas que cambiaron o fallaron

Esto es esencial para pipelines de larga ejecución donde las fallas pueden ocurrir tarde en la ejecución.

### 5.2. Probar resume con molkart

Ejecute el mismo comando nuevamente, pero agregue `-resume`:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results \
  -resume
```

Debería ver una salida como: <!-- TODO: full output -->

```console
executor >  local (0)
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)        [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)                    [100%] 2 of 2, cached: 2 ✔
[7f/8g9h0i] NFCORE_MOLKART:MOLKART:CREATE_STACK (mem_only)             [100%] 1 of 1, cached: 1 ✔
[9h/0i1j2k] NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only) [100%] 1 of 1, cached: 1 ✔
[2k/3l4m5n] NFCORE_MOLKART:MOLKART:CELLPOSE (mem_only)                 [100%] 1 of 1, cached: 1 ✔
...
```

Note `cached: 2` o `cached: 1` para cada proceso - ¡nada fue re-ejecutado!

### 5.3. Cuándo resume es útil

Resume es particularmente valioso cuando:

- Un pipeline falla debido a límites de recursos (memoria insuficiente, límite de tiempo excedido)
- Necesita modificar procesos posteriores sin re-ejecutar pasos anteriores
- Su conexión de red se interrumpe durante la descarga de datos
- Desea agregar salidas adicionales sin rehacer el cómputo

!!! Warning

    Resume solo funciona si no ha cambiado los datos de entrada, el código del pipeline o los parámetros.
    Si cambia cualquiera de estos, Nextflow correctamente re-ejecutará las tareas afectadas.

### Conclusión

Sabe cómo usar `-resume` para re-ejecutar pipelines eficientemente sin repetir tareas exitosas.

### ¿Qué sigue?

Ahora que puede ejecutar nf-core/molkart con datos de prueba, está listo para aprender cómo configurarlo para sus propios conjuntos de datos.
