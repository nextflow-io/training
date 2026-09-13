# Parte 1: Adaptación al entorno de cómputo

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En [Nextflow Run](../nextflow_run/index.md), configuró las entradas, los parámetros y las salidas de un pipeline.
Este curso cubre la otra mitad del panorama: adaptar la ejecución de un pipeline a cualquier entorno de cómputo en el que se ejecute, sin modificar el código del workflow.

!!! example "Escenario"

    Desarrolló y probó su pipeline en su laptop usando Docker.
    Ahora necesita entregarlo: un colaborador solo tiene Conda configurado, y el clúster HPC de su institución espera que los trabajos pasen por su propio planificador con sus propios límites de recursos.
    Nada de eso debería requerir reescribir el pipeline en sí.

El mismo código del pipeline puede ejecutarse en todos estos lugares, porque nada de eso está integrado en el workflow.
El empaquetado de software, la plataforma de ejecución y la asignación de recursos se controlan mediante configuración, superpuesta al código, y eso es lo que cubre este curso: cómo adaptar el mismo pipeline a un nuevo entorno cambiando la configuración, no el código.

---

## 1. Seleccionar una tecnología de empaquetado de software

En [Nextflow Run](../nextflow_run/index.md), vio un perfil `conda` ya configurado en `nextflow.config` como alternativa a Docker.
Aquí construirá ese mismo cambio usted mismo, y verá qué se necesita para que un proceso sea realmente utilizable con Conda.

### 1.1. Deshabilitar Docker y habilitar Conda

Cambie `docker.enabled` a `false` y agregue una directiva que habilite Conda.

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Esto permite que Nextflow cree y use entornos Conda para cualquier proceso que tenga un paquete Conda especificado.
El proceso `cowpy` aún no tiene uno, así que agreguemos uno, completamente desde la configuración.

### 1.2. Agregar un paquete Conda mediante configuración

Una directiva `conda` puede establecerse en la definición del proceso en sí, de la misma manera en que `container` ya está en `modules/cowpy.nf`, pero no es obligatorio: `withName` le permite establecerla desde la configuración, con alcance limitado solo al proceso `cowpy`.

=== "Después"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Esto no reemplaza la directiva `container` que ya está en el código del pipeline, sino que agrega una alternativa junto a ella, sin tocar ese código en absoluto.

!!! tip "Consejo"

    La búsqueda en [Seqera Containers](https://seqera.io/containers/) es una forma conveniente de buscar el URI del paquete Conda para una herramienta determinada, incluso si no planea construir un contenedor a partir de él.

### 1.3. Ejecutar el workflow para verificar que puede usar Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/execution-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Esto produce la misma salida que ejecutar con Docker, aunque los mecanismos son diferentes en segundo plano: Nextflow recupera el paquete Conda y construye un entorno a partir de él, en lugar de descargar una imagen de contenedor.

!!! info "Info"

    Construir un nuevo entorno Conda puede tardar un poco más que descargar un contenedor la primera vez, pero el paquete utilizado aquí es pequeño, por lo que debería ser rápido.

Ahora vuelva a Docker para el resto de este curso.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Combinar Docker y Conda"

    Dado que estas configuraciones tienen alcance por proceso, puede combinarlas: algunos procesos usan Docker, otros usan Conda, dependiendo de lo que esté disponible para cada herramienta.
    Si tanto una directiva `container` (en el código del pipeline) como una directiva `conda` (aquí, desde la configuración) están establecidas para el mismo proceso y ambos sistemas de empaquetado están habilitados, Nextflow prioriza los contenedores.

### Conclusión

Sabe cómo configurar qué tecnología de empaquetado de software debe usar un proceso, y cómo cambiar entre Docker y Conda.

### ¿Qué sigue?

Aprenda cómo cambiar la plataforma de ejecución que Nextflow usa para ejecutar sus tareas.

---

## 2. Seleccionar una plataforma de ejecución

Cada pipeline que ha ejecutado hasta ahora ha usado el executor local: cada tarea se ejecuta en la misma máquina que Nextflow.
Nextflow verifica los CPUs y la memoria disponibles, y retiene las tareas hasta que haya suficientes recursos libres.

El executor local es conveniente, pero no escala más allá de una sola máquina.
Nextflow admite [muchos otros backends de ejecución](https://nextflow.io/docs/latest/executor.html), incluidos planificadores HPC (Slurm, LSF, SGE, PBS y otros) y plataformas en la nube (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes y más).

### 2.1. Apuntar a un backend diferente

El executor se establece mediante una directiva de proceso llamada `executor`.
Por defecto es `local`, por lo que lo siguiente está implícito:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Para apuntar a un backend diferente, establezca la directiva al executor que desee.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Advertencia"

    El entorno de capacitación no está conectado a un clúster HPC, por lo que esto no es algo que pueda ejecutar aquí.

### 2.2. La sintaxis específica del backend se abstrae

La mayoría de las plataformas HPC requieren que los envíos de trabajos especifiquen solicitudes de recursos, como CPUs, memoria y un nombre de cola, usando su propia sintaxis.
La misma solicitud de 8 CPUs y 4 GB de RAM en una cola llamada `my-science-work` se ve completamente diferente dependiendo del planificador.

??? abstract "Ejemplos"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow abstrae todo esto: usted especifica propiedades estandarizadas como `cpus`, `memory` y `queue` una sola vez (consulte las [directivas de proceso](https://nextflow.io/docs/latest/reference/process.html#process-directives) para la lista completa), y Nextflow las traduce a los scripts específicos del backend correspondiente en tiempo de ejecución.

### 2.3. Ver qué ejecuta realmente Nextflow

Esa traducción no es solo una conveniencia del archivo de configuración: está respaldada por algo concreto que puede inspeccionar ahora mismo, incluso con el executor local.
En [Nextflow Run, sección 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory), miró dentro de un directorio de tarea bajo `work/` y encontró `.command.sh`, el comando exacto que ejecutó Nextflow.
Ese mismo directorio también contiene un archivo que aún no examinó: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Salida del comando (extracto)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` es el script real que Nextflow entrega para su ejecución.
Envuelve `.command.sh` con todo lo necesario para ejecutarlo realmente: configuración del entorno, staging de entradas/salidas, y reporte del resultado de vuelta a Nextflow.
Con el executor `local`, Nextflow simplemente ejecuta este script en la misma máquina.

Esto es exactamente lo que cambia cuando establece un `executor` diferente.
Para un planificador HPC como Slurm o PBS, Nextflow genera ese mismo tipo de script envolvente, agrega el encabezado específico del planificador que vio en [2.2](#22-backend-specific-syntax-is-abstracted-away) (traducido desde sus configuraciones de `cpus`, `memory` y `queue`), y entrega el resultado al comando de envío propio de ese planificador, por ejemplo `sbatch` para Slurm.
A partir de ahí, Nextflow consulta al planificador el estado del trabajo en lugar de observar un proceso local directamente.
Los backends de lotes en la nube funcionan de manera un poco diferente, ya que se manejan mediante llamadas a API en lugar de un comando de envío, pero la misma idea subyacente aplica: el mismo script de tarea se ejecuta, solo cambia cómo se lanza y se rastrea.

### Conclusión

Sabe cómo cambiar el executor para apuntar a diferentes infraestructuras de cómputo, que Nextflow abstrae la sintaxis de envío específica del backend, y qué sucede realmente en segundo plano cuando una tarea se ejecuta en un backend diferente.

### ¿Qué sigue?

Continúe con la [Parte 2](./02_resources_and_retries.md), donde aprenderá cómo perfilar y asignar recursos de cómputo, y manejar fallos de tareas con reintentos.

---

## Resumen

En esta parte aprendió a:

- Cambiar la tecnología de empaquetado de software entre Docker y Conda
- Agregar una directiva `conda` a una definición de proceso
- Cambiar la plataforma de ejecución con la directiva `executor`
- Inspeccionar qué genera y ejecuta realmente Nextflow para una tarea, y cómo cambia eso entre distintos executors
