# Parte 4: Configuración

En las Partes 1-3, aprendimos cómo ejecutar Nextflow, ejecutar un pipeline de nf-core y gestionar entradas con archivos de parámetros y samplesheets.
Ahora exploraremos cómo configurar pipelines para diferentes entornos de cómputo usando **archivos de configuración** y **perfiles**.

## Objetivos de aprendizaje

Al finalizar esta parte, podrás:

- Entender cómo Nextflow resuelve la configuración desde múltiples fuentes
- Usar perfiles integrados de nf-core para contenedores y pruebas
- Crear perfiles personalizados para diferentes entornos de cómputo
- Personalizar solicitudes de recursos usando etiquetas de proceso
- Gestionar límites de recursos en entornos restringidos
- Inspeccionar la configuración resuelta con `nextflow config`

---

## 1. Entendiendo la configuración de Nextflow

### 1.1. ¿Qué es un archivo de configuración?

Nextflow usa archivos de configuración para separar la **lógica del workflow** (qué hacer) de la **configuración de ejecución** (cómo y dónde hacerlo).

Los archivos de configuración controlan:

- Motores de contenedores (Docker, Singularity, Conda)
- Recursos de cómputo (CPUs, memoria, tiempo)
- Plataformas de ejecución (local, HPC, nube)
- Parámetros del pipeline

### 1.2. Precedencia de configuración

Nextflow carga la configuración desde múltiples fuentes, donde las fuentes posteriores sobrescriben a las anteriores:

1. **Configuración del pipeline**: `nextflow.config` en el repositorio del pipeline
2. **Configuración del directorio**: `nextflow.config` en tu directorio de trabajo actual
3. **Configuración del usuario**: `~/.nextflow/config`
4. **Línea de comandos**: Parámetros y opciones pasados directamente

Este enfoque por capas te permite mantener valores predeterminados en el pipeline, sobrescribir con configuraciones específicas del usuario y hacer ajustes rápidos en la línea de comandos.

### 1.3. Nuestra configuración actual

Veamos la configuración que hemos estado usando:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Comentemos o cambiemos la línea `docker.enabled = true` de la Parte 2, y descubramos cómo podemos lograr el mismo resultado usando un perfil en molkart.

---

## 2. Usando perfiles

### 2.1. ¿Qué son los perfiles?

Los perfiles son conjuntos nombrados de configuración que pueden activarse con la bandera `-profile` mediante el comando `nextflow run`.
Facilitan el cambio entre diferentes escenarios de cómputo sin editar archivos de configuración.

Todos los pipelines de nf-core vienen con varios perfiles predeterminados que podemos usar.

### 2.2. Inspeccionando perfiles integrados

Inspeccionémoslos en el archivo `molkart/nextflow.config` asociado con el código base del pipeline:

```bash
code molkart/nextflow.config
```

Busca el bloque `profiles`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Perfiles de contenedores comunes:

- `docker`: Usa contenedores Docker (más común para desarrollo local)
- `singularity`: Usa Singularity/Apptainer (común en HPC)
- `conda`: Usa entornos Conda
- `apptainer`: Usa contenedores Apptainer

### 2.3. Re-ejecutando con perfiles en lugar de nextflow.config

Ahora que hemos deshabilitado la configuración de docker en nuestro archivo local `nextflow.config` y entendemos los perfiles, volvamos a ejecutar el pipeline usando la bandera `-profile`.

Anteriormente en la Parte 3, creamos un archivo `params.yaml` con nuestros parámetros personalizados.
Ahora podemos combinarlo con el perfil Docker integrado:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Desglosemos qué hace cada bandera:

- `-profile docker`: Activa el perfil Docker del `nextflow.config` de molkart, que establece `docker.enabled = true`
- `-params-file params.yaml`: Carga todos los parámetros del pipeline desde nuestro archivo YAML
- `-resume`: Reutiliza resultados en caché de ejecuciones previas

Debido a que estamos usando `-resume`, Nextflow verificará si algo cambió desde la última ejecución.
Si los parámetros, entradas y código son los mismos, todas las tareas se recuperarán del caché y el pipeline se completará casi instantáneamente.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

¡Observa que todos los procesos muestran `cached: 2` o `cached: 1` - nada fue re-ejecutado!

### 2.4. Perfiles de prueba

Los perfiles de prueba proporcionan formas rápidas de especificar parámetros de entrada predeterminados y archivos de datos para permitirte verificar que el pipeline funciona.
Los pipelines de nf-core siempre incluirán al menos dos perfiles de prueba:

- `test`: Conjunto de datos pequeño con parámetros rápidos para pruebas rápidas
- `test_full`: Prueba más completa con datos más grandes

Veamos más de cerca el perfil `test` en molkart que se incluye usando la directiva `includeConfig`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Esto significa que cada vez que ejecutamos el pipeline con `-profile test`, Nextflow cargará la configuración desde `conf/test.config`.

```groovy title="molkart/conf/test.config (excerpt)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Observa que este perfil contiene los mismos parámetros que usamos en nuestro archivo `params.yaml` anteriormente.

Puedes activar múltiples perfiles separándolos con comas.
Usemos eso para probar nuestro pipeline sin necesitar nuestro archivo de parámetros:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Esto combina:

- `docker`: Habilita contenedores Docker
- `test`: Usa conjunto de datos y parámetros de prueba

Los perfiles se aplican de izquierda a derecha, por lo que los perfiles posteriores sobrescriben a los anteriores si establecen los mismos valores.

### Conclusión

Los pipelines de nf-core vienen con perfiles integrados para contenedores, pruebas y entornos especiales.
Puedes combinar múltiples perfiles para construir la configuración que necesitas.

### ¿Qué sigue?

Aprende cómo crear tus propios perfiles personalizados para diferentes entornos de cómputo.

---

## 3. Creando perfiles personalizados

### 3.1. Crear perfiles para cambiar entre desarrollo local y ejecución en HPC

Creemos perfiles personalizados para dos escenarios:

1. Desarrollo local con Docker
2. HPC universitario con planificador Slurm y Singularity

Agrega lo siguiente a tu `nextflow.config`:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Ahora puedes cambiar entre entornos fácilmente:

```bash
# Para desarrollo local
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Para HPC (cuando esté disponible)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! Note "Nota"

    No podemos probar el perfil HPC en este entorno de capacitación ya que no tenemos acceso a un planificador Slurm.
    Pero esto muestra cómo lo configurarías para uso en el mundo real.

### 3.2. Usar `nextflow config` para inspeccionar la configuración

El comando `nextflow config` muestra la configuración completamente resuelta sin ejecutar el pipeline.

Ver la configuración predeterminada:

```bash
nextflow config ./molkart
```

Ver la configuración con un perfil específico:

```bash
nextflow config -profile local_dev ./molkart
```

Esto es extremadamente útil para:

- Depurar problemas de configuración
- Entender qué valores se usarán realmente
- Verificar cómo interactúan múltiples perfiles

### Conclusión

Los perfiles personalizados te permiten cambiar entre diferentes entornos de cómputo con una sola bandera de línea de comandos.
Usa `nextflow config` para inspeccionar la configuración resuelta antes de ejecutar.

### ¿Qué sigue?

Aprende cómo personalizar solicitudes de recursos para procesos individuales usando el sistema de etiquetas de proceso de nf-core.

---

## 4. Personalizando solicitudes de recursos

### 4.1. Entendiendo las etiquetas de proceso en pipelines de nf-core

Para simplificar, los pipelines de nf-core usan [**etiquetas de proceso**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) para estandarizar la asignación de recursos en todos los pipelines.
Cada proceso está etiquetado con una etiqueta como `process_low`, `process_medium` o `process_high` para describir requisitos de recursos de cómputo bajos, medios o altos, respectivamente.
Estas etiquetas se convierten en solicitudes de recursos específicas en uno de los archivos de configuración ubicados en el directorio `conf/` del pipeline.

```groovy title="molkart/conf/base.config (excerpt)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Observa el multiplicador `task.attempt` - esto permite que los reintentos subsecuentes de tareas soliciten más recursos, si el pipeline está configurado con `process.maxRetries > 1`.

### 4.2. Sobrescribiendo recursos para procesos específicos

Para un control más detallado, apunta a procesos individuales por nombre:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Si intentamos ejecutar este pipeline con la sobrescritura anterior, el proceso `CELLPOSE` solicitará 16 CPUs y 32 GB de memoria en lugar del valor predeterminado definido por su etiqueta.
Esto causará que el pipeline falle en nuestro entorno actual ya que no tenemos tanta RAM disponible.
Aprenderemos cómo prevenir estos tipos de fallas en la siguiente sección.

!!! Tip "Consejo"

    Para encontrar nombres de procesos, revisa la salida de ejecución del pipeline o verifica `.nextflow.log`.
    Los nombres de procesos siguen el patrón `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Conclusión

Los pipelines de nf-core usan etiquetas de proceso para estandarizar la asignación de recursos.
Puedes sobrescribir recursos por etiqueta (afecta múltiples procesos) o por nombre (afecta un proceso específico).

### ¿Qué sigue?

Aprende cómo gestionar límites de recursos en entornos restringidos como GitHub Codespaces.

---

## 5. Gestionando recursos en entornos restringidos

### 5.1. El problema de los límites de recursos

Si intentáramos ejecutar molkart con un proceso solicitando 16 CPUs y 32 GB de memoria (como se mostró en la sección 4.2), fallaría en nuestro entorno actual porque no tenemos tantos recursos disponibles.
En un entorno de clúster con nodos más grandes, tales solicitudes se enviarían al planificador.

En entornos restringidos como GitHub Codespaces, sin límites, Nextflow se negaría a ejecutar procesos que excedan los recursos disponibles.

### 5.2. Estableciendo límites de recursos

La directiva `resourceLimits` limita las solicitudes de recursos a valores especificados:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Esto le dice a Nextflow: "Si algún proceso solicita más de 2 CPUs o 7 GB de memoria, limítalo a estos valores en su lugar."

### 5.3. Agregando límites de recursos a perfiles personalizados

Actualiza tus perfiles personalizados para incluir límites apropiados:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! Warning "Advertencia"

    Establecer límites de recursos demasiado bajos puede causar que los procesos fallen o se ejecuten lentamente.
    El pipeline puede necesitar usar algoritmos menos intensivos en memoria o procesar datos en fragmentos más pequeños.

### Conclusión

Usa `resourceLimits` para ejecutar pipelines en entornos con recursos restringidos limitando las solicitudes de recursos de los procesos.
Diferentes perfiles pueden tener diferentes límites apropiados para su entorno.

### ¿Qué sigue?

¡Has completado la capacitación principal de Nextflow para Bioimagen!

---

## Conclusión

Ahora entiendes cómo configurar pipelines de Nextflow para diferentes entornos de cómputo.

Habilidades clave que has aprendido:

- **Precedencia de configuración**: Cómo Nextflow resuelve configuraciones desde múltiples fuentes
- **Perfiles de nf-core**: Usar perfiles integrados para contenedores, pruebas y utilidades
- **Perfiles personalizados**: Crear tus propios perfiles para diferentes entornos
- **Etiquetas de proceso**: Entender y sobrescribir solicitudes de recursos por etiqueta
- **Límites de recursos**: Gestionar entornos restringidos con `resourceLimits`
- **Inspección de configuración**: Usar `nextflow config` para depurar y verificar configuraciones

Estas habilidades de configuración son transferibles a cualquier pipeline de Nextflow y te ayudarán a ejecutar workflows eficientemente en máquinas locales, clústeres HPC y plataformas en la nube.

### ¿Qué sigue?

¡Felicidades por completar el curso de Nextflow para Bioimagen!

Próximos pasos:

- Completa la encuesta del curso para proporcionar retroalimentación
- Consulta [Hello Nextflow](../hello_nextflow/index.md) para aprender más sobre el desarrollo de workflows
- Explora [Hello nf-core](../hello_nf-core/index.md) para profundizar en las herramientas de nf-core
- Navega otros cursos en las [colecciones de capacitación](../training_collections/index.md)
