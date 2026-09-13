# Parte 3: Usar perfiles para cambiar configuraciones

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


En la [Parte 1](./01_packaging_and_execution.md) y la [Parte 2](./02_resources_and_retries.md), acumuló algunas opciones de configuración: empaquetado de software, plataforma de ejecución y asignación de recursos.
En la práctica, con frecuencia querrá cambiar entre conjuntos completos de estas opciones dependiendo de dónde esté ejecutando, por ejemplo, una laptop para desarrollo y un clúster HPC para producción.

Nextflow le permite configurar cualquier número de [perfiles](https://nextflow.io/docs/latest/config.html#profiles) que describen diferentes configuraciones, y seleccionar uno (o varios) en tiempo de ejecución con una sola bandera.

Ya usó uno: el perfil `test` de [Nextflow Run](../nextflow_run/index.md) reemplaza los parámetros de entrada con un conjunto pequeño y bien definido.
Ahora creará sus propios perfiles de infraestructura y los combinará con él.

---

## 1. Crear perfiles para diferentes entornos

### 1.1. Configurar los perfiles

Agregue dos perfiles a `nextflow.config`: uno para ejecutar en una laptop regular con Docker, y otro para un clúster HPC universitario con un planificador Slurm y Conda.

=== "Después"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="35"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

El perfil `univ_hpc` también establece límites de recursos, ya que esto generalmente es requerido en infraestructura HPC compartida.

### 1.2. Ejecutar el workflow con un perfil

Seleccione un perfil en tiempo de ejecución con `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Advertencia"

    El perfil `univ_hpc` no se ejecutará en el entorno de capacitación, ya que no hay un planificador Slurm disponible.

Si encuentra otras configuraciones que siempre van juntas, agréguelas al perfil correspondiente.
También puede crear perfiles adicionales para agrupar cualquier otra combinación que necesite.

### 1.3. Ejecutar con múltiples perfiles

Los perfiles no son mutuamente excluyentes.
Puede activar varios a la vez con `-profile <perfil1>,<perfil2>`.
Combine `my_laptop` con el perfil `test` que ya conoce de Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Los nombres de los archivos individuales recogen correctamente `batch = 'test'` del perfil `test` (`COLLECTED-test-output.txt`, y así sucesivamente).

Si combina perfiles que establecen la misma opción, Nextflow resuelve el conflicto usando el valor que lee en último lugar, es decir, el que aparece más tarde en el archivo.
Si las configuraciones en conflicto provienen de fuentes de configuración completamente diferentes, se aplica el [orden de precedencia](https://www.nextflow.io/docs/latest/config.html) estándar.

### Conclusión

Sabe cómo definir perfiles que agrupan configuración específica de infraestructura, seleccionar uno en tiempo de ejecución con `-profile`, combinar múltiples perfiles en una sola ejecución, y cómo Nextflow resuelve los conflictos cuando más de un perfil establece la misma opción.

### ¿Qué sigue?

Aprenda cómo inspeccionar la configuración completamente resuelta antes de ejecutar cualquier cosa.

---

## 2. Inspeccionar la configuración resuelta

Ya usó `nextflow config -profile test` en [Nextflow Run](../nextflow_run/02_configure_pipeline.md) para verificar a qué se resuelve un solo perfil.
Ese comando se vuelve especialmente útil una vez que está combinando múltiples perfiles: como acaba de ver, cuando dos perfiles establecen la misma opción, puede ser complicado determinar manualmente qué valor prevalece.
El comando `nextflow config` resuelve todo eso por usted, sin ejecutar el pipeline.

### 2.1. Resolver la configuración predeterminada

```bash
nextflow config
```

??? success "Salida del comando"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

Esto es exactamente lo que se aplicaría si ejecutara el pipeline sin banderas adicionales.

### 2.2. Resolver la configuración con perfiles activados

Agregue los mismos perfiles que usaría para una ejecución real.

```bash
nextflow config -profile my_laptop,test
```

??? success "Salida del comando"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

Comparar los dos confirma qué cambió: `params.batch`, `params.character` y `process.executor` reflejan los perfiles `my_laptop,test`.
Esto resulta especialmente valioso para pipelines con muchas capas de configuración, donde determinar manualmente las configuraciones resueltas sería tedioso y propenso a errores.

### Conclusión

Sabe cómo usar `nextflow config` para inspeccionar la configuración completamente resuelta para cualquier combinación de perfiles, antes de ejecutar cualquier cosa.

### ¿Qué sigue?

Ha cubierto los aspectos esenciales de la configuración de pipelines en Nextflow.
Consulte el [Resumen del curso](next_steps.md) para saber hacia dónde ir desde aquí.

---

## Resumen

En esta parte aprendió a:

- Definir perfiles que agrupan configuración específica de infraestructura
- Combinar múltiples perfiles en una sola ejecución, y entender cómo se resuelven los conflictos entre ellos
- Usar `nextflow config` para inspeccionar la configuración completamente resuelta
