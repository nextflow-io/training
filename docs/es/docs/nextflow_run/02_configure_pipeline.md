# Parte 2: Configurar el pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la [Parte 1](./01_run_nextflow.md), ejecutó un pipeline completo de múltiples pasos que procesa múltiples entradas en paralelo usando contenedores.
Ahora vamos a ver cómo configurar el comportamiento del pipeline usando `nextflow.config`: primero examinando el archivo de configuración que ya le proporcionamos, luego explorando un par de otras formas de suministrar configuración, y finalmente controlando cómo y dónde se publican las salidas.

---

## 1. Examinar el archivo de configuración principal

Nextflow recoge automáticamente `nextflow.config` del directorio de trabajo y aplica su configuración a cada ejecución.

Le proporcionamos un archivo de configuración que cubre cuatro áreas: empaquetado de software, configuración de procesos, parámetros del pipeline y perfiles de ejecución.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Empaquetado de software
     */
    docker.enabled = true

    /*
     * Configuración de procesos
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Parámetros del pipeline
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Perfiles
     */
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

Revisemos cada uno, y luego pongamos los perfiles en práctica ejecutando el pipeline con uno de ellos.

!!! note "Nota"

    Esta configuración cubre la ejecución local en una sola máquina.
    Nextflow también admite planificadores HPC (SLURM, PBS, LSF) y ejecutores en la nube (AWS Batch, Google Cloud Batch, Azure Batch), todos configurados a través del mismo mecanismo de `nextflow.config`.
    Consulte la [Parte 1: Adaptarse a su entorno de cómputo](../execution_config/01_packaging_and_execution.md) en el curso [Execution Config](../execution_config/index.md) para un recorrido completo de estas opciones.

### 1.1. Empaquetado de software

El empaquetado de software es la forma en que Nextflow suministra las herramientas reales que sus procesos necesitan, ya sea una imagen de contenedor, un entorno Conda u otra cosa.

```groovy title="nextflow.config" linenums="1"
/*
 * Empaquetado de software
 */
docker.enabled = true
```

Esta línea habilita Docker para cada proceso.
Cualquier proceso que declare una directiva `container` se ejecuta dentro de la imagen especificada.

### 1.2. Configuración de procesos

Recuerde que un proceso es un paso individual en su pipeline, como `sayHello` o `cowpy`.
Nextflow le permite configurar varias cosas sobre cómo se ejecuta cada uno: cuánta CPU y memoria obtiene, qué contenedor o entorno Conda usa, y más.

```groovy title="nextflow.config" linenums="6"
/*
 * Configuración de procesos
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Esto limita cada proceso a una sola CPU y 1 GB de memoria.

Nextflow también le permite establecer valores diferentes para procesos individuales con nombre o grupos de procesos; aprenderá cómo en la [Parte 2: Gestionar recursos de cómputo y fallos](../execution_config/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) del curso [Execution Config](../execution_config/index.md).

### 1.3. Parámetros del pipeline

Los parámetros son las entradas de línea de comandos del pipeline, los mismos indicadores `--input`, `--batch` y `--character` que ya ha estado configurando directamente en la línea de comandos.
Establecer valores predeterminados aquí significa que no tiene que escribirlos cada vez, aunque como verá más adelante en esta parte, hay un par de otras formas de suministrarlos también.

```groovy title="nextflow.config" linenums="14"
/*
 * Parámetros del pipeline
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Estos valores predeterminados se activan cuando un parámetro no se suministra en la línea de comandos, por lo que ejecutar `nextflow run main.nf` sin indicadores sigue funcionando.

### 1.4. Perfiles

Los perfiles le permiten agrupar un conjunto de configuraciones bajo un solo nombre, para poder cambiar entre configuraciones completas con un solo indicador en lugar de cambiar valores manualmente cada vez.

```groovy title="nextflow.config" linenums="23"
/*
 * Perfiles
 */
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

El perfil `test` sobreescribe tres parámetros para ejecutar el pipeline con un conjunto de entradas pequeño y bien definido; cada pipeline de nf-core incluye uno de estos para validación rápida, y es una convención que vale la pena seguir en sus propios pipelines también.

El perfil `conda` cambia el empaquetado de software de Docker a Conda.

Se activa un perfil pasando `-profile <nombre>` en la línea de comandos.

Pongamos el perfil `test` en práctica.

```bash
nextflow run main.nf -profile test
```

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

El pipeline se ejecuta con `batch = 'test'` y `character = 'tux'`.
Revise `results/test/`: el nombre del lote ahora forma parte de la ruta del directorio, y el arte ASCII muestra el pingüino tux en lugar de un pavo.

!!! note "Nota"

    Puede activar varios perfiles a la vez, y usar `nextflow config -profile <nombre>,<nombre>` para ver el resultado completamente resuelto antes de ejecutar nada.
    La combinación de perfiles y cómo Nextflow resuelve los conflictos entre ellos se trata en profundidad en la [Parte 3: Usar perfiles para cambiar configuraciones](../execution_config/03_profiles.md) del curso [Execution Config](../execution_config/index.md).

### Conclusión

Sabe qué hacen los elementos más comunes de un archivo `nextflow.config` y cómo activar un perfil.

### ¿Qué sigue?

Aprenda un par de otras formas de suministrar valores de configuración sin modificar el archivo `nextflow.config` principal, útiles para configurar ejecuciones individuales y para compartir un conjunto exacto de configuraciones con otra persona.

---

## 2. Proporcionar configuración mediante archivos suplementarios

Establecer valores predeterminados en `nextflow.config` funciona bien para valores que rara vez cambian.
Nextflow también le ofrece dos mecanismos más específicos: un archivo de configuración específico para una ejecución, para adaptar la ejecución a un entorno particular, y un archivo de parámetros para compartir un conjunto exacto de valores de entrada con un colaborador.

### 2.1. Usar un archivo de configuración específico para una ejecución

Suponga que está moviendo el pipeline a una máquina que no tiene Docker y quiere darle a cada proceso más recursos con los que trabajar.
Cree un nuevo archivo de configuración con solo las sobreescrituras que necesita:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Páselo junto con su pipeline principal con `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow combina `custom.config` sobre el `nextflow.config` propio del pipeline, por lo que cada proceso ahora obtiene 2 CPUs y 2 GB de memoria en lugar de los valores predeterminados, y se ejecuta a través de Conda en lugar de Docker.
`cowpy` es el único proceso con un paquete Conda declarado junto a su contenedor, por lo que es el que verá que Nextflow construye un entorno para él:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Un archivo pequeño que solo sobreescribe la asignación de recursos y el empaquetado, sin tocar los parámetros del pipeline, es exactamente el patrón que los pipelines de nf-core esperan de las configuraciones institucionales.
Explore el repositorio [nf-core/configs](https://github.com/nf-core/configs) para ver ejemplos del mundo real.

Esto le da una forma desechable de adaptar un pipeline a un nuevo entorno sin tocar su configuración normal.

### 2.2. Usar un archivo de parámetros

Suponga que en cambio necesita compartir un conjunto exacto de parámetros de ejecución con un colaborador, o registrarlos para una publicación.

Nextflow le permite suministrar [archivos de parámetros](https://nextflow.io/docs/latest/config.html#parameter-file) en formato YAML o JSON, que son una forma más sencilla de distribuir un conjunto exacto y reproducible de valores.

Ya se proporciona un archivo de parámetros llamado `test-params.yaml` en su directorio de trabajo:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

La sintaxis usa dos puntos (`:`) en lugar de los signos de igual (`=`) usados en `nextflow.config`, ya que este archivo es YAML simple en lugar de Groovy.

!!! info "Info"

    También se proporciona una versión JSON, `test-params.json`. Siéntase libre de probarla por su cuenta; la sintaxis para pasarla es idéntica.

Pase el archivo con `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Contenido del archivo"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Un archivo de parámetros es especialmente valioso cuando un pipeline tiene más de un puñado de parámetros: le permite suministrarlos todos a la vez, sin una línea de comandos extensa ni ningún cambio en el script del workflow, y es fácil de distribuir junto con sus resultados.

### Conclusión

Sabe dos formas más de suministrar configuración: un archivo de configuración específico para una ejecución, para adaptar la ejecución a un nuevo entorno, y un archivo de parámetros para compartir valores de entrada exactos y reproducibles.

### ¿Qué sigue?

Aprenda cómo controlar cómo y dónde se publican las salidas de su pipeline.

---

## 3. Gestionar las salidas del pipeline

El autor de un pipeline decide cómo se organizan las salidas en el código, pero no necesita tocar ese código para controlar dónde terminan o cómo llegan allí.
Nextflow le ofrece formas de hacerlo a nivel de configuración: establecer un directorio base de salida y elegir si los archivos se copian o se enlazan simbólicamente.

### 3.1. Personalizar el directorio de salida

De forma predeterminada, Nextflow publica las salidas en `results/`.
Apúntelo a otro lugar con `-output-dir` (o su forma abreviada, `-o`):

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Contenido del directorio"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Las salidas ahora se ubican en `outputs/batch/` en lugar del valor predeterminado `results/batch/`.
El código propio del pipeline sigue decidiendo la estructura dentro de ese directorio base, como los subdirectorios `batch/` e `intermediates/`; `-output-dir` solo controla dónde comienza esa estructura.

`-output-dir` es en realidad solo un atajo de línea de comandos para la opción de configuración `outputDir`, por lo que puede ir en cualquier lugar donde pueda ir la configuración: directamente en `nextflow.config`, dentro de un perfil, o en un archivo de superposición `-c` como el que usó anteriormente en esta parte.
Por ejemplo, este fragmento muestra la misma configuración colocada directamente en `nextflow.config` en lugar de pasarla en la línea de comandos:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Consulte [Configuration file](https://nextflow.io/docs/latest/config.html) en la referencia de Nextflow para la lista completa de lugares donde puede vivir una opción de configuración como esta.

### 3.2. Elegir cómo se publican las salidas

De forma predeterminada, Nextflow publica las salidas como enlaces simbólicos que apuntan a las ubicaciones de las salidas en `work/`, no como copias reales:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Los autores del pipeline pueden establecer el 'modo de publicación' en `'copy'` o `'move'` para cada proceso individual en el código del workflow.
Típicamente hacen esto para las salidas finales del pipeline, mientras dejan el comportamiento predeterminado `'symlink'` para los archivos intermedios que pueden eliminarse una vez que el pipeline completo se ha ejecutado.

Eso evita duplicar datos en disco, pero significa que no puede eliminar los directorios de tareas en `work/` sin romper el enlace, perdiendo la capacidad de usar `-resume`.
Si desea que todos los archivos de salida se copien correctamente, establezca [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) en `'copy'` en la configuración de su pipeline. (A diferencia de `-output-dir`, no hay un indicador de línea de comandos para esto; es solo de configuración.)

Intente establecerlo en `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Luego ejecute el pipeline, cambiando el nombre del lote para que pueda ver la diferencia en las salidas:

```bash
nextflow run main.nf --batch withmode
```

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Revise uno de los archivos de salida como antes:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Ahora es un archivo real e independiente que seguirá disponible incluso si `work/` se limpia.

!!! warning "Advertencia"

    La configuración `workflow.output.mode` solo establece un valor predeterminado para las salidas que aún no tienen un modo establecido en el código del pipeline.
    No puede sobreescribir un modo que el autor haya codificado de forma fija, sin importar lo que establezca.

### Conclusión

Sabe cómo personalizar el directorio base de salida y elegir entre salidas copiadas y enlazadas simbólicamente, todo sin tocar el código del pipeline.

### ¿Qué sigue?

Continúe con la [Parte 3](./03_manage_executions.md), donde aprenderá cómo inspeccionar el historial de ejecuciones pasadas, generar informes de ejecución y limpiar directorios de trabajo antiguos.

---

## Resumen

En esta parte aprendió a:

- Configurar el comportamiento del pipeline usando `nextflow.config` y perfiles
- Suministrar configuración mediante un archivo de configuración específico para una ejecución o un archivo de parámetros
- Personalizar el directorio de salida y elegir entre salidas copiadas y enlazadas simbólicamente
