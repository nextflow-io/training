# Parte 2: Gestionar recursos de cómputo y fallos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


En la [Parte 1](./01_packaging_and_execution.md), adaptó dónde y cómo se ejecutan las tareas de un pipeline.
Aquí aprenderá a controlar cuántos recursos de cómputo recibe cada tarea y qué sucede cuando una tarea falla a pesar de su mejor estimación de asignación.

---

## 1. Controlar las asignaciones de recursos de cómputo

Por defecto, Nextflow asigna un solo CPU a cada proceso mediante la directiva `cpus`, y no impone un límite de memoria a menos que usted establezca uno:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Ya sabe por [Nextflow Run](../nextflow_run/index.md) que la configuración de este pipeline establece `memory` en 1 GB para todos los procesos.
Pero, ¿cómo sabe qué valores usar realmente para sus propios pipelines?

### 1.1. Generar un informe de utilización de recursos

Ya generó un informe de ejecución con `-with-report` en [Nextflow Run](../nextflow_run/02_configure_pipeline.md).
Ese mismo informe es la forma de averiguar cuántos CPU y memoria necesitan realmente sus procesos: ejecute el workflow con algunas asignaciones predeterminadas, registre el uso real y luego ajuste a partir de ahí.

```bash
nextflow run main.nf -with-report report-config-1.html
```

El informe es un archivo HTML que puede abrir en un navegador.
Desglosa el tiempo de ejecución y la utilización de recursos por proceso, incluyendo qué porcentaje de los recursos asignados se usó realmente.
Esto es lo que muestra para `cowpy` con los valores predeterminados actuales (1 CPU, 1 GB de memoria):

| Métrica                  | Valor  |
| ------------------------ | ------ |
| Uso de CPU               | 116%   |
| Memoria máxima utilizada | 6.4 MB |
| Memoria asignada         | 1 GB   |

`cowpy` usa bastante menos del 1% de su asignación de 1 GB; el `%cpu` por encima del 100% simplemente significa que usa brevemente más de la capacidad de un CPU dentro del contenedor, en ráfagas cortas.

Consulte [Reports](https://nextflow.io/docs/latest/reports.html) para ver la lista completa de funcionalidades disponibles.

### 1.2. Establecer asignaciones de recursos para un proceso específico

El informe anterior muestra que `cowpy` está cómodamente dentro de su asignación actual, pero suponga que desea darle más margen de todas formas, por ejemplo porque espera entradas más grandes en producción.
Puede anular los valores predeterminados para un solo proceso con `withName`.

=== "Después"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Con esto en su lugar, cada proceso solicita 1 GB de memoria y un solo CPU, excepto `cowpy`, que solicita 2 GB y 2 CPUs (además de la configuración de `conda` de la [Parte 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Si su máquina tiene pocos CPUs y asigna un número alto por proceso, las llamadas a tareas pueden quedar en cola una detrás de otra, ya que Nextflow no solicitará más CPUs de los disponibles.

Ejecútelo de nuevo con un nombre de archivo de informe diferente para poder comparar antes y después.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Comparando los dos informes para `cowpy`:

| Métrica                  | Antes (1 CPU, 1 GB) | Después (2 CPUs, 2 GB) |
| ------------------------ | ------------------- | ---------------------- |
| Memoria máxima utilizada | 6.4 MB              | 6.4 MB                 |
| Uso de CPU               | 116%                | 118%                   |

Duplicar la asignación no cambió el uso real en absoluto, lo que indica que el 1 GB / 1 CPU original ya era generoso para esta carga de trabajo de prueba.
En un pipeline real que procese datos no triviales, esperaría que los números difieran significativamente entre procesos, que es exactamente por qué se hace un perfil antes de decidir qué asignar, en lugar de adivinar.

### 1.3. Agregar límites de recursos

Dependiendo de su infraestructura de cómputo, puede haber restricciones estrictas sobre lo que puede solicitar, por ejemplo un límite a nivel de clúster.
La directiva `resourceLimits` le permite establecer esos límites:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traduce estos valores a lo que el executor de destino espera.
Si un proceso solicita más que el límite, la solicitud se reduce en lugar de rechazarse.

!!! warning "Advertencia"

    Esto no es algo que pueda ejecutar en el entorno de capacitación, ya que requiere infraestructura HPC para tener efecto.

??? info "Configuraciones de referencia institucionales"

    El proyecto nf-core mantiene una [colección de archivos de configuración](https://nf-co.re/configs/) compartidos por instituciones de todo el mundo, que cubren una amplia variedad de executors HPC y en la nube.
    Son un punto de partida útil independientemente de si su propia institución está entre ellas.

### Conclusión

Sabe cómo generar un informe de perfil para evaluar la utilización de recursos, anular las asignaciones de recursos para un proceso específico y limitar las asignaciones con `resourceLimits`.

### ¿Qué sigue?

Aprenda cómo hacer que un pipeline se recupere automáticamente cuando una tarea falla, independientemente de si su estimación de asignación de recursos fue correcta.

---

## 2. Manejar fallos de tareas con reintentos

El perfil le indica lo que un proceso necesita la mayor parte del tiempo, pero las cargas de trabajo reales varían: una asignación que es cómoda para la mayoría de las entradas puede ser insuficiente para una entrada inusualmente grande, y las estimaciones simplemente pueden estar equivocadas.
En lugar de dejar que una sola tarea fallida detenga toda la ejecución, Nextflow puede reintentar una tarea fallida automáticamente, opcionalmente dándole más recursos en cada intento.

### 2.1. Reintentar una tarea fallida automáticamente

Para ver esto en acción, establezca deliberadamente la asignación de memoria de `cowpy` por debajo de lo que realmente necesita: recuerde de la [sección 1.1](#11-generate-a-resource-utilization-report) que alcanza un pico de alrededor de 6.4 MB, por lo que 6 MB debería ser justo por debajo de lo suficiente.

=== "Después"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` le indica a Nextflow qué hacer cuando una tarea falla: `'retry'` reenvía la tarea en lugar de detener todo el pipeline.
`maxRetries` limita cuántos intentos adicionales tiene antes de que Nextflow se rinda.

```bash
nextflow run main.nf
```

??? failure "Salida del comando (abreviada)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/config-exec/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

El código de salida 137 es la señal estándar para una terminación por falta de memoria: el contenedor no tenía suficiente memoria para ejecutar `cowpy` en absoluto.
Nextflow reintentó la tarea dos veces, tres intentos en total, coincidiendo con `maxRetries = 2`.
Como la asignación de memoria nunca cambió entre intentos, cada intento chocó contra la misma barrera; una vez agotados los reintentos, Nextflow reporta el fallo completo y detiene el pipeline, saliendo con un estado distinto de cero.

Reintentar por sí solo no soluciona nada si la causa subyacente no cambia entre intentos.

### 2.2. Aumentar los recursos en cada reintento

Dentro de una directiva de proceso, `task.attempt` contiene el número del intento actual, comenzando en 1.
Puede usarlo en un closure para escalar una asignación de recursos hacia arriba con cada reintento.

=== "Después"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Ejecute el workflow de nuevo:

```bash
nextflow run main.nf
```

??? success "Salida del comando (abreviada)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

El primer intento sigue fallando con 6 MB, pero el reintento se ejecuta con 12 MB (`6.MB * 2`) y tiene éxito, y el pipeline se completa con todas las salidas publicadas.

!!! warning "Advertencia"

    La salida de la consola todavía incluye una línea `NOTE:` que reporta el primer intento fallido, aunque el pipeline en su conjunto haya tenido éxito: Nextflow registra cada reintento individualmente, pero un fallo reintentado no afecta el resultado general.
    Verifique el resumen `Outputs:`, o el estado de salida del comando, para confirmar si la ejecución realmente tuvo éxito.

Consulte [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) en la documentación de Nextflow para patrones de reintento más avanzados, incluyendo el escalado basado en el error específico que ocurrió.

### Conclusión

Sabe cómo hacer que un pipeline reintente automáticamente las tareas fallidas y cómo escalar las asignaciones de recursos con cada reintento usando `task.attempt`.

### ¿Qué sigue?

Continúe con la [Parte 3](./03_profiles.md), donde aprenderá cómo agrupar configuraciones como esta en perfiles intercambiables.

---

## Resumen

En esta parte aprendió a:

- Generar un informe de perfil de recursos y establecer asignaciones de recursos por proceso
- Limitar las solicitudes de recursos con `resourceLimits`
- Reintentar automáticamente una tarea fallida con `errorStrategy` y `maxRetries`
- Escalar una asignación de recursos hacia arriba con cada reintento usando `task.attempt`
