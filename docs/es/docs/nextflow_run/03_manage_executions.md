# Parte 3: Gestionar ejecuciones de workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A medida que ejecuta y vuelve a ejecutar pipelines, acumula historial de ejecuciones y directorios `work/` antiguos.
En la [Parte 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) ya utilizó `-resume` para omitir el trabajo que ya estaba hecho.
Aquí aprenderá a generar reportes sobre una ejecución, inspeccionar el historial de ejecuciones pasadas con [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), y eliminar directorios de trabajo antiguos que ya no necesita con [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

---

## 1. Generar reportes del pipeline

Nextflow puede generar varios tipos de reportes sobre una ejecución, cada uno añadido con su propio indicador `-with-*`: un reporte de ejecución (`-with-report`), una línea de tiempo de ejecución (`-with-timeline`), un archivo de traza de tareas (`-with-trace`), y un diagrama del workflow (`-with-dag`).
Aquí generaremos los dos primeros; consulte [Execution reports](https://nextflow.io/docs/latest/reports.html) en la referencia de Nextflow para el resto.

### 1.1. Generar un reporte de ejecución

Agregue `-with-report` a cualquier comando `nextflow run` para generar un reporte HTML después de que el pipeline finalice:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

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

Nextflow escribe el reporte en un archivo llamado `report-<timestamp>.html` en el directorio de trabajo.
Ábralo en un navegador para ver un resumen de la ejecución, una tabla de cada tarea con su estado y tiempo de ejecución, y gráficos de uso de recursos desglosados por proceso.

La pestaña **Tasks** lista cada tarea que ejecutó el pipeline, con el nombre del proceso, el estado y el uso de recursos:

![Tabla de tareas del reporte de ejecución](img/execution_report_tasks.png)

El reporte es especialmente útil cuando un pipeline tarda más de lo esperado o una tarea falla: la tabla de tareas muestra exactamente dónde se invirtió el tiempo y qué tareas tuvieron éxito o fallaron.

### 1.2. Generar una línea de tiempo de ejecución

Agregue `-with-timeline` a una ejecución para obtener una vista estilo diagrama de Gantt de cuándo se ejecutó cada tarea:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

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

Nextflow escribe la línea de tiempo en un archivo llamado `timeline-<timestamp>.html`.
Ábralo en un navegador para ver una barra por cada tarea, posicionada y dimensionada según cuándo se ejecutó y cuánto tiempo tomó:

![Línea de tiempo de ejecución](img/execution_timeline.png)

La línea de tiempo hace visible de un vistazo la forma de expansión y contracción de la [Parte 1](./01_run_nextflow.md#31-run-the-workflow): las tres tareas `sayHello` se ejecutan en paralelo, luego las tres tareas `convertToUpper`, y después `collectGreetings` y `cowpy` se ejecutan una tras otra, ya que cada una depende de todo lo anterior.

### Conclusión

Sabe cómo generar un reporte de ejecución HTML con `-with-report` y una línea de tiempo de ejecución con `-with-timeline`, y dónde buscar los otros tipos de reportes que admite Nextflow.

### ¿Qué sigue?

Aprenda a inspeccionar el historial de ejecuciones pasadas.

---

## 2. Inspeccionar el registro de ejecuciones pasadas

Ya sea que esté desarrollando un pipeline o ejecutándolo en producción, en algún momento necesitará consultar información sobre ejecuciones pasadas.

### 2.1. El archivo de historial

Cada vez que lanza un workflow de Nextflow, se escribe una línea en un archivo de registro llamado `history`, dentro de un directorio oculto llamado `.nextflow` en el directorio de trabajo actual.

??? abstract "Contenido del archivo"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Cada línea le proporciona la marca de tiempo, la duración, el nombre de la ejecución, el estado, el ID de revisión, el ID de sesión y la línea de comando completa de una ejecución lanzada desde este directorio.

Observe las dos últimas líneas: son dos invocaciones separadas (una normal y una con `-resume`) del mismo comando exacto, y comparten el mismo ID de sesión.
El ID de sesión solo cambia cuando lanza una ejecución genuinamente nueva; usar `-resume` lo mantiene, que es como Nextflow sabe qué caché reutilizar.

### 2.2. Usar `nextflow log` para una vista más amigable

Leer el archivo de historial sin procesar funciona, pero `nextflow log` formatea la misma información con un encabezado:

```bash
nextflow log
```

??? success "Salida del comando"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow agrupa la información de caché que utiliza para `-resume` bajo `.nextflow/cache`, indexada por ID de sesión.
Por eso, buscar el nombre de ejecución o el ID de sesión correcto aquí es el primer paso cuando necesita investigar o limpiar una ejecución pasada.

### Conclusión

Sabe dónde Nextflow registra el historial de ejecuciones pasadas y cómo inspeccionarlo con `nextflow log`.

### ¿Qué sigue?

Aprenda a eliminar directorios de trabajo antiguos que ya no necesita.

---

## 3. Eliminar directorios de trabajo antiguos

Cada ejecución deja sus directorios de tareas en `work/`, incluso después de haber copiado las salidas que le interesan a `results/`.
Si ejecuta suficientes pipelines durante el desarrollo, esos subdirectorios se acumulan, por lo que Nextflow proporciona `nextflow clean` para eliminar los que ya no necesita.

### 3.1. Determinar los criterios de eliminación

`nextflow clean` admite varias formas de seleccionar qué eliminar; consulte la [documentación de referencia](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para la lista completa.
Aquí eliminará todo lo correspondiente a ejecuciones anteriores a una ejecución determinada, usando su nombre de ejecución.

Busque la ejecución más reciente que desea conservar usando `nextflow log`; en el [ejemplo de la sección 2.2](#22-use-nextflow-log-for-a-friendlier-view) esa es `elegant_panini`, la última ejecución normal antes de la ejecución con `-resume`.
El nombre de ejecución es la cadena de dos partes generada automáticamente que aparece en la línea de consola `Launching (...)`, o en la columna `RUN NAME` de `nextflow log`.

### 3.2. Realizar una ejecución de prueba

Agregue `-n` primero para verificar qué eliminaría un comando determinado sin eliminar nada realmente:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Salida del comando"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

Son 16 directorios de tareas: las 8 tareas de la ejecución `turkey` más las 8 de la ejecución `tux`, exactamente las que cabría esperar para dos ejecuciones completas de este pipeline de cuatro procesos.
La ejecución `elegant_panini` en sí misma, y las tareas en caché que la ejecución con `-resume` reutilizó de ella, se dejan intactas.

Su salida listará nombres de directorios diferentes, y la cantidad de líneas que obtenga depende de cuántas ejecuciones haya realizado; si no ve ninguna línea, es posible que el nombre de ejecución no coincida con ninguno en su registro, o que no haya nada que eliminar antes de él.

### 3.3. Proceder con la eliminación

Una vez que la ejecución de prueba se vea correcta, vuelva a ejecutar el mismo comando con `-f` en lugar de `-n`:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Salida del comando"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` vacía los directorios de tareas pero deja en su lugar los directorios padre de dos caracteres (como `e5/`).

!!! warning "Advertencia"

    Eliminar directorios de trabajo de ejecuciones pasadas los elimina del caché de Nextflow y borra cualquier salida almacenada únicamente en ellos.
    Esto rompe la capacidad de Nextflow para reanudar la ejecución sin volver a ejecutar los procesos correspondientes, así que solo limpie las ejecuciones de las que esté seguro de que no necesitará reanudar.
    Esta es también la razón por la que vale la pena publicar todo lo que le importe en `results/` con `mode 'copy'` en lugar de depender del directorio `work/` o de un modo de publicación con `symlink`.

### Conclusión

Sabe cómo eliminar directorios de trabajo antiguos con `nextflow clean`, y por qué hacerlo implica perder la capacidad de reanudar desde esas ejecuciones.

### ¿Qué sigue?

Aprenda a ejecutar pipelines directamente desde repositorios remotos como GitHub en la [Parte 4](./04_remote_repositories.md).

---

## Resumen

En esta parte aprendió a:

- Generar un reporte de ejecución HTML con `-with-report` y una línea de tiempo de ejecución con `-with-timeline`
- Inspeccionar el historial de ejecuciones pasadas con `nextflow log`
- Eliminar directorios de trabajo antiguos con `nextflow clean`, y comprender la compensación con resume que esto conlleva
