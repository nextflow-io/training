# Parte 4: Ejecutar pipelines remotos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hasta ahora, ha ejecutado scripts de workflow almacenados localmente.
En la práctica, con frecuencia querrá ejecutar pipelines publicados en repositorios remotos, como GitHub, sin necesidad de descargarlos usted mismo.

Nextflow hace esto sencillo: puede ejecutar cualquier pipeline directamente desde la URL de un repositorio Git.

---

## 1. Ejecutar un pipeline desde GitHub

La sintaxis básica para ejecutar un pipeline remoto es `nextflow run <repository>`, donde `<repository>` puede ser una ruta de repositorio de GitHub como `nextflow-io/hello`, una URL completa, o una ruta a GitLab, Bitbucket u otro servicio de alojamiento Git.

### 1.1. Lanzar el pipeline

Ejecute el pipeline de demostración oficial "hello" de Nextflow.
Este es un pipeline diferente, mucho más simple que el que ha estado ejecutando en este curso: es anterior al pipeline "Hello" utilizado a lo largo de esta capacitación, y simplemente imprime un saludo para cada uno de algunos idiomas predefinidos, así que no espere la entrada CSV ni el arte ASCII al que está acostumbrado.

```bash
nextflow run nextflow-io/hello
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Encontrar dónde se almacena en caché el pipeline

La primera vez que ejecuta un pipeline remoto, Nextflow lo descarga y lo almacena en caché localmente.
Las ejecuciones posteriores reutilizan la versión en caché a menos que solicite explícitamente una actualización.

De forma predeterminada, Nextflow guarda los pipelines descargados en `$NXF_HOME/assets`.
Para encontrar dónde quedó un pipeline específico y qué revisiones están disponibles, consulte a Nextflow directamente:

```bash
nextflow info nextflow-io/hello
```

??? success "Salida del comando"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow marca cada revisión que ya ha descargado localmente con `>`; el resto están disponibles pero aún no se han obtenido en una copia de trabajo.

También puede listar todos los pipelines que ha descargado hasta ahora con `nextflow list`:

```bash
nextflow list
```

??? success "Salida del comando"

    ```console
    nextflow-io/hello
    ```

El curso [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) cubre este mecanismo de caché con mayor profundidad, incluyendo cómo explorar el código fuente de un pipeline descargado.

### Conclusión

Sabe cómo ejecutar un pipeline directamente desde un repositorio de GitHub sin descargarlo usted mismo, y dónde encontrarlo localmente después.

### ¿Qué sigue?

Aprenda cómo fijar una versión específica de un pipeline remoto para garantizar la reproducibilidad.

---

## 2. Especificar una versión para reproducibilidad

De forma predeterminada, Nextflow ejecuta la última revisión de la rama predeterminada.
Puede fijar una versión (tag), rama o commit específico usando el indicador `-r`.

### 2.1. Fijar una revisión específica

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow obtiene esta revisión la primera vez que la solicita, de ahí las líneas `Pulling` y `downloaded from`; solicitar la misma revisión nuevamente más adelante omite directamente el paso de descarga y va a `Launching`.
Fijar una revisión exacta es esencial para la reproducibilidad.
Garantiza que usted y sus colaboradores ejecuten exactamente el mismo código de pipeline, independientemente de los cambios que haya habido en el repositorio desde entonces.

### 2.2. Las revisiones se aplican solo por invocación

Fijar una revisión con `-r` solo afecta la ejecución en la que lo especifica: no cambia lo que usa una ejecución posterior con `nextflow run` sin argumentos adicionales.
Intente ejecutar el pipeline nuevamente sin `-r`:

```bash
nextflow run nextflow-io/hello
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Aunque la ejecución anterior fijó explícitamente `v1.3`, esta ejecución vuelve directamente a la rama predeterminada (`master`).
Nextflow mantiene una copia de trabajo local separada para cada revisión que ha utilizado, que es lo que muestran los marcadores `>` en `nextflow info`, pero nunca recuerda cuál ejecutó por última vez.
Puede encontrar el nombre de la rama predeterminada de un pipeline ejecutando `nextflow info <pipeline>`; es la que está marcada con `(default)`.
La reproducibilidad depende completamente de usted: siempre pase `-r` explícitamente cuando sea importante, en lugar de asumir que una revisión fijada en una ejecución anterior sigue aplicándose.

### Conclusión

Sabe cómo fijar un pipeline remoto a una versión, rama o commit específico para una ejecución reproducible, y que la fijación se aplica solo a esa invocación, no a ejecuciones posteriores.

### ¿Qué sigue?

Ha cubierto los fundamentos de la ejecución y gestión de pipelines de Nextflow.
Consulte el [Resumen del curso](next_steps.md) para saber hacia dónde continuar desde aquí.

---

## Resumen

En esta parte aprendió a:

- Ejecutar un pipeline directamente desde un repositorio de GitHub sin descargarlo
- Fijar un pipeline remoto a una revisión específica para garantizar la reproducibilidad, y comprender que la fijación se aplica solo a esa invocación
