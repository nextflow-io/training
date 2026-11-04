# Parte 1: Ejecutar Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta parte, presentamos los conceptos fundamentales para ejecutar pipelines de Nextflow.
Comenzamos con un workflow simple de Hello World y luego avanzamos hacia un pipeline completo de múltiples pasos que procesa varias entradas en paralelo usando contenedores.

---

## 1. Hello World

El workflow `1-hello.nf` recibe un saludo mediante un argumento de línea de comandos y lo escribe en un archivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Lanzar el workflow

Ejecute el siguiente comando en su terminal.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Salida del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

La línea clave en la salida es la línea de estado del proceso:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Esto nos indica que el proceso `sayHello` se ejecutó correctamente una vez.
El prefijo `[6d/740edd]` es una ruta truncada al directorio de trabajo de la tarea — más información sobre esto a continuación.
El bloque `Outputs:` que sigue lista todos los archivos que el pipeline publicó, etiquetados según el bloque `output` que se describe en [1.4](#14-optional-code-walkthrough) más abajo.

### 1.2. Encontrar la salida

Este workflow está configurado para publicar su salida en un directorio `results`.
Después de ejecutarlo, debería encontrar la salida allí:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Abra el archivo para confirmar que contiene `Hello World!`.

### 1.3. Explorar el directorio `work/`

En segundo plano, Nextflow crea un directorio de tarea único para cada llamada a un proceso dentro de un directorio llamado `work/`.
El hash que aparece en la salida de la consola (`[6d/740edd]`) es la ruta a ese directorio.

```bash
ls work/6d/740edd*
```

Dentro encontrará el archivo de salida junto con varios archivos de registro ocultos:

- **`.command.sh`**: el comando exacto que ejecutó Nextflow
- **`.command.out`** / **`.command.err`**: stdout y stderr del proceso
- **`.command.log`**: salida de registro combinada
- **`.exitcode`**: el código de salida del proceso

El archivo `.command.sh` es especialmente útil al depurar — muestra exactamente qué se ejecutó.

### 1.4. Opcional: Revisión del código

Entender el código no es esencial si solo quiere ejecutar pipelines, pero si tiene curiosidad, vale la pena echarle un vistazo.

??? optional "Haga clic para explorar el código asociado a este ejercicio"

    Abramos `1-hello.nf` y veamos sus componentes principales.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Parámetros del pipeline
     */
    params {
        input: String
    }

    workflow {

        main:
        // emite un saludo
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Vemos lo siguiente:

    - una declaración `include` que apunta a un módulo de `process`
    - un bloque `params` que define los parámetros del pipeline
    - un bloque `workflow` que describe el trabajo a realizar
    - un bloque `output` que describe qué hacer con las salidas

    Veamos cada uno por separado.

    ### El módulo `process`

    La declaración `include` le indica a Nextflow que cargue algo llamado `sayHello` desde un archivo de código separado.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    En ese archivo, encontramos la definición de un proceso llamado `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

    Un **process** define un único paso en el pipeline.
    Declara sus entradas, salidas y el script a ejecutar.
    El calificador `val` significa que la entrada es un valor simple (string, número, etc.).
    El calificador `path` significa que la salida es una ruta de archivo.

    Es posible escribir la definición del proceso en el archivo principal del workflow, pero mantenerlos en archivos de módulo separados los hace reutilizables: el mismo módulo puede ser importado por múltiples scripts de workflow.

    ### El bloque `params`

    El bloque `params` declara los parámetros de línea de comandos que acepta el workflow:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Cualquier parámetro declarado aquí estará disponible en la línea de comandos con doble guion (`--input`).
    Los tipos admitidos incluyen `String`, `Integer`, `Float`, `Boolean` y `Path`.

    !!! tip "Consejo"

        Los parámetros del workflow siempre usan dos guiones (`--input`) para distinguirlos de las propias opciones de la CLI de Nextflow, que usan un guion (por ejemplo, `-resume`).

    ### El bloque `workflow`

    El bloque **workflow** define la lógica del flujo de datos: qué procesos ejecutar y en qué orden.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // emite un saludo
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Aquí solo se llama a un proceso, por lo que es muy simple; cubriremos ejemplos más realistas más adelante.

    La sección `main:` llama al proceso `sayHello` con el valor de `--input`.
    La sección `publish:` lista qué salidas deben copiarse al directorio de resultados.

    ### El bloque `output`

    El bloque `output` al final del archivo especifica la ruta de destino y el modo de copia.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Cada entrada con nombre corresponde a una etiqueta `publish:` en el workflow y la mapea a un subdirectorio dentro de `results/`.

### Conclusión

Sabe cómo ejecutar un pipeline de Nextflow y encontrar sus salidas, y sabe que el trabajo se ejecuta en directorios de tarea dentro de `work/`.

### ¿Qué sigue?

Descubra cómo Nextflow maneja múltiples entradas de manera eficiente.

---

## 2. Procesar múltiples entradas

Los pipelines del mundo real típicamente procesan muchos datos, no solo uno.
El workflow `2-inputs.nf` lee desde un archivo CSV y ejecuta `sayHello` una vez por fila, en paralelo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Ejecutemos el workflow primero y luego veremos qué mecanismo usa Nextflow para manejar estas múltiples entradas.

### 2.1. Ejecutar el workflow

Ejecute el siguiente comando en su terminal.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Salida del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

El `3 of 3` nos indica que el proceso `sayHello` fue llamado tres veces, una por cada fila del CSV.

En el directorio `results`, ahora debería ver tres archivos de salida, uno por saludo:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Abra cualquiera de los archivos de salida para confirmar que cada uno contiene un saludo.

La salida condensada anterior muestra una única línea de resumen para `sayHello`, pero Nextflow en realidad lanzó tres ejecuciones de tarea separadas, una por fila del CSV, y las ejecutó en paralelo tan pronto como su máquina tuvo los recursos para hacerlo.

Al igual que la tarea individual que exploró en [1.3](#13-explore-the-work-directory), cada una de estas tres ejecuciones obtiene su propio directorio de tarea dentro de `work/`, completamente aislado de las demás:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Cada `.command.sh` solo contiene el comando para ese saludo en particular:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Este aislamiento es lo que hace que la ejecución en paralelo sea segura: tres tareas ejecutándose al mismo tiempo nunca comparten un directorio de trabajo, por lo que nada de lo que escribe una tarea puede colisionar con o sobrescribir lo que escribe otra tarea, incluso si producen archivos con el mismo nombre.
También es por eso que `-resume` (que se cubre a continuación) puede almacenar en caché y reutilizar tareas individuales de forma independiente: las entradas, salidas y registros de cada tarea viven completamente dentro de su propio directorio, sin nada compartido entre tareas que pueda desincronizarse.

### 2.2. Ejecutar el workflow nuevamente con `-ansi-log false`

Por defecto, Nextflow condensa la salida en una única línea de resumen por proceso.
Para ver cada llamada a proceso listada individualmente, agregue `-ansi-log false`:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Salida del comando"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Esto muestra las tres llamadas a proceso y el subdirectorio de trabajo único creado para cada una.

### 2.3. Usar `-resume` para omitir el trabajo completado

Ahora cambie al archivo de entrada extendido, que agrega dos saludos más, y añada `-resume` a la línea de comandos:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Salida del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow ejecutó solo las dos nuevas entradas.
Los tres saludos procesados en la ejecución anterior fueron almacenados en caché y reutilizados automáticamente.

Esto también funciona para omitir la ejecución de procesos para pasos que ya se completaron exitosamente en un pipeline de múltiples pasos.
Por ejemplo, si una ejecución del pipeline fue interrumpida por un error del sistema, o si agregó nuevos pasos a un pipeline en desarrollo.

La capacidad de `-resume` es especialmente valiosa en pipelines largos donde recuperarse de un fallo puede ahorrar tiempo y recursos críticos.

### 2.4. Opcional: Revisión del código

Entender el código no es esencial si solo quiere ejecutar pipelines, pero si tiene curiosidad, vale la pena echarle un vistazo.

??? optional "Haga clic para explorar el código asociado a este ejercicio"

    El cambio clave en `2-inputs.nf` está en la sección `main:` del workflow:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // crea un canal para las entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emite un saludo
        sayHello(greeting_ch)
    ```

    Lo que ve aquí se llama un **canal**: una construcción de cola que maneja los datos de entrada de una manera que facilita la paralelización de operaciones.

    - `channel.fromPath(params.input)` crea un canal a partir de la ruta de archivo proporcionada con `--input`
    - `.splitCsv()` analiza el CSV en filas
    - `#!groovy .map { line -> line[0] }` extrae la primera columna de cada fila

    El resultado es un canal que contiene `Hello`, `Bonjour` y `Hola`.
    Cuando se pasa a `sayHello(greeting_ch)`, Nextflow llama automáticamente al proceso una vez por elemento, ejecutándolos en paralelo cuando los recursos lo permiten.

### Conclusión

Sabe cómo procesar múltiples entradas desde un archivo CSV en paralelo, y cómo usar `-resume` para evitar repetir el trabajo completado.

### ¿Qué sigue?

Aprenda cómo un pipeline completo de múltiples pasos encadena procesos usando canales, y cómo usar contenedores para gestionar las herramientas de análisis y sus dependencias.

---

## 3. Ejecutar un pipeline de múltiples pasos

Hasta ahora ha ejecutado un único proceso y luego lo ha ejecutado múltiples veces en paralelo sobre un conjunto de entradas.
Los pipelines reales suelen ir más lejos: encadenan varios procesos, alimentando la salida de uno hacia el siguiente, y frecuentemente dependen de más de una pieza de software en el camino.
El workflow `main.nf` combina ambas cosas en un pipeline completo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Cada saludo de entrada fluye a través de los cuatro pasos: `sayHello` lo escribe en un archivo, `convertToUpper` convierte el texto a mayúsculas, `collectGreetings` fusiona todos los resultados en un archivo, y `cowpy` genera arte ASCII a partir de la salida fusionada usando una herramienta en contenedor.
Nextflow conecta estos pasos con canales: la salida de un proceso se convierte en la entrada del siguiente, por lo que toda la cadena se ejecuta automáticamente a medida que los datos están disponibles, sin que usted tenga que orquestar cada paso manualmente.

Tenga en cuenta que este workflow usa módulos: cada proceso está definido en su propio archivo dentro de `modules/`, y `main.nf` los importa con declaraciones `include` en lugar de definirlos directamente.
Esto hace que cada proceso sea reutilizable en múltiples workflows sin duplicar código. Para obtener más información, consulte la sección de exploración de código más abajo.

### 3.1. Ejecutar el workflow

Ejecute el siguiente comando en su terminal.

```bash
nextflow run main.nf --input data/greetings.csv
```

El parámetro `character` tiene como valor predeterminado `turkey` en `nextflow.config`, por lo que el arte ASCII usa un pavo a menos que lo cambie (pruebe agregando `--character tux`).

??? success "Salida del comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

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

Se ejecutaron cuatro procesos, pero no el mismo número de veces.
`sayHello` y `convertToUpper` se ejecutaron una vez por entrada (3 of 3): cada saludo necesita ser escrito y convertido a mayúsculas por separado.
`collectGreetings` y `cowpy` se ejecutaron solo una vez (1 of 1): fusionar los saludos y generar el arte ASCII solo tiene sentido una vez que todos los resultados individuales están listos.
Esta forma de expansión y luego contracción — varias tareas en paralelo que alimentan a un número menor de tareas posteriores — es común en los pipelines reales.

Nextflow no espera a que un paso completo termine antes de comenzar el siguiente.
Tan pronto como una salida de `sayHello` está lista, la tarea correspondiente de `convertToUpper` puede comenzar, por lo que las tareas de diferentes procesos se ejecutan de forma concurrente en lugar de en lotes estrictos.
`collectGreetings` y `cowpy` sí tienen que esperar, ya que cada uno depende de que todos los resultados anteriores estén disponibles primero.

El directorio `results` refleja esa contracción, más lo que el autor del pipeline eligió publicar y dónde: recuerde el bloque `output` de la revisión de código en 1.4, que es lo que define esta estructura.

```console title="results/"
results
└── batch
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

El directorio de nivel superior recibe el nombre del parámetro `batch`, que tiene como valor predeterminado `batch`; verá que cambia en ejercicios posteriores.

Revise `cowpy-COLLECTED-batch-output.txt` para ver el archivo de arte ASCII.

??? abstract "Contenido del archivo"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
     ---------
      \                                  ,+*^^*+___+++_
       \                           ,*^^^^              )
        \                       _+*                     ^**+_
         \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                 {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
               {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
               U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
             (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
               (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                 ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                         ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Al igual que en [2.1](#21-run-the-workflow), cada una de estas 8 ejecuciones de tarea, en los cuatro procesos, obtiene su propio directorio dentro de `work/`, completamente aislado de las demás.
`collectGreetings` es un buen ejemplo de por qué eso importa: depende de las salidas de las tres tareas de `convertToUpper`, que viven en tres directorios de tarea diferentes, por lo que Nextflow crea enlaces simbólicos a esos archivos dentro del propio directorio de `collectGreetings` en lugar de leer directamente desde los directorios de las tareas anteriores:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Cada tarea solo ve los archivos específicos que necesita, sin importar de dónde vengan, y nunca el contenido interno del directorio de otra tarea.
A lo largo de todo un pipeline, ese mismo aislamiento que vio con un único proceso en [2.1](#21-run-the-workflow) es lo que permite a Nextflow ejecutar cada tarea de cada proceso de forma concurrente y segura.

!!! note "Nota"

    El paso `cowpy` se ejecuta dentro de un contenedor Docker en lugar de depender del software instalado localmente.
    Un contenedor empaqueta una aplicación junto con todo lo que necesita para ejecutarse, por lo que no tiene que instalar y gestionar dependencias usted mismo, y el pipeline se comporta de la misma manera en cualquier máquina que pueda ejecutar el contenedor.
    Nextflow también admite Conda como alternativa a los contenedores; consulte la [Parte 2](./02_configure_pipeline.md) para saber cómo cambiar entre ellos.

### 3.2. Opcional: Revisión del código

Entender el código no es esencial si solo quiere ejecutar pipelines, pero si tiene curiosidad, vale la pena echarle un vistazo.

??? optional "Haga clic para explorar el código asociado a este ejercicio"

    ### Cómo fluyen los datos de un paso al siguiente

    Cada proceso pasa su canal de salida al siguiente:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // crea un canal para las entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    El patrón `processName.out` hace referencia al canal de salida de un proceso.

    El operador `.collect()` reúne todas las salidas individuales de `convertToUpper` en un único elemento de canal antes de pasarlas a `collectGreetings`.

    ### Usar módulos de proceso

    `main.nf` no define ningún código de proceso directamente.
    En cambio, importa cada proceso desde su propio archivo dentro de `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Cada archivo de módulo contiene una única definición de proceso, estructurada de la misma manera que el módulo `sayHello` en [1.4](#14-optional-code-walkthrough).
    Mantener los procesos en archivos separados los hace reutilizables en múltiples workflows sin duplicar código.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Usar software en contenedor

    El proceso `cowpy` se ejecuta dentro de un contenedor Docker especificado en su archivo de módulo:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

    Nextflow descarga automáticamente la imagen, ejecuta el script dentro del contenedor y limpia después.
    Docker está habilitado para este proyecto en `nextflow.config`:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Esta única línea habilita Docker para cualquier proceso del pipeline que tenga un contenedor especificado.

### Conclusión

Ha ejecutado un pipeline completo de múltiples pasos que procesa varias entradas en paralelo usando una herramienta en contenedor.

### ¿Qué sigue?

Continúe con la [Parte 2](./02_configure_pipeline.md), donde aprenderá cómo configurar el comportamiento del pipeline usando `nextflow.config`.

---

## Resumen

En esta parte aprendió a:

- Ejecutar un workflow de Nextflow y encontrar sus salidas
- Explorar el directorio `work/` y sus archivos de registro
- Procesar múltiples entradas desde un archivo CSV en paralelo
- Usar `-resume` para omitir el trabajo completado al agregar nuevas entradas
- Ejecutar un pipeline de múltiples pasos que usa una herramienta en contenedor
