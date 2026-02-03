# Parte 2: Reescribir Hello para nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta segunda parte del curso de entrenamiento Hello nf-core, te mostramos cómo crear una versión compatible con nf-core del pipeline producido por el curso para principiantes [Hello Nextflow](../hello_nextflow/index.md).

Habrás notado en la primera sección del entrenamiento que los pipelines de nf-core siguen una estructura bastante elaborada con muchos archivos accesorios.
Crear todo eso desde cero sería muy tedioso, por lo que la comunidad nf-core ha desarrollado herramientas para hacerlo desde una plantilla, facilitando el proceso inicial.

Vamos a mostrarte cómo usar estas herramientas para crear una estructura base de pipeline, y luego adaptar el código de pipeline 'regular' existente sobre esta estructura nf-core.

Si no estás familiarizado con el pipeline Hello o necesitas un recordatorio, consulta [esta página de información](../info/hello_pipeline.md).

---

## 1. Crear un nuevo proyecto de pipeline

Primero, creamos la estructura base para el nuevo pipeline.

!!! note "Nota"

    Asegúrate de estar en el directorio `hello-nf-core` en tu terminal.

### 1.1. Ejecutar la herramienta de creación de pipeline basada en plantilla

Comencemos creando un nuevo pipeline con el comando `nf-core pipelines create`.
Esto creará una estructura base de pipeline usando la plantilla base de nf-core, personalizada con un nombre de pipeline, descripción y autor.

```bash
nf-core pipelines create
```

Ejecutar este comando abrirá una Interfaz de Usuario de Texto (TUI) para la creación del pipeline:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Esta TUI te pedirá que proporciones información básica sobre tu pipeline y te ofrecerá opciones de funcionalidades para incluir o excluir en tu estructura base.

- En la pantalla de bienvenida, haz clic en **Let's go!**.
- En la pantalla `Choose pipeline type`, haz clic en **Custom**.
- Ingresa los detalles de tu pipeline como se indica a continuación (reemplazando `< TU NOMBRE >` con tu propio nombre), luego haz clic en **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < TU NOMBRE >
```

- En la pantalla Template features, establece `Toggle all features` en **off**, luego **habilita** selectivamente las siguientes opciones. Verifica tus selecciones y haz clic en **Continue**.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- En la pantalla `Final details`, haz clic en **Finish**. Espera a que se cree el pipeline, luego haz clic en **Continue**.
- En la pantalla Create GitHub repository, haz clic en **Finish without creating a repo**. Esto mostrará instrucciones para crear un repositorio de GitHub más adelante. Ignóralas y haz clic en **Close**.

Una vez que la TUI se cierre, deberías ver la siguiente salida en la consola.

??? success "Salida del comando"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

No hay una confirmación explícita en la salida de consola de que la creación del pipeline funcionó, pero deberías ver un nuevo directorio llamado `core-hello`.

Visualiza el contenido del nuevo directorio para ver cuánto trabajo te ahorraste usando la plantilla.

```bash
tree core-hello
```

??? abstract "Contenido del directorio"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

¡Son muchos archivos!

Con suerte reconocerás muchos de ellos como los mismos que encontramos cuando exploramos la estructura del pipeline `nf-core/demo`.
Pero no te preocupes si todavía te sientes un poco perdido; recorreremos juntos las partes importantes en el transcurso de este entrenamiento.

!!! note "Nota"

    Una diferencia importante comparada con el pipeline `nf-core/demo` que examinamos en la primera parte de este entrenamiento es que no hay un directorio `modules`.
    Esto es porque no elegimos incluir ninguno de los módulos predeterminados de nf-core.

### 1.2. Probar que la estructura base es funcional

Lo creas o no, aunque aún no has agregado módulos para que haga trabajo real, la estructura base del pipeline puede ejecutarse usando el perfil de prueba, de la misma manera que ejecutamos el pipeline `nf-core/demo`.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Esto te muestra que toda la estructura básica está en su lugar.
Entonces, ¿dónde están las salidas? ¿Hay alguna?

De hecho, se creó un nuevo directorio de resultados llamado `core-hello-results` que contiene los reportes de ejecución estándar:

```bash
tree core-hello-results
```

??? abstract "Contenido del directorio"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Puedes echar un vistazo a los reportes para ver qué se ejecutó, y la respuesta es: ¡nada en absoluto!

![reporte de línea de tiempo de ejecución vacío](./img/execution_timeline_empty.png)

Echemos un vistazo a lo que realmente hay en el código.

### 1.3. Examinar el workflow placeholder

Si miras dentro del archivo `main.nf`, verás que importa un workflow llamado `HELLO` desde `workflows/hello`.

Esto es equivalente al workflow `workflows/demo.nf` que encontramos en la Parte 1, y sirve como workflow placeholder para nuestro workflow de interés, con algo de funcionalidad nf-core ya en su lugar.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // canal: samplesheet leído desde --input
    main:

    ch_versions = channel.empty()

    //
    // Recopilar y guardar versiones de software
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // canal: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

Comparado con un workflow básico de Nextflow como el desarrollado en [Hello Nextflow](../hello_nextflow/index.md), notarás algunas cosas que son nuevas aquí (líneas destacadas arriba):

- El bloque workflow tiene un nombre
- Las entradas del workflow se declaran usando la palabra clave `take:` y la construcción del canal se mueve al workflow padre
- El contenido del workflow se coloca dentro de un bloque `main:`
- Las salidas se declaran usando la palabra clave `emit:`

Estas son características opcionales de Nextflow que hacen que el workflow sea **componible**, lo que significa que puede ser llamado desde dentro de otro workflow.

!!! note "Workflows componibles en profundidad"

    La [Side Quest Workflows de Workflows](../side_quests/workflows_of_workflows.md) explora la composición de workflows con mucha más profundidad, incluyendo cómo componer múltiples workflows juntos y gestionar flujos de datos complejos entre ellos. Estamos introduciendo la componibilidad aquí porque es un requisito fundamental de la arquitectura de plantilla nf-core, que utiliza workflows anidados para organizar la inicialización del pipeline, el workflow de análisis principal y las tareas de finalización en componentes separados y reutilizables.

Vamos a necesitar conectar la lógica relevante de nuestro workflow de interés en esa estructura.
El primer paso para eso es hacer que nuestro workflow original sea componible.

### Conclusión

Ahora sabes cómo crear una estructura base de pipeline usando las herramientas nf-core.

### ¿Qué sigue?

Aprende cómo hacer que un workflow simple sea componible como preludio para hacerlo compatible con nf-core.

---

## 2. Hacer el workflow original de Hello Nextflow componible

Ahora es momento de trabajar integrando nuestro workflow en la estructura base nf-core.
Como recordatorio, estamos trabajando con el workflow presentado en nuestro curso de entrenamiento [Hello Nextflow](../hello_nextflow/index.md).

!!! tip "Consejo"

    Si no estás familiarizado con ese pipeline o necesitas un recordatorio, consulta [El pipeline Hello](../info/hello_pipeline.md).

Te proporcionamos una copia limpia y completamente funcional del workflow Hello Nextflow completado en el directorio `original-hello` junto con sus módulos y el archivo CSV predeterminado que espera usar como entrada.

```bash
tree original-hello/
```

??? abstract "Contenido del directorio"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

Siéntete libre de ejecutarlo para verificar que funciona:

```bash
nextflow run original-hello/hello.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Abramos el archivo de workflow `hello.nf` para inspeccionar el código, que se muestra completo a continuación (sin contar los procesos, que están en módulos):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Incluir módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // crear un canal para entradas desde un archivo CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emitir un saludo
  sayHello(greeting_ch)

  // convertir el saludo a mayúsculas
  convertToUpper(sayHello.out)

  // recopilar todos los saludos en un archivo
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // generar arte ASCII de los saludos con cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Como puedes ver, este workflow fue escrito como un workflow simple sin nombre que puede ejecutarse por sí solo.
Para hacerlo ejecutable desde dentro de un workflow padre como requiere la plantilla nf-core, necesitamos hacerlo **componible**.

Recorramos los cambios necesarios uno por uno.

### 2.1. Nombrar el workflow

Primero, démosle un nombre al workflow para poder referirnos a él desde un workflow padre.

=== "Después"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

Las mismas convenciones se aplican a los nombres de workflow que a los nombres de módulos.

### 2.2. Reemplazar la construcción de canal con `take`

Ahora, reemplaza la construcción del canal con una simple declaración `take` que declare las entradas esperadas.

=== "Después"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // channel of greetings
        greeting_ch
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // crear un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Esto deja los detalles de cómo se proporcionan las entradas al workflow padre.

Ya que estamos en eso, también podemos comentar la línea `params.greeting = 'greetings.csv'`

=== "Después"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Nota"

    Si tienes instalada la extensión del servidor de lenguaje Nextflow, el verificador de sintaxis iluminará tu código con garabatos rojos.
    Eso es porque si pones una declaración `take:`, también debes tener un `main:`.

    Agregaremos eso en el siguiente paso.

### 2.3. Anteponer las operaciones del workflow con la declaración `main`

A continuación, agrega una declaración `main` antes del resto de las operaciones llamadas en el cuerpo del workflow.

=== "Después"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // emitir un saludo
        sayHello(greeting_ch)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // recopilar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emitir un saludo
        sayHello(greeting_ch)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // recopilar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Esto básicamente dice 'esto es lo que este workflow _hace_'.

### 2.4. Agregar declaración `emit`

Finalmente, agrega una declaración `emit` que declare cuáles son las salidas finales del workflow.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Esta es una adición completamente nueva al código comparado con el workflow original.

### 2.5. Resumen de los cambios completados

Si has hecho todos los cambios como se describe, tu workflow ahora debería verse así:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Incluir módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // channel of greetings
    greeting_ch

    main:

    // emitir un saludo
    sayHello(greeting_ch)

    // convertir el saludo a mayúsculas
    convertToUpper(sayHello.out)

    // recopilar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // generar arte ASCII de los saludos con cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Esto describe todo lo que Nextflow necesita EXCEPTO qué alimentar en el canal de entrada.
Eso se va a definir en el workflow padre, también llamado workflow de **punto de entrada**.

### 2.6. Hacer un workflow de punto de entrada ficticio

Antes de integrar nuestro workflow componible en la compleja estructura base nf-core, verifiquemos que funcione correctamente.
Podemos hacer un workflow de punto de entrada ficticio simple para probar el workflow componible de forma aislada.

Crea un archivo en blanco llamado `main.nf` en el mismo directorio `original-hello`.

```bash
touch original-hello/main.nf
```

Copia el siguiente código en el archivo `main.nf`.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// importar el código del workflow desde el archivo hello.nf
include { HELLO } from './hello.nf'

// declarar parámetro de entrada
params.greeting = 'greetings.csv'

workflow {
  // crear un canal para entradas desde un archivo CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // llamar al workflow importado con el canal de saludos
  HELLO(greeting_ch)

  // ver las salidas emitidas por el workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

Hay dos observaciones importantes que hacer aquí:

- La sintaxis para llamar al workflow importado es esencialmente la misma que la sintaxis para llamar módulos.
- Todo lo que está relacionado con llevar las entradas al workflow (parámetro de entrada y construcción del canal) ahora se declara en este workflow padre.

!!! note "Nota"

    Nombrar el archivo de workflow de punto de entrada `main.nf` es una convención, no un requisito.

    Si sigues esta convención, puedes omitir especificar el nombre del archivo del workflow en tu comando `nextflow run`.
    Nextflow buscará automáticamente un archivo llamado `main.nf` en el directorio de ejecución.

    Sin embargo, puedes nombrar el archivo de workflow de punto de entrada de otra manera si lo prefieres.
    En ese caso, asegúrate de especificar el nombre del archivo del workflow en tu comando `nextflow run`.

### 2.7. Probar que el workflow se ejecuta

Finalmente tenemos todas las piezas que necesitamos para verificar que el workflow componible funciona.
¡Ejecutémoslo!

```bash
nextflow run ./original-hello
```

Aquí ves la ventaja de usar la convención de nomenclatura `main.nf`.
Si hubiéramos nombrado el workflow de punto de entrada `algo_mas.nf`, habríamos tenido que hacer `nextflow run original-hello/algo_mas.nf`.

Si hiciste todos los cambios correctamente, esto debería ejecutarse hasta completarse.

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Esto significa que hemos actualizado exitosamente nuestro workflow HELLO para que sea componible.

### Conclusión

Sabes cómo hacer un workflow componible dándole un nombre y agregando declaraciones `take`, `main` y `emit`, y cómo llamarlo desde un workflow de punto de entrada.

### ¿Qué sigue?

Aprende cómo injertar un workflow componible básico en la estructura base nf-core.

---

## 3. Ajustar la lógica del workflow actualizado en el workflow placeholder

Ahora que hemos verificado que nuestro workflow componible funciona correctamente, volvamos a la estructura base del pipeline nf-core que creamos en la sección 1.
Queremos integrar el workflow componible que acabamos de desarrollar en la estructura de plantilla nf-core, para que el resultado final se vea algo así.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Entonces, ¿cómo hacemos que eso suceda? Echemos un vistazo al contenido actual del workflow `HELLO` en `core-hello/workflows/hello.nf` (la estructura base nf-core).

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // canal: samplesheet leído desde --input
    main:

    ch_versions = channel.empty()

    //
    // Recopilar y guardar versiones de software
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // canal: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

En general, este código hace muy poco aparte de algunas tareas de mantenimiento que tienen que ver con capturar la versión de cualquier herramienta de software que se ejecute en el pipeline.

Necesitamos agregar el código relevante del workflow componible original que desarrollamos en la sección 2.

Vamos a abordar esto en las siguientes etapas:

1. Copiar los módulos y configurar las importaciones de módulos
2. Dejar la declaración `take` como está
3. Agregar la lógica del workflow al bloque `main`
4. Actualizar el bloque `emit`

!!! note "Nota"

    Vamos a ignorar la captura de versión por este primer paso y veremos cómo conectar eso en una parte posterior de este entrenamiento.

### 3.1. Copiar los módulos y configurar las importaciones de módulos

Los cuatro procesos de nuestro workflow Hello Nextflow están almacenados como módulos en `original-hello/modules/`.
Necesitamos copiar esos módulos en la estructura del proyecto nf-core (bajo `core-hello/modules/local/`) y agregar declaraciones de importación al archivo de workflow nf-core.

Primero copiemos los archivos de módulos de `original-hello/` a `core-hello/`:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Ahora deberías ver el directorio de módulos listado bajo `core-hello/`.

```bash
tree core-hello/modules
```

??? abstract "Contenido del directorio"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Ahora configuremos las declaraciones de importación de módulos.

Estas eran las declaraciones de importación en el workflow `original-hello/hello.nf`:

```groovy title="original-hello/hello.nf" linenums="9"
// Incluir módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Abre el archivo `core-hello/workflows/hello.nf` y transpón esas declaraciones de importación en él como se muestra a continuación.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Dos observaciones más interesantes aquí:

- Hemos adaptado el formato de las declaraciones de importación para seguir la convención de estilo nf-core.
- Hemos actualizado las rutas relativas a los módulos para reflejar que ahora están almacenados en un nivel diferente de anidamiento.

### 3.2. Dejar la declaración `take` como está

El proyecto nf-core tiene mucha funcionalidad preconstruida alrededor del concepto de samplesheet, que típicamente es un archivo CSV que contiene datos en columnas.
Como eso es esencialmente lo que es nuestro archivo `greetings.csv`, mantendremos la declaración `take` actual como está, y simplemente actualizaremos el nombre del canal de entrada en el siguiente paso.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // canal: samplesheet leído desde --input
```

El manejo de entrada se hará antes de este workflow (no en este archivo de código).

### 3.3. Agregar la lógica del workflow al bloque `main`

Ahora que nuestros módulos están disponibles para el workflow, podemos conectar la lógica del workflow en el bloque `main`.

Como recordatorio, este es el código relevante en el workflow original, que no cambió mucho cuando lo hicimos componible (solo agregamos la línea `main:`):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // emitir un saludo
    sayHello(greeting_ch)

    // convertir el saludo a mayúsculas
    convertToUpper(sayHello.out)

    // recopilar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // generar arte ASCII de los saludos con cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

Necesitamos copiar el código que viene después de `main:` en la nueva versión del workflow.

Ya hay algo de código allí que tiene que ver con capturar las versiones de las herramientas que ejecuta el workflow. Vamos a dejarlo solo por ahora (nos ocuparemos de las versiones de las herramientas más tarde).
Mantendremos la inicialización `ch_versions = channel.empty()` en la parte superior, luego insertaremos nuestra lógica de workflow, manteniendo el código de recopilación de versiones al final.
Este orden tiene sentido porque en un pipeline real, los procesos emitirían información de versión que se agregaría al canal `ch_versions` a medida que el workflow se ejecuta.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // canal: samplesheet leído desde --input

        main:

        ch_versions = Channel.empty()

        // emitir un saludo
        sayHello(greeting_ch)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // recopilar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Recopilar y guardar versiones de software
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // canal: [ path(versions.yml) ]

    }
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // canal: samplesheet leído desde --input
        main:

        ch_versions = Channel.empty()

        //
        // Recopilar y guardar versiones de software
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // canal: [ path(versions.yml) ]

    }
    ```

Notarás que también agregamos una línea en blanco antes de `main:` para hacer el código más legible.

Esto se ve bien, pero todavía necesitamos actualizar el nombre del canal que estamos pasando al proceso `sayHello()` de `greeting_ch` a `ch_samplesheet` como se muestra a continuación, para coincidir con lo que está escrito bajo la palabra clave `take:`.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir un saludo (updated to use the nf-core convention for samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir un saludo
        sayHello(greeting_ch)
    ```

Ahora la lógica del workflow está correctamente conectada.

### 3.4. Actualizar el bloque `emit`

Finalmente, necesitamos actualizar el bloque `emit` para incluir la declaración de las salidas finales del workflow.

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // canal: [ path(versions.yml) ]
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // canal: [ path(versions.yml) ]
    ```

Esto concluye las modificaciones que necesitamos hacer al workflow HELLO en sí mismo.
En este punto, hemos logrado la estructura general de código que nos propusimos implementar.

### Conclusión

Sabes cómo ajustar las piezas centrales de un workflow componible en un workflow placeholder nf-core.

### ¿Qué sigue?

Aprende cómo adaptar cómo se manejan las entradas en la estructura base del pipeline nf-core.

---

## 4. Adaptar el manejo de entradas

Ahora que hemos integrado exitosamente nuestra lógica de workflow en la estructura base nf-core, necesitamos abordar una pieza más crítica: asegurar que nuestros datos de entrada se procesen correctamente.
La plantilla nf-core viene con un manejo de entrada sofisticado diseñado para conjuntos de datos genómicos complejos, por lo que necesitamos adaptarlo para que funcione con nuestro archivo `greetings.csv` más simple.

### 4.1. Identificar dónde se manejan las entradas

El primer paso es averiguar dónde se realiza el manejo de entrada.

Puede que recuerdes que cuando reescribimos el workflow Hello Nextflow para que fuera componible, movimos la declaración del parámetro de entrada un nivel arriba, en el workflow de punto de entrada `main.nf`.
Así que echemos un vistazo al workflow de punto de entrada `main.nf` de nivel superior que se creó como parte de la estructura base del pipeline:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Ejecutar el pipeline de análisis principal según el tipo de entrada
//
workflow CORE_HELLO {

    take:
    samplesheet // canal: samplesheet leído desde --input

    main:

    //
    // WORKFLOW: Ejecutar pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Ejecutar tareas de inicialización
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Ejecutar workflow principal
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Ejecutar tareas de finalización
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

El proyecto nf-core hace un uso intensivo de subworkflows anidados, por lo que esta parte puede ser un poco confusa al principio.

Lo que importa aquí es que hay dos workflows definidos:

- `CORE_HELLO` es un envoltorio delgado para ejecutar el workflow HELLO que acabamos de terminar de adaptar en `core-hello/workflows/hello.nf`.
- Un workflow sin nombre que llama a `CORE_HELLO` así como a otros dos subworkflows, `PIPELINE_INITIALISATION` y `PIPELINE_COMPLETION`.

Aquí hay un diagrama de cómo se relacionan entre sí:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Importante, no podemos encontrar ningún código construyendo un canal de entrada en este nivel, solo referencias a un samplesheet proporcionado a través del parámetro `--input`.

Un poco de investigación revela que el manejo de entrada lo realiza el subworkflow `PIPELINE_INITIALISATION`, apropiadamente, que se importa de `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`.

Si abrimos ese archivo y nos desplazamos hacia abajo, llegamos a este fragmento de código:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

Esta es la fábrica de canales que analiza el samplesheet y lo pasa en una forma que está lista para ser consumida por el workflow HELLO.

!!! note "Nota"

    La sintaxis anterior es un poco diferente de lo que hemos usado anteriormente, pero básicamente esto:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    es equivalente a esto:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Este código involucra algunos pasos de análisis y validación que son altamente específicos del samplesheet de ejemplo incluido con la plantilla del pipeline nf-core, que al momento de escribir esto es muy específico del dominio y no es adecuado para nuestro proyecto de pipeline simple.

### 4.2. Reemplazar el código del canal de entrada de la plantilla

La buena noticia es que las necesidades de nuestro pipeline son mucho más simples, por lo que podemos reemplazar todo eso por el código de construcción de canal que desarrollamos en el workflow Hello Nextflow original.

Como recordatorio, así se veía la construcción del canal (como se ve en el directorio de soluciones):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // crear un canal para entradas desde un archivo CSV
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Así que solo necesitamos conectar eso en el workflow de inicialización, con cambios menores: actualizamos el nombre del canal de `greeting_ch` a `ch_samplesheet`, y el nombre del parámetro de `params.greeting` a `params.input` (ver línea destacada).

=== "Después"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Create channel from input file provided through params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Create channel from input file provided through params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Eso completa los cambios que necesitamos hacer para que el procesamiento de entrada funcione.

En su forma actual, esto no nos permitirá aprovechar las capacidades integradas de nf-core para la validación de esquema, pero podemos agregar eso más tarde.
Por ahora, nos estamos enfocando en mantenerlo lo más simple posible para llegar a algo que podamos ejecutar exitosamente en datos de prueba.

### 4.3. Actualizar el perfil de prueba

Hablando de datos de prueba y parámetros, actualicemos el perfil de prueba para este pipeline para usar el mini-samplesheet `greetings.csv` en lugar del samplesheet de ejemplo proporcionado en la plantilla.

Bajo `core-hello/conf`, encontramos dos perfiles de prueba de plantilla: `test.config` y `test_full.config`, que están destinados a probar una muestra pequeña de datos y una de tamaño completo.
Dado el propósito de nuestro pipeline, realmente no hay sentido en configurar un perfil de prueba de tamaño completo, así que siéntete libre de ignorar o eliminar `test_full.config`.
Nos vamos a enfocar en configurar `test.config` para que se ejecute en nuestro archivo `greetings.csv` con algunos parámetros predeterminados.

#### 4.3.1. Copiar el archivo `greetings.csv`

Primero necesitamos copiar el archivo `greetings.csv` a un lugar apropiado en nuestro proyecto de pipeline.
Típicamente los archivos de prueba pequeños se almacenan en el directorio `assets`, así que copiemos el archivo desde nuestro directorio de trabajo.

```bash
cp greetings.csv core-hello/assets/.
```

Ahora el archivo `greetings.csv` está listo para ser usado como entrada de prueba.

#### 4.3.2. Actualizar el archivo `test.config`

Ahora podemos actualizar el archivo `test.config` de la siguiente manera:

=== "Después"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Datos de entrada
        input  = "${projectDir}/assets/greetings.csv"

        // Other parameters
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Antes"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Datos de entrada
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Puntos clave:

- **Usando `${projectDir}`**: Esta es una variable implícita de Nextflow que apunta al directorio donde se encuentra el script de workflow principal (la raíz del pipeline). Usarla asegura que la ruta funcione independientemente de dónde se ejecute el pipeline.
- **Rutas absolutas**: Al usar `${projectDir}`, creamos una ruta absoluta, lo cual es importante para datos de prueba que se envían con el pipeline.
- **Ubicación de datos de prueba**: Los pipelines nf-core típicamente almacenan datos de prueba en el directorio `assets/` dentro del repositorio del pipeline para archivos de prueba pequeños, o referencian conjuntos de datos de prueba externos para archivos más grandes.

Y ya que estamos en eso, ajustemos los límites de recursos predeterminados para asegurar que esto se ejecutará en máquinas muy básicas (como las VMs mínimas en Github Codespaces):

=== "Después"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Antes"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Esto completa las modificaciones de código que necesitamos hacer.

### 4.4. Ejecutar el pipeline con el perfil de prueba

Eso fue mucho, ¡pero finalmente podemos intentar ejecutar el pipeline!
Ten en cuenta que tenemos que agregar `--validate_params false` a la línea de comando porque no configuramos la validación todavía (eso vendrá más tarde).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Si has hecho todas las modificaciones correctamente, debería ejecutarse hasta completarse.

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Como puedes ver, esto produjo el resumen típico de nf-core al inicio gracias al subworkflow de inicialización, y las líneas para cada módulo ahora muestran los nombres completos PIPELINE:WORKFLOW:module.

### 4.5. Encontrar las salidas del pipeline

La pregunta ahora es: ¿dónde están las salidas del pipeline?
Y la respuesta es bastante interesante: ahora hay dos lugares diferentes donde buscar los resultados.

Como puede que recuerdes de antes, nuestra primera ejecución del workflow recién creado produjo un directorio llamado `core-hello-results/` que contenía varios reportes de ejecución y metadatos.

```bash
tree core-hello-results
```

??? abstract "Contenido del directorio"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

Ves que obtuvimos otro conjunto de reportes de ejecución además de los que obtuvimos de la primera ejecución, cuando el workflow era solo un placeholder.
Esta vez ves todas las tareas que se ejecutaron como se esperaba.

![reporte de línea de tiempo de ejecución para el pipeline Hello](./img/execution_timeline_hello.png)

!!! note "Nota"

    Una vez más, las tareas no se ejecutaron en paralelo porque estamos ejecutando en una máquina minimalista en Github Codespaces.
    Para ver que se ejecuten en paralelo, intenta aumentar la asignación de CPU de tu codespace y los límites de recursos en la configuración de prueba.

Eso es genial, ¡pero nuestros resultados reales del pipeline no están ahí!

Aquí está lo que pasó: no cambiamos nada de los módulos en sí, por lo que las salidas manejadas por las directivas `publishDir` a nivel de módulo todavía van a un directorio `results` como se especificó en el pipeline original.

```bash
tree results
```

??? abstract "Contenido del directorio"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, ahí están, mezclados con las salidas de ejecuciones anteriores del pipeline Hello original.

Si queremos que estén organizados ordenadamente como lo estaban las salidas del pipeline demo, necesitaremos cambiar cómo configuramos las salidas para ser publicadas.
Te mostraremos cómo hacer eso más adelante en este curso de entrenamiento.

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

¡Y ahí lo tienes! Puede parecer mucho trabajo para lograr el mismo resultado que el pipeline original, pero obtienes todos esos reportes hermosos generados automáticamente, y ahora tienes una base sólida para aprovechar las características adicionales de nf-core, incluyendo la validación de entrada y algunas capacidades ingeniosas de manejo de metadatos que cubriremos en una sección posterior.

---

### Conclusión

Sabes cómo convertir un pipeline regular de Nextflow en un pipeline de estilo nf-core usando la plantilla nf-core.
Como parte de eso, aprendiste cómo hacer un workflow componible, y cómo identificar los elementos de la plantilla nf-core que más comúnmente necesitan ser adaptados al desarrollar un pipeline de estilo nf-core personalizado.

### ¿Qué sigue?

¡Toma un descanso, fue un trabajo duro! Cuando estés listo, continúa con [Parte 3: Usar un módulo nf-core](./03_use_module.md) para aprender cómo aprovechar módulos mantenidos por la comunidad del repositorio nf-core/modules.
