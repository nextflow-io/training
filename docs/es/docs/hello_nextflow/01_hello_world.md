# Parte 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea [la lista de reproducción completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/01_hello_world.md).
///
-->

En esta primera parte del curso de entrenamiento Hello Nextflow, nos introducimos al tema con un ejemplo muy básico de Hello World independiente del dominio, que iremos construyendo progresivamente para demostrar el uso de la lógica y componentes fundamentales de Nextflow.

??? info "¿Qué es un ejemplo Hello World?"

    Un "Hello World!" es un ejemplo minimalista destinado a demostrar la sintaxis básica y la estructura de un lenguaje de programación o framework de software.
    El ejemplo típicamente consiste en imprimir la frase "Hello, World!" al dispositivo de salida, como la consola o terminal, o escribirla en un archivo.

---

## 0. Calentamiento: Ejecutar un ejemplo Hello World directamente

Demostremos esto con un comando simple que ejecutamos directamente en el terminal, para mostrar lo que hace antes de envolverlo en Nextflow.

!!! tip "Consejo"

    Recuerde que ahora debería estar dentro del directorio `hello-nextflow/` como se describe en la página de [Comenzando](00_orientation.md).

### 0.1. Hacer que el terminal diga hola

Ejecute el siguiente comando en su terminal.

```bash
echo 'Hello World!'
```

??? success "Salida del comando"

    ```console
    Hello World!
    ```

Esto muestra el texto 'Hello World' directamente en el terminal.

### 0.2. Escribir la salida a un archivo

Ejecutar pipelines principalmente implica leer datos de archivos y escribir resultados en otros archivos, así que modifiquemos el comando para escribir la salida de texto en un archivo para hacer el ejemplo un poco más relevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Salida del comando"

    ```console

    ```

Esto no muestra nada en el terminal.

### 0.3. Encontrar la salida

El texto 'Hello World' ahora debería estar en el archivo de salida que especificamos, llamado `output.txt`.
Puede abrirlo en el explorador de archivos o desde la línea de comandos usando la utilidad `cat`, por ejemplo.

??? abstract "Contenido del archivo"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Esto es lo que vamos a intentar replicar con nuestro primer flujo de trabajo de Nextflow.

### Conclusión

Ahora sabe cómo ejecutar un comando simple en el terminal que produce algún texto, y opcionalmente, cómo hacer que escriba la salida en un archivo.

### ¿Qué sigue?

Descubra cómo se vería esto escrito como un flujo de trabajo de Nextflow.

---

## 1. Examinar el script y ejecutarlo

Le proporcionamos un script de flujo de trabajo completamente funcional, aunque minimalista, llamado `hello-world.nf` que hace lo mismo que antes (escribir 'Hello World!') pero con Nextflow.

Para comenzar, abramos el script del flujo de trabajo para que pueda tener una idea de cómo está estructurado.
Luego lo ejecutaremos y buscaremos sus salidas.

### 1.1. Examinar el código

Encontrará el script `hello-world.nf` en su directorio actual, que debería ser `hello-nextflow`. Ábralo en el panel del editor.

??? full-code "Archivo de código completo"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usar echo para imprimir 'Hello World!' a un archivo
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // emitir un saludo
        sayHello()
    }
    ```

Un script de flujo de trabajo de Nextflow típicamente incluye una o más definiciones de **process** y el **workflow** en sí, además de algunos bloques opcionales (no presentes aquí) que introduciremos más adelante.

Cada **process** describe qué operación(es) debe realizar el paso correspondiente en el pipeline, mientras que el **workflow** describe la lógica de flujo de datos que conecta los diversos pasos.

Vamos a examinar más de cerca el bloque **process** primero, luego veremos el bloque **workflow**.

#### 1.1.1. La definición de `process`

El primer bloque de código describe un **process**.

La definición del proceso comienza con la palabra clave `process`, seguida del nombre del proceso y finalmente el cuerpo del proceso delimitado por llaves.
El cuerpo del proceso debe contener un bloque script que especifica el comando a ejecutar, que puede ser cualquier cosa que podría ejecutar en una terminal de línea de comandos.

```groovy title="hello-world.nf" linenums="3"
/*
* Usar echo para imprimir 'Hello World!' a un archivo
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Aquí tenemos un **process** llamado `sayHello` que escribe su **output** en un archivo llamado `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Esta es una definición de proceso muy mínima que solo contiene una definición de `output` y el `script` a ejecutar.

La definición de `output` incluye el calificador `path`, que le dice a Nextflow que esto debe manejarse como una ruta (incluye tanto rutas de directorio como archivos).
Otro calificador común es `val`.

Importante: la definición de salida no _determina_ qué salida se creará.
Simplemente _declara_ cuál es la salida esperada, para que Nextflow pueda buscarla una vez que la ejecución esté completa.
Esto es necesario para verificar que el comando se ejecutó correctamente y para pasar la salida a procesos posteriores si es necesario. La salida producida que no coincida con lo declarado en el bloque de salida no se pasará a procesos posteriores.

!!! warning "Advertencia"

    Este ejemplo es frágil porque codificamos el nombre del archivo de salida en dos lugares separados (los bloques script y output).
    Si cambiamos uno pero no el otro, el script fallará.
    Más adelante, aprenderá formas de usar variables para mitigar este problema.

En un pipeline del mundo real, un proceso generalmente contiene bloques adicionales como directivas y entradas, que introduciremos en breve.

#### 1.1.2. La definición de `workflow`

El segundo bloque de código describe el **workflow** en sí.
La definición del flujo de trabajo comienza con la palabra clave `workflow`, seguida de un nombre opcional, luego el cuerpo del flujo de trabajo delimitado por llaves.

Aquí tenemos un **workflow** que consiste en un bloque `main:` (que dice 'este es el cuerpo principal del flujo de trabajo') que contiene una llamada al proceso `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emitir un saludo
    sayHello()
}
```

Esta es una definición de **workflow** muy mínima.
En un pipeline del mundo real, el flujo de trabajo típicamente contiene múltiples llamadas a **processes** conectados por **channels**, y los procesos esperan una o más **input(s)** variables.

Aprenderá cómo agregar entradas variables más adelante en este módulo de entrenamiento; y aprenderá cómo agregar más procesos y conectarlos mediante canales en la Parte 3 de este curso.

!!! tip "Consejo"

    Técnicamente, la línea `main:` no es requerida para flujos de trabajo simples como este, por lo que puede encontrar flujos de trabajo que no la tienen.
    Pero la necesitaremos para aprovechar las salidas a nivel de flujo de trabajo, así que es mejor incluirla desde el principio.

### 1.2. Ejecutar el flujo de trabajo

Mirar código no es ni de lejos tan divertido como ejecutarlo, así que probemos esto en la práctica.

#### 1.2.1. Iniciar el flujo de trabajo y monitorear la ejecución

En el terminal, ejecute el siguiente comando:

```bash
nextflow run hello-world.nf
```

??? success "Salida del comando"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Si la salida de su consola se parece a eso, ¡felicitaciones, acaba de ejecutar su primer flujo de trabajo de Nextflow!

La salida más importante aquí es la última línea, que está resaltada en la salida anterior:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Esto nos dice que el proceso `sayHello` se ejecutó exitosamente una vez (`1 of 1 ✔`).

Importante: esta línea también le dice dónde encontrar la salida de la llamada al proceso `sayHello`.
Veamos eso ahora.

#### 1.2.2. Encontrar la salida y los logs en el directorio `work`

Cuando ejecuta Nextflow por primera vez en un directorio dado, crea un directorio llamado `work` donde escribirá todos los archivos (y cualquier enlace simbólico) generados en el curso de la ejecución.

Dentro del directorio `work`, Nextflow organiza las salidas y logs por llamada de proceso.
Para cada llamada de proceso, Nextflow crea un subdirectorio anidado, nombrado con un hash para hacerlo único, donde preparará todas las entradas necesarias (usando enlaces simbólicos por defecto), escribirá archivos auxiliares, y escribirá logs y cualquier salida del proceso.

La ruta a ese subdirectorio se muestra en forma truncada entre corchetes en la salida de la consola.
Mirando lo que obtuvimos para la ejecución mostrada arriba, la línea de log de la consola para el proceso sayHello comienza con `[65/7be2fa]`. Eso corresponde a la siguiente ruta de directorio: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Veamos qué hay ahí.

??? abstract "Contenido del directorio"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "¿No ve lo mismo?"

    Los nombres exactos de los subdirectorios serán diferentes en su sistema.

    Si navega por los contenidos del subdirectorio de tarea en el explorador de archivos de VSCode, verá todos los archivos de inmediato.
    Sin embargo, los archivos de log están configurados para ser invisibles en el terminal, así que si quiere usar `ls` o `tree` para verlos, necesitará establecer la opción relevante para mostrar archivos invisibles.

    ```bash
    tree -a work
    ```

Lo primero que quiere ver es la salida real del flujo de trabajo, es decir, el archivo `output.txt` producido por el proceso `sayHello`.
Ábralo y encontrará el saludo `Hello World!`, que era el objetivo de nuestro flujo de trabajo minimalista.

??? abstract "Contenido del archivo"

    ```console title="output.txt"
    Hello World!
    ```

¡Funcionó!

Es cierto que puede parecer mucho código envolvente para un resultado tan pequeño, pero el valor de todo ese código envolvente se volverá más obvio una vez que empecemos a leer archivos de entrada y encadenar múltiples pasos.

Dicho esto, también veamos los otros archivos en ese directorio. Esos son archivos auxiliares y de log producidos por Nextflow como parte de la ejecución de la tarea.

- **`.command.begin`**: Metadatos relacionados con el inicio de la ejecución de la llamada al proceso
- **`.command.err`**: Mensajes de error (`stderr`) emitidos por la llamada al proceso
- **`.command.log`**: Salida de log completa emitida por la llamada al proceso
- **`.command.out`**: Salida regular (`stdout`) de la llamada al proceso
- **`.command.run`**: Script completo ejecutado por Nextflow para ejecutar la llamada al proceso
- **`.command.sh`**: El comando que fue realmente ejecutado por la llamada al proceso
- **`.exitcode`**: El código de salida resultante del comando

El archivo `.command.sh` es especialmente útil porque le dice el comando principal que Nextflow ejecutó, sin incluir toda la contabilidad y configuración de tarea/entorno.

??? abstract "Contenido del archivo"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Esto coincide con lo que ejecutamos manualmente antes.

En este caso es muy directo porque el comando del proceso estaba codificado de forma fija, pero más adelante en el curso verá comandos de proceso que involucran cierta interpolación de variables.
Eso hace especialmente valioso poder ver exactamente cómo Nextflow interpretó el código y qué comando se produjo cuando está solucionando problemas de una ejecución fallida.

### 1.3. Ejecutar el flujo de trabajo de nuevo

Intente volver a ejecutar el flujo de trabajo varias veces, luego mire los directorios de tareas bajo `work/`.

??? abstract "Contenido del directorio"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Ve que se ha creado un nuevo subdirectorio con un conjunto completo de archivos de salida y log para cada ejecución.
Esto le muestra que ejecutar el mismo flujo de trabajo varias veces no sobrescribirá los resultados de ejecuciones anteriores.

### Conclusión

Sabe cómo descifrar un script simple de Nextflow, ejecutarlo y encontrar la salida y los archivos de log relevantes en el directorio de trabajo.

### ¿Qué sigue?

Aprenda a publicar las salidas del flujo de trabajo en una ubicación más conveniente.

---

## 2. Publicar salidas

Como acaba de aprender, la salida producida por nuestro pipeline está enterrada en un directorio de trabajo varios niveles de profundidad.
Esto se hace a propósito; Nextflow tiene el control de este directorio y no debemos interactuar con él.
Sin embargo, eso hace inconveniente recuperar las salidas que nos importan.

Afortunadamente, Nextflow proporciona una forma de publicar salidas en un directorio designado usando [definiciones de salida a nivel de flujo de trabajo](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Uso básico

Esto va a involucrar dos nuevas piezas de código:

1. Un bloque `publish:` dentro del cuerpo del `workflow`, declarando las salidas del proceso.
2. Un bloque `output` en el script especificando opciones de salida como modo y ubicación.

#### 2.1.1. Declarar la salida del proceso `sayHello`

Necesitamos agregar un bloque `publish:` al cuerpo del flujo de trabajo (el mismo tipo de elemento de código que el bloque `main:`) y listar la salida del proceso `sayHello()`.

En el archivo de script del flujo de trabajo `hello-world.nf`, agregue las siguientes líneas de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emitir un saludo
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emitir un saludo
        sayHello()
    }
    ```

Ve que podemos referirnos a la salida del proceso simplemente haciendo `sayHello().out`, y asignarle un nombre arbitrario, `first_output`.

#### 2.1.2. Agregar un bloque `output:` al script

Ahora solo necesitamos agregar el bloque `output:` donde se especificará la ruta del directorio de salida. Note que este nuevo bloque se ubica **fuera** y **debajo** del bloque `workflow` dentro del script.

En el archivo de script del flujo de trabajo `hello-world.nf`, agregue las siguientes líneas de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emitir un saludo
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emitir un saludo
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Podemos usar esto para asignar rutas específicas a cualquier salida de proceso declarada en el bloque `workflow`.
Más adelante, aprenderá sobre formas de generar estructuras de directorio de salida sofisticadas, pero por ahora, simplemente estamos codificando una ruta mínima para simplificar.

#### 2.1.3. Ejecutar el flujo de trabajo

Ahora ejecute el script de flujo de trabajo modificado:

```bash
nextflow run hello-world.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

La salida del terminal debería parecer familiar. Externamente, nada ha cambiado.

Sin embargo, revise su explorador de archivos: esta vez, Nextflow ha creado un nuevo directorio llamado `results/`.

??? abstract "Contenido del directorio"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

Dentro del directorio `results`, encontramos un enlace simbólico al `output.txt` producido en el directorio de trabajo por el comando que acabamos de ejecutar.

Esto nos permite recuperar fácilmente los archivos de salida sin tener que buscar en el subdirectorio de trabajo.

### 2.2. Establecer una ubicación personalizada

Tener una ubicación predeterminada es genial, pero puede querer personalizar dónde se guardan los resultados y cómo se organizan.

Por ejemplo, puede querer organizar sus salidas en subdirectorios.
La forma más simple de hacer eso es asignar una ruta de salida específica por salida.

#### 2.2.1. Modificar la ruta de salida

Una vez más, modificar el comportamiento de publicación para una salida específica es realmente sencillo.
Para establecer una ubicación personalizada, simplemente edite el `path` de acuerdo:

=== "Después"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Como esto se establece a nivel de la salida individual, puede especificar diferentes ubicaciones y subdirectorios para adaptarse a sus necesidades.

#### 2.2.2. Ejecutar el flujo de trabajo de nuevo

Probémoslo.

```bash
nextflow run hello-world.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Esta vez el resultado se escribe bajo el subdirectorio especificado.

??? abstract "Contenido del directorio"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Ve que el resultado de la ejecución anterior todavía está ahí.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Puede usar tantos niveles de anidamiento como desee.
También es posible usar el nombre del proceso u otras variables para nombrar los directorios usados para organizar los resultados, y es posible cambiar el nombre predeterminado del directorio de salida de nivel superior (que es controlado por la variable especial `outputDir`).
Cubriremos estas opciones en entrenamientos posteriores.

### 2.3. Establecer el modo de publicación a copia

Por defecto, las salidas se publican como enlaces simbólicos desde el directorio `work`.
Eso significa que solo hay un único archivo en el sistema de archivos.

Esto es genial cuando está tratando con archivos muy grandes, para los cuales no quiere almacenar múltiples copias.
Sin embargo, si elimina el directorio de trabajo en algún momento (cubriremos las operaciones de limpieza en breve), perderá acceso al archivo.
Así que necesita tener un plan para guardar copias de cualquier archivo importante en un lugar seguro.

Una opción fácil es cambiar el modo de publicación a copia para las salidas que le importan.

#### 2.3.1. Agregar la directiva de modo

Esta parte es realmente sencilla.
Simplemente agregue `mode 'copy'` a la definición de salida a nivel de flujo de trabajo relevante:

=== "Después"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Esto establece el modo de publicación para esa salida específica.

#### 2.3.2. Ejecutar el flujo de trabajo de nuevo

Probémoslo.

```bash
nextflow run hello-world.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Esta vez, si mira los resultados, el archivo es una copia real en lugar de solo un enlace simbólico.

??? abstract "Contenido del directorio"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Como esto también se establece a nivel de la salida individual, le permite establecer el modo de publicación de manera granular.
Esto será especialmente útil más adelante cuando pasemos a pipelines de múltiples pasos, donde puede querer solo copiar las salidas finales y dejar las salidas intermedias como enlaces simbólicos, por ejemplo.

Como se señaló antes, hay otras opciones más sofisticadas para controlar cómo se publican las salidas.
Le mostraremos cómo usarlas a su debido tiempo en su viaje con Nextflow.

### 2.4. Nota sobre directivas `publishDir` a nivel de proceso

Hasta hace muy poco, la forma establecida de publicar salidas era hacerlo a nivel de cada proceso individual usando una directiva `publishDir`.

Para lograr lo que acabamos de hacer para las salidas del proceso `sayHello`, habríamos agregado la siguiente línea a la definición del proceso:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Todavía encontrará este patrón de código en todos los pipelines más antiguos de Nextflow y módulos de proceso, por lo que es importante estar al tanto de él.
Sin embargo, no recomendamos usarlo en ningún trabajo nuevo ya que eventualmente se deshabilitará en futuras versiones del lenguaje Nextflow.

### Conclusión

Sabe cómo publicar salidas del flujo de trabajo en una ubicación más conveniente.

### ¿Qué sigue?

Aprenda a proporcionar una entrada variable a través de un parámetro de línea de comandos y utilizar valores predeterminados de manera efectiva.

---

## 3. Usar una entrada variable pasada en la línea de comandos

En su estado actual, nuestro flujo de trabajo usa un saludo codificado en el comando del proceso.
Queremos agregar algo de flexibilidad usando una variable de entrada, para poder cambiar más fácilmente el saludo en tiempo de ejecución.

Esto requiere que hagamos tres conjuntos de cambios en nuestro script:

1. Cambiar el proceso para esperar una entrada variable
2. Configurar un parámetro de línea de comandos para capturar la entrada del usuario
3. Pasar la entrada al proceso en el cuerpo del flujo de trabajo

Hagamos estos cambios uno a la vez.

### 3.1. Cambiar el proceso `sayHello` para esperar una entrada variable

Necesitamos editar la definición del proceso para (1) aceptar una variable de entrada y (2) usar esa variable en la línea de comandos.

#### 3.1.1. Agregar un bloque de entrada a la definición del proceso

Primero, adaptemos la definición del proceso para aceptar una entrada llamada `greeting`.

En el bloque del proceso, haga el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

La variable `greeting` está prefijada con `val` para decirle a Nextflow que es un valor (no una ruta).

#### 3.1.2. Editar el comando del proceso para usar la variable de entrada

Ahora intercambiamos el valor original codificado por el valor de la variable de entrada que esperamos recibir.

En el bloque del proceso, haga el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

El símbolo `$` y las llaves (`{ }`) le dicen a Nextflow que este es un nombre de variable que necesita ser reemplazado con el valor de entrada real (=interpolado).

!!! tip "Consejo"

    Las llaves (`{ }`) eran técnicamente opcionales en versiones anteriores de Nextflow, así que puede ver flujos de trabajo más antiguos donde esto está escrito como `echo '$greeting' > output.txt`.

Ahora que el proceso `sayHello()` está listo para aceptar una entrada variable, necesitamos una forma de proporcionar un valor de entrada a la llamada del proceso a nivel de flujo de trabajo.

### 3.2. Configurar un parámetro de línea de comandos para capturar la entrada del usuario

Podríamos simplemente codificar una entrada directamente haciendo la llamada al proceso `sayHello('Hello World!')`.
Sin embargo, cuando estamos haciendo trabajo real con nuestro flujo de trabajo, querremos poder controlar sus entradas desde la línea de comandos.

Buenas noticias: Nextflow tiene un sistema de parámetros de flujo de trabajo incorporado llamado `params`, que facilita declarar y usar parámetros CLI.

La sintaxis general es declarar `params.<nombre_parametro>` para decirle a Nextflow que espere un parámetro `--<nombre_parametro>` en la línea de comandos.

Aquí, queremos crear un parámetro llamado `--input`, así que necesitamos declarar `params.input` en algún lugar del flujo de trabajo.
En principio podemos escribirlo en cualquier lugar; pero como vamos a querer dárselo a la llamada del proceso `sayHello()`, podemos conectarlo directamente escribiendo `sayHello(params.input)`.

En el bloque del flujo de trabajo, haga el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emitir un saludo
    sayHello(params.input)
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emitir un saludo
    sayHello()
    ```

Esto le dice a Nextflow que ejecute el proceso `sayHello` con el valor proporcionado a través del parámetro `--input`.

En efecto, hemos logrado los pasos (2) y (3) descritos al inicio de la sección de una sola vez.

### 3.3. Ejecutar el comando del flujo de trabajo

¡Ejecutémoslo!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Si hizo todas estas ediciones correctamente, debería obtener otra ejecución exitosa.

Asegúrese de abrir el archivo de salida para verificar que ahora tiene la nueva versión del saludo.

??? abstract "Contenido del archivo"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

¡Voilà!

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Note cómo la nueva ejecución ha sobrescrito el archivo de salida publicado en el directorio `results`.
Sin embargo, los resultados de las ejecuciones anteriores todavía se conservan en los directorios de tareas bajo `work`.

!!! tip "Consejo"

    Puede distinguir fácilmente los parámetros a nivel de Nextflow de los parámetros a nivel de pipeline.

    - Los parámetros que se aplican a un pipeline siempre llevan un doble guión (`--`).
    - Los parámetros que modifican una configuración de Nextflow, _por ejemplo_ la función `-resume` que usamos antes, llevan un solo guión (`-`).

### 3.4. Usar valores predeterminados para parámetros de línea de comandos

Ok, eso fue conveniente, pero en muchos casos, tiene sentido proporcionar un valor predeterminado para un parámetro dado para no tener que especificarlo en cada ejecución.

#### 3.4.1. Establecer un valor predeterminado para el parámetro CLI

Démosle al parámetro `input` un valor predeterminado declarándolo antes de la definición del flujo de trabajo.

```groovy title="hello-world.nf" linenums="20"
/*
 * Parámetros del pipeline
 */
params {
    input: String = 'Holà mundo!'
}
```

Como ve, podemos especificar el tipo de entrada que el flujo de trabajo espera (Nextflow 25.10.2 y posterior).
La sintaxis es `nombre: Tipo = valor_predeterminado`.
Los tipos soportados incluyen `String`, `Integer`, `Float`, `Boolean` y `Path`.

!!! info "Información"

    En flujos de trabajo más antiguos, puede ver todo ese bloque `params` escrito simplemente como `input = 'Holà mundo!'`.

A medida que agregue más parámetros a su pipeline, debería agregarlos todos a este bloque, ya sea que necesite o no darles un valor predeterminado.
Esto facilitará encontrar todos los parámetros configurables de un vistazo.

#### 3.4.2. Ejecutar el flujo de trabajo de nuevo sin especificar el parámetro

Ahora que tiene un valor predeterminado establecido, puede ejecutar el flujo de trabajo de nuevo sin tener que especificar un valor en la línea de comandos.

```bash
nextflow run hello-world.nf
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

La salida estará en el mismo lugar que antes, pero los contenidos deberían actualizarse con el nuevo texto.

??? abstract "Contenido del archivo"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow usó el valor predeterminado del parámetro greeting para crear la salida.

#### 3.4.3. Sobrescribir el valor predeterminado

Si proporciona el parámetro en la línea de comandos, el valor CLI sobrescribirá el valor predeterminado.

Pruébelo:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Una vez más, debería encontrar la salida actualizada correspondiente en su directorio de resultados.

??? abstract "Contenido del archivo"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Nota"

    En Nextflow, hay múltiples lugares donde puede especificar valores para parámetros.
    Si el mismo parámetro se establece con diferentes valores en múltiples lugares, Nextflow determinará qué valor usar basándose en el orden de precedencia que se describe [aquí](https://www.nextflow.io/docs/latest/config.html).

    Cubriremos esto en más detalle en la Parte 6 (Configuration).

### Conclusión

Sabe cómo usar una entrada variable simple proporcionada en tiempo de ejecución a través de un parámetro de línea de comandos, así como configurar, usar y sobrescribir valores predeterminados.

### ¿Qué sigue?

Aprenda a gestionar ejecuciones de manera más conveniente.

---

## 4. Gestionar ejecuciones del flujo de trabajo

Saber cómo iniciar flujos de trabajo y recuperar salidas es genial, pero encontrará rápidamente que hay algunos otros aspectos de la gestión de flujos de trabajo que harán su vida más fácil, especialmente si está desarrollando sus propios flujos de trabajo.

Aquí le mostramos cómo usar la función `resume` cuando necesita volver a iniciar el mismo flujo de trabajo, cómo inspeccionar el log de ejecuciones pasadas con `nextflow log`, y cómo eliminar directorios de trabajo antiguos con `nextflow clean`.

<!-- ¿Alguna otra opción interesante que deberíamos incluir? Se agregó log -->

### 4.1. Volver a iniciar un flujo de trabajo con `-resume`

A veces, querrá volver a ejecutar un pipeline que ya ha iniciado previamente sin rehacer ningún paso que ya se completó exitosamente.

Nextflow tiene una opción llamada `-resume` que le permite hacer esto.
Específicamente, en este modo, cualquier proceso que ya se haya ejecutado con exactamente el mismo código, configuración y entradas se omitirá.
Esto significa que Nextflow solo ejecutará procesos que haya agregado o modificado desde la última ejecución, o a los cuales está proporcionando nuevas configuraciones o entradas.

Hay dos ventajas clave de hacer esto:

- Si está en medio del desarrollo de su pipeline, puede iterar más rápidamente ya que solo tiene que ejecutar el(los) proceso(s) en los que está trabajando activamente para probar sus cambios.
- Si está ejecutando un pipeline en producción y algo sale mal, en muchos casos puede solucionar el problema y volver a iniciar el pipeline, y reanudará la ejecución desde el punto de fallo, lo que puede ahorrarle mucho tiempo y cómputo.

Para usarlo, simplemente agregue `-resume` a su comando y ejecútelo:

```bash
nextflow run hello-world.nf -resume
```

??? success "Salida del comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

La salida de la consola debería parecer familiar, pero hay algo que es un poco diferente comparado con antes.

Busque la parte `cached:` que se ha agregado en la línea de estado del proceso (línea 5), que significa que Nextflow ha reconocido que ya hizo este trabajo y simplemente reutilizó el resultado de la ejecución exitosa anterior.

También puede ver que el hash del subdirectorio de trabajo es el mismo que en la ejecución anterior.
Nextflow literalmente le está señalando la ejecución anterior y diciendo "Ya hice eso allí."

!!! tip "Consejo"

    Cuando vuelve a ejecutar un pipeline con `resume`, Nextflow no sobrescribe ningún archivo publicado fuera del directorio de trabajo por ninguna ejecución que se ejecutó exitosamente anteriormente.

### 4.2. Inspeccionar el log de ejecuciones pasadas

Ya sea que esté desarrollando un nuevo pipeline o ejecutando pipelines en producción, en algún momento probablemente necesitará buscar información sobre ejecuciones pasadas.
Aquí está cómo hacerlo.

Cada vez que inicia un flujo de trabajo de Nextflow, se escribe una línea en un archivo de log llamado `history`, bajo un directorio oculto llamado `.nextflow` en el directorio de trabajo actual.

??? abstract "Contenido del archivo"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Este archivo le da la marca de tiempo, nombre de ejecución, estado, ID de revisión, ID de sesión y línea de comando completa para cada ejecución de Nextflow que se ha iniciado desde dentro del directorio de trabajo actual.

Una forma más conveniente de acceder a esta información es usar el comando `nextflow log`.

```bash
nextflow log
```

??? success "Salida del comando"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Esto mostrará los contenidos del archivo de log en el terminal, aumentado con una línea de encabezado.

Notará que el ID de sesión cambia cada vez que ejecuta un nuevo comando `nextflow run`, EXCEPTO si está usando la opción `-resume`.
En ese caso, el ID de sesión permanece igual.

Nextflow usa el ID de sesión para agrupar información de caché de ejecución bajo el directorio `cache`, también ubicado bajo `.nextflow`.

### 4.3. Eliminar directorios de trabajo antiguos

Durante el proceso de desarrollo, típicamente ejecutará su borrador de pipeline un gran número de veces, lo que puede llevar a una acumulación de muchos archivos en muchos subdirectorios.

Afortunadamente Nextflow incluye un útil subcomando `clean` que puede eliminar automáticamente los subdirectorios de trabajo para ejecuciones pasadas que ya no le importan.

#### 4.3.1. Determinar criterios de eliminación

Hay múltiples [opciones](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para determinar qué eliminar.

Aquí le mostramos un ejemplo que elimina todos los subdirectorios de ejecuciones antes de una ejecución dada, especificada usando su nombre de ejecución.

Busque la ejecución exitosa más reciente donde no usó `-resume`; en nuestro caso el nombre de ejecución fue `golden_cantor`.

El nombre de ejecución es la cadena de dos partes generada por la máquina que se muestra entre corchetes en la línea de salida de consola `Launching (...)`.
También puede usar el log de Nextflow para buscar una ejecución basada en su marca de tiempo y/o línea de comando.

#### 4.3.2. Hacer una ejecución de prueba

Primero usamos el flag de ejecución de prueba `-n` para verificar qué se eliminará dado el comando:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Salida del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Su salida tendrá diferentes nombres de directorio de tarea y puede tener un número diferente de líneas, pero debería parecer similar al ejemplo.

Si no ve ninguna línea de salida, o no proporcionó un nombre de ejecución válido o no hay ejecuciones pasadas para eliminar. Asegúrese de cambiar `golden_cantor` en el comando de ejemplo por el correspondiente nombre de ejecución más reciente en su log.

#### 4.3.3. Proceder con la eliminación

Si la salida parece como se esperaba y quiere proceder con la eliminación, vuelva a ejecutar el comando con el flag `-f` en lugar de `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Salida del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

La salida debería ser similar a antes, pero ahora diciendo 'Removed' en lugar de 'Would remove'.
Note que esto no elimina los subdirectorios de dos caracteres (como `a3/` arriba) pero sí vacía sus contenidos.

!!! Warning "Advertencia"

    Eliminar subdirectorios de trabajo de ejecuciones pasadas los elimina del caché de Nextflow y elimina cualquier salida que se almacenó en esos directorios.
    Eso significa que rompe la capacidad de Nextflow de reanudar la ejecución sin volver a ejecutar los procesos correspondientes.

    ¡Usted es responsable de guardar cualquier salida que le importe o en la que planee confiar! Esa es la razón principal por la que preferimos usar el modo `copy` en lugar del modo `symlink` para la directiva `publish`.

### Conclusión

Sabe cómo publicar salidas en un directorio específico, volver a iniciar un pipeline sin repetir pasos que ya se ejecutaron de manera idéntica, y usar el comando `nextflow clean` para limpiar directorios de trabajo antiguos.

Más generalmente, sabe cómo interpretar un flujo de trabajo simple de Nextflow, gestionar su ejecución y recuperar salidas.

### ¿Qué sigue?

¡Tómese un pequeño descanso, se lo ha ganado!

Cuando esté listo, continúe a [**Parte 2: Hello Channels**](./02_hello_channels.md) para aprender cómo usar canales para alimentar entradas a su flujo de trabajo, lo que le permitirá aprovechar el paralelismo de flujo de datos integrado de Nextflow y otras características poderosas.

---

## Cuestionario

<quiz>
¿Cuáles son los componentes mínimos requeridos de un proceso de Nextflow?
- [ ] Solo bloques de entrada y salida
- [x] Bloques de salida y script
- [ ] Bloques de entrada, salida y script
- [ ] Solo un bloque de script

Más información: [1.1.1. La definición de process](#111-la-definicion-de-process)
</quiz>

<quiz>
¿Cuál es el propósito del bloque output en un proceso?
- [ ] Imprimir resultados en la consola
- [ ] Guardar archivos en el directorio de trabajo
- [x] Declarar las salidas esperadas del proceso
- [ ] Definir variables de entorno

Más información: [1.1.1. La definición de process](#111-la-definicion-de-process)
</quiz>

<quiz>
¿Qué comando se usa para ejecutar un flujo de trabajo de Nextflow?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Mirando el directorio de trabajo de una tarea, ¿qué archivo contiene el comando real que se ejecutó?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Más información: [1.2.2. Encontrar la salida y los logs en el directorio `work`](#122-encontrar-la-salida-y-los-logs-en-el-directorio-work)
</quiz>

<quiz>
¿Qué hace el flag `-resume`?
- [ ] Reinicia el flujo de trabajo desde el principio
- [ ] Pausa el flujo de trabajo
- [x] Omite procesos que ya se completaron exitosamente
- [ ] Crea una copia de seguridad del flujo de trabajo

Más información: [4.1. Volver a iniciar un flujo de trabajo con `-resume`](#41-volver-a-iniciar-un-flujo-de-trabajo-con--resume)
</quiz>

<quiz>
¿Cuál es el modo predeterminado para publicar salidas del flujo de trabajo?
- [ ] Copiar archivos al directorio de salida
- [x] Crear enlaces simbólicos en el directorio de salida
- [ ] Mover archivos al directorio de salida
- [ ] Comprimir archivos en el directorio de salida

Más información: [2.3. Establecer el modo de publicación a copia](#23-establecer-el-modo-de-publicacion-a-copia)
</quiz>

<quiz>
¿Cómo se pasa un valor de parámetro a un flujo de trabajo de Nextflow desde la línea de comandos?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Más información: [3.2. Configurar un parámetro de línea de comandos para capturar la entrada del usuario](#32-configurar-un-parametro-de-linea-de-comandos-para-capturar-la-entrada-del-usuario)
</quiz>

<quiz>
¿Cómo se referencia una variable dentro de un bloque script de Nextflow?
- [ ] Usar sintaxis `%variable%`
- [x] Usar sintaxis `#!groovy ${variable}`
- [ ] Usar sintaxis `{{variable}}`
- [ ] Usar sintaxis `[variable]`
</quiz>
