# Parte 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vea [la lista de reproducción completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/01_hello_world.md).
///

En esta primera parte del curso de capacitación Hello Nextflow, comenzamos con un ejemplo muy básico de Hello World independiente del dominio, que desarrollaremos progresivamente para demostrar el uso de la lógica y los componentes fundamentales de Nextflow.

??? info "¿Qué es un ejemplo Hello World?"

    Un "Hello World!" es un ejemplo minimalista diseñado para demostrar la sintaxis básica y la estructura de un lenguaje de programación o framework de software.
    El ejemplo típicamente consiste en imprimir la frase "Hello, World!" en el dispositivo de salida, como la consola o terminal, o escribirla en un archivo.

---

## 0. Calentamiento: Ejecute un ejemplo Hello World directamente

Demostremos esto con un comando simple que ejecutamos directamente en la terminal, para mostrar qué hace antes de envolverlo en Nextflow.

!!! tip "Consejo"

    Recuerde que ahora debe estar dentro del directorio `hello-nextflow/` como se describe en la página [Primeros pasos](00_orientation.md).

### 0.1. Haga que la terminal diga hola

Ejecute el siguiente comando en su terminal.

```bash
echo 'Hello World!'
```

??? success "Salida del comando"

    ```console
    Hello World!
    ```

Esto muestra el texto 'Hello World' directamente en la terminal.

### 0.2. Escriba la salida en un archivo

Ejecutar pipelines implica principalmente leer datos de archivos y escribir resultados en otros archivos, así que modifiquemos el comando para escribir la salida de texto en un archivo para hacer el ejemplo un poco más relevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Salida del comando"

    ```console

    ```

Esto no muestra nada en la terminal.

### 0.3. Encuentre la salida

El texto 'Hello World' ahora debería estar en el archivo de salida que especificamos, llamado `output.txt`.
Puede abrirlo en el explorador de archivos o desde la línea de comandos usando la utilidad `cat`, por ejemplo.

??? abstract "Contenido del archivo"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Esto es lo que vamos a intentar replicar con nuestro primer workflow de Nextflow.

### Conclusión

Ahora sabe cómo ejecutar un comando simple en la terminal que genera texto y, opcionalmente, cómo hacer que escriba la salida en un archivo.

### ¿Qué sigue?

Descubra cómo se vería esto escrito como un workflow de Nextflow.

---

## 1. Examine el script y ejecútelo

Le proporcionamos un script de workflow completamente funcional, aunque minimalista, llamado `hello-world.nf` que hace lo mismo que antes (escribir 'Hello World!') pero con Nextflow.

Para comenzar, abramos el script del workflow para que pueda tener una idea de cómo está estructurado.
Luego lo ejecutaremos y buscaremos sus salidas.

### 1.1. Examine el código

Encontrará el script `hello-world.nf` en su directorio actual, que debería ser `hello-nextflow`. Ábralo en el panel del editor.

??? full-code "Archivo de código completo"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo para imprimir 'Hello World!' en un archivo
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
        // emite un saludo
        sayHello()
    }
    ```

Un script de workflow de Nextflow típicamente incluye una o más definiciones de [**process**](https://nextflow.io/docs/latest/process.html) y el [**workflow**](https://nextflow.io/docs/latest/workflow.html) en sí, además de algunos bloques opcionales (no presentes aquí) que presentaremos más adelante.

Cada **process** describe qué operación(es) debe realizar el paso correspondiente en el pipeline, mientras que el **workflow** describe la lógica de flujo de datos que conecta los diversos pasos.

Vamos a examinar más de cerca el bloque **process** primero, luego veremos el bloque **workflow**.

#### 1.1.1. La definición del `process`

El primer bloque de código describe un **process**.

La definición del proceso comienza con la palabra clave `process`, seguida del nombre del proceso y finalmente el cuerpo del proceso delimitado por llaves.
El cuerpo del proceso debe contener un bloque script que especifica el comando a ejecutar, que puede ser cualquier cosa que pueda ejecutar en una terminal de línea de comandos.

```groovy title="hello-world.nf" linenums="3"
/*
* Use echo para imprimir 'Hello World!' en un archivo
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

Es importante destacar que la definición de salida no _determina_ qué salida se creará.
Simplemente _declara_ cuál es la salida esperada, para que Nextflow pueda buscarla una vez que se complete la ejecución.
Esto es necesario para verificar que el comando se ejecutó correctamente y para pasar la salida a procesos posteriores si es necesario. La salida producida que no coincida con lo declarado en el bloque de salida no se pasará a procesos posteriores.

!!! warning "Advertencia"

    Este ejemplo es frágil porque codificamos el nombre del archivo de salida en dos lugares separados (los bloques script y output).
    Si cambiamos uno pero no el otro, el script fallará.
    Más adelante, aprenderá formas de usar variables para mitigar este problema.

En un pipeline del mundo real, un proceso generalmente contiene bloques adicionales como directivas y entradas, que presentaremos en un momento.

#### 1.1.2. La definición del `workflow`

El segundo bloque de código describe el **workflow** en sí.
La definición del workflow comienza con la palabra clave `workflow`, seguida de un nombre opcional, luego el cuerpo del workflow delimitado por llaves.

Aquí tenemos un **workflow** que consiste en un bloque `main:` (que dice 'este es el cuerpo principal del workflow') que contiene una llamada al proceso `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // emite un saludo
    sayHello()
}
```

Esta es una definición de **workflow** muy mínima.
En un pipeline del mundo real, el workflow típicamente contiene múltiples llamadas a **procesos** conectados por **canales**, y los procesos esperan una o más **entrada(s)** variables.

Aprenderá cómo agregar entradas variables más adelante en este módulo de capacitación; y aprenderá cómo agregar más procesos y conectarlos mediante canales en la Parte 3 de este curso.

!!! tip "Consejo"

    Técnicamente, la línea `main:` no es necesaria para workflows simples como este, por lo que puede encontrar workflows que no la tienen.
    Pero la necesitaremos para aprovechar las salidas a nivel de workflow, así que es mejor incluirla desde el principio.

### 1.2. Ejecute el workflow

Mirar código no es tan divertido como ejecutarlo, así que probemos esto en la práctica.

#### 1.2.1. Lance el workflow y monitoree la ejecución

En la terminal, ejecute el siguiente comando:

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

Si la salida de su consola se ve algo así, entonces felicidades, ¡acaba de ejecutar su primer workflow de Nextflow!

La salida más importante aquí es la última línea, que está resaltada en la salida anterior:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Esto nos dice que el proceso `sayHello` se ejecutó exitosamente una vez (`1 of 1 ✔`).

Es importante destacar que esta línea también le indica dónde encontrar la salida de la llamada al proceso `sayHello`.
Veamos eso ahora.

#### 1.2.2. Encuentre la salida y los logs en el directorio `work`

Cuando ejecuta Nextflow por primera vez en un directorio dado, crea un directorio llamado `work` donde escribirá todos los archivos (y cualquier enlace simbólico) generados en el curso de la ejecución.

Dentro del directorio `work`, Nextflow organiza las salidas y logs por llamada de proceso.
Para cada llamada de proceso, Nextflow crea un subdirectorio anidado, nombrado con un hash para hacerlo único, donde preparará todas las entradas necesarias (usando enlaces simbólicos por defecto), escribirá archivos auxiliares y escribirá logs y cualquier salida del proceso.

La ruta a ese subdirectorio se muestra en forma truncada entre corchetes en la salida de la consola.
Mirando lo que obtuvimos para la ejecución mostrada arriba, la línea de log de la consola para el proceso sayHello comienza con `[65/7be2fa]`. Eso corresponde a la siguiente ruta de directorio: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Echemos un vistazo a lo que hay allí.

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

    Si navega por el contenido del subdirectorio de la tarea en el explorador de archivos de VSCode, verá todos los archivos de inmediato.
    Sin embargo, los archivos de log están configurados para ser invisibles en la terminal, por lo que si desea usar `ls` o `tree` para verlos, deberá establecer la opción relevante para mostrar archivos invisibles.

    ```bash
    tree -a work
    ```

Lo primero que querrá ver es la salida real del workflow, es decir, el archivo `output.txt` producido por el proceso `sayHello`.
Ábralo y encontrará el saludo `Hello World!`, que era el objetivo de nuestro workflow minimalista.

??? abstract "Contenido del archivo"

    ```console title="output.txt"
    Hello World!
    ```

¡Funcionó!

Ciertamente, puede parecer mucho código envolvente para un resultado tan pequeño, pero el valor de todo ese código envolvente se volverá más obvio una vez que comencemos a leer archivos de entrada y encadenar múltiples pasos.

Dicho esto, veamos también los otros archivos en ese directorio. Esos son archivos auxiliares y de log producidos por Nextflow como parte de la ejecución de la tarea.

- **`.command.begin`**: Metadatos relacionados con el inicio de la ejecución de la llamada al proceso
- **`.command.err`**: Mensajes de error (`stderr`) emitidos por la llamada al proceso
- **`.command.log`**: Salida de log completa emitida por la llamada al proceso
- **`.command.out`**: Salida regular (`stdout`) de la llamada al proceso
- **`.command.run`**: Script completo ejecutado por Nextflow para ejecutar la llamada al proceso
- **`.command.sh`**: El comando que realmente fue ejecutado por la llamada al proceso
- **`.exitcode`**: El código de salida resultante del comando

El archivo `.command.sh` es especialmente útil porque le dice el comando principal que Nextflow ejecutó, sin incluir toda la contabilidad y configuración de tarea/entorno.

??? abstract "Contenido del archivo"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Esto coincide con lo que ejecutamos anteriormente manualmente.

En este caso es muy sencillo porque el comando del proceso estaba codificado, pero más adelante en el curso verá comandos de proceso que involucran alguna interpolación de variables.
Eso hace que sea especialmente valioso poder ver exactamente cómo Nextflow interpretó el código y qué comando se produjo cuando está solucionando problemas de una ejecución fallida.

### 1.3. Ejecute el workflow nuevamente

Intente volver a ejecutar el workflow varias veces, luego mire los directorios de tareas bajo `work/`.

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

Verá que se ha creado un nuevo subdirectorio con un conjunto completo de archivos de salida y log para cada ejecución.
Esto le muestra que ejecutar el mismo workflow varias veces no sobrescribirá los resultados de ejecuciones anteriores.

### Conclusión

Sabe cómo descifrar un script simple de Nextflow, ejecutarlo y encontrar la salida y los archivos de log relevantes en el directorio work.

### ¿Qué sigue?

Aprenda cómo publicar las salidas del workflow en una ubicación más conveniente.

---

## 2. Publique salidas

Como acaba de aprender, la salida producida por nuestro pipeline está enterrada en un directorio de trabajo varias capas más abajo.
Esto se hace a propósito; Nextflow tiene el control de este directorio y no se supone que interactuemos con él.
Sin embargo, eso hace que sea inconveniente recuperar salidas que nos interesan.

Afortunadamente, Nextflow proporciona una forma de publicar salidas en un directorio designado usando [definiciones de salida de workflow](https://nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Uso básico

Esto va a involucrar dos nuevas piezas de código:

1. Un bloque `publish:` dentro del cuerpo del `workflow`, declarando salidas de proceso.
2. Un bloque `output` en el script especificando opciones de salida como modo y ubicación.

#### 2.1.1. Declare la salida del proceso `sayHello`

Necesitamos agregar un bloque `publish:` al cuerpo del workflow (el mismo tipo de elemento de código que el bloque `main:`) y listar la salida del proceso `sayHello()`.

En el archivo de script del workflow `hello-world.nf`, agregue las siguientes líneas de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // emite un saludo
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // emite un saludo
        sayHello()
    }
    ```

Verá que podemos referirnos a la salida del proceso simplemente haciendo `sayHello().out`, y asignarle un nombre arbitrario, `first_output`.

#### 2.1.2. Agregue un bloque `output:` al script

Ahora solo necesitamos agregar el bloque `output:` donde se especificará la ruta del directorio de salida. Tenga en cuenta que este nuevo bloque se encuentra **fuera** y **debajo** del bloque `workflow` dentro del script.

En el archivo de script del workflow `hello-world.nf`, agregue las siguientes líneas de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // emite un saludo
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
        // emite un saludo
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Podemos usar esto para asignar rutas específicas a cualquier salida de proceso declarada en el bloque `workflow`.
Más adelante, aprenderá sobre formas de generar estructuras de directorios de salida sofisticadas, pero por ahora, solo estamos codificando una ruta mínima por simplicidad.

#### 2.1.3. Ejecute el workflow

Ahora ejecute el script de workflow modificado:

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

La salida de la terminal debería verse familiar. Externamente, nada ha cambiado.

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

Dentro del directorio `results`, encontramos un enlace simbólico al `output.txt` producido en el directorio work por el comando que acabamos de ejecutar.

Esto nos permite recuperar fácilmente archivos de salida sin tener que buscar en el subdirectorio work.

### 2.2. Establezca una ubicación personalizada

Tener una ubicación predeterminada es excelente, pero es posible que desee personalizar dónde se guardan los resultados y cómo se organizan.

Por ejemplo, es posible que desee organizar sus salidas en subdirectorios.
La forma más simple de hacer eso es asignar una ruta de salida específica por salida.

#### 2.2.1. Modifique la ruta de salida

Una vez más, modificar el comportamiento de publicación para una salida específica es realmente sencillo.
Para establecer una ubicación personalizada, simplemente edite el `path` en consecuencia:

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

Dado que esto se establece a nivel de la salida individual, puede especificar diferentes ubicaciones y subdirectorios según sus necesidades.

#### 2.2.2. Ejecute el workflow nuevamente

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

Verá que el resultado de la ejecución anterior todavía está allí.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Puede usar tantos niveles de anidamiento como desee.
También es posible usar el nombre del proceso u otras variables para nombrar los directorios utilizados para organizar resultados, y es posible cambiar el nombre predeterminado del directorio de salida de nivel superior (que está controlado por el flag CLI `-o` o la variable de configuración `outputDir`).
Cubriremos estas opciones más adelante en la capacitación.

### 2.3. Establezca el modo de publicación en copy

Por defecto, las salidas se publican como enlaces simbólicos desde el directorio `work`.
Eso significa que solo hay un único archivo en el sistema de archivos.

Esto es excelente cuando está tratando con archivos muy grandes, para los cuales no desea almacenar múltiples copias.
Sin embargo, si elimina el directorio work en algún momento (cubriremos las operaciones de limpieza en breve), perderá el acceso al archivo.
Por lo tanto, necesita tener un plan para guardar copias de cualquier archivo importante en un lugar seguro.

Una opción fácil es cambiar el modo de publicación a copy para las salidas que le interesan.

#### 2.3.1. Agregue la directiva mode

Esta parte es realmente sencilla.
Simplemente agregue `mode 'copy'` a la definición de salida a nivel de workflow relevante:

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

#### 2.3.2. Ejecute el workflow nuevamente

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

Esta vez, si observa los resultados, el archivo es una copia adecuada en lugar de solo un enlace simbólico.

??? abstract "Contenido del directorio"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Dado que esto también se establece a nivel de la salida individual, le permite establecer el modo de publicación de manera granular.
Esto será especialmente útil más adelante cuando pasemos a pipelines de múltiples pasos, donde es posible que desee copiar solo las salidas finales y dejar las salidas intermedias como enlaces simbólicos, por ejemplo.

Como se señaló anteriormente, hay otras opciones más sofisticadas para controlar cómo se publican las salidas.
Le mostraremos cómo usarlas a su debido tiempo en su viaje con Nextflow.

### 2.4. Nota sobre las directivas `publishDir` a nivel de proceso

Hasta hace muy poco, la forma establecida de publicar salidas era hacerlo a nivel de cada proceso individual usando una directiva `publishDir`.

Para lograr lo que acabamos de hacer para las salidas del proceso `sayHello`, habríamos agregado en su lugar la siguiente línea a la definición del proceso:

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

Todavía encontrará este patrón de código en todas partes en pipelines de Nextflow más antiguos y módulos de proceso, por lo que es importante estar al tanto de ello.
Sin embargo, no recomendamos usarlo en ningún trabajo nuevo, ya que eventualmente no se permitirá en futuras versiones del lenguaje Nextflow.

### Conclusión

Sabe cómo publicar salidas de workflow en una ubicación más conveniente.

### ¿Qué sigue?

Aprenda a proporcionar una entrada variable a través de un parámetro de línea de comandos y utilizar valores predeterminados de manera efectiva.

---

## 3. Use una entrada variable pasada en la línea de comandos

En su estado actual, nuestro workflow usa un saludo codificado en el comando del proceso.
Queremos agregar algo de flexibilidad usando una variable de entrada, para que podamos cambiar más fácilmente el saludo en tiempo de ejecución.

Esto requiere que hagamos tres conjuntos de cambios en nuestro script:

1. Cambiar el proceso para esperar una entrada variable
2. Configurar un parámetro de línea de comandos para capturar la entrada del usuario
3. Pasar la entrada al proceso en el cuerpo del workflow

Hagamos estos cambios uno a la vez.

### 3.1. Cambie el proceso `sayHello` para esperar una entrada variable

Necesitamos editar la definición del proceso para (1) aceptar una variable de entrada y (2) usar esa variable en la línea de comandos.

#### 3.1.1. Agregue un bloque input a la definición del proceso

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

La variable `greeting` tiene el prefijo `val` para decirle a Nextflow que es un valor (no una ruta).

#### 3.1.2. Edite el comando del proceso para usar la variable de entrada

Ahora intercambiamos el valor codificado original por el valor de la variable de entrada que esperamos recibir.

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

    Las llaves (`{ }`) eran técnicamente opcionales en versiones anteriores de Nextflow, por lo que puede ver workflows más antiguos donde esto está escrito como `echo '$greeting' > output.txt`.

Ahora que el proceso `sayHello()` está listo para aceptar una entrada variable, necesitamos una forma de proporcionar un valor de entrada a la llamada del proceso a nivel de workflow.

### 3.2. Configure un parámetro de línea de comandos para capturar la entrada del usuario

Podríamos simplemente codificar una entrada directamente haciendo la llamada al proceso `sayHello('Hello World!')`.
Sin embargo, cuando estemos haciendo trabajo real con nuestro workflow, vamos a querer poder controlar sus entradas desde la línea de comandos, para que podamos hacer algo como esto:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Afortunadamente, Nextflow tiene un sistema de parámetros de workflow incorporado llamado [`params`](https://nextflow.io/docs/latest/config.html#params) que facilita declarar y usar parámetros CLI.

La sintaxis general es declarar `params.<nombre_parámetro>` para decirle a Nextflow que espere un parámetro `--<nombre_parámetro>` en la línea de comandos.

Aquí, queremos crear un parámetro llamado `--input`, por lo que necesitamos declarar `params.input` en algún lugar del workflow.
En principio podemos escribirlo en cualquier lugar; pero como vamos a querer dárselo a la llamada del proceso `sayHello()`, podemos conectarlo allí directamente escribiendo `sayHello(params.input)`.

En el bloque del workflow, haga el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emite un saludo
    sayHello(params.input)
    ```

=== "Antes"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // emite un saludo
    sayHello()
    ```

Esto le dice a Nextflow que ejecute el proceso `sayHello` con el valor proporcionado a través del parámetro `--input`.

En efecto, hemos logrado los pasos (2) y (3) descritos al inicio de la sección de una sola vez.

### 3.3. Ejecute el comando del workflow

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

¡Et voilà!

Observe cómo la nueva ejecución ha sobrescrito el archivo de salida publicado en el directorio `results`.
Sin embargo, los resultados de las ejecuciones anteriores todavía se conservan en los directorios de tareas bajo `work`.

!!! tip "Consejo"

    Puede distinguir fácilmente los parámetros a nivel de Nextflow de los parámetros a nivel de pipeline.

    - Los parámetros que se aplican a un pipeline siempre llevan un guion doble (`--`).
    - Los parámetros que modifican una configuración de Nextflow, _p. ej._ la función `-resume` que usamos anteriormente, llevan un solo guion (`-`).

### 3.4. Use valores predeterminados para parámetros de línea de comandos

Ok, eso fue conveniente, pero en muchos casos, tiene sentido proporcionar un valor predeterminado para un parámetro dado para que no tenga que especificarlo en cada ejecución.

#### 3.4.1. Establezca un valor predeterminado para el parámetro CLI

Démosle al parámetro `input` un valor predeterminado declarándolo antes de la definición del workflow.

```groovy title="hello-world.nf" linenums="20"
/*
 * Parámetros del pipeline
 */
params {
    input: String = 'Holà mundo!'
}
```

Como puede ver, podemos especificar el tipo de entrada que el workflow espera (Nextflow 25.10.2 y posteriores).
La sintaxis es `nombre: Tipo = valor_predeterminado`.
Los tipos soportados incluyen `String`, `Integer`, `Float`, `Boolean` y `Path`.

!!! info "Info"

    En workflows más antiguos, puede ver que todo ese bloque `params` está escrito como solo `input = 'Holà mundo!'`.

A medida que agregue más parámetros a su pipeline, debe agregarlos todos a este bloque, ya sea que necesite darles un valor predeterminado o no.
Esto facilitará encontrar todos los parámetros configurables de un vistazo.

#### 3.4.2. Ejecute el workflow nuevamente sin especificar el parámetro

Ahora que tiene un valor predeterminado establecido, puede ejecutar el workflow nuevamente sin tener que especificar un valor en la línea de comandos.

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

La salida estará en el mismo lugar que anteriormente, pero el contenido debería actualizarse con el nuevo texto.

??? abstract "Contenido del archivo"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow usó el valor predeterminado del parámetro greeting para crear la salida.

#### 3.4.3. Sobrescriba el valor predeterminado

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
    Si el mismo parámetro se establece en valores diferentes en múltiples lugares, Nextflow determinará qué valor usar según el orden de precedencia que se describe [aquí](https://www.nextflow.io/docs/latest/config.html).

    Cubriremos esto con más detalle en la Parte 6 (Configuración).

### Conclusión

Sabe cómo usar una entrada variable simple proporcionada en tiempo de ejecución a través de un parámetro de línea de comandos, así como configurar, usar y sobrescribir valores predeterminados.

### ¿Qué sigue?

Aprenda cómo administrar ejecuciones de manera más conveniente.

---

## 4. Administre ejecuciones de workflow

Saber cómo lanzar workflows y recuperar salidas es excelente, pero rápidamente encontrará que hay algunos otros aspectos de la administración de workflows que harán su vida más fácil, especialmente si está desarrollando sus propios workflows.

Aquí le mostramos cómo usar la función [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) para cuando necesite volver a lanzar el mismo workflow, cómo inspeccionar el log de ejecuciones pasadas con [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log), y cómo eliminar directorios work más antiguos con [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

### 4.1. Vuelva a lanzar un workflow con `-resume`

A veces, va a querer volver a ejecutar un pipeline que ya ha lanzado anteriormente sin rehacer ningún paso que ya se completó exitosamente.

Nextflow tiene una opción llamada [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) que le permite hacer esto.
Específicamente, en este modo, cualquier proceso que ya se haya ejecutado con exactamente el mismo código, configuraciones y entradas se omitirá.
Esto significa que Nextflow solo ejecutará procesos que haya agregado o modificado desde la última ejecución, o a los que esté proporcionando nuevas configuraciones o entradas.

Hay dos ventajas clave al hacer esto:

- Si está en medio del desarrollo de su pipeline, puede iterar más rápidamente ya que solo tiene que ejecutar el(los) proceso(s) en el(los) que está trabajando activamente para probar sus cambios.
- Si está ejecutando un pipeline en producción y algo sale mal, en muchos casos puede solucionar el problema y volver a lanzar el pipeline, y se reanudará desde el punto de falla, lo que puede ahorrarle mucho tiempo y cómputo.

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

La salida de la consola debería verse familiar, pero hay una cosa que es un poco diferente en comparación con antes.

Busque la parte `cached:` que se ha agregado en la línea de estado del proceso (línea 5), lo que significa que Nextflow ha reconocido que ya ha hecho este trabajo y simplemente reutilizó el resultado de la ejecución exitosa anterior.

También puede ver que el hash del subdirectorio work es el mismo que en la ejecución anterior.
Nextflow literalmente le está señalando la ejecución anterior y diciendo "Ya hice eso allá".

!!! tip "Consejo"

    Cuando vuelve a ejecutar un pipeline con `resume`, Nextflow no sobrescribe ningún archivo publicado fuera del directorio work por ninguna ejecución que se ejecutó exitosamente anteriormente.

### 4.2. Inspeccione el log de ejecuciones pasadas

Ya sea que esté desarrollando un nuevo pipeline o ejecutando pipelines en producción, en algún momento probablemente necesitará buscar información sobre ejecuciones pasadas.
Aquí le mostramos cómo hacerlo.

Cada vez que lanza un workflow de nextflow, se escribe una línea en un archivo de log llamado `history`, bajo un directorio oculto llamado `.nextflow` en el directorio de trabajo actual.

??? abstract "Contenido del archivo"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Este archivo le proporciona la marca de tiempo, nombre de ejecución, estado, ID de revisión, ID de sesión y línea de comandos completa para cada ejecución de Nextflow que se ha lanzado desde dentro del directorio de trabajo actual.

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

Esto mostrará el contenido del archivo de log en la terminal, aumentado con una línea de encabezado.

Notará que el ID de sesión cambia cada vez que ejecuta un nuevo comando `nextflow run`, EXCEPTO si está usando la opción `-resume`.
En ese caso, el ID de sesión permanece igual.

Nextflow usa el ID de sesión para agrupar información de caché de ejecución bajo el directorio `cache`, también ubicado bajo `.nextflow`.

### 4.3. Elimine directorios work más antiguos

Durante el proceso de desarrollo, típicamente ejecutará su borrador de pipeline un gran número de veces, lo que puede llevar a una acumulación de muchos archivos en muchos subdirectorios.

Afortunadamente, Nextflow incluye un útil subcomando `clean` que puede eliminar automáticamente los subdirectorios work de ejecuciones pasadas que ya no le interesan.

#### 4.3.1. Determine los criterios de eliminación

Hay múltiples [opciones](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para determinar qué eliminar.

Aquí le mostramos un ejemplo que elimina todos los subdirectorios de ejecuciones anteriores a una ejecución dada, especificada usando su nombre de ejecución.

Busque la ejecución exitosa más reciente donde no usó `-resume`; en nuestro caso el nombre de ejecución fue `golden_cantor`.

El nombre de ejecución es la cadena de dos partes generada por la máquina que se muestra entre corchetes en la línea de salida de consola `Launching (...)`.
También puede usar el log de Nextflow para buscar una ejecución según su marca de tiempo y/o línea de comandos.

#### 4.3.2. Haga una ejecución de prueba

Primero usamos el flag de ejecución de prueba `-n` para verificar qué se eliminará dado el comando:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Salida del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Su salida tendrá diferentes nombres de directorio de tareas y puede tener un número diferente de líneas, pero debería verse similar al ejemplo.

Si no ve ninguna línea de salida, o no proporcionó un nombre de ejecución válido o no hay ejecuciones pasadas para eliminar. Asegúrese de cambiar `golden_cantor` en el comando de ejemplo por el nombre de ejecución más reciente correspondiente en su log.

#### 4.3.3. Proceda con la eliminación

Si la salida se ve como se esperaba y desea proceder con la eliminación, vuelva a ejecutar el comando con el flag `-f` en lugar de `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Salida del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

La salida debería ser similar a la anterior, pero ahora diciendo 'Removed' en lugar de 'Would remove'.
Tenga en cuenta que esto no elimina los subdirectorios de dos caracteres (como `a3/` arriba) pero sí vacía su contenido.

!!! Warning "Advertencia"

    Eliminar subdirectorios work de ejecuciones pasadas los elimina del caché de Nextflow y elimina cualquier salida que se almacenó en esos directorios.
    Eso significa que rompe la capacidad de Nextflow de reanudar la ejecución sin volver a ejecutar los procesos correspondientes.

    ¡Usted es responsable de guardar cualquier salida que le importe o en la que planee confiar! Esa es la razón principal por la que preferimos usar el modo `copy` en lugar del modo `symlink` para la directiva `publish`.

### Conclusión

Sabe cómo publicar salidas en un directorio específico, volver a lanzar un pipeline sin repetir pasos que ya se ejecutaron de manera idéntica, y usar el comando `nextflow clean` para limpiar directorios work antiguos.

Más generalmente, sabe cómo interpretar un workflow simple de Nextflow, administrar su ejecución y recuperar salidas.

### ¿Qué sigue?

Tómese un pequeño descanso, ¡se lo ha ganado!

Cuando esté listo, pase a [**Parte 2: Hello Channels**](./02_hello_channels.md) para aprender cómo usar canales para alimentar entradas en su workflow, lo que le permitirá aprovechar el paralelismo de flujo de datos incorporado de Nextflow y otras características poderosas.

---

## Quiz

<quiz>
¿Cuáles son los componentes mínimos requeridos de un proceso de Nextflow?
- [ ] Solo bloques input y output
- [x] Bloques output y script
- [ ] Bloques input, output y script
- [ ] Solo un bloque script

Aprenda más: [1.1.1. La definición del proceso](#111-the-process-definition)
</quiz>

<quiz>
¿Cuál es el propósito del bloque output en un proceso?
- [ ] Imprimir resultados en la consola
- [ ] Guardar archivos en el directorio work
- [x] Declarar las salidas esperadas del proceso
- [ ] Definir variables de entorno

Aprenda más: [1.1.1. La definición del proceso](#111-the-process-definition)
</quiz>

<quiz>
¿Qué comando se usa para ejecutar un workflow de Nextflow?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Mirando el directorio work de una tarea, ¿qué archivo contiene el comando real que se ejecutó?

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

Aprenda más: [1.2.2. Encuentre la salida y los logs en el directorio `work`](#122-find-the-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
¿Qué hace el flag `-resume`?
- [ ] Reinicia el workflow desde el principio
- [ ] Pausa el workflow
- [x] Omite procesos que ya se completaron exitosamente
- [ ] Crea una copia de seguridad del workflow

Aprenda más: [4.1. Vuelva a lanzar un workflow con `-resume`](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
¿Cuál es el modo predeterminado para publicar salidas de workflow?
- [ ] Copiar archivos al directorio de salida
- [x] Crear enlaces simbólicos en el directorio de salida
- [ ] Mover archivos al directorio de salida
- [ ] Comprimir archivos en el directorio de salida

Aprenda más: [2.3. Establezca el modo de publicación en copy](#23-set-the-publish-mode-to-copy)
</quiz>

<quiz>
¿Cómo se pasa un valor de parámetro a un workflow de Nextflow desde la línea de comandos?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Aprenda más: [3.2. Configure un parámetro de línea de comandos para capturar la entrada del usuario](#32-set-up-a-command-line-parameter-to-capture-user-input)
</quiz>

<quiz>
¿Cómo se hace referencia a una variable dentro de un bloque script de Nextflow?
- [ ] Usar sintaxis `%variable%`
- [x] Usar sintaxis `#!groovy ${variable}`
- [ ] Usar sintaxis `{{variable}}`
- [ ] Usar sintaxis `[variable]`
</quiz>
