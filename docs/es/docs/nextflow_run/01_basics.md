# Parte 1: Ejecutar operaciones básicas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta primera parte del curso de capacitación Nextflow Run, nos introducimos en el tema con un ejemplo muy básico de Hello World independiente del dominio, que usaremos para demostrar operaciones esenciales y señalar los componentes de código de Nextflow correspondientes.

??? info "¿Qué es un ejemplo Hello World?"

    Un "Hello World!" es un ejemplo minimalista que pretende demostrar la sintaxis y estructura básica de un lenguaje de programación o framework de software.
    El ejemplo típicamente consiste en imprimir la frase "Hello, World!" al dispositivo de salida, como la consola o terminal, o escribirla en un archivo.

---

## 1. Ejecutar un Hello World directamente

Demostremos este concepto con un comando simple que ejecutamos directamente en la terminal, para mostrar lo que hace antes de envolverlo en Nextflow.

!!! tip "Consejo"

    Recuerde que ahora debería estar dentro del directorio `nextflow-run/` como se describe en la página de [Primeros pasos](00_orientation.md).

### 1.1. Hacer que la terminal diga hola

Ejecute el siguiente comando en su terminal.

```bash
echo 'Hello World!'
```

??? success "Salida del comando"

    ```console
    Hello World!
    ```

Esto muestra el texto 'Hello World' directamente en la terminal.

### 1.2. Escribir la salida a un archivo

Ejecutar pipelines principalmente implica leer datos de archivos y escribir resultados en otros archivos, así que modifiquemos el comando para escribir la salida de texto a un archivo para hacer el ejemplo un poco más relevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Salida del comando"

    ```console

    ```

Esto no muestra nada en la terminal.

### 1.3. Encontrar la salida

El texto 'Hello World' ahora debería estar en el archivo de salida que especificamos, llamado `output.txt`.
Puede abrirlo en el explorador de archivos o desde la línea de comandos usando la utilidad `cat`, por ejemplo.

??? abstract "Contenido del archivo"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Esto es lo que vamos a intentar replicar con nuestro primer workflow de Nextflow.

### Conclusión

Ahora sabe cómo ejecutar un comando simple en la terminal que produce algún texto, y opcionalmente, cómo hacer que escriba la salida a un archivo.

### ¿Qué sigue?

Descubra qué se necesita para ejecutar un workflow de Nextflow que logre el mismo resultado.

---

## 2. Ejecutar el workflow

Le proporcionamos un script de workflow llamado `1-hello.nf` que toma un saludo de entrada a través de un argumento de línea de comandos llamado `--input` y produce un archivo de texto que contiene ese saludo.

No vamos a mirar el código todavía; primero veamos cómo se ve ejecutarlo.

### 2.1. Lanzar el workflow y monitorear la ejecución

En la terminal, ejecute el siguiente comando:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Salida del comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Si su salida de consola se ve algo así, ¡felicidades, acaba de ejecutar su primer workflow de Nextflow!

La salida más importante aquí es la última línea, que está resaltada en la salida anterior:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Esto nos dice que el process `sayHello` fue ejecutado exitosamente una vez (`1 of 1 ✔`).

Genial, pero puede estar preguntándose: ¿dónde está la salida?

### 2.2. Encontrar el archivo de salida en el directorio `results`

Este workflow está configurado para publicar su salida en un directorio de resultados.
Si mira su directorio actual, verá que cuando ejecutó el workflow, Nextflow creó un nuevo directorio llamado `results`, así como un subdirectorio llamado `1-hello` bajo él, que contiene un archivo llamado `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Abra el archivo; el contenido debería coincidir con la cadena que especificó en la línea de comandos.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

¡Genial, nuestro workflow hizo lo que se suponía que debía hacer!

### 2.3. Guardar los resultados en un directorio diferente

Por defecto, Nextflow guardará las salidas del pipeline en un directorio llamado `results` en su ruta actual.
Para cambiar dónde se publican sus archivos, use la bandera CLI `-output-dir` (o `-o` para abreviar)

!!! danger "Advertencia"

    ¡Note que `--input` tiene dos guiones y `-output-dir` tiene uno!
    Esto es porque `--input` es un _parámetro_ del pipeline y `-output-dir` es una bandera CLI central de Nextflow.
    Más sobre esto más adelante.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Debería ver que sus salidas ahora se publican en un directorio llamado `hello_results` en lugar de `results`:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Los archivos dentro de este directorio son exactamente los mismos que antes, solo que el directorio de nivel superior es diferente.
Sin embargo, tenga en cuenta que en ambos casos el resultado 'publicado' es una copia (o en algunos casos un enlace simbólico) de la salida real producida por Nextflow cuando ejecutó el workflow.

Ahora vamos a echar un vistazo bajo el capó para ver dónde ejecutó Nextflow realmente el trabajo.

!!! warning "Advertencia"

    No todos los workflows estarán configurados para publicar salidas en un directorio de resultados, y/o los nombres y estructura de directorios pueden ser diferentes.
    Un poco más adelante en esta sección, le mostraremos cómo averiguar dónde se especifica este comportamiento.

### 2.4. Encontrar la salida original y los registros en el directorio `work/`

Cuando ejecuta un workflow, Nextflow crea un 'directorio de tarea' distinto para cada invocación individual de cada process en el workflow (=cada paso en el pipeline).
Para cada uno, preparará las entradas necesarias, ejecutará la(s) instrucción(es) relevante(s) y escribirá las salidas y archivos de registro dentro de ese único directorio, que se nombra automáticamente usando un hash para hacerlo único.

Todos estos directorios de tareas vivirán bajo un directorio llamado `work` dentro de su directorio actual (donde está ejecutando el comando).

Eso puede sonar confuso, así que veamos cómo se ve en la práctica.

Volviendo a la salida de consola del workflow que ejecutamos antes, teníamos esta línea:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

¿Ve cómo la línea comienza con `[a3/1e1535]`?
Esa es una forma truncada de la ruta del directorio de tarea para esa llamada de process, y le dice dónde encontrar la salida de la llamada del process `sayHello` dentro de la ruta del directorio `work/`.

Puede encontrar la ruta completa escribiendo el siguiente comando (reemplazando `a3/1e1535` con lo que ve en su propia terminal) y presionando la tecla tab para autocompletar la ruta o agregando un asterisco:

```bash
ls work/a3/1e1535*
```

Esto debería producir la ruta completa del directorio: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Veamos qué hay ahí.

??? abstract "Contenidos del directorio"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
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
    Sin embargo, los archivos de registro están configurados para ser invisibles en la terminal, así que si quiere usar `ls` o `tree` para verlos, necesitará establecer la opción relevante para mostrar archivos invisibles.

    ```bash
    tree -a work
    ```

Hay dos conjuntos de directorios en `work/`, de las dos ejecuciones diferentes del pipeline que hemos hecho.
Cada ejecución de tarea obtiene su propio directorio aislado para trabajar.
En este caso el pipeline hizo lo mismo ambas veces, así que los contenidos de cada directorio de tarea son idénticos

Debería reconocer inmediatamente el archivo `output.txt`, que es de hecho la salida original del process `sayHello` que se publicó en el directorio `results`.
Si lo abre, encontrará el saludo `Hello World!` nuevamente.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

¿Entonces qué pasa con todos esos otros archivos?

Estos son los archivos auxiliares y de registro que Nextflow escribió como parte de la ejecución de la tarea:

- **`.command.begin`**: Archivo centinela creado tan pronto como se lanza la tarea.
- **`.command.err`**: Mensajes de error (`stderr`) emitidos por la llamada del process
- **`.command.log`**: Salida de registro completa emitida por la llamada del process
- **`.command.out`**: Salida regular (`stdout`) por la llamada del process
- **`.command.run`**: Script completo ejecutado por Nextflow para ejecutar la llamada del process
- **`.command.sh`**: El comando que realmente ejecutó la llamada del process
- **`.exitcode`**: El código de salida resultante del comando

El archivo `.command.sh` es especialmente útil porque le muestra el comando principal que Nextflow ejecutó, sin incluir toda la contabilidad y configuración de tarea/entorno.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Esto confirma que el workflow compuso el mismo comando que ejecutamos directamente en la línea de comandos anteriormente.

Cuando algo sale mal y necesita solucionar lo que sucedió, puede ser útil mirar el script `command.sh` para verificar exactamente qué comando compuso Nextflow basándose en las instrucciones del workflow, interpolación de variables y demás.

### 2.5. Re-ejecutar el workflow con diferentes saludos

Intente re-ejecutar el workflow algunas veces con diferentes valores para el argumento `--input`, luego mire los directorios de tarea.

??? abstract "Contenidos del directorio"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Puede ver que se ha creado un nuevo subdirectorio con un conjunto completo de archivos de salida y registro para cada ejecución.

En contraste, si mira el directorio `results`, todavía hay solo un conjunto de resultados, y el contenido del archivo de salida corresponde a lo que ejecutó por última vez.

??? abstract "Contenidos del directorio"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Esto le muestra que los resultados publicados serán sobrescritos por ejecuciones posteriores, mientras que los directorios de tarea bajo `work/` se preservan.

### Conclusión

Sabe cómo ejecutar un script simple de Nextflow, monitorear su ejecución y encontrar sus salidas.

### ¿Qué sigue?

Aprenda cómo leer un script básico de Nextflow e identificar cómo sus componentes se relacionan con su funcionalidad.

---

## 3. Examinar el script inicial del workflow Hello World

Lo que hicimos allí fue básicamente tratar el script del workflow como una caja negra.
Ahora que hemos visto lo que hace, abramos la caja y miremos adentro.

Nuestro objetivo aquí no es memorizar la sintaxis del código de Nextflow, sino formar alguna intuición básica de cuáles son los componentes principales y cómo están organizados.

### 3.1. Examinar la estructura general del código

Encontrará el script `1-hello.nf` en su directorio actual, que debería ser `nextflow-run`. Ábralo en el panel del editor.

??? full-code "Archivo de código completo"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usar echo para imprimir 'Hello World!' a un archivo
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // emitir un saludo
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

Un script de workflow de Nextflow típicamente incluye una o más definiciones de **process**, el **workflow** en sí, y algunos bloques opcionales como **params** y **output**.

Cada **process** describe qué operación(es) debería realizar el paso correspondiente en el pipeline, mientras que el **workflow** describe la lógica de flujo de datos que conecta los varios pasos.

Echemos un vistazo más de cerca al bloque **process** primero, luego veremos el bloque **workflow**.

### 3.2. La definición del `process`

El primer bloque de código describe un [**process**](https://nextflow.io/docs/latest/process.html).
La definición del process comienza con la palabra clave `process`, seguida del nombre del process y finalmente el cuerpo del process delimitado por llaves.
El cuerpo del process debe contener un bloque script que especifica el comando a ejecutar, que puede ser cualquier cosa que pueda ejecutar en una terminal de línea de comandos.

```groovy title="1-hello.nf" linenums="3"
/*
* Usar echo para imprimir un saludo a un archivo
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

Aquí tenemos un **process** llamado `sayHello` que toma una variable de **entrada** llamada `greeting` y escribe su **salida** a un archivo llamado `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Esta es una definición de process muy mínima que solo contiene una definición de `input`, una definición de `output` y el `script` a ejecutar.

La definición de `input` incluye el calificador `val`, que le dice a Nextflow que espere un valor de algún tipo (puede ser una cadena, un número, lo que sea).

La definición de `output` incluye el calificador `path`, que le dice a Nextflow que esto debe manejarse como una ruta (incluye tanto rutas de directorio como archivos).

### 3.3. La definición del `workflow`

El segundo bloque de código describe el [**workflow**](https://nextflow.io/docs/latest/workflow.html) en sí.
La definición del workflow comienza con la palabra clave `workflow`, seguida de un nombre opcional, luego el cuerpo del workflow delimitado por llaves.

Aquí tenemos un **workflow** que consiste en un bloque `main:` y un bloque `publish:`.
El bloque `main:` es el cuerpo principal del workflow y el bloque `publish:` lista las salidas que deben publicarse en el directorio `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emitir un saludo
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

En este caso el bloque `main:` contiene una llamada al process `sayHello` y le da una entrada llamada `params.input` para usar como el saludo.

Como discutiremos con más detalle en un momento, `params.input` contiene el valor que dimos al parámetro `--input` en nuestra línea de comandos.

El bloque `publish:` lista la salida de la llamada del process `sayHello()`, a la cual se refiere como `sayHello.out` y le da el nombre `first_output` (esto puede ser cualquier cosa que el autor del workflow quiera).

Esta es una definición de **workflow** muy mínima.
En un pipeline del mundo real, el workflow típicamente contiene múltiples llamadas a **processes** conectados por **channels**, y puede haber valores predeterminados configurados para las entradas variables.

Llegaremos a eso en la Parte 2 del curso.
Por ahora, echemos un vistazo más de cerca a cómo nuestro workflow está manejando entradas y salidas.

### 3.4. El sistema `params` de parámetros de línea de comandos

El `params.input` que proporcionamos a la llamada del process `sayHello()` es un código elegante de Nextflow y vale la pena dedicarle un minuto extra.

Como se mencionó anteriormente, así es como pasamos el valor del parámetro de línea de comandos `--input` a la llamada del process `sayHello()`.
De hecho, simplemente declarar `params.someParameterName` es suficiente para dar al workflow un parámetro llamado `--someParameterName` desde la línea de comandos.

Aquí hemos formalizado esa declaración de parámetro configurando un bloque `params` que especifica el tipo de entrada que espera el workflow (Nextflow 25.10.2 y posterior).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

Los tipos soportados incluyen `String`, `Integer`, `Float`, `Boolean` y `Path`.
Para más información, consulte [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) en la documentación de referencia de Nextflow.

!!! tip "Consejo"

    Recuerde que los parámetros del _workflow_ declarados usando el sistema `params` siempre llevan dos guiones en la línea de comandos (`--`).
    Esto los distingue de las banderas CLI _a nivel de Nextflow_, que solo llevan un guión (`-`).

### 3.5. La directiva `publish`

En el otro extremo del workflow, ya hemos echado un vistazo al bloque `publish:`.
Esa es una mitad del sistema de manejo de salidas; la otra mitad es el bloque `output` ubicado abajo.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Esto especifica que la salida `first_output` listada en el bloque `publish:` debe copiarse a un subdirectorio llamado `1-hello` bajo el directorio de salida `results` predeterminado.

La línea `mode 'copy'` anula el comportamiento predeterminado del sistema, que es hacer un enlace simbólico (o symlink) al archivo original en el directorio `work/` en lugar de una copia propiamente dicha.

Hay más opciones que las mostradas aquí para controlar el comportamiento de publicación; cubriremos algunas más adelante.
También verá que cuando un workflow genera múltiples salidas, cada una se lista de esta manera en el bloque `output`.

Para más información, consulte [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) en la documentación de referencia de Nextflow.

??? info "Sintaxis antigua para publicar salidas usando `publishDir`"

    Hasta muy recientemente, la forma establecida de publicar salidas era hacerlo a nivel de cada process individual usando una directiva `publishDir`.

    Todavía encontrará este patrón de código en todas partes en pipelines de Nextflow más antiguos y módulos de process, por lo que es importante estar al tanto de ello.

    En lugar de tener un bloque `publish:` en el workflow y un bloque `output` en el nivel superior, vería una línea `publishDir` en la definición del process `sayHello`:

    ```groovy title="Ejemplo de sintaxis" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Sin embargo, no recomendamos usar esto en ningún trabajo nuevo ya que eventualmente será prohibido en futuras versiones del lenguaje Nextflow.

### Conclusión

Ahora sabe cómo está estructurado un workflow simple de Nextflow, y cómo los componentes básicos se relacionan con su funcionalidad.

### ¿Qué sigue?

Aprenda a gestionar las ejecuciones de su workflow de manera conveniente.

---

## 4. Gestionar ejecuciones de workflow

Saber cómo lanzar workflows y recuperar salidas es genial, pero rápidamente encontrará que hay algunos otros aspectos de la gestión de workflows que harán su vida más fácil.

Aquí le mostramos cómo aprovechar la función `resume` para cuando necesite re-lanzar el mismo workflow, cómo inspeccionar los registros de ejecución con `nextflow log`, y cómo eliminar directorios de trabajo antiguos con `nextflow clean`.

### 4.1. Re-lanzar un workflow con `-resume`

A veces, va a querer re-ejecutar un pipeline que ya lanzó anteriormente sin rehacer ningún trabajo que ya se completó exitosamente.

Nextflow tiene una opción llamada `-resume` que le permite hacer esto.
Específicamente, en este modo, cualquier process que ya se haya ejecutado con exactamente el mismo código, configuración y entradas será omitido.
Esto significa que Nextflow solo ejecutará los processes que haya agregado o modificado desde la última ejecución, o a los que esté proporcionando nuevas configuraciones o entradas.

Hay dos ventajas clave de hacer esto:

- Si está en medio del desarrollo de un pipeline, puede iterar más rápidamente ya que solo tiene que ejecutar el o los process(es) en los que está trabajando activamente para probar sus cambios.
- Si está ejecutando un pipeline en producción y algo sale mal, en muchos casos puede arreglar el problema y relanzar el pipeline, y reanudará la ejecución desde el punto de falla, lo que puede ahorrarle mucho tiempo y cómputo.

Para usarlo, simplemente agregue `-resume` a su comando y ejecútelo:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Salida del comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

La salida de consola debería verse familiar, pero hay una cosa que es un poco diferente comparado con antes.

Busque el bit `cached:` que se ha agregado en la línea de estado del process (línea 5), lo que significa que Nextflow ha reconocido que ya hizo este trabajo y simplemente reutilizó el resultado de la ejecución exitosa anterior.

También puede ver que el hash del subdirectorio de trabajo es el mismo que en la ejecución anterior.
Nextflow está literalmente señalándole la ejecución anterior y diciendo "Ya hice eso ahí".

!!! tip "Consejo"

    Cuando re-ejecuta un pipeline con `resume`, Nextflow no sobrescribe ningún archivo publicado fuera del directorio de trabajo por cualquier ejecución que se ejecutó exitosamente anteriormente.

    Para más información, consulte [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) en la documentación de referencia de Nextflow.

### 4.2. Inspeccionar el registro de ejecuciones pasadas

Cada vez que lanza un workflow de Nextflow, se escribe una línea en un archivo de registro llamado `history`, bajo un directorio oculto llamado `.nextflow` en el directorio de trabajo actual.

??? abstract "Contenido del archivo"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Este archivo le da la marca de tiempo, nombre de ejecución, estado, ID de revisión, ID de sesión y línea de comandos completa para cada ejecución de Nextflow que se ha lanzado desde el directorio de trabajo actual.

Una forma más conveniente de acceder a esta información es usar el comando [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

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

Esto mostrará el contenido del archivo de registro en la terminal, aumentado con una línea de encabezado.

Notará que el ID de sesión cambia cada vez que ejecuta un nuevo comando `nextflow run`, EXCEPTO si está usando la opción `-resume`.
En ese caso, el ID de sesión permanece igual.

Nextflow usa el ID de sesión para agrupar información de caché de ejecución bajo el directorio `cache`, también ubicado bajo `.nextflow`.

### 4.3. Eliminar directorios de trabajo antiguos

Si ejecuta muchos pipelines, puede terminar acumulando muchos archivos a través de muchos subdirectorios.
Dado que los subdirectorios se nombran aleatoriamente, es difícil saber por sus nombres cuáles son ejecuciones más antiguas vs. más recientes.

Afortunadamente Nextflow incluye un comando útil llamado [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) que puede eliminar automáticamente los subdirectorios de trabajo de ejecuciones pasadas que ya no le importan.

#### 4.3.1. Determinar criterios de eliminación

Hay múltiples opciones para determinar qué eliminar, que puede explorar en la documentación vinculada arriba.
Aquí le mostramos un ejemplo que elimina todos los subdirectorios de ejecuciones anteriores a una ejecución dada, especificada usando su nombre de ejecución.

Busque la ejecución exitosa más reciente donde no usó `-resume`; en nuestro caso el nombre de ejecución fue `backstabbing_swartz`.

El nombre de ejecución es la cadena de dos partes generada por la máquina mostrada entre corchetes en la línea de salida de consola `Launching (...)`.
También puede usar el registro de Nextflow para buscar una ejecución basándose en su marca de tiempo y/o línea de comandos.

#### 4.3.2. Hacer una ejecución de prueba

Primero usamos la bandera de ejecución de prueba `-n` para verificar qué se eliminará dado el comando:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Salida del comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Su salida tendrá diferentes nombres de directorio de tarea y puede tener un número diferente de líneas, pero debería verse similar al ejemplo.

Si no ve ninguna línea de salida, ya sea que no proporcionó un nombre de ejecución válido o no hay ejecuciones pasadas para eliminar. Asegúrese de cambiar `backstabbing_swartz` en el comando de ejemplo a lo que sea el nombre de ejecución más reciente correspondiente en su registro.

#### 4.3.3. Proceder con la eliminación

Si la salida se ve como se esperaba y quiere proceder con la eliminación, re-ejecute el comando con la bandera `-f` en lugar de `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Salida del comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

La salida debería ser similar a antes, pero ahora diciendo 'Removed' en lugar de 'Would remove'.
Note que esto no elimina los subdirectorios de dos caracteres (como `eb/` arriba) pero sí vacía su contenido.

!!! warning "Advertencia"

    Eliminar subdirectorios de trabajo de ejecuciones pasadas los elimina de la caché de Nextflow y elimina cualquier salida que se almacenó en esos directorios.
    Eso significa que rompe la capacidad de Nextflow de reanudar la ejecución sin re-ejecutar los processes correspondientes.

    ¡Usted es responsable de guardar cualquier salida que le importe! Esa es la razón principal por la que preferimos usar el modo `copy` en lugar del modo `symlink` para la directiva `publish`.

### Conclusión

Sabe cómo relanzar un pipeline sin repetir pasos que ya se ejecutaron de manera idéntica, inspeccionar el registro de ejecución, y usar el comando `nextflow clean` para limpiar directorios de trabajo antiguos.

### ¿Qué sigue?

¡Tome un pequeño descanso! Acaba de absorber los bloques de construcción de la sintaxis de Nextflow e instrucciones básicas de uso.

En la próxima sección de esta capacitación, vamos a ver cuatro versiones sucesivamente más realistas del pipeline Hello World que demostrarán cómo Nextflow le permite procesar múltiples entradas eficientemente, ejecutar workflows compuestos de múltiples pasos conectados, aprovechar componentes de código modulares, y utilizar contenedores para mayor reproducibilidad y portabilidad.

---

## Cuestionario

<quiz>
En la línea de salida de consola `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, ¿qué representa `[a3/7be2fa]`?
- [ ] El número de versión del process
- [ ] Un identificador único de ejecución
- [x] La ruta truncada al directorio de trabajo de la tarea
- [ ] La suma de verificación del archivo de salida

Más información: [2.4. Encontrar la salida original y los registros en el directorio `work/`](#24-encontrar-la-salida-original-y-los-registros-en-el-directorio-work)
</quiz>

<quiz>
¿Cuál es el propósito del archivo `.command.sh` en un directorio de tarea?
- [ ] Almacena la configuración de la tarea
- [x] Muestra el comando real que fue ejecutado por el process
- [ ] Contiene mensajes de error de tareas fallidas
- [ ] Lista los archivos de entrada preparados para la tarea

Más información: [2.4. Encontrar la salida original y los registros en el directorio `work/`](#24-encontrar-la-salida-original-y-los-registros-en-el-directorio-work)
</quiz>

<quiz>
¿Qué sucede con los resultados publicados cuando re-ejecuta un workflow sin `-resume`?
- [ ] Se preservan en directorios separados con marca de tiempo
- [x] Son sobrescritos por la nueva ejecución
- [ ] Nextflow previene la sobrescritura y falla
- [ ] Son respaldados automáticamente

Más información: [2.5. Re-ejecutar el workflow con diferentes saludos](#25-re-ejecutar-el-workflow-con-diferentes-saludos)
</quiz>

<quiz>
¿Qué indica esta salida de consola?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] La tarea falló y fue omitida
- [ ] La tarea está esperando en una cola
- [x] Nextflow reutilizó resultados de una ejecución idéntica anterior
- [ ] La tarea fue cancelada manualmente

Más información: [4.1. Re-lanzar un workflow con `-resume`](#41-re-lanzar-un-workflow-con--resume)
</quiz>

<quiz>
¿Dónde almacena Nextflow el historial de ejecución que muestra el comando `nextflow log`?
- [ ] En el directorio de resultados
- [ ] En el directorio de trabajo
- [x] En el archivo `.nextflow/history`
- [ ] En `nextflow.config`

Más información: [4.2. Inspeccionar el registro de ejecuciones pasadas](#42-inspeccionar-el-registro-de-ejecuciones-pasadas)
</quiz>

<quiz>
¿Cuál es el propósito del bloque `params` en un archivo de workflow?
- [ ] Definir los requisitos de recursos del process
- [ ] Configurar el executor
- [x] Declarar y tipificar los parámetros de entrada del workflow
- [ ] Especificar opciones de publicación de salida

Más información: [3.4. El sistema params de parámetros de línea de comandos](#34-el-sistema-params-de-parámetros-de-línea-de-comandos)
</quiz>

<quiz>
En el bloque `output` del workflow, ¿qué hace `mode 'copy'`?
- [ ] Crea un respaldo del directorio de trabajo
- [x] Hace una copia completa de archivos en lugar de enlaces simbólicos
- [ ] Copia el script del workflow a los resultados
- [ ] Habilita la copia incremental de archivos

Más información: [3.5. La directiva publish](#35-la-directiva-publish)
</quiz>

<quiz>
¿Cuál es la bandera recomendada para usar con el comando `nextflow clean` antes de realmente eliminar archivos?
- [x] `-n` (ejecución de prueba) para previsualizar lo que se eliminaría
- [ ] `-v` (verbose) para ver salida detallada
- [ ] `-a` (all) para seleccionar todos los directorios
- [ ] `-q` (quiet) para suprimir advertencias

Más información: [4.3. Eliminar directorios de trabajo antiguos](#43-eliminar-directorios-de-trabajo-antiguos)
</quiz>
