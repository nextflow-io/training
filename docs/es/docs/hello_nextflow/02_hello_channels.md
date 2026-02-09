# Parte 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulte [la lista de reproducción completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) en el canal de YouTube de Nextflow.

:green_book: La transcripción del video está disponible [aquí](./transcripts/02_hello_channels.md).
///

En la Parte 1 de este curso (Hello World), le mostramos cómo proporcionar una entrada variable a un proceso proporcionando la entrada directamente en la llamada al proceso: `sayHello(params.input)`.
Ese fue un enfoque deliberadamente simplificado.
En la práctica, ese enfoque tiene limitaciones importantes; es decir, solo funciona para casos muy simples donde solo queremos ejecutar el proceso una vez, con un solo valor.
En la mayoría de los casos de uso de workflow realistas, queremos procesar múltiples valores (datos experimentales para múltiples muestras, por ejemplo), por lo que necesitamos una forma más sofisticada de manejar las entradas.

Para eso están los [**canales**](https://nextflow.io/docs/latest/channel.html) de Nextflow.
Los canales son colas diseñadas para manejar entradas de manera eficiente y transportarlas de un paso a otro en workflows de múltiples pasos, al tiempo que proporcionan paralelismo integrado y muchos beneficios adicionales.

En esta parte del curso, aprenderá a usar un canal para manejar múltiples entradas de una variedad de fuentes diferentes.
También aprenderá a usar [**operadores**](https://nextflow.io/docs/latest/reference/operator.html) para transformar el contenido del canal según sea necesario.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado la Parte 1 del curso [Hello Nextflow](./index.md), pero si se siente cómodo con los conceptos básicos cubiertos en esa sección, puede comenzar desde aquí sin hacer nada especial.

---

## 0. Calentamiento: Ejecute `hello-channels.nf`

Vamos a usar el script de workflow `hello-channels.nf` como punto de partida.
Es equivalente al script producido al trabajar en la Parte 1 de este curso de capacitación, excepto que hemos cambiado el destino de salida:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Solo para asegurarnos de que todo funciona, ejecute el script una vez antes de hacer cualquier cambio:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Como antes, encontrará el archivo de salida llamado `output.txt` en el directorio `results/hello_channels` (como se especifica en el bloque `output` del script de workflow, mostrado arriba).

??? abstract "Contenido del directorio"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contenido del archivo"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Si eso funcionó para usted, está listo para aprender sobre canales.

---

## 1. Proporcione entradas variables a través de un canal explícitamente

Vamos a crear un **canal** para pasar la entrada variable al proceso `sayHello()` en lugar de depender del manejo implícito, que tiene ciertas limitaciones.

### 1.1. Cree un canal de entrada

Hay una variedad de [**constructores de canales**](https://nextflow.io/docs/latest/reference/channel.html) que podemos usar para configurar un canal.
Para mantener las cosas simples por ahora, vamos a usar el constructor de canales más básico, llamado [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of), que creará un canal que contiene un solo valor.
Funcionalmente, esto será similar a cómo lo teníamos configurado antes, pero en lugar de que Nextflow cree un canal implícitamente, ahora lo estamos haciendo explícitamente.

Esta es la línea de código que vamos a usar:

```console title="Sintaxis"
greeting_ch = channel.of('Hello Channels!')
```

Esto crea un canal llamado `greeting_ch` usando el constructor de canales `channel.of()`, que configura un canal de cola simple, y carga la cadena `'Hello Channels!'` para usar como valor de saludo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Nota"

    Estamos volviendo temporalmente a cadenas codificadas en lugar de usar un parámetro CLI por razones de legibilidad. Volveremos a usar parámetros CLI una vez que hayamos cubierto lo que está sucediendo a nivel del canal.

En el bloque workflow, agregue el código del constructor de canales:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite un saludo
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // emite un saludo
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Esto aún no es funcional ya que todavía no hemos cambiado la entrada a la llamada del proceso.

### 1.2. Agregue el canal como entrada a la llamada del proceso

Ahora necesitamos conectar nuestro canal recién creado a la llamada del proceso `sayHello()`, reemplazando el parámetro CLI que estábamos proporcionando directamente antes.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite un saludo
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Esto le dice a Nextflow que ejecute el proceso `sayHello` con el contenido del canal `greeting_ch`.

Ahora nuestro workflow es completamente funcional; es el equivalente explícito de escribir `sayHello('Hello Channels!')`.

### 1.3. Ejecute el workflow

¡Ejecutémoslo!

```bash
nextflow run hello-channels.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Si realizó ambas ediciones correctamente, debería obtener una ejecución exitosa.
Puede verificar el directorio de resultados para asegurarse de que el resultado sigue siendo el mismo que antes.

??? abstract "Contenido del archivo"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Así que hemos aumentado la flexibilidad de nuestro workflow mientras logramos el mismo resultado final.
Esto puede parecer que estamos escribiendo más código sin ningún beneficio tangible, pero el valor quedará claro tan pronto como comencemos a manejar más entradas.

Como adelanto de eso, veamos una cosa más antes de continuar: un pequeño pero conveniente beneficio de usar un canal explícito para administrar la entrada de datos.

### 1.4. Use `view()` para inspeccionar el contenido del canal

Los canales de Nextflow están construidos de una manera que nos permite operar en su contenido usando operadores, que cubriremos en detalle más adelante en este capítulo.

Por ahora, solo vamos a mostrarle cómo usar un operador súper simple llamado [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) para inspeccionar el contenido de un canal.
Puede pensar en `view()` como una herramienta de depuración, como una declaración `print()` en Python, o su equivalente en otros lenguajes.

Agregue esta pequeña línea al bloque workflow:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

La cantidad exacta de espacios no importa siempre que sea un múltiplo de 4; solo estamos tratando de alinear el inicio de la declaración `.view()` con la parte `.of()` de la construcción del canal.

Ahora ejecute el workflow nuevamente:

```bash
nextflow run hello-channels.nf
```

??? success "Salida del comando"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Como puede ver, esto muestra el contenido del canal en la consola.
Aquí solo tenemos un elemento, pero cuando comencemos a cargar múltiples valores en el canal en la siguiente sección, verá que esto está configurado para mostrar un elemento por línea.

### Conclusión

Sabe cómo usar un constructor de canales básico para proporcionar una entrada a un proceso.

### ¿Qué sigue?

Aprenda a usar canales para hacer que el workflow itere sobre múltiples valores de entrada.

---

## 2. Modifique el workflow para ejecutarse con múltiples valores de entrada

Los workflows típicamente se ejecutan en lotes de entradas que están destinadas a procesarse en masa, por lo que queremos actualizar el workflow para aceptar múltiples valores de entrada.

### 2.1. Cargue múltiples saludos en el canal de entrada

Convenientemente, el constructor de canales `channel.of()` que hemos estado usando está muy contento de aceptar más de un valor, por lo que no necesitamos modificar eso en absoluto.
Simplemente podemos cargar múltiples valores en el canal.

Hagámoslos `'Hello'`, `'Bonjour'` y `'Holà'`.

#### 2.1.1. Agregue más saludos

Antes del bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // crea un canal para las entradas
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // crea un canal para las entradas
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

La documentación nos dice que esto debería funcionar. ¿Puede ser realmente tan simple?

#### 2.1.2. Ejecute el comando y observe la salida del registro

Intentémoslo.

```bash
nextflow run hello-channels.nf
```

??? success "Salida del comando"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Ciertamente parece haberse ejecutado bien.
El monitor de ejecución muestra que se realizaron `3 de 3` llamadas para el proceso `sayHello`, y vemos los tres saludos enumerados por la declaración `view()`, uno por línea como se prometió.

Sin embargo, todavía solo hay una salida en el directorio de resultados:

??? abstract "Contenido del directorio"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Contenido del archivo"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Debería ver uno de los tres saludos allí, aunque el que obtuvo podría ser diferente de lo que se muestra aquí.
¿Puede pensar por qué podría ser eso?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_En el diagrama, el canal está representado en verde, y el orden de los elementos está representado como canicas en una tubería: el primero cargado está a la derecha, luego el segundo en el medio, luego el tercero está a la izquierda._

Mirando hacia atrás en el monitor de ejecución, nos dio solo una ruta de subdirectorio (`f4/c9962c`).
Echemos un vistazo allí.

??? abstract "Contenido del directorio"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contenido del archivo"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

¡Ese ni siquiera es el mismo saludo que obtuvimos en el directorio de resultados! ¿Qué está pasando?

En este punto, necesitamos decirle que, por defecto, el sistema de registro ANSI escribe el registro de múltiples llamadas al mismo proceso en la misma línea.
Entonces, el estado de las tres llamadas al proceso sayHello() están aterrizando en el mismo lugar.

Afortunadamente, podemos deshabilitar ese comportamiento para ver la lista completa de llamadas a procesos.

#### 2.1.3. Ejecute el comando nuevamente con la opción `-ansi-log false`

Para expandir el registro y mostrar una línea por llamada de proceso, agregue `-ansi-log false` al comando.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Salida del comando"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Esta vez vemos las tres ejecuciones de procesos y sus subdirectorios de trabajo asociados listados en la salida.

Eso es mucho mejor, al menos para un workflow simple.
Para un workflow complejo, o un gran número de entradas, tener la lista completa en la terminal sería un poco abrumador.
Es por eso que `-ansi-log false` no es el comportamiento predeterminado.

!!! tip "Consejo"

    La forma en que se informa el estado es un poco diferente entre los dos modos de registro.
    En el modo condensado, Nextflow informa si las llamadas se completaron con éxito o no.
    En este modo expandido, solo informa que fueron enviadas.

De todos modos, ahora que tenemos los subdirectorios de cada llamada de proceso, podemos buscar sus registros y salidas.

??? abstract "Contenido del directorio"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Contenido del archivo"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Esto muestra que los tres procesos se ejecutaron con éxito (¡bien!).

Dicho esto, todavía tenemos el problema de que solo hay un archivo de salida en el directorio de resultados.

Puede recordar que codificamos el nombre del archivo de salida para el proceso `sayHello`, por lo que las tres llamadas produjeron un archivo llamado `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Mientras los archivos de salida permanezcan en los subdirectorios de trabajo, aislados de los otros procesos, eso está bien.
Pero cuando se publican en el mismo directorio de resultados, el que se copió primero es sobrescrito por el siguiente, y así sucesivamente.

### 2.2. Asegúrese de que los nombres de los archivos de salida sean únicos

Podemos continuar publicando todas las salidas en el mismo directorio de resultados, pero necesitamos asegurarnos de que tendrán nombres únicos.
Específicamente, necesitamos modificar el primer proceso para generar un nombre de archivo dinámicamente para que los nombres de archivo finales sean únicos.

Entonces, ¿cómo hacemos que los nombres de archivo sean únicos?
Una forma común de hacerlo es usar algún dato único de metadatos de las entradas (recibidas del canal de entrada) como parte del nombre del archivo de salida.
Aquí, por conveniencia, solo usaremos el saludo en sí, ya que es solo una cadena corta, y lo antepondremos al nombre base del archivo de salida.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Construya un nombre de archivo de salida dinámico

En el bloque de proceso, realice los siguientes cambios de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

Asegúrese de reemplazar `output.txt` tanto en la definición de salida como en el bloque de comando `script:`.

!!! tip "Consejo"

    En la definición de salida, DEBE usar comillas dobles alrededor de la expresión del nombre de archivo de salida (NO comillas simples), de lo contrario fallará.

Esto debería producir un nombre de archivo de salida único cada vez que se llame al proceso, para que pueda distinguirse de las salidas de otras llamadas al mismo proceso en el directorio de salida.

#### 2.2.2. Ejecute el workflow

Ejecutémoslo. Tenga en cuenta que estamos de vuelta a ejecutar con la configuración de registro ANSI predeterminada.

```bash
nextflow run hello-channels.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Volviendo a la vista de resumen, la salida se resume en una línea nuevamente.
Eche un vistazo al directorio `results` para ver si todos los saludos de salida están allí.

??? abstract "Contenido del directorio"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

¡Sí! Y cada uno tiene el contenido esperado.

??? abstract "Contenido del archivo"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

¡Éxito! Ahora podemos agregar tantos saludos como queramos sin preocuparnos de que los archivos de salida se sobrescriban.

!!! tip "Consejo"

    En la práctica, nombrar archivos basándose en los datos de entrada en sí es casi siempre poco práctico.
    La mejor manera de generar nombres de archivo dinámicos es pasar metadatos a un proceso junto con los archivos de entrada.
    Los metadatos generalmente se proporcionan a través de una 'hoja de muestra' o equivalentes.
    Aprenderá cómo hacer eso más adelante en su capacitación de Nextflow (consulte [Misión secundaria de metadatos](../side_quests/metadata.md)).

### Conclusión

Sabe cómo alimentar múltiples elementos de entrada a través de un canal.

### ¿Qué sigue?

Aprenda a usar un operador para transformar el contenido de un canal.

---

## 3. Proporcione múltiples entradas a través de un array

Acabamos de mostrarle cómo manejar múltiples elementos de entrada que estaban codificados directamente en el constructor de canales.
¿Qué pasa si quisiéramos proporcionar esas múltiples entradas de una manera diferente?

Por ejemplo, imagine que configuramos una variable de entrada que contiene un array de elementos como este:

`greetings_array = ['Hello','Bonjour','Holà']`

¿Podemos cargar eso en nuestro canal de salida y esperar que funcione?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Averigüémoslo.

### 3.1. Proporcione un array de valores como entrada al canal

El sentido común sugiere que deberíamos poder simplemente pasar un array de valores en lugar de un solo valor.
Intentémoslo; necesitaremos configurar la variable de entrada y cargarla en el constructor de canales.

#### 3.1.1. Configure la variable de entrada

Tomemos la variable `greetings_array` que acabamos de imaginar y hagámosla realidad agregándola al bloque workflow:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Esto aún no es funcional, solo hemos agregado una declaración para el array.

#### 3.1.2. Establezca el array de saludos como entrada al constructor de canales

Ahora vamos a reemplazar los valores `'Hello','Bonjour','Holà'` actualmente codificados en el constructor de canales con el `greetings_array` que acabamos de crear.

En el bloque workflow, realice el siguiente cambio:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Esto debería ser funcional ahora.

#### 3.1.3. Ejecute el workflow

Intentemos ejecutarlo:

```bash
nextflow run hello-channels.nf
```

??? failure "Salida del comando"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

¡Oh no! ¡Hay un error!

Observe la salida de `view()` y los mensajes de error.

Parece que Nextflow intentó ejecutar una sola llamada de proceso, usando `[Hello, Bonjour, Holà]` como un solo valor de cadena, en lugar de usar las tres cadenas en el array como valores separados.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Entonces es el 'empaquetado' lo que está causando el problema.
¿Cómo hacemos que Nextflow desempaquete el array y cargue las cadenas individuales en el canal?

### 3.2. Use un operador para transformar el contenido del canal

Aquí es donde entran en juego los [**operadores**](https://nextflow.io/docs/latest/reference/operator.html).
Ya ha usado el operador `.view()`, que solo mira lo que hay allí.
Ahora vamos a ver operadores que nos permiten actuar sobre el contenido de un canal.

Si hojea la [lista de operadores](https://nextflow.io/docs/latest/reference/operator.html) en la documentación de Nextflow, encontrará [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten), que hace exactamente lo que necesitamos: desempaquetar el contenido de un array y emitirlos como elementos individuales.

#### 3.2.1. Agregue el operador `flatten()`

Para aplicar el operador `flatten()` a nuestro canal de entrada, lo agregamos a la declaración del constructor de canales.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Aquí agregamos el operador en la siguiente línea para mayor legibilidad, pero puede agregar operadores en la misma línea que el constructor de canales si lo prefiere, así:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Refine las declaraciones `view()`

Podríamos ejecutar esto de inmediato para probar si funciona, pero mientras estamos en ello, vamos a refinar cómo inspeccionamos el contenido del canal.

Queremos poder contrastar cómo se ve el contenido antes y después de que se aplique el operador `flatten()`, por lo que vamos a agregar un segundo, Y vamos a agregar un poco de código para que se etiqueten más claramente en la salida.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Antes de flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Después de flatten: $greeting" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Verá que hemos agregado una segunda declaración `.view`, y para cada una de ellas, hemos reemplazado los paréntesis vacíos (`()`) con llaves que contienen algo de código, como `{ greeting -> "Antes de flatten: $greeting" }`.

Estos se llaman _closures_. El código que contienen se ejecutará para cada elemento en el canal.
Definimos una variable temporal para el valor interno, aquí llamada `greeting` (pero podría ser cualquier nombre arbitrario), que solo se usa dentro del alcance de ese closure.

En este ejemplo, `$greeting` representa cada elemento individual cargado en el canal.
Esto resultará en una salida de consola claramente etiquetada.

!!! info

    En algunos pipelines puede ver una variable especial llamada `$it` usada dentro de closures de operadores.
    Esta es una variable _implícita_ que permite un acceso abreviado a la variable interna,
    sin necesidad de definirla con un `->`.

    Preferimos ser explícitos para ayudar a la claridad del código, por lo que la sintaxis `$it` está desaconsejada y se eliminará lentamente del lenguaje Nextflow.

#### 3.2.3. Ejecute el workflow

Finalmente, ¡puede intentar ejecutar el workflow nuevamente!

```bash
nextflow run hello-channels.nf
```

??? success "Salida del comando"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Antes de flatten: [Hello, Bonjour, Holà]
    Después de flatten: Hello
    Después de flatten: Bonjour
    Después de flatten: Holà
    ```

Esta vez funciona Y nos da información adicional sobre cómo se ve el contenido del canal antes y después de ejecutar el operador `flatten()`.

- Una sola declaración `Antes de flatten:` porque en ese punto el canal contiene un elemento, el array original.
- Tres declaraciones separadas `Después de flatten:`, una para cada saludo, que ahora son elementos individuales en el canal.

Importante, esto significa que cada elemento ahora puede ser procesado por separado por el workflow.

!!! tip "Consejo"

    Técnicamente es posible lograr los mismos resultados usando un constructor de canales diferente, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), que incluye un paso de mapeo implícito en su operación.
    Aquí elegimos no usar eso para demostrar el uso de un operador en un caso de uso simple.

### Conclusión

Sabe cómo usar un operador como `flatten()` para transformar el contenido de un canal, y cómo usar el operador `view()` para inspeccionar el contenido del canal antes y después de aplicar un operador.

### ¿Qué sigue?

Aprenda a hacer que el workflow tome un archivo como su fuente de valores de entrada.

---

## 4. Lea valores de entrada desde un archivo CSV

Realísticamente, rara vez o nunca vamos a comenzar desde un array de valores.
Lo más probable es que tengamos uno o más archivos que contengan los datos que deben procesarse, en algún tipo de formato estructurado.

Hemos preparado un archivo CSV llamado `greetings.csv` que contiene varios saludos de entrada, imitando el tipo de datos en columnas que podría querer procesar en un análisis de datos real, almacenado en `data/`.
(Los números no son significativos, solo están ahí con fines ilustrativos.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Nuestra próxima tarea es adaptar nuestro workflow para leer los valores de este archivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Veamos cómo podemos hacer que eso suceda.

### 4.1. Modifique el script para esperar un archivo CSV como fuente de saludos

Para comenzar, vamos a necesitar hacer dos cambios clave en el script:

- Cambiar el parámetro de entrada para que apunte al archivo CSV
- Cambiar el constructor de canales a uno diseñado para manejar un archivo

#### 4.1.1. Cambie el parámetro de entrada para que apunte al archivo CSV

¿Recuerda el parámetro `params.input` que configuramos en la Parte 1?
Vamos a actualizarlo para que apunte al archivo CSV que contiene nuestros saludos.

Realice la siguiente edición a la declaración del parámetro:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Parámetros del pipeline
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Parámetros del pipeline
     */
    input: String = 'Holà mundo!'
    ```

Esto asume que el archivo está ubicado junto con el código del workflow.
Aprenderá cómo manejar otras ubicaciones de datos más adelante en su viaje con Nextflow.

#### 4.1.2. Cambie a un constructor de canales diseñado para manejar un archivo

Dado que ahora queremos usar un archivo en lugar de cadenas simples como entrada, no podemos usar el constructor de canales `channel.of()` de antes.
Necesitamos cambiar a usar un nuevo constructor de canales, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath), que tiene alguna funcionalidad incorporada para manejar rutas de archivos.

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Antes de flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Después de flatten: $greeting" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // declara un array de saludos de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // crea un canal para las entradas
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Antes de flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Después de flatten: $greeting" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Notará que cambiamos la entrada del canal de vuelta a `param.input`, y eliminamos la declaración `greetings_array` ya que ya no la necesitaremos.
También hemos comentado el `flatten()` y la segunda declaración `view()`.

#### 4.1.3. Ejecute el workflow

Intentemos ejecutar el workflow con el nuevo constructor de canales y el archivo de entrada.

```bash
nextflow run hello-channels.nf
```

??? failure "Salida del comando"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Antes de flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh no, no funciona. Eche un vistazo al inicio de la salida de la consola y el mensaje de error.
La parte `Command executed:` es especialmente útil aquí.

Esto puede parecer un poco familiar.
Parece que Nextflow intentó ejecutar una sola llamada de proceso usando la ruta del archivo en sí como un valor de cadena.
Entonces ha resuelto la ruta del archivo correctamente, pero en realidad no analizó su contenido, que es lo que queríamos.

¿Cómo hacemos que Nextflow abra el archivo y cargue su contenido en el canal?

¡Parece que necesitamos otro [operador](https://nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Use el operador `splitCsv()` para analizar el archivo

Mirando nuevamente la lista de operadores, encontramos [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), que está diseñado para analizar y dividir texto con formato CSV.

#### 4.2.1. Aplique `splitCsv()` al canal

Para aplicar el operador, lo agregamos a la línea del constructor de canales como antes.

En el bloque workflow, realice el siguiente cambio de código para reemplazar `flatten()` con `splitcsv()` (sin comentar):

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Antes de splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Después de splitCsv: $csv" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Antes de flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Después de flatten: $greeting" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Como puede ver, también actualizamos las declaraciones `view()` de antes/después.
Técnicamente podríamos haber usado el mismo nombre de variable (`greeting`) pero lo actualizamos a algo más apropiado (`csv`) para hacer el código más legible para otros.

#### 4.2.2. Ejecute el workflow nuevamente

Intentemos ejecutar el workflow con la lógica de análisis CSV agregada.

```bash
nextflow run hello-channels.nf
```

??? failure "Salida del comando"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Antes de splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Después de splitCsv: [Hello, English, 123]
    Después de splitCsv: [Bonjour, French, 456]
    Después de splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Curiosamente, esto también falla, pero con un error diferente.
Esta vez Nextflow ha analizado el contenido del archivo (¡bien!) pero ha cargado cada fila como un array, y cada array es un elemento en el canal.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Necesitamos decirle que solo tome la primera columna en cada fila.
Entonces, ¿cómo desempaquetamos esto?

Anteriormente usamos `flatten()` para desempaquetar el contenido de un canal, pero eso no funcionaría aquí porque flatten desempaqueta _todo_ (siéntase libre de intentarlo si quiere verlo por sí mismo).

En su lugar, usaremos otro operador llamado `map()` que es realmente útil y aparece mucho en los pipelines de Nextflow.

### 4.3. Use el operador `map()` para extraer los saludos

El operador [`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) es una pequeña herramienta muy útil que nos permite hacer todo tipo de mapeos al contenido de un canal.

En este caso, vamos a usarlo para extraer ese elemento que queremos de cada fila en nuestro archivo de datos.
Así es como se ve la sintaxis:

```groovy title="Sintaxis"
.map { row -> row[0] }
```

Esto significa 'para cada fila en el canal, tome el elemento 0 (primero) que contiene'.

Entonces apliquemos eso a nuestro análisis CSV.

#### 4.3.1. Aplique `map()` al canal

En el bloque workflow, realice el siguiente cambio de código:

=== "Después"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Antes de splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Después de splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "Después de map: $csv" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // crea un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Antes de splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Después de splitCsv: $csv" }
        // emite un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Verá que agregamos otra llamada `view()` para confirmar que el operador hace lo que esperamos.

#### 4.3.2. Ejecute el workflow

Ejecutemos esto una vez más:

```bash
nextflow run hello-channels.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Antes de splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Después de splitCsv: [Hello, English, 123]
    Después de splitCsv: [Bonjour, French, 456]
    Después de splitCsv: [Holà, Spanish, 789]
    Después de map: Hello
    Después de map: Bonjour
    Después de map: Holà
    ```

Esta vez debería ejecutarse sin error.

Mirando la salida de las declaraciones `view()`, verá lo siguiente:

- Una sola declaración `Antes de splitCsv:`: en ese punto el canal contiene un elemento, la ruta del archivo original.
- Tres declaraciones separadas `Después de splitCsv:`: una para cada saludo, pero cada una está contenida dentro de un array que corresponde a esa línea en el archivo.
- Tres declaraciones separadas `Después de map:`: una para cada saludo, que ahora son elementos individuales en el canal.

_Tenga en cuenta que las líneas pueden aparecer en un orden diferente en su salida._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

También puede mirar los archivos de salida para verificar que cada saludo fue extraído correctamente y procesado a través del workflow.

Hemos logrado el mismo resultado que antes, pero ahora tenemos mucha más flexibilidad para agregar más elementos al canal de saludos que queremos procesar modificando un archivo de entrada, sin modificar ningún código.
Aprenderá enfoques más sofisticados para manejar entradas complejas en una capacitación posterior.

### Conclusión

Sabe cómo usar el constructor de canales `.fromPath()` y los operadores `splitCsv()` y `map()` para leer un archivo de valores de entrada y manejarlos apropiadamente.

Más generalmente, tiene una comprensión básica de cómo Nextflow usa **canales** para administrar entradas a procesos y **operadores** para transformar su contenido.
También ha visto cómo los canales manejan la ejecución paralela implícitamente.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### ¿Qué sigue?

¡Tome un gran descanso, trabajó duro en este!

Cuando esté listo, continúe con [**Parte 3: Hello Workflow**](./03_hello_workflow.md) para aprender cómo agregar más pasos y conectarlos juntos en un workflow adecuado.

---

## Cuestionario

<quiz>
¿Qué es un canal en Nextflow?
- [ ] Una especificación de ruta de archivo
- [ ] Una definición de proceso
- [x] Una estructura tipo cola para pasar datos entre procesos
- [ ] Una configuración

Más información: [1.1. Cree un canal de entrada](#11-cree-un-canal-de-entrada)
</quiz>

<quiz>
¿Qué mostrará este código?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (una sola lista)
- [x] Cada elemento en una línea separada: `Hello`, `Bonjour`, `Hola`
- [ ] Nada (los canales no imprimen por defecto)
- [ ] Un error (sintaxis inválida)

Más información: [1.1. Cree un canal de entrada](#11-cree-un-canal-de-entrada)
</quiz>

<quiz>
Cuando un canal contiene múltiples valores, ¿cómo maneja Nextflow la ejecución del proceso?
- [ ] El proceso se ejecuta una vez con todos los valores
- [x] El proceso se ejecuta una vez por cada valor en el canal
- [ ] El proceso se ejecuta solo con el primer valor
- [ ] El proceso se ejecuta solo con el último valor

Más información: [2. Modifique el workflow para ejecutarse con múltiples valores de entrada](#2-modifique-el-workflow-para-ejecutarse-con-multiples-valores-de-entrada)
</quiz>

<quiz>
¿Qué hace el operador `flatten()`?
- [ ] Combina múltiples canales en uno
- [ ] Ordena elementos del canal
- [x] Desempaqueta arrays en elementos individuales
- [ ] Elimina elementos duplicados

Más información: [3.2.1. Agregue el operador `flatten()`](#321-agregue-el-operador-flatten)
</quiz>

<quiz>
¿Cuál es el propósito del operador `view()`?
- [ ] Filtrar el contenido del canal
- [ ] Transformar elementos del canal
- [x] Inspeccionar y depurar el contenido del canal
- [ ] Guardar el contenido del canal en un archivo

Más información: [1.4. Use `view()` para inspeccionar el contenido del canal](#14-use-view-para-inspeccionar-el-contenido-del-canal)
</quiz>

<quiz>
¿Qué hace `splitCsv()`?
- [ ] Crea un archivo CSV a partir del contenido del canal
- [ ] Divide una cadena por comas
- [x] Analiza un archivo CSV en arrays que representan cada fila
- [ ] Combina múltiples archivos CSV

Más información: [4.2. Use el operador `splitCsv()` para analizar el archivo](#42-use-el-operador-splitcsv-para-analizar-el-archivo)
</quiz>

<quiz>
¿Cuál es el propósito del operador `map()`?
- [ ] Filtrar elementos de un canal
- [ ] Combinar múltiples canales
- [x] Transformar cada elemento en un canal
- [ ] Contar elementos en un canal

Más información: [4.3. Use el operador `map()` para extraer los saludos](#43-use-el-operador-map-para-extraer-los-saludos)
</quiz>

<quiz>
¿Por qué es importante usar nombres de archivo de salida dinámicos al procesar múltiples entradas?
- [ ] Para mejorar el rendimiento
- [ ] Para reducir el espacio en disco
- [x] Para evitar que los archivos de salida se sobrescriban entre sí
- [ ] Para habilitar la funcionalidad de reanudación

Más información: [2.2. Asegúrese de que los nombres de los archivos de salida sean únicos](#22-asegurese-de-que-los-nombres-de-los-archivos-de-salida-sean-unicos)
</quiz>
