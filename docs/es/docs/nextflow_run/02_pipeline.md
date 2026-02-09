# Parte 2: Ejecutar pipelines reales

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la Parte 1 de este curso (Ejecutar Operaciones Básicas), comenzamos con un workflow de ejemplo que tenía solo características mínimas para mantener baja la complejidad del código.
Por ejemplo, `1-hello.nf` usaba un parámetro de línea de comandos (`--input`) para proporcionar un solo valor a la vez.

Sin embargo, la mayoría de los pipelines del mundo real usan características más sofisticadas para permitir el procesamiento eficiente de grandes cantidades de datos a escala, y aplican múltiples pasos de procesamiento encadenados por lógica a veces compleja.

En esta parte de la capacitación, demostramos características clave de pipelines del mundo real probando versiones expandidas del pipeline original Hello World.

## 1. Procesar datos de entrada desde un archivo

En un pipeline del mundo real, típicamente queremos procesar múltiples puntos de datos (o series de datos) contenidos en uno o más archivos de entrada.
Y donde sea posible, queremos ejecutar el procesamiento de datos independientes en paralelo, para acortar el tiempo de espera del análisis.

Para demostrar cómo Nextflow hace esto, hemos preparado un archivo CSV llamado `greetings.csv` que contiene varios saludos de entrada, imitando el tipo de datos en columnas que podrías querer procesar en un análisis de datos real.
Nota que los números no son significativos, están ahí solo con fines ilustrativos.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

También hemos escrito una versión mejorada del workflow original, ahora llamado `2a-inputs.nf`, que leerá el archivo CSV, extraerá los saludos y escribirá cada uno de ellos en un archivo separado.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Ejecutemos primero el workflow, y luego veremos el código Nextflow relevante.

### 1.1. Ejecutar el workflow

Ejecuta el siguiente comando en tu terminal.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Emocionantemente, esto parece indicar que se realizaron '3 de 3' llamadas para el proceso, lo cual es alentador, ya que había tres filas de datos en el CSV que proporcionamos como entrada.
Esto sugiere que el proceso `sayHello()` fue llamado tres veces, una vez en cada fila de entrada.

### 1.2. Encontrar las salidas publicadas en el directorio `results`

Veamos el directorio 'results' para ver si nuestro workflow todavía está escribiendo una copia de nuestras salidas allí.

??? abstract "Contenido del directorio"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

¡Sí! Vemos un nuevo directorio llamado `2a-inputs` con tres archivos de salida con diferentes nombres, convenientemente.

Puedes abrir cada uno de ellos para verificar que contienen la cadena de saludo apropiada.

??? abstract "Contenido de los archivos"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Esto confirma que cada saludo en el archivo de entrada ha sido procesado apropiadamente.

### 1.3. Encontrar las salidas originales y los logs

Puede que hayas notado que la salida de consola anterior se refería a solo un directorio de tarea.
¿Significa eso que las tres llamadas a `sayHello()` fueron ejecutadas dentro de ese único directorio de tarea?

#### 1.3.1. Examinar el directorio de tarea dado en la terminal

Echemos un vistazo dentro de ese directorio de tarea `8e/0eb066`.

??? abstract "Contenido del directorio"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Solo encontramos la salida correspondiente a uno de los saludos (así como los archivos accesorios si habilitamos la visualización de archivos ocultos).

Entonces, ¿qué está pasando aquí?

Por defecto, el sistema de logging ANSI escribe la información de estado para todas las llamadas al mismo proceso en la misma línea.
Como resultado, solo nos mostró una de las tres rutas de directorio de tarea (`8e/0eb066`) en la salida de consola.
Hay otras dos que no están listadas allí.

#### 1.3.2. Hacer que la terminal muestre más detalles

Podemos modificar el comportamiento del logging para ver la lista completa de llamadas a procesos agregando `-ansi-log false` al comando de la siguiente manera:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Salida del comando"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Esta vez vemos las tres ejecuciones de proceso y sus subdirectorios de trabajo asociados listados en la salida.
Deshabilitar el logging ANSI también evitó que Nextflow usara colores en la salida de terminal.

Nota que la forma en que se reporta el estado es un poco diferente entre los dos modos de logging.
En el modo condensado, Nextflow reporta si las llamadas se completaron exitosamente o no.
En este modo expandido, solo reporta que fueron enviadas.

Esto confirma que el proceso `sayHello()` se llama tres veces, y se crea un directorio de tarea separado para cada uno.

Si miramos dentro de cada uno de los directorios de tarea listados allí, podemos verificar que cada uno corresponde a uno de los saludos.

??? abstract "Contenido del directorio"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Esto confirma que cada llamada a proceso se ejecuta aisladamente de todas las demás.
Eso tiene muchas ventajas, incluyendo evitar colisiones si el proceso produce archivos intermedios con nombres no únicos.

!!! tip

    Para un workflow complejo, o un gran número de entradas, tener la lista completa en la terminal podría ser un poco abrumador, por lo que la gente normalmente no usa `-ansi-log false` en el uso rutinario.

### 1.4. Examinar el código del workflow

Entonces esta versión del workflow es capaz de leer un archivo CSV de entradas, procesar las entradas por separado y nombrar las salidas de manera única.

Echemos un vistazo a qué hace esto posible en el código del workflow.

??? full-code "Archivo de código completo"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // crear un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emitir un saludo
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Una vez más, no necesitas memorizar la sintaxis del código, pero es bueno aprender a reconocer componentes clave del workflow que proporcionan funcionalidad importante.

#### 1.4.1. Cargar los datos de entrada desde el CSV

Esta es la parte más interesante: ¿cómo cambiamos de tomar un solo valor desde la línea de comandos, a tomar un archivo CSV, parsearlo y procesar los saludos individuales que contiene?

En Nextflow, hacemos eso con un [**channel**](https://nextflow.io/docs/latest/channel.html): una construcción de cola diseñada para manejar entradas eficientemente y transportarlas de un paso a otro en workflows de múltiples pasos, mientras proporciona paralelismo incorporado y muchos beneficios adicionales.

Vamos a desglosarlo.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // crear un canal para entradas desde un archivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emitir un saludo
    sayHello(greeting_ch)
```

Este código crea un canal llamado `greeting_ch` que lee el archivo CSV, lo parsea y extrae la primera columna de cada fila.
El resultado es un canal que contiene `Hello`, `Bonjour` y `Holà`.

??? tip "¿Cómo funciona esto?"

    Aquí está lo que esa línea significa en lenguaje simple:

    - `channel.fromPath` es una **channel factory** que crea un canal desde ruta(s) de archivo
    - `(params.input)` especifica que la ruta de archivo es proporcionada por `--input` en la línea de comandos

    En otras palabras, esa línea le dice a Nextflow: toma la ruta de archivo dada con `--input` y prepárate para tratar su contenido como datos de entrada.

    Luego las siguientes dos líneas aplican **operators** que hacen el parseo real del archivo y la carga de los datos en la estructura de datos apropiada:

    - `.splitCsv()` le dice a Nextflow que parsee el archivo CSV en un array que representa filas y columnas
    - `.map { line -> line[0] }` le dice a Nextflow que tome solo el elemento en la primera columna de cada fila

    Entonces en la práctica, comenzando desde el siguiente archivo CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Hemos transformado eso en un array que se ve así:

    ```txt title="Contenido del array"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    Y luego hemos tomado el primer elemento de cada una de las tres filas y los hemos cargado en un canal de Nextflow que ahora contiene: `Hello`, `Bonjour` y `Holà`.

    Si quieres entender los canales y operadores en profundidad, incluyendo cómo escribirlos tú mismo, consulta [Hello Nextflow Parte 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Llamar al proceso en cada saludo

A continuación, en la última línea del bloque `main:` del workflow, proporcionamos el canal `greeting_ch` cargado como entrada al proceso `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // crear un canal para entradas desde un archivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emitir un saludo
    sayHello(greeting_ch)
```

Esto le dice a Nextflow que ejecute el proceso individualmente en cada elemento del canal, _es decir_, en cada saludo.
Y porque Nextflow es inteligente así, ejecutará estas llamadas a proceso en paralelo si es posible, dependiendo de la infraestructura de cómputo disponible.

Así es como puedes lograr procesamiento eficiente y escalable de muchos datos (muchas muestras, o puntos de datos, cualquiera que sea tu unidad de investigación) con comparativamente muy poco código.

#### 1.4.3. Cómo se nombran las salidas

Finalmente, vale la pena echar un vistazo rápido al código del proceso para ver cómo logramos que los archivos de salida se nombren de manera única.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

Ves que, comparado con la versión de este proceso en `1-hello.nf`, la declaración de salida y la parte relevante del comando han cambiado para incluir el valor del saludo en el nombre del archivo de salida.

Esta es una forma de asegurar que los nombres de archivo de salida no colisionen cuando se publiquen en el directorio de resultados común.

¡Y ese es el único cambio que hemos tenido que hacer dentro de la declaración del proceso!

### Conclusión

Entiendes a un nivel básico cómo los canales y operadores nos permiten procesar múltiples entradas eficientemente.

### ¿Qué sigue?

Descubre cómo se construyen los workflows de múltiples pasos y cómo operan.

---

## 2. Ejecutar workflows de múltiples pasos

La mayoría de los workflows del mundo real involucran más de un paso.
Construyamos sobre lo que acabamos de aprender sobre canales, y veamos cómo Nextflow usa canales y operadores para conectar procesos juntos en un workflow de múltiples pasos.

Para ese fin, te proporcionamos un workflow de ejemplo que encadena tres pasos separados y demuestra lo siguiente:

1. Hacer que los datos fluyan de un proceso al siguiente
2. Recolectar salidas de múltiples llamadas a proceso en una sola llamada a proceso

Específicamente, hicimos una versión expandida del workflow llamada `2b-multistep.nf` que toma cada saludo de entrada, lo convierte a mayúsculas, luego recolecta todos los saludos en mayúsculas en un solo archivo de salida.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Como antes, ejecutaremos primero el workflow y luego veremos el código para ver qué es nuevo.

### 2.1. Ejecutar el workflow

Ejecuta el siguiente comando en tu terminal:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Salida del comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Ves que como se prometió, se ejecutaron múltiples pasos como parte del workflow; los primeros dos (`sayHello` y `convertToUpper`) presumiblemente se ejecutaron en cada saludo individual, y el tercero (`collectGreetings`) se habrá ejecutado solo una vez, en las salidas de las tres llamadas a `convertToUpper`.

### 2.2. Encontrar las salidas

Verifiquemos que eso es de hecho lo que sucedió echando un vistazo al directorio `results`.

??? abstract "Contenido del directorio"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Como puedes ver, tenemos un nuevo directorio llamado `2b-multistep`, y contiene bastantes más archivos que antes.
Algunos de los archivos han sido agrupados en un subdirectorio llamado `intermediates`, mientras que dos archivos están ubicados en el nivel superior.

Esos dos son los resultados finales del workflow de múltiples pasos.
Tómate un minuto para ver los nombres de archivo y verificar su contenido para confirmar que son lo que esperas.

??? abstract "Contenido de los archivos"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

El primero contiene nuestros tres saludos, en mayúsculas y recolectados de vuelta en un solo archivo como se prometió.
El segundo es un archivo de reporte que resume información sobre la ejecución.

### 2.3. Examinar el código

Veamos el código e identifiquemos los patrones clave para workflows de múltiples pasos.

??? full-code "Archivo de código completo"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // crear un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emitir un saludo
        sayHello(greeting_ch)
        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)
        // recolectar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Hay mucho sucediendo ahí, pero la diferencia más obvia comparada con la versión anterior del workflow es que ahora hay múltiples definiciones de proceso, y correspondientemente, varias llamadas a proceso en el bloque workflow.

Echemos un vistazo más de cerca y veamos si podemos identificar las piezas más interesantes.

#### 2.3.1. Visualizar la estructura del workflow

Si estás usando VSCode con la extensión de Nextflow, puedes obtener un diagrama útil de cómo los procesos están conectados haciendo clic en el pequeño enlace `DAG preview` mostrado justo encima del bloque workflow en cualquier script de Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Esto te da una buena visión general de cómo los procesos están conectados y qué producen.

Ves que además del proceso original `sayHello`, ahora también tenemos `convertToUpper` y `collectGreetings`, que coinciden con los nombres de los procesos que vimos en la salida de consola.
Las dos nuevas definiciones de proceso están estructuradas de la misma manera que el proceso `sayHello`, excepto que `collectGreetings` toma un parámetro de entrada adicional llamado `batch` y produce dos salidas.

No entraremos en detalle en el código de cada uno, pero si tienes curiosidad, puedes buscar los detalles en [Parte 2 de Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Por ahora, profundicemos en cómo los procesos están conectados entre sí.

#### 2.3.2. Cómo se conectan los procesos

Lo realmente interesante para ver aquí es cómo las llamadas a proceso están encadenadas juntas en el bloque `main:` del workflow.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // crear un canal para entradas desde un archivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emitir un saludo
    sayHello(greeting_ch)
    // convertir el saludo a mayúsculas
    convertToUpper(sayHello.out)
    // recolectar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Puedes ver que la primera llamada a proceso, `sayHello(greeting_ch)`, no ha cambiado.
Luego la siguiente llamada a proceso, a `convertToUpper`, se refiere a la salida de `sayHello` como `sayHello.out`.

El patrón es simple: `processName.out` se refiere al canal de salida de un proceso, que puede ser pasado directamente al siguiente proceso.
Así es como transportamos datos de un paso al siguiente en Nextflow.

#### 2.3.3. Un proceso puede tomar múltiples entradas

La tercera llamada a proceso, a `collectGreetings`, es un poco diferente.

```groovy title="2b-multistep.nf" linenums="77"
    // recolectar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Ves que esta llamada recibe dos entradas, `convertToUpper.out.collect()` y `params.batch`.
Ignorando la parte `.collect()` por ahora, podemos generalizar esto como `collectGreetings(input1, input2)`.

Eso coincide con las dos declaraciones de entrada en el módulo de proceso:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Cuando Nextflow parsea esto, asignará la primera entrada en la llamada a `path input_files`, y la segunda a `val batch_name`.

Entonces ahora sabes que un proceso puede tomar múltiples entradas, y cómo se ve la llamada en el bloque workflow.

Ahora echemos un vistazo más de cerca a esa primera entrada, `convertToUpper.out.collect()`.

#### 2.3.4. Qué hace `collect()` en la llamada a `collectGreetings`

Para pasar la salida de `sayHello` a `convertToUpper`, simplemente nos referimos al canal de salida de `sayHello` como `sayHello.out`. Pero para el siguiente paso, estamos viendo una referencia a `convertToUpper.out.collect()`.

¿Qué es esta parte de `collect()` y qué hace?

Es un operador, por supuesto. Al igual que los operadores `splitCsv` y `map` que encontramos antes.
Esta vez el operador se llama `collect`, y se aplica al canal de salida producido por `convertToUpper`.

El operador `collect` se usa para recolectar las salidas de múltiples llamadas al mismo proceso y empaquetarlas en un solo elemento de canal.

En el contexto de este workflow, está tomando los tres saludos en mayúsculas en el canal `convertToUpper.out` (que son tres elementos de canal separados, y normalmente serían manejados en llamadas separadas por el siguiente proceso) y empaquetándolos en un solo elemento.
Así es como obtenemos todos los saludos de vuelta en el mismo archivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

En contraste, si no aplicáramos `collect()` a la salida de `convertToUpper()` antes de alimentarla a `collectGreetings()`, Nextflow simplemente ejecutaría `collectGreetings()` independientemente en cada saludo, lo que no lograría nuestro objetivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Hay muchos otros [operators](https://nextflow.io/docs/latest/reference/operator.html) disponibles para aplicar transformaciones al contenido de los canales entre llamadas a proceso.

Esto da a los desarrolladores de pipelines mucha flexibilidad para personalizar la lógica de flujo de su pipeline.
La desventaja es que a veces puede hacer más difícil descifrar qué está haciendo el pipeline.

#### 2.3.5. Un parámetro de entrada puede tener un valor predeterminado

Puede que hayas notado que `collectGreetings` toma una segunda entrada, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // recolectar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Esto pasa un parámetro CLI llamado `--batch` al workflow.
Sin embargo, cuando lanzamos el workflow antes, no especificamos un parámetro `--batch`.

¿Qué está pasando ahí?
Echa un vistazo al bloque `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Hay un valor predeterminado configurado en el workflow, por lo que no tenemos que proporcionarlo.
Pero si proporcionamos uno en la línea de comandos, el valor que especifiquemos se usará en lugar del predeterminado.

Pruébalo:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Salida del comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Deberías ver nuevas salidas finales nombradas con tu nombre de lote personalizado.

??? abstract "Contenido del directorio"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Este es un aspecto de la configuración de entrada, que cubriremos con más detalle en la Parte 3, pero por ahora lo importante es saber que los parámetros de entrada pueden tener valores predeterminados.

#### 2.3.6. Un proceso puede producir múltiples salidas

En la definición del proceso `collectGreetings`, vemos las siguientes declaraciones de salida:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Que luego se refieren por el nombre dado con `emit:` en el bloque `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Esto hace fácil pasar salidas específicas individualmente a otros procesos en el workflow, en combinación con varios operadores.

#### 2.3.7. Las salidas publicadas pueden ser organizadas

En el bloque `output`, hemos usado rutas personalizadas para agrupar resultados intermedios con el fin de hacer más fácil seleccionar solo las salidas finales del workflow.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Hay formas más sofisticadas de organizar salidas publicadas; tocaremos algunas en la parte sobre configuración.

!!! tip "¿Quieres aprender más sobre construir workflows?"

    Para cobertura detallada de construir workflows de múltiples pasos, consulta [Hello Nextflow Parte 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Conclusión

Entiendes a un nivel básico cómo se construyen los workflows de múltiples pasos usando canales y operadores y cómo operan.
También has visto que los procesos pueden tomar múltiples entradas y producir múltiples salidas, y que estas pueden ser publicadas de manera estructurada.

### ¿Qué sigue?

Aprende cómo los pipelines de Nextflow pueden ser modularizados para promover la reutilización de código y la mantenibilidad.

---

## 3. Ejecutar pipelines modularizados

Hasta ahora, todos los workflows que hemos visto han consistido en un solo archivo de workflow que contiene todo el código relevante.

Sin embargo, los pipelines del mundo real típicamente se benefician de ser _modularizados_, lo que significa que el código se divide en diferentes archivos.
Esto puede hacer su desarrollo y mantenimiento más eficiente y sostenible.

Aquí vamos a demostrar la forma más común de modularidad de código en Nextflow, que es el uso de **modules**.

En Nextflow, un [**module**](https://nextflow.io/docs/latest/module.html) es una definición de proceso única que está encapsulada por sí misma en un archivo de código independiente.
Para usar un módulo en un workflow, solo agregas una declaración de importación de una sola línea a tu archivo de código de workflow; luego puedes integrar el proceso en el workflow de la misma manera que normalmente lo harías.
Eso hace posible reutilizar definiciones de proceso en múltiples workflows sin producir múltiples copias del código.

Hasta ahora hemos estado ejecutando workflows que tenían todos sus procesos incluidos en un archivo de código monolítico.
Ahora vamos a ver cómo se ve cuando los procesos están almacenados en módulos individuales.

Por supuesto, hemos preparado nuevamente un workflow adecuado para propósitos de demostración, llamado `2c-modules.nf`, junto con un conjunto de módulos ubicados en el directorio `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Contenido del directorio"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Ves que hay cuatro archivos Nextflow, cada uno nombrado según uno de los procesos.
Puedes ignorar el archivo `cowpy.nf` por ahora; llegaremos a ese más tarde.

### 3.1. Examinar el código

Esta vez vamos a ver el código primero.
Comienza abriendo el archivo de workflow `2c-modules.nf`.

??? full-code "Archivo de código completo"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // crear un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emitir un saludo
        sayHello(greeting_ch)
        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)
        // recolectar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Ves que la lógica del workflow es exactamente la misma que en la versión anterior del workflow.
Sin embargo, el código del proceso ha desaparecido del archivo de workflow, y en su lugar hay declaraciones `include` apuntando a archivos separados bajo `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Incluir módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Abre uno de esos archivos y encontrarás el código para el proceso correspondiente.

??? full-code "Archivo de código completo"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

Como puedes ver, el código del proceso no ha cambiado; solo ha sido copiado en un archivo de módulo individual en lugar de estar en el archivo de workflow principal.
Lo mismo aplica a los otros dos procesos.

Entonces veamos cómo se ve ejecutar esta nueva versión.

### 3.2. Ejecutar el workflow

Ejecuta este comando en tu terminal, con la bandera `-resume`:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Notarás que las ejecuciones de proceso todas se almacenaron en caché exitosamente, lo que significa que Nextflow reconoció que ya ha hecho el trabajo solicitado, aunque el código ha sido dividido y el archivo de workflow principal ha sido renombrado.

Nada de eso le importa a Nextflow; lo que importa es el script de trabajo que se genera una vez que todo el código ha sido reunido y evaluado.

!!! tip

    También es posible encapsular una sección de un workflow como un 'subworkflow' que puede ser importado en un pipeline más grande, pero eso está fuera del alcance de este curso.

    Puedes aprender más sobre desarrollar workflows componibles en la Misión Secundaria sobre [Workflows de Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Conclusión

Sabes cómo los procesos pueden ser almacenados en módulos independientes para promover la reutilización de código y mejorar la mantenibilidad.

### ¿Qué sigue?

Aprende a usar contenedores para gestionar dependencias de software.

---

## 4. Usar software en contenedores

Hasta ahora los workflows que hemos estado usando como ejemplos solo necesitaban ejecutar operaciones de procesamiento de texto muy básicas usando herramientas UNIX disponibles en nuestro entorno.

Sin embargo, los pipelines del mundo real típicamente requieren herramientas y paquetes especializados que no están incluidos por defecto en la mayoría de los entornos.
Usualmente, necesitarías instalar estas herramientas, gestionar sus dependencias y resolver cualquier conflicto.

Todo eso es muy tedioso y molesto.
Una forma mucho mejor de abordar este problema es usar **containers**.

Un **container** es una unidad de software ligera, independiente y ejecutable creada a partir de una **image** de contenedor que incluye todo lo necesario para ejecutar una aplicación incluyendo código, bibliotecas del sistema y configuraciones.

!!! Tip

    Enseñamos esto usando la tecnología [Docker](https://www.docker.com/get-started/), pero Nextflow soporta varias otras tecnologías de contenedores también.
    Puedes aprender más sobre el soporte de Nextflow para contenedores [aquí](https://nextflow.io/docs/latest/container.html).

### 4.1. Usar un contenedor directamente

Primero, intentemos interactuar con un contenedor directamente.
Esto ayudará a solidificar tu comprensión de qué son los contenedores antes de que comencemos a usarlos en Nextflow.

#### 4.1.1. Descargar la imagen del contenedor

Para usar un contenedor, usualmente descargas o "pulls" una imagen de contenedor desde un registro de contenedores, y luego ejecutas la imagen de contenedor para crear una instancia de contenedor.

La sintaxis general es la siguiente:

```bash title="Sintaxis"
docker pull '<container>'
```

- `docker pull` es la instrucción al sistema de contenedores para descargar una imagen de contenedor desde un repositorio.
- `'<container>'` es la dirección URI de la imagen de contenedor.

Como ejemplo, descarguemos una imagen de contenedor que contiene [cowpy](https://github.com/jeffbuttars/cowpy), una implementación en python de una herramienta llamada `cowsay` que genera arte ASCII para mostrar entradas de texto arbitrarias de una manera divertida.

Hay varios repositorios donde puedes encontrar contenedores publicados.
Usamos el servicio [Seqera Containers](https://seqera.io/containers/) para generar esta imagen de contenedor Docker desde el paquete Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Ejecuta el comando pull completo:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Salida del comando"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Esto le dice al sistema que descargue la imagen especificada.
Una vez que la descarga está completa, tienes una copia local de la imagen de contenedor.

#### 4.1.2. Iniciar el contenedor

Los contenedores pueden ejecutarse como un comando único, pero también puedes usarlos interactivamente, lo que te da un prompt de shell dentro del contenedor y te permite jugar con el comando.

La sintaxis general es la siguiente:

```bash title="Sintaxis"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` es la instrucción al sistema de contenedores para iniciar una instancia de contenedor desde una imagen de contenedor y ejecutar un comando en ella.
- `--rm` le dice al sistema que apague la instancia de contenedor después de que el comando se haya completado.

Completamente ensamblado, el comando de ejecución de contenedor se ve así:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Ejecuta ese comando, y deberías ver tu prompt cambiar a algo como `(base) root@b645838b3314:/tmp#`, lo que indica que ahora estás dentro del contenedor.

Puedes verificar esto ejecutando `ls` para listar el contenido del directorio:

```bash
ls /
```

??? success "Salida del comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Ves que el sistema de archivos dentro del contenedor es diferente del sistema de archivos en tu sistema host.

!!! Tip

    Cuando ejecutas un contenedor, está aislado del sistema host por defecto.
    Esto significa que el contenedor no puede acceder a ningún archivo en el sistema host a menos que explícitamente le permitas hacerlo especificando que quieres montar un volumen como parte del comando `docker run` usando la siguiente sintaxis:

    ```bash title="Sintaxis"
    -v <outside_path>:<inside_path>
    ```

    Esto efectivamente establece un túnel a través de la pared del contenedor que puedes usar para acceder a esa parte de tu sistema de archivos.

    Esto se cubre con más detalle en [Parte 5 de Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Ejecutar la herramienta `cowpy`

Desde dentro del contenedor, puedes ejecutar el comando `cowpy` directamente.

```bash
cowpy "Hello Containers"
```

??? success "Salida del comando"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Esto produce arte ASCII del personaje de vaca predeterminado (o 'cowacter') con una burbuja de diálogo que contiene el texto que especificamos.

Ahora que has probado el uso básico, puedes intentar darle algunos parámetros.
Por ejemplo, la documentación de la herramienta dice que podemos establecer el personaje con `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Salida del comando"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Esta vez la salida de arte ASCII muestra el pingüino de Linux, Tux, porque especificamos el parámetro `-c tux`.

Como estás dentro del contenedor, puedes ejecutar el comando cowpy tantas veces como quieras, variando los parámetros de entrada, sin tener que preocuparte por instalar ninguna biblioteca en tu sistema mismo.

??? tip "Otros personajes disponibles"

    Usa la bandera '-c' para elegir un personaje diferente, incluyendo:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Siéntete libre de jugar con esto.
Cuando hayas terminado, sal del contenedor usando el comando `exit`:

```bash
exit
```

Te encontrarás de vuelta en tu shell normal.

### 4.2. Usar un contenedor en un workflow

Cuando ejecutamos un pipeline, queremos poder decirle a Nextflow qué contenedor usar en cada paso, e importantemente, queremos que maneje todo ese trabajo que acabamos de hacer: descargar el contenedor, iniciarlo, ejecutar el comando y desmontar el contenedor cuando haya terminado.

Buenas noticias: eso es exactamente lo que Nextflow va a hacer por nosotros.
Solo necesitamos especificar un contenedor para cada proceso.

Para demostrar cómo funciona esto, hicimos otra versión de nuestro workflow que ejecuta `cowpy` en el archivo de saludos recolectados producido en el tercer paso.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Esto debería producir un archivo que contiene el arte ASCII con los tres saludos en la burbuja de diálogo.

#### 4.2.1. Examinar el código

El workflow es muy similar al anterior, más el paso extra para ejecutar `cowpy`.

??? full-code "Archivo de código completo"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Incluir módulos
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // crear un canal para entradas desde un archivo CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emitir un saludo
        sayHello(greeting_ch)
        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)
        // recolectar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Ves que este workflow importa un proceso `cowpy` desde un archivo de módulo, y lo llama en la salida de la llamada a `collectGreetings()`, más un parámetro de entrada llamado `params.character`.

```groovy title="2d-container.nf" linenums="31"
// generar arte ASCII de los saludos con cowpy
cowpy(collectGreetings.out.outfile, params.character)
```

El proceso `cowpy`, que envuelve el comando cowpy para generar arte ASCII, está definido en el módulo `cowpy.nf`.

??? full-code "Archivo de código completo"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
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

El proceso `cowpy` requiere dos entradas: la ruta a un archivo de entrada que contiene el texto para poner en la burbuja de diálogo (`input_file`), y un valor para la variable character.

Importantemente, también incluye la línea `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, que apunta al URI del contenedor que usamos antes.

#### 4.2.2. Verificar que Docker está habilitado en la configuración

Vamos a anticipar ligeramente la Parte 3 de este curso de capacitación introduciendo el archivo de configuración `nextflow.config`, que es una de las principales formas que Nextflow ofrece para configurar la ejecución del workflow.
Cuando un archivo llamado `nextflow.config` está presente en el directorio actual, Nextflow lo cargará automáticamente y aplicará cualquier configuración que contenga.

Para ese fin, incluimos un archivo `nextflow.config` con una sola línea de código que habilita Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Esta configuración le dice a Nextflow que use Docker para cualquier proceso que especifique un contenedor compatible.

!!! tip

    Es técnicamente posible habilitar la ejecución de Docker desde la línea de comandos, por ejecución, usando el parámetro `-with-docker <container>`.
    Sin embargo, eso solo nos permite especificar un contenedor para todo el workflow, mientras que el enfoque que acabamos de mostrarte nos permite especificar un contenedor diferente por proceso.
    Este último es mucho mejor para modularidad, mantenimiento de código y reproducibilidad.

#### 4.2.3. Ejecutar el workflow

Solo para recapitular, esto es lo que estamos a punto de ejecutar:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

¿Crees que va a funcionar?

Ejecutemos el workflow con la bandera `-resume`, y especifiquemos que queremos que el personaje sea el pavo.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

Los primeros tres pasos se almacenaron en caché ya que ya los hemos ejecutado antes, pero el proceso `cowpy` es nuevo así que ese realmente se ejecuta.

Puedes encontrar la salida del paso `cowpy` en el directorio `results`.

??? abstract "Contenido del archivo"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
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

Ves que el personaje está diciendo todos los saludos, ya que se ejecutó en el archivo de saludos en mayúsculas recolectados.

Más al punto, pudimos ejecutar esto como parte de nuestro pipeline sin tener que hacer una instalación apropiada de cowpy y todas sus dependencias.
Y ahora podemos compartir el pipeline con colaboradores y hacer que lo ejecuten en su infraestructura sin que ellos necesiten instalar nada tampoco, aparte de Docker o una de sus alternativas (como Singularity/Apptainer) como se mencionó anteriormente.

#### 4.2.4. Inspeccionar cómo Nextflow lanzó la tarea en contenedor

Como coda final a esta sección, echemos un vistazo al subdirectorio de trabajo para una de las llamadas al proceso `cowpy` para obtener un poco más de información sobre cómo Nextflow trabaja con contenedores bajo el capó.

Verifica la salida de tu comando `nextflow run` para encontrar la ruta al subdirectorio de trabajo para el proceso `cowpy`.
Mirando lo que obtuvimos para la ejecución mostrada arriba, la línea de log de consola para el proceso `cowpy` comienza con `[7f/caf718]`.
Eso corresponde a la siguiente ruta de directorio truncada: `work/7f/caf718`.

En ese directorio, encontrarás el archivo `.command.run` que contiene todos los comandos que Nextflow ejecutó en tu nombre en el curso de ejecutar el pipeline.

??? abstract "Contenido del archivo"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Si buscas `nxf_launch` en este archivo, deberías ver algo como esto:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Este comando de lanzamiento muestra que Nextflow está usando un comando `docker run` muy similar para lanzar la llamada al proceso como lo hicimos cuando lo ejecutamos manualmente.
También monta el subdirectorio de trabajo correspondiente en el contenedor, establece el directorio de trabajo dentro del contenedor en consecuencia, y ejecuta nuestro script bash con plantilla en el archivo `.command.sh`.

¡Esto confirma que todo el trabajo duro que tuvimos que hacer manualmente en la sección anterior ahora es hecho por nosotros por Nextflow!

### Conclusión

Entiendes qué papel juegan los contenedores en la gestión de versiones de herramientas de software y en asegurar la reproducibilidad.

Más generalmente, tienes una comprensión básica de cuáles son los componentes centrales de los pipelines de Nextflow del mundo real y cómo están organizados.
Conoces los fundamentos de cómo Nextflow puede procesar múltiples entradas eficientemente, ejecutar workflows compuestos de múltiples pasos conectados juntos, aprovechar componentes de código modulares, y utilizar contenedores para mayor reproducibilidad y portabilidad.

### ¿Qué sigue?

¡Toma otro descanso! Esa fue una gran pila de información sobre cómo funcionan los pipelines de Nextflow.

En la última sección de esta capacitación, vamos a profundizar más en el tema de configuración.
Aprenderás cómo configurar la ejecución de tu pipeline para ajustarse a tu infraestructura así como gestionar la configuración de entradas y parámetros.

---

## Quiz

<quiz>
¿Por qué Nextflow crea un directorio de tarea separado para cada llamada a proceso?
- [ ] Para mejorar la velocidad de ejecución
- [ ] Para reducir el uso de memoria
- [x] Para aislar ejecuciones y evitar colisiones entre salidas
- [ ] Para habilitar compresión de archivos en paralelo

Aprende más: [1.3. Encontrar las salidas originales y los logs](#13-encontrar-las-salidas-originales-y-los-logs)
</quiz>

<quiz>
¿Qué hace la opción `-ansi-log false` al ejecutar un workflow?
- [ ] Deshabilita toda la salida de consola
- [x] Elimina el color de la salida
- [x] Muestra todas las rutas de directorio de tarea en lugar de condensarlas en una línea
- [ ] Habilita el modo de depuración verbose

Aprende más: [1.3.2. Hacer que la terminal muestre más detalles](#132-hacer-que-la-terminal-muestre-más-detalles)

También puedes usar cualquiera de las siguientes variables de entorno si prefieres este estilo:

```bash
export NXF_ANSI_LOG=0
# o
export NO_COLOR=1
```

</quiz>

<quiz>
En el código `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, ¿qué hace `#!groovy .map { line -> line[0] }`?
- [ ] Filtra líneas vacías
- [ ] Ordena las líneas alfabéticamente
- [x] Extrae la primera columna de cada fila CSV
- [ ] Cuenta el número de líneas

Aprende más: [1.4.1. Cargar los datos de entrada desde el CSV](#141-cargar-los-datos-de-entrada-desde-el-csv)
</quiz>

<quiz>
¿Por qué es importante incluir el valor de entrada en los nombres de archivo de salida (por ejemplo, `#!groovy "${greeting}-output.txt"`)?
- [ ] Para mejorar la velocidad de procesamiento
- [ ] Para habilitar la funcionalidad de resume
- [x] Para evitar que los archivos de salida se sobrescriban entre sí al procesar múltiples entradas
- [ ] Para hacer los archivos más fáciles de comprimir

Aprende más: [1.4.3. Cómo se nombran las salidas](#143-cómo-se-nombran-las-salidas)
</quiz>

<quiz>
¿Cuál es el propósito de la declaración `include` en un workflow modularizado?
- [ ] Para copiar código de proceso en el archivo de workflow
- [x] Para importar una definición de proceso desde un archivo de módulo externo
- [ ] Para incluir configuraciones
- [ ] Para agregar comentarios de documentación

Aprende más: [3. Ejecutar pipelines modularizados](#3-ejecutar-pipelines-modularizados)
</quiz>

<quiz>
Cuando modularizas un workflow y lo ejecutas con `-resume`, ¿qué sucede?
- [ ] El almacenamiento en caché se deshabilita para procesos modulares
- [ ] Todas las tareas deben ser re-ejecutadas
- [x] El almacenamiento en caché funciona normalmente basado en los scripts de trabajo generados
- [ ] Solo el archivo de workflow principal se almacena en caché

Aprende más: [3.2. Ejecutar el workflow](#32-ejecutar-el-workflow)
</quiz>

<quiz>
¿Qué especifica la directiva `container` en una definición de proceso?
- [ ] El directorio de trabajo para el proceso
- [ ] La asignación máxima de memoria
- [x] El URI de la imagen de contenedor a usar para ejecutar el proceso
- [ ] El formato de archivo de salida

Aprende más: [4.2. Usar un contenedor en un workflow](#42-usar-un-contenedor-en-un-workflow)
</quiz>

<quiz>
En el archivo `.command.run`, ¿qué contiene la función `nxf_launch`?
- [ ] La información de versión de Nextflow
- [ ] Los parámetros del workflow
- [x] El comando `docker run` con montajes de volumen y configuraciones de contenedor
- [ ] Las declaraciones de entrada del proceso

Aprende más: [4.2.4. Inspeccionar cómo Nextflow lanzó la tarea en contenedor](#424-inspeccionar-cómo-nextflow-lanzó-la-tarea-en-contenedor)
</quiz>

<quiz>
¿Qué maneja automáticamente Nextflow al ejecutar un proceso en contenedor? (Selecciona todas las que apliquen)
- [x] Descargar la imagen de contenedor si es necesario
- [x] Montar el directorio de trabajo en el contenedor
- [x] Ejecutar el script del proceso dentro del contenedor
- [x] Limpiar la instancia de contenedor después de la ejecución

Aprende más: [4. Usar software en contenedores](#4-usar-software-en-contenedores)
</quiz>
