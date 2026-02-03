# Parte 2: Ejecutar pipelines reales

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la Parte 1 de este curso (Ejecutar operaciones básicas), comenzamos con un workflow de ejemplo que solo tenía características mínimas para mantener baja la complejidad del código.
Por ejemplo, `1-hello.nf` usaba un parámetro de línea de comandos (`--input`) para proporcionar un único valor a la vez.

Sin embargo, la mayoría de los pipelines del mundo real usan características más sofisticadas para permitir el procesamiento eficiente de grandes cantidades de datos a escala, y aplican múltiples pasos de procesamiento encadenados por una lógica a veces compleja.

En esta parte del entrenamiento, demostramos características clave de pipelines del mundo real probando versiones expandidas del pipeline Hello World original.

## 1. Procesar datos de entrada desde un archivo

En un pipeline del mundo real, típicamente queremos procesar múltiples puntos de datos (o series de datos) contenidos en uno o más archivos de entrada.
Y siempre que sea posible, queremos ejecutar el procesamiento de datos independientes en paralelo, para acortar el tiempo de espera para el análisis.

Para demostrar cómo hace esto Nextflow, hemos preparado un archivo CSV llamado `greetings.csv` que contiene varios saludos de entrada, imitando el tipo de datos columnares que podría querer procesar en un análisis de datos real.
Note que los números no tienen significado, están ahí solo con propósitos ilustrativos.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

También hemos escrito una versión mejorada del workflow original, ahora llamada `2a-inputs.nf`, que leerá el archivo CSV, extraerá los saludos y escribirá cada uno de ellos en un archivo separado.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Ejecutemos el workflow primero, y veremos el código de Nextflow relevante después.

### 1.1. Ejecutar el workflow

Ejecute el siguiente comando en su terminal.

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

Emocionantemente, esto parece indicar que se hicieron '3 of 3' llamadas para el process, lo cual es alentador, ya que había tres filas de datos en el CSV que proporcionamos como entrada.
Esto sugiere que el process `sayHello()` fue llamado tres veces, una vez en cada fila de entrada.

### 1.2. Encontrar las salidas publicadas en el directorio `results`

Veamos el directorio 'results' para ver si nuestro workflow todavía está escribiendo una copia de nuestras salidas allí.

??? abstract "Contenidos del directorio"

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

Puede abrir cada uno de ellos para asegurarse de que contienen la cadena de saludo apropiada.

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

### 1.3. Encontrar las salidas y registros originales

Puede haber notado que la salida de consola anterior se refería a solo un directorio de tarea.
¿Significa eso que las tres llamadas a `sayHello()` fueron ejecutadas dentro de ese único directorio de tarea?

#### 1.3.1. Examinar el directorio de tarea dado en la terminal

Echemos un vistazo dentro de ese directorio de tarea `8e/0eb066`.

??? abstract "Contenidos del directorio"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Solo encontramos la salida correspondiente a uno de los saludos (así como los archivos accesorios si habilitamos la visualización de archivos ocultos).

Entonces, ¿qué está pasando aquí?

Por defecto, el sistema de registro ANSI escribe la información de estado para todas las llamadas al mismo process en la misma línea.
Como resultado, solo nos mostró una de las tres rutas de directorio de tarea (`8e/0eb066`) en la salida de consola.
Hay otras dos que no están listadas allí.

#### 1.3.2. Hacer que la terminal muestre más detalles

Podemos modificar el comportamiento de registro para ver la lista completa de llamadas de process agregando `-ansi-log false` al comando de la siguiente manera:

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

Esta vez vemos las tres ejecuciones de process y sus subdirectorios de trabajo asociados listados en la salida.
Deshabilitar el registro ANSI también evitó que Nextflow usara colores en la salida de terminal.

Note que la forma en que se reporta el estado es un poco diferente entre los dos modos de registro.
En el modo condensado, Nextflow reporta si las llamadas se completaron exitosamente o no.
En este modo expandido, solo reporta que fueron enviadas.

Esto confirma que el process `sayHello()` es llamado tres veces, y se crea un directorio de tarea separado para cada uno.

Si miramos dentro de cada uno de los directorios de tarea listados allí, podemos verificar que cada uno corresponde a uno de los saludos.

??? abstract "Contenidos del directorio"

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

Esto confirma que cada llamada de process se ejecuta aislada de todas las demás.
Eso tiene muchas ventajas, incluyendo evitar colisiones si el process produce algún archivo intermedio con nombres no únicos.

!!! tip "Consejo"

    Para un workflow complejo, o un gran número de entradas, tener la lista completa en la terminal podría volverse un poco abrumador, así que las personas normalmente no usan `-ansi-log false` en uso rutinario.

### 1.4. Examinar el código del workflow

Así que esta versión del workflow es capaz de leer un archivo CSV de entradas, procesar las entradas por separado, y nombrar las salidas de manera única.

Veamos qué hace posible eso en el código del workflow.

??? full-code "Archivo de código completo"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Usar echo para imprimir 'Hello World!' a un archivo
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

Una vez más, no necesita memorizar la sintaxis del código, pero es bueno aprender a reconocer componentes clave del workflow que proporcionan funcionalidad importante.

#### 1.4.1. Cargar los datos de entrada desde el CSV

Esta es la parte más interesante: ¿cómo cambiamos de tomar un único valor desde la línea de comandos, a tomar un archivo CSV, analizarlo y procesar los saludos individuales que contiene?

En Nextflow, hacemos eso con un [**channel**](https://nextflow.io/docs/latest/channel.html): una estructura de cola diseñada para manejar entradas eficientemente y transportarlas de un paso a otro en workflows de múltiples pasos, mientras proporciona paralelismo incorporado y muchos beneficios adicionales.

Desglosémoslo.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // crear un canal para entradas desde un archivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emitir un saludo
    sayHello(greeting_ch)
```

Este código crea un channel llamado `greeting_ch` que lee el archivo CSV, lo analiza y extrae la primera columna de cada fila.
El resultado es un channel que contiene `Hello`, `Bonjour` y `Holà`.

??? tip "¿Cómo funciona esto?"

    Esto es lo que significa esa línea en español sencillo:

    - `channel.fromPath` es una **fábrica de canales** que crea un channel a partir de ruta(s) de archivo
    - `(params.input)` especifica que la ruta del archivo se proporciona mediante `--input` en la línea de comandos

    En otras palabras, esa línea le dice a Nextflow: toma la ruta del archivo dada con `--input` y prepárate para tratar su contenido como datos de entrada.

    Luego las siguientes dos líneas aplican **operadores** que hacen el análisis real del archivo y la carga de los datos en la estructura de datos apropiada:

    - `.splitCsv()` le dice a Nextflow que analice el archivo CSV en un array representando filas y columnas
    - `.map { line -> line[0] }` le dice a Nextflow que tome solo el elemento en la primera columna de cada fila

    Así que en la práctica, comenzando desde el siguiente archivo CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Hemos transformado eso en un array que se ve así:

    ```txt title="Contenido del array"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    Y luego hemos tomado el primer elemento de cada una de las tres filas y los hemos cargado en un channel de Nextflow que ahora contiene: `Hello`, `Bonjour` y `Holà`.

    Si quiere entender los channels y operadores en profundidad, incluyendo cómo escribirlos usted mismo, vea [Hello Nextflow Part 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Llamar al process en cada saludo

Luego, en la última línea del bloque `main:` del workflow, proporcionamos el channel cargado `greeting_ch` como entrada al process `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // crear un canal para entradas desde un archivo CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emitir un saludo
    sayHello(greeting_ch)
```

Esto le dice a Nextflow que ejecute el process individualmente en cada elemento del channel, _es decir_, en cada saludo.
Y como Nextflow es inteligente, ejecutará estas llamadas de process en paralelo si es posible, dependiendo de la infraestructura de computación disponible.

Así es como puede lograr un procesamiento eficiente y escalable de muchos datos (muchas muestras, o puntos de datos, lo que sea su unidad de investigación) con comparativamente muy poco código.

#### 1.4.3. Cómo se nombran las salidas

Finalmente, vale la pena echar un vistazo rápido al código del process para ver cómo hacemos que los archivos de salida se nombren de manera única.

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

Puede ver que, comparado con la versión de este process en `1-hello.nf`, la declaración de salida y la parte relevante del comando han cambiado para incluir el valor del saludo en el nombre del archivo de salida.

Esta es una forma de asegurar que los nombres de los archivos de salida no colisionen cuando se publican en el directorio de resultados común.

¡Y ese es el único cambio que hemos tenido que hacer dentro de la declaración del process!

### Conclusión

Comprende a un nivel básico cómo los channels y operadores nos permiten procesar múltiples entradas eficientemente.

### ¿Qué sigue?

Descubra cómo se construyen los workflows de múltiples pasos y cómo operan.

---

## 2. Ejecutar workflows de múltiples pasos

La mayoría de los workflows del mundo real involucran más de un paso.
Construyamos sobre lo que acabamos de aprender sobre channels, y veamos cómo Nextflow usa channels y operadores para conectar processes juntos en un workflow de múltiples pasos.

Para ese fin, le proporcionamos un workflow de ejemplo que encadena tres pasos separados y demuestra lo siguiente:

1. Hacer que los datos fluyan de un process al siguiente
2. Recolectar salidas de múltiples llamadas de process en una sola llamada de process

Específicamente, hicimos una versión expandida del workflow llamada `2b-multistep.nf` que toma cada saludo de entrada, lo convierte a mayúsculas, luego recolecta todos los saludos en mayúsculas en un único archivo de salida.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Como anteriormente, ejecutaremos el workflow primero y luego veremos el código para ver qué hay de nuevo.

### 2.1. Ejecutar el workflow

Ejecute el siguiente comando en su terminal:

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

Puede ver que como se prometió, se ejecutaron múltiples pasos como parte del workflow; los primeros dos (`sayHello` y `convertToUpper`) presumiblemente se ejecutaron en cada saludo individual, y el tercero (`collectGreetings`) se habrá ejecutado solo una vez, en las salidas de las tres llamadas `convertToUpper`.

### 2.2. Encontrar las salidas

Verifiquemos que eso es de hecho lo que sucedió echando un vistazo en el directorio `results`.

??? abstract "Contenidos del directorio"

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

Como puede ver, tenemos un nuevo directorio llamado `2b-multistep`, y contiene bastantes más archivos que antes.
Algunos de los archivos han sido agrupados en un subdirectorio llamado `intermediates`, mientras que dos archivos están ubicados en el nivel superior.

Esos dos son los resultados finales del workflow de múltiples pasos.
Tómese un minuto para mirar los nombres de los archivos y verificar su contenido para confirmar que son lo que espera.

??? abstract "Contenido de los archivos"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

El primero contiene nuestros tres saludos, en mayúsculas y recolectados de vuelta en un único archivo como se prometió.
El segundo es un archivo de informe que resume alguna información sobre la ejecución.

### 2.3. Examinar el código

Veamos el código e identifiquemos los patrones clave para workflows de múltiples pasos.

??? full-code "Archivo de código completo"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Usar echo para imprimir 'Hello World!' a un archivo
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
    * Usar una herramienta de reemplazo de texto para convertir el saludo a mayúsculas
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
    * Recopilar saludos en mayúsculas en un único archivo de salida
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
        // recopilar todos los saludos en un archivo
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

Hay mucho sucediendo allí, pero la diferencia más obvia comparada con la versión anterior del workflow es que ahora hay múltiples definiciones de process, y correspondientemente, varias llamadas de process en el bloque workflow.

Echemos un vistazo más de cerca y veamos si podemos identificar las piezas más interesantes.

#### 2.3.1. Visualizar la estructura del workflow

Si está usando VSCode con la extensión de Nextflow, puede obtener un diagrama útil de cómo los processes están conectados haciendo clic en el pequeño enlace `DAG preview` que se muestra justo encima del bloque workflow en cualquier script de Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Esto le da una buena visión general de cómo los processes están conectados y qué producen.

Puede ver que además del process original `sayHello`, ahora también tenemos `convertToUpper` y `collectGreetings`, que coinciden con los nombres de los processes que vimos en la salida de consola.
Las dos nuevas definiciones de process están estructuradas de la misma manera que el process `sayHello`, excepto que `collectGreetings` toma un parámetro de entrada adicional llamado `batch` y produce dos salidas.

No entraremos en el código para cada uno en detalle, pero si tiene curiosidad, puede buscar los detalles en [Part 2 of Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Por ahora, profundicemos en cómo los processes están conectados entre sí.

#### 2.3.2. Cómo los processes están conectados

Lo realmente interesante de ver aquí es cómo las llamadas de process están encadenadas juntas en el bloque `main:` del workflow.

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
    // recopilar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Puede ver que la primera llamada de process, `sayHello(greeting_ch)`, no ha cambiado.
Luego la siguiente llamada de process, a `convertToUpper`, se refiere a la salida de `sayHello` como `sayHello.out`.

El patrón es simple: `processName.out` se refiere al channel de salida de un process, que puede pasarse directamente al siguiente process.
Así es como transportamos datos de un paso al siguiente en Nextflow.

#### 2.3.3. Un process puede tomar múltiples entradas

La tercera llamada de process, a `collectGreetings`, es un poco diferente.

```groovy title="2b-multistep.nf" linenums="77"
    // recopilar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Puede ver que esta llamada recibe dos entradas, `convertToUpper.out.collect()` y `params.batch`.
Ignorando la parte `.collect()` por ahora, podemos generalizar esto como `collectGreetings(input1, input2)`.

Eso coincide con las dos declaraciones de entrada en el módulo del process:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Cuando Nextflow analiza esto, asignará la primera entrada en la llamada a `path input_files`, y la segunda a `val batch_name`.

Así que ahora sabe que un process puede tomar múltiples entradas, y cómo se ve la llamada en el bloque workflow.

Ahora echemos un vistazo más de cerca a esa primera entrada, `convertToUpper.out.collect()`.

#### 2.3.4. Qué hace `collect()` en la llamada `collectGreetings`

Para pasar la salida de `sayHello` a `convertToUpper`, simplemente nos referimos al channel de salida de `sayHello` como `sayHello.out`. Pero para el siguiente paso, estamos viendo una referencia a `convertToUpper.out.collect()`.

¿Qué es esta parte de `collect()` y qué hace?

Es un operador, por supuesto. Al igual que los operadores `splitCsv` y `map` que encontramos antes.
Esta vez el operador se llama `collect`, y se aplica al channel de salida producido por `convertToUpper`.

El operador `collect` se usa para recolectar las salidas de múltiples llamadas al mismo process y empaquetarlas en un único elemento de channel.

En el contexto de este workflow, está tomando los tres saludos en mayúsculas en el channel `convertToUpper.out` (que son tres elementos de channel separados, y normalmente serían manejados en llamadas separadas por el siguiente process) y empaquetándolos en un único elemento.
Así es como obtenemos todos los saludos de vuelta en el mismo archivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

En contraste, si no aplicáramos `collect()` a la salida de `convertToUpper()` antes de alimentarla a `collectGreetings()`, Nextflow simplemente ejecutaría `collectGreetings()` independientemente en cada saludo, lo cual no lograría nuestro objetivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Hay muchos otros [operadores](https://nextflow.io/docs/latest/reference/operator.html) disponibles para aplicar transformaciones al contenido de los channels entre llamadas de process.

Esto le da a los desarrolladores de pipelines mucha flexibilidad para personalizar la lógica de flujo de su pipeline.
La desventaja es que a veces puede hacer más difícil descifrar lo que está haciendo el pipeline.

#### 2.3.5. Un parámetro de entrada puede tener un valor predeterminado

Puede haber notado que `collectGreetings` toma una segunda entrada, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // recopilar todos los saludos en un archivo
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Esto pasa un parámetro CLI llamado `--batch` al workflow.
Sin embargo, cuando lanzamos el workflow antes, no especificamos un parámetro `--batch`.

¿Qué está pasando ahí?
Eche un vistazo al bloque `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Hay un valor predeterminado configurado en el workflow, así que no tenemos que proporcionarlo.
Pero si proporcionamos uno en la línea de comandos, el valor que especificamos se usará en lugar del predeterminado.

Inténtelo:

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

Debería ver nuevas salidas finales nombradas con su nombre de lote personalizado.

??? abstract "Contenidos del directorio"

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

#### 2.3.6. Un process puede producir múltiples salidas

En la definición del process `collectGreetings`, vemos las siguientes declaraciones de salida:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Las cuales luego se refieren por el nombre dado con `emit:` en el bloque `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Esto facilita pasar salidas específicas individualmente a otros processes en el workflow, en combinación con varios operadores.

#### 2.3.7. Las salidas publicadas pueden organizarse

En el bloque `output`, hemos usado rutas personalizadas para agrupar resultados intermedios para facilitar identificar solo las salidas finales del workflow.

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

Hay formas más sofisticadas de organizar las salidas publicadas; tocaremos algunas en la parte sobre configuración.

!!! tip "¿Quiere aprender más sobre construir workflows?"

    Para cobertura detallada de la construcción de workflows de múltiples pasos, vea [Hello Nextflow Part 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Conclusión

Comprende a un nivel básico cómo se construyen los workflows de múltiples pasos usando channels y operadores y cómo operan.
También ha visto que los processes pueden tomar múltiples entradas y producir múltiples salidas, y que estas pueden publicarse de manera estructurada.

### ¿Qué sigue?

Aprenda cómo los pipelines de Nextflow pueden modularizarse para promover la reutilización de código y la mantenibilidad.

---

## 3. Ejecutar pipelines modularizados

Hasta ahora, todos los workflows que hemos visto han consistido en un único archivo de workflow que contiene todo el código relevante.

Sin embargo, los pipelines del mundo real típicamente se benefician de ser _modularizados_, lo que significa que el código se divide en diferentes archivos.
Esto puede hacer su desarrollo y mantenimiento más eficiente y sostenible.

Aquí vamos a demostrar la forma más común de modularidad de código en Nextflow, que es el uso de **módulos**.

En Nextflow, un [**módulo**](https://nextflow.io/docs/latest/module.html) es una definición de process única que se encapsula por sí misma en un archivo de código independiente.
Para usar un módulo en un workflow, simplemente agrega una declaración de importación de una sola línea a su archivo de código de workflow; luego puede integrar el process en el workflow de la misma manera que normalmente lo haría.
Eso hace posible reutilizar definiciones de process en múltiples workflows sin producir múltiples copias del código.

Hasta ahora hemos estado ejecutando workflows que tenían todos sus processes incluidos en un archivo de código monolítico.
Ahora vamos a ver cómo se ve cuando los processes se almacenan en módulos individuales.

Por supuesto, hemos preparado nuevamente un workflow adecuado para propósitos de demostración, llamado `2c-modules.nf`, junto con un conjunto de módulos ubicados en el directorio `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Contenidos del directorio"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Puede ver que hay cuatro archivos de Nextflow, cada uno nombrado según uno de los processes.
Puede ignorar el archivo `cowpy.nf` por ahora; llegaremos a ese más tarde.

### 3.1. Examinar el código

Esta vez vamos a ver el código primero.
Comience abriendo el archivo de workflow `2c-modules.nf`.

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
        // recopilar todos los saludos en un archivo
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

Puede ver que la lógica del workflow es exactamente la misma que en la versión anterior del workflow.
Sin embargo, el código del process ya no está en el archivo del workflow, y en su lugar hay declaraciones `include` que apuntan a archivos separados bajo `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Incluir módulos
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Abra uno de esos archivos y encontrará el código para el process correspondiente.

??? full-code "Archivo de código completo"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usar echo para imprimir 'Hello World!' a un archivo
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

Como puede ver, el código del process no ha cambiado; simplemente se ha copiado en un archivo de módulo individual en lugar de estar en el archivo principal del workflow.
Lo mismo aplica a los otros dos processes.

Así que veamos cómo se ve ejecutar esta nueva versión.

### 3.2. Ejecutar el workflow

Ejecute este comando en su terminal, con la bandera `-resume`:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
    [95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
    [5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
    ```

Notará que las ejecuciones de process todas se almacenaron en caché exitosamente, lo que significa que Nextflow reconoció que ya ha hecho el trabajo solicitado, aunque el código se haya dividido y el archivo principal del workflow se haya renombrado.

Nada de eso importa a Nextflow; lo que importa es el script de trabajo que se genera una vez que todo el código ha sido juntado y evaluado.

!!! tip "Consejo"

    También es posible encapsular una sección de un workflow como un 'subworkflow' que puede importarse en un pipeline más grande, pero eso está fuera del alcance de este curso.

    Puede aprender más sobre desarrollar workflows componibles en el Side Quest sobre [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Conclusión

Sabe cómo los processes pueden almacenarse en módulos independientes para promover la reutilización de código y mejorar la mantenibilidad.

### ¿Qué sigue?

Aprenda a usar contenedores para gestionar dependencias de software.

---

## 4. Usar software en contenedores

Hasta ahora los workflows que hemos usado como ejemplos solo necesitaban ejecutar operaciones de procesamiento de texto muy básicas usando herramientas UNIX disponibles en nuestro entorno.

Sin embargo, los pipelines del mundo real típicamente requieren herramientas y paquetes especializados que no están incluidos por defecto en la mayoría de los entornos.
Usualmente, necesitaría instalar estas herramientas, gestionar sus dependencias y resolver cualquier conflicto.

Todo eso es muy tedioso y molesto.
Una manera mucho mejor de abordar este problema es usar **contenedores**.

Un **contenedor** es una unidad de software ligera, independiente y ejecutable creada a partir de una **imagen** de contenedor que incluye todo lo necesario para ejecutar una aplicación incluyendo código, bibliotecas del sistema y configuraciones.

!!! Tip "Consejo"

    Enseñamos esto usando la tecnología [Docker](https://www.docker.com/get-started/), pero Nextflow soporta varias otras tecnologías de contenedores también.
    Puede aprender más sobre el soporte de Nextflow para contenedores [aquí](https://nextflow.io/docs/latest/container.html).

### 4.1. Usar un contenedor directamente

Primero, intentemos interactuar con un contenedor directamente.
Esto ayudará a solidificar su comprensión de qué son los contenedores antes de que comencemos a usarlos en Nextflow.

#### 4.1.1. Descargar la imagen del contenedor

Para usar un contenedor, usualmente descarga o "pull" una imagen de contenedor de un registro de contenedores, y luego ejecuta la imagen del contenedor para crear una instancia de contenedor.

La sintaxis general es la siguiente:

```bash title="Sintaxis"
docker pull '<container>'
```

- `docker pull` es la instrucción al sistema de contenedores para descargar una imagen de contenedor de un repositorio.
- `'<container>'` es la dirección URI de la imagen del contenedor.

Como ejemplo, descarguemos una imagen de contenedor que contiene [cowpy](https://github.com/jeffbuttars/cowpy), una implementación en python de una herramienta llamada `cowsay` que genera arte ASCII para mostrar entradas de texto arbitrarias de una manera divertida.

Hay varios repositorios donde puede encontrar contenedores publicados.
Usamos el servicio [Seqera Containers](https://seqera.io/containers/) para generar esta imagen de contenedor Docker del paquete Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Ejecute el comando completo de descarga:

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
Una vez que la descarga esté completa, tiene una copia local de la imagen del contenedor.

#### 4.1.2. Iniciar el contenedor

Los contenedores pueden ejecutarse como un comando único, pero también puede usarlos interactivamente, lo que le da un prompt de shell dentro del contenedor y le permite jugar con el comando.

La sintaxis general es la siguiente:

```bash title="Sintaxis"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` es la instrucción al sistema de contenedores para iniciar una instancia de contenedor a partir de una imagen de contenedor y ejecutar un comando en él.
- `--rm` le dice al sistema que apague la instancia del contenedor después de que el comando se haya completado.

Completamente ensamblado, el comando de ejecución del contenedor se ve así:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Ejecute ese comando, y debería ver que su prompt cambia a algo como `(base) root@b645838b3314:/tmp#`, lo que indica que ahora está dentro del contenedor.

Puede verificar esto ejecutando `ls` para listar contenidos del directorio:

```bash
ls /
```

??? success "Salida del comando"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Puede ver que el sistema de archivos dentro del contenedor es diferente del sistema de archivos en su sistema host.

!!! Tip "Consejo"

    Cuando ejecuta un contenedor, está aislado del sistema host por defecto.
    Esto significa que el contenedor no puede acceder a ningún archivo en el sistema host a menos que explícitamente lo permita especificando que quiere montar un volumen como parte del comando `docker run` usando la siguiente sintaxis:

    ```bash title="Sintaxis"
    -v <outside_path>:<inside_path>
    ```

    Esto efectivamente establece un túnel a través de la pared del contenedor que puede usar para acceder a esa parte de su sistema de archivos.

    Esto se cubre con más detalle en [Part 5 of Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Ejecutar la herramienta `cowpy`

Desde dentro del contenedor, puede ejecutar el comando `cowpy` directamente.

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

Esto produce arte ASCII del personaje de vaca predeterminado (o 'cowacter') con un globo de diálogo que contiene el texto que especificamos.

Ahora que ha probado el uso básico, puede intentar darle algunos parámetros.
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

Ya que está dentro del contenedor, puede ejecutar el comando cowpy tantas veces como quiera, variando los parámetros de entrada, sin tener que preocuparse por instalar ninguna biblioteca en su sistema mismo.

??? tip "Otros personajes disponibles"

    Use la bandera '-c' para elegir un personaje diferente, incluyendo:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Siéntase libre de jugar con esto.
Cuando haya terminado, salga del contenedor usando el comando `exit`:

```bash
exit
```

Se encontrará de vuelta en su shell normal.

### 4.2. Usar un contenedor en un workflow

Cuando ejecutamos un pipeline, queremos poder decirle a Nextflow qué contenedor usar en cada paso, e importantemente, queremos que maneje todo ese trabajo que acabamos de hacer: descargar el contenedor, iniciarlo, ejecutar el comando y destruir el contenedor cuando haya terminado.

Buenas noticias: eso es exactamente lo que Nextflow va a hacer por nosotros.
Solo necesitamos especificar un contenedor para cada process.

Para demostrar cómo funciona esto, hicimos otra versión de nuestro workflow que ejecuta `cowpy` en el archivo de saludos recolectados producido en el tercer paso.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Esto debería producir un archivo que contiene el arte ASCII con los tres saludos en el globo de diálogo.

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
        // recopilar todos los saludos en un archivo
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

Puede ver que este workflow importa un process `cowpy` de un archivo de módulo, y lo llama en la salida de la llamada `collectGreetings()`, más un parámetro de entrada llamado `params.character`.

```groovy title="2d-container.nf" linenums="25"
// generar arte ASCII con cowpy
cowpy(collectGreetings.out, params.character)
```

El process `cowpy`, que envuelve el comando cowpy para generar arte ASCII, está definido en el módulo `cowpy.nf`.

??? full-code "Archivo de código completo"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generar arte ASCII con cowpy (https://github.com/jeffbuttars/cowpy)
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

El process `cowpy` requiere dos entradas: la ruta a un archivo de entrada que contiene el texto para poner en el globo de diálogo (`input_file`), y un valor para la variable de personaje.

Importantemente, también incluye la línea `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, que apunta a la URI del contenedor que usamos antes.

#### 4.2.2. Verificar que Docker esté habilitado en la configuración

Vamos a anticipar ligeramente la Parte 3 de este curso de entrenamiento introduciendo el archivo de configuración `nextflow.config`, que es una de las principales formas que ofrece Nextflow para configurar la ejecución de workflows.
Cuando un archivo llamado `nextflow.config` está presente en el directorio actual, Nextflow lo cargará automáticamente y aplicará cualquier configuración que contenga.

Para ese fin, incluimos un archivo `nextflow.config` con una sola línea de código que habilita Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Esta configuración le dice a Nextflow que use Docker para cualquier process que especifique un contenedor compatible.

!!! tip "Consejo"

    Es técnicamente posible habilitar la ejecución de Docker desde la línea de comandos, por ejecución, usando el parámetro `-with-docker <container>`.
    Sin embargo, eso solo nos permite especificar un contenedor para todo el workflow, mientras que el enfoque que acabamos de mostrarle le permite especificar un contenedor diferente por process.
    Este último es mucho mejor para modularidad, mantenimiento de código y reproducibilidad.

#### 4.2.3. Ejecutar el workflow

Solo para recapitular, esto es lo que estamos a punto de ejecutar:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

¿Cree que va a funcionar?

Ejecutemos el workflow con la bandera `-resume`, y especifiquemos que queremos que el personaje sea el pavo (turkey).

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

Los primeros tres pasos se almacenaron en caché ya que los hemos ejecutado antes, pero el process `cowpy` es nuevo así que realmente se ejecuta.

Puede encontrar la salida del paso `cowpy` en el directorio `results`.

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

Puede ver que el personaje está diciendo todos los saludos, ya que se ejecutó en el archivo de saludos recolectados en mayúsculas.

Más al punto, pudimos ejecutar esto como parte de nuestro pipeline sin tener que hacer una instalación propiamente dicha de cowpy y todas sus dependencias.
Y ahora podemos compartir el pipeline con colaboradores y hacer que lo ejecuten en su infraestructura sin que ellos necesiten instalar nada tampoco, aparte de Docker o una de sus alternativas (como Singularity/Apptainer) como se mencionó anteriormente.

#### 4.2.4. Inspeccionar cómo Nextflow lanzó la tarea en contenedor

Como coda final a esta sección, echemos un vistazo al subdirectorio de trabajo para una de las llamadas del process `cowpy` para obtener un poco más de información sobre cómo trabaja Nextflow con contenedores bajo el capó.

Verifique la salida de su comando `nextflow run` para encontrar la ruta al subdirectorio de trabajo para el process `cowpy`.
Mirando lo que obtuvimos para la ejecución mostrada arriba, la línea de registro de consola para el process `cowpy` comienza con `[7f/caf718]`.
Eso corresponde a la siguiente ruta de directorio truncada: `work/7f/caf718`.

En ese directorio, encontrará el archivo `.command.run` que contiene todos los comandos que Nextflow ejecutó en su nombre en el curso de ejecutar el pipeline.

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
    ...
    ```

Si busca `nxf_launch` en este archivo, debería ver algo como esto:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Este comando de lanzamiento muestra que Nextflow está usando un comando `docker run` muy similar para lanzar la llamada del process como hicimos cuando lo ejecutamos manualmente.
También monta el subdirectorio de trabajo correspondiente en el contenedor, establece el directorio de trabajo dentro del contenedor en consecuencia, y ejecuta nuestro script bash con plantilla en el archivo `.command.sh`.

¡Esto confirma que todo el trabajo duro que tuvimos que hacer manualmente en la sección anterior ahora lo hace Nextflow por nosotros!

### Conclusión

Comprende qué papel juegan los contenedores en la gestión de versiones de herramientas de software y asegurar la reproducibilidad.

Más generalmente, tiene una comprensión básica de cuáles son los componentes principales de los pipelines de Nextflow del mundo real y cómo están organizados.
Conoce los fundamentos de cómo Nextflow puede procesar múltiples entradas eficientemente, ejecutar workflows compuestos de múltiples pasos conectados juntos, aprovechar componentes de código modulares, y utilizar contenedores para mayor reproducibilidad y portabilidad.

### ¿Qué sigue?

¡Tome otro descanso! Esa fue una gran cantidad de información sobre cómo funcionan los pipelines de Nextflow.

En la última sección de este entrenamiento, vamos a profundizar más en el tema de la configuración.
Aprenderá cómo configurar la ejecución de su pipeline para adaptarlo a su infraestructura así como gestionar la configuración de entradas y parámetros.

---

## Cuestionario

<quiz>
¿Por qué Nextflow crea un directorio de tarea separado para cada llamada de process?
- [ ] Para mejorar la velocidad de ejecución
- [ ] Para reducir el uso de memoria
- [x] Para aislar ejecuciones y evitar colisiones entre salidas
- [ ] Para habilitar la compresión paralela de archivos

Más información: [1.3. Encontrar las salidas y registros originales](#13-encontrar-las-salidas-y-registros-originales)
</quiz>

<quiz>
¿Qué hace la opción `-ansi-log false` al ejecutar un workflow?
- [ ] Deshabilita toda la salida de consola
- [x] Elimina el color de la salida
- [x] Muestra todas las rutas de directorio de tarea en lugar de condensarlas en una línea
- [ ] Habilita el modo de depuración detallado

Más información: [1.3.2. Hacer que la terminal muestre más detalles](#132-hacer-que-la-terminal-muestre-más-detalles)

También puede usar cualquiera de las siguientes variables de entorno si prefiere este estilo:

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

Más información: [1.4.1. Cargar los datos de entrada desde el CSV](#141-cargar-los-datos-de-entrada-desde-el-csv)
</quiz>

<quiz>
¿Por qué es importante incluir el valor de entrada en los nombres de archivos de salida (por ejemplo, `#!groovy "${greeting}-output.txt"`)?
- [ ] Para mejorar la velocidad de procesamiento
- [ ] Para habilitar la funcionalidad resume
- [x] Para evitar que los archivos de salida se sobrescriban entre sí al procesar múltiples entradas
- [ ] Para facilitar la compresión de archivos

Más información: [1.4.3. Cómo se nombran las salidas](#143-cómo-se-nombran-las-salidas)
</quiz>

<quiz>
¿Cuál es el propósito de la declaración `include` en un workflow modularizado?
- [ ] Copiar código de process al archivo de workflow
- [x] Importar una definición de process de un archivo de módulo externo
- [ ] Incluir configuraciones
- [ ] Agregar comentarios de documentación

Más información: [3. Ejecutar pipelines modularizados](#3-ejecutar-pipelines-modularizados)
</quiz>

<quiz>
Cuando modulariza un workflow y lo ejecuta con `-resume`, ¿qué sucede?
- [ ] El almacenamiento en caché está deshabilitado para processes modulares
- [ ] Todas las tareas deben re-ejecutarse
- [x] El almacenamiento en caché funciona normalmente basándose en los scripts de trabajo generados
- [ ] Solo el archivo principal del workflow se almacena en caché

Más información: [3.2. Ejecutar el workflow](#32-ejecutar-el-workflow)
</quiz>

<quiz>
¿Qué especifica la directiva `container` en una definición de process?
- [ ] El directorio de trabajo para el process
- [ ] La asignación máxima de memoria
- [x] La URI de la imagen del contenedor a usar para ejecutar el process
- [ ] El formato del archivo de salida

Más información: [4.2. Usar un contenedor en un workflow](#42-usar-un-contenedor-en-un-workflow)
</quiz>

<quiz>
En el archivo `.command.run`, ¿qué contiene la función `nxf_launch`?
- [ ] La información de versión de Nextflow
- [ ] Los parámetros del workflow
- [x] El comando `docker run` con montajes de volumen y configuraciones de contenedor
- [ ] Las declaraciones de entrada del process

Más información: [4.2.4. Inspeccionar cómo Nextflow lanzó la tarea en contenedor](#424-inspeccionar-cómo-nextflow-lanzó-la-tarea-en-contenedor)
</quiz>

<quiz>
¿Qué maneja automáticamente Nextflow al ejecutar un process en contenedor? (Seleccione todas las que apliquen)
- [x] Descargar la imagen del contenedor si es necesario
- [x] Montar el directorio de trabajo en el contenedor
- [x] Ejecutar el script del process dentro del contenedor
- [x] Limpiar la instancia del contenedor después de la ejecución

Más información: [4. Usar software en contenedores](#4-usar-software-en-contenedores)
</quiz>
