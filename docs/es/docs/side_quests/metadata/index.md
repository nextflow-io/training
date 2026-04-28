# Metadatos y Meta Maps

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En cualquier análisis científico, rara vez trabajamos solo con los archivos de datos sin procesar.
Cada archivo viene con su propia información adicional: qué es, de dónde proviene y qué lo hace especial.
Esta información extra es lo que llamamos metadatos.

Los metadatos son datos que describen otros datos.
Los metadatos rastrean detalles importantes sobre los archivos y las condiciones experimentales, y ayudan a adaptar los análisis a las características únicas de cada conjunto de datos.

Piénsalo como el catálogo de una biblioteca: mientras que los libros contienen el contenido real (datos sin procesar), las fichas del catálogo proporcionan información esencial sobre cada libro—cuándo fue publicado, quién lo escribió, dónde encontrarlo (metadatos).
En los pipelines de Nextflow, los metadatos se pueden usar para:

- Rastrear información específica de cada archivo a lo largo del workflow
- Configurar procesos según las características de los archivos
- Agrupar archivos relacionados para análisis conjunto

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo manejar metadatos en workflows.
Comenzando con una hoja de datos simple (a menudo llamada samplesheet en bioinformática) que contiene información básica de los archivos, aprenderás a:

- Leer y analizar metadatos de archivos CSV
- Entender por qué la interfaz "meta map + archivo de datos" es una convención ampliamente utilizada
- Agregar nuevos campos de metadatos durante la ejecución del workflow
- Usar metadatos para personalizar el comportamiento de los procesos y organizar las salidas

Estas habilidades te ayudarán a construir pipelines más robustos y flexibles que puedan manejar relaciones complejas entre archivos y requisitos de procesamiento.

### Requisitos previos

Antes de comenzar esta misión secundaria, debes:

- Haber completado el tutorial [Hello Nextflow](../../hello_nextflow/index.md) o un curso equivalente para principiantes.
- Sentirte cómodo/a usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores)

---

## 0. Primeros pasos

#### Abre el codespace de capacitación

Si aún no lo has hecho, asegúrate de abrir el entorno de capacitación como se describe en la [Configuración del entorno](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Muévete al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos de este tutorial.

```bash
cd side-quests/metadata
```

Puedes configurar VSCode para que se enfoque en este directorio:

```bash
code .
```

El editor se abre con el directorio del proyecto en foco.

#### Revisa los materiales

Encontrarás un archivo principal del workflow y un directorio `data` que contiene una hoja de datos y algunos archivos de datos.

??? abstract "Contenido del directorio"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

El workflow en el archivo `main.nf` es un esqueleto que irás expandiendo gradualmente hasta convertirlo en un workflow completamente funcional.

La hoja de datos lista las rutas a los archivos de datos y algunos metadatos asociados, organizados en 3 columnas:

- `id`: autoexplicativo, un ID asignado al archivo
- `character`: un nombre de personaje, que usaremos más adelante para dibujar diferentes criaturas
- `data`: rutas a archivos `.txt` que contienen saludos en diferentes idiomas

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Cada archivo de datos contiene texto de saludo en uno de cinco idiomas (fr: francés, de: alemán, es: español, it: italiano, en: inglés).

Vamos a usar una herramienta llamada [`COWPY`](https://github.com/jeffbuttars/cowpy) para generar arte ASCII de cada personaje pronunciando su saludo grabado.

??? info "¿Qué hace `COWPY`?"

    `COWPY` es una herramienta de línea de comandos que genera arte ASCII para mostrar entradas de texto arbitrarias de una manera divertida.
    Es una implementación en Python de la clásica herramienta [cowsay](https://en.wikipedia.org/wiki/Cowsay) de Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Opcionalmente, puedes seleccionar un personaje (o 'cowacter') para usar en lugar de la vaca predeterminada.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
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

Además, usaremos una herramienta de análisis de idiomas llamada `langid` para identificar qué idioma habla cada personaje y organizar las salidas del pipeline en consecuencia.

#### Revisa la tarea

Tu desafío es escribir un workflow de Nextflow que:

1. **Genere arte ASCII** de cada personaje
2. **Organice** las salidas por familia lingüística (lenguas germánicas vs. románicas)

Esto representa un patrón típico de workflow donde los metadatos específicos de cada archivo impulsan las decisiones de procesamiento; exactamente el tipo de problema que los meta maps resuelven de manera elegante.

#### Lista de verificación de preparación

¿Crees que estás listo/a para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo correctamente
- [ ] Entiendo la tarea

Si puedes marcar todas las casillas, estás listo/a para continuar.

---

## 1. Opciones básicas para cargar y usar metadatos

Abre el archivo del workflow `main.nf` para examinar el esqueleto del workflow que te damos como punto de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

El operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) lee cada fila del archivo como un elemento del canal.
Este es el mismo enfoque que usamos para cargar datos CSV en Hello Nextflow, nuestro curso para principiantes.
Consulta [esta sección](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) si necesitas recordar cómo funciona.

Con `header: true`, la primera fila se trata como encabezados de columna, por lo que cada elemento se convierte en un map de pares clave-valor indexados por nombre de columna.

Ten en cuenta que como aún no estamos ejecutando ningún proceso sobre los datos, los bloques `publish` y `output` son solo esqueletos.

### 1.1. Ejecutar el workflow

Ejecuta el workflow para ver cómo está estructurado el contenido del canal una vez que todo está cargado:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Como puedes ver, el operador ha construido un map de pares clave-valor para cada fila del archivo CSV, con los encabezados de columna como claves para los valores correspondientes.

Cada entrada del map corresponde a una columna en nuestra hoja de datos:

- `id`
- `character`
- `recording`

Esto facilita el acceso a campos específicos de cada fila.
Por ejemplo, podríamos acceder al ID del archivo con `id` o a la ruta del archivo txt con `recording`.

??? info "(Opcional) Más sobre los maps de Groovy"

    En Groovy, el lenguaje de programación sobre el que está construido Nextflow, un map es una estructura de datos de clave-valor similar a los diccionarios en Python, los objetos en JavaScript o los hashes en Ruby.

    Aquí hay un script ejecutable que muestra cómo puedes definir un map y acceder a su contenido en la práctica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crea un map simple
    def my_map = [id:'sampleA', character:'squirrel']

    // Imprime el map completo
    println "map: ${my_map}"

    // Accede a valores individuales usando notación de punto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Aunque no tiene un bloque `workflow` propiamente dicho, Nextflow puede ejecutar esto como si fuera un workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Y esto es lo que puedes esperar ver en la salida:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Seleccionar un campo específico con `map`

Vamos a usar el operador `map` para iterar sobre cada elemento de un canal y seleccionar únicamente el campo `character`, al que podemos acceder por nombre usando notación de punto.

#### 1.2.1. Agregar la operación map

Para acceder a la columna `character`, agrega la operación `map` antes de la operación `.view()` de la siguiente manera:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Esta forma de acceder a un campo específico se explica con más detalle en [esta sección](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) de Hello Nextflow, si necesitas recordar cómo funciona.

#### 1.2.2. Ejecutar el workflow

Ejecuta el workflow para verificar que puedes ver los nombres de personajes extraídos.

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Esto muestra que podemos acceder a los valores de la columna `character` para cada fila.

Ahora hagamos algo con estos datos: usar los campos `character` y `recording` juntos para generar arte ASCII usando `COWPY`.

### 1.3. Emitir sub-canales con `multiMap`

Te proporcionamos un módulo de proceso preescrito llamado `COWPY`, así que primero necesitas examinar los requisitos de entrada del proceso.

Puedes abrir el archivo para ver cómo luce el proceso:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Genera arte ASCII con cowpy
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Como puedes ver, el proceso toma dos entradas separadas: un archivo de grabación y un nombre de personaje.
Es importante notar que tenemos valores para ambos, pero actualmente están agrupados dentro de cada elemento del canal.

Una forma de extraer múltiples campos en canales separados es el operador [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap), que divide un canal en múltiples sub-canales con nombre en una sola operación.

#### 1.3.1. Agregar la operación multiMap

Reemplaza la operación `map` con `multiMap`:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

El bloque `multiMap` define dos sub-canales con nombre (`file` y `character`) a partir de cada fila, a los que podemos acceder como `ch_datasheet.file` y `ch_datasheet.character`.

#### 1.3.2. Llamar a COWPY con los sub-canales

Ahora incluye el proceso `COWPY` y pásale cada sub-canal como argumento separado:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Esto nos permite pasar los dos campos por separado como requiere `COWPY`.

#### 1.3.3. Configurar la publicación de salidas

Finalmente, agrega la salida de `COWPY` al bloque `publish:`:

=== "Después"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Esto nos permitirá ver fácilmente las salidas producidas por el workflow.

#### 1.3.4. Ejecutar el workflow

Ejecuta el workflow para verificar que `COWPY` se ejecuta con las entradas que proporcionamos:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Como puedes ver, `COWPY` se ejecutó en cada archivo usando el personaje correcto para cada uno.

??? abstract "Contenido del directorio de resultados"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Contenido de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Este enfoque funciona, pero tiene una limitación: tuvimos que dividir el canal en dos sub-canales separados.
Si quisiéramos pasar más campos al proceso, tendríamos que dividirlos en más sub-canales.
Eso podría volverse tedioso y desordenado.

Buenas noticias: hay una forma más sencilla de hacer esto.

### 1.4. Agrupar todo como una única entrada al proceso

En lugar de dividir los campos en canales separados, podemos actualizar el proceso para recibir todas las entradas como una única tupla, lo que simplifica la llamada al proceso.

#### 1.4.1. Actualizar el proceso COWPY

Actualiza `COWPY` para aceptar una tupla que corresponda a los tres elementos de cada fila:

=== "Después"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera arte ASCII con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Genera arte ASCII con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
        """
    }
    ```

Ahora el proceso toma una única entrada que contiene todo lo que podríamos querer pasarle.

#### 1.4.2. Usar `map()` para crear la tupla de entrada

Aún necesitamos usar una operación de mapeo para enumerar los elementos que queremos pasar en la tupla al proceso:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Quizás te preguntes por qué no podemos simplemente pasar el map de Groovy completo que proviene de `splitCsv` tal como está.
Es porque necesitamos indicarle a Nextflow explícitamente que el archivo de grabación debe manejarse como una ruta (es decir, debe ser preparado correctamente).
Eso ocurre a nivel de la interfaz de entrada de `COWPY`, donde el elemento `recording` se designa explícitamente como un `path`.

#### 1.4.3. Actualizar la llamada al proceso

Finalmente, reemplazamos las dos entradas separadas en la llamada al proceso con la única tupla que acabamos de crear:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Esto simplifica un poco la llamada al proceso.

#### 1.4.4. Ejecutar el workflow

Ejecuta el workflow para verificar que `COWPY` aún puede procesar los datos correctamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

La salida son los mismos siete archivos `cowpy-*.txt` que antes, ahora producidos con una llamada más simple a `COWPY`.

??? abstract "Contenido del directorio de resultados"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Contenido de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Esto es una ligera mejora sobre el enfoque con `multiMap`.
Pero aún tuvimos que desempaquetar el map de Groovy original para crear la tupla de entrada, y hay un acoplamiento estrecho entre el proceso y la hoja de datos: la definición de entrada de `COWPY` ahora referencia directamente los nombres de columna `id`, `character` y `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Si un colaborador usa una hoja de datos con una estructura diferente —con columnas adicionales o en un orden distinto— este proceso no funcionará sin modificaciones.
Esto hace que el proceso sea frágil, porque su estructura de entrada está ligada a la composición exacta de la hoja de datos.

Para resolver esto, necesitamos una forma de pasar todos los metadatos como un paquete sin codificar su estructura exacta en la interfaz del proceso.

### 1.5. Usar una interfaz de meta map + archivo

La solución es separar dos preocupaciones distintas en el canal: los **metadatos sobre una muestra** y el **archivo de datos** en sí.
Al agrupar todos los metadatos en un único map —el "meta map"— obtenemos una tupla consistente de dos elementos independientemente de cuántas columnas de metadatos contenga la hoja de datos:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Agregar o eliminar columnas de la hoja de datos cambia lo que hay dentro de `meta`, pero la forma de la tupla `[meta, file]` permanece constante.
Los procesos que aceptan esta estructura no necesitan saber ni preocuparse por cuántos campos de metadatos existen.

#### 1.5.1. Reorganizar el contenido de la tupla en un meta map

Reestructuremos la operación `map` para producir una tupla `[meta, file]`:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Se actualizará en el siguiente paso

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Notarás que también agregamos una declaración `view()`, comentamos la llamada a `COWPY` y reemplazamos `COWPY.out` con `channel.empty()` porque la definición de entrada del proceso aún no coincide con la nueva estructura.

#### 1.5.2. Ejecutar el workflow para inspeccionar el contenido reorganizado

Ejecuta el workflow para ver la nueva forma del canal:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Cada elemento del canal es ahora una tupla de dos elementos: el meta map primero, el archivo segundo.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Si más adelante agregamos una columna `language` a la hoja de datos, estará disponible como `meta.language` sin necesidad de cambiar la definición de entrada del proceso.

#### 1.5.3. Actualizar el proceso `COWPY` para usar el meta map

Actualiza `COWPY` para aceptar la estructura de tupla `[meta, file]`:

=== "Después"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera arte ASCII con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera arte ASCII con cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Dentro del bloque script, `meta.character` accede al campo `character` del meta map.
Cualquier campo del meta map es accesible de la misma manera.

#### 1.5.4. Actualizar la llamada al proceso

Restaura la llamada a `COWPY` y conecta su salida para publicación:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Se actualizará en el siguiente paso

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

También hemos restaurado la publicación de salidas.

#### 1.5.5. Ejecutar el workflow

Ejecuta el workflow para verificar que todo funciona:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

El directorio de resultados ahora contiene los archivos de arte ASCII.

??? abstract "Contenido del directorio"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Contenido de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

El proceso ahora recibe todos los metadatos como un paquete a través de `meta`, usa lo que necesita (`meta.character`) e ignora el resto.

Esta es la interfaz estándar utilizada por todos los módulos de [nf-core](https://nf-co.re/).
El patrón `tuple val(meta), path(file)` aparece de manera consistente en toda la biblioteca de módulos de nf-core, razón por la cual los workflows que adoptan esta convención pueden incorporar módulos de nf-core con mínima fricción.

### Conclusión

En esta sección, has aprendido:

- **Cómo leer hojas de datos:** Usar `splitCsv` para analizar archivos CSV con información de encabezado
- **Por qué existe la convención del meta map:** Separar los metadatos de los archivos de datos en tuplas `[meta, file]` mantiene la estructura del canal estable a medida que evoluciona la hoja de datos
- **Cómo usar los campos del meta map dentro de un proceso:** Cualquier campo del meta map es accesible mediante notación de punto en el bloque script

---

## 2. Manipulaciones adicionales de metadatos

Ahora que la interfaz del meta map está en su lugar, podemos enriquecerla a medida que los datos fluyen por el pipeline.

Vamos a usar una herramienta llamada [`langid`](https://github.com/saffsd/langid.py) para identificar el idioma en cada archivo de grabación.
Dado un fragmento de texto, produce una predicción del idioma y una puntuación de probabilidad en `stdout`.

### 2.1. Agregar un paso de identificación de idioma

Te proporcionamos un módulo de proceso preescrito llamado `IDENTIFY_LANGUAGE` que envuelve la herramienta `langid`.

Abre el archivo del módulo para examinar su código:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
// Usa langid para predecir el idioma de cada archivo de entrada
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

La definición de entrada usa la misma estructura `tuple val(meta), path(file)` que acabamos de construir en la Sección 1, por lo que `ch_datasheet` puede alimentar directamente este proceso sin ninguna adaptación.

La salida agrega `stdout` como tercer elemento: esto captura la predicción del idioma que `langid` imprime en la consola.
El comando `sed` elimina la puntuación de probabilidad y el salto de línea final, dejando solo el código de idioma de dos letras.

#### 2.1.1. Agregar una llamada a `IDENTIFY_LANGUAGE`

Incluye el módulo del proceso `IDENTIFY_LANGUAGE` y llámalo en el canal de la hoja de datos:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

La salida principal de este proceso es solo una cadena de texto, por lo que no hay archivos de salida que publicar.
En cambio, usamos `IDENTIFY_LANGUAGE.out.view()` para ver los resultados de la operación.

#### 2.1.2. Ejecutar el workflow

Ejecuta el workflow para producir la identificación de idioma, usando `-resume` para evitar volver a ejecutar las tareas de `COWPY`:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Ahora tenemos una predicción del idioma para cada archivo del conjunto de datos.

Ten en cuenta que la tupla de salida está compuesta por `[meta, file, lang_id]`, lo que significa que el meta map y el archivo se transportan junto con el nuevo resultado.

!!! note "Nota"

    Este patrón de mantener el meta map asociado con los resultados facilita la unión de resultados entre canales más adelante.
    No puedes confiar en el orden de los elementos en los canales para asociar los datos correctamente.
    En cambio, debes usar claves.
    Los meta maps proporcionan una estructura ideal para este propósito.

    Este caso de uso se explora en detalle en la misión secundaria [Splitting & Grouping](../splitting_and_grouping/index.md).

### 2.2. Ampliar los metadatos con las salidas del proceso

La predicción del idioma es en sí misma un metadato sobre los datos del archivo.
En lugar de mantenerla como un elemento separado, incorporémosla de vuelta al meta map.

#### 2.2.1. Crear un nuevo meta map ampliado

Podemos crear un nuevo meta map para reemplazar el original usando el operador `+` de Groovy:

=== "Después"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

El núcleo de esta operación es `#!groovy meta + [lang: lang_id]`.

Ese código esencialmente crea un map temporal con un único par clave-valor que contiene el código de idioma (`[lang: lang_id]`), luego usa el operador `+` de Groovy para combinarlo con el map `meta` original que contiene los metadatos preexistentes, produciendo un nuevo meta map ampliado.

Para una explicación más detallada, consulta el recuadro a continuación.

??? info "Creación del nuevo meta map usando el operador `+`"

    **Primero, necesitas saber que podemos fusionar el contenido de dos maps usando el operador de Groovy `+`.**

    Supongamos que tenemos los siguientes maps:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Podemos fusionarlos así:

    ```groovy
    new_map = map1 + map2
    ```

    El contenido de `new_map` será:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    ¡Genial!

    **¿Pero qué pasa si necesitas agregar un campo que aún no forma parte de un map?**

    Supongamos que comienzas de nuevo desde `map1`, pero la predicción del idioma no está en su propio map (no hay `map2`).
    En cambio, está almacenada en una variable llamada `lang_id`, y sabes que quieres guardar su valor (`'fr'`) con la clave `lang`.

    En realidad puedes hacer lo siguiente:

    ```groovy
    new_map = map1 + [lang: lang_id]
    ```

    Aquí, `[lang: lang_id]` crea un nuevo map sin nombre al vuelo, y `map1 + ` fusiona `map1` con el nuevo map sin nombre, produciendo el mismo contenido de `new_map` que antes.

    ¡Muy ingenioso, ¿verdad?!

    **Ahora transpongamos eso al contexto de una operación `channel.map()` de Nextflow.**

    El código se convierte en:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Esto hace lo siguiente:

    - `#!groovy map1, lang_id ->` toma los dos elementos de la tupla
    - `#!groovy map1 + [lang: lang_id]` crea el nuevo map como se detalla arriba

    La salida es un único map sin nombre con el mismo contenido que `new_map` en nuestro ejemplo anterior.
    Así que hemos transformado efectivamente:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    en:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Con suerte puedes ver que si cambiamos `map1` por `meta`, eso es básicamente todo lo que necesitamos para agregar la predicción del idioma a nuestro meta map en nuestro workflow.

    ¡Excepto por una cosa!

    En el caso de nuestro workflow, **también necesitamos tener en cuenta la presencia del objeto `file` en la tupla**, que está compuesta por `meta, file, lang_id`.

    Entonces el código aquí se convertiría en:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Si te cuesta seguir por qué el `file` parece moverse en la operación `map`, imagina que en lugar de `#!groovy [meta + [lang: lang_id], file]`, esa línea dice `[new_map, file]`.
    Esto debería dejar más claro que simplemente estamos dejando el `file` en su lugar original en segunda posición en la tupla. Solo hemos tomado el valor `new_info` y lo hemos incorporado al map que está en primera posición.

    **¡Y esto nos lleva de vuelta a la estructura de canal `tuple val(meta), path(file)`!**

#### 2.2.2. Ejecutar el workflow

Una vez que estés seguro/a de entender lo que hace el código, ejecuta el workflow para ver si funcionó:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

¡Sí, eso es correcto!
Hemos reorganizado ordenadamente la salida del proceso de `meta, file, lang_id` de modo que `lang_id` ahora es una de las claves en el meta map, y las tuplas del canal encajan nuevamente en el modelo `meta, file`.

!!! tip "Eliminar claves de un meta map"

    Puedes eliminar una clave de un meta map usando el método [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) de Groovy, que devuelve un nuevo map que contiene solo las claves que especifiques:

    ```groovy
    meta.subMap(['id', 'character'])  // devuelve un map con solo 'id' y 'character'
    ```

    Esto es útil cuando un proceso o módulo posterior no necesita todos los campos que se han acumulado en el meta map.

### 2.3. Asignar un grupo lingüístico usando condicionales

Con la predicción del idioma en el meta map, podemos derivar más metadatos a partir de ella.
Los idiomas en nuestro conjunto de datos pertenecen a dos familias: germánicas (inglés, alemán) y románicas (francés, español, italiano).
Agregar un campo `lang_group` hará que esa clasificación esté disponible más adelante en el pipeline.

#### 2.3.1. Agregar una operación `map` con la lógica condicional

Vamos a usar una segunda operación `map` con lógica condicional para asignar la familia lingüística:

```groovy
.map { meta, file ->

    // la lógica condicional que define lang_group va aquí

    [meta + [lang_group: lang_group], file]
}
```

Aquí está la lógica a aplicar:

- Comenzar con `lang_group = 'unknown'` como valor predeterminado.
- Si `meta.lang` es `'de'` o `'en'`, establecer `lang_group` en `'germanic'`.
- De lo contrario, si `meta.lang` está en `['fr', 'es', 'it']`, establecer `lang_group` en `'romance'`.

!!! tip "Consejo"

    Puedes acceder al valor de `lang` dentro de la operación map con `meta.lang`.

Realiza los siguientes cambios al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        ch_languages.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Puntos clave:

- `def lang_group = "unknown"` inicializa la variable con un valor predeterminado seguro.
- La estructura `if / else if` maneja las dos familias lingüísticas; cualquier otro caso permanece como `'unknown'`.
- `#!groovy .set { ch_languages }` le da un nombre al canal resultante para usarlo en el siguiente paso.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Ejecutar el workflow

Ejecuta el workflow para verificar que funciona:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

El meta map ahora contiene cuatro campos: `id`, `character`, `lang` y `lang_group`.
La estructura del canal sigue siendo `[meta, file]`.

### 2.4. Usar metadatos para nombrar y organizar las salidas

Con `lang` y `lang_group` ahora disponibles en el meta map, podemos usarlos para agregar un código de idioma a los nombres de los archivos de salida y organizarlos en subdirectorios por familia lingüística.

Esto requiere tres cambios: actualizar el proceso `COWPY` para renombrar su salida e incluir `meta` en lo que emite, actualizar la llamada a `COWPY` para ejecutarse en `ch_languages`, y actualizar el bloque de salida para especificar la ruta del subdirectorio.

#### 2.4.1. Actualizar el proceso `COWPY`

Renombra el archivo de salida usando el código de idioma del meta map, y agrega `meta` a la salida para que el bloque de salida pueda acceder a `lang_group` para el enrutamiento a subdirectorios:

=== "Después"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Esto muestra cómo podemos aprovechar otros campos de metadatos para personalizar el comportamiento de un proceso, sin necesidad de modificar la definición de entrada en absoluto.

#### 2.4.2. Actualizar la llamada a `COWPY` para ejecutarse en `ch_languages`

Reemplaza `COWPY(ch_datasheet)` con `COWPY(ch_languages)`:

=== "Después"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

También eliminamos la línea `ch_languages.view()` ya que no necesitamos inspeccionar el contenido del canal más.

#### 2.4.3. Actualizar el bloque de salida

Agrega un closure `path` al bloque `output {}` para enrutar cada archivo a su subdirectorio de grupo lingüístico:

=== "Después"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Esto muestra cómo podemos usar metadatos para organizar las salidas con gran flexibilidad.

#### 2.4.4. Ejecutar el pipeline completo

Elimina los resultados anteriores y ejecuta el pipeline completo:

```bash
rm -r results
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

El directorio de resultados ahora está organizado por familia lingüística, con cada archivo nombrado según el idioma detectado:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

El closure `path` en el bloque `output {}` recibe cada tupla `[meta, file]` y devuelve `meta.lang_group` como nombre del subdirectorio.
El nombre del archivo en sí proviene de lo que el proceso produce (`#!groovy "${meta.lang}-${input_file}"`).
Ambas piezas de metadatos (código de idioma y grupo lingüístico) provienen del meta map enriquecido construido en esta sección.

### Conclusión

En esta sección, has aprendido:

- **Cómo ampliar el meta map con salidas de procesos:** Agregar nuevas claves con `#!groovy meta + [clave: valor]` mantiene intacta la estructura del canal `[meta, file]` mientras enriquece los metadatos.
- **Cómo derivar metadatos a partir de metadatos:** La lógica condicional dentro de una operación `map` puede calcular nuevos campos a partir de los existentes.
- **Cómo usar metadatos para organizar las salidas:** El closure `path` en el bloque `output {}` puede leer del meta map para enrutar archivos a subdirectorios.

---

## 3. Consideraciones de robustez

Cuando los valores de los metadatos impulsan el comportamiento de los procesos, los datos faltantes o incompletos pueden causar problemas difíciles de diagnosticar.
Aquí te explicamos qué esperar y cómo manejarlo.

### 3.1. Qué sucede cuando falta un campo de metadatos requerido

El valor `character` es necesario para que el proceso `COWPY` produzca un resultado válido.
El modo de fallo depende de si la columna existe en la hoja de datos pero está vacía, o si está completamente ausente.

#### 3.1.1. La columna existe pero un valor está vacío

Supongamos que una entrada en la hoja de datos tiene el campo `character` en blanco:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La clave `character` se crea para todas las entradas cuando se analiza la hoja de datos, pero `meta.character` para `sampleA` será una cadena vacía.
Cuando Nextflow sustituye `#!groovy ${meta.character}` en el comando, la herramienta `COWPY` recibe un argumento vacío para `-c` y falla:

??? failure "Salida del comando"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

El mensaje de error (`expected one argument`) apunta al indicador `-c` vacío.
Revisar el archivo `.command.sh` en el directorio de trabajo confirma que el comando se ejecutó con un valor vacío.

#### 3.1.2. La columna no existe en la hoja de datos

Si la columna `character` está completamente ausente:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La clave `character` nunca se crea en el meta map.
Cuando el script del proceso evalúa `#!groovy ${meta.character}`, la clave faltante devuelve `null`, y Nextflow literalmente sustituye la cadena `null` en el comando:

??? failure "Salida del comando"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

El `cowpy -c null` en el comando ejecutado es la pista diagnóstica.

### 3.2. Estrategias para manejar metadatos faltantes

Hay dos enfoques complementarios para hacer los workflows más robustos ante metadatos faltantes.

**1. Validación de entrada**

La solución más confiable es validar la hoja de datos antes de que comience cualquier procesamiento, para que los problemas se detecten temprano con un mensaje de error claro en lugar de aparecer como un fallo críptico del proceso a mitad de la ejecución.
La capacitación [Hello nf-core](../../hello_nf-core/05_input_validation.md) cubre cómo agregar validación de entrada usando el plugin nf-schema. <!-- TODO (future) pending a proper Validation side quest -->

**2. Entradas explícitas del proceso para valores requeridos**

Si quieres que la interfaz del proceso comunique que un valor particular es obligatorio, considera extraerlo del meta map como una entrada explícita:

=== "Definición del proceso"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Llamada en el workflow"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Este enfoque hace que `character` sea una parte visible y requerida del contrato del proceso.
Cualquier persona que lea el módulo puede ver inmediatamente que se debe proporcionar un valor de personaje.
Si el campo está ausente, el workflow falla claramente a nivel del canal antes de que el proceso siquiera se ejecute.

Esto destaca un principio de diseño útil:

**Usa el meta map para información opcional o descriptiva; extrae los valores requeridos como entradas explícitas.**

El meta map mantiene las estructuras de canal limpias y estables, pero para los valores que son genuinamente requeridos por un proceso, exponerlos como entradas con nombre mejora la claridad y hace que el módulo sea más fácil de usar correctamente en otros contextos.

### Conclusión

En esta sección, has visto:

- **Cómo se manifiestan los metadatos faltantes:** Un campo vacío produce un argumento vacío; un campo ausente produce `null` sustituido literalmente en el comando.
- **Dos estrategias complementarias:** Validación de entrada para detectar problemas temprano, y entradas explícitas del proceso para comunicar los requisitos claramente.

---

## Resumen

En esta misión secundaria, has explorado cómo trabajar eficazmente con metadatos en los workflows de Nextflow.

El patrón de tupla "meta map + archivo de datos" es una convención fundamental en Nextflow, que ofrece varias ventajas sobre pasar metadatos como valores individuales:

- La estructura del canal permanece estable a medida que evoluciona la hoja de datos
- El comportamiento del proceso puede personalizarse por muestra sin codificar directamente los nombres de los campos
- Los metadatos están disponibles a lo largo del pipeline para nombrar, agrupar y organizar las salidas
- Los módulos escritos para esta interfaz son intercambiables, incluidos los módulos de nf-core

### Patrones clave

1.  **Lectura y estructuración de metadatos:** Analizar una hoja de datos CSV y crear un meta map.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Expansión de metadatos durante el workflow:** Agregar nuevas claves a partir de salidas de procesos o lógica derivada.

    ```groovy
    // A partir de una salida de proceso
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // A partir de lógica condicional
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Uso de metadatos dentro de un proceso:** Acceder a cualquier campo mediante notación de punto en el bloque script.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organización de salidas por valor de metadatos:** Usar el closure `path` en el bloque `output {}`.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Recursos adicionales

- [operador map](https://www.nextflow.io/docs/latest/operator.html#map)
- [operador multiMap](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [calificador de salida stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## ¿Qué sigue?

Regresa al [menú de misiones secundarias](../index.md) o haz clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
