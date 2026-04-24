# Metadatos y meta maps

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
- Crear y manipular meta maps
- Agregar nuevos campos de metadatos durante la ejecución del workflow
- Usar metadatos para personalizar el comportamiento de los procesos

Estas habilidades te ayudarán a construir pipelines más robustos y flexibles que puedan manejar relaciones complejas entre archivos y requisitos de procesamiento.

### Requisitos previos

Antes de comenzar esta misión secundaria, debes:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso equivalente para principiantes.
- Sentirte cómodo/a usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores)

---

## 0. Primeros pasos

#### Abre el codespace de capacitación

Si aún no lo has hecho, asegúrate de abrir el entorno de capacitación como se describe en la [Configuración del entorno](../envsetup/index.md).

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

También te proporcionaremos una herramienta de análisis de idiomas en contenedor llamada `langid`.

#### Revisa la tarea

Tu desafío es escribir un workflow de Nextflow que:

1. **Identifique** el idioma en cada archivo automáticamente
2. **Agrupe** los archivos por familia lingüística (lenguas germánicas vs. románicas)
3. **Personalice** el procesamiento de cada archivo según su idioma y metadatos
4. **Organice** las salidas por grupo lingüístico

Esto representa un patrón típico de workflow donde los metadatos específicos de cada archivo impulsan las decisiones de procesamiento; exactamente el tipo de problema que los meta maps resuelven de manera elegante.

#### Lista de verificación de preparación

¿Crees que estás listo/a para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo correctamente
- [ ] Entiendo la tarea

Si puedes marcar todas las casillas, estás listo/a para continuar.

---

## 1. Cargar metadatos desde una hoja de datos

Abre el archivo del workflow `main.nf` para examinar el esqueleto del workflow que te damos como punto de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Puedes ver que hemos configurado una fábrica de canales básica para cargar la hoja de datos de ejemplo como un archivo, pero eso aún no leerá el contenido del archivo.
Empecemos agregando eso.

### 1.1. Leer el contenido con `splitCsv`

Necesitamos elegir un operador que analice el contenido del archivo de manera apropiada con el mínimo esfuerzo de nuestra parte.
Como nuestra hoja de datos está en formato CSV, este es un trabajo para el operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), que carga cada fila del archivo como un elemento en el canal.

Realiza los siguientes cambios para agregar una operación `splitCsv()` al código de construcción del canal, más una operación `view()` para verificar que el contenido del archivo se está cargando correctamente en el canal.

=== "Después"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Ten en cuenta que estamos usando la opción `header: true` para indicarle a Nextflow que lea la primera fila del archivo CSV como la fila de encabezado.

¿Veamos qué sale de eso?
Ejecuta el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Podemos ver que el operador ha construido un map de pares clave-valor para cada fila del archivo CSV, con los encabezados de columna como claves para los valores correspondientes.

Cada entrada del map corresponde a una columna en nuestra hoja de datos:

- `id`
- `character`
- `recording`

¡Excelente! Esto facilita el acceso a campos específicos de cada archivo.
Por ejemplo, podríamos acceder al ID del archivo con `id` o a la ruta del archivo txt con `recording`.

??? info "(Opcional) Más sobre los maps"

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
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Seleccionar campos específicos con `map`

Supongamos que queremos acceder a la columna `character` de la hoja de datos e imprimirla.
Podemos usar el operador `map` de Nextflow para iterar sobre cada elemento de nuestro canal y seleccionar específicamente la entrada `character` del objeto map.

Realiza las siguientes ediciones al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Ahora ejecuta el workflow nuevamente:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

¡Éxito! Hemos aprovechado la estructura del map derivada de nuestra hoja de datos para acceder a los valores de columnas individuales para cada fila.

Ahora que hemos leído exitosamente la hoja de datos y tenemos acceso a los datos de cada fila, podemos comenzar a implementar la lógica de nuestro pipeline.

### 1.3. Organizar los metadatos en un 'meta map'

En el estado actual del workflow, los archivos de entrada (bajo la clave `recording`) y los metadatos asociados (`id`, `character`) están todos al mismo nivel, como si estuvieran todos en una misma bolsa.
La consecuencia práctica es que cada proceso que consuma este canal necesitaría estar configurado con esta estructura en mente:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Eso está bien siempre que el número de columnas en la hoja de datos no cambie.
Sin embargo, si agregas aunque sea una sola columna a la hoja de datos, la forma del canal ya no coincidirá con lo que el proceso espera, y el workflow producirá errores.
También hace que el proceso sea difícil de compartir con otros que puedan tener datos de entrada ligeramente diferentes, y podrías terminar teniendo que codificar variables directamente en el proceso que no son necesarias para el bloque de script.

Para evitar este problema, necesitamos encontrar una manera de mantener la estructura del canal consistente independientemente de cuántas columnas contenga la hoja de datos.

Podemos hacer eso recopilando todos los metadatos en un elemento dentro de la tupla, al que llamaremos el mapa de metadatos, o más simplemente 'meta map'.

Realiza las siguientes ediciones a la operación `map`:

=== "Después"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Hemos reestructurado los elementos de nuestro canal en una tupla que consta de dos elementos: el meta map y el objeto de archivo correspondiente.

Ejecutemos el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Ahora, cada elemento en el canal contiene primero el meta map y segundo el objeto de archivo correspondiente:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Como resultado, agregar más columnas en la hoja de datos hará que haya más metadatos disponibles en el map `meta`, pero no cambiará la forma del canal.
Esto nos permite escribir procesos que consuman el canal sin tener que codificar directamente los elementos de metadatos en la especificación de entrada:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Esta es una convención ampliamente utilizada para organizar metadatos en los workflows de Nextflow.

### Conclusión

En esta sección, has aprendido:

- **Por qué los metadatos son importantes:** Mantener los metadatos junto con tus datos preserva información importante sobre los archivos a lo largo del workflow.
- **Cómo leer hojas de datos:** Usar `splitCsv` para leer archivos CSV con información de encabezado y transformar filas en datos estructurados
- **Cómo crear un meta map:** Separar los metadatos de los datos de archivo usando la estructura de tupla `[ [id:valor, ...], archivo ]`

---

## 2. Manipulación de metadatos

¡Ahora que tenemos nuestros metadatos cargados, hagamos algo con ellos!

Vamos a usar una herramienta llamada [`langid`](https://github.com/saffsd/langid.py) para identificar el idioma contenido en el archivo de grabación de cada criatura.
La herramienta viene preentrenada en un conjunto de idiomas, y dado un fragmento de texto, producirá una predicción del idioma y una puntuación de probabilidad asociada, ambas en `stdout`.

### 2.1. Importar el proceso y examinar el código

Te proporcionamos un módulo de proceso preescrito llamado `IDENTIFY_LANGUAGE` que envuelve la herramienta `langid`, por lo que solo necesitas agregar una declaración include antes del bloque del workflow.

Realiza la siguiente edición al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Puedes abrir el archivo del módulo para examinar su código:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

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

Como puedes ver, la definición de entrada usa la misma estructura `tuple val(meta), path(file)` que acabamos de aplicar a nuestro canal de entrada.

La definición de salida está estructurada como una tupla con una estructura similar a la de la entrada, excepto que también contiene `stdout` como tercer elemento.
Este patrón `tuple val(meta), path(file), <salida>` mantiene los metadatos asociados tanto con los datos de entrada como con las salidas a medida que fluyen por el pipeline.

Ten en cuenta que estamos usando el calificador de salida [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) de Nextflow porque la herramienta imprime su salida directamente en la consola en lugar de escribir un archivo; y usamos `sed` en la línea de comandos para eliminar la puntuación de probabilidad, limpiar la cadena eliminando los caracteres de nueva línea y devolver solo la predicción del idioma.

### 2.2. Agregar una llamada a `IDENTIFY_LANGUAGE`

Ahora que el proceso está disponible para el workflow, podemos agregar una llamada al proceso `IDENTIFY_LANGUAGE` para ejecutarlo en el canal de datos.

Realiza las siguientes ediciones al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Ten en cuenta que hemos eliminado la operación `.view()` original en la construcción del canal.

Ahora podemos ejecutar el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

¡Excelente! Ahora tenemos una predicción del idioma que habla cada personaje.

Y como se señaló anteriormente, también hemos incluido el archivo de entrada y el meta map en la salida, lo que significa que ambos permanecen asociados con la nueva información que acabamos de producir.
Esto resultará útil en el siguiente paso.

!!! note "Nota"

    De manera más general, este patrón de mantener el meta map asociado con los resultados facilita la asociación de resultados relacionados que comparten los mismos identificadores.

    Como ya habrás aprendido, no puedes confiar en el orden de los elementos en los canales para hacer coincidir resultados entre ellos.
    En cambio, debes usar claves para asociar los datos correctamente, y los meta maps proporcionan una estructura ideal para este propósito.

    Exploramos este caso de uso en detalle en la misión secundaria [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Ampliar los metadatos con las salidas del proceso

Dado que los resultados que acabamos de producir son en sí mismos una forma de metadatos sobre el contenido de los archivos, sería útil agregarlos a nuestro meta map.

Sin embargo, no queremos modificar el meta map existente en su lugar.
Desde un punto de vista técnico, es _posible_ hacerlo, pero no es seguro.

Por eso, en su lugar, crearemos un nuevo meta map que contenga el contenido del meta map existente más un nuevo par clave-valor `lang: lang_id` que almacene la nueva información, usando el operador `+` (una característica de Groovy).
Y lo combinaremos con una operación [`map`](https://www.nextflow.io/docs/latest/operator.html#map) para reemplazar el map antiguo con el nuevo.

Aquí están las ediciones que necesitas hacer al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Si aún no estás familiarizado/a con el operador `+`, o si esto parece confuso, tómate unos minutos para revisar la explicación detallada a continuación.

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
    new_map = [map1 + [lang: lang_id]]
    ```

    Aquí, `[lang: new_info]` crea un nuevo map sin nombre al vuelo, y `map1 + ` fusiona `map1` con el nuevo map sin nombre, produciendo el mismo contenido de `new_map` que antes.

    ¡Muy ingenioso, ¿verdad?!

    **Ahora transpongamos eso al contexto de una operación `channel.map()` de Nextflow.**

    El código se convierte en:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Esto hace lo siguiente:

    - `map1, lang_id ->` toma los dos elementos de la tupla
    - `[map1 + [lang: lang_id]]` crea el nuevo map como se detalla arriba

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

    Si te cuesta seguir por qué el `file` parece moverse en la operación `map`, imagina que en lugar de `[meta + [lang: lang_id], file]`, esa línea dice `[new_map, file]`.
    Esto debería dejar más claro que simplemente estamos dejando el `file` en su lugar original en segunda posición en la tupla. Solo hemos tomado el valor `new_info` y lo hemos incorporado al map que está en primera posición.

    **¡Y esto nos lleva de vuelta a la estructura de canal `tuple val(meta), path(file)`!**

Una vez que estés seguro/a de entender lo que hace este código, ejecuta el workflow para ver si funcionó:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

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

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Asignar un grupo lingüístico usando condicionales

Ahora que tenemos nuestras predicciones de idioma, usemos la información para asignar nuevas agrupaciones.

En nuestros datos de ejemplo, los idiomas que usan nuestros personajes se pueden agrupar en lenguas germánicas (inglés, alemán) y lenguas románicas (francés, español, italiano).
Podría ser útil tener esa clasificación disponible en algún punto posterior del pipeline, así que agreguemos esa información en el meta map.

Y, buenas noticias, ¡este es otro caso que se presta perfectamente para usar el operador `map`!

Específicamente, vamos a definir una variable llamada `lang_group`, usar algo de lógica condicional simple para determinar qué valor asignar al `lang_group` para cada pieza de datos.

La sintaxis general se verá así:

```groovy
.map { meta, file ->

    // la lógica condicional que define lang_group va aquí

    [meta + [lang_group: lang_group], file]
}
```

Puedes ver que esto es muy similar a la operación de fusión de maps al vuelo que usamos en el paso anterior.
Solo necesitamos escribir las declaraciones condicionales.

Aquí está la lógica condicional que queremos aplicar:

- Definir una variable llamada `lang_group` con valor predeterminado `'unknown'`.
- Si `lang` es alemán (`'de'`) o inglés (`'en'`), cambiar `lang_group` a `germanic`.
- De lo contrario, si `lang` está incluido en una lista que contiene francés (`'fr'`), español (`'es'`) e italiano (`'it'`), cambiar `lang_group` a `romance`.

Intenta escribirlo tú mismo/a si ya sabes cómo escribir declaraciones condicionales en Nextflow.

!!! tip "Consejo"

    Puedes acceder al valor de `lang` dentro de la operación map con `meta.lang`.

Deberías terminar haciendo los siguientes cambios al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
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
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Ejecuta langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Aquí están los puntos clave:

- Usamos `def lang_group = "unknown"` para crear la variable `lang_group` con valor predeterminado establecido en `unknown`.
- Usamos una estructura `if {} else if {}` para la lógica condicional, con pruebas alternativas `.equals()` para los dos idiomas germánicos, y una prueba de existencia en una lista para los tres idiomas románicos.
- Usamos la operación de fusión `meta + [lang_group:lang_group]` como antes para generar el meta map actualizado.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Una vez que todo eso tenga sentido, ejecuta el workflow nuevamente para ver el resultado:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Como puedes ver, los elementos del canal mantienen su estructura `[meta, file]`, pero el meta map ahora incluye esta nueva clasificación.

### Conclusión

En esta sección, has aprendido a:

- **Aplicar metadatos de entrada a canales de salida**: Copiar metadatos de esta manera nos permite asociar resultados más adelante basándonos en el contenido de los metadatos.
- **Crear claves personalizadas**: Creaste dos nuevas claves en tu meta map, fusionándolas con `meta + [nueva_clave:valor]` en el meta map existente. Una basada en un valor calculado de un proceso, y otra basada en una condición que estableciste en el operador `map`.

Estas habilidades te permiten asociar metadatos nuevos y existentes con archivos a medida que avanzas por tu pipeline.
Incluso si no estás usando metadatos como parte de un proceso, mantener el meta map asociado con los datos de esta manera facilita mantener toda la información relevante junta.

---

## 3. Usar información del meta map en un proceso

Ahora que sabes cómo crear y actualizar el meta map, podemos llegar a la parte realmente divertida: usar los metadatos en un proceso.

Más específicamente, vamos a agregar un segundo paso a nuestro workflow para dibujar cada animal como arte ASCII y hacer que diga el texto grabado en un globo de diálogo.
Vamos a hacer esto usando una herramienta llamada [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "¿Qué hace `cowpy`?"

    `cowpy` es una herramienta de línea de comandos que genera arte ASCII para mostrar entradas de texto arbitrarias de una manera divertida.
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

Si trabajaste en el curso Hello Nextflow, ya has visto esta herramienta en acción.
Si no, no te preocupes; cubriremos todo lo que necesitas saber a medida que avancemos.

### 3.1. Importar el proceso y examinar el código

Te proporcionamos un módulo de proceso preescrito llamado `COWPY` que envuelve la herramienta `cowpy`, por lo que solo necesitas agregar una declaración include antes del bloque del workflow.

Realiza la siguiente edición al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Puedes abrir el archivo del módulo para examinar su código:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Genera arte ASCII con cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

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

Como puedes ver, este proceso está actualmente diseñado para tomar un archivo de entrada (que contiene el texto a mostrar) y un valor que especifica el personaje que debe dibujarse en arte ASCII, generalmente proporcionado a nivel del workflow mediante un parámetro de línea de comandos.

### 3.2. Pasar un campo del meta map como entrada

Cuando usamos la herramienta `cowpy` en el curso Hello Nextflow, usamos un parámetro de línea de comandos para determinar qué personaje usar para dibujar la imagen final.
Eso tenía sentido, porque solo generábamos una imagen por ejecución del pipeline.

Sin embargo, en este tutorial, queremos generar una imagen apropiada para cada sujeto que estamos procesando, por lo que usar un parámetro de línea de comandos sería demasiado limitante.

Buenas noticias: tenemos una columna `character` en nuestra hoja de datos y, por lo tanto, en nuestro meta map.
Usemos eso para establecer el personaje que el proceso debe usar para cada entrada.

Para ello, necesitaremos hacer tres cosas:

1. Dar un nombre al canal de salida que proviene del proceso anterior para poder operar sobre él de manera más conveniente.
2. Determinar cómo acceder a la información de interés
3. Agregar una llamada al segundo proceso y proporcionar la información de manera apropiada.

¡Comencemos!

#### 3.2.1. Nombrar el canal de salida anterior

Aplicamos las manipulaciones anteriores directamente en el canal de salida del primer proceso, `IDENTIFY_LANGUAGE.out`.
Para poder pasar el contenido de ese canal al siguiente proceso (y hacerlo de una manera clara y fácil de leer) queremos darle su propio nombre, `ch_languages`.

Podemos hacer eso usando el operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

En el workflow principal, reemplaza el operador `.view()` con `.set { ch_languages }`, y agrega una línea para verificar que podemos referirnos al canal por nombre.

=== "Después"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
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

        // Temporal: inspeccionar ch_languages
        ch_languages.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
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
            .view()
    ```

Ejecutemos esto:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Esto confirma que ahora podemos referirnos al canal por nombre.

#### 3.2.2. Acceder al archivo y a los metadatos del personaje

Sabemos por el código del módulo que el proceso `COWPY` espera recibir un archivo de texto y un valor `character`.
Para escribir la llamada al proceso `COWPY`, solo necesitamos saber cómo extraer el objeto de archivo correspondiente y los metadatos de cada elemento en el canal.

Como suele ser el caso, la forma más sencilla de hacerlo es usar una operación `map`.

Nuestro canal contiene tuplas estructuradas como `[meta, file]`, por lo que podemos acceder al objeto `file` directamente, y podemos acceder al valor `character` almacenado dentro del meta map refiriéndonos a él como `meta.character`.

En el workflow principal, realiza los siguientes cambios de código:

=== "Después"

    ```groovy title="main.nf" linenums="34"
        // Temporal: acceder al archivo y al personaje
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporal: inspeccionar ch_languages
        ch_languages.view()
    ```

Ten en cuenta que estamos usando closures (como `{ file -> "File: " + file }`) para hacer que la salida de las operaciones `.view` sea más legible.

Ejecutemos esto:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Las rutas de archivos y los valores de personajes pueden aparecer en un orden diferente en tu salida._

Esto confirma que podemos acceder al archivo y al personaje para cada elemento en el canal.

#### 3.2.3. Llamar al proceso `COWPY`

Ahora juntemos todo y llamemos al proceso `COWPY` en el canal `ch_languages`.

En el workflow principal, realiza los siguientes cambios de código:

=== "Después"

    ```groovy title="main.nf" linenums="34"
        // Ejecuta cowpy para generar arte ASCII
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporal: acceder al archivo y al personaje
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Verás que simplemente copiamos las dos operaciones map (sin las declaraciones `.view()`) como entradas a la llamada del proceso.
¡Solo asegúrate de no olvidar la coma entre ellas!

Es un poco torpe, pero veremos cómo mejorarlo en la siguiente sección.

Ejecutemos esto:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Si miras en el directorio de resultados, deberías ver los archivos individuales que contienen el arte ASCII de cada saludo pronunciado por el personaje correspondiente.

??? abstract "Contenido del directorio y ejemplo de archivo"

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

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Esto muestra que pudimos usar la información en el meta map para parametrizar el comando en el segundo paso del pipeline.

Sin embargo, como se señaló anteriormente, parte del código involucrado fue un poco torpe, ya que tuvimos que desempaquetar los metadatos mientras aún estábamos en el contexto del cuerpo del workflow.
Ese enfoque funciona bien para usar un pequeño número de campos del meta map, pero escalaría mal si quisiéramos usar muchos más.

Existe otro operador llamado `multiMap()` que nos permite agilizar esto un poco, pero incluso así no es ideal.

??? info "(Opcional) Versión alternativa con `multiMap()`"

    En caso de que te lo estés preguntando, no podíamos simplemente escribir una única operación `map()` que produzca tanto el `file` como el `character`, porque eso los devolvería como una tupla.
    Tuvimos que escribir dos operaciones `map()` separadas para pasar los elementos `file` y `character` al proceso por separado.

    Técnicamente hay otra manera de hacer esto a través de una única operación de mapeo, usando el operador `multiMap()`, que es capaz de emitir múltiples canales.
    Por ejemplo, podrías reemplazar la llamada a `COWPY` anterior con el siguiente código:

    === "Después"

        ```groovy title="main.nf" linenums="34"
            // Ejecuta cowpy para generar arte ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Antes"

        ```groovy title="main.nf" linenums="34"
            // Ejecuta cowpy para generar arte ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Esto produce exactamente el mismo resultado.

En cualquier caso, es incómodo tener que hacer algo de desempaquetado a nivel del workflow.

Sería mejor si pudiéramos pasar el meta map completo al proceso y seleccionar lo que necesitamos una vez allí.

### 3.3. Pasar y usar el meta map completo

El objetivo del meta map es, después de todo, pasar todos los metadatos juntos como un paquete.
La única razón por la que no pudimos hacer eso anteriormente es que el proceso no está configurado para aceptar un meta map.
Pero como controlamos el código del proceso, podemos cambiarlo.

Modifiquemos el proceso `COWPY` para aceptar la estructura de tupla `[meta, file]` que usamos en el primer proceso para poder agilizar el workflow.

Para ello, necesitaremos hacer tres cosas:

1. Modificar las definiciones de entrada del módulo del proceso `COWPY`
2. Actualizar el comando del proceso para usar el meta map
3. Actualizar la llamada al proceso en el cuerpo del workflow

¿Listo/a? ¡Vamos!

#### 3.3.1. Modificar la entrada del módulo `COWPY`

Realiza las siguientes ediciones al archivo del módulo `cowpy.nf`:

=== "Después"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Antes"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Esto nos permite usar la estructura de tupla `[meta, file]` que cubrimos anteriormente en el tutorial.

Ten en cuenta que no actualizamos la definición de salida del proceso para producir el meta map, con el fin de mantener el tutorial simplificado, pero siéntete libre de hacerlo tú mismo/a como ejercicio siguiendo el modelo del proceso `IDENTIFY_LANGUAGE`.

#### 3.3.2. Actualizar el comando para usar el campo del meta map

El meta map completo ahora está disponible dentro del proceso, por lo que podemos referirnos a la información que contiene directamente desde dentro del bloque de comandos.

Realiza las siguientes ediciones al archivo del módulo `cowpy.nf`:

=== "Después"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Antes"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Hemos reemplazado la referencia al valor `character` que anteriormente se pasaba como entrada independiente con el valor almacenado en el meta map, al que nos referimos usando `meta.character`.

Ahora actualicemos la llamada al proceso en consecuencia.

#### 3.3.3. Actualizar la llamada al proceso y ejecutarlo

El proceso ahora espera que su entrada use la estructura de tupla `[meta, file]`, que es lo que produce el proceso anterior, por lo que simplemente podemos pasar el canal `ch_languages` completo al proceso `COWPY`.

Realiza las siguientes ediciones al workflow principal:

=== "Después"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Ejecuta cowpy para generar arte ASCII
    COWPY(ch_languages)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Ejecuta cowpy para generar arte ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

¡Eso simplifica significativamente la llamada!

Eliminemos los resultados de la ejecución anterior y ejecutémoslo:

```bash
rm -r results
nextflow run main.nf
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Si miras en el directorio de resultados, deberías ver las mismas salidas que antes, _es decir_, archivos individuales que contienen el arte ASCII de cada saludo pronunciado por el personaje correspondiente.

??? abstract "Contenido del directorio"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Así que esto produce los mismos resultados que antes con un código más simple.

Por supuesto, esto asume que puedes modificar el código del proceso.
En algunos casos, es posible que tengas que depender de procesos existentes que no tienes libertad de modificar, lo que limita tus opciones.
La buena noticia, si planeas usar módulos del proyecto [nf-core](https://nf-co.re/), es que los módulos de nf-core están todos configurados para usar la estructura de tupla `[meta, file]` como estándar.

### 3.4. Solución de problemas con entradas requeridas faltantes

El valor `character` es necesario para que el proceso `COWPY` se ejecute correctamente.
Si no establecemos un valor predeterminado en un archivo de configuración, DEBEMOS proporcionar un valor en la hoja de datos.

**¿Qué sucede si no lo hacemos?**
Depende de lo que contenga la hoja de datos de entrada y qué versión del workflow estemos ejecutando.

#### 3.4.1. La columna character existe pero está vacía

Supongamos que eliminamos el valor del personaje para una de las entradas en nuestra hoja de datos para simular un error de recopilación de datos:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Para cualquiera de las versiones del workflow que hemos usado anteriormente, la clave `character` se creará para todas las entradas cuando se lea la hoja de datos, pero para `sampleA` el valor será una cadena vacía.

Esto causará un error.

??? failure "Salida del comando"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

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

Cuando Nextflow ejecuta la línea de comandos `cowpy` para esa muestra, `${meta.character}` se rellena con una cadena vacía en la línea de comandos de `cowpy`, por lo que la herramienta `cowpy` lanza un error indicando que no se proporcionó ningún valor para el argumento `-c`.

#### 3.4.2. La columna character no existe en la hoja de datos

Ahora supongamos que eliminamos la columna `character` completamente de nuestra hoja de datos:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

En este caso, la clave `character` no se creará en absoluto cuando se lea la hoja de datos.

##### 3.4.2.1. Valor accedido a nivel del workflow

Si estamos usando la versión del código que escribimos en la sección 3.2, Nextflow intentará acceder a la clave `character` en el meta map ANTES de llamar al proceso `COWPY`.

No encontrará ningún elemento que coincida con la instrucción, por lo que no ejecutará `COWPY` en absoluto.

??? success "Salida del comando"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

¡Desde la perspectiva de Nextflow, este workflow se ejecutó correctamente!
Sin embargo, no se producirá ninguna de las salidas que queremos.

##### 3.4.2.2. Valor accedido a nivel del proceso

Si estamos usando la versión de la sección 3.3, Nextflow pasará el meta map completo al proceso `COWPY` e intentará ejecutar el comando.

Esto causará un error, pero diferente al del primer caso.

??? failure "Salida del comando"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

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

Esto ocurre porque `meta.character` no existe, por lo que nuestro intento de acceder a él devuelve `null`. Como resultado, Nextflow literalmente inserta `null` en la línea de comandos, que por supuesto no es reconocido por la herramienta `cowpy`.

#### 3.4.3. Soluciones

Además de proporcionar un valor predeterminado como parte de la configuración del workflow, hay dos cosas que podemos hacer para manejar esto de manera más robusta:

1. Implementar validación de entrada en tu workflow para asegurarte de que la hoja de datos contenga toda la información requerida. Puedes encontrar una [introducción a la validación de entrada](../hello_nf-core/05_input_validation.md) en el curso de capacitación Hello nf-core. <!-- TODO (future) pending a proper Validation side quest -->

2. Si quieres asegurarte de que cualquier persona que use tu módulo de proceso pueda identificar inmediatamente las entradas requeridas, también puedes hacer que la propiedad de metadatos requerida sea una entrada explícita.

Aquí hay un ejemplo de cómo funcionaría eso.

Primero, a nivel del proceso, actualiza la definición de entrada de la siguiente manera:

=== "Después"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Antes"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Luego, a nivel del workflow, usa una operación de mapeo para extraer la propiedad `character` de los metadatos y convertirla en un componente explícito de la tupla de entrada:

=== "Después"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Este enfoque tiene la ventaja de mostrar explícitamente que `character` es requerido, y hace que el proceso sea más fácil de reutilizar en otros contextos.

Esto destaca un principio de diseño importante:

**Usa el meta map para información opcional y descriptiva, pero extrae los valores requeridos como entradas explícitas.**

El meta map es excelente para mantener las estructuras de canal limpias y evitar estructuras de canal arbitrarias, pero para los elementos obligatorios que se referencian directamente en un proceso, extraerlos como entradas explícitas crea un código más robusto y mantenible.

### Conclusión

En esta sección, has aprendido a utilizar metadatos para personalizar la ejecución de un proceso, accediendo a ellos ya sea a nivel del workflow o a nivel del proceso.

---

## Ejercicio suplementario

Si deseas practicar el uso de información del meta map desde dentro de un proceso, intenta usar otras piezas de información del meta map como `lang` y `lang_group` para personalizar cómo se nombran y/o organizan las salidas.

Por ejemplo, intenta modificar el código para producir este resultado:

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

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Resumen

En esta misión secundaria, has explorado cómo trabajar eficazmente con metadatos en los workflows de Nextflow.

Este patrón de mantener los metadatos explícitos y adjuntos a los datos es una práctica recomendada fundamental en Nextflow, que ofrece varias ventajas sobre la codificación directa de información de archivos:

- Los metadatos de los archivos permanecen asociados con los archivos a lo largo del workflow
- El comportamiento del proceso puede personalizarse por archivo
- La organización de las salidas puede reflejar los metadatos de los archivos
- La información de los archivos puede expandirse durante la ejecución del pipeline

Aplicar este patrón en tu propio trabajo te permitirá construir workflows de bioinformática robustos y mantenibles.

### Patrones clave

1.  **Lectura y estructuración de metadatos:** Leer archivos CSV y crear mapas de metadatos organizados que permanezcan asociados con tus archivos de datos.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Expansión de metadatos durante el workflow:** Agregar nueva información a tus metadatos a medida que avanza tu pipeline, añadiendo salidas de procesos y derivando valores mediante lógica condicional.

    - Agregar nuevas claves basadas en la salida del proceso

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Agregar nuevas claves usando una cláusula condicional

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Personalización del comportamiento del proceso:** Usar metadatos dentro del proceso.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Recursos adicionales

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## ¿Qué sigue?

Regresa al [menú de misiones secundarias](../) o haz clic en el botón en la parte inferior derecha de la página para continuar con el siguiente tema de la lista.
