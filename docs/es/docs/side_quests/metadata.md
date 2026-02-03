# Metadatos y mapas de metadatos

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En cualquier análisis científico, rara vez trabajamos solo con los archivos de datos en bruto.
Cada archivo viene con su propia información adicional: qué es, de dónde proviene y qué lo hace especial.
Esta información extra es lo que llamamos metadatos.

Los metadatos son datos que describen otros datos.
Los metadatos rastrean detalles importantes sobre archivos y condiciones experimentales, y ayudan a adaptar los análisis a las características únicas de cada conjunto de datos.

Piense en ello como un catálogo de biblioteca: mientras que los libros contienen el contenido real (datos en bruto), las fichas del catálogo proporcionan información esencial sobre cada libro—cuándo fue publicado, quién lo escribió, dónde encontrarlo (metadatos).
En los pipelines de Nextflow, los metadatos se pueden usar para:

- Rastrear información específica de archivos a lo largo del workflow
- Configurar procesos basándose en las características de los archivos
- Agrupar archivos relacionados para análisis conjunto

### Objetivos de aprendizaje

En esta misión secundaria, exploraremos cómo manejar metadatos en workflows.
Comenzando con una hoja de datos simple (a menudo llamada samplesheet en bioinformática) que contiene información básica de archivos, aprenderá a:

- Leer y procesar metadatos de archivos desde archivos CSV
- Crear y manipular mapas de metadatos
- Agregar nuevos campos de metadatos durante la ejecución del workflow
- Usar metadatos para personalizar el comportamiento de los procesos

Estas habilidades lo ayudarán a construir pipelines más robustos y flexibles que puedan manejar relaciones complejas de archivos y requisitos de procesamiento.

### Requisitos previos

Antes de embarcarse en esta misión secundaria, debería:

- Haber completado el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curso para principiantes equivalente.
- Sentirse cómodo usando conceptos y mecanismos básicos de Nextflow (procesos, canales, operadores)

---

## 0. Primeros pasos

#### Abrir el codespace de entrenamiento

Si aún no lo ha hecho, asegúrese de abrir el entorno de entrenamiento como se describe en [Configuración del Entorno](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moverse al directorio del proyecto

Vamos a movernos al directorio donde se encuentran los archivos para este tutorial.

```bash
cd side-quests/metadata
```

Puede configurar VSCode para enfocarse en este directorio:

```bash
code .
```

#### Revisar los materiales

Encontrará un archivo de workflow principal y un directorio `data` que contiene una hoja de datos y un puñado de archivos de datos.

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

El workflow en el archivo `main.nf` es un esqueleto que expandirá gradualmente en un workflow completamente funcional.

La hoja de datos enumera las rutas a los archivos de datos y algunos metadatos asociados, organizados en 3 columnas:

- `id`: autoexplicativo, un ID asignado al archivo
- `character`: un nombre de personaje, que usaremos más tarde para dibujar diferentes criaturas
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

También le proporcionaremos una herramienta de análisis de lenguaje en contenedor llamada `langid`.

#### Revisar la asignación

Su desafío es escribir un workflow de Nextflow que:

1. **Identifique** el idioma en cada archivo automáticamente
2. **Agrupe** archivos por familia de idiomas (idiomas germánicos vs románicos)
3. **Personalice** el procesamiento para cada archivo basándose en su idioma y metadatos
4. **Organice** las salidas por grupo de idiomas

Esto representa un patrón típico de workflow donde los metadatos específicos de archivos impulsan las decisiones de procesamiento; exactamente el tipo de problema que los mapas de metadatos resuelven elegantemente.

#### Lista de verificación de preparación

¿Cree que está listo para comenzar?

- [ ] Entiendo el objetivo de este curso y sus requisitos previos
- [ ] Mi codespace está en funcionamiento
- [ ] He configurado mi directorio de trabajo apropiadamente
- [ ] Entiendo la asignación

Si puede marcar todas las casillas, está listo para comenzar.

---

## 1. Cargar metadatos desde una hoja de datos

Abra el archivo de workflow `main.nf` para examinar el esqueleto de workflow que le proporcionamos como punto de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Puede ver que hemos configurado un factory de canal básico para cargar la hoja de datos de ejemplo como un archivo, pero eso aún no leerá el contenido del archivo.
Comencemos agregando eso.

### 1.1. Leer el contenido con `splitCsv`

Necesitamos elegir un operador que procese el contenido del archivo apropiadamente con el mínimo esfuerzo de nuestra parte.
Como nuestra hoja de datos está en formato CSV, este es un trabajo para el operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), que carga cada fila en el archivo como un elemento en el canal.

Realice los siguientes cambios para agregar una operación `splitCsv()` al código de construcción del canal, más una operación `view()` para verificar que el contenido del archivo se esté cargando en el canal correctamente.

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

Tenga en cuenta que estamos usando la opción `header: true` para indicarle a Nextflow que lea la primera fila del archivo CSV como la fila de encabezado.

Veamos qué sale de eso, ¿de acuerdo?
Ejecute el workflow:

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

Podemos ver que el operador ha construido un mapa de pares clave-valor para cada fila en el archivo CSV, con los encabezados de columna como claves para los valores correspondientes.

Cada entrada de mapa corresponde a una columna en nuestra hoja de datos:

- `id`
- `character`
- `recording`

¡Esto es genial! Facilita el acceso a campos específicos de cada archivo.
Por ejemplo, podríamos acceder al ID del archivo con `id` o a la ruta del archivo txt con `recording`.

??? info "(Opcional) Más sobre mapas"

    En Groovy, el lenguaje de programación sobre el que se construye Nextflow, un mapa es una estructura de datos clave-valor similar a los diccionarios en Python, objetos en JavaScript o hashes en Ruby.

    Aquí hay un script ejecutable que muestra cómo puede definir un mapa y acceder a su contenido en la práctica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crear un mapa simple
    def my_map = [id:'sampleA', character:'squirrel']

    // Imprimir todo el mapa
    println "map: ${my_map}"

    // Acceder a valores individuales usando notación de punto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Aunque no tiene un bloque `workflow` apropiado, Nextflow puede ejecutar esto como si fuera un workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Y aquí está lo que puede esperar ver en la salida:

    ```console title="Salida"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Seleccionar campos específicos con `map`

Digamos que queremos acceder a la columna `character` de la hoja de datos e imprimirla.
Podemos usar el operador `map` de Nextflow para iterar sobre cada elemento en nuestro canal y específicamente seleccionar la entrada `character` del objeto mapa.

Realice las siguientes ediciones al workflow:

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

Ahora ejecute el workflow nuevamente:

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

¡Éxito! Hemos aprovechado la estructura de mapa derivada de nuestra hoja de datos para acceder a los valores de columnas individuales para cada fila.

Ahora que hemos leído exitosamente la hoja de datos y tenemos acceso a los datos en cada fila, podemos comenzar a implementar nuestra lógica de pipeline.

### 1.3. Organizar los metadatos en un 'mapa de metadatos'

En el estado actual del workflow, los archivos de entrada (bajo la clave `recording`) y los metadatos asociados (`id`, `character`) están todos al mismo nivel, como si estuvieran todos en una gran bolsa.
La consecuencia práctica es que cada proceso que consuma este canal necesitaría configurarse con esta estructura en mente:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Eso está bien siempre que el número de columnas en la hoja de datos no cambie.
Sin embargo, si agrega incluso solo una columna a la hoja de datos, la forma del canal ya no coincidirá con lo que el proceso espera, y el workflow producirá errores.
También hace que el proceso sea difícil de compartir con otros que podrían tener datos de entrada ligeramente diferentes, y podría terminar teniendo que codificar variables en el proceso que no son necesarias para el bloque de script.

Para evitar este problema, necesitamos encontrar una forma de mantener la estructura del canal consistente independientemente de cuántas columnas contenga esa hoja de datos.

Podemos hacer eso recopilando todos los metadatos en un elemento dentro de la tupla, que llamaremos el mapa de metadatos, o más simplemente 'mapa de metadatos'.

Realice las siguientes ediciones a la operación `map`:

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

Hemos reestructurado nuestros elementos de canal en una tupla que consta de dos elementos, el mapa de metadatos y el objeto archivo correspondiente.

Ejecutemos el workflow:

```bash
nextflow run main.nf
```

??? success "Salida del comando"

    ```console title="Ver mapa de metadatos"
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

Ahora, cada elemento en el canal contiene primero el mapa de metadatos y segundo el objeto archivo correspondiente:

```console title="Ejemplo de estructura de salida"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Como resultado, agregar más columnas en la hoja de datos hará que más metadatos estén disponibles en el mapa `meta`, pero no cambiará la forma del canal.
Esto nos permite escribir procesos que consuman el canal sin tener que codificar los elementos de metadatos en la especificación de entrada:

```groovy title="Ejemplo de sintaxis"
    input:
    tuple val(meta), file(recording)
```

Esta es una convención ampliamente utilizada para organizar metadatos en workflows de Nextflow.

### Conclusión

En esta sección, ha aprendido:

- **Por qué los metadatos son importantes:** Mantener los metadatos con sus datos preserva información importante de archivos a lo largo del workflow.
- **Cómo leer hojas de datos:** Usar `splitCsv` para leer archivos CSV con información de encabezado y transformar filas en datos estructurados
- **Cómo crear un mapa de metadatos:** Separar metadatos de datos de archivos usando la estructura de tupla `[ [id:value, ...], file ]`

---

## 2. Manipulando metadatos

Ahora que tenemos nuestros metadatos cargados, ¡hagamos algo con ellos!

Vamos a usar una herramienta llamada [`langid`](https://github.com/saffsd/langid.py) para identificar el idioma contenido en cada archivo de grabación de la criatura.
La herramienta viene pre-entrenada con un conjunto de idiomas, y dado un fragmento de texto, producirá una predicción de idioma y una puntuación de probabilidad asociada, ambas a `stdout`.

### 2.1. Importar el proceso y examinar el código

Le proporcionamos un módulo de proceso pre-escrito llamado `IDENTIFY_LANGUAGE` que envuelve la herramienta `langid`, por lo que solo necesita agregar una declaración de inclusión antes del bloque de workflow.

Realice la siguiente edición al workflow:

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

Puede abrir el archivo del módulo para examinar su código:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Usar langid para predecir el idioma de cada archivo de entrada
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

Como puede ver, la definición de entrada usa la misma estructura `tuple val(meta), path(file)` que acabamos de aplicar a nuestro canal de entrada.

La definición de salida está estructurada como una tupla con una estructura similar a la entrada, excepto que también contiene `stdout` como tercer elemento.
Este patrón `tuple val(meta), path(file), <output>` mantiene los metadatos asociados tanto con los datos de entrada como con las salidas mientras fluye a través del pipeline.

Tenga en cuenta que estamos usando el calificador de salida [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) de Nextflow aquí porque la herramienta imprime su salida directamente a la consola en lugar de escribir un archivo; y usamos `sed` en la línea de comando para eliminar la puntuación de probabilidad, limpiar la cadena eliminando caracteres de nueva línea y devolver solo la predicción de idioma.

### 2.2. Agregar una llamada a `IDENTIFY_LANGUAGE`

Ahora que el proceso está disponible para el workflow, podemos agregar una llamada al proceso `IDENTIFY_LANGUAGE` para ejecutarlo en el canal de datos.

Realice las siguientes ediciones al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Ejecutar langid para identificar el idioma de cada saludo
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

Tenga en cuenta que hemos eliminado la operación `.view()` original en la construcción del canal.

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

¡Excelente! Ahora tenemos una predicción de qué idioma habla cada personaje.

Y como se señaló anteriormente, también hemos incluido el archivo de entrada y el mapa de metadatos en la salida, lo que significa que ambos permanecen asociados con la nueva información que acabamos de producir.
Esto resultará útil en el siguiente paso.

!!! note

    De manera más general, este patrón de mantener el mapa de metadatos asociado con los resultados facilita la asociación de resultados relacionados que comparten los mismos identificadores.

    Como ya habrá aprendido, no puede confiar en el orden de los elementos en los canales para emparejar resultados entre ellos.
    En su lugar, debe usar claves para asociar datos correctamente, y los mapas de metadatos proporcionan una estructura ideal para este propósito.

    Exploramos este caso de uso en detalle en la misión secundaria [Dividir y Agrupar](./splitting_and_grouping.md).

### 2.3. Aumentar metadatos con salidas de procesos

Dado que los resultados que acabamos de producir son en sí mismos una forma de metadatos sobre el contenido de los archivos, sería útil agregarlos a nuestro mapa de metadatos.

Sin embargo, no queremos modificar el mapa de metadatos existente en su lugar.
Desde un punto de vista técnico, es _posible_ hacer eso, pero no es seguro.

Así que en su lugar, crearemos un nuevo mapa de metadatos que contenga el contenido del mapa de metadatos existente más un nuevo par clave-valor `lang: lang_id` que contenga la nueva información, usando el operador `+` (una característica de Groovy).
Y combinaremos esto con una operación [`map`](https://www.nextflow.io/docs/latest/operator.html#map) para reemplazar el mapa antiguo con el nuevo.

Aquí están las ediciones que necesita hacer al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Ejecutar langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Ejecutar langid para identificar el idioma de cada saludo
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Si aún no está familiarizado con el operador `+`, o si esto parece confuso, tome unos minutos para revisar la explicación detallada a continuación.

??? info "Creación del nuevo mapa de metadatos usando el operador `+`"

    **Primero, necesita saber que podemos fusionar el contenido de dos mapas usando el operador `+` de Groovy.**

    Digamos que tenemos los siguientes mapas:

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

    **¿Pero qué pasa si necesita agregar un campo que aún no es parte de un mapa?**

    Digamos que comienza nuevamente desde `map1`, pero la predicción de idioma no está en su propio mapa (no hay `map2`).
    En su lugar, se mantiene en una variable llamada `lang_id`, y sabe que quiere almacenar su valor (`'fr'`) con la clave `lang`.

    En realidad puede hacer lo siguiente:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Aquí, `[lang: new_info]` crea un nuevo mapa sin nombre sobre la marcha, y `map1 + ` fusiona `map1` con el nuevo mapa sin nombre, produciendo el mismo contenido `new_map` que antes.

    Ordenado, ¿verdad?

    **Ahora transpongamos eso al contexto de una operación `channel.map()` de Nextflow.**

    El código se convierte en:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Esto hace lo siguiente:

    - `map1, lang_id ->` toma los dos elementos en la tupla
    - `[map1 + [lang: lang_id]]` crea el nuevo mapa como se detalló arriba

    La salida es un solo mapa sin nombre con el mismo contenido que `new_map` en nuestro ejemplo anterior.
    Entonces efectivamente hemos transformado:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    en:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Con suerte puede ver que si cambiamos `map1` a `meta`, eso es básicamente todo lo que necesitamos para agregar la predicción de idioma a nuestro mapa de metadatos en nuestro workflow.

    ¡Excepto por una cosa!

    En el caso de nuestro workflow, **también necesitamos tener en cuenta la presencia del objeto `file` en la tupla**, que está compuesta de `meta, file, lang_id`.

    Entonces el código aquí se convertiría en:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Si está teniendo dificultades para entender por qué el `file` parece estar moviéndose en la operación `map`, imagine que en lugar de `[meta + [lang: lang_id], file]`, esa línea lee `[new_map, file]`.
    Esto debería dejar más claro que simplemente estamos dejando el `file` en su lugar original en segunda posición en la tupla. Simplemente hemos tomado el valor `new_info` y lo hemos incorporado al mapa que está en primera posición.

    **¡Y esto nos lleva de vuelta a la estructura de canal `tuple val(meta), path(file)`!**

Una vez que esté seguro de entender qué está haciendo este código, ejecute el workflow para ver si funcionó:

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

Sí, ¡eso lo confirma!
Hemos reorganizado ordenadamente la salida del proceso de `meta, file, lang_id` para que `lang_id` sea ahora una de las claves en el mapa de metadatos, y las tuplas del canal se ajusten al modelo `meta, file` una vez más.

### 2.4. Asignar un grupo de idiomas usando condicionales

Ahora que tenemos nuestras predicciones de idioma, usemos la información para asignar algunas agrupaciones nuevas.

En nuestros datos de ejemplo, los idiomas usados por nuestros personajes pueden agruparse en idiomas germánicos (inglés, alemán) e idiomas románicos (francés, español, italiano).
Podría ser útil tener esa clasificación fácilmente disponible en algún lugar más adelante en el pipeline, así que agreguemos esa información en el mapa de metadatos.

Y, buenas noticias, ¡este es otro caso que se presta perfectamente para usar el operador `map`!

Específicamente, vamos a definir una variable llamada `lang_group`, usar algo de lógica condicional simple para determinar qué valor asignar al `lang_group` para cada pieza de datos.

La sintaxis general se verá así:

```groovy
.map { meta, file ->

    // la lógica condicional que define lang_group va aquí

    [meta + [lang_group: lang_group], file]
}
```

Puede ver que esto es muy similar a la operación de fusión de mapa sobre la marcha que usamos en el paso anterior.
Solo necesitamos escribir las declaraciones condicionales.

Aquí está la lógica condicional que queremos aplicar:

- Definir una variable llamada `lang_group` con valor predeterminado `'unknown'`.
- Si `lang` es alemán (`'de'`) o inglés (`'en'`), cambiar `lang_group` a `germanic`.
- De lo contrario, si `lang` está incluido en una lista que contiene francés (`'fr'`), español (`'es'`) e italiano (`'it'`), cambiar `lang_group` a `romance`.

Intente escribirlo usted mismo si ya sabe cómo escribir declaraciones condicionales en Nextflow.

!!! tip

    Puede acceder al valor de `lang` dentro de la operación map con `meta.lang`.

Debería terminar haciendo los siguientes cambios al workflow:

=== "Después"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Ejecutar langid para identificar el idioma de cada saludo
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
        // Ejecutar langid para identificar el idioma de cada saludo
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
- Usamos la operación de fusión `meta + [lang_group:lang_group]` como anteriormente para generar el mapa de metadatos actualizado.

Una vez que todo tenga sentido, ejecute el workflow nuevamente para ver el resultado:

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

Como puede ver, los elementos del canal mantienen su estructura `[meta, file]`, pero el mapa de metadatos ahora incluye esta nueva clasificación.

### Conclusión

En esta sección, ha aprendido cómo:

- **Aplicar metadatos de entrada a canales de salida**: Copiar metadatos de esta manera nos permite asociar resultados más adelante basándonos en el contenido de los metadatos.
- **Crear claves personalizadas**: Creó dos nuevas claves en su mapa de metadatos, fusionándolas con `meta + [new_key:value]` en el mapa de metadatos existente. Una basada en un valor calculado de un proceso, y una basada en una condición que estableció en el operador `map`.

Estos le permiten asociar metadatos nuevos y existentes con archivos a medida que avanza a través de su pipeline.
Incluso si no está usando metadatos como parte de un proceso, mantener el mapa de metadatos asociado con los datos de esta manera facilita mantener toda la información relevante junta.

---

## 3. Usando información del mapa de metadatos en un proceso

Ahora que sabe cómo crear y actualizar el mapa de metadatos, podemos llegar a la parte realmente divertida: usar realmente los metadatos en un proceso.

Más específicamente, vamos a agregar un segundo paso a nuestro workflow para dibujar cada animal como arte ASCII y hacer que diga el texto grabado en una burbuja de diálogo.
Vamos a hacer esto usando una herramienta llamada [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "¿Qué hace `cowpy`?"

    `cowpy` es una herramienta de línea de comandos que genera arte ASCII para mostrar entradas de texto arbitrarias de una manera divertida.
    Es una implementación en Python de la herramienta clásica [cowsay](https://en.wikipedia.org/wiki/Cowsay) de Tony Monroe.

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

    Opcionalmente, puede seleccionar un personaje (o 'cowacter') para usar en lugar de la vaca predeterminada.

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

Si trabajó en el curso Hello Nextflow, ya ha visto esta herramienta en acción.
Si no, no se preocupe; cubriremos todo lo que necesita saber a medida que avancemos.

### 3.1. Importar el proceso y examinar el código

Le proporcionamos un módulo de proceso pre-escrito llamado `COWPY` que envuelve la herramienta `cowpy`, por lo que solo necesita agregar una declaración de inclusión antes del bloque de workflow.

Realice la siguiente edición al workflow:

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

Puede abrir el archivo del módulo para examinar su código:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generar arte ASCII con cowpy
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

Como puede ver, este proceso está actualmente diseñado para tomar un archivo de entrada (que contiene el texto a mostrar) y un valor que especifica el personaje que debe dibujarse en arte ASCII, generalmente proporcionado a nivel de workflow por un parámetro de línea de comandos.

### 3.2. Pasar un campo del mapa de metadatos como entrada

Cuando usamos la herramienta `cowpy` en el curso Hello Nextflow, usamos un parámetro de línea de comandos para determinar qué personaje usar para dibujar la imagen final.
Eso tenía sentido, porque solo estábamos generando una imagen por ejecución del pipeline.

Sin embargo, en este tutorial, queremos generar una imagen apropiada para cada sujeto que estamos procesando, por lo que usar un parámetro de línea de comandos sería demasiado limitante.

Buenas noticias: tenemos una columna `character` en nuestra hoja de datos y por lo tanto, en nuestro mapa de metadatos.
Usemos eso para establecer el personaje que el proceso debería usar para cada entrada.

Para ese fin, necesitaremos hacer tres cosas:

1. Dar un nombre al canal de salida que sale del proceso anterior para que podamos operar en él más convenientemente.
2. Determinar cómo acceder a la información de interés
3. Agregar una llamada al segundo proceso y alimentar la información apropiadamente.

Comencemos.

#### 3.2.1. Nombrar el canal de salida anterior

Aplicamos las manipulaciones anteriores directamente en el canal de salida del primer proceso, `IDENTIFY_LANGUAGE.out`.
Para alimentar el contenido de ese canal al siguiente proceso (y hacerlo de una manera que sea clara y fácil de leer) queremos darle su propio nombre, `ch_languages`.

Podemos hacer eso usando el operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

En el workflow principal, reemplace el operador `.view()` con `.set { ch_languages }`, y agregue una línea probando que podemos referirnos al canal por nombre.

=== "Después"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Ejecutar langid para identificar el idioma de cada saludo
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

        // Temporal: mirar dentro de ch_languages
        ch_languages.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Ejecutar langid para identificar el idioma de cada saludo
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

#### 3.2.2. Acceder al archivo y metadatos del personaje

Sabemos por mirar el código del módulo que el proceso `COWPY` espera que se le proporcione un archivo de texto y un valor `character`.
Para escribir la llamada al proceso `COWPY`, solo necesitamos saber cómo extraer el objeto archivo correspondiente y los metadatos de cada elemento en el canal.

Como suele ser el caso, la forma más simple de hacer eso es usar una operación `map`.

Nuestro canal contiene tuplas estructuradas como `[meta, file]`, por lo que podemos acceder al objeto `file` directamente, y podemos acceder al valor `character` almacenado dentro del mapa de metadatos refiriéndonos a él como `meta.character`.

En el workflow principal, realice los siguientes cambios de código:

=== "Después"

    ```groovy title="main.nf" linenums="34"
        // Temporal: acceder al archivo y personaje
        ch_languages.map { meta, file -> file }.view { file -> "Archivo: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Personaje: " + character }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporal: mirar dentro de ch_languages
        ch_languages.view()
    ```

Tenga en cuenta que estamos usando closures (como `{ file -> "Archivo: " + file }`) para hacer la salida de las operaciones `.view` más legible.

Ejecutemos esto:

```bash
nextflow run main.nf -resume
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Personaje: squirrel
    Archivo: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    Archivo: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Personaje: tux
    Personaje: turkey
    Archivo: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    Archivo: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Personaje: sheep
    Personaje: moose
    Personaje: stegosaurus
    Archivo: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    Archivo: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    Archivo: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Personaje: turtle
    ```

_Las rutas de archivo y valores de personaje pueden aparecer en un orden diferente en su salida._

Esto confirma que podemos acceder al archivo y al personaje para cada elemento en el canal.

#### 3.2.3. Llamar al proceso `COWPY`

Ahora pongamos todo junto y llamemos realmente al proceso `COWPY` en el canal `ch_languages`.

En el workflow principal, realice los siguientes cambios de código:

=== "Después"

    ```groovy title="main.nf" linenums="34"
        // Ejecutar cowpy para generar arte ASCII
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporal: acceder al archivo y personaje
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Ve que simplemente copiamos las dos operaciones map (menos las declaraciones `.view()`) como las entradas para la llamada al proceso.
¡Solo asegúrese de no olvidar la coma entre ellas!

Es un poco torpe, pero veremos cómo mejorar eso en la siguiente sección.

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

Si mira en el directorio de resultados, debería ver los archivos individuales que contienen el arte ASCII de cada saludo hablado por el personaje correspondiente.

??? abstract "Directorio y contenido de archivo de ejemplo"

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

Esto muestra que pudimos usar la información en el mapa de metadatos para parametrizar el comando en el segundo paso del pipeline.

Sin embargo, como se señaló anteriormente, parte del código involucrado fue un poco torpe, ya que tuvimos que desempaquetar metadatos mientras todavía estábamos en el contexto del cuerpo del workflow.
Ese enfoque funciona bien para usar un pequeño número de campos del mapa de metadatos, pero escalaría pobremente si quisiéramos usar muchos más.

Hay otro operador llamado `multiMap()` que nos permite optimizar esto un poco, pero incluso entonces no es ideal.

??? info "(Opcional) Versión alternativa con `multiMap()`"

    En caso de que se lo esté preguntando, no podíamos simplemente escribir una sola operación `map()` que produjera tanto el `file` como el `character`, porque eso los devolvería como una tupla.
    Tuvimos que escribir dos operaciones `map()` separadas para alimentar los elementos `file` y `character` al proceso por separado.

    Técnicamente hay otra forma de hacer esto a través de una sola operación de mapeo, usando el operador `multiMap()`, que es capaz de emitir múltiples canales.
    Por ejemplo, podría reemplazar la llamada a `COWPY` anterior con el siguiente código:

    === "Después"

        ```groovy title="main.nf" linenums="34"
            // Ejecutar cowpy para generar arte ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Antes"

        ```groovy title="main.nf" linenums="34"
            // Ejecutar cowpy para generar arte ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Esto produce exactamente el mismo resultado.

En cualquier caso, es incómodo que tengamos que hacer algo de desempaquetado a nivel de workflow.

Sería mejor si pudiéramos alimentar todo el mapa de metadatos en el proceso y seleccionar lo que necesitamos una vez allí.

### 3.3. Pasar y usar todo el mapa de metadatos

El punto del mapa de metadatos es después de todo pasar todos los metadatos juntos como un paquete.
La única razón por la que no pudimos hacer eso anteriormente es que el proceso no está configurado para aceptar un mapa de metadatos.
Pero como controlamos el código del proceso, podemos cambiar eso.

Modifiquemos el proceso `COWPY` para aceptar la estructura de tupla `[meta, file]` que usamos en el primer proceso para poder optimizar el workflow.

Para ese fin, necesitaremos hacer tres cosas:

1. Modificar las definiciones de entrada del módulo de proceso `COWPY`
2. Actualizar el comando del proceso para usar el mapa de metadatos
3. Actualizar la llamada al proceso en el cuerpo del workflow

¿Listo? ¡Vamos!

#### 3.3.1. Modificar la entrada del módulo `COWPY`

Realice las siguientes ediciones al archivo de módulo `cowpy.nf`:

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

Tenga en cuenta que no actualizamos la definición de salida del proceso para producir el mapa de metadatos, con el fin de mantener el tutorial simplificado, pero siéntase libre de hacerlo usted mismo como ejercicio siguiendo el modelo del proceso `IDENTIFY_LANGUAGE`.

#### 3.3.2. Actualizar el comando para usar el campo del mapa de metadatos

Todo el mapa de metadatos ahora está disponible dentro del proceso, por lo que podemos referirnos a la información que contiene directamente desde dentro del bloque de comando.

Realice las siguientes ediciones al archivo de módulo `cowpy.nf`:

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

Hemos reemplazado la referencia al valor `character` previamente pasado como una entrada independiente con el valor contenido en el mapa de metadatos, al que nos referimos usando `meta.character`.

Ahora actualicemos la llamada al proceso en consecuencia.

#### 3.3.3. Actualizar la llamada al proceso y ejecutarlo

El proceso ahora espera que su entrada use la estructura de tupla `[meta, file]`, que es lo que produce el proceso anterior, por lo que simplemente podemos pasar todo el canal `ch_languages` al proceso `COWPY`.

Realice las siguientes ediciones al workflow principal:

=== "Después"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Ejecutar cowpy para generar arte ASCII
    COWPY(ch_languages)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Ejecutar cowpy para generar arte ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
