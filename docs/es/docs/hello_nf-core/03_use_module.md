# Parte 3: Usar un módulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta tercera parte del curso de entrenamiento Hello nf-core, le mostramos cómo encontrar, instalar y usar un módulo nf-core existente en su pipeline.

Uno de los grandes beneficios de trabajar con nf-core es la capacidad de aprovechar módulos preconstruidos y probados del repositorio [nf-core/modules](https://github.com/nf-core/modules).
En lugar de escribir cada proceso desde cero, puede instalar y usar módulos mantenidos por la comunidad que siguen las mejores prácticas.

Para demostrar cómo funciona esto, reemplazaremos el módulo personalizado `collectGreetings` con el módulo `cat/cat` de nf-core/modules en el pipeline `core-hello`.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que ha completado la [Parte 2: Reescribir Hello para nf-core](./02_rewrite_hello.md) y tiene un pipeline `core-hello` funcional.

    Si no completó la Parte 2 o desea comenzar desde cero para esta parte, puede usar la solución `core-hello-part2` como punto de partida.
    Ejecute este comando desde el directorio `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Esto le proporciona un pipeline nf-core completamente funcional listo para agregar módulos.
    Puede probar que se ejecuta correctamente ejecutando el siguiente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Encontrar e instalar un módulo nf-core adecuado

Primero, aprendamos cómo encontrar un módulo nf-core existente e instalarlo en nuestro pipeline.

Nuestro objetivo será reemplazar el proceso `collectGreetings`, que utiliza el comando Unix `cat` para concatenar múltiples archivos de saludo en uno.
Concatenar archivos es una operación muy común, por lo que es razonable pensar que podría haber ya un módulo en nf-core diseñado para ese propósito.

Comencemos.

### 1.1. Explorar los módulos disponibles en el sitio web de nf-core

El proyecto nf-core mantiene un catálogo centralizado de módulos en [https://nf-co.re/modules](https://nf-co.re/modules).

Navegue a la página de módulos en su navegador web y use la barra de búsqueda para buscar 'concatenate'.

![module search results](../../../en/docs/hello_nf-core/img/module-search-results.png)

Como puede ver, hay bastantes resultados, muchos de ellos módulos diseñados para concatenar tipos de archivos muy específicos.
Entre ellos, debería ver uno llamado `cat_cat` que es de propósito general.

!!! note "Convención de nomenclatura de módulos"

    El guion bajo (`_`) se usa como sustituto del carácter de barra diagonal (`/`) en los nombres de los módulos.

    Los módulos nf-core siguen la convención de nomenclatura `software/comando` cuando una herramienta proporciona múltiples comandos, como `samtools/view` (paquete samtools, comando view) o `gatk/haplotypecaller` (paquete GATK, comando HaplotypeCaller).
    Para herramientas que proporcionan solo un comando principal, los módulos usan un solo nivel como `fastqc` o `multiqc`.

Haga clic en el recuadro del módulo `cat_cat` para ver la documentación del módulo.

La página del módulo muestra:

- Una breve descripción: "A module for concatenation of gzipped or uncompressed files"
- Comando de instalación: `nf-core modules install cat/cat`
- Estructura de canales de entrada y salida
- Parámetros disponibles

### 1.2. Listar los módulos disponibles desde la línea de comandos

Alternativamente, también puede buscar módulos directamente desde la línea de comandos usando las herramientas nf-core.

```bash
nf-core modules list remote
```

Esto mostrará una lista de todos los módulos disponibles en el repositorio nf-core/modules, aunque es un poco menos conveniente si no conoce ya el nombre del módulo que está buscando.
Sin embargo, si lo conoce, puede canalizar la lista a `grep` para encontrar módulos específicos:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Salida del comando"

    ```console
    │ cat/cat
    ```

Tenga en cuenta que el enfoque de `grep` solo mostrará resultados con el término de búsqueda en su nombre, lo que no funcionaría para `cat_cat`.

### 1.3. Obtener información detallada sobre el módulo

Para ver información detallada sobre un módulo específico desde la línea de comandos, use el comando `info`:

```bash
nf-core modules info cat/cat
```

Esto muestra documentación sobre el módulo, incluyendo sus entradas, salidas e información básica de uso.

??? success "Salida del comando"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Esta es exactamente la misma información que puede encontrar en el sitio web.

### 1.4. Instalar el módulo cat/cat

Ahora que hemos encontrado el módulo que queremos, necesitamos agregarlo al código fuente de nuestro pipeline.

La buena noticia es que el proyecto nf-core incluye algunas herramientas para facilitar esta parte.
Específicamente, el comando `nf-core modules install` permite automatizar la recuperación del código y hacerlo disponible para su proyecto en un solo paso.

Navegue al directorio de su pipeline y ejecute el comando de instalación:

```bash
cd core-hello
nf-core modules install cat/cat
```

La herramienta puede primero solicitarle que especifique un tipo de repositorio.
(Si no lo hace, salte a "Finalmente, la herramienta procederá a instalar el módulo.")

??? success "Salida del comando"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Si es así, presione enter para aceptar la respuesta predeterminada (`Pipeline`) y continuar.

Luego, la herramienta ofrecerá modificar la configuración de su proyecto para evitar esta solicitud en el futuro.

??? success "Salida del comando"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

¡Más vale aprovechar esta práctica herramienta!
Presione enter para aceptar la respuesta predeterminada (sí).

Finalmente, la herramienta procederá a instalar el módulo.

??? success "Salida del comando"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

El comando automáticamente:

- Descarga los archivos del módulo a `modules/nf-core/cat/cat/`
- Actualiza `modules.json` para rastrear el módulo instalado
- Le proporciona la declaración `include` correcta para usar en su workflow

!!! tip

    Asegúrese siempre de que su directorio de trabajo actual sea la raíz de su proyecto de pipeline antes de ejecutar el comando de instalación del módulo.

Verifiquemos que el módulo se instaló correctamente:

```bash
tree -L 4 modules
```

??? abstract "Contenidos del directorio"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

También puede verificar la instalación pidiendo a la utilidad nf-core que liste los módulos instalados localmente:

```bash
nf-core modules list local
```

??? success "Salida del comando"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Esto confirma que el módulo `cat/cat` ahora es parte del código fuente de su proyecto.

Sin embargo, para realmente usar el nuevo módulo, necesitamos importarlo en nuestro pipeline.

### 1.5. Actualizar las importaciones de módulos

Reemplacemos la declaración `include` para el módulo `collectGreetings` con la de `CAT_CAT` en la sección de importaciones del workflow `workflows/hello.nf`.

Como recordatorio, la herramienta de instalación del módulo nos dio la declaración exacta a usar:

```groovy title="Declaración de importación producida por el comando install"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Tenga en cuenta que la convención nf-core es usar mayúsculas para los nombres de los módulos al importarlos.

Abra [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) y realice la siguiente sustitución:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
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

Observe cómo la ruta para el módulo nf-core difiere de los módulos locales:

- **Módulo nf-core**: `'../modules/nf-core/cat/cat/main'` (referencia a `main.nf`)
- **Módulo local**: `'../modules/local/collectGreetings.nf'` (referencia a archivo único)

El módulo ahora está disponible para el workflow, así que todo lo que necesitamos hacer es intercambiar la llamada a `collectGreetings` para usar `CAT_CAT`. ¿Correcto?

No tan rápido.

En este punto, podría sentirse tentado a sumergirse y comenzar a editar código, pero vale la pena tomarse un momento para examinar cuidadosamente qué espera el nuevo módulo y qué produce.

Vamos a abordar eso como una sección separada porque involucra un nuevo mecanismo que no hemos cubierto todavía: los mapas de metadatos.

!!! note

    Opcionalmente puede eliminar el archivo `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Sin embargo, podría querer mantenerlo como referencia para comprender las diferencias entre los módulos locales y los de nf-core.

### Conclusión

Sabe cómo encontrar un módulo nf-core y hacerlo disponible para su proyecto.

### ¿Qué sigue?

Evaluar qué requiere un nuevo módulo e identificar cualquier cambio importante necesario para integrarlo en un pipeline.

---

## 2. Evaluar los requisitos del nuevo módulo

Específicamente, necesitamos examinar la **interfaz** del módulo, es decir, sus definiciones de entrada y salida, y compararla con la interfaz del módulo que buscamos reemplazar.
Esto nos permitirá determinar si podemos simplemente tratar el nuevo módulo como un reemplazo directo o si necesitaremos adaptar parte del cableado.

Idealmente, esto es algo que debería hacer _antes_ de instalar el módulo, pero bueno, más vale tarde que nunca.
(Para su información, existe un comando `uninstall` para deshacerse de los módulos que decide que ya no quiere.)

!!! note

    El proceso CAT_CAT incluye un manejo bastante ingenioso de diferentes tipos de compresión, extensiones de archivo, etc., que no son estrictamente relevantes para lo que estamos tratando de mostrarle aquí, así que ignoraremos la mayor parte y nos enfocaremos solo en las partes que son importantes.

### 2.1. Comparar las interfaces de los dos módulos

Como recordatorio, esto es cómo se ve la interfaz de nuestro módulo `collectGreetings`:

```groovy title="modules/local/collectGreetings.nf (extracto)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

El módulo `collectGreetings` toma dos entradas:

- `input_files` contiene uno o más archivos de entrada para procesar;
- `batch_name` es un valor que usamos para asignar un nombre específico de ejecución al archivo de salida, que es una forma de metadatos.

Al completarse, `collectGreetings` produce una sola ruta de archivo, emitida con la etiqueta `outfile`.

En comparación, la interfaz del módulo `cat/cat` es más compleja:

```groovy title="modules/nf-core/cat/cat/main.nf (extracto)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

El módulo CAT_CAT toma una sola entrada, pero esa entrada es una tupla que contiene dos cosas:

- `meta` es una estructura que contiene metadatos, llamada metamapa;
- `files_in` contiene uno o más archivos de entrada para procesar, equivalente a `input_files` de `collectGreetings`.

Al completarse, CAT_CAT entrega sus salidas en dos partes:

- Otra tupla que contiene el metamapa y el archivo de salida concatenado, emitido con la etiqueta `file_out`;
- Un archivo `versions.yml` que captura información sobre la versión del software que se usó, emitido con la etiqueta `versions`.

Observe también que, por defecto, el archivo de salida se nombrará en base a un identificador que es parte de los metadatos (código no mostrado aquí).

Esto puede parecer mucho para llevar un registro solo mirando el código, así que aquí hay un diagrama para ayudarlo a visualizar cómo todo encaja.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Puede ver que los dos módulos tienen requisitos de entrada similares en términos de contenido (un conjunto de archivos de entrada más algunos metadatos) pero expectativas muy diferentes sobre cómo se empaqueta ese contenido.
Ignorando el archivo de versiones por ahora, su salida principal también es equivalente (un archivo concatenado), excepto que CAT_CAT también emite el metamapa junto con el archivo de salida.

Las diferencias de empaquetado serán bastante fáciles de manejar, como verá en un momento.
Sin embargo, para entender la parte del metamapa, necesitamos presentarle algo de contexto adicional.

### 2.2. Comprender los metamapas

Acabamos de decirle que el módulo CAT_CAT espera un mapa de metadatos como parte de su tupla de entrada.
Tomemos unos minutos para examinar más de cerca qué es eso.

El **mapa de metadatos**, a menudo referido como **metamapa** para abreviar, es un mapa estilo Groovy que contiene información sobre unidades de datos.
En el contexto de los pipelines Nextflow, las unidades de datos pueden ser lo que usted quiera: muestras individuales, lotes de muestras o conjuntos de datos completos.

Por convención, un metamapa nf-core se llama `meta` y contiene el campo requerido `id`, que se usa para nombrar salidas y rastrear unidades de datos.

Por ejemplo, un mapa de metadatos típico podría verse así:

```groovy title="Ejemplo de metamapa a nivel de muestra"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

O en un caso donde los metadatos se adjuntan a nivel de lote:

```groovy title="Ejemplo de metamapa a nivel de lote"
[id: 'batch1', date: '25.10.01']
```

Ahora pongamos esto en el contexto del proceso `CAT_CAT`, que espera que los archivos de entrada se empaqueten en una tupla con un metamapa, y también produce el metamapa como parte de la tupla de salida.

```groovy title="modules/nf-core/cat/cat/main.nf (extracto)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Como resultado, cada unidad de datos viaja a través del pipeline con los metadatos relevantes adjuntos.
Los procesos subsecuentes pueden entonces acceder fácilmente a esos metadatos también.

¿Recuerda cómo le dijimos que el archivo producido por `CAT_CAT` se nombrará en base a un identificador que es parte de los metadatos?
Este es el código relevante:

```groovy title="modules/nf-core/cat/cat/main.nf (extracto)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Esto se traduce aproximadamente de la siguiente manera: si se proporciona un `prefix` a través del sistema de parámetros externos de tarea (`task.ext`), úselo para nombrar el archivo de salida; de lo contrario, cree uno usando `${meta.id}`, que corresponde al campo `id` en el metamapa.

Puede imaginar el canal de entrada llegando a este módulo con contenidos como este:

```groovy title="Ejemplo de contenidos del canal de entrada"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Luego, los contenidos del canal de salida saliendo así:

```groovy title="Ejemplo de contenidos del canal de salida"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Como se mencionó anteriormente, la configuración de entrada `tuple val(meta), path(files_in)` es un patrón estándar usado en todos los módulos nf-core.

Con suerte, puede comenzar a ver qué tan útil puede ser esto.
No solo le permite nombrar salidas basadas en metadatos, sino que también puede hacer cosas como usarlo para aplicar diferentes valores de parámetros, y en combinación con operadores específicos, incluso puede agrupar, ordenar o filtrar datos a medida que fluyen a través del pipeline.

!!! note "Aprenda más sobre metadatos"

    Para una introducción completa sobre cómo trabajar con metadatos en workflows Nextflow, incluyendo cómo leer metadatos de hojas de muestras y usarlos para personalizar el procesamiento, consulte la [Metadatos en workflows](../side_quests/metadata) misión secundaria.

### 2.3. Resumir los cambios a realizar

Basándonos en lo que hemos revisado, estos son los cambios principales que necesitamos hacer a nuestro pipeline para utilizar el módulo `cat/cat`:

- Crear un metamapa que contenga el nombre del lote;
- Empaquetar el metamapa en una tupla con el conjunto de archivos de entrada a concatenar (provenientes de `convertToUpper`);
- Cambiar la llamada de `collectGreetings()` a `CAT_CAT`;
- Extraer el archivo de salida de la tupla producida por el proceso `CAT_CAT` antes de pasarlo a `cowpy`.

¡Eso debería funcionar! Ahora que tenemos un plan, estamos listos para sumergirnos.

### Conclusión

Sabe cómo evaluar la interfaz de entrada y salida de un nuevo módulo para identificar sus requisitos, y ha aprendido cómo los metamapas son usados por los pipelines nf-core para mantener los metadatos estrechamente asociados con los datos a medida que fluyen a través de un pipeline.

### ¿Qué sigue?

Integrar el nuevo módulo en un workflow.

---

## 3. Integrar CAT_CAT en el workflow `hello.nf`

Ahora que sabe todo sobre los metamapas (o suficiente para los propósitos de este curso, de todos modos), es hora de implementar realmente los cambios que describimos anteriormente.

Por claridad, desglosaremos esto y cubriremos cada paso por separado.

!!! note

    Todos los cambios mostrados a continuación se realizan en la lógica del workflow en el bloque `main` en el archivo de workflow `core-hello/workflows/hello.nf`.

### 3.1. Crear un mapa de metadatos

Primero, necesitamos crear un mapa de metadatos para `CAT_CAT`, teniendo en cuenta que los módulos nf-core requieren que el metamapa tenga al menos un campo `id`.

Dado que no necesitamos otros metadatos, podemos mantenerlo simple y usar algo como esto:

```groovy title="Ejemplo de sintaxis"
def cat_meta = [id: 'test']
```

Excepto que no queremos codificar el valor de `id`; queremos usar el valor del parámetro `params.batch`.
Entonces el código se convierte en:

```groovy title="Ejemplo de sintaxis"
def cat_meta = [id: params.batch]
```

Sí, es literalmente así de simple crear un metamapa básico.

Agreguemos estas líneas después de la llamada a `convertToUpper`, eliminando la llamada a `collectGreetings`:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // recopilar todos los saludos en un archivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Esto crea un mapa de metadatos simple donde el `id` se establece en nuestro nombre de lote (que será `test` cuando se use el perfil de prueba).

### 3.2. Crear un canal con tuplas de metadatos

A continuación, transforme el canal de archivos en un canal de tuplas que contengan metadatos y archivos:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // crear un canal con metadatos y archivos en formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

La línea que hemos agregado logra dos cosas:

- `.collect()` reúne todos los archivos de la salida de `convertToUpper` en una sola lista
- `.map { files -> tuple(cat_meta, files) }` crea una tupla de `[metadatos, archivos]` en el formato que `CAT_CAT` espera

Eso es todo lo que necesitamos hacer para configurar la tupla de entrada para `CAT_CAT`.

### 3.3. Llamar al módulo CAT_CAT

Ahora llame a `CAT_CAT` en el canal recién creado:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // crear un canal con metadatos y archivos en formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenar archivos usando el módulo nf-core cat/cat
        CAT_CAT(ch_for_cat)

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // crear un canal con metadatos y archivos en formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Esto completa la parte más complicada de esta sustitución, pero aún no hemos terminado: todavía necesitamos actualizar cómo pasamos la salida concatenada al proceso `cowpy`.

### 3.4. Extraer el archivo de salida de la tupla para `cowpy`

Anteriormente, el proceso `collectGreetings` simplemente producía un archivo que podíamos pasar a `cowpy` directamente.
Sin embargo, el proceso `CAT_CAT` produce una tupla que incluye el metamapa además del archivo de salida.

Dado que `cowpy` aún no acepta tuplas de metadatos (arreglaremos esto en la siguiente parte del curso), necesitamos extraer el archivo de salida de la tupla producida por `CAT_CAT` antes de pasarlo a `cowpy`:

=== "Después"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // crear un canal con metadatos y archivos en formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenar los saludos
        CAT_CAT(ch_for_cat)

        // extraer el archivo de la tupla ya que cowpy aún no usa metadatos
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generar arte ASCII de los saludos con cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emitir un saludo
        sayHello(ch_samplesheet)

        // convertir el saludo a mayúsculas
        convertToUpper(sayHello.out)

        // crear mapa de metadatos con el nombre del lote como ID
        def cat_meta = [ id: params.batch ]

        // crear un canal con metadatos y archivos en formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenar los saludos
        CAT_CAT(ch_for_cat)

        // generar arte ASCII de los saludos con cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

La operación `.map{ meta, file -> file }` extrae el archivo de la tupla `[metadatos, archivo]` producida por `CAT_CAT` en un nuevo canal, `ch_for_cowpy`.

Luego, es solo cuestión de pasar `ch_for_cowpy` a `cowpy` en lugar de `collectGreetings.out.outfile` en esa última línea.

!!! note

    En la siguiente parte del curso, actualizaremos `cowpy` para trabajar con tuplas de metadatos directamente, por lo que este paso de extracción ya no será necesario.

### 3.5. Probar el workflow

Probemos que el workflow funciona con el módulo `cat/cat` recién integrado:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Esto debería ejecutarse razonablemente rápido.

??? success "Salida del comando"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Observe que `CAT_CAT` ahora aparece en la lista de ejecución de procesos en lugar de `collectGreetings`.

¡Y eso es todo! Ahora estamos usando un módulo robusto mantenido por la comunidad en lugar de código personalizado de grado prototipo para ese paso en el pipeline.

### Conclusión

Ahora sabe cómo:

- Encontrar e instalar módulos nf-core
- Evaluar los requisitos de un módulo nf-core
- Crear un mapa de metadatos simple para usar con un módulo nf-core
- Integrar un módulo nf-core en su workflow

### ¿Qué sigue?

Aprenda a adaptar sus módulos locales para seguir las convenciones nf-core.
También le mostraremos cómo crear nuevos módulos nf-core a partir de una plantilla usando las herramientas nf-core.
