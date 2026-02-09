# Parte 5: Validación de entradas

En esta quinta parte del curso de capacitación Hello nf-core, te mostramos cómo usar el plugin nf-schema para validar las entradas y parámetros del pipeline.

??? info "Cómo comenzar desde esta sección"

    Esta sección asume que has completado la [Parte 4: Crear un módulo nf-core](./04_make_module.md) y has actualizado el módulo del proceso `COWPY` a los estándares nf-core en tu pipeline.

    Si no completaste la Parte 4 o quieres comenzar de nuevo para esta parte, puedes usar la solución `core-hello-part4` como punto de partida.
    Ejecuta estos comandos desde dentro del directorio `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part4 core-hello
    cd core-hello
    ```

    Esto te proporciona un pipeline con el módulo `COWPY` ya actualizado para seguir los estándares nf-core.
    Puedes verificar que se ejecuta correctamente ejecutando el siguiente comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 0. Calentamiento: Un poco de contexto

### 0.1. Por qué importa la validación

Imagina ejecutar tu pipeline durante dos horas, solo para que falle porque un usuario proporcionó un archivo con la extensión incorrecta. O pasar horas depurando errores crípticos, solo para descubrir que un parámetro estaba mal escrito. Sin validación de entradas, estos escenarios son comunes.

Considera este ejemplo:

```console title="Without validation"
$ nextflow run my-pipeline --input data.txt --output results

...2 hours later...

ERROR ~ No such file: 'data.fq.gz'
  Expected FASTQ format but received TXT
```

El pipeline aceptó entradas inválidas y se ejecutó durante horas antes de fallar. Con validación adecuada:

```console title="With validation"
$ nextflow run my-pipeline --input data.txt --output results

ERROR ~ Validation of pipeline parameters failed!

 * --input (data.txt): File extension '.txt' does not match required pattern '.fq.gz' or '.fastq.gz'
 * --output: required parameter is missing (expected: --outdir)

Pipeline failed before execution - please fix the errors above
```

El pipeline falla inmediatamente con mensajes de error claros y accionables. Esto ahorra tiempo, recursos de cómputo y frustración.

### 0.2. El plugin nf-schema

El [plugin nf-schema](https://nextflow-io.github.io/nf-schema/latest/) es un plugin de Nextflow que proporciona capacidades de validación integrales para pipelines de Nextflow.
Aunque nf-schema funciona con cualquier workflow de Nextflow, es la solución de validación estándar para todos los pipelines nf-core.

nf-schema proporciona varias funciones clave:

- **Validación de parámetros**: Valida los parámetros del pipeline contra `nextflow_schema.json`
- **Validación de hojas de muestras**: Valida archivos de entrada contra `assets/schema_input.json`
- **Conversión de canales**: Convierte hojas de muestras validadas a canales de Nextflow
- **Generación de texto de ayuda**: Genera automáticamente la salida de `--help` a partir de las definiciones del schema
- **Resumen de parámetros**: Muestra qué parámetros difieren de los valores predeterminados

nf-schema es el sucesor del plugin nf-validation obsoleto y utiliza el estándar [JSON Schema Draft 2020-12](https://json-schema.org/) para la validación.

??? info "¿Qué son los plugins de Nextflow?"

    Los plugins son extensiones que agregan nueva funcionalidad al lenguaje Nextflow en sí. Se instalan a través de un bloque `plugins{}` en `nextflow.config` y pueden proporcionar:

    - Nuevas funciones y clases que se pueden importar (como `samplesheetToList`)
    - Nuevas características DSL y operadores
    - Integración con servicios externos

    El plugin nf-schema se especifica en `nextflow.config`:

    ```groovy
    plugins {
        id 'nf-schema@2.1.1'
    }
    ```

    Una vez instalado, puedes importar funciones de plugins usando la sintaxis `include { functionName } from 'plugin/plugin-name'`.

### 0.3. Dos archivos schema para dos tipos de validación

Un pipeline nf-core utilizará dos archivos schema separados, que corresponden a dos tipos de validación:

| Archivo Schema             | Propósito                  | Valida                                               |
| -------------------------- | -------------------------- | ---------------------------------------------------- |
| `nextflow_schema.json`     | Validación de parámetros   | Flags de línea de comandos: `--input`, `--outdir`, `--batch` |
| `assets/schema_input.json` | Validación de datos de entrada | Contenidos de hojas de muestras y archivos de entrada |

Ambos schemas usan el formato JSON Schema, un estándar ampliamente adoptado para describir y validar estructuras de datos.

**La validación de parámetros** valida parámetros de línea de comandos (flags como `--outdir`, `--batch`, `--input`):

- Verifica tipos, rangos y formatos de parámetros
- Asegura que se proporcionen los parámetros requeridos
- Valida que las rutas de archivos existan
- Definida en `nextflow_schema.json`

**La validación de datos de entrada** valida la estructura de hojas de muestras y archivos manifest (archivos CSV/TSV que describen tus datos):

- Verifica la estructura de columnas y tipos de datos
- Valida que las rutas de archivos referenciadas en la hoja de muestras existan
- Asegura que los campos requeridos estén presentes
- Definida en `assets/schema_input.json`

!!! warning "Lo que la validación de datos de entrada NO hace"

    La validación de datos de entrada verifica la estructura de *archivos manifest* (hojas de muestras, archivos CSV), NO los contenidos de tus archivos de datos reales (FASTQ, BAM, VCF, etc.).

    Para datos a gran escala, validar los contenidos de archivos (como verificar la integridad de BAM) debería ocurrir en procesos del pipeline ejecutándose en nodos de trabajo, no durante la etapa de validación en la máquina de orquestación.

### 0.4. ¿Cuándo debería ocurrir la validación?

```mermaid
graph LR
    A[User runs pipeline] --> B[Parameter validation]
    B -->|✓ Valid| C[Input data validation]
    B -->|✗ Invalid| D[Error: Fix parameters]
    C -->|✓ Valid| E[Pipeline executes]
    C -->|✗ Invalid| F[Error: Fix input data]
```

La validación debería ocurrir **antes** de que se ejecuten los procesos del pipeline, para proporcionar retroalimentación rápida y prevenir tiempo de cómputo desperdiciado.

Ahora apliquemos estos principios en la práctica, comenzando con la validación de parámetros.

---

## 1. Validación de parámetros (nextflow_schema.json)

Comencemos agregando validación de parámetros a nuestro pipeline. Esto valida flags de línea de comandos como `--input`, `--outdir` y `--batch`.

### 1.1. Configurar la validación para omitir la validación de archivos de entrada

La plantilla de pipeline nf-core viene con nf-schema ya instalado y configurado:

- El plugin nf-schema se instala a través del bloque `plugins{}` en `nextflow.config`
- La validación de parámetros está habilitada por defecto a través de `params.validate_params = true`
- La validación se realiza mediante el subworkflow `UTILS_NFSCHEMA_PLUGIN` durante la inicialización del pipeline

El comportamiento de validación se controla a través del scope `validation{}` en `nextflow.config`.

Como trabajaremos primero en la validación de parámetros (esta sección) y no configuraremos el schema de datos de entrada hasta la sección 2, necesitamos decirle temporalmente a nf-schema que omita la validación de los contenidos del archivo del parámetro `input`.

Abre `nextflow.config` y encuentra el bloque `validation` (alrededor de la línea 246). Agrega `ignoreParams` para omitir la validación de archivos de entrada:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="3" linenums="246"
    validation {
        defaultIgnoreParams = ["genomes"]
        ignoreParams = ['input']
        monochromeLogs = params.monochrome_logs
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="246"
    validation {
        defaultIgnoreParams = ["genomes"]
        monochromeLogs = params.monochrome_logs
    }
    ```

Esta configuración le dice a nf-schema que:

- **`defaultIgnoreParams`**: Omita la validación de parámetros complejos como `genomes` (establecido por los desarrolladores de la plantilla)
- **`ignoreParams`**: Omita la validación de los contenidos del archivo del parámetro `input` (temporal; volveremos a habilitar esto en la sección 2)
- **`monochromeLogs`**: Deshabilite la salida coloreada en mensajes de validación cuando se establezca en `true` (controlado por `params.monochrome_logs`)

!!! note "¿Por qué ignorar el parámetro input?"

    El parámetro `input` en `nextflow_schema.json` tiene `"schema": "assets/schema_input.json"` que le dice a nf-schema que valide los *contenidos* del archivo CSV de entrada contra ese schema.
    Como aún no hemos configurado ese schema, temporalmente ignoramos esta validación.
    Eliminaremos esta configuración en la sección 2 después de configurar el schema de datos de entrada.

### 1.2. Examinar el schema de parámetros

Veamos una sección del archivo `nextflow_schema.json` que vino con nuestra plantilla de pipeline:

```bash
grep -A 25 '"input_output_options"' nextflow_schema.json
```

El schema de parámetros está organizado en grupos. Aquí está el grupo `input_output_options`:

```json title="core-hello/nextflow_schema.json (excerpt)" linenums="8"
        "input_output_options": {
            "title": "Input/output options",
            "type": "object",
            "fa_icon": "fas fa-terminal",
            "description": "Define where the pipeline should find input data and save output data.",
            "required": ["input", "outdir"],
            "properties": {
                "input": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "schema": "assets/schema_input.json",
                    "mimetype": "text/csv",
                    "pattern": "^\\S+\\.csv$",
                    "description": "Path to comma-separated file containing information about the samples in the experiment.",
                    "help_text": "You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.",
                    "fa_icon": "fas fa-file-csv"
                },
                "outdir": {
                    "type": "string",
                    "format": "directory-path",
                    "description": "The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
                    "fa_icon": "fas fa-folder-open"
                }
            }
        },
```

Cada entrada descrita aquí tiene las siguientes propiedades clave que pueden ser validadas:

- **`type`**: Tipo de dato (string, integer, boolean, number)
- **`format`**: Formatos especiales como `file-path` o `directory-path`
- **`exists`**: Para rutas de archivos, verifica si el archivo existe
- **`pattern`**: Expresión regular que el valor debe cumplir
- **`required`**: Array de nombres de parámetros que deben proporcionarse
- **`mimetype`**: Mimetype de archivo esperado para validación

Si tienes buen ojo, podrías notar que el parámetro de entrada `batch` que hemos estado usando aún no está definido en el schema.
Lo vamos a agregar en la siguiente sección.

??? info "¿De dónde vienen los parámetros del schema?"

    La validación del schema usa `nextflow.config` como base para las definiciones de parámetros.
    Los parámetros declarados en otros lugares de tus scripts de workflow (como en `main.nf` o archivos de módulos) **no** son recogidos automáticamente por el validador del schema.

    Esto significa que siempre debes declarar los parámetros de tu pipeline en `nextflow.config`, y luego definir sus reglas de validación en `nextflow_schema.json`.

### 1.3. Agregar el parámetro batch

Aunque el schema es un archivo JSON que puede editarse manualmente, **la edición manual es propensa a errores y no se recomienda**.
En su lugar, nf-core proporciona una herramienta GUI interactiva que maneja la sintaxis JSON Schema por ti y valida tus cambios:

```bash
nf-core pipelines schema build
```

Deberías ver algo como esto:

```console
                                      ,--./,-.
      ___     __   __   __   ___     /,-._.--\
|\ | |__  __ /  ` /  \ |__) |__         }  {
| \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                      `._,._,'

nf-core/tools version 3.4.1 - https://nf-co.re

INFO     [✓] Default parameters match schema validation
INFO     [✓] Pipeline schema looks valid (found 17 params)
INFO     Writing schema with 17 params: 'nextflow_schema.json'
🚀  Launch web builder for customisation and editing? [y/n]:
```

Escribe `y` y presiona Enter para lanzar la interfaz web interactiva.

Tu navegador se abrirá mostrando el constructor de schema de parámetros:

![Interfaz del constructor de schema](./img/schema_build.png)

Para agregar el parámetro `batch`:

1. Haz clic en el botón **"Add parameter"** en la parte superior
2. Usa el asa de arrastre (⋮⋮) para mover el nuevo parámetro hacia arriba en el grupo "Input/output options", debajo del parámetro `input`
3. Completa los detalles del parámetro:
   - **ID**: `batch`
   - **Description**: `Name for this batch of greetings`
   - **Type**: `string`
   - **Required**: marca la casilla
   - Opcionalmente, selecciona un ícono del selector de íconos (ej., `fas fa-layer-group`)

![Agregando el parámetro batch](./img/schema_add.png)

Cuando termines, haz clic en el botón **"Finished"** en la parte superior derecha.

De vuelta en tu terminal, verás:

```console
INFO     Writing schema with 18 params: 'nextflow_schema.json'
⣾ Use ctrl+c to stop waiting and force exit.
```

Presiona `Ctrl+C` para salir del constructor de schema.

La herramienta ahora ha actualizado tu archivo `nextflow_schema.json` con el nuevo parámetro `batch`, manejando toda la sintaxis JSON Schema correctamente.

### 1.4. Verificar los cambios

```bash
grep -A 25 '"input_output_options"' nextflow_schema.json
```

```json title="core-hello/nextflow_schema.json (excerpt)" linenums="8" hl_lines="19-23"
    "input_output_options": {
      "title": "Input/output options",
      "type": "object",
      "fa_icon": "fas fa-terminal",
      "description": "Define where the pipeline should find input data and save output data.",
      "required": ["input", "outdir", "batch"],
      "properties": {
        "input": {
          "type": "string",
          "format": "file-path",
          "exists": true,
          "schema": "assets/schema_input.json",
          "mimetype": "text/csv",
          "pattern": "^\\S+\\.csv$",
          "description": "Path to comma-separated file containing information about the samples in the experiment.",
          "help_text": "You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.",
          "fa_icon": "fas fa-file-csv"
        },
        "batch": {
          "type": "string",
          "description": "Name for this batch of greetings",
          "fa_icon": "fas fa-layer-group"
        },
```

Deberías ver que el parámetro `batch` ha sido agregado al schema con el campo "required" ahora mostrando `["input", "outdir", "batch"]`.

### 1.5. Probar la validación de parámetros

Ahora probemos que la validación de parámetros funciona correctamente.

Primero, intenta ejecutar sin el parámetro requerido `input`:

```bash
nextflow run . --outdir test-results -profile docker
```

??? warning "Salida del comando"

    ```console
    ERROR ~ Validation of pipeline parameters failed!

    -- Check '.nextflow.log' file for details
    The following invalid input values have been detected:

    * Missing required parameter(s): input, batch
    ```

¡Perfecto! La validación detecta el parámetro requerido faltante antes de que se ejecute el pipeline.

Ahora intenta con un conjunto válido de parámetros:

```bash
nextflow run . --input assets/greetings.csv --outdir results --batch my-batch -profile test,docker
```

??? success "Salida del comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [peaceful_wozniak] DSL2 - revision: b9e9b3b8de

    executor >  local (8)
    [de/a1b2c3] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [4f/d5e6f7] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [8a/b9c0d1] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [e2/f3a4b5] CORE_HELLO:HELLO:COWPY (test)       | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

El pipeline debería ejecutarse exitosamente, y el parámetro `batch` ahora está validado.

### Conclusión

Has aprendido cómo usar la herramienta interactiva `nf-core pipelines schema build` para agregar parámetros a `nextflow_schema.json` y has visto la validación de parámetros en acción.
La interfaz web maneja toda la sintaxis JSON Schema por ti, facilitando la gestión de schemas de parámetros complejos sin edición manual de JSON propensa a errores.

### ¿Qué sigue?

Ahora que la validación de parámetros está funcionando, agreguemos validación para los contenidos del archivo de datos de entrada.

---

## 2. Validación de datos de entrada (schema_input.json)

Vamos a agregar validación para los contenidos de nuestro archivo CSV de entrada.
Mientras que la validación de parámetros verifica flags de línea de comandos, la validación de datos de entrada asegura que los datos dentro del archivo CSV estén estructurados correctamente.

### 2.1. Entender el formato de greetings.csv

Recordemos cómo se ve nuestra entrada:

```bash
cat assets/greetings.csv
```

```csv title="assets/greetings.csv"
Hello,en,87
Bonjour,fr,96
Holà,es,98
```

Este es un CSV simple con:

- Tres columnas (sin encabezado)
- En cada línea: un saludo, un idioma y un puntaje
- Las primeras dos columnas son cadenas de texto sin requisitos de formato especiales
- La tercera columna es un entero

Para nuestro pipeline, solo se requiere la primera columna.

### 2.2. Diseñar la estructura del schema

Para nuestro caso de uso, queremos:

1. Aceptar entrada CSV con al menos una columna
2. Tratar el primer elemento de cada fila como una cadena de saludo
3. Asegurar que los saludos no estén vacíos y no comiencen con espacios en blanco
4. Asegurar que el campo de idioma coincida con uno de los códigos de idioma soportados (en, fr, es, it, de)
5. Asegurar que el campo de puntaje sea un entero con un valor entre 0 y 100

Estructuraremos esto como un array de objetos, donde cada objeto tiene al menos un campo `greeting`.

### 2.3. Actualizar el archivo schema

La plantilla de pipeline nf-core incluye un `assets/schema_input.json` predeterminado diseñado para datos de secuenciación paired-end.
Necesitamos reemplazarlo con un schema más simple para nuestro caso de uso de saludos.

Abre `assets/schema_input.json` y reemplaza las secciones `properties` y `required`:

=== "Después"

    ```json title="assets/schema_input.json" linenums="1" hl_lines="10-25 27"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/core/hello/main/assets/schema_input.json",
        "title": "core/hello pipeline - params.input schema",
        "description": "Schema for the greetings file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "greeting": {
                    "type": "string",
                    "pattern": "^\\S.*$",
                    "errorMessage": "Greeting must be provided and cannot be empty or start with whitespace"
                },
                "language": {
                    "type": "string",
                    "enum": ["en", "fr", "es", "it", "de"],
                    "errorMessage": "Language must be one of: en, fr, es, it, de"
                },
                "score": {
                    "type": "integer",
                    "minimum": 0,
                    "maximum": 100,
                    "errorMessage": "Score must be an integer with a value between 0 and 100"
                }
            },
            "required": ["greeting"]
        }
    }
    ```

=== "Antes"

    ```json title="assets/schema_input.json" linenums="1" hl_lines="10-29 31"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/core/hello/main/assets/schema_input.json",
        "title": "core/hello pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Los cambios clave:

- **`description`**: Actualizada para mencionar "greetings file"
- **`properties`**: Reemplazadas `sample`, `fastq_1` y `fastq_2` con `greeting`, `language` y `score`
  - **`type:`** Forzar ya sea string (`greeting`, `language`) o integer (`score`)
  - **`pattern: "^\\S.*$"`**: El saludo debe comenzar con un carácter que no sea espacio en blanco (pero puede contener espacios después de eso)
  - **`"enum": ["en", "fr", "es", "it", "de"]`**: El código de idioma debe estar en el conjunto soportado
  - **`"minimum": 0` y `"maximum": 100`**: El valor del puntaje debe estar entre 0 y 100
  - **`errorMessage`**: Mensaje de error personalizado mostrado si la validación falla
- **`required`**: Cambiado de `["sample", "fastq_1"]` a `["greeting"]`

### 2.4. Agregar un encabezado al archivo greetings.csv

Cuando nf-schema lee un archivo CSV, espera que la primera fila contenga encabezados de columna que coincidan con los nombres de campo en el schema.

Para nuestro caso simple, necesitamos agregar una línea de encabezado a nuestro archivo de saludos:

=== "Después"

    ```csv title="assets/greetings.csv" linenums="1" hl_lines="1"
    greeting,language,score
    Hello,en,87
    Bonjour,fr,96
    Holà,es,98
    ```

=== "Antes"

    ```csv title="assets/greetings.csv" linenums="1"
    Hello,en,87
    Bonjour,fr,96
    Holà,es,98
    ```

Ahora el archivo CSV tiene una línea de encabezado que coincide con los nombres de campo en nuestro schema.

El paso final es implementar la validación en el código del pipeline usando `samplesheetToList`.

### 2.5. Implementar la validación en el pipeline

Ahora necesitamos reemplazar nuestro análisis simple de CSV con la función `samplesheetToList` de nf-schema, que validará y analizará la hoja de muestras.

La función `samplesheetToList`:

1. Lee la hoja de muestras de entrada (CSV, TSV, JSON o YAML)
2. La valida contra el schema JSON proporcionado
3. Devuelve una lista de Groovy donde cada entrada corresponde a una fila
4. Lanza mensajes de error útiles si la validación falla

Actualicemos el código de manejo de entrada:

Abre `subworkflows/local/utils_nfcore_hello_pipeline/main.nf` y localiza la sección donde creamos el canal de entrada (alrededor de la línea 80).

Necesitamos:

1. Usar la función `samplesheetToList` (ya importada en la plantilla)
2. Validar y analizar la entrada
3. Extraer solo las cadenas de saludo para nuestro workflow

Primero, nota que la función `samplesheetToList` ya está importada en la parte superior del archivo (la plantilla nf-core incluye esto por defecto):

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="1" hl_lines="13"
//
// Subworkflow with functionality specific to the core/hello pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
```

Ahora actualiza el código de creación del canal:

=== "Después"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="80" hl_lines="4"
        //
        // Create channel from input file provided through params.input
        //
        ch_samplesheet = channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="80" hl_lines="4 5"
        //
        // Create channel from input file provided through params.input
        //
        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Desglosemos lo que cambió:

1. **`samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")`**: Valida el archivo de entrada contra nuestro schema y devuelve una lista
2. **`Channel.fromList(...)`**: Convierte la lista en un canal de Nextflow

Esto completa la implementación de validación de datos de entrada usando `samplesheetToList` y schemas JSON.

Ahora que hemos configurado el schema de datos de entrada, podemos eliminar la configuración de ignorar temporal que agregamos anteriormente.

### 2.6. Volver a habilitar la validación de entrada

Abre `nextflow.config` y elimina la línea `ignoreParams` del bloque `validation`:

=== "Después"

    ```groovy title="nextflow.config" linenums="246"
    validation {
        defaultIgnoreParams = ["genomes"]
        monochromeLogs = params.monochrome_logs
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" hl_lines="3" linenums="246"
    validation {
        defaultIgnoreParams = ["genomes"]
        ignoreParams = ['input']
        monochromeLogs = params.monochrome_logs
    }
    ```

Ahora nf-schema validará tanto los tipos de parámetros COMO los contenidos del archivo de entrada.

### 2.7. Probar la validación de entrada

Verifiquemos que nuestra validación funciona probando tanto entradas válidas como inválidas.

#### 2.7.1. Probar con entrada válida

Primero, confirma que el pipeline se ejecuta exitosamente con entrada válida.
¡Nota que ya no necesitamos `--validate_params false` ya que la validación está funcionando!

```bash
nextflow run . --outdir core-hello-results -profile test,docker
```

??? success "Salida del comando"

    ```console
    ------------------------------------------------------
    WARN: The following invalid input values have been detected:

    * --character: tux


    executor >  local (8)
    [c1/39f64a] CORE_HELLO:HELLO:sayHello (1)       | 3 of 3 ✔
    [44/c3fb82] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [62/80fab2] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [e1/4db4fd] CORE_HELLO:HELLO:COWPY (test)       | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

¡Excelente! El pipeline se ejecuta exitosamente y la validación pasa silenciosamente.
La advertencia sobre `--character` es solo informativa ya que no está definida en el schema.
Si quieres, ¡usa lo que has aprendido para agregar validación para ese parámetro también!

#### 2.7.2. Probar con entrada inválida

Pasar la validación siempre es una buena sensación, pero asegurémonos de que la validación realmente detectará errores.

Para crear un archivo de prueba con un nombre de columna inválido, comienza haciendo una copia del archivo `greetings.csv`:

```bash
cp assets/greetings.csv assets/invalid_greetings.csv
```

Ahora abre el archivo y cambia el nombre de la primera columna, en la línea de encabezado, de `greeting` a `message`:

=== "Después"

    ```csv title="tmp_invalid_greetings.csv" hl_lines="1" linenums="1"
    message,language,score
    Hello,en,87
    Bonjour,fr,96
    Holà,es,98
    ```

=== "Antes"

    ```csv title="tmp_invalid_greetings.csv" hl_lines="1" linenums="1"
    greeting,language,score
    Hello,en,87
    Bonjour,fr,96
    Holà,es,98
    ```

Esto no coincide con nuestro schema, por lo que la validación debería lanzar un error.

Intenta ejecutar el pipeline con esta entrada inválida:

```bash
nextflow run . --input assets/invalid_greetings.csv --outdir test-results -profile docker
```

??? failure "Salida del comando"

    ```console
    N E X T F L O W   ~  version 24.10.4

    Launching `./main.nf` [trusting_ochoa] DSL2 - revision: b9e9b3b8de

    Input/output options
      input              : assets/invalid_greetings.csv
      outdir             : test-results

    Generic options
      trace_report_suffix: 2025-01-27_03-16-04

    Core Nextflow options
      runName            : trusting_ochoa
      containerEngine    : docker
      launchDir          : /workspace/hello-nf-core
      workDir            : /workspace/hello-nf-core/work
      projectDir         : /workspace/hello-nf-core
      userName           : user
      profile            : docker
      configFiles        : /workspace/hello-nf-core/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    ERROR ~ Validation of pipeline parameters failed!

     -- Check '.nextflow.log' file for details
    The following invalid input values have been detected:

    * Missing required parameter(s): batch
    * --input (assets/invalid_greetings.csv): Validation of file failed:
        -> Entry 1: Missing required field(s): greeting
        -> Entry 2: Missing required field(s): greeting
        -> Entry 3: Missing required field(s): greeting

     -- Check script 'subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
    ```

¡Perfecto! La validación detectó el error y proporcionó un mensaje de error claro y útil señalando:

- Qué archivo falló la validación
- Qué entrada (fila 1, la primera fila de datos) tiene el problema
- Cuál es el problema específico (campo requerido faltante `greeting`)

La validación del schema asegura que los archivos de entrada tengan la estructura correcta antes de que se ejecute el pipeline, ahorrando tiempo y previniendo errores confusos más adelante en la ejecución.

Si quieres practicar esto, siéntete libre de crear otros archivos de entrada de saludos que violen el schema de otras formas divertidas.

### Conclusión

Has implementado y probado tanto la validación de parámetros como la validación de datos de entrada. Tu pipeline ahora valida las entradas antes de la ejecución, proporcionando retroalimentación rápida y mensajes de error claros.

!!! tip "Lectura adicional"

    Para aprender más sobre características y patrones de validación avanzados, consulta la [documentación de nf-schema](https://nextflow-io.github.io/nf-schema/latest/). El comando `nf-core pipelines schema build` proporciona una GUI interactiva para gestionar schemas complejos.

### ¿Qué sigue?

¡Has completado las cinco partes del curso de capacitación Hello nf-core!

Continúa al [Resumen](summary.md) para reflexionar sobre lo que has construido y aprendido.
