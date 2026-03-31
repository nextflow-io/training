# Parte 1: Conceptos Básicos de Plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En esta sección, aprenderá cómo los plugins extienden Nextflow y luego probará tres plugins diferentes para verlos en acción.

---

## 1. Cómo funcionan los plugins

Los plugins extienden Nextflow a través de varios tipos de extensiones:

| Tipo de extensión    | Qué hace                                                  | Ejemplo                      |
| -------------------- | --------------------------------------------------------- | ---------------------------- |
| Funciones            | Agrega funciones personalizadas invocables desde workflows | `samplesheetToList()`        |
| Monitores de workflow | Responden a eventos como la finalización de tareas        | Logging personalizado, alertas de Slack |
| Executors            | Agregan backends de ejecución de tareas                   | AWS Batch, Kubernetes        |
| Filesystems          | Agregan backends de almacenamiento                        | S3, Azure Blob               |

Las funciones y los monitores de workflow (llamados "trace observers" en la API de Nextflow) son los tipos más comunes para los autores de plugins.
Los executors y filesystems son creados típicamente por proveedores de plataformas.

Los siguientes ejercicios muestran plugins de funciones y un plugin observer, para que pueda ver ambos tipos en acción.

---

## 2. Usar plugins de funciones

Los plugins de funciones agregan funciones invocables que se importan en los workflows.
Probará dos: nf-hello (un ejemplo simple) y nf-schema (un plugin real ampliamente utilizado).
Ambos ejercicios modifican el mismo pipeline `hello.nf`, para que pueda ver cómo los plugins mejoran un workflow existente.

### 2.1. nf-hello: reemplazar código escrito a mano

El plugin [nf-hello](https://github.com/nextflow-io/nf-hello) proporciona una función `randomString` que genera cadenas de texto aleatorias.
El pipeline ya define su propia versión en línea de esta función, que reemplazará con la del plugin.

#### 2.1.1. Ver el punto de partida

Examine el pipeline:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Genera una cadena alfanumérica aleatoria
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

El pipeline define su propia función `randomString` en línea y luego la usa para agregar un ID aleatorio a cada saludo.

Ejecútelo:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

El orden de la salida y las cadenas aleatorias serán diferentes, y si ejecuta el script nuevamente obtendrá un conjunto distinto de saludos aleatorios.

#### 2.1.2. Configurar el plugin

Reemplace la función en línea con una del plugin. Agregue esto a su `nextflow.config`:

```groovy title="nextflow.config"
// Configuración para los ejercicios de desarrollo de plugins
plugins {
    id 'nf-hello@0.5.0'
}
```

Los plugins se declaran en `nextflow.config` usando el bloque `plugins {}`.
Nextflow los descarga automáticamente desde el [Nextflow Plugin Registry](https://registry.nextflow.io/), un repositorio central de plugins oficiales y de la comunidad.

#### 2.1.3. Usar la función del plugin

Reemplace la función `randomString` en línea con la versión del plugin:

=== "Después"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Antes"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Genera una cadena alfanumérica aleatoria
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

La sentencia `include` importa `randomString` desde una librería probada, testeada y mantenida por un grupo más amplio de colaboradores que pueden detectar y corregir errores.
En lugar de que cada pipeline mantenga su propia copia de la función, todos los pipelines que usan el plugin obtienen la misma implementación verificada.
Esto reduce el código duplicado y la carga de mantenimiento que conlleva.
La sintaxis `#!groovy include { function } from 'plugin/plugin-id'` es el mismo `include` que se usa para los módulos de Nextflow, con el prefijo `plugin/`.
Puede ver el [código fuente de `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) en el repositorio nf-hello en GitHub.

#### 2.1.4. Ejecutarlo

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Sus cadenas aleatorias serán diferentes.)

La salida sigue teniendo sufijos aleatorios, pero ahora `randomString` proviene del plugin nf-hello en lugar del código en línea.
Los mensajes "Pipeline is starting!" y "Pipeline complete!" son nuevos.
Provienen del componente observer del plugin, que explorará en la Parte 5.

Nextflow descarga los plugins automáticamente la primera vez que se usan, por lo que cualquier pipeline que declare `nf-hello@0.5.0` obtiene exactamente la misma función `randomString` testeada sin necesidad de copiar código entre proyectos.

Ahora ha visto los tres pasos para usar un plugin de funciones: declararlo en `nextflow.config`, importar la función con `include` y llamarla en su workflow.
El siguiente ejercicio aplica estos mismos pasos a un plugin del mundo real.

### 2.2. nf-schema: análisis de CSV con validación

El plugin [nf-schema](https://github.com/nextflow-io/nf-schema) es uno de los plugins de Nextflow más ampliamente utilizados.
Proporciona `samplesheetToList`, una función que analiza archivos CSV/TSV usando un esquema JSON que define las columnas y tipos esperados.

El pipeline actualmente lee `greetings.csv` usando `splitCsv` y un `map` manual, pero nf-schema puede reemplazar esto con un análisis validado basado en esquemas.
Ya se proporciona un archivo de esquema JSON (`greetings_schema.json`) en el directorio del ejercicio.

??? info "¿Qué es un esquema?"

    Un esquema es una descripción formal de cómo deben verse los datos válidos.
    Define cosas como qué columnas se esperan, qué tipo debe tener cada valor (string, número, etc.) y qué campos son obligatorios.

    Piénselo como un contrato: si los datos de entrada no coinciden con el esquema, la herramienta puede detectar el problema de forma temprana en lugar de dejar que cause errores confusos más adelante en el pipeline.

#### 2.2.1. Examinar el esquema

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

El esquema define dos columnas (`greeting` y `language`) y marca `greeting` como obligatoria.
Si alguien pasa un CSV sin la columna `greeting`, nf-schema detecta el error antes de que el pipeline se ejecute.

#### 2.2.2. Agregar nf-schema a la configuración

Actualice `nextflow.config` para incluir ambos plugins:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. Actualizar hello.nf para usar samplesheetToList

Reemplace la entrada con `splitCsv` por `samplesheetToList`:

=== "Después"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Antes"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

El código personalizado de análisis con `splitCsv` y `map` se reemplaza con `samplesheetToList`, una función probada y testeada que también valida el samplesheet contra el esquema antes de que el pipeline se ejecute.
Esto reduce la carga de mantenimiento de la lógica de análisis escrita a mano, al mismo tiempo que mejora la experiencia para los usuarios del pipeline, quienes reciben mensajes de error claros cuando su entrada no coincide con el formato esperado.
Cada fila se convierte en una lista de valores en orden de columna, por lo que `row[0]` es el saludo y `row[1]` es el idioma.

#### 2.2.4. Ejecutarlo

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Sus cadenas aleatorias serán diferentes.)

La salida es la misma, pero ahora el esquema valida la estructura del CSV antes de que el pipeline se ejecute.
En pipelines reales con sample sheets complejos y muchas columnas, este tipo de validación previene errores que `splitCsv` + `map` manual pasaría por alto.

#### 2.2.5. Ver la validación en acción

Para ver qué detecta la validación de esquemas, intente introducir errores en `greetings.csv`.

Cambie el nombre de la columna obligatoria `greeting` a `message`:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Ejecute el pipeline:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

El pipeline se niega a ejecutarse porque el esquema requiere una columna `greeting` y no puede encontrarla.

Ahora restaure la columna obligatoria pero cambie el nombre de la columna opcional `language` a `lang`:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Esta vez el pipeline se ejecuta, pero muestra una advertencia:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Las columnas obligatorias generan errores graves; las columnas opcionales generan advertencias.
Este es el tipo de retroalimentación temprana que ahorra tiempo de depuración en pipelines reales con docenas de columnas.

#### 2.2.6. Configurar el comportamiento de validación

La advertencia sobre `lang` es útil, pero puede controlar su severidad a través de la configuración.
Los plugins pueden incluir sus propios ámbitos de configuración que controlan su comportamiento.
El plugin nf-schema incluye el ámbito de configuración `validation`; al modificar los ajustes aquí puede cambiar cómo se comporta nf-schema.

Agregue un bloque `validation` a `nextflow.config` para que los encabezados no reconocidos generen un error en lugar de una advertencia:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Ejecute el pipeline nuevamente con la misma columna `lang` todavía en su lugar:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

El pipeline ahora falla en lugar de advertir.
El código del pipeline no cambió; solo lo hizo la configuración.

Restaure `greetings.csv` a su estado original y elimine el bloque `validation` antes de continuar:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

Tanto nf-hello como nf-schema son plugins de funciones: proporcionan funciones que se importan con `include` y se llaman en el código del workflow.
El siguiente ejercicio muestra un tipo diferente de plugin que funciona sin ninguna sentencia `include`.

---

## 3. Usar un plugin observer: nf-co2footprint

No todos los plugins proporcionan funciones para importar.
El plugin [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) usa un **trace observer** para monitorear el uso de recursos de su pipeline y estimar su huella de carbono.
No necesita cambiar ningún código del pipeline; simplemente agréguelo a la configuración.

### 3.1. Agregar nf-co2footprint a la configuración

Actualice `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Ejecutar el pipeline

```bash
nextflow run hello.nf
```

El plugin produce varios mensajes INFO y WARN durante la ejecución.
Estos son normales para un ejemplo pequeño que se ejecuta en una máquina local:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Las advertencias sobre la zona, el executor, el modelo de CPU y la memoria aparecen porque el plugin no puede detectar todos los detalles del hardware en un entorno de capacitación local.
En un entorno de producción (por ejemplo, un clúster HPC o la nube), estos valores estarían disponibles y las estimaciones serían más precisas.

Al final, busque una línea como:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Sus números serán diferentes.)

### 3.3. Ver el reporte

El plugin genera archivos de salida en su directorio de trabajo:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Examine el resumen:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(Sus números serán diferentes.)

La primera sección muestra las cifras brutas de energía y emisiones.
La sección "Which equals" pone esos números en perspectiva convirtiéndolos a equivalentes familiares.
El resumen también incluye una sección que lista las opciones de configuración del plugin y una cita al artículo de investigación [Green Algorithms](https://doi.org/10.1002/advs.202100707) en el que se basa el método de cálculo.

### 3.4. Configurar el plugin

La advertencia "Target zone null" de la sección 3.2 apareció porque el plugin no tenía ninguna ubicación configurada.
El plugin nf-co2footprint define un ámbito de configuración `co2footprint` donde puede establecer su ubicación geográfica.

Agregue un bloque `co2footprint` a `nextflow.config`:

=== "Después"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Consejo"

    Use el código de su propio país si lo prefiere (por ejemplo, `'US'`, `'DE'`, `'FR'`).

Ejecute el pipeline:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

La advertencia de zona desapareció.
El plugin ahora usa la intensidad de carbono específica de GB (163.92 gCO₂eq/kWh) en lugar del valor global de respaldo (480.0 gCO₂eq/kWh).

!!! note "Nota"

    También puede ver un mensaje `WARN: Unrecognized config option 'co2footprint.location'`.
    Esto es cosmético y puede ignorarse sin problema; el plugin sigue leyendo el valor correctamente.

En la Parte 6, creará un ámbito de configuración para su propio plugin.

Este plugin funciona completamente a través del mecanismo observer, conectándose a los eventos del ciclo de vida del workflow para recopilar métricas de recursos y generar su reporte cuando el pipeline finaliza.

Ahora ha probado plugins de funciones (importados con `include`) y un plugin observer (activado solo a través de la configuración).
Estos son los dos tipos de extensión más comunes, pero como muestra la tabla en la sección 1, los plugins también pueden agregar executors y filesystems.

---

## 4. Descubrir plugins

El [Nextflow Plugin Registry](https://registry.nextflow.io/) es el centro principal para encontrar plugins disponibles.

![La página del plugin nf-hello en registry.nextflow.io](img/plugin-registry-nf-hello.png)

Cada página de plugin muestra su descripción, versiones disponibles, instrucciones de instalación y enlaces a la documentación.

---

## 5. Prepararse para el desarrollo de plugins

Las siguientes secciones (Partes 2-6) usan un archivo de pipeline separado, `greet.nf`, que depende de nf-schema pero no de nf-hello ni de nf-co2footprint.

Actualice `nextflow.config` para mantener solo nf-schema:

```groovy title="nextflow.config"
// Configuración para los ejercicios de desarrollo de plugins
plugins {
    id 'nf-schema@2.6.1'
}
```

Elimine los archivos de salida de co2footprint:

```bash
rm -f co2footprint_*
```

El archivo `hello.nf` conserva el trabajo de la Parte 1 como referencia; a partir de ahora, trabajará con `greet.nf`.

---

## Conclusión

Utilizó tres plugins diferentes:

- **nf-hello**: Un plugin de funciones que proporciona `randomString`, importado con `include`
- **nf-schema**: Un plugin de funciones que proporciona `samplesheetToList` para el análisis de CSV con validación de esquemas
- **nf-co2footprint**: Un plugin observer que monitorea el uso de recursos automáticamente, sin necesidad de `include`

Patrones clave:

- Los plugins se declaran en `nextflow.config` con `#!groovy plugins { id 'plugin-name@version' }`
- Los plugins de funciones requieren `#!groovy include { function } from 'plugin/plugin-id'`
- Los plugins observer funcionan automáticamente una vez declarados en la configuración
- Los plugins pueden definir ámbitos de configuración (por ejemplo, `#!groovy validation {}`, `#!groovy co2footprint {}`) para personalizar el comportamiento
- El [Nextflow Plugin Registry](https://registry.nextflow.io/) lista los plugins disponibles

---

## ¿Qué sigue?

Las siguientes secciones le muestran cómo construir su propio plugin.
Si no está interesado en el desarrollo de plugins, puede detenerse aquí o saltar directamente al [Resumen](summary.md).

[Continuar a la Parte 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
