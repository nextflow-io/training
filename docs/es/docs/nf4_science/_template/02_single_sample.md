# Parte 2: Procesamiento de una sola muestra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la Parte 1, probaste los comandos de {TOOL_A} y {TOOL_B} manualmente en sus respectivos contenedores.
Ahora vamos a envolver esos mismos comandos en un workflow de Nextflow.

## Asignación

En esta parte del curso, vamos a desarrollar un workflow que hace lo siguiente:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Esto replica los pasos de la Parte 1, donde ejecutaste estos comandos manualmente en sus contenedores.

Como punto de partida, te proporcionamos un archivo de workflow, `{DOMAIN_DIR}.nf`, que describe las partes principales del workflow, así como dos archivos de módulo, {TOOL_A_MODULE}.nf y {TOOL_B_MODULE}.nf, que describen la estructura de los módulos.
Estos archivos no son funcionales; su propósito es solo servir como esqueletos para que los completes con las partes interesantes del código.

## Plan de la lección

Para hacer el proceso de desarrollo más educativo, lo hemos dividido en {N} pasos:

1. **Escribe un workflow de una sola etapa que ejecute {TOOL_A_ACTION}.**
   Esto cubre la creación de un módulo, su importación y su llamada en un workflow.
2. **Agrega un segundo proceso para ejecutar {TOOL_B_ACTION}.**
   Esto introduce el encadenamiento de salidas de procesos a entradas y el manejo de archivos accesorios.
3. **Adapta el workflow para ejecutarse en un lote de muestras.**
   Esto cubre la ejecución paralela e introduce tuplas para mantener juntos los archivos asociados.
4. **Haz que el workflow acepte una hoja de muestras que contenga un lote de archivos de entrada.**
   Esto demuestra un patrón común para proporcionar entradas en lote.

Cada paso se enfoca en un aspecto específico del desarrollo de workflows.

---

## 1. Escribe un workflow de una sola etapa que ejecute {TOOL_A_ACTION}

Este primer paso se enfoca en lo básico: cargar {PRIMARY_INPUT_TYPE} y {TOOL_A_OUTPUT_DESCRIPTION}.

Recuerda el comando `{TOOL_A_COMMAND_NAME}` de la [Parte 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

El comando toma {INPUT_DESCRIPTION} y produce {OUTPUT_DESCRIPTION}.
El URI del contenedor era `{TOOL_A_CONTAINER_URI}`.

Vamos a tomar esta información y envolverla en Nextflow en tres etapas:

1. Configurar la entrada
2. Escribir el proceso y llamarlo en el workflow
3. Configurar el manejo de la salida

### 1.1. Configura la entrada

Necesitamos declarar un parámetro de entrada, crear un perfil de prueba para proporcionar un valor predeterminado conveniente y crear un canal de entrada.

#### 1.1.1. Agrega una declaración de parámetro de entrada

En el archivo principal del workflow `{DOMAIN_DIR}.nf`, bajo la sección `Pipeline parameters`, declara un parámetro CLI llamado `{PRIMARY_PARAM_NAME}`.

=== "Después"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Antes"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Eso configura el parámetro CLI, pero no queremos escribir la ruta del archivo cada vez que ejecutemos el workflow durante el desarrollo.
Hay múltiples opciones para proporcionar un valor predeterminado; aquí usamos un perfil de prueba.

#### 1.1.2. Crea un perfil de prueba con un valor predeterminado en `nextflow.config`

Un perfil de prueba proporciona valores predeterminados convenientes para probar un workflow sin especificar entradas en la línea de comandos.
Esta es una convención común en el ecosistema de Nextflow (consulta [Hello Config](../../hello_nextflow/06_hello_config.md) para más detalles).

Agrega un bloque `profiles` a `nextflow.config` con un perfil `test` que establezca el parámetro `{PRIMARY_PARAM_NAME}` a uno de los {PRIMARY_INPUT_TYPE}s de prueba.

=== "Después"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aquí, estamos usando `${projectDir}`, una variable integrada de Nextflow que apunta al directorio donde se encuentra el script del workflow.
Esto facilita la referencia a archivos de datos y otros recursos sin codificar rutas absolutas.

#### 1.1.3. Configura el canal de entrada

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Escribe el módulo de {TOOL_A_NAME}

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Importa y llama al módulo en el workflow

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Ejecuta el workflow

En este punto, tenemos un workflow de un paso que debería ser completamente funcional.

Podemos ejecutarlo con `-profile test` para usar el valor predeterminado configurado en el perfil de prueba y evitar tener que escribir la ruta en la línea de comandos.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Salida del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusión

Sabes cómo crear un módulo que contiene un proceso, importarlo en un workflow, llamarlo con un canal de entrada y publicar los resultados.

### ¿Qué sigue?

Agrega un segundo proceso para encadenar pasos de análisis adicionales.

---

## 2. Agrega un segundo proceso para ejecutar {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Recuerda el comando `{TOOL_B_COMMAND_NAME}` de la [Parte 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Escribe el módulo de {TOOL_B_NAME}

{MODULE_INSTRUCTIONS}

### 2.2. Importa y llama al módulo en el workflow

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Ejecuta el workflow

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Salida del comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusión

Sabes cómo encadenar salidas de procesos a entradas y manejar archivos accesorios en el workflow.

### ¿Qué sigue?

Escala el workflow para procesar múltiples muestras en paralelo.

---

## 3. Adapta el workflow para ejecutarse en un lote de muestras

Hasta ahora, el workflow procesa una sola muestra.
Para manejar múltiples muestras, necesitamos modificar cómo se proporcionan las entradas y aprovechar el paradigma de flujo de datos de Nextflow para paralelizar la ejecución.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Ejecuta el workflow en múltiples muestras

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Salida del comando"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Conclusión

Sabes cómo aprovechar el paradigma de flujo de datos de Nextflow para paralelizar el procesamiento por muestra en múltiples muestras de entrada.

### ¿Qué sigue?

En la [Parte 3](03_multi_sample.md), agregarás agregación de múltiples muestras para combinar resultados de todas las muestras.
