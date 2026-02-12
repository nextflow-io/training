# Parte 3: Agregación de múltiples muestras

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En la Parte 2, construiste un pipeline de procesamiento por muestra que manejaba cada muestra de forma independiente.
Ahora vamos a extenderlo para implementar {AGGREGATION_METHOD} de múltiples muestras, como se cubrió en la [Parte 1](01_method.md).

## Asignación

En esta parte del curso, vamos a extender el workflow para hacer lo siguiente:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Esta parte se construye directamente sobre el workflow producido en la Parte 2.

??? info "Cómo comenzar desde esta sección"

    Esta sección del curso asume que has completado la [Parte 2: Procesamiento de una sola muestra](./02_single_sample.md) y tienes un pipeline `{DOMAIN_DIR}.nf` funcional.

    Si no completaste la Parte 2 o quieres comenzar de nuevo para esta parte, puedes usar la solución de la Parte 2 como punto de partida.
    Ejecuta estos comandos desde dentro del directorio `nf4-science/{DOMAIN_DIR}/`:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Esto te proporciona un workflow completo de procesamiento de una sola muestra.
    Puedes verificar que se ejecuta correctamente ejecutando el siguiente comando:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Plan de la lección

Hemos dividido esto en dos pasos:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Esto cubre la actualización de comandos y salidas de procesos.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Esto introduce el operador `collect()` {AND_OTHER_CONCEPTS}.

!!! note "Nota"

     Asegúrate de estar en el directorio de trabajo correcto:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Recuerda el comando modificado de la [Parte 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Después"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Antes"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Ejecuta el workflow para verificar la modificación

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Salida del comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusión

Sabes cómo modificar comandos y salidas de procesos para adaptar el comportamiento del workflow.

### ¿Qué sigue?

Agrega el paso de agregación de múltiples muestras.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Escribe el módulo de agregación

{MODULE_INSTRUCTIONS}

### 2.2. Recopila las salidas por muestra y envíalas al proceso de agregación

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Ejecuta el workflow completo

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Salida del comando"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Conclusión

Tienes un pipeline completo que procesa muestras individualmente y agrega resultados de todas las muestras.
Sabes cómo usar operadores de canal como `collect()` para agregar salidas por muestra para análisis de múltiples muestras.

### ¿Qué sigue?

¡Felicidades por completar este curso! Dirígete al [resumen del curso](next_steps.md) para revisar lo que has aprendido y explorar los próximos pasos.
