# Parte 1: Descripción general del método y pruebas manuales

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![Descripción general del pipeline](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Métodos

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Antes de comenzar a escribir cualquier código de workflow, vamos a probar los comandos manualmente con algunos datos de prueba.

### Conjunto de datos

Proporcionamos los siguientes datos y recursos relacionados:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Software

Las herramientas principales involucradas son [{TOOL_A}]({TOOL_A_URL}) y [{TOOL_B}]({TOOL_B_URL}).

Estas herramientas no están instaladas en el entorno de GitHub Codespaces, por lo que las usaremos a través de contenedores (consulte [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

    Asegúrese de estar en el directorio `nf4-science/{DOMAIN_DIR}` para que la última parte de la ruta que se muestra al escribir `pwd` sea `{DOMAIN_DIR}`.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

En esta sección probamos los comandos que conforman el enfoque de procesamiento de muestra única.
Estos son los comandos que envolveremos en un workflow de Nextflow en la Parte 2 de este curso.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Comenzamos probando los comandos en una sola muestra.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Descargar el contenedor

Ejecute el comando `docker pull` para descargar la imagen del contenedor:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Salida del comando"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Iniciar el contenedor de forma interactiva

Inicie el contenedor y monte el directorio `data` para que las herramientas puedan acceder a los archivos de entrada:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

Su prompt cambia para indicar que está dentro del contenedor.

#### 1.1.3. Ejecutar el comando

```bash
{TOOL_A_COMMAND}
```

??? success "Salida del comando"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Descargar el contenedor

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Salida del comando"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Iniciar el contenedor de forma interactiva

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Ejecutar el comando

```bash
{TOOL_B_COMMAND}
```

??? success "Salida del comando"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Salir del contenedor

Para salir del contenedor, escriba `exit`.

```bash
exit
```

Su prompt debería volver a la normalidad.
Esto concluye la prueba de procesamiento de muestra única.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

En esta sección probamos los comandos adicionales necesarios para el procesamiento de múltiples muestras.
Estos son los comandos que envolveremos en un workflow de Nextflow en la Parte 3 de este curso.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Conclusión

Ahora sabe cómo probar los comandos de {TOOL_A} y {TOOL_B} en sus respectivos contenedores, incluyendo cómo {MULTI_SAMPLE_SUMMARY}.

### ¿Qué sigue?

Aprenda cómo envolver esos mismos comandos en workflows que usan contenedores para ejecutar el trabajo.
