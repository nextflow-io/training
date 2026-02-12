# Part 1: Visió general del mètode i proves manual

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![Visió general del pipeline](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Mètodes

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Abans de començar a escriure cap codi de workflow, provarem les comandes manualment amb algunes dades de prova.

### Conjunt de dades

Proporcionem les dades i recursos relacionats següents:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Programari

Les eines principals implicades són [{TOOL_A}]({TOOL_A_URL}) i [{TOOL_B}]({TOOL_B_URL}).

Aquestes eines no estan instal·lades a l'entorn de GitHub Codespaces, així que les utilitzarem mitjançant contenidors (vegeu [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Nota"

    Assegureu-vos que esteu al directori `nf4-science/{DOMAIN_DIR}` de manera que l'última part del camí que es mostra quan escriviu `pwd` sigui `{DOMAIN_DIR}`.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

En aquesta secció provem les comandes que componen l'enfocament de processament d'una sola mostra.
Aquestes són les comandes que encapsularem en un workflow de Nextflow a la Part 2 d'aquest curs.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Comencem provant les comandes amb només una mostra.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Descarregueu el contenidor

Executeu la comanda `docker pull` per descarregar la imatge del contenidor:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Inicieu el contenidor de manera interactiva

Inicieu el contenidor i munteu el directori `data` perquè les eines puguin accedir als fitxers d'entrada:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

El vostre indicador canvia per indicar que esteu dins del contenidor.

#### 1.1.3. Executeu la comanda

```bash
{TOOL_A_COMMAND}
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Sortiu del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre indicador hauria de tornar a la normalitat.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Descarregueu el contenidor

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Inicieu el contenidor de manera interactiva

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Executeu la comanda

```bash
{TOOL_B_COMMAND}
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Sortiu del contenidor

Per sortir del contenidor, escriviu `exit`.

```bash
exit
```

El vostre indicador hauria de tornar a la normalitat.
Això conclou la prova de processament d'una sola mostra.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

En aquesta secció provem les comandes addicionals necessàries per al processament de múltiples mostres.
Aquestes són les comandes que encapsularem en un workflow de Nextflow a la Part 3 d'aquest curs.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Conclusió

Sabeu com provar les comandes de {TOOL_A} i {TOOL_B} als seus respectius contenidors, incloent com {MULTI_SAMPLE_SUMMARY}.

### Què segueix?

Apreneu com encapsular aquestes mateixes comandes en workflows que utilitzen contenidors per executar el treball.
