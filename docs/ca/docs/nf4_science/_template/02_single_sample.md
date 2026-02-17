# Part 2: Processament d'una sola mostra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la Part 1, vau provar les comandes de {TOOL_A} i {TOOL_B} manualment als seus respectius contenidors.
Ara encapsularem aquestes mateixes comandes en un workflow de Nextflow.

## Assignació

En aquesta part del curs, desenvoluparem un workflow que fa el següent:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Això replica els passos de la Part 1, on vau executar aquestes comandes manualment als seus contenidors.

Com a punt de partida, us proporcionem un fitxer de workflow, `{DOMAIN_DIR}.nf`, que esbossa les parts principals del workflow, així com dos fitxers de mòdul, {TOOL_A_MODULE}.nf i {TOOL_B_MODULE}.nf, que esbossen l'estructura dels mòduls.
Aquests fitxers no són funcionals; el seu propòsit és només servir com a plantilles perquè les ompliu amb les parts interessants del codi.

## Pla de la lliçó

Per tal de fer el procés de desenvolupament més educatiu, ho hem dividit en {N} passos:

1. **Escriure un workflow d'una sola etapa que executi {TOOL_A_ACTION}.**
   Això cobreix la creació d'un mòdul, la seva importació i la seva crida en un workflow.
2. **Afegir un segon procés per executar {TOOL_B_ACTION}.**
   Això introdueix l'encadenament de sortides de processos a entrades i la gestió de fitxers accessoris.
3. **Adaptar el workflow per executar-se en un lot de mostres.**
   Això cobreix l'execució paral·lela i introdueix les tuples per mantenir els fitxers associats junts.
4. **Fer que el workflow accepti un full de mostres que contingui un lot de fitxers d'entrada.**
   Això demostra un patró comú per proporcionar entrades en bloc.

Cada pas se centra en un aspecte específic del desenvolupament de workflows.

---

## 1. Escriure un workflow d'una sola etapa que executi {TOOL_A_ACTION}

Aquest primer pas se centra en els conceptes bàsics: carregar {PRIMARY_INPUT_TYPE} i {TOOL_A_OUTPUT_DESCRIPTION}.

Recordeu la comanda `{TOOL_A_COMMAND_NAME}` de la [Part 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

La comanda pren {INPUT_DESCRIPTION} i produeix {OUTPUT_DESCRIPTION}.
L'URI del contenidor era `{TOOL_A_CONTAINER_URI}`.

Encapsularem aquesta informació en Nextflow en tres etapes:

1. Configurar l'entrada
2. Escriure el procés i cridar-lo al workflow
3. Configurar la gestió de la sortida

### 1.1. Configurar l'entrada

Hem de declarar un paràmetre d'entrada, crear un perfil de prova per proporcionar un valor per defecte convenient i crear un canal d'entrada.

#### 1.1.1. Afegir una declaració de paràmetre d'entrada

Al fitxer principal del workflow `{DOMAIN_DIR}.nf`, sota la secció `Pipeline parameters`, declareu un paràmetre CLI anomenat `{PRIMARY_PARAM_NAME}`.

=== "Després"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Abans"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Això configura el paràmetre CLI, però no volem escriure el camí del fitxer cada vegada que executem el workflow durant el desenvolupament.
Hi ha múltiples opcions per proporcionar un valor per defecte; aquí utilitzem un perfil de prova.

#### 1.1.2. Crear un perfil de prova amb un valor per defecte a `nextflow.config`

Un perfil de prova proporciona valors per defecte convenients per provar un workflow sense especificar entrades a la línia de comandes.
Aquesta és una convenció comuna a l'ecosistema Nextflow (vegeu [Hello Config](../../hello_nextflow/06_hello_config.md) per a més detalls).

Afegiu un bloc `profiles` a `nextflow.config` amb un perfil `test` que estableixi el paràmetre `{PRIMARY_PARAM_NAME}` a un dels {PRIMARY_INPUT_TYPE}s de prova.

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aquí, estem utilitzant `${projectDir}`, una variable integrada de Nextflow que apunta al directori on es troba l'script del workflow.
Això facilita la referència a fitxers de dades i altres recursos sense codificar camins absoluts.

#### 1.1.3. Configurar el canal d'entrada

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Escriure el mòdul {TOOL_A_NAME}

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Importar i cridar el mòdul al workflow

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Executar el workflow

En aquest punt, tenim un workflow d'un sol pas que hauria de ser completament funcional.

Podem executar-lo amb `-profile test` per utilitzar el valor per defecte configurat al perfil de prova i evitar haver d'escriure el camí a la línia de comandes.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusió

Sabeu com crear un mòdul que contingui un procés, importar-lo a un workflow, cridar-lo amb un canal d'entrada i publicar els resultats.

### Què segueix?

Afegir un segon procés per encadenar passos d'anàlisi addicionals.

---

## 2. Afegir un segon procés per executar {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Recordeu la comanda `{TOOL_B_COMMAND_NAME}` de la [Part 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Escriure el mòdul {TOOL_B_NAME}

{MODULE_INSTRUCTIONS}

### 2.2. Importar i cridar el mòdul al workflow

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Executar el workflow

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Conclusió

Sabeu com encadenar sortides de processos a entrades i gestionar fitxers accessoris al workflow.

### Què segueix?

Escalar el workflow per processar múltiples mostres en paral·lel.

---

## 3. Adaptar el workflow per executar-se en un lot de mostres

Fins ara, el workflow processa una sola mostra.
Per gestionar múltiples mostres, hem de modificar com es proporcionen les entrades i aprofitar el paradigma de flux de dades de Nextflow per paral·lelitzar l'execució.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Executar el workflow amb múltiples mostres

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Conclusió

Sabeu com aprofitar el paradigma de flux de dades de Nextflow per paral·lelitzar el processament per mostra a través de múltiples mostres d'entrada.

### Què segueix?

A la [Part 3](03_multi_sample.md), afegireu l'agregació de múltiples mostres per combinar resultats de totes les mostres.
