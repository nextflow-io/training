# Workflows de Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Quan esteu desenvolupant un pipeline, sovint us trobeu creant seqüències similars de processos per a diferents tipus de dades o passos d'anàlisi. Podríeu acabar copiant i enganxant aquestes seqüències de processos, cosa que porta a codi duplicat difícil de mantenir; o podríeu crear un workflow massiu que sigui difícil d'entendre i modificar.

Una de les característiques més potents de Nextflow és la seva capacitat per compondre pipelines complexos a partir de mòduls de workflow més petits i reutilitzables. Aquest enfocament modular fa que els pipelines siguin més fàcils de desenvolupar, provar i mantenir.

### Objectius d'aprenentatge

En aquesta missió secundària, explorarem com desenvolupar mòduls de workflow que es puguin provar i utilitzar per separat, compondre aquests mòduls en un pipeline més gran i gestionar el flux de dades entre mòduls.

Al final d'aquesta missió secundària, sereu capaços de:

- Descompondre pipelines complexos en unitats lògiques i reutilitzables
- Provar cada mòdul de workflow de manera independent
- Combinar workflows per crear nous pipelines
- Compartir mòduls de workflow comuns entre diferents pipelines
- Fer el vostre codi més mantenible i fàcil d'entendre

Aquestes habilitats us ajudaran a construir pipelines complexos mantenint una estructura de codi neta i mantenible.

### Prerequisits

Abans d'abordar aquesta missió secundària hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Estar còmodes utilitzant conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors, mòduls)

---

## 0. Primers passos

#### Obriu el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Movem-nos al directori on es troben els fitxers d'aquest tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Podeu configurar VSCode perquè es centri en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un directori `modules` que conté diverses definicions de processos que es basen en el que vau aprendre a 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### Reviseu l'assignació

El vostre repte és assemblar aquests mòduls en dos workflows separats que després composarem en un workflow principal:

- Un `GREETING_WORKFLOW` que valida noms, crea salutacions i afegeix marques de temps
- Un `TRANSFORM_WORKFLOW` que converteix text a majúscules i l'inverteix

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Llista de verificació de preparació

Creieu que esteu preparats per submergir-vos?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu codespace està en funcionament
- [ ] He configurat el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu a punt per començar.

---

## 1. Creeu el Greeting Workflow

Comencem creant un workflow que validi noms i generi salutacions amb marca de temps.

### 1.1. Creeu l'estructura del workflow

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Afegiu el codi del primer (sub)workflow

Afegiu aquest codi a `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Encadena processos: valida -> crea salutació -> afegeix marca de temps
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Aquest és un workflow complet, amb una estructura similar a les que vau veure al tutorial 'Hello Nextflow', que podem provar de manera independent. Provem-ho ara:

```bash
nextflow run workflows/greeting.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Això funciona com s'esperava, però per fer-lo composable hi ha algunes coses que hem de canviar.

### 1.3. Feu el workflow composable

Els workflows composables tenen algunes diferències respecte als que vau veure al tutorial 'Hello Nextflow':

- El bloc del workflow ha de tenir un nom
- Les entrades es declaren utilitzant la paraula clau `take:`
- El contingut del workflow es col·loca dins del bloc `main:`
- Les sortides es declaren utilitzant la paraula clau `emit:`

Actualitzem el greeting workflow perquè coincideixi amb aquesta estructura. Canvieu el codi pel següent:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Canal d'entrada amb noms

    main:
        // Encadena processos: valida -> crea salutació -> afegeix marca de temps
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Salutacions originals
        timestamped = timestamped_ch  // Salutacions amb marca de temps
}
```

Podeu veure que el workflow ara té un nom i té un bloc `take:` i `emit:`, i aquestes són les connexions que utilitzarem per compondre un workflow de nivell superior.
El contingut del workflow també es col·loca dins del bloc `main:`. Noteu també que hem eliminat la declaració del canal d'entrada `names_ch`, ja que ara es passa com a argument al workflow.

Provem el workflow de nou per veure si funciona com s'esperava:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Això us informa d'un altre concepte nou, un 'entry workflow'. L'entry workflow és el workflow que s'invoca quan executeu un script de Nextflow. Per defecte, Nextflow utilitzarà un workflow sense nom com a entry workflow, quan estigui present, i això és el que heu estat fent fins ara, amb blocs de workflow que comencen així:

```groovy title="hello.nf" linenums="1"
workflow {
```

Però el nostre greeting workflow no té un workflow sense nom, sinó que tenim un workflow amb nom:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Per això Nextflow va llançar un error i no va fer el que volíem.

No hem afegit la sintaxi `take:`/`emit:` per poder cridar el workflow directament - ho hem fet per poder compondre'l amb altres workflows. La solució és crear un script principal amb un entry workflow sense nom que importi i cridi el nostre workflow amb nom.

### 1.4. Creeu i proveu el workflow principal

Ara crearem un workflow principal que importi i utilitzi el workflow `greeting`.

Creeu `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Noteu que l'entrada del workflow en aquest fitxer no té nom, i és perquè l'utilitzarem com a entry workflow.

Executeu-lo i observeu la sortida:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Funciona! Hem embolcallat el greeting workflow amb nom en un workflow principal amb un bloc `workflow` d'entrada sense nom. El workflow principal utilitza el workflow `GREETING_WORKFLOW` gairebé (no exactament) com un procés, i passa el canal `names` com a argument.

### Conclusió

En aquesta secció, heu après diversos conceptes importants:

- **Workflows amb nom**: Crear un workflow amb nom (`GREETING_WORKFLOW`) que es pot importar i reutilitzar
- **Interfícies de workflow**: Definir entrades clares amb `take:` i sortides amb `emit:` per crear un workflow composable
- **Punts d'entrada**: Entendre que Nextflow necessita un entry workflow sense nom per executar un script
- **Composició de workflows**: Importar i utilitzar un workflow amb nom dins d'un altre workflow
- **Espais de noms de workflow**: Accedir a les sortides del workflow utilitzant l'espai de noms `.out` (`GREETING_WORKFLOW.out.greetings`)

Ara teniu un greeting workflow funcional que:

- Rep un canal de noms com a entrada
- Valida cada nom
- Crea una salutació per a cada nom vàlid
- Afegeix marques de temps a les salutacions
- Exposa tant les salutacions originals com les salutacions amb marca de temps com a sortides

Aquest enfocament modular us permet provar el greeting workflow de manera independent o utilitzar-lo com a component en pipelines més grans.

---

## 2. Afegiu el Transform Workflow

Ara creem un workflow que apliqui transformacions de text a les salutacions.

### 2.1. Creeu el fitxer del workflow

```bash
touch workflows/transform.nf
```

### 2.2. Afegiu el codi del workflow

Afegiu aquest codi a `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Canal d'entrada amb missatges

    main:
        // Aplica les transformacions en seqüència
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Salutacions en majúscules
        reversed = reversed_ch  // Salutacions en majúscules invertides
}
```

No repetirem l'explicació de la sintaxi composable aquí, però noteu que el workflow amb nom es torna a declarar amb un bloc `take:` i `emit:`, i el contingut del workflow es col·loca dins del bloc `main:`.

### 2.3. Actualitzeu el workflow principal

Actualitzeu `main.nf` per utilitzar ambdós workflows:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Executa el greeting workflow
    GREETING_WORKFLOW(names)

    // Executa el transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Mostra els resultats
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Executeu el pipeline complet:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Si observeu un d'aquests fitxers invertits, veureu que és la versió en majúscules de la salutació invertida:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Conclusió

Ara hauríeu de tenir un pipeline complet que:

- Processa noms a través del greeting workflow
- Alimenta les salutacions amb marca de temps al transform workflow
- Produeix versions en majúscules i invertides de les salutacions

---

## Resum

En aquesta missió secundària, hem explorat el potent concepte de composició de workflows a Nextflow, que ens permet construir pipelines complexos a partir de components més petits i reutilitzables.

Aquest enfocament modular ofereix diversos avantatges respecte als pipelines monolítics:

- Cada workflow es pot desenvolupar, provar i depurar de manera independent
- Els workflows es poden reutilitzar en diferents pipelines
- L'estructura general del pipeline es torna més llegible i mantenible
- Els canvis en un workflow no afecten necessàriament els altres si les interfícies es mantenen consistents
- Els punts d'entrada es poden configurar per executar diferents parts del vostre pipeline segons sigui necessari

_És important tenir en compte, però, que tot i que cridar workflows és una mica com cridar processos, en realitat no és el mateix. No podeu, per exemple, executar un workflow N vegades cridant-lo amb un canal de mida N - hauríeu de passar un canal de mida N al workflow i iterar internament._

Aplicar aquestes tècniques en el vostre propi treball us permetrà construir pipelines de Nextflow més sofisticats que puguin gestionar tasques complexes de bioinformàtica mantenint-se mantenibles i escalables.

### Patrons clau

1.  **Estructura del workflow**: Hem definit entrades i sortides clares per a cada workflow utilitzant la sintaxi `take:` i `emit:`, creant interfícies ben definides entre components, i hem embolcallat la lògica del workflow dins del bloc `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Els canals d'entrada es declaren aquí
            input_ch

        main:
            // La lògica del workflow va aquí
            // Aquí és on es criden els processos i es manipulen els canals
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Els canals de sortida es declaren aquí
            output_ch = result_ch
    }
    ```

2.  **Importacions de workflow:** Hem construït dos mòduls de workflow independents i els hem importat en un pipeline principal amb sentències include.

    - Incloure un sol workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Incloure múltiples workflows

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Incloure amb àlies per evitar conflictes de noms

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Punts d'entrada**: Nextflow requereix un entry workflow sense nom per saber on iniciar l'execució. Aquest entry workflow crida els vostres workflows amb nom.

    - Workflow sense nom (punt d'entrada)

    ```groovy
    workflow {
        // Aquest és el punt d'entrada quan s'executa l'script
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow amb nom (cridat des de l'entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Ha de ser cridat des de l'entry workflow
    }
    ```

4.  **Gestió del flux de dades:** Hem après com accedir a les sortides del workflow utilitzant la notació d'espai de noms (`WORKFLOW_NAME.out.channel_name`) i passar-les a altres workflows.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Recursos addicionals

- [Documentació de Workflow de Nextflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Referència d'operadors de canal](https://www.nextflow.io/docs/latest/operator.html)
- [Documentació d'estratègia d'errors](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Què segueix?

Torneu al [menú de missions secundàries](../) o feu clic al botó a la part inferior dreta de la pàgina per continuar amb el tema següent de la llista.
