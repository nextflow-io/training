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

- Haver completat el tutorial [Hello Nextflow](../../hello_nextflow/index.md) o un curs equivalent per a principiants.
- Estar còmodes utilitzant conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors, mòduls)

---

## 0. Primers passos

#### Obriu el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../../envsetup/index.md).

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

L'editor s'obre amb el directori del projecte en primer pla.

#### Reviseu els materials

Trobareu un directori `modules` amb definicions de processos, un directori `workflows` amb dos scripts de workflow escrits prèviament, i un fitxer `main.nf` que anireu actualitzant progressivament:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

El directori `modules/` conté les definicions individuals de processos, i el directori `workflows/` conté els dos scripts de workflow escrits prèviament amb els quals treballareu en aquesta missió secundària.

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

## 1. Afegiu el greeting workflow al pipeline

El greeting workflow valida noms i genera salutacions amb marca de temps.

### 1.1. Reviseu i executeu el greeting workflow

Obriu `workflows/greeting.nf` i examineu el codi:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Encadena processos: valida -> crea salutació -> afegeix marca de temps
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Aquest és un workflow complet i autònom amb la mateixa estructura que vau veure al tutorial 'Hello Nextflow'.
Té els noms d'entrada codificats, encadena tres processos i publica dues sortides.

Executeu-lo per verificar que tot funciona:

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

Per fer-lo composable amb altres workflows, cal canviar algunes coses.

### 1.2. Feu el workflow composable

Per fer un workflow composable, cal canviar quatre coses:
el workflow rep un nom, les entrades es mouen a un bloc `take:`, les sortides es mouen a un bloc `emit:`,
i els blocs autònoms `publish:`/`output {}` s'eliminen (pertanyen al entry workflow).

Anem veient aquests canvis un per un.

#### 1.2.1. Doneu nom al workflow

Doneu un nom al workflow perquè es pugui importar des d'un workflow pare.

=== "Després"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Abans"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Amb un nom, el workflow es pot importar en altres scripts.

#### 1.2.2. Declareu les entrades amb `take:`

Substituïu la declaració del canal codificada per un bloc `take:` que declari quines entrades espera el workflow.
El bloc `take:` va abans de `main:`, i s'elimina la línia `names_ch = channel.of(...)`.

=== "Després"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Canal d'entrada amb noms

        main:
        // Encadena processos: valida -> crea salutació -> afegeix marca de temps
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Abans"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Encadena processos: valida -> crea salutació -> afegeix marca de temps
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

El bloc `take:` declara el canal només pel nom — els detalls del que hi entra els definirà el workflow pare.

#### 1.2.3. Declareu les sortides amb `emit:`

Substituïu la secció `publish:` i elimineu el bloc `output {}`, reemplaçant-los per un bloc `emit:` que nomeni les sortides.

=== "Després"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Salutacions originals
        timestamped = timestamped_ch // Salutacions amb marca de temps
    }
    ```

=== "Abans"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

El bloc `emit:` exposa sortides amb nom a les quals els workflows pare poden accedir mitjançant `GREETING_WORKFLOW.out.greetings` i `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Verifiqueu el resultat i proveu-lo

Després dels tres canvis, el fitxer complet hauria de tenir aquest aspecte:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Canal d'entrada amb noms

    main:
    // Encadena processos: valida -> crea salutació -> afegeix marca de temps
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Salutacions originals
    timestamped = timestamped_ch // Salutacions amb marca de temps
}
```

Ara proveu d'executar-lo directament:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Això introdueix un concepte clau: el **entry workflow**.
Nextflow utilitza un bloc `workflow {}` sense nom com a punt d'entrada quan executeu un script directament.
`GREETING_WORKFLOW` té nom, de manera que Nextflow no sap com executar-lo per si sol.

Això és intencionat — els workflows composables estan dissenyats per ser cridats des d'un entry workflow, no per executar-se directament.
La solució és un entry workflow a `main.nf` que importi i cridi `GREETING_WORKFLOW`.

### 1.3. Actualitzeu i proveu el workflow principal

Ara actualitzem el workflow principal per cridar el greeting workflow.

#### 1.3.1. Incloeu el greeting workflow i crideu-lo

Afegiu la sentència `include`, actualitzeu el cos del workflow per cridar `GREETING_WORKFLOW` i substituïu el marcador de posició `channel.empty()` a `publish:`:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Executa el greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

El entry workflow es manté sense nom perquè Nextflow l'utilitzi com a punt d'entrada del pipeline.

#### 1.3.2. Actualitzeu el bloc de sortida

Afegiu una directiva `path` per dirigir les salutacions publicades a un subdirectori `greetings/`:

=== "Després"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Executeu el workflow

Executeu el workflow per comprovar que funciona:

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
    ```

??? abstract "Contingut del directori"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "Contingut del fitxer"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Els fitxers de salutació es publiquen a `results/greetings/`.
El workflow principal crida `GREETING_WORKFLOW` i connecta la seva sortida directament a la secció `publish:`.

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

## 2. Afegiu el transform workflow al pipeline

El transform workflow aplica transformacions de text a les salutacions amb marca de temps.

### 2.1. Reviseu i executeu el workflow

Obriu `workflows/transform.nf` i examineu el codi:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Aplica les transformacions en seqüència
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Aquest workflow autònom llegeix els fitxers de salutació amb marca de temps del directori `results/` produït per `greeting.nf`, els converteix a majúscules i després inverteix el text.

Executeu-lo per verificar que funciona amb els resultats del greeting de la secció 1.1:

```bash
nextflow run workflows/transform.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

Per fer-lo composable amb `GREETING_WORKFLOW`, cal aplicar els mateixos tres canvis de la secció 1.2.

### 2.2. Feu-lo composable

Apliqueu els mateixos tres canvis que a la secció 1.2: doneu nom al workflow, substituïu l'entrada codificada per `take:`, i substituïu `publish:`/`output {}` per `emit:`.

El fitxer acabat hauria de tenir aquest aspecte:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Canal d'entrada amb missatges

    main:
    // Aplica les transformacions en seqüència
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Salutacions en majúscules
    reversed = reversed_ch // Salutacions en majúscules invertides
}
```

El transform workflow ara és composable i està llest per ser importat al workflow principal.

### 2.3. Actualitzeu i proveu el workflow principal

Ara actualitzem el workflow principal per cridar el transform workflow.

#### 2.3.1. Incloeu el transform workflow i crideu-lo

Afegiu la sentència include, una crida a `TRANSFORM_WORKFLOW` encadenada a les salutacions amb marca de temps, i les dues noves entrades a `publish:`:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Executa el greeting workflow
        GREETING_WORKFLOW(names)

        // Executa el transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Executa el greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Això executarà el transform workflow sobre les salutacions amb marca de temps.

#### 2.3.2. Actualitzeu el bloc de sortida

Afegiu les entrades `upper` i `reversed` al bloc `output {}`, cadascuna amb una directiva `path` per al seu subdirectori:

=== "Després"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Això publicarà les sortides finals als directoris corresponents.

#### 2.3.3. Executeu el pipeline complet

Executeu el pipeline per comprovar que tot funciona:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Contingut del directori"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "Contingut del fitxer"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

El pipeline funciona de principi a fi: la salutació s'ha convertit a majúscules i s'ha invertit.

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

És important tenir en compte que, tot i que cridar workflows és una mica com cridar processos, en realitat no és el mateix. No podeu, per exemple, executar un workflow N vegades cridant-lo amb un canal de mida N — hauríeu de passar un canal de mida N al workflow i iterar internament.

Aplicar aquestes tècniques en el vostre propi treball us permetrà construir pipelines de Nextflow més sofisticats que puguin gestionar tasques complexes de processament de dades mantenint-se mantenibles i escalables.

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

2.  **Importacions de workflow:** Hem construït dos mòduls de workflow independents i els hem importat en un pipeline principal amb sentències `include`.

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

Torneu al [menú de missions secundàries](../index.md) o feu clic al botó a la part inferior dreta de la pàgina per continuar amb el tema següent de la llista.
