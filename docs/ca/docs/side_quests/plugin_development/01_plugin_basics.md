# Part 1: Conceptes bàsics dels plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta secció, aprendreu com els plugins amplien Nextflow i, a continuació, provareu tres plugins diferents per veure'ls en acció.

---

## 1. Com funcionen els plugins

Els plugins amplien Nextflow mitjançant diversos tipus d'extensió:

| Tipus d'extensió    | Què fa                                                  | Exemple                      |
| ------------------- | ------------------------------------------------------- | ---------------------------- |
| Funcions            | Afegeix funcions personalitzades invocables des de workflows | `samplesheetToList()`        |
| Monitors de workflow | Responen a esdeveniments com la finalització de tasques | Registre personalitzat, alertes de Slack |
| Executors           | Afegeix backends d'execució de tasques                  | AWS Batch, Kubernetes        |
| Sistemes de fitxers | Afegeix backends d'emmagatzematge                       | S3, Azure Blob               |

Les funcions i els monitors de workflow (anomenats "trace observers" a l'API de Nextflow) són els tipus més habituals per als autors de plugins.
Els executors i els sistemes de fitxers els creen habitualment els proveïdors de plataformes.

Els exercicis següents us mostren plugins de funcions i un plugin d'observador, perquè pugueu veure ambdós tipus en acció.

---

## 2. Ús de plugins de funcions

Els plugins de funcions afegeixen funcions invocables que importeu als vostres workflows.
En provareu dos: nf-hello (un exemple senzill) i nf-schema (un plugin real d'ús generalitzat).
Tots dos exercicis modifiquen el mateix pipeline `hello.nf`, de manera que podeu veure com els plugins milloren un workflow existent.

### 2.1. nf-hello: substituir codi escrit manualment

El plugin [nf-hello](https://github.com/nextflow-io/nf-hello) proporciona una funció `randomString` que genera cadenes de text aleatòries.
El pipeline ja defineix la seva pròpia versió en línia d'aquesta funció, que substituireu per la del plugin.

#### 2.1.1. Veure el punt de partida

Examineu el pipeline:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Genera una cadena alfanumèrica aleatòria
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

El pipeline defineix la seva pròpia funció `randomString` en línia i l'utilitza per afegir un identificador aleatori a cada salutació.

Executeu-lo:

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

L'ordre de la vostra sortida i les cadenes aleatòries seran diferents, i si torneu a executar l'script obtindreu un conjunt diferent de salutacions aleatòries.

#### 2.1.2. Configurar el plugin

Substituïu la funció en línia per una del plugin. Afegiu això al vostre `nextflow.config`:

```groovy title="nextflow.config"
// Configuració per als exercicis de desenvolupament de plugins
plugins {
    id 'nf-hello@0.5.0'
}
```

Els plugins es declaren a `nextflow.config` mitjançant el bloc `plugins {}`.
Nextflow els descarrega automàticament des del [Nextflow Plugin Registry](https://registry.nextflow.io/), un repositori central de plugins de la comunitat i oficials.

#### 2.1.3. Usar la funció del plugin

Substituïu la funció `randomString` en línia per la versió del plugin:

=== "Després"

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

=== "Abans"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Genera una cadena alfanumèrica aleatòria
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

La instrucció `include` importa `randomString` d'una biblioteca provada, testada i mantinguda per un grup més ampli de col·laboradors que poden detectar i corregir errors.
En lloc que cada pipeline mantingui la seva pròpia còpia de la funció, tots els pipelines que utilitzen el plugin obtenen la mateixa implementació verificada.
Això redueix el codi duplicat i la càrrega de manteniment que comporta.
La sintaxi `#!groovy include { function } from 'plugin/plugin-id'` és el mateix `include` que s'utilitza per als mòduls de Nextflow, amb el prefix `plugin/`.
Podeu consultar el [codi font de `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) al repositori nf-hello a GitHub.

#### 2.1.4. Executar-lo

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

(Les vostres cadenes aleatòries seran diferents.)

La sortida encara té sufixos aleatoris, però ara `randomString` prové del plugin nf-hello en lloc del codi en línia.
Els missatges "Pipeline is starting!" i "Pipeline complete!" són nous.
Provenen del component observador del plugin, que explorareu a la Part 5.

Nextflow descarrega els plugins automàticament la primera vegada que s'utilitzen, de manera que qualsevol pipeline que declari `nf-hello@0.5.0` obté exactament la mateixa funció `randomString` testada sense haver de copiar codi entre projectes.

Ara ja heu vist els tres passos per usar un plugin de funcions: declarar-lo a `nextflow.config`, importar la funció amb `include` i cridar-la al vostre workflow.
El proper exercici aplica aquests mateixos passos a un plugin del món real.

### 2.2. nf-schema: anàlisi de CSV amb validació

El plugin [nf-schema](https://github.com/nextflow-io/nf-schema) és un dels plugins de Nextflow més utilitzats.
Proporciona `samplesheetToList`, una funció que analitza fitxers CSV/TSV mitjançant un esquema JSON que defineix les columnes i els tipus esperats.

El pipeline llegeix actualment `greetings.csv` amb `splitCsv` i un `map` manual, però nf-schema pot substituir-ho per una anàlisi validada i basada en esquemes.
Ja hi ha un fitxer d'esquema JSON (`greetings_schema.json`) al directori de l'exercici.

??? info "Què és un esquema?"

    Un esquema és una descripció formal de com han de ser les dades vàlides.
    Defineix coses com quines columnes s'esperen, quin tipus ha de tenir cada valor (cadena de text, número, etc.) i quins camps són obligatoris.

    Penseu-hi com un contracte: si les dades d'entrada no coincideixen amb l'esquema, l'eina pot detectar el problema aviat en lloc de deixar que causi errors confusos més endavant al pipeline.

#### 2.2.1. Examinar l'esquema

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

L'esquema defineix dues columnes (`greeting` i `language`) i marca `greeting` com a obligatòria.
Si algú passa un CSV sense la columna `greeting`, nf-schema detecta l'error abans que s'executi el pipeline.

#### 2.2.2. Afegir nf-schema a la configuració

Actualitzeu `nextflow.config` per incloure tots dos plugins:

=== "Després"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. Actualitzar hello.nf per usar samplesheetToList

Substituïu l'entrada `splitCsv` per `samplesheetToList`:

=== "Després"

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

=== "Abans"

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

El codi d'anàlisi personalitzat amb `splitCsv` i `map` es substitueix per `samplesheetToList`, una funció provada i testada que també valida el samplesheet contra l'esquema abans que s'executi el pipeline.
Això redueix la càrrega de manteniment de la lògica d'anàlisi escrita manualment i millora l'experiència dels usuaris del pipeline, que reben missatges d'error clars quan la seva entrada no coincideix amb el format esperat.
Cada fila es converteix en una llista de valors en ordre de columna, de manera que `row[0]` és la salutació i `row[1]` és l'idioma.

#### 2.2.4. Executar-lo

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

(Les vostres cadenes aleatòries seran diferents.)

La sortida és la mateixa, però ara l'esquema valida l'estructura del CSV abans que s'executi el pipeline.
En pipelines reals amb fulls de mostres complexos i moltes columnes, aquest tipus de validació evita errors que `splitCsv` + `map` manual passaria per alt.

#### 2.2.5. Veure la validació en acció

Per veure què detecta la validació d'esquemes, proveu d'introduir errors a `greetings.csv`.

Canvieu el nom de la columna obligatòria `greeting` per `message`:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Executeu el pipeline:

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

El pipeline es nega a executar-se perquè l'esquema requereix una columna `greeting` i no en troba cap.

Ara restaureu la columna obligatòria però canvieu el nom de la columna opcional `language` per `lang`:

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

Aquesta vegada el pipeline s'executa, però mostra un avís:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Les columnes obligatòries causen errors greus; les columnes opcionals causen avisos.
Aquest és el tipus de retroalimentació primerenca que estalvia temps de depuració en pipelines reals amb desenes de columnes.

#### 2.2.6. Configurar el comportament de la validació

L'avís sobre `lang` és útil, però podeu controlar la seva gravetat mitjançant la configuració.
Els plugins poden incloure els seus propis àmbits de configuració que controlen el seu comportament.
El plugin nf-schema inclou l'àmbit de configuració `validation`; modificant els paràmetres aquí podeu canviar com es comporta nf-schema.

Afegiu un bloc `validation` a `nextflow.config` perquè les capçaleres no reconegudes causin un error en lloc d'un avís:

=== "Després"

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

=== "Abans"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Torneu a executar el pipeline amb la mateixa columna `lang` encara present:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

El pipeline ara falla en lloc d'avisar.
El codi del pipeline no ha canviat; només ho ha fet la configuració.

Restaureu `greetings.csv` al seu estat original i elimineu el bloc `validation` abans de continuar:

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

Tant nf-hello com nf-schema són plugins de funcions: proporcionen funcions que importeu amb `include` i crideu al codi del vostre workflow.
El proper exercici mostra un tipus diferent de plugin que funciona sense cap instrucció `include`.

---

## 3. Ús d'un plugin observador: nf-co2footprint

No tots els plugins proporcionen funcions per importar.
El plugin [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) utilitza un **trace observer** per monitorar l'ús de recursos del vostre pipeline i estimar la seva petjada de carboni.
No cal canviar cap codi del pipeline; simplement afegiu-lo a la configuració.

### 3.1. Afegir nf-co2footprint a la configuració

Actualitzeu `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Executar el pipeline

```bash
nextflow run hello.nf
```

El plugin produeix diversos missatges INFO i WARN durant l'execució.
Aquests són normals per a un exemple petit que s'executa en una màquina local:

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

Els avisos sobre la zona, l'executor, el model de CPU i la memòria apareixen perquè el plugin no pot detectar tots els detalls del maquinari d'un entorn de formació local.
En un entorn de producció (per exemple, un clúster HPC o al núvol), aquests valors estarien disponibles i les estimacions serien més precises.

Al final, busqueu una línia com:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Els vostres números seran diferents.)

### 3.3. Veure l'informe

El plugin genera fitxers de sortida al vostre directori de treball:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Examineu el resum:

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

(Els vostres números seran diferents.)

La primera secció mostra les xifres brutes d'energia i emissions.
La secció "Which equals" posa aquests números en perspectiva convertint-los en equivalents familiars.
El resum també inclou una secció que llista les opcions de configuració del plugin i una citació a l'article de recerca [Green Algorithms](https://doi.org/10.1002/advs.202100707) en el qual es basa el mètode de càlcul.

### 3.4. Configurar el plugin

L'avís "Target zone null" de la secció 3.2 va aparèixer perquè el plugin no tenia cap ubicació configurada.
El plugin nf-co2footprint defineix un àmbit de configuració `co2footprint` on podeu establir la vostra ubicació geogràfica.

Afegiu un bloc `co2footprint` a `nextflow.config`:

=== "Després"

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

=== "Abans"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Consell"

    Useu el codi del vostre propi país si ho preferiu (per exemple, `'US'`, `'DE'`, `'FR'`).

Executeu el pipeline:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

L'avís de zona ha desaparegut.
El plugin ara utilitza la intensitat de carboni específica de GB (163.92 gCO₂eq/kWh) en lloc del valor global de reserva (480.0 gCO₂eq/kWh).

!!! note "Nota"

    També podeu veure un missatge `WARN: Unrecognized config option 'co2footprint.location'`.
    Això és cosmètic i es pot ignorar sense problemes; el plugin llegeix el valor correctament igualment.

A la Part 6, creareu un àmbit de configuració per al vostre propi plugin.

Aquest plugin funciona completament mitjançant el mecanisme d'observador, connectant-se als esdeveniments del cicle de vida del workflow per recollir mètriques de recursos i generar el seu informe quan el pipeline finalitza.

Ara ja heu provat plugins de funcions (importats amb `include`) i un plugin observador (activat només mitjançant la configuració).
Aquests són els dos tipus d'extensió més habituals, però tal com mostra la taula de la secció 1, els plugins també poden afegir executors i sistemes de fitxers.

---

## 4. Descobrir plugins

El [Nextflow Plugin Registry](https://registry.nextflow.io/) és el centre principal per trobar plugins disponibles.

![La pàgina del plugin nf-hello a registry.nextflow.io](img/plugin-registry-nf-hello.png)

Cada pàgina de plugin mostra la seva descripció, les versions disponibles, les instruccions d'instal·lació i enllaços a la documentació.

---

## 5. Preparar-se per al desenvolupament de plugins

Les seccions següents (Parts 2-6) utilitzen un fitxer de pipeline separat, `greet.nf`, que depèn de nf-schema però no de nf-hello ni de nf-co2footprint.

Actualitzeu `nextflow.config` per conservar només nf-schema:

```groovy title="nextflow.config"
// Configuració per als exercicis de desenvolupament de plugins
plugins {
    id 'nf-schema@2.6.1'
}
```

Elimineu els fitxers de sortida de co2footprint:

```bash
rm -f co2footprint_*
```

El fitxer `hello.nf` conserva el vostre treball de la Part 1 com a referència; a partir d'ara, treballareu amb `greet.nf`.

---

## Conclusió

Heu utilitzat tres plugins diferents:

- **nf-hello**: Un plugin de funcions que proporciona `randomString`, importat amb `include`
- **nf-schema**: Un plugin de funcions que proporciona `samplesheetToList` per a l'anàlisi de CSV amb validació d'esquemes
- **nf-co2footprint**: Un plugin observador que monitoritza l'ús de recursos automàticament, sense necessitat de cap `include`

Patrons clau:

- Els plugins es declaren a `nextflow.config` amb `#!groovy plugins { id 'plugin-name@version' }`
- Els plugins de funcions requereixen `#!groovy include { function } from 'plugin/plugin-id'`
- Els plugins observadors funcionen automàticament un cop declarats a la configuració
- Els plugins poden definir àmbits de configuració (per exemple, `#!groovy validation {}`, `#!groovy co2footprint {}`) per personalitzar el comportament
- El [Nextflow Plugin Registry](https://registry.nextflow.io/) llista els plugins disponibles

---

## Què segueix?

Les seccions següents us mostren com construir el vostre propi plugin.
Si no esteu interessats en el desenvolupament de plugins, podeu aturar-vos aquí o saltar directament al [Resum](summary.md).

[Continua a la Part 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
