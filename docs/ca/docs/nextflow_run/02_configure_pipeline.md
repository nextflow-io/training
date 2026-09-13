# Part 2: Configurar el pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la [Part 1](./01_run_nextflow.md), heu executat un pipeline complet de múltiples passos que processa múltiples entrades en paral·lel utilitzant contenidors.
Ara veurem com configurar el comportament del pipeline mitjançant `nextflow.config`: primer examinant el fitxer de configuració que ja us hem proporcionat, després explorant un parell d'altres maneres de subministrar configuració, i finalment controlant com i on es publiquen les sortides.

---

## 1. Examinar el fitxer de configuració principal

Nextflow recull automàticament `nextflow.config` del directori de treball i aplica la seva configuració a cada execució.

Us proporcionem un fitxer de configuració que cobreix quatre àrees: empaquetament de programari, configuració de processos, paràmetres del pipeline i perfils d'execució.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Empaquetament de programari
     */
    docker.enabled = true

    /*
     * Configuració de processos
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Paràmetres del pipeline
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Perfils
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Repassem cadascun d'ells i, a continuació, posem els perfils en pràctica executant el pipeline amb un d'ells.

!!! note "Nota"

    Aquesta configuració cobreix l'execució local en una sola màquina.
    Nextflow també admet planificadors HPC (SLURM, PBS, LSF) i executors al núvol (AWS Batch, Google Cloud Batch, Azure Batch), tots configurats mitjançant el mateix mecanisme `nextflow.config`.
    Consulteu la [Part 1: Adapt to your compute environment](../execution_config/01_packaging_and_execution.md) del curs [Execution Config](../execution_config/index.md) per a una guia completa d'aquestes opcions.

### 1.1. Empaquetament de programari

L'empaquetament de programari és la manera com Nextflow subministra les eines reals que necessiten els vostres processos, ja sigui una imatge de contenidor, un entorn Conda o qualsevol altra cosa.

```groovy title="nextflow.config" linenums="1"
/*
 * Empaquetament de programari
 */
docker.enabled = true
```

Aquesta línia habilita Docker per a cada procés.
Qualsevol procés que declari una directiva `container` s'executa dins de la imatge especificada.

### 1.2. Configuració de processos

Recordeu que un procés és un pas únic del vostre pipeline, com ara `sayHello` o `cowpy`.
Nextflow us permet configurar diverses coses sobre com s'executa cadascun: quanta CPU i memòria obté, quin contenidor o entorn Conda utilitza, i més.

```groovy title="nextflow.config" linenums="6"
/*
 * Configuració de processos
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Això limita cada procés a una sola CPU i 1 GB de memòria.

Nextflow també us permet establir valors diferents per a processos individuals amb nom o grups de processos; aprendreu com fer-ho a la [Part 2: Manage compute resources and failures](../execution_config/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) del curs [Execution Config](../execution_config/index.md).

### 1.3. Paràmetres del pipeline

Els paràmetres són les entrades de línia de comandes del pipeline, els mateixos indicadors `--input`, `--batch` i `--character` que ja heu estat establint directament a la línia de comandes.
Establir valors per defecte aquí significa que no els heu d'escriure cada vegada, tot i que, com veureu més endavant en aquesta part, hi ha un parell d'altres maneres de subministrar-los.

```groovy title="nextflow.config" linenums="14"
/*
 * Paràmetres del pipeline
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Aquests valors per defecte s'activen sempre que no es proporciona un paràmetre a la línia de comandes, de manera que executar `nextflow run main.nf` sense cap indicador continua funcionant.

### 1.4. Perfils

Els perfils us permeten agrupar un conjunt de configuracions sota un sol nom, de manera que podeu canviar entre configuracions completes amb un sol indicador en lloc de canviar els valors manualment cada vegada.

```groovy title="nextflow.config" linenums="23"
/*
 * Perfils
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

El perfil `test` sobreescriu tres paràmetres per executar el pipeline amb un conjunt d'entrades petit i ben definit; cada pipeline nf-core n'inclou un per a la validació ràpida, i és una convenció que val la pena seguir en els vostres propis pipelines.

El perfil `conda` canvia l'empaquetament de programari de Docker a Conda.

Activeu un perfil passant `-profile <nom>` a la línia de comandes.

Posem el perfil `test` en pràctica.

```bash
nextflow run main.nf -profile test
```

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

El pipeline s'executa amb `batch = 'test'` i `character = 'tux'`.
Comproveu `results/test/`: el nom del lot ara forma part del camí del directori, i l'art ASCII mostra el pingüí tux en lloc d'un gall dindi.

!!! note "Nota"

    Podeu activar diversos perfils alhora i utilitzar `nextflow config -profile <nom>,<nom>` per veure el resultat completament resolt abans d'executar res.
    La combinació de perfils i com Nextflow resol els conflictes entre ells es tracta en profunditat a la [Part 3: Use profiles to switch configurations](../execution_config/03_profiles.md) del curs [Execution Config](../execution_config/index.md).

### Conclusió

Sabeu per a què serveixen els elements més comuns d'un fitxer `nextflow.config` i com activar un perfil.

### Què segueix?

Apreneu un parell d'altres maneres de subministrar valors de configuració sense modificar el fitxer `nextflow.config` principal, útils per configurar execucions individuals i per compartir un conjunt exacte de configuracions amb una altra persona.

---

## 2. Proporcionar configuració mitjançant fitxers suplementaris

Establir valors per defecte a `nextflow.config` funciona bé per a valors que rarament canvien.
Nextflow també us ofereix dos mecanismes més específics: un fitxer de configuració específic per a una execució per adaptar-la a un entorn particular, i un fitxer de paràmetres per compartir un conjunt exacte de valors d'entrada amb un col·laborador.

### 2.1. Utilitzar un fitxer de configuració específic per a una execució

Suposem que esteu movent el pipeline a una màquina que no té Docker i voleu donar a cada procés més recursos per treballar.
Creeu un nou fitxer de configuració amb només les sobreescriptures que necessiteu:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Passeu-lo juntament amb el vostre pipeline principal amb `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow fusiona `custom.config` per sobre del `nextflow.config` propi del pipeline, de manera que cada procés ara obté 2 CPUs i 2 GB de memòria en lloc dels valors per defecte, i s'executa amb Conda en lloc de Docker.
`cowpy` és l'únic procés amb un paquet Conda declarat juntament amb el seu contenidor, de manera que és el que veureu que Nextflow construeix realment un entorn per a ell:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Un fitxer petit que només sobreescriu l'assignació de recursos i l'empaquetament, sense tocar els paràmetres del pipeline, és exactament el patró que els pipelines nf-core esperen de les configuracions institucionals.
Consulteu el repositori [nf-core/configs](https://github.com/nf-core/configs) per a exemples del món real.

Això us ofereix una manera desechable d'adaptar un pipeline a un nou entorn sense tocar la vostra configuració habitual.

### 2.2. Utilitzar un fitxer de paràmetres

Suposem que necessiteu compartir un conjunt exacte de paràmetres d'execució amb un col·laborador, o registrar-los per a una publicació.

Nextflow us permet subministrar [fitxers de paràmetres](https://nextflow.io/docs/latest/config.html#parameter-file) en format YAML o JSON, que són una manera més senzilla de distribuir un conjunt exacte i reproduïble de valors.

Un fitxer de paràmetres anomenat `test-params.yaml` ja es troba al vostre directori de treball:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

La sintaxi utilitza dos punts (`:`) en lloc dels signes d'igual (`=`) que s'utilitzen a `nextflow.config`, ja que aquest fitxer és YAML simple en lloc de Groovy.

!!! info "Info"

    També es proporciona una versió JSON, `test-params.json`. Podeu provar-la pel vostre compte; la sintaxi per passar-la és idèntica.

Passeu el fitxer amb `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Contingut del fitxer"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Un fitxer de paràmetres és especialment valuós quan un pipeline té més d'uns quants paràmetres: us permet subministrar-los tots alhora, sense una línia de comandes extensa ni cap canvi a l'script del workflow, i és fàcil de distribuir juntament amb els vostres resultats.

### Conclusió

Coneixeu dues maneres més de subministrar configuració: un fitxer de configuració específic per a una execució per adaptar-la a un nou entorn, i un fitxer de paràmetres per compartir valors d'entrada exactes i reproduïbles.

### Què segueix?

Apreneu a controlar com i on es publiquen les sortides del vostre pipeline.

---

## 3. Gestionar les sortides del pipeline

L'autor d'un pipeline decideix com s'organitzen les sortides en el codi, però no cal que toqueu aquest codi per controlar on acaben o com hi arriben.
Nextflow us ofereix maneres de fer-ho a nivell de configuració: establir un directori de sortida base i triar si els fitxers es copien o s'enllacen simbòlicament.

### 3.1. Personalitzar el directori de sortida

Per defecte, Nextflow publica les sortides sota `results/`.
Apunteu-lo a un altre lloc amb `-output-dir` (o la seva forma abreujada, `-o`):

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Contingut del directori"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Les sortides ara es troben sota `outputs/batch/` en lloc del valor per defecte integrat `results/batch/`.
El codi propi del pipeline continua decidint l'estructura dins d'aquest directori base, com ara els subdirectoris `batch/` i `intermediates/`; `-output-dir` només controla on comença aquesta estructura.

`-output-dir` és realment una drecera de línia de comandes per a l'opció de configuració `outputDir`, de manera que pot anar a qualsevol lloc on pugui anar la configuració: directament a `nextflow.config`, dins d'un perfil, o en un fitxer de superposició `-c` com el que heu utilitzat abans en aquesta part.
Per exemple, aquest fragment mostra la mateixa configuració col·locada directament a `nextflow.config` en lloc de passar-la a la línia de comandes:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Consulteu [Configuration file](https://nextflow.io/docs/latest/config.html) a la referència de Nextflow per a la llista completa de llocs on pot viure una opció de configuració com aquesta.

### 3.2. Triar com es publiquen les sortides

Per defecte, Nextflow publica les sortides com a enllaços simbòlics que apunten a les ubicacions de les sortides sota `work/`, no còpies reals:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Els autors del pipeline poden establir el 'mode de publicació' a `'copy'` o `'move'` per a cada procés individual en el codi del workflow.
Normalment ho fan per a les sortides finals del pipeline, mentre deixen el comportament per defecte `'symlink'` establert per als fitxers intermedis que es poden eliminar un cop s'ha executat el pipeline complet.

Això evita duplicar dades al disc, però significa que no podeu eliminar els directoris de tasques sota `work/` sense trencar l'enllaç, perdent la capacitat d'utilitzar `-resume`.
Si voleu que tots els fitxers de sortida es copiïn correctament, establiu [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) a `'copy'` a la configuració del vostre pipeline. (A diferència de `-output-dir`, no hi ha cap indicador de línia de comandes per a això; és només de configuració.)

Proveu d'establir-ho a `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Abans"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

A continuació, executeu el pipeline canviant el nom del lot per poder veure la diferència a les sortides:

```bash
nextflow run main.nf --batch withmode
```

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Examineu un dels fitxers de sortida com abans:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Ara és un fitxer real i independent que continuarà disponible fins i tot si `work/` es neteja.

!!! warning "Advertència"

    La configuració `workflow.output.mode` només omple un valor per defecte per a les sortides que encara no tenen un mode establert al codi del pipeline.
    No pot sobreescriure un mode que l'autor ha codificat de manera fixa, independentment del que establiu.

### Conclusió

Sabeu com personalitzar el directori de sortida base i triar entre sortides copiades i enllaçades simbòlicament, tot sense tocar el codi del pipeline.

### Què segueix?

Continueu a la [Part 3](./03_manage_executions.md), on aprendreu a inspeccionar l'historial d'execucions passades, generar informes d'execució i netejar directoris de treball antics.

---

## Resum

En aquesta part heu après a:

- Configurar el comportament del pipeline mitjançant `nextflow.config` i perfils
- Subministrar configuració mitjançant un fitxer de configuració específic per a una execució o un fitxer de paràmetres
- Personalitzar el directori de sortida i triar entre sortides copiades i enllaçades simbòlicament
