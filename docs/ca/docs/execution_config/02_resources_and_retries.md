# Part 2: Gestionar recursos de còmput i fallades

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la [Part 1](./01_packaging_and_execution.md), heu adaptat on i com s'executen les tasques d'un pipeline.
Aquí adaptareu quant còmput rep cada tasca i què passa quan una tasca falla malgrat la vostra millor estimació d'assignació.

---

## 1. Controlar les assignacions de recursos de còmput

Per defecte, Nextflow assigna una sola CPU a cada procés mitjançant la directiva `cpus`, i no imposa un límit de memòria tret que en definiu un:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

Ja sabeu de [Nextflow Run](../nextflow_run/index.md) que la configuració d'aquest pipeline defineix `memory` a 1 GB per a tots els processos.
Però com sabeu quins valors heu d'utilitzar realment per als vostres propis pipelines?

### 1.1. Generar un informe d'utilització de recursos

Ja heu generat un informe d'execució amb `-with-report` a [Nextflow Run](../nextflow_run/02_configure_pipeline.md).
Aquest mateix informe és com descobriu quanta CPU i memòria necessiten realment els vostres processos: executeu el workflow amb algunes assignacions per defecte, registreu l'ús real i ajusteu a partir d'aquí.

```bash
nextflow run main.nf -with-report report-config-1.html
```

L'informe és un fitxer HTML que podeu obrir en un navegador.
Desglossa el temps d'execució i la utilització de recursos per procés, incloent-hi quin percentatge dels recursos assignats s'ha utilitzat realment.
Aquí teniu el que mostra per a `cowpy` amb els valors per defecte actuals (1 CPU, 1 GB de memòria):

| Mètrica              | Valor  |
| -------------------- | ------ |
| Ús de CPU            | 116%   |
| Memòria màxima usada | 6.4 MB |
| Memòria assignada    | 1 GB   |

`cowpy` utilitza molt menys de l'1% de la seva assignació d'1 GB; el `%cpu` per sobre del 100% simplement significa que utilitza breument més d'una CPU de processament dins del contenidor, en ràfegues curtes.

Consulteu [Reports](https://nextflow.io/docs/latest/reports.html) per a la llista completa de funcionalitats disponibles.

### 1.2. Definir assignacions de recursos per a un procés específic

L'informe anterior mostra que `cowpy` es troba còmodament dins de la seva assignació actual, però suposem que voleu donar-li més marge igualment, per exemple perquè espereu entrades més grans en producció.
Podeu sobreescriure els valors per defecte per a un sol procés amb `withName`.

=== "Després"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Amb això al seu lloc, cada procés sol·licita 1 GB de memòria i una sola CPU, excepte `cowpy`, que sol·licita 2 GB i 2 CPUs (a més de la configuració de `conda` de la [Part 1](./01_packaging_and_execution.md)).

!!! info "Info"

    Si la vostra màquina té poques CPUs i n'assigneu un nombre elevat per procés, les crides a les tasques poden posar-se en cua les unes darrere les altres, ja que Nextflow no sol·licitarà més CPUs de les disponibles.

Executeu-lo de nou amb un nom de fitxer d'informe diferent, per poder comparar abans i després.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

Comparant els dos informes per a `cowpy`:

| Mètrica              | Abans (1 CPU, 1 GB) | Després (2 CPUs, 2 GB) |
| -------------------- | ------------------- | ---------------------- |
| Memòria màxima usada | 6.4 MB              | 6.4 MB                 |
| Ús de CPU            | 116%                | 118%                   |

Doblar l'assignació no va canviar gens l'ús real, la qual cosa us indica que l'original d'1 GB / 1 CPU ja era generós per a aquesta càrrega de treball de prova.
En un pipeline real que processa dades no trivials, esperaríeu que els números diferissin de manera significativa entre processos, que és exactament per això que feu el perfil abans de decidir què assignar, en lloc d'endevinar.

### 1.3. Afegir límits de recursos

Depenent de la vostra infraestructura de còmput, pot haver-hi restriccions estrictes sobre el que podeu sol·licitar, per exemple un límit a tot el clúster.
La directiva `resourceLimits` us permet definir aquests límits:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow tradueix aquests límits al que espera l'executor de destinació.
Si un procés sol·licita més del límit, la sol·licitud es limita en lloc de rebutjar-se.

!!! warning "Advertència"

    Això no és quelcom que pugueu executar a l'entorn de formació, ja que requereix infraestructura HPC per tenir efecte.

??? info "Configuracions de referència institucionals"

    El projecte nf-core manté una [col·lecció de fitxers de configuració](https://nf-co.re/configs/) compartits per institucions de tot el món, que cobreixen una àmplia gamma d'executors HPC i al núvol.
    Són un bon punt de partida tant si la vostra institució hi és com si no.

### Conclusió

Sabeu com generar un informe de perfil per avaluar la utilització de recursos, sobreescriure les assignacions de recursos per a un procés específic i limitar les assignacions amb `resourceLimits`.

### Què segueix?

Apreneu a fer que un pipeline es recuperi automàticament quan una tasca falla, tant si la vostra estimació d'assignació de recursos era correcta com si no.

---

## 2. Gestionar les fallades de tasques amb reintents

El perfil us indica el que necessita un procés la majoria de les vegades, però les càrregues de treball reals varien: una assignació còmoda per a la majoria d'entrades pot ser massa justa per a una d'inusualment gran, i les estimacions simplement poden ser incorrectes.
En lloc de deixar que una sola tasca fallida ensorri tota l'execució, Nextflow pot reintentar una tasca fallida automàticament, opcionalment donant-li més recursos en cada intent.

### 2.1. Reintentar una tasca fallida automàticament

Per veure-ho en acció, definiu deliberadament l'assignació de memòria de `cowpy` per sota del que realment necessita: recordeu de la [1.1](#11-generate-a-resource-utilization-report) que arriba a un màxim d'uns 6.4 MB, de manera que 6 MB hauria de ser just per sota del necessari.

=== "Després"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy` indica a Nextflow què ha de fer quan una tasca falla: `'retry'` resubmet la tasca en lloc d'aturar tot el pipeline.
`maxRetries` limita quants intents addicionals té abans que Nextflow es rendeixi.

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda (abreujada)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/execution-config/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

El codi de sortida 137 és el senyal estàndard per a una terminació per falta de memòria: el contenidor no tenia prou memòria per executar `cowpy` en absolut.
Nextflow va reintentar la tasca dues vegades, tres intents en total, coincidint amb `maxRetries = 2`.
Com que l'assignació de memòria no va canviar entre intents, cada intent va topar amb el mateix obstacle; un cop esgotats els reintents, Nextflow informa de la fallada completament i atura el pipeline, sortint amb un estat no zero.

Reintentar per si sol no soluciona res si la causa subjacent no canvia entre intents.

### 2.2. Augmentar els recursos en cada reintent

Dins d'una directiva de procés, `task.attempt` conté el número de l'intent actual, començant per 1.
Podeu utilitzar-lo en una closure per escalar una assignació de recursos cap amunt amb cada reintent.

=== "Després"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

Executeu el workflow de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda (abreujada)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

El primer intent falla igualment amb 6 MB, però el reintent s'executa amb 12 MB (`6.MB * 2`) i té èxit, i el pipeline es completa amb totes les sortides publicades.

!!! warning "Advertència"

    La sortida de la consola encara inclou una línia `NOTE:` que informa del primer intent fallat, tot i que el pipeline en conjunt ha tingut èxit: Nextflow registra cada reintent individualment, però una fallada reintentada no afecta el resultat global.
    Comproveu el resum `Outputs:`, o l'estat de sortida de la comanda, per confirmar si l'execució ha tingut èxit realment.

Consulteu [Dynamic computing resources](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) a la documentació de Nextflow per a patrons de reintent més avançats, incloent-hi l'escalat basat en quin error específic s'ha produït.

### Conclusió

Sabeu com fer que un pipeline reintenti automàticament les tasques fallides i com escalar les assignacions de recursos amb cada reintent utilitzant `task.attempt`.

### Què segueix?

Continueu amb la [Part 3](./03_profiles.md), on aprendreu a agrupar configuració com aquesta en perfils commutables.

---

## Resum

En aquesta part heu après a:

- Generar un informe de perfil de recursos i definir assignacions de recursos per procés
- Limitar les sol·licituds de recursos amb `resourceLimits`
- Reintentar automàticament una tasca fallida amb `errorStrategy` i `maxRetries`
- Escalar una assignació de recursos cap amunt amb cada reintent utilitzant `task.attempt`
