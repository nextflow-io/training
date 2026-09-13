# Part 3: Utilitzar perfils per canviar configuracions

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


A través de la [Part 1](./01_packaging_and_execution.md) i la [Part 2](./02_resources_and_retries.md), heu acumulat diverses opcions de configuració: empaquetament de programari, plataforma d'execució i assignació de recursos.
A la pràctica, sovint voldreu canviar entre conjunts sencers d'aquestes opcions depenent d'on esteu executant, per exemple un ordinador portàtil per al desenvolupament i un clúster HPC per a producció.

Nextflow us permet configurar qualsevol nombre de [perfils](https://nextflow.io/docs/latest/config.html#profiles) que descriuen configuracions diferents, i seleccionar-ne un (o diversos) en temps d'execució amb un únic indicador.

Ja n'heu utilitzat un: el perfil `test` de [Nextflow Run](../nextflow_run/index.md) sobreescriu els paràmetres d'entrada amb un conjunt petit i ben definit.
Ara creareu els vostres propis perfils d'infraestructura i els combinareu amb aquest.

---

## 1. Crear perfils per a entorns diferents

### 1.1. Configurar els perfils

Afegiu dos perfils a `nextflow.config`: un per executar en un ordinador portàtil normal amb Docker, i un per a un clúster HPC universitari amb un planificador Slurm i Conda.

=== "Després"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
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
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="35"
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

El perfil `univ_hpc` també estableix límits de recursos, ja que això és habitualment necessari en infraestructures HPC compartides.

### 1.2. Executar el workflow amb un perfil

Seleccioneu un perfil en temps d'execució amb `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Advertència"

    El perfil `univ_hpc` no s'executarà en l'entorn de formació, ja que no hi ha cap planificador Slurm disponible.

Si trobeu altres configuracions que sempre van juntes, afegiu-les al perfil corresponent.
També podeu crear perfils addicionals per agrupar qualsevol altra combinació que necessiteu.

### 1.3. Executar amb múltiples perfils

Els perfils no són mútuament excloents.
Podeu activar-ne diversos alhora amb `-profile <perfil1>,<perfil2>`.
Combineu `my_laptop` amb el perfil `test` que ja coneixeu de Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Els noms dels fitxers individuals recullen correctament `batch = 'test'` del perfil `test` (`COLLECTED-test-output.txt`, i així successivament).

Si combineu perfils que estableixen la mateixa opció, Nextflow resol el conflicte utilitzant el valor que llegeix en darrer lloc, és a dir, el que apareix més tard al fitxer.
Si les configuracions en conflicte provenen de fonts de configuració completament diferents, s'aplica l'[ordre de precedència](https://www.nextflow.io/docs/latest/config.html) estàndard.

### Conclusió

Sabeu com definir perfils que agrupen configuració específica d'infraestructura, seleccionar-ne un en temps d'execució amb `-profile`, combinar múltiples perfils en una sola execució, i com Nextflow resol els conflictes quan més d'un perfil estableix la mateixa opció.

### Què segueix?

Apreneu a inspeccionar la configuració completament resolta abans d'executar res.

---

## 2. Inspeccionar la configuració resolta

Ja heu utilitzat `nextflow config -profile test` a [Nextflow Run](../nextflow_run/02_configure_pipeline.md) per comprovar a què es resol un únic perfil.
Aquesta comanda esdevé especialment útil quan combineu múltiples perfils: com acabeu de veure, quan dos perfils estableixen la mateixa opció, pot ser complicat esbrinar manualment quin valor guanya realment.
La comanda `nextflow config` ho resol tot per vosaltres, sense executar el pipeline.

### 2.1. Resoldre la configuració per defecte

```bash
nextflow config
```

??? success "Sortida de la comanda"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

Això és exactament el que s'aplicaria si executéssiu el pipeline sense cap indicador addicional.

### 2.2. Resoldre la configuració amb perfils activats

Afegiu els mateixos perfils que utilitzaríeu per a una execució real.

```bash
nextflow config -profile my_laptop,test
```

??? success "Sortida de la comanda"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

Comparar els dos resultats confirma què ha canviat: `params.batch`, `params.character` i `process.executor` reflecteixen tots els perfils `my_laptop,test`.
Això és especialment valuós per a pipelines amb moltes capes de configuració, on esbrinar manualment les configuracions resoltes seria tediós i propens a errors.

### Conclusió

Sabeu com utilitzar `nextflow config` per inspeccionar la configuració completament resolta per a qualsevol combinació de perfils, abans d'executar res.

### Què segueix?

Heu cobert els aspectes essencials de la configuració de pipelines Nextflow.
Consulteu el [Resum del curs](next_steps.md) per saber on anar a partir d'aquí.

---

## Resum

En aquesta part heu après a:

- Definir perfils que agrupen configuració específica d'infraestructura
- Combinar múltiples perfils en una sola execució, i entendre com es resolen els conflictes entre ells
- Utilitzar `nextflow config` per inspeccionar la configuració completament resolta
