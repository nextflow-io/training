# Part 1: Adaptació a l'entorn de còmput

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A [Nextflow Run](../nextflow_run/index.md), heu configurat les entrades, els paràmetres i les sortides d'un pipeline.
Aquest curs cobreix l'altra meitat: adaptar l'execució d'un pipeline a qualsevol entorn de còmput on s'hagi d'executar, sense modificar el codi del workflow.

!!! example "Escenari"

    Heu desenvolupat i provat el vostre pipeline al vostre ordinador portàtil amb Docker.
    Ara cal transferir-lo: un col·laborador només té Conda configurat, i el clúster HPC de la vostra institució espera que les tasques passin pel seu propi planificador amb els seus propis límits de recursos.
    Res d'això hauria de requerir reescriure el pipeline.

El mateix codi del pipeline pot executar-se en tots aquests llocs, perquè res d'això està integrat al workflow.
L'empaquetament de programari, la plataforma d'execució i l'assignació de recursos es controlen mitjançant la configuració, en capes sobre el codi, i això és el que cobreix aquest curs: com adaptar el mateix pipeline a un nou entorn canviant la configuració, no el codi.

---

## 1. Seleccioneu una tecnologia d'empaquetament de programari

A [Nextflow Run](../nextflow_run/index.md), heu vist un perfil `conda` ja configurat a `nextflow.config` com a alternativa a Docker.
Aquí construireu vosaltres mateixos aquest mateix canvi, i veureu què cal fer perquè un procés sigui realment utilitzable amb Conda.

### 1.1. Desactiveu Docker i activeu Conda

Canvieu `docker.enabled` a `false` i afegiu una directiva que activi Conda.

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Això permet a Nextflow crear i utilitzar entorns Conda per a qualsevol procés que tingui un paquet Conda especificat.
El procés `cowpy` encara no en té cap, així que n'afegirem un, completament des de la configuració.

### 1.2. Afegiu un paquet Conda mitjançant la configuració

Una directiva `conda` es pot definir a la pròpia definició del procés, de la mateixa manera que `container` ja ho és a `modules/cowpy.nf`, però no és obligatori: `withName` us permet definir-la des de la configuració, limitada només al procés `cowpy`.

=== "Després"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Això no substitueix la directiva `container` que ja hi ha al codi del pipeline, sinó que afegeix una alternativa al seu costat, sense tocar gens aquest codi.

!!! tip "Consell"

    La cerca de [Seqera Containers](https://seqera.io/containers/) és una manera convenient de cercar l'URI del paquet Conda per a una eina determinada, fins i tot si no teniu previst construir-ne un contenidor.

### 1.3. Executeu el workflow per verificar que pot utilitzar Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/execution-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Això produeix la mateixa sortida que executar amb Docker, tot i que la mecànica és diferent entre bastidors: Nextflow recupera el paquet Conda i construeix un entorn a partir d'ell, en lloc de descarregar una imatge de contenidor.

!!! info "Info"

    Construir un nou entorn Conda pot trigar una mica més que descarregar un contenidor la primera vegada, però el paquet utilitzat aquí és petit, de manera que hauria de ser ràpid.

Ara torneu a Docker per a la resta d'aquest curs.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Combinar Docker i Conda"

    Com que aquests paràmetres estan delimitats per procés, podeu combinar-los: alguns processos utilitzen Docker, d'altres utilitzen Conda, depenent del que estigui disponible per a cada eina.
    Si tant una directiva `container` (al codi del pipeline) com una directiva `conda` (aquí, des de la configuració) estan definides per al mateix procés i tots dos sistemes d'empaquetament estan activats, Nextflow prioritza els contenidors.

### Conclusió

Sabeu com configurar quina tecnologia d'empaquetament de programari ha d'utilitzar un procés, i com canviar entre Docker i Conda.

### Què segueix?

Apreneu a canviar la plataforma d'execució que utilitza Nextflow per executar les vostres tasques.

---

## 2. Seleccioneu una plataforma d'execució

Tots els pipelines que heu executat fins ara han utilitzat l'executor local: cada tasca s'executa a la mateixa màquina que Nextflow.
Nextflow comprova les CPU i la memòria disponibles, i reté les tasques fins que hi hagi prou recursos lliures.

L'executor local és convenient, però no escala més enllà d'una sola màquina.
Nextflow admet [molts altres backends d'execució](https://nextflow.io/docs/latest/executor.html), incloent planificadors HPC (Slurm, LSF, SGE, PBS i d'altres) i plataformes al núvol (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i més).

### 2.1. Apunteu a un backend diferent

L'executor es defineix amb una directiva de procés anomenada `executor`.
Per defecte és `local`, de manera que el següent està implícit:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Per apuntar a un backend diferent, definiu la directiva a l'executor que voleu.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Advertència"

    L'entorn de formació no està connectat a un clúster HPC, de manera que això no és quelcom que pugueu executar aquí.

### 2.2. La sintaxi específica del backend s'abstreu

La majoria de plataformes HPC requereixen que les trameses de tasques especifiquin sol·licituds de recursos, com ara CPU, memòria i un nom de cua, utilitzant la seva pròpia sintaxi.
La mateixa sol·licitud de 8 CPU i 4 GB de RAM en una cua anomenada `my-science-work` té un aspecte completament diferent depenent del planificador.

??? abstract "Exemples"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Nextflow abstreu tot això: especifiqueu propietats estandarditzades com ara `cpus`, `memory` i `queue` una sola vegada (vegeu les [directives de procés](https://nextflow.io/docs/latest/reference/process.html#process-directives) per a la llista completa), i Nextflow les tradueix als scripts específics del backend corresponent en temps d'execució.

### 2.3. Veieu què executa realment Nextflow

Aquesta traducció no és només una comoditat del fitxer de configuració: està recolzada per quelcom concret que podeu inspeccionar ara mateix, fins i tot amb l'executor local.
A [Nextflow Run, secció 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory), heu mirat dins d'un directori de tasca sota `work/` i heu trobat `.command.sh`, la comanda exacta que Nextflow va executar.
Aquest mateix directori també conté un fitxer que encara no heu vist: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Sortida de la comanda (extracte)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` és l'script real que Nextflow lliura per a l'execució.
Embolcalla `.command.sh` amb tot el necessari per executar-lo realment: configuració de l'entorn, staging d'entrades/sortides i notificació del resultat a Nextflow.
Amb l'executor `local`, Nextflow simplement executa aquest script a la mateixa màquina.

Això és exactament el que canvia quan definiu un `executor` diferent.
Per a un planificador HPC com Slurm o PBS, Nextflow genera el mateix tipus d'script embolcall, afegeix la capçalera específica del planificador que heu vist a [2.2](#22-backend-specific-syntax-is-abstracted-away) (traduïda des dels vostres paràmetres `cpus`, `memory` i `queue`), i lliura el resultat a la comanda de tramesa pròpia d'aquell planificador, per exemple `sbatch` per a Slurm.
A partir d'aquí, Nextflow consulta l'estat de les tasques al planificador en lloc de supervisar directament un procés local.
Els backends de cloud batch funcionen una mica diferent, ja que es gestionen mitjançant crides a l'API en lloc d'una comanda de tramesa, però la mateixa idea subjacent s'aplica: el mateix script de tasca s'executa, només canvia com es llança i es fa el seguiment.

### Conclusió

Sabeu com canviar l'executor per apuntar a una infraestructura de còmput diferent, que Nextflow abstreu la sintaxi de tramesa específica del backend, i què passa realment entre bastidors quan una tasca s'executa en un backend diferent.

### Què segueix?

Continueu amb la [Part 2](./02_resources_and_retries.md), on aprendreu a perfilar i assignar recursos de còmput, i a gestionar els errors de les tasques amb reintents.

---

## Resum

En aquesta part heu après a:

- Canviar la tecnologia d'empaquetament de programari entre Docker i Conda
- Afegir una directiva `conda` a una definició de procés
- Canviar la plataforma d'execució amb la directiva `executor`
- Inspeccionar el que Nextflow genera i executa realment per a una tasca, i com canvia entre executors
