# Part 1: Executar un pipeline de demostració

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta primera part del curs Use nf-core, us mostrem com trobar un pipeline d'nf-core i provar-lo utilitzant el seu perfil de prova integrat.

Farem servir un pipeline anomenat nf-core/demo que és mantingut pel projecte nf-core com a part del seu inventari de pipelines per a demostracions i formació.

Assegureu-vos que el vostre directori de treball estigui configurat a `nfcore-use/` tal com s'indica a la pàgina de [Primers passos](./00_orientation.md).

---

## 1. Trobar i recuperar el pipeline nf-core/demo

Comencem localitzant el pipeline nf-core/demo al lloc web del projecte a [nf-co.re](https://nf-co.re), que centralitza tota la informació com ara: documentació general i articles d'ajuda, documentació per a cadascun dels pipelines, entrades de blog, anuncis d'esdeveniments i molt més.

### 1.1. Trobar el pipeline al lloc web

Al vostre navegador web, aneu a [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) i escriviu `demo` a la barra de cerca.

![resultats de la cerca](./img/search-results.png)

Feu clic al nom del pipeline, `demo`, per accedir a la pàgina de documentació del pipeline.

Cada pipeline publicat té una pàgina dedicada que inclou les seccions de documentació següents:

- **Introduction:** Una introducció i visió general del pipeline
- **Usage:** Descripcions de com executar el pipeline
- **Parameters:** Paràmetres del pipeline agrupats amb descripcions
- **Output:** Descripcions i exemples dels fitxers de sortida esperats
- **Results:** Fitxers de sortida d'exemple generats a partir del conjunt de dades de prova complet
- **Releases & Statistics:** Historial de versions del pipeline i estadístiques

Sempre que estigueu considerant adoptar un nou pipeline, hauríeu de llegir la documentació del pipeline amb atenció primer per entendre què fa i com s'ha de configurar abans d'intentar executar-lo.

Feu-hi un cop d'ull ara i vegeu si podeu esbrinar:

- Quines eines executarà el pipeline (Consulteu la pestanya: `Introduction`)
- Quines entrades i paràmetres accepta o requereix el pipeline (Consulteu la pestanya: `Parameters`)
- Quines són les sortides produïdes pel pipeline (Consulteu la pestanya: `Output`)

#### 1.1.1. Visió general del pipeline

La pestanya `Introduction` proporciona una visió general del pipeline, incloent-hi una representació visual (anomenada mapa de metro) i una llista d'eines que s'executen com a part del pipeline.

![mapa de metro del pipeline](./img/nf-core-demo-subway-cropped.png)

1. Control de qualitat de lectures ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Retallada d'adaptadors i qualitat ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Presentació del control de qualitat per a lectures en brut ([MULTIQC](http://multiqc.info/))
4. Generació d'un missatge de text divertit d'una vaca ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Exemple de línia de comandes

La documentació també proporciona un fitxer d'entrada d'exemple (comentat més endavant) i una línia de comandes d'exemple.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Notareu que la comanda d'exemple NO especifica un fitxer de workflow, només la referència al repositori del pipeline, `nf-core/demo`.

Quan s'invoca d'aquesta manera, Nextflow assumirà que el codi està organitzat d'una determinada manera.
Recuperem el codi per poder examinar aquesta estructura.

### 1.2. Recuperar el codi del pipeline

Un cop hem determinat que el pipeline sembla adequat per als nostres propòsits, provem-lo.
Afortunadament, Nextflow facilita la recuperació de pipelines de repositoris correctament formatats sense haver de descarregar res manualment.

#### 1.2.1. Utilitzar `nextflow pull`

Tornem al terminal i executem el següent:

```bash
nextflow pull nf-core/demo
```

??? success "Sortida de la comanda"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow fa un `pull` del codi del pipeline, és a dir, descarrega el repositori complet al vostre disc local.

Per ser clars, podeu fer això amb qualsevol pipeline de Nextflow que estigui configurat adequadament a GitHub, no només amb els pipelines d'nf-core.
No obstant això, nf-core és la col·lecció de codi obert de pipelines de Nextflow més gran.

#### 1.2.2. Utilitzar `nextflow list`

Podeu demanar a Nextflow que us doni una llista dels pipelines que heu recuperat d'aquesta manera:

```bash
nextflow list
```

??? success "Sortida de la comanda"

    ```console
    nf-core/demo
    ```

Podeu provar de descarregar uns quants pipelines més per veure com apareixen llistats quan en teniu més d'un.

#### 1.2.3. Trobar on s'ha descarregat el pipeline

Notareu que els fitxers no es troben al vostre directori de treball actual.
Per defecte, Nextflow desa els pipelines descarregats a `$NXF_HOME/assets`.

Per trobar on viu un pipeline específic, pregunteu-ho directament a Nextflow:

```bash
nextflow info nf-core/demo
```

??? success "Sortida de la comanda"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Info"

    El camí complet pot diferir al vostre sistema si no esteu utilitzant el nostre entorn de formació.

Nextflow manté el codi font descarregat intencionadament "fora del camí" seguint el principi que aquests pipelines s'han d'utilitzar més com a biblioteques que com a codi amb el qual interactuaríeu directament.

Internament, Nextflow emmagatzema cada pipeline descarregat com un repositori git a `$NXF_HOME/assets/.repos/`, i extreu el codi per a cada revisió en un subdirectori `clones/<commit>/`.
Com que `.repos` és un directori ocult, un simple `tree -L 2 $NXF_HOME/assets/` semblarà buit.

#### 1.2.4. Crear un enllaç simbòlic per accedir fàcilment al codi font

No examinarem el codi en detall, però fem-hi un cop d'ull ràpid per tenir una idea de com és l'organització general.

Per facilitar la navegació pel codi font del pipeline, creeu un enllaç simbòlic que apunti a la còpia extreta del pipeline:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Això crea una drecera perquè pugueu explorar el codi amb `tree -L 2 pipelines/nf-core/demo` o obrir fitxers directament.

#### 1.2.5. Visió general de l'organització del codi

Podeu utilitzar `tree` o l'explorador de fitxers per trobar i obrir el directori `nf-core/demo`.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Contingut del directori"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows

    7 directories, 12 files
    ```

Com podeu veure, hi ha moltes coses, la majoria de les quals no cal que us preocupeu.

Breument, observem que al nivell superior podeu trobar un fitxer README amb informació resumida, així com fitxers accessoris que resumeixen informació del projecte com ara llicències, directrius de contribució, citació i codi de conducta.
La documentació detallada del pipeline es troba al directori `docs`.
Tot aquest contingut s'utilitza per generar les pàgines web del lloc web d'nf-core de manera programàtica, de manera que sempre estan actualitzades amb el codi.

Per a la resta, podem distingir tres grups funcionals de fitxers de codi:

1. Components del codi del pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuració del pipeline
3. Paràmetres / entrades i validació del pipeline

No repassarem els components del codi del pipeline en aquesta part del curs, però sí que tractarem elements de configuració i validació que probablement us seran rellevants com a usuaris finals de pipelines d'nf-core.

!!! tip "Consell"

    També podeu navegar pel codi font de qualsevol pipeline d'nf-core a GitHub, per exemple [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Tots els pipelines d'nf-core segueixen el mateix disseny de directoris, de manera que un cop coneixeu l'estructura, podeu trobar fitxers de configuració, mòduls i workflows per a qualsevol pipeline de la mateixa manera.

Ara, a executar el pipeline!

### Conclusio

Ara sabeu com trobar un pipeline al lloc web d'nf-core i recuperar una còpia local del codi font.

### Què segueix?

Apreneu com provar un pipeline d'nf-core amb el mínim esforç.

---

## 2. Provar el pipeline amb el seu perfil de prova

Convenientment, tots els pipelines d'nf-core inclouen un perfil de prova.
Es tracta d'un conjunt mínim de paràmetres de configuració perquè el pipeline s'executi utilitzant un petit conjunt de dades de prova allotjat al repositori [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
És una manera excel·lent de provar ràpidament un pipeline a petita escala.

!!! tip "Consell"

    El sistema de perfils de configuració de Nextflow us permet canviar fàcilment entre diferents motors de contenidors o entorns d'execució.
    Per a més detalls, consulteu [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar el perfil de prova

És una bona pràctica comprovar què especifica el perfil de prova d'un pipeline abans d'executar-lo.
El perfil `test` per a `nf-core/demo` es troba al fitxer de configuració `conf/test.config`.
Podeu trobar-lo localment dins del codi font del pipeline que `nextflow pull` va descarregar, a través de l'enllaç simbòlic `pipelines` creat a la secció 1.2.4:

```bash
code pipelines/nf-core/demo/conf/test.config
```

Aquí teniu el contingut d'aquest fitxer:

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Dades d'entrada
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Notareu de seguida que el bloc de comentaris a la part superior inclou un exemple d'ús que mostra com executar el pipeline amb aquest perfil de prova.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Les úniques coses que hem de proporcionar són les que es mostren entre claudàtors angulars a la comanda d'exemple: `<docker/singularity>` i `<OUTDIR>`.

Com a recordatori, `<docker/singularity>` fa referència a l'elecció del sistema de contenidors. Tots els pipelines d'nf-core estan dissenyats per ser utilitzables amb contenidors (Docker, Singularity, etc.) per garantir la reproductibilitat i eliminar problemes d'instal·lació de programari.
Per tant, haurem d'especificar si volem utilitzar Docker o Singularity per provar el pipeline.

La part `--outdir <OUTDIR>` fa referència al directori on Nextflow escriurà les sortides del pipeline.
Hem de proporcionar-li un nom, que podem inventar-nos.
Si no existeix ja, Nextflow el crearà per nosaltres en temps d'execució.

Continuant amb la secció posterior al bloc de comentaris, el perfil de prova ens mostra què s'ha preconfigurar per a les proves: el més destacat és que el paràmetre `input` ja està configurat per apuntar a un conjunt de dades de prova, de manera que no cal que proporcionem les nostres pròpies dades.
Si seguiu l'enllaç a l'entrada preconfigurada, veureu que és un fitxer CSV que conté identificadors de mostres i camins de fitxers per a diverses mostres experimentals.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Això s'anomena full de mostres (samplesheet), i és la forma d'entrada més habitual als pipelines d'nf-core.
No us preocupeu si no esteu familiaritzats amb els formats i tipus de dades, no és important per al que segueix.

Ara tenim tot el que necessitem per provar el pipeline.

### 2.2. Executar el pipeline

Tal com s'ha indicat anteriorment, podem utilitzar la comanda de prova d'exemple gairebé tal com és; només cal especificar quin empaquetament de programari volem utilitzar i quin nom donar al directori de sortida.
Aquí farem servir Docker com a sistema de contenidors i `demo-results`, respectivament.

Amb això, podem executar la comanda de prova:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Si la vostra sortida coincideix amb aquesta, felicitats! Acabeu d'executar el vostre primer pipeline d'nf-core.

Notareu que hi ha molta més sortida de consola que quan executeu un pipeline bàsic de Nextflow.
Hi ha una capçalera que inclou un resum de la versió del pipeline, les entrades i sortides, i alguns elements de configuració.

!!! info "Info"

    La vostra sortida mostrarà marques de temps, noms d'execució i camins de fitxers diferents, però l'estructura general i l'execució dels processos haurien de ser similars.

Fixeu-vos en la línia prop de la part superior de la sortida:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Això us indica quina revisió del pipeline s'ha utilitzat.
Com que no hem especificat cap versió, Nextflow ha utilitzat l'últim commit a `master`.
Per a execucions reproduïbles, hauríeu de fixar una versió específica utilitzant l'indicador `-r`:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

Això garanteix que s'utilitzi el mateix codi del pipeline cada vegada, independentment de nous commits o versions.
Per a aquesta formació ometem `-r` per simplicitat, però en producció sempre hauríeu d'especificar-lo.

Continuant amb la sortida d'execució, fem un cop d'ull a les línies que ens indiquen quins processos s'han executat:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Això ens indica que s'han executat quatre processos, corresponents a les quatre eines que es mostren a la pàgina de documentació del pipeline al lloc web d'nf-core: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` i `COWPY`.

Els noms complets dels processos tal com es mostren aquí, com ara `NFCORE_DEMO:DEMO:MULTIQC`, són més llargs del que potser heu vist al material introductori de Hello Nextflow.
Aquests inclouen els noms dels seus workflows pare i reflecteixen la modularitat del codi del pipeline.
Si voleu aprendre a desenvolupar pipelines d'estil nf-core vosaltres mateixos, consulteu el curs [Build with nf-core](../nfcore_build/index.md).

### 2.3. Examinar les sortides del pipeline

Finalment, fem un cop d'ull al directori `demo-results` produït pel pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contingut del directori"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Pot semblar molt.
Per obtenir més informació sobre les sortides del pipeline `nf-core/demo`, consulteu la seva [pàgina de documentació](https://nf-co.re/demo/1.2.0/docs/output/).

En aquesta etapa, el que és important observar és que els resultats estan organitzats per mòdul, i a més hi ha un directori anomenat `pipeline_info` que conté diversos informes amb marca de temps sobre l'execució del pipeline.

Per exemple, el fitxer `execution_timeline_*` us mostra quins processos s'han executat, en quin ordre i quant de temps han trigat:

![informe de la línia de temps d'execució](./img/execution_timeline.png)

!!! info "Info"

    Aquí les tasques no s'han executat en paral·lel perquè estem executant en una màquina minimalista a Github Codespaces.
    Per veure-les executar-se en paral·lel, proveu d'augmentar l'assignació de CPU del vostre codespace i els límits de recursos a la configuració de prova.

Aquests informes es generen automàticament per a tots els pipelines d'nf-core.

### Conclusio

Sabeu com executar un pipeline d'nf-core utilitzant el seu perfil de prova integrat i on trobar les seves sortides.

### Què segueix?

Aneu a la [Part 2](./02_configure_execution.md), on aprendreu a configurar l'execució del pipeline.

---

## Resum

En aquesta part heu après a:

- Trobar i recuperar un pipeline d'nf-core i examinar la seva estructura de codi
- Executar un pipeline utilitzant el seu perfil de prova integrat
