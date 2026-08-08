# Part 1: Executar un pipeline de demostració

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta primera part del curs de formació Build with nf-core, us mostrem com trobar i provar un pipeline nf-core, configurar i personalitzar la seva execució per a les vostres necessitats, i entendre com la validació d'entrada protegeix contra errors comuns.

Utilitzarem un pipeline anomenat nf-core/demo que és mantingut pel projecte nf-core com a part del seu inventari de pipelines per a finalitats de demostració i formació.

Assegureu-vos que el vostre directori de treball està configurat a `hello-nf-core/` tal com s'indica a la pàgina [Primers passos](./00_orientation.md).

---

## 1. Trobar i recuperar el pipeline nf-core/demo

Comencem localitzant el pipeline nf-core/demo al lloc web del projecte a [nf-co.re](https://nf-co.re), que centralitza tota la informació com ara: documentació general i articles d'ajuda, documentació per a cadascun dels pipelines, entrades de blog, anuncis d'esdeveniments, etc.

### 1.1. Trobar el pipeline al lloc web

Al vostre navegador web, aneu a [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) i escriviu `demo` a la barra de cerca.

![resultats de cerca](./img/search-results.png)

Feu clic al nom del pipeline, `demo`, per accedir a la pàgina de documentació del pipeline.

Cada pipeline publicat té una pàgina dedicada que inclou les següents seccions de documentació:

- **Introduction:** Una introducció i visió general del pipeline
- **Usage:** Descripcions de com executar el pipeline
- **Parameters:** Paràmetres del pipeline agrupats amb descripcions
- **Output:** Descripcions i exemples dels fitxers de sortida esperats
- **Results:** Exemples de fitxers de sortida generats a partir del conjunt de dades de prova complet
- **Releases & Statistics:** Historial de versions del pipeline i estadístiques

Sempre que estigueu considerant adoptar un nou pipeline, hauríeu de llegir la documentació del pipeline amb atenció primer per entendre què fa i com s'ha de configurar abans d'intentar executar-lo.

Doneu-hi una ullada ara i veieu si podeu esbrinar:

- Quines eines executarà el pipeline (Consulteu la pestanya: `Introduction`)
- Quines entrades i paràmetres accepta o requereix el pipeline (Consulteu la pestanya: `Parameters`)
- Quines són les sortides produïdes pel pipeline (Consulteu la pestanya: `Output`)

#### 1.1.1. Visió general del pipeline

La pestanya `Introduction` proporciona una visió general del pipeline, incloent una representació visual (anomenada mapa de metro) i una llista d'eines que s'executen com a part del pipeline.

![mapa de metro del pipeline](./img/nf-core-demo-subway-cropped.png)

1. Read QC ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Present QC for raw reads ([MULTIQC](http://multiqc.info/))
4. Generate a lighthearted text message from a cow ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Exemple de línia de comandes

La documentació també proporciona un exemple de fitxer d'entrada (que es discuteix més endavant) i un exemple de línia de comandes.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Notareu que l'exemple de comanda NO especifica un fitxer de workflow, només la referència al repositori del pipeline, `nf-core/demo`.

Quan s'invoca d'aquesta manera, Nextflow assumirà que el codi està organitzat d'una manera determinada.
Recuperem el codi perquè puguem examinar aquesta estructura.

### 1.2. Recuperar el codi del pipeline

Un cop hem determinat que el pipeline sembla adequat per als nostres propòsits, provem-lo.
Afortunadament, Nextflow facilita la recuperació de pipelines des de repositoris correctament formatats sense haver de descarregar res manualment.

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

Nextflow fa un `pull` del codi del pipeline, és a dir, descarrega el repositori complet a la vostra unitat local.

Per ser clars, podeu fer això amb qualsevol pipeline Nextflow que estigui configurat adequadament a GitHub, no només pipelines nf-core.
No obstant això, nf-core és la col·lecció de codi obert més gran de pipelines Nextflow.

#### 1.2.2. Utilitzar `nextflow list`

Podeu fer que Nextflow us doni una llista dels pipelines que heu recuperat d'aquesta manera:

```bash
nextflow list
```

??? success "Sortida de la comanda"

    ```console
    nf-core/demo
    ```

Podeu provar de recuperar alguns altres pipelines per veure com apareixen llistats quan en teniu més d'un.

#### 1.2.3. Trobar on s'ha descarregat el pipeline

Notareu que els fitxers no estan al vostre directori de treball actual.
Per defecte, Nextflow desa els pipelines recuperats a `$NXF_HOME/assets`.

Per saber on es troba un pipeline específic, pregunteu-ho directament a Nextflow:

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
      1.0.0 [t]
      1.0.1 [t]
      1.0.2 [t]
      1.1.0 [t]
    > 1.2.0 [t]
    ```

!!! info "Info"

    El camí complet pot diferir al vostre sistema si no esteu utilitzant el nostre entorn de formació.

Nextflow manté el codi font descarregat intencionadament 'fora del camí' amb el principi que aquests pipelines s'haurien d'utilitzar més com a biblioteques que com a codi amb el qual interactuaríeu directament.

Internament, Nextflow emmagatzema cada pipeline recuperat com un repositori git a `$NXF_HOME/assets/.repos/`, i extreu el codi de cada revisió a un subdirectori `clones/<commit>/`.
Com que `.repos` és un directori ocult, un simple `tree -L 2 $NXF_HOME/assets/` semblarà buit.

#### 1.2.4. Crear un enllaç simbòlic per accedir fàcilment al codi font

No examinarem el codi en detall, però fem-hi una ullada ràpida per tenir una idea de com és l'organització general.

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

Com podeu veure, hi ha molt en marxa allà, però la majoria no us hauria de preocupar.

Breument, observem que al nivell superior podeu trobar un fitxer README amb informació resumida, així com fitxers accessoris que resumeixen informació del projecte com ara llicència, directrius de contribució, citació i codi de conducta.
La documentació detallada del pipeline es troba al directori `docs`.
Tot aquest contingut s'utilitza per generar les pàgines web al lloc web nf-core de manera programàtica, de manera que sempre estan actualitzades amb el codi.

Per a la resta, podem distingir tres grups funcionals de fitxers de codi:

1. Components del codi del pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuració del pipeline
3. Paràmetres del pipeline / entrades i validació

No repassarem els components del codi del pipeline en aquesta part del curs, però sí que tractarem elements de configuració i validació que probablement us seran rellevants com a usuaris finals de pipelines nf-core.

!!! tip "Consell"

    També podeu navegar pel codi font de qualsevol pipeline nf-core a GitHub, per exemple [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Tots els pipelines nf-core segueixen el mateix disseny de directoris, de manera que un cop coneixeu l'estructura, podeu trobar fitxers de configuració, mòduls i workflows per a qualsevol pipeline de la mateixa manera.

Però de moment, anem a executar el pipeline!

### Conclusió

Ara sabeu com trobar un pipeline a través del lloc web nf-core i recuperar una còpia local del codi font.

### Què segueix?

Apreneu com provar un pipeline nf-core amb un esforç mínim.

---

## 2. Provar el pipeline amb el seu perfil de prova

Convenientment, cada pipeline nf-core ve amb un perfil de prova.
Aquest és un conjunt mínim de paràmetres de configuració perquè el pipeline s'executi utilitzant un petit conjunt de dades de prova allotjat al repositori [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
És una manera excel·lent de provar ràpidament un pipeline a petita escala.

!!! tip "Consell"

    El sistema de perfils de configuració de Nextflow us permet canviar fàcilment entre diferents motors de contenidors o entorns d'execució.
    Per a més detalls, consulteu [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar el perfil de prova

És una bona pràctica comprovar què especifica el perfil de prova d'un pipeline abans d'executar-lo.
El perfil `test` per a `nf-core/demo` es troba al fitxer de configuració `conf/test.config`.
Podeu trobar-lo localment dins del codi font del pipeline que `nextflow pull` ha descarregat, mitjançant l'enllaç simbòlic `pipelines` creat a la secció 1.2.4:

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

Notareu immediatament que el bloc de comentaris a la part superior inclou un exemple d'ús que mostra com executar el pipeline amb aquest perfil de prova.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Les úniques coses que hem de proporcionar són el que es mostra entre claudàtors a l'exemple de comanda: `<docker/singularity>` i `<OUTDIR>`.

Com a recordatori, `<docker/singularity>` es refereix a l'elecció del sistema de contenidors. Tots els pipelines nf-core estan dissenyats per ser utilitzables amb contenidors (Docker, Singularity, etc.) per garantir la reproducibilitat i eliminar problemes d'instal·lació de programari.
Així que haurem d'especificar si volem utilitzar Docker o Singularity per provar el pipeline.

La part `--outdir <OUTDIR>` es refereix al directori on Nextflow escriurà les sortides del pipeline.
Hem de proporcionar un nom per a ell, que podem inventar.
Si encara no existeix, Nextflow el crearà per a nosaltres en temps d'execució.

Passant a la secció després del bloc de comentaris, el perfil de prova ens mostra el que s'ha preconfigurat per a les proves: més notablement, el paràmetre `input` ja està configurat per apuntar a un conjunt de dades de prova, així que no necessitem proporcionar les nostres pròpies dades.
Si seguiu l'enllaç a l'entrada preconfigurada, veureu que és un fitxer csv que conté identificadors de mostra i camins de fitxer per a diverses mostres experimentals.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Això s'anomena samplesheet, i és la forma més comuna d'entrada als pipelines nf-core.
No us preocupeu si no esteu familiaritzats amb els formats i tipus de dades, no és important per al que segueix.

Ara tenim tot el que necessitem per provar el pipeline.

### 2.2. Executar el pipeline

Tal com s'ha indicat anteriorment, podem utilitzar l'exemple de comanda de prova gairebé tal com és; només hem d'especificar quin sistema d'empaquetament de programari volem utilitzar i quin nom donar al directori de sortida.
Aquí utilitzarem Docker per al sistema de contenidors i `demo-results`, respectivament.

Amb això, podem executar la comanda de prova:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
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
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : docker,test
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
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     [100%] 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) [100%] 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   [100%] 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          [100%] 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Si la vostra sortida coincideix amb això, felicitats! Acabeu d'executar el vostre primer pipeline nf-core.

Notareu que hi ha molta més sortida a la consola que quan executeu un pipeline Nextflow bàsic.
Hi ha una capçalera que inclou un resum de la versió del pipeline, entrades i sortides, i alguns elements de configuració.

!!! info "Info"

    La vostra sortida mostrarà diferents marques de temps, noms d'execució i camins de fitxer, però l'estructura general i l'execució del procés haurien de ser similars.

Fixeu-vos en la línia prop de la part superior de la sortida:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Això us indica quina revisió del pipeline s'ha utilitzat.
Com que no hem especificat cap versió, Nextflow ha utilitzat el darrer commit a `master`.
Per a execucions reproduïbles, hauríeu de fixar una versió específica amb el flag `-r`:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile docker,test --outdir demo-results
```

Això garanteix que s'utilitzi sempre el mateix codi del pipeline, independentment de nous commits o versions.
En aquesta formació ometem `-r` per simplicitat, però en producció sempre hauríeu d'especificar-lo.

Passant a la sortida d'execució, donem una ullada a les línies que ens diuen quins processos s'han executat:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     [100%] 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) [100%] 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   [100%] 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          [100%] 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Això ens diu que s'han executat quatre processos, corresponents a les quatre eines mostrades a la pàgina de documentació del pipeline al lloc web nf-core: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` i `COWPY`.

Els noms complets dels processos tal com es mostren aquí, com ara `NFCORE_DEMO:DEMO:MULTIQC`, són més llargs del que potser heu vist al material introductori Hello Nextflow.
Aquests inclouen els noms dels seus workflows pare i reflecteixen la modularitat del codi del pipeline.
Entrarem en més detall sobre això a la Part 2 d'aquest curs.

### 2.3. Examinar les sortides del pipeline

Finalment, donem una ullada al directori `demo-results` produït pel pipeline.

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

Això pot semblar molt.
Per aprendre més sobre les sortides del pipeline `nf-core/demo`, consulteu la seva [pàgina de documentació](https://nf-co.re/demo/1.2.0/docs/output/).

En aquesta etapa, el que és important observar és que els resultats estan organitzats per mòdul, i hi ha addicionalment un directori anomenat `pipeline_info` que conté diversos informes amb marca de temps sobre l'execució del pipeline.

Per exemple, el fitxer `execution_timeline_*` us mostra quins processos s'han executat, en quin ordre i quant de temps han trigat a executar-se:

![informe de línia de temps d'execució](./img/execution_timeline.png)

!!! info "Info"

    Aquí les tasques no s'han executat en paral·lel perquè estem executant en una màquina minimalista a Github Codespaces.
    Per veure-les executar-se en paral·lel, proveu d'augmentar l'assignació de CPU del vostre codespace i els límits de recursos a la configuració de prova.

Aquests informes es generen automàticament per a tots els pipelines nf-core.

### Conclusió

Sabeu com executar un pipeline nf-core utilitzant el seu perfil de prova integrat i on trobar les seves sortides.

### Què segueix?

Apreneu com configurar el pipeline per personalitzar la seva execució.

---

## 3. Configurar l'execució del pipeline

Tal com s'explica a [Hello Config](../hello_nextflow/06_hello_config.md), volem poder canviar les dades sobre les quals s'executarà el nostre pipeline i com s'executarà sense modificar el codi del pipeline en si.
Per a aquest fi, Nextflow admet múltiples maneres de controlar la configuració del pipeline, cosa que pot resultar una mica aclaparadora.

El projecte nf-core especifica convencions per organitzar els elements de configuració, distingint dos tipus de configuració al nivell superior: **paràmetres del pipeline** i **configuració** en sentit estricte.

- Els **paràmetres del pipeline** (establerts mitjançant el sistema `params`) típicament inclouen coses com fitxers d'entrada, flags de comportament d'eines i paràmetres d'anàlisi.
- La **configuració** en sentit estricte fa referència a la logística de com s'executa el pipeline, és a dir, l'executor, les assignacions de recursos informàtics, etc.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/params_vs_config.excalidraw.svg"
</figure>

Comencem tractant els paràmetres del pipeline i després veurem la configuració en sentit estricte.

### 3.1. Paràmetres del pipeline

Per a tots els pipelines nf-core, podeu obtenir una llista completa dels paràmetres del pipeline directament des de la línia de comandes utilitzant el flag `--help`, que és en si mateix un paràmetre del pipeline.

#### 3.1.1. Obtenir la llista de paràmetres amb `--help`

Executeu la comanda d'ajuda per al pipeline de demostració:

```bash
nextflow run nf-core/demo --help
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.

      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
    !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Com podeu veure, la sortida agrupa els paràmetres en categories (opcions d'entrada/sortida, opcions de genoma de referència, etc.) amb tipus i descripcions per a cadascun.

Aquesta categorització ve determinada per un fitxer d'esquema, que es tracta més endavant.
En pipelines Nextflow simples, `--help` només funciona si el desenvolupador l'ha implementat manualment.

!!! tip "Consell"

    Utilitzeu `--help --show_hidden` per veure paràmetres addicionals que estan ocults per defecte, com ara `--publish_dir_mode` o `--monochrome_logs`.

#### 3.1.2. Establir valors de paràmetres

Tal com es tracta a [Hello Config](../hello_nextflow/06_hello_config.md), podeu establir valors de paràmetres a la línia de comandes amb `--nom_parametre` o recollir un conjunt de paràmetres en un fitxer YAML i passar-lo amb `-params-file`.
Tots dos enfocaments funcionen de la mateixa manera amb els pipelines nf-core.

Per exemple, per ometre el pas de retallada, volem establir el paràmetre booleà `skip_trim` a `true`.
Al vostre directori de treball hi ha un fitxer de paràmetres anomenat `my_params.yml` amb aquest valor ja configurat:

```yaml title="my_params.yml"
skip_trim: true
```

Passeu-lo amb `-params-file`:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


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
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) [100%] 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               [100%] 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      [100%] 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El procés `SEQTK_TRIM` ja no apareix a la sortida.

!!! warning "Advertència: limitacions importants sobre les entrades de paràmetres"

    **Establir paràmetres booleans a la línia de comandes**

    A partir de la versió 26.04 de Nextflow, tots els valors subministrats a la línia de comandes es tracten com a strings.
    Per a un paràmetre booleà com `skip_trim`, passar-lo com a flag simple (`--skip_trim`) o com `--skip_trim true` s'avalua com el **string** `"true"`, cosa que fa fallar la validació de l'esquema:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Per establir un paràmetre booleà a un valor genuí `true`/`false`, utilitzeu un `-params-file` tal com es mostra més amunt, o establiu-lo en un fitxer de configuració.
    Els paràmetres de tipus string, integer i file-path no es veuen afectats i es poden continuar establint directament a la línia de comandes.
    Aquest curs utilitza aquest patró per a tots els paràmetres booleans.

    **Utilitzar fitxers de configuració personalitzats**

    Tot i que tècnicament és possible establir paràmetres del pipeline en un fitxer de configuració personalitzat passat amb `-c`, és possible que no sobreescrigui els valors per defecte ja establerts al `nextflow.config` propi del pipeline, depenent de les regles de precedència de configuració de Nextflow.
    Utilitzar `--nom_parametre` a la línia de comandes o `-params-file` és més fiable, ja que aquests sempre tenen prioritat.

    Com a regla general: si apareix a la sortida de `--help`, establiu-lo mitjançant la línia de comandes o un fitxer de paràmetres en lloc d'un fitxer de configuració.

#### 3.1.3. Validació de paràmetres

Curiositat: la comanda `--help` funciona per a tots els pipelines nf-core perquè el projecte nf-core requereix que els desenvolupadors defineixin formalment tots els paràmetres del pipeline en un fitxer d'esquema JSON (`nextflow_schema.json`).
Aquest esquema registra el tipus, la descripció, el valor per defecte i l'agrupació de cada paràmetre.

A més de generar la sortida de `--help`, el fitxer d'esquema també permet la validació automatitzada en el moment del llançament.
Això significa que Nextflow pot comprovar que cada paràmetre que passeu existeix i té un valor adequat (del tipus adequat, dins del rang de valors permesos, etc.).

Ho tractem amb més detall a [Part 5: Input Validation](05_input_validation.md), però ja podeu veure-ho en acció donant al pipeline de demostració alguna entrada de paràmetres no vàlida.

##### 3.1.3.1. Paràmetres no reconeguts

Proveu de passar un paràmetre que no existeix:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --foobar "invalid"
```

La sortida de la consola inclou un avís:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

El pipeline continua executant-se, però l'avís us alerta immediatament que `--foobar` no és un paràmetre reconegut.
Això vol cridar la vostra atenció sobre errors tipogràfics que no trenquen l'execució, com ara `--outDir` en lloc de `--outdir`, cosa que us pot ajudar a evitar malgastar temps i recursos de còmput.

##### 3.1.3.2. Valors de paràmetres no vàlids

La validació també comprova els **valors** dels paràmetres.
El paràmetre `--skip_trim` és un flag booleà, de manera que passar un valor de tipus string fa que el pipeline falli immediatament:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]
```

El pipeline s'atura abans que s'executi cap procés, estalviant-vos una execució fallida o incorrecta.
Tal com s'indica a la secció 3.1.2, els paràmetres booleans s'han d'establir a un valor genuí `true`/`false` en un fitxer de paràmetres en lloc de passar-los a la línia de comandes, ja que els valors de la línia de comandes es tracten com a strings.

#### 3.1.4. Validació d'entrada

La mateixa lògica de validació també es pot utilitzar per comprovar la validesa dels fitxers d'entrada.
Per exemple, si un pipeline espera un samplesheet com a entrada de dades principal (que és el cas de molts, si no la majoria, dels pipelines nf-core), el desenvolupador pot proporcionar un esquema d'entrada (diferent de l'esquema de paràmetres) que descrigui com s'ha d'estructurar el fitxer d'entrada.

Llavors, en temps d'execució, Nextflow pot comprovar que el fitxer d'entrada proporcionat és vàlid.

També ho tractem amb més detall a [Part 5: Input Validation](05_input_validation.md), però ja podeu veure-ho en acció donant al pipeline de demostració un samplesheet d'entrada no vàlid.

El pipeline `nf-core/demo` espera un fitxer CSV amb les columnes `sample`, `fastq_1` i `fastq_2`.
Això es defineix en un fitxer d'esquema (`assets/schema_input.json`) que especifica l'estructura esperada, els tipus de columnes i les restriccions.

??? abstract "Fitxer d'esquema per a les entrades"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

L'esquema especifica que `sample` i `fastq_1` són obligatoris, mentre que `fastq_2` és opcional (admetent tant dades paired-end com single-end).
Els camins de fitxer es validen per existència i patró d'extensió.

Per demostrar-ho, al vostre directori de treball hi ha un samplesheet mal format anomenat `malformed_samplesheet.csv`:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

A aquest samplesheet li falta la columna obligatòria `fastq_1` i té un camí de fitxer inexistent a `fastq_2`.

Executeu el pipeline de demostració utilitzant `malformed_samplesheet.csv` com a entrada:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory
       '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces
       and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1
```

Com podeu veure, el pipeline falla immediatament i informa de **tots** els errors de validació alhora.
nf-schema no s'atura al primer error: recull tots els problemes i els llista junts, de manera que podeu corregir-ho tot d'una vegada en lloc de descobrir els problemes un per un.

Cada error identifica l'entrada i el camp exactes que han causat el problema, de manera que podeu corregir el vostre samplesheet i tornar a llançar el pipeline amb la confiança que no fallarà en algun punt posterior quan Nextflow intenti accedir al camí del fitxer.

Per als desenvolupadors, tot això es tracta amb més detall a la [Part 5](./05_input_validation.md) d'aquest curs.

### 3.2. Configuració

La configuració en sentit estricte controla **com** s'executa el pipeline: assignació de recursos, arguments específics d'eines, on s'executen les tasques i quin sistema d'empaquetament de programari s'utilitza.

Els pipelines nf-core inclouen configuració per defecte a `nextflow.config` i al directori `conf/`.
Abans de sobreescriure res, és útil saber on es troben els valors per defecte.

Ja heu vist a la secció 2.1 que el codi font del pipeline es troba a `$NXF_HOME/assets`.
Utilitzant l'enllaç simbòlic `pipelines` de la secció 1.2.4, llisteu els fitxers de configuració per veure què hi ha disponible:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/nfcore_config_files.excalidraw.svg"
</figure>

Els fitxers de configuració més importants són:

- **`conf/base.config`**: Defineix etiquetes de recursos (`process_low`, `process_medium`, `process_high`) que assignen CPUs, memòria i temps als processos. Quan veieu un procés que utilitza més recursos dels esperats, aquí és d'on provenen aquests valors per defecte.
- **`conf/modules.config`**: Estableix arguments d'eines per procés (`ext.args`) i configuració de publicació de sortides (`publishDir`). Obriu aquest fitxer per veure quins arguments rep cada eina per defecte.
- **`conf/test.config`**: El perfil de prova que heu utilitzat a la secció 2.1, que limita els recursos mitjançant `resourceLimits` i estableix un samplesheet de prova. S'activa amb `-profile test`.
  També hi ha un `conf/test_full.config` per executar amb un conjunt de dades de prova de mida completa, útil per a benchmarking.

El `nextflow.config` central carrega tots els anteriors i estableix els valors per defecte adequats per a tot.

Si voleu modificar qualsevol dels paràmetres especificats en aquests fitxers, no modifiqueu cap d'ells directament.
En canvi, creeu el vostre propi fitxer de configuració i passeu-lo amb `-c`.
Els valors que especifiqueu sobreescriuran els valors per defecte establerts en aquells altres fitxers.

Practiquem-ho.

#### 3.2.1. Personalitzar els recursos dels processos i els arguments de les eines

Els mòduls nf-core admeten dos tipus comuns de sobreescriptura de configuració: **assignació de recursos** (CPUs, memòria, temps) i **arguments d'eines** mitjançant `ext.args`.

Moltes eines de línia de comandes tenen arguments que no s'utilitzen prou freqüentment com per exposar-los com a paràmetres del pipeline.
La convenció `ext.args` us permet passar aquests arguments a l'eina subjacent mitjançant un fitxer de configuració.

El fitxer `custom.config` proporcionat al vostre directori de treball demostra tots dos tipus de sobreescriptura:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

El primer bloc sobreescriu l'assignació de recursos de `FASTQC`.
Per defecte, `FASTQC` utilitza l'etiqueta `process_medium` de `base.config`, que assigna 6 CPUs i 36 GB de memòria; aquí ho limitem a 2 CPUs i 4 GB.

El segon bloc passa un argument addicional a `SEQTK_TRIM` mitjançant `ext.args`.
El flag `-b 5` indica a `seqtk trimfq` que retalli 5 bases del principi de cada lectura a més de la retallada per qualitat.

Executeu el pipeline amb aquesta configuració:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-custom -c custom.config
```

??? success "Sortida de la comanda"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El flag `-c` afegeix la vostra configuració a sobre de la configuració integrada del pipeline.

Per verificar que la sobreescriptura de `ext.args` ha tingut efecte, trobeu el hash del directori de treball de `SEQTK_TRIM` a la sortida de l'execució (per exemple, `work/17/428668...`) i comproveu el fitxer `.command.sh` que hi ha dins:

```bash
cat work/17/428668/.command.sh
```

??? success "Sortida de la comanda"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Hauríeu de veure `-b 5` a la comanda `seqtk trimfq`.

Una cosa important a saber sobre `ext.args`: si un mòdul ja té un valor per defecte establert, el vostre valor el **reemplaçarà completament** en lloc d'afegir-s'hi.
Per exemple, `FASTQC` té `ext.args = '--quiet'` establert per defecte a `conf/modules.config`:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Si establiu `ext.args = '--kmers 8'` per a `FASTQC`, el flag `--quiet` ja no s'aplicarà.
Per mantenir tots dos, establiu `ext.args = '--quiet --kmers 8'`.

Sempre hauríeu de comprovar la configuració per defecte d'un mòdul abans de sobreescriure `ext.args`.

### Conclusió

Sabeu com obtenir ajuda d'un pipeline nf-core, establir paràmetres i entendre com es validen, i personalitzar la configuració mitjançant fitxers de configuració.

### Què segueix?

Si simplement voleu executar pipelines nf-core, ja heu acabat!

Si voleu aprendre a desenvolupar els vostres propis pipelines seguint els estàndards nf-core, feu una pausa i passeu a la Part 2 quan estigueu preparats. Aprendreu a crear el vostre propi pipeline compatible amb nf-core utilitzant les eines basades en la plantilla nf-core.
