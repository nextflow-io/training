# Part 1: Executar un pipeline de demostració

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta primera part del curs de formació Hello nf-core, us mostrem com trobar i provar un pipeline nf-core, entendre com s'organitza el codi i reconèixer com difereix del codi Nextflow simple tal com es mostra a [Hello Nextflow](../hello_nextflow/index.md).

Utilitzarem un pipeline anomenat nf-core/demo que és mantingut pel projecte nf-core com a part del seu inventari de pipelines per demostrar l'estructura del codi i les operacions de les eines.

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

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

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

Tornem al terminal i executem el següent:

```bash
nextflow pull nf-core/demo
```

??? success "Sortida de la comanda"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow fa un `pull` del codi del pipeline, és a dir, descarrega el repositori complet a la vostra unitat local.

Per ser clars, podeu fer això amb qualsevol pipeline Nextflow que estigui configurat adequadament a GitHub, no només amb pipelines nf-core.
No obstant això, nf-core és la col·lecció de codi obert més gran de pipelines Nextflow.

Podeu fer que Nextflow us doni una llista dels pipelines que heu recuperat d'aquesta manera:

```bash
nextflow list
```

??? success "Sortida de la comanda"

    ```console
    nf-core/demo
    ```

Notareu que els fitxers no estan al vostre directori de treball actual.
Per defecte, Nextflow els desa a `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Directory contents"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note "Nota"

    El camí complet pot diferir al vostre sistema si no esteu utilitzant el nostre entorn de formació.

Nextflow manté el codi font descarregat intencionadament 'fora del camí' amb el principi que aquests pipelines s'haurien d'utilitzar més com a biblioteques que com a codi amb el qual interactuaríeu directament.

No obstant això, per als propòsits d'aquesta formació, volem poder explorar i veure què hi ha allà.
Així que per facilitar-ho, creem un enllaç simbòlic a aquesta ubicació des del nostre directori de treball actual.

```bash
ln -s $NXF_HOME/assets pipelines
```

Això crea una drecera que facilita l'exploració del codi que acabem de descarregar.

```bash
tree -L 2 pipelines
```

```console title="Directory contents"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Ara podem mirar més fàcilment el codi font segons sigui necessari.

Però primer, provem d'executar el nostre primer pipeline nf-core!

### Conclusió

Ara sabeu com trobar un pipeline a través del lloc web nf-core i recuperar una còpia local del codi font.

### Què segueix?

Apreneu com provar un pipeline nf-core amb un esforç mínim.

---

## 2. Provar el pipeline amb el seu perfil de prova

Convenientment, cada pipeline nf-core ve amb un perfil de prova.
Aquest és un conjunt mínim de paràmetres de configuració perquè el pipeline s'executi utilitzant un petit conjunt de dades de prova allotjat al repositori [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
És una manera excel·lent de provar ràpidament un pipeline a petita escala.

!!! note "Nota"

    El sistema de perfils de configuració de Nextflow us permet canviar fàcilment entre diferents motors de contenidors o entorns d'execució.
    Per a més detalls, consulteu [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar el perfil de prova

És una bona pràctica comprovar què especifica el perfil de prova d'un pipeline abans d'executar-lo.
El perfil `test` per a `nf-core/demo` es troba al fitxer de configuració `conf/test.config` i es mostra a continuació.

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
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

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

!!! note "Nota"

    No us preocupeu si no esteu familiaritzats amb els formats i tipus de dades, no és important per al que segueix.

Així que això confirma que tenim tot el que necessitem per provar el pipeline.

### 2.2. Executar el pipeline

Decidim utilitzar Docker per al sistema de contenidors i `demo-results` com a directori de sortida, i estem preparats per executar la comanda de prova:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Si la vostra sortida coincideix amb això, felicitats! Acabeu d'executar el vostre primer pipeline nf-core.

Notareu que hi ha molta més sortida a la consola que quan executeu un pipeline Nextflow bàsic.
Hi ha una capçalera que inclou un resum de la versió del pipeline, entrades i sortides, i alguns elements de configuració.

!!! note "Nota"

    La vostra sortida mostrarà diferents marques de temps, noms d'execució i camins de fitxer, però l'estructura general i l'execució del procés haurien de ser similars.

Passant a la sortida d'execució, donem una ullada a les línies que ens diuen quins processos s'han executat:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Això ens diu que s'han executat tres processos, corresponents a les tres eines mostrades a la pàgina de documentació del pipeline al lloc web nf-core: FASTQC, SEQTK_TRIM i MULTIQC.

Els noms complets dels processos tal com es mostren aquí, com ara `NFCORE_DEMO:DEMO:MULTIQC`, són més llargs del que potser heu vist al material introductori Hello Nextflow.
Aquests inclouen els noms dels seus workflows pare i reflecteixen la modularitat del codi del pipeline.
Entrarem en més detall sobre això d'aquí a una estona.

### 2.3. Examinar les sortides del pipeline

Finalment, donem una ullada al directori `demo-results` produït pel pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Contingut del directori"

    ```console
    demo-results
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
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

Això pot semblar molt.
Per aprendre més sobre les sortides del pipeline `nf-core/demo`, consulteu la seva [pàgina de documentació](https://nf-co.re/demo/1.0.2/docs/output/).

En aquesta etapa, el que és important observar és que els resultats estan organitzats per mòdul, i hi ha addicionalment un directori anomenat `pipeline_info` que conté diversos informes amb marca de temps sobre l'execució del pipeline.

Per exemple, el fitxer `execution_timeline_*` us mostra quins processos s'han executat, en quin ordre i quant de temps han trigat a executar-se:

![informe de línia de temps d'execució](./img/execution_timeline.png)

!!! note "Nota"

    Aquí les tasques no s'han executat en paral·lel perquè estem executant en una màquina minimalista a Github Codespaces.
    Per veure-les executar-se en paral·lel, proveu d'augmentar l'assignació de CPU del vostre codespace i els límits de recursos a la configuració de prova.

Aquests informes es generen automàticament per a tots els pipelines nf-core.

### Conclusió

Sabeu com executar un pipeline nf-core utilitzant el seu perfil de prova integrat i on trobar les seves sortides.

### Què segueix?

Apreneu com s'organitza el codi del pipeline.

---

## 3. Examinar l'estructura del codi del pipeline

Ara que hem executat amb èxit el pipeline com a usuaris, canviem la nostra perspectiva per veure com estan estructurats internament els pipelines nf-core.

El projecte nf-core aplica directrius estrictes sobre com s'estructuren els pipelines, i sobre com s'organitza, configura i documenta el codi.
Entendre com s'organitza tot això és el primer pas cap al desenvolupament dels vostres propis pipelines compatibles amb nf-core, que abordarem a la Part 2 d'aquest curs.

Donem una ullada a com s'organitza el codi del pipeline al repositori `nf-core/demo`, utilitzant l'enllaç simbòlic `pipelines` que hem creat anteriorment.

Podeu utilitzar `tree` o utilitzar l'explorador de fitxers per trobar i obrir el directori `nf-core/demo`.

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
    ```

Hi ha molt en marxa allà, així que abordarem això pas a pas.

Primer, observem que al nivell superior, podeu trobar un fitxer README amb informació resumida, així com fitxers accessoris que resumeixen informació del projecte com ara llicència, directrius de contribució, citació i codi de conducta.
La documentació detallada del pipeline es troba al directori `docs`.
Tot aquest contingut s'utilitza per generar les pàgines web al lloc web nf-core de manera programàtica, així que sempre estan actualitzades amb el codi.

Ara, per a la resta, dividirem la nostra exploració en tres etapes:

1. Components del codi del pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuració del pipeline
3. Entrades i validació

Comencem amb els components del codi del pipeline.
Ens centrarem en la jerarquia de fitxers i l'organització estructural, en lloc d'aprofundir en el codi dins dels fitxers individuals.

### 3.1. Components del codi del pipeline

L'organització estàndard del codi del pipeline nf-core segueix una estructura modular que està dissenyada per maximitzar la reutilització del codi, tal com es va introduir a [Hello Modules](../hello_nextflow/04_hello_modules.md), Part 4 del curs [Hello Nextflow](../hello_nextflow/index.md), encara que a l'estil nf-core, això s'implementa amb una mica de complexitat addicional.
Específicament, els pipelines nf-core fan un ús abundant de subworkflows, és a dir, scripts de workflow que són importats per un workflow pare.

Això pot sonar una mica abstracte, així que donem una ullada a com s'utilitza això a la pràctica al pipeline `nf-core/demo`.

!!! note "Nota"

    No repassarem el codi real de _com_ es connecten aquests components modulars, perquè hi ha una certa complexitat addicional associada amb l'ús de subworkflows que pot ser confusa, i entendre això no és necessari en aquesta etapa de la formació.
    Per ara, ens centrarem en l'organització general i la lògica.

#### 3.1.1. Visió general

Així és com es veuen les relacions entre els components de codi rellevants per al pipeline `nf-core/demo`:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

Hi ha un script anomenat _punt d'entrada_ anomenat `main.nf`, que actua com a embolcall per a dos tipus de workflows niats: el workflow que conté la lògica d'anàlisi real, ubicat a `workflows/` i anomenat `demo.nf`, i un conjunt de workflows de manteniment ubicats a `subworkflows/`.
El workflow `demo.nf` crida **mòduls** ubicats a `modules/`; aquests contenen els **processos** que realitzaran els passos d'anàlisi reals.

!!! note "Nota"

    Els subworkflows no es limiten a funcions de manteniment, i poden fer ús de mòduls de procés.

    El pipeline `nf-core/demo` mostrat aquí resulta ser del costat més simple de l'espectre, però altres pipelines nf-core (com ara `nf-core/rnaseq`) utilitzen subworkflows que estan involucrats en l'anàlisi real.

Ara, revisem aquests components per torn.

#### 3.1.2. L'script de punt d'entrada: `main.nf`

L'script `main.nf` és el punt d'entrada des del qual Nextflow comença quan executem `nextflow run nf-core/demo`.
Això significa que quan executeu `nextflow run nf-core/demo` per executar el pipeline, Nextflow troba i executa automàticament l'script `main.nf`.
Això funciona per a qualsevol pipeline Nextflow que segueixi aquesta nomenclatura i estructura convencionals, no només per als pipelines nf-core.

L'ús d'un script de punt d'entrada facilita l'execució de subworkflows de 'manteniment' estandarditzats abans i després que s'executi l'script d'anàlisi real.
Repassarem aquests després d'haver revisat el workflow d'anàlisi real i els seus mòduls.

#### 3.1.3. L'script d'anàlisi: `workflows/demo.nf`

El workflow `workflows/demo.nf` és on s'emmagatzema la lògica central del pipeline.
Està estructurat de manera molt similar a un workflow Nextflow normal, excepte que està dissenyat per ser cridat des d'un workflow pare, cosa que requereix algunes característiques addicionals.
Cobrirem les diferències rellevants a la següent part d'aquest curs, quan abordem la conversió del simple pipeline Hello de Hello Nextflow a una forma compatible amb nf-core.

El workflow `demo.nf` crida **mòduls** ubicats a `modules/`, que revisarem a continuació.

!!! note "Nota"

    Alguns workflows d'anàlisi nf-core mostren nivells addicionals de niament cridant subworkflows de nivell inferior.
    Això s'utilitza principalment per embolcallar dos o més mòduls que s'utilitzen comunament junts en segments de pipeline fàcilment reutilitzables.
    Podeu veure alguns exemples navegant pels [subworkflows nf-core](https://nf-co.re/subworkflows/) disponibles al lloc web nf-core.

    Quan l'script d'anàlisi utilitza subworkflows, aquests s'emmagatzemen al directori `subworkflows/`.

#### 3.1.4. Els mòduls

Els mòduls són on viu el codi del procés, tal com es descriu a la [Part 4 del curs de formació Hello Nextflow](../hello_nextflow/04_hello_modules.md).

Al projecte nf-core, els mòduls s'organitzen utilitzant una estructura niada de múltiples nivells que reflecteix tant el seu origen com el seu contingut.
Al nivell superior, els mòduls es diferencien com a `nf-core` o `local` (no part del projecte nf-core), i després es col·loquen en un directori anomenat segons l'eina o eines que embolcallen.
Si l'eina pertany a un toolkit (és a dir, un paquet que conté múltiples eines) llavors hi ha un nivell de directori intermedi anomenat segons el toolkit.

Podeu veure això aplicat a la pràctica als mòduls del pipeline `nf-core/demo`:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Contingut del directori"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Aquí veieu que els mòduls `fastqc` i `multiqc` es troben al nivell superior dins dels mòduls `nf-core`, mentre que el mòdul `trim` es troba sota el toolkit al qual pertany, `seqtk`.
En aquest cas no hi ha mòduls `local`.

El fitxer de codi del mòdul que descriu el procés sempre s'anomena `main.nf`, i està acompanyat de proves i fitxers `.yml` que ignorarem per ara.

Presos conjuntament, el workflow de punt d'entrada, el workflow d'anàlisi i els mòduls són suficients per executar les parts 'interessants' del pipeline.
No obstant això, sabem que també hi ha subworkflows de manteniment allà, així que mirem-los ara.

#### 3.1.5. Els subworkflows de manteniment

Com els mòduls, els subworkflows es diferencien en directoris `local` i `nf-core`, i cada subworkflow té la seva pròpia estructura de directori niat amb el seu propi script `main.nf`, proves i fitxer `.yml`.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Contingut del directori"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Com s'ha indicat anteriorment, el pipeline `nf-core/demo` no inclou cap subworkflow específic d'anàlisi, així que tots els subworkflows que veiem aquí són els anomenats workflows de 'manteniment' o 'utilitat', tal com es denota pel prefix `utils_` als seus noms.
Aquests subworkflows són el que produeix la capçalera nf-core elegant a la sortida de la consola, entre altres funcions accessòries.

!!! tip "Consell"

    A part del seu patró de nomenclatura, una altra indicació que aquests subworkflows no realitzen cap funció realment relacionada amb l'anàlisi és que no criden cap procés en absolut.

Això completa el resum dels components de codi bàsics que constitueixen el pipeline `nf-core/demo`.
Ara donem una ullada als elements restants que hauríeu de conèixer una mica abans d'endinsar-vos en el desenvolupament: configuració del pipeline i validació d'entrada.

### 3.2. Configuració del pipeline

Heu après anteriorment que Nextflow ofereix moltes opcions per configurar l'execució del pipeline, ja sigui en termes d'entrades i paràmetres, recursos informàtics i altres aspectes de l'orquestració.
El projecte nf-core aplica directrius altament estandarditzades per a la configuració del pipeline que tenen com a objectiu construir sobre les opcions de personalització flexibles de Nextflow d'una manera que proporcioni una major coherència i mantenibilitat entre pipelines.

El fitxer de configuració central `nextflow.config` s'utilitza per establir valors per defecte per a paràmetres i altres opcions de configuració.
La majoria d'aquestes opcions de configuració s'apliquen per defecte mentre que d'altres (per exemple, perfils de dependències de programari) s'inclouen com a perfils opcionals.

Hi ha diversos fitxers de configuració addicionals que s'emmagatzemen a la carpeta `conf` i que es poden afegir a la configuració per defecte o opcionalment com a perfils:

- `base.config`: Un fitxer de configuració 'en blanc', apropiat per a l'ús general a la majoria d'entorns informàtics d'alt rendiment. Això defineix grups amplis d'ús de recursos, per exemple, que són convenients d'aplicar als mòduls.
- `modules.config`: Directives i arguments de mòdul addicionals.
- `test.config`: Un perfil per executar el pipeline amb dades de prova mínimes, que hem utilitzat quan hem executat el pipeline de demostració.
- `test_full.config`: Un perfil per executar el pipeline amb un conjunt de dades de prova de mida completa.

Tocarem alguns d'aquests fitxers més endavant al curs.

### 3.3. Entrades i validació

Com hem observat anteriorment, quan hem examinat el perfil de prova del pipeline `nf-core/demo`, està dissenyat per prendre com a entrada un samplesheet que conté camins de fitxer i identificadors de mostra.
Els camins de fitxer enllaçats a dades reals ubicades al repositori `nf-core/test-datasets`.

També es proporciona un exemple de samplesheet al directori `assets`, encara que els camins en aquest no són reals.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Aquest samplesheet en particular és bastant simple, però alguns pipelines s'executen amb samplesheets que són més complexos, amb molt més metadades associades amb les entrades primàries.

Malauradament, com que aquests fitxers poden ser difícils de comprovar a ull, el format inadequat de les dades d'entrada és una font molt comuna de fallades del pipeline.
Un problema relacionat és quan els paràmetres es proporcionen incorrectament.

La solució a aquests problemes és executar comprovacions de validació automatitzades en tots els fitxers d'entrada per assegurar que contenen els tipus d'informació esperats, formatats correctament, i en paràmetres per assegurar que són del tipus esperat.
Això s'anomena validació d'entrada, i idealment s'hauria de fer _abans_ d'intentar executar un pipeline, en lloc d'esperar que el pipeline falli per esbrinar que hi havia un problema amb les entrades.

Igual que per a la configuració, el projecte nf-core té opinions molt fortes sobre la validació d'entrada, i recomana l'ús del [plugin nf-schema](https://nextflow-io.github.io/nf-schema/latest/), un plugin Nextflow que proporciona capacitats de validació completes per als pipelines Nextflow.

Cobrirem aquest tema amb més detall a la Part 5 d'aquest curs.
Per ara, només sigueu conscients que hi ha dos fitxers JSON proporcionats per a aquest propòsit, `nextflow_schema.json` i `assets/schema_input.json`.

El `nextflow_schema.json` és un fitxer utilitzat per emmagatzemar informació sobre els paràmetres del pipeline incloent tipus, descripció i text d'ajuda en un format llegible per màquina.
Això s'utilitza per a diversos propòsits, incloent validació automatitzada de paràmetres, generació de text d'ajuda i renderització de formularis de paràmetres interactius en interfícies d'usuari.

El `schema_input.json` és un fitxer utilitzat per definir l'estructura del samplesheet d'entrada.
Cada columna pot tenir un tipus, patró, descripció i text d'ajuda en un format llegible per màquina.
L'esquema s'utilitza per a diversos propòsits, incloent validació automatitzada i proporcionar missatges d'error útils.

### Conclusió

Sabeu quins són els components principals d'un pipeline nf-core i com s'organitza el codi; on es troben els elements principals de configuració; i sou conscients de per a què serveix la validació d'entrada.

### Què segueix?

Feu una pausa! Això ha estat molt. Quan us sentiu refrescats i preparats, passeu a la següent secció per aplicar el que heu après per escriure un pipeline compatible amb nf-core.

!!! tip "Consell"

    Si voleu aprendre com compondre workflows amb subworkflows abans de passar a la següent part, consulteu la [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest.
