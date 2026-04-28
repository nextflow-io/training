# Part 1: Executar un pipeline de demostració

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta primera part del curs de formació Hello nf-core, us mostrem com trobar i provar un pipeline nf-core, configurar i personalitzar la seva execució per a les vostres necessitats, i entendre com la validació d'entrada protegeix contra errors comuns.

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

#### 1.2.1. Utilitzar `nextflow pull`

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

#### 1.2.3. Trobar els vostres pipelines a `$NXF_HOME/assets/`

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

#### 1.2.4. Crear un enllaç simbòlic per accedir fàcilment al codi font

No examinarem el codi en detall, però fem-hi una ullada ràpida per tenir una idea de com és l'organització general.

Per facilitar la navegació pel codi font del pipeline, creeu un enllaç simbòlic al directori d'assets:

```bash
ln -s $NXF_HOME/assets pipelines
```

Això crea una drecera perquè pugueu explorar el codi amb `tree -L 2 pipelines` o obrir fitxers directament.

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

!!! note "Nota"

    El sistema de perfils de configuració de Nextflow us permet canviar fàcilment entre diferents motors de contenidors o entorns d'execució.
    Per a més detalls, consulteu [Hello Nextflow Part 6: Configuration](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar el perfil de prova

És una bona pràctica comprovar què especifica el perfil de prova d'un pipeline abans d'executar-lo.
El perfil `test` per a `nf-core/demo` es troba al fitxer de configuració `conf/test.config`.
Podeu trobar-lo localment dins del codi font del pipeline que `nextflow pull` ha descarregat:

```bash
code $NXF_HOME/assets/nf-core/demo/conf/test.config
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
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Dades d'entrada
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
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
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

Fixeu-vos en la línia prop de la part superior de la sortida:

```console
Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]
```

Això us indica quina revisió del pipeline s'ha utilitzat.
Com que no hem especificat cap versió, Nextflow ha utilitzat el darrer commit a `master`.
Per a execucions reproduïbles, hauríeu de fixar una versió específica amb el flag `-r`:

```bash
nextflow run nf-core/demo -r 1.1.0 -profile docker,test --outdir demo-results
```

Això garanteix que s'utilitzi sempre el mateix codi del pipeline, independentment de nous commits o versions.
En aquesta formació ometem `-r` per simplicitat, però en producció sempre hauríeu d'especificar-lo.

Passant a la sortida d'execució, donem una ullada a les línies que ens diuen quins processos s'han executat:

```console
executor >  local (7)
[ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Això ens diu que s'han executat tres processos, corresponents a les tres eines mostrades a la pàgina de documentació del pipeline al lloc web nf-core: FASTQC, SEQTK_TRIM i MULTIQC.

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
Per aprendre més sobre les sortides del pipeline `nf-core/demo`, consulteu la seva [pàgina de documentació](https://nf-co.re/demo/1.1.0/docs/output/).

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
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [run_name] DSL2 - revision: 45904cb9d1 [master]

    ----------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
    ----------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

    Input/output options
      --input                       [string]           Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string]           The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string]           Email address for completion summary.
      --multiqc_title               [string]           MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string]           Name of iGenomes reference.
      --fasta                       [string]           Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean]          Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]           Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string]  Display the help message.
      --help_full                   [boolean]          Display the full detailed help message.
      --show_hidden                 [boolean]          Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 20 param(s), use the `--show_hidden` parameter to show them !!
    ----------------------------------------------------

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

Per exemple, per ometre el pas de retallada:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-notrim --skip_trim
```

??? success "Sortida de la comanda"

    ```console
    executor >  local (4)
    [3f/a82c91] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [7d/c5e014] NFCORE_DEMO:DEMO:MULTIQC             | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El procés `SEQTK_TRIM` ja no apareix a la sortida.

!!! info "Info"

    Tot i que tècnicament és possible establir paràmetres del pipeline en un fitxer de configuració personalitzat passat amb `-c`, és possible que no sobreescrigui els valors per defecte ja establerts al `nextflow.config` propi del pipeline, depenent de les regles de precedència de configuració de Nextflow.
    Utilitzar `--nom_parametre` a la línia de comandes o `-params-file` és més fiable, ja que aquests sempre tenen prioritat.

    **Com a regla general:** si apareix a la sortida de `--help`, establiu-lo mitjançant la línia de comandes o un fitxer de paràmetres en lloc d'un fitxer de configuració.

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
Això detecta errors tipogràfics com `--outDir` en lloc de `--outdir` abans que malgasteu temps de còmput preguntant-vos per què la sortida ha anat al lloc equivocat.

##### 3.1.3.2. Valors de paràmetres no vàlids

La validació també comprova els **valors** dels paràmetres.
El paràmetre `--skip_trim` és un flag booleà, de manera que passar un valor de tipus string fa que el pipeline falli immediatament:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]
```

El pipeline s'atura abans que s'executi cap procés, estalviant-vos una execució fallida o incorrecta.
Els paràmetres booleans s'han de passar com a flags (`--skip_trim`) sense cap valor, o establir-se a `true`/`false` en un fitxer de paràmetres.

#### 3.1.4. Validació d'entrada

La mateixa lògica de validació també es pot utilitzar per comprovar la validesa dels fitxers d'entrada.
Per exemple, si un pipeline espera un samplesheet com a entrada de dades principal (que és el cas de molts, si no la majoria, dels pipelines nf-core), el desenvolupador pot proporcionar un esquema d'entrada (diferent de l'esquema de paràmetres) que descrigui com s'ha d'estructurar el fitxer d'entrada.

Llavors, en temps d'execució, Nextflow pot comprovar que el fitxer d'entrada proporcionat és vàlid.

També ho tractem amb més detall a [Part 5: Input Validation](05_input_validation.md), però ja podeu veure-ho en acció donant al pipeline de demostració un samplesheet d'entrada no vàlid.

El pipeline `nf-core/demo` espera un fitxer CSV amb les columnes `sample`, `fastq_1` i `fastq_2`.
Això es defineix en un fitxer d'esquema (`assets/schema_input.json`) que especifica l'estructura esperada, els tipus de columnes i les restriccions.

??? abstract "assets/schema_input.json"

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

##### 3.1.4.1. Crear un samplesheet no vàlid

Creeu un samplesheet amb una columna que falta i un camí de fitxer inexistent:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

A aquest samplesheet li falta la columna obligatòria `fastq_1` i té un camí de fitxer inexistent a `fastq_2`.
Tots dos problemes produiran errors de validació al pas següent.

##### 3.1.4.2. Executar el pipeline de demostració amb el samplesheet no vàlid

Executeu el pipeline de demostració utilitzant `malformed_samplesheet.csv` com a entrada.

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

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
Llisteu els fitxers de configuració per veure què hi ha disponible:

```bash
ls $NXF_HOME/assets/nf-core/demo/conf/
```

```console
base.config  igenomes.config  igenomes_ignored.config  modules.config  test.config  test_full.config
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

Fem alguns exercicis per practicar-ho.

#### 3.2.1. Canviar l'assignació de recursos per a un procés

El pipeline de demostració assigna recursos utilitzant etiquetes definides a `base.config`.
Per exemple, `FASTQC` utilitza l'etiqueta `process_medium`, que assigna 6 CPUs i 36 GB de memòria.

El perfil de prova limita els recursos mitjançant `resourceLimits`, però també podeu sobreescriure els recursos per a processos específics.

Creeu un fitxer anomenat `custom.config`:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
}
```

Executeu el pipeline amb la vostra configuració personalitzada:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-custom -c custom.config
```

??? success "Sortida de la comanda"

    ```console
    executor >  local (7)
    [2a/f17b3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [9c/e4d028] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [5b/a93c71] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El flag `-c` afegeix la vostra configuració a sobre de la configuració integrada del pipeline.

#### 3.2.2. Establir valors d'arguments d'eines amb `ext.args`

Moltes eines de línia de comandes tenen arguments que no són obligatoris i, per tant, no es configuren com a paràmetres del pipeline tret que s'utilitzin molt freqüentment.
Per a aquests arguments d'eines, els mòduls nf-core utilitzen una convenció de Nextflow anomenada `ext.args` per passar arguments a l'eina subjacent mitjançant un fitxer de configuració.

Per exemple, afegim un argument de retallada al mòdul `SEQTK_TRIM` utilitzant `ext.args`.

##### 3.2.2.1. Actualitzar la configuració personalitzada

Actualitzeu el vostre `custom.config`:

```groovy title="custom.config" linenums="1" hl_lines="6 7 8"
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

Això indica a `seqtk trimfq` que retalli 5 bases del principi de cada lectura a més de la retallada per qualitat.

##### 3.2.2.2. Executar el pipeline

Executeu el pipeline de nou amb aquesta configuració per veure l'efecte:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-extargs -c custom.config
```

??? success "Sortida de la comanda"

    ```console
    executor >  local (7)
    [1e/b7a392] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [ab/cd1234] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [4f/c8d105] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Per verificar que l'argument s'ha aplicat, trobeu el hash del directori de treball de `SEQTK_TRIM` a la sortida de l'execució (per exemple, `work/ab/cd1234...`) i comproveu el fitxer `.command.sh` que hi ha dins:

```bash
cat work/ab/cd1234/.command.sh
```

??? success "Sortida de la comanda"

    ```console
    #!/usr/bin/env bash
    ...
    seqtk trimfq -b 5 SAMPLE3_SE.fastq.gz | gzip -c > SAMPLE3_SE.trimmed.fastq.gz
    ```

Hauríeu de veure `-b 5` a la comanda `seqtk trimfq`, confirmant que la sobreescriptura de `ext.args` ha tingut efecte.

##### 3.2.2.3. Sobreescriure valors per defecte

Alguns mòduls ja tenen `ext.args` establert per defecte.
Per exemple, el mòdul `FASTQC` està configurat amb `ext.args = '--quiet'` per defecte (definit a `conf/modules.config`).

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}"
        ]
    }
```

Si proporcioneu un valor per a `ext.args` mitjançant un fitxer de configuració personalitzat, aquest valor reemplaçarà completament el valor per defecte establert per a aquell procés.

Així, per exemple, si el valor per defecte era `'--quiet'` i establiu `ext.args = '--kmers 8'`, el flag `--quiet` ja no s'aplicarà.
Per mantenir tots dos, establiu `ext.args = '--quiet --kmers 8'`.

Això significa que sou responsables de comprovar quina és la configuració per defecte de les eines a les quals voleu proporcionar valors d'arguments amb `ext.args`.

### Conclusió

Sabeu com obtenir ajuda d'un pipeline nf-core, establir paràmetres i entendre com es validen, i personalitzar la configuració mitjançant fitxers de configuració.

### Què segueix?

Feu una pausa! Quan estigueu preparats, passeu a la Part 2, on creareu el vostre propi pipeline compatible amb nf-core des de zero.
