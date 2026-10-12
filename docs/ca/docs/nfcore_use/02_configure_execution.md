# Part 2: Configurar l'execució del pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la [Part 1](./01_run_demo.md), heu trobat i executat el pipeline nf-core/demo utilitzant el seu perfil de prova.
Ara veurem com configurar l'execució del pipeline: establir paràmetres, entendre la validació i personalitzar l'assignació de recursos i els arguments de les eines.

Tal com s'explica a [Hello Config](../hello_nextflow/06_hello_config.md), volem poder canviar les dades sobre les quals s'executarà el nostre pipeline i com s'executarà sense modificar el codi del pipeline en si.
Per a això, Nextflow admet múltiples maneres de controlar la configuració del pipeline, la qual cosa pot resultar una mica aclaparadora.

El projecte nf-core especifica convencions per organitzar els elements de configuració, distingint dos tipus de configuració al nivell superior: **paràmetres del pipeline** i **configuració** en sentit estricte.

- **Paràmetres del pipeline** (establerts mitjançant el sistema `params`) inclouen típicament coses com fitxers d'entrada, indicadors de comportament de les eines i paràmetres d'anàlisi.
- **Configuració** en sentit estricte fa referència a la logística de com s'executa el pipeline, és a dir, l'executor, l'assignació de recursos de càlcul, etc.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Comencem per abordar els paràmetres del pipeline i, a continuació, veurem la configuració en sentit estricte.

---

## 1. Paràmetres del pipeline

Per a tots els pipelines nf-core, podeu obtenir una llista completa dels paràmetres del pipeline directament des de la línia de comandes utilitzant l'indicador `--help`, que és en si mateix un paràmetre del pipeline.

### 1.1. Obtenir la llista de paràmetres amb `--help`

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

Com podeu veure, la sortida agrupa els paràmetres en categories (opcions d'entrada/sortida, opcions del genoma de referència, etc.) amb tipus i descripcions per a cadascun.

Aquesta categorització està determinada per un fitxer d'esquema, que es tracta més endavant.
En pipelines Nextflow simples, `--help` només funciona si el desenvolupador l'ha implementat manualment.

!!! tip "Consell"

    Utilitzeu `--help --show_hidden` per veure paràmetres addicionals que estan ocults per defecte, com ara `--publish_dir_mode` o `--monochrome_logs`.

### 1.2. Establir valors de paràmetres

Tal com es tracta a [Hello Config](../hello_nextflow/06_hello_config.md), podeu establir valors de paràmetres a la línia de comandes amb `--nom_param` o recollir un conjunt de paràmetres en un fitxer YAML i passar-lo amb `-params-file`.
Tots dos enfocaments funcionen de la mateixa manera amb els pipelines nf-core.

Per exemple, per ometre el pas de retallada, volem establir el paràmetre booleà `skip_trim` a `true`.
Un fitxer de paràmetres anomenat `my_params.yml` es proporciona al vostre directori de treball amb aquest valor ja establert:

```yaml title="my_params.yml"
skip_trim: true
```

Passeu-lo amb `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
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

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

El procés `SEQTK_TRIM` ja no apareix a la sortida.

!!! warning "Limitacions importants sobre les entrades de paràmetres"

    **Establir paràmetres booleans a la línia de comandes**

    A partir de la versió 26.04 de Nextflow, tots els valors subministrats a la línia de comandes es tipifiquen com a strings.
    Per a un paràmetre booleà com `skip_trim`, passar-lo com a indicador simple (`--skip_trim`) o com a `--skip_trim true` s'avalua com el **string** `"true"`, la qual cosa falla en la validació de l'esquema:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Per establir un paràmetre booleà a un valor genuí `true`/`false`, utilitzeu un `-params-file` tal com es mostra més amunt, o establiu-lo en un fitxer de configuració.
    Els paràmetres de tipus string, integer i ruta de fitxer no es veuen afectats i es poden continuar establint directament a la línia de comandes.
    Aquest curs utilitza aquest patró al llarg de tot el material per als paràmetres booleans.

    **Ús de fitxers de configuració personalitzats**

    Tot i que és tècnicament possible establir paràmetres del pipeline en un fitxer de configuració personalitzat passat amb `-c`, és possible que no sobreescrigui els valors per defecte ja establerts al `nextflow.config` propi del pipeline, depenent de les regles de precedència de configuració de Nextflow.
    Utilitzar `--nom_param` a la línia de comandes o `-params-file` és més fiable, ja que aquests sempre tenen precedència.

    Com a regla general: si apareix a la sortida de `--help`, establiu-lo mitjançant la línia de comandes o un fitxer de paràmetres en lloc d'un fitxer de configuració.

### 1.3. Validació de paràmetres

Curiositat: la comanda `--help` funciona per a tots els pipelines nf-core perquè el projecte nf-core requereix que els desenvolupadors defineixin formalment tots els paràmetres del pipeline en un fitxer d'esquema JSON (`nextflow_schema.json`).
Aquest esquema registra el tipus, la descripció, el valor per defecte i l'agrupació de cada paràmetre.

A més de potenciar la sortida de `--help`, el fitxer d'esquema també permet la validació automatitzada en el moment del llançament.
Això significa que Nextflow pot comprovar que cada paràmetre que passeu existeix i se li ha donat un valor adequat (del tipus adequat, dins del rang de valors permès, etc.).

Ho tractem amb més detall a la [secció de validació d'entrades](../nfcore_build/04_input_validation.md), però ja podeu veure-ho en acció donant al pipeline de demostració alguna entrada de paràmetres no vàlida.

#### 1.3.1. Paràmetres no reconeguts

Proveu de passar un paràmetre que no existeix:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

La sortida de la consola inclou un avís:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

El pipeline continua executant-se, però l'avís us alerta immediatament que `--foobar` no és un paràmetre reconegut.
Això té com a objectiu cridar la vostra atenció sobre errors tipogràfics no crítics, com ara l'ús de `--outDir` en lloc de `--outdir`, la qual cosa us pot ajudar a evitar malgastar temps i recursos de càlcul.

#### 1.3.2. Valors de paràmetres no vàlids

La validació també comprova els **valors** dels paràmetres.
El paràmetre `--skip_trim` és un indicador booleà, de manera que passar un valor de tipus string fa que el pipeline falli immediatament:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

El pipeline s'atura abans que s'executi cap procés, estalviant-vos una execució fallida o incorrecta.
Tal com s'indica a [1.2](#12-set-parameter-values), els paràmetres booleans s'han d'establir a un valor genuí `true`/`false` en un fitxer de paràmetres en lloc de passar-los a la línia de comandes, ja que els valors de la línia de comandes es tipifiquen com a strings.

### 1.4. Validació d'entrades

La mateixa lògica de validació també es pot utilitzar per comprovar la validesa dels fitxers d'entrada.
Per exemple, si un pipeline espera un full de mostres com a entrada de dades principal (que és el cas de molts, si no de la majoria, dels pipelines nf-core), el desenvolupador pot proporcionar un esquema d'entrada (diferent de l'esquema de paràmetres) que descrigui com ha d'estar estructurat el fitxer d'entrada.

Aleshores, en temps d'execució, Nextflow pot comprovar que el fitxer d'entrada proporcionat és vàlid.

També ho tractem amb més detall a la [secció de validació d'entrades](../nfcore_build/04_input_validation.md), però ja podeu veure-ho en acció donant al pipeline de demostració un full de mostres d'entrada no vàlid.

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

L'esquema especifica que `sample` i `fastq_1` són obligatoris, mentre que `fastq_2` és opcional (admetent dades tant de parells com de lectura única).
Les rutes de fitxer es validen per a l'existència i el patró d'extensió.

Per demostrar-ho, proporcionem un full de mostres mal format anomenat `malformed_samplesheet.csv` al vostre directori de treball:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

A aquest full de mostres li falta la columna obligatòria `fastq_1` i té una ruta de fitxer inexistent a `fastq_2`.

Executeu el pipeline de demostració utilitzant `malformed_samplesheet.csv` com a entrada:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Com podeu veure, el pipeline falla immediatament i informa de **tots** els errors de validació alhora.
nf-schema no s'atura al primer error — recull tots els problemes i els llista junts, de manera que podeu corregir-ho tot d'una vegada en lloc de descobrir els problemes un per un.

Cada error identifica l'entrada i el camp exactes que han causat el problema, de manera que podeu corregir el vostre full de mostres i tornar a llançar el pipeline amb la confiança que no fallarà en algun punt posterior quan Nextflow intenti accedir a la ruta del fitxer.

Per als desenvolupadors, tot això es tracta amb més detall a la [Part 4 de Build with nf-core](../nfcore_build/04_input_validation.md).

### Conclusió

Sabeu com obtenir una llista completa dels paràmetres d'un pipeline amb `--help`, establir-los mitjançant la línia de comandes o un fitxer de paràmetres, i com Nextflow valida tant els valors dels paràmetres com els fitxers d'entrada respecte als esquemes del pipeline.

### Què segueix?

Apreneu sobre l'altre tipus de configuració: com s'executa el pipeline, incloent-hi l'assignació de recursos i els arguments de les eines.

---

## 2. Configuració

La configuració en sentit estricte controla **com** s'executa el pipeline: l'assignació de recursos, els arguments específics de les eines, on s'executen les tasques i quin sistema d'empaquetament de programari s'utilitza.

Els pipelines nf-core inclouen la configuració per defecte a `nextflow.config` i al directori `conf/`.
Abans de sobreescriure res, és útil saber on viuen els valors per defecte.

### 2.1. Explorar els fitxers de configuració

Ja heu vist a la [Part 1](./01_run_demo.md) que el codi font del pipeline es troba a `$NXF_HOME/assets`.
Utilitzant l'enllaç simbòlic `pipelines` que heu creat a la [Part 1](./01_run_demo.md), llisteu els fitxers de configuració per veure què hi ha disponible:

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
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

Els fitxers de configuració més importants són:

- **`conf/base.config`**: Defineix etiquetes de recursos (`process_low`, `process_medium`, `process_high`) que assignen CPUs, memòria i temps als processos. Quan veieu un procés que utilitza més recursos dels esperats, aquí és d'on provenen aquests valors per defecte.
- **`conf/modules.config`**: Estableix els arguments de les eines per procés (`ext.args`) i la configuració de publicació de sortides (`publishDir`). Obriu aquest fitxer per veure quins arguments rep cada eina per defecte.
- **`conf/test.config`**: El perfil de prova que heu utilitzat a la [Part 1](./01_run_demo.md), que limita els recursos mitjançant `resourceLimits` i estableix un full de mostres de prova. S'activa amb `-profile test`.
  També hi ha un `conf/test_full.config` per executar amb un conjunt de dades de prova de mida completa, útil per a la comparació de rendiment.

El `nextflow.config` central carrega tots els anteriors i estableix els valors per defecte adequats per a tot.

Si voleu modificar qualsevol dels paràmetres especificats en aquests fitxers, no els modifiqueu directament.
En lloc d'això, creeu el vostre propi fitxer de configuració i passeu-lo amb `-c`.
Els valors que especifiqueu sobreescriuran els valors per defecte establerts en aquells altres fitxers.

Provem-ho a la pràctica.

### 2.2. Personalitzar els recursos dels processos i els arguments de les eines

Els mòduls nf-core admeten dos tipus comuns de sobreescriptura de configuració: **assignació de recursos** (CPUs, memòria, temps) i **arguments de les eines** mitjançant `ext.args`.

Moltes eines de línia de comandes tenen arguments que no s'utilitzen prou habitualment com per exposar-los com a paràmetres del pipeline.
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
L'indicador `-b 5` indica a `seqtk trimfq` que retalli 5 bases del principi de cada lectura a més de la retallada per qualitat.

Executeu el pipeline amb aquesta configuració:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
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

L'indicador `-c` afegeix la vostra configuració a sobre de la configuració integrada del pipeline.

Per verificar que la sobreescriptura de `ext.args` ha tingut efecte, trobeu el hash del directori de treball de `SEQTK_TRIM` a la sortida de l'execució (p. ex. `work/17/428668...`) i comproveu el fitxer `.command.sh` que hi ha dins:

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

Una cosa important que cal saber sobre `ext.args`: si un mòdul ja té un valor per defecte establert, el vostre valor el **substituirà completament** en lloc d'afegir-s'hi.
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

Si establiu `ext.args = '--kmers 8'` per a `FASTQC`, l'indicador `--quiet` ja no s'aplicarà.
Per mantenir tots dos, establiu `ext.args = '--quiet --kmers 8'`.

Sempre hauríeu de comprovar la configuració per defecte d'un mòdul abans de sobreescriure `ext.args`.

### Conclusió

Sabeu on viuen els valors per defecte de configuració dels pipelines nf-core i com sobreescriure les assignacions de recursos i els arguments de les eines amb un fitxer de configuració personalitzat.

### Què segueix?

Aneu a la [Part 3](./03_run_production_pipeline.md), on aplicareu el que heu après a un pipeline de producció real.

---

## Resum

En aquesta part heu après a:

- Obtenir ajuda, establir paràmetres i entendre la validació de paràmetres i d'entrades
- Personalitzar l'assignació de recursos i els arguments de les eines mitjançant fitxers de configuració
