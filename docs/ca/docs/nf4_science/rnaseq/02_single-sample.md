# Part 2: Implementació d'una sola mostra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta part del curs, escriurem el workflow més senzill possible que encapsula totes les comandes que vam executar a la Part 1 per automatitzar-ne l'execució, i processarem només una mostra cada vegada.

!!! warning "Prerequisit"

    Heu de completar la [Part 1: Visió general del mètode](./01_method.md) abans de començar aquesta lliçó.
    Concretament, treballar la secció 1.2.3 crea el fitxer d'índex del genoma (`data/genome_index.tar.gz`) necessari per al pas d'alineament d'aquesta lliçó.

## Assignació

En aquesta part del curs, desenvoluparem un workflow que fa el següent:

1. Executar control de qualitat (FastQC) sobre les lectures d'entrada
2. Retallar adaptadors i executar QC post-retallat (Trim Galore)
3. Alinear les lectures retallades a un genoma de referència (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

Això automatitza els passos de la primera secció de la [Part 1: Visió general del mètode](./01_method.md#1-single-sample-processing), on vau executar aquestes comandes manualment dins dels seus contenidors.

Com a punt de partida, us proporcionem un fitxer de workflow, `rnaseq.nf`, que esbossa les parts principals del workflow, així com quatre fitxers de mòdul al directori `modules/` (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf` i `multiqc.nf`) que esbossen l'estructura de cada procés.

??? full-code "Fitxers d'estructura"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Module INCLUDE statements

    /*
     * Pipeline parameters
     */

    // Primary input

    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }

    output {
        // Configure publish targets
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

Aquests fitxers no són funcionals; el seu propòsit és només servir com a estructures perquè les ompliu amb les parts interessants del codi.

## Pla de la lliçó

Per tal de fer el procés de desenvolupament més educatiu, ho hem dividit en tres etapes:

1. **Escriure un workflow d'una sola etapa que executi el pas inicial de QC.**
   Això cobreix la configuració d'un parametre CLI, la creació d'un canal d'entrada, l'escriptura d'un mòdul de procés i la configuració de la publicació de sortides.
2. **Afegir el retallat d'adaptadors i QC post-retallat.**
   Això introdueix l'encadenament de processos connectant la sortida d'un procés amb l'entrada d'un altre.
3. **Afegir l'alineament al genoma de referència.**
   Això cobreix la gestió d'entrades de referència addicionals i el treball amb arxius comprimits.

Cada pas se centra en un aspecte específic del desenvolupament de workflows.

!!! tip "Consell"

     Assegureu-vos que esteu al directori de treball correcte:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Escriure un workflow d'una sola etapa que executi el QC inicial

Aquest primer pas se centra en els conceptes bàsics: carregar un fitxer FASTQ i executar-hi control de qualitat.

Recordeu la comanda `fastqc` de la [Part 1](01_method.md):

```bash
fastqc <reads>
```

La comanda pren un fitxer FASTQ com a entrada i produeix un informe de control de qualitat com a arxiu `.zip` i un resum `.html`.
L'URI del contenidor era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Prendrem aquesta informació i l'encapsularem en Nextflow en tres etapes:

1. Configurar l'entrada
2. Escriure el procés de QC i cridar-lo al workflow
3. Configurar la gestió de sortides

### 1.1. Configurar l'entrada

Hem de declarar un parametre d'entrada, crear un perfil de prova per proporcionar un valor per defecte convenient, i crear un canal d'entrada.

#### 1.1.1. Afegir una declaració de parametre d'entrada

A `rnaseq.nf`, sota la secció `Pipeline parameters`, declareu un parametre anomenat `reads` amb el tipus `Path`.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        input: Path
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Això configura el parametre CLI, però no volem escriure el camí del fitxer cada vegada que executem el workflow durant el desenvolupament.
Hi ha múltiples opcions per proporcionar un valor per defecte; aquí utilitzem un perfil de prova.

#### 1.1.2. Crear un perfil de prova amb un valor per defecte a `nextflow.config`

Un perfil de prova proporciona valors per defecte convenients per provar un workflow sense especificar entrades a la línia de comandes.
Aquesta és una convenció comuna a l'ecosistema Nextflow (consulteu [Hello Config](../../hello_nextflow/06_hello_config.md) per a més detalls).

Afegiu un bloc `profiles` a `nextflow.config` amb un perfil `test` que estableixi el parametre `reads` a un dels fitxers FASTQ de prova.

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aquí, estem utilitzant `#!groovy ${projectDir}`, una variable integrada de Nextflow que apunta al directori on es troba l'script del workflow.
Això facilita la referència a fitxers de dades i altres recursos sense codificar camins absoluts.

El parametre ara té un valor per defecte convenient. A continuació, hem de crear un canal a partir d'ell.

#### 1.1.3. Configurar el canal d'entrada

Al bloc workflow, creeu un canal d'entrada a partir del valor del parametre utilitzant la factoria de canals `.fromPath` (com s'utilitza a [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Després"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Create input channel

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

A continuació, haurem de crear el procés per executar QC sobre aquesta entrada.

### 1.2. Escriure el procés de QC i cridar-lo al workflow

Hem d'omplir la definició del procés al fitxer de mòdul, importar-lo al workflow utilitzant una declaració include, i cridar-lo sobre l'entrada.

#### 1.2.1. Omplir el mòdul per al procés de QC

Obriu `modules/fastqc.nf` i examineu l'esquema de la definició del procés.
Hauríeu de reconèixer els elements estructurals principals; si no, considereu llegir [Hello Nextflow](../../hello_nextflow/01_hello_world.md) per refrescar-vos.

Aneu endavant i ompliu la definició del procés per vosaltres mateixos utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball amb la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Després"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Executa FastQC sobre les lectures d'entrada
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

L'accessor `simpleName` elimina totes les extensions del nom del fitxer, així que `ENCSR000COQ1_1.fastq.gz` es converteix en `ENCSR000COQ1_1`.
Utilitzem la sintaxi `emit:` per assignar noms a cada canal de sortida, cosa que serà útil per connectar sortides al bloc publish.

Un cop hàgiu completat això, el procés està complet.
Per utilitzar-lo al workflow, haureu d'importar el mòdul i afegir una crida al procés.

#### 1.2.2. Incloure el mòdul

A `rnaseq.nf`, afegiu una declaració `include` per fer el procés disponible al workflow:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    ```

El procés ara està disponible a l'àmbit del workflow.

#### 1.2.3. Cridar el procés de QC sobre l'entrada

Afegiu una crida a `FASTQC` al bloc workflow, passant el canal d'entrada com a argument.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Control de qualitat inicial
        FASTQC(read_ch)

        publish:
        // Declare outputs to publish
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Call processes

        publish:
        // Declare outputs to publish
    }
    ```

El workflow ara carrega l'entrada i executa el procés de QC sobre ella.
A continuació, hem de configurar com es publica la sortida.

### 1.3. Configurar la gestió de sortides

Hem de declarar quines sortides de procés publicar i especificar on han d'anar.

#### 1.3.1. Declarar sortides a la secció `publish:`

La secció `publish:` dins del bloc workflow declara quines sortides de procés s'han de publicar.
Assigneu les sortides de `FASTQC` a objectius amb nom.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Declare outputs to publish
    }
    ```

A continuació, haurem de dir a Nextflow on posar les sortides publicades.

#### 1.3.2. Configurar els objectius de sortida al bloc `output {}`

El bloc `output {}` es troba fora del workflow i especifica on es publica cada objectiu amb nom.
Configureu ambdós objectius per publicar-se en un subdirectori `fastqc/`.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configure publish targets
    }
    ```

!!! note "Nota"

    Per defecte, Nextflow publica fitxers de sortida com a enllaços simbòlics, cosa que evita duplicacions innecessàries.
    Tot i que els fitxers de dades que estem utilitzant aquí són molt petits, en genòmica poden arribar a ser molt grans.
    Els enllaços simbòlics es trenquen quan netegeu el vostre directori `work`, així que per a workflows de producció potser voldreu sobreescriure el mode de publicació per defecte a `'copy'`.

### 1.4. Executar el workflow

En aquest punt, tenim un workflow de QC d'un sol pas que hauria de ser completament funcional.

Executem amb `-profile test` per utilitzar el valor per defecte configurat al perfil de prova, evitant la necessitat d'escriure el camí a la línia de comandes.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

Això hauria d'executar-se molt ràpidament si heu treballat la Part 1 i ja heu descarregat el contenidor.
Si l'heu omès, Nextflow descarregarà el contenidor per vosaltres; no heu de fer res perquè passi, però potser haureu d'esperar fins a un minut.

Podeu comprovar les sortides al directori de resultats.

```bash
ls results/fastqc
```

```console title="Output"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

Els informes de QC per a la mostra ara estan publicats al subdirectori `fastqc/`.

### Conclusió

Sabeu com crear un mòdul que conté un procés, importar-lo a un workflow, cridar-lo amb un canal d'entrada i publicar els resultats utilitzant el bloc output a nivell de workflow.

### Què segueix?

Afegir el retallat d'adaptadors amb QC post-retallat com a segon pas al workflow.

---

## 2. Afegir retallat d'adaptadors i QC post-retallat

Ara que tenim el QC inicial en funcionament, podem afegir el pas de retallat d'adaptadors amb el seu QC post-retallat integrat.

Recordeu la comanda `trim_galore` de la [Part 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

La comanda retalla adaptadors d'un fitxer FASTQ i executa FastQC sobre la sortida retallada.
Produeix lectures retallades, un informe de retallat i informes FastQC per a les lectures retallades.
L'URI del contenidor era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Només hem d'escriure la definició del procés, importar-lo, cridar-lo al workflow i actualitzar la gestió de sortides.

### 2.1. Escriure el procés de retallat i cridar-lo al workflow

Com abans, hem d'omplir la definició del procés, importar el mòdul i afegir la crida al procés.

#### 2.1.1. Omplir el mòdul per al procés de retallat

Obriu `modules/trim_galore.nf` i examineu l'esquema de la definició del procés.

Aneu endavant i ompliu la definició del procés per vosaltres mateixos utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball amb la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Després"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Retalla adaptadors i executa QC post-retallat
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    }
    ```

Aquest procés té tres sortides amb nom: les lectures retallades que alimenten el pas d'alineament, l'informe de retallat i els informes FastQC post-retallat.
La bandera `--fastqc` indica a Trim Galore que executi automàticament FastQC sobre la sortida retallada.

#### 2.1.2. Incloure el mòdul

Actualitzeu `rnaseq.nf` per importar el nou mòdul:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    ```

A continuació, afegirem la crida al procés al workflow.

#### 2.1.3. Cridar el procés de retallat sobre l'entrada

Afegiu la crida al procés al bloc workflow:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Control de qualitat inicial
        FASTQC(read_ch)

        // Retallat d'adaptadors i QC post-retallat
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Control de qualitat inicial
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

El procés de retallat ara està connectat al workflow.

### 2.2. Actualitzar la gestió de sortides

Hem d'afegir les sortides de retallat a la declaració publish i configurar on van.

#### 2.2.1. Afegir objectius publish per a les sortides de retallat

Afegiu les sortides de retallat a la secció `publish:`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

A continuació, haurem de dir a Nextflow on posar aquestes sortides.

#### 2.2.2. Configurar els nous objectius de sortida

Afegiu entrades per als objectius de retallat al bloc `output {}`, publicant-los en un subdirectori `trimming/`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

La configuració de sortida està completa.

### 2.3. Executar el workflow

El workflow ara inclou tant el QC inicial com el retallat d'adaptadors.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

Això també hauria d'executar-se molt ràpidament, ja que estem executant sobre un fitxer d'entrada tan petit.

Podeu trobar les sortides de retallat al directori de resultats.

```bash
ls results/trimming
```

```console title="Output"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

Les sortides de retallat i els informes de QC post-retallat ara estan al subdirectori `trimming/`.

### Conclusió

Sabeu com afegir un segon pas de processament que s'executa independentment sobre la mateixa entrada, produint múltiples sortides amb nom.

### Què segueix?

Afegir el pas d'alineament que s'encadena a partir de la sortida de lectures retallades.

---

## 3. Afegir alineament al genoma de referència

Finalment podem afegir el pas d'alineament del genoma utilitzant HISAT2.

Recordeu la comanda d'alineament de la [Part 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

La comanda alinea lectures a un genoma de referència i converteix la sortida a format BAM.
Requereix un arxiu d'índex del genoma preconstruït i produeix un fitxer BAM i un registre de resum d'alineament.
L'URI del contenidor era `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

Aquest procés requereix una entrada addicional (l'arxiu d'índex del genoma), així que primer hem de configurar-ho, després escriure i connectar el procés.

### 3.1. Configurar les entrades

Hem de declarar un parametre per a l'arxiu d'índex del genoma.

#### 3.1.1. Afegir un parametre per a l'índex del genoma

Afegiu una declaració de parametre per a l'arxiu d'índex del genoma a `rnaseq.nf`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Primary input
        input: Path

        // Arxiu del genoma de referència
        hisat2_index_zip: Path
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Primary input
        input: Path
    }
    ```

#### 3.1.2. Afegir el valor per defecte de l'índex del genoma al perfil de prova

Tal com vam fer per a `reads` a la secció 1.1.2, afegiu un valor per defecte per a l'índex del genoma al perfil de prova a `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

El parametre està llest; ara podem crear el procés d'alineament.

### 3.2. Escriure el procés d'alineament i cridar-lo al workflow

Com abans, hem d'omplir la definició del procés, importar el mòdul i afegir la crida al procés.

#### 3.2.1. Omplir el mòdul per al procés d'alineament

Obriu `modules/hisat2_align.nf` i examineu l'esquema de la definició del procés.

Aneu endavant i ompliu la definició del procés per vosaltres mateixos utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball amb la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Després"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Alinea lectures a un genoma de referència
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    }
    ```

Aquest procés pren dues entrades: les lectures i l'arxiu d'índex del genoma.
El bloc script primer extreu l'índex de l'arxiu, després executa l'alineament HISAT2 canalitzat a `samtools view` per convertir la sortida a format BAM.
L'accessor `simpleName` sobre `index_zip` extreu el nom base de l'arxiu (`genome_index`) per utilitzar-lo com a prefix de l'índex.

#### 3.2.2. Incloure el mòdul

Actualitzeu `rnaseq.nf` per importar el nou mòdul:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="3"
    // Module INCLUDE statements
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

A continuació, afegirem la crida al procés al workflow.

#### 3.2.3. Cridar el procés d'alineament

Les lectures retallades estan al canal `TRIM_GALORE.out.trimmed_reads` produït pel pas anterior.
Utilitzem `#!groovy file(params.hisat2_index_zip)` per proporcionar l'arxiu d'índex del genoma.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Control de qualitat inicial
        FASTQC(read_ch)

        // Retallat d'adaptadors i QC post-retallat
        TRIM_GALORE(read_ch)

        // Alineament a un genoma de referència
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Crea un canal d'entrada a partir d'un camí de fitxer
        read_ch = channel.fromPath(params.input)

        // Control de qualitat inicial
        FASTQC(read_ch)

        // Retallat d'adaptadors i QC post-retallat
        TRIM_GALORE(read_ch)
    ```

El procés d'alineament ara està connectat al workflow.

### 3.3. Actualitzar la gestió de sortides

Hem d'afegir les sortides d'alineament a la declaració publish i configurar on van.

#### 3.3.1. Afegir objectius publish per a les sortides d'alineament

Afegiu les sortides d'alineament a la secció `publish:`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

A continuació, haurem de dir a Nextflow on posar aquestes sortides.

#### 3.3.2. Configurar els nous objectius de sortida

Afegiu entrades per als objectius d'alineament al bloc `output {}`, publicant-los en un subdirectori `align/`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

La configuració de sortida està completa.

### 3.4. Executar el workflow

El workflow ara inclou els tres passos de processament: QC, retallat i alineament.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

Podeu trobar les sortides d'alineament al directori de resultats.

```bash
ls results/align
```

```console title="Output"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Això completa el processament bàsic que hem d'aplicar a cada mostra.

_Afegirem l'agregació d'informes MultiQC a la Part 3, després d'haver fet que el workflow accepti múltiples mostres alhora._

---

### Conclusió

Sabeu com encapsular tots els passos bàsics per processar mostres RNAseq d'extrem únic individualment.

### Què segueix?

Feu una pausa! Ha estat molt.

Quan us sentiu refrescats, aneu a la [Part 3](./03_multi-sample.md), on aprendreu com modificar el workflow per processar múltiples mostres en paral·lel, agregar informes de QC a través de tots els passos per a totes les mostres, i habilitar l'execució del workflow sobre dades RNAseq d'extrem aparellat.
