# Part 3: Implementació multi-mostra amb dades paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Anteriorment, heu construït un pipeline de detecció de variants per mostra que processava les dades de cada mostra de manera independent.
En aquesta part del curs, portarem el nostre workflow senzill al següent nivell convertint-lo en una potent eina d'automatització per lots capaç de gestionar un nombre arbitrari de mostres.
I mentre hi som, també l'actualitzarem perquè esperi dades paired-end, que són més comunes en estudis recents.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat la [Part 1: Visió general del mètode](./01_method.md), la [Part 2: Implementació d'una sola mostra](./02_single-sample.md) i teniu un pipeline `rnaseq.nf` funcional amb els fitxers de mòdul completats.

    Si no heu completat la Part 2 o voleu començar de nou per a aquesta part, podeu utilitzar la solució de la Part 2 com a punt de partida.
    Executeu aquestes comandes des de dins del directori `nf4-science/rnaseq/`:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Això us proporciona un workflow complet de processament d'una sola mostra.
    Podeu comprovar que s'executa correctament:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Assignació

En aquesta part del curs, ampliarem el workflow per fer el següent:

1. Llegir informació de mostres des d'un full de càlcul CSV
2. Executar QC per mostra, retallat i alineament en totes les mostres en paral·lel
3. Agregar tots els informes de QC en un informe MultiQC complet

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Això automatitza els passos de la segona secció de la [Part 1: Visió general del mètode](./01_method.md#2-multi-sample-qc-aggregation), on vau executar aquestes comandes manualment en els seus contenidors.

## Pla de la lliçó

Ho hem dividit en tres etapes:

1. **Fer que el workflow accepti múltiples mostres d'entrada.**
   Això cobreix el canvi d'un sol camí de fitxer a un full de càlcul CSV, analitzar-lo amb `splitCsv()` i executar tots els processos existents en múltiples mostres.
2. **Afegir generació d'informes de QC complets.**
   Això introdueix l'operador `collect()` per agregar sortides entre mostres, i afegeix un procés MultiQC per produir un informe combinat.
3. **Canviar a dades RNAseq paired-end.**
   Això cobreix l'adaptació de processos per a entrades paired-end (utilitzant tuples), la creació de mòduls paired-end i la configuració d'un perfil de prova separat.

Això implementa el mètode descrit a la [Part 1: Visió general del mètode](./01_method.md) (segona secció que cobreix el cas d'ús multi-mostra) i es basa directament en el workflow produït per la Part 2.

!!! tip "Consell"

     Assegureu-vos que esteu al directori de treball correcte:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Fer que el workflow accepti múltiples mostres d'entrada

Per executar-se en múltiples mostres, hem de canviar com gestionem l'entrada: en lloc de proporcionar un sol camí de fitxer, llegirem la informació de les mostres des d'un fitxer CSV.

Proporcionem un fitxer CSV que conté IDs de mostra i camins de fitxers FASTQ al directori `data/`.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Aquest fitxer CSV inclou una línia de capçalera que anomena les columnes.

Tingueu en compte que aquestes són encara dades de lectura single-end.

!!! warning "Advertència"

    Els camins de fitxer al CSV són camins absoluts que han de coincidir amb el vostre entorn.
    Si no esteu executant això a l'entorn de formació que proporcionem, haureu d'actualitzar els camins perquè coincideixin amb el vostre sistema.

### 1.1. Canviar l'entrada primària a un CSV de camins de fitxer al perfil de prova

Primer, hem d'actualitzar el perfil de prova a `nextflow.config` per proporcionar el camí del fitxer CSV en lloc del camí FASTQ únic.

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
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
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

A continuació, haurem d'actualitzar la creació del canal per llegir des d'aquest CSV.

### 1.2. Actualitzar la factory del canal per analitzar l'entrada CSV

Hem de carregar els continguts del fitxer al canal en lloc de només el camí del fitxer.

Podem fer-ho utilitzant el mateix patró que vam utilitzar a la [Part 2 de Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): aplicant l'operador [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) per analitzar el fitxer, després una operació `map` per extreure el camí del fitxer FASTQ de cada fila.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Crea el canal d'entrada des dels continguts d'un fitxer CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Crea el canal d'entrada des d'un camí de fitxer
        read_ch = channel.fromPath(params.input)
    ```

Una cosa que és nova en comparació amb el que vau trobar al curs Hello Nextflow és que aquest CSV té una línia de capçalera, així que afegim `#!groovy header: true` a la crida `splitCsv()`.
Això ens permet referenciar columnes per nom a l'operació `map`: `#!groovy row.fastq_path` extreu el camí del fitxer de la columna `fastq_path` de cada fila.

La gestió d'entrada està actualitzada i el workflow està llest per provar.

### 1.3. Executar el workflow

El workflow ara llegeix informació de mostres des d'un fitxer CSV i processa totes les mostres en paral·lel.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Aquesta vegada cada pas s'executa 6 vegades, una per cada mostra al fitxer CSV.

Això és tot el que va fer falta per aconseguir que el workflow s'executés en múltiples fitxers.
Nextflow gestiona tot el paral·lelisme per nosaltres.

### Conclusió

Sabeu com canviar d'una entrada d'un sol fitxer a una entrada multi-mostra basada en CSV que Nextflow processa en paral·lel.

### Què segueix?

Afegir un pas d'agregació d'informes de QC que combini mètriques de totes les mostres.

---

## 2. Agregar mètriques de QC de preprocessament en un sol informe MultiQC

Tot això produeix molts informes de QC, i no volem haver de buscar entre informes individuals.
Aquest és el punt perfecte per posar un pas d'agregació d'informes MultiQC.

Recordeu la comanda `multiqc` de la [Part 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

La comanda escaneja el directori actual per trobar fitxers de sortida de QC reconeguts i els agrega en un sol informe HTML.
L'URI del contenidor era `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Hem de configurar un paràmetre addicional, preparar les entrades, escriure el procés, connectar-lo i actualitzar la gestió de sortides.

### 2.1. Configurar les entrades

El procés MultiQC necessita un paràmetre de nom d'informe i les sortides de QC recollides de tots els passos anteriors agrupades.

#### 2.1.1. Afegir un paràmetre `report_id`

Afegiu un paràmetre per anomenar l'informe de sortida.

=== "Després"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Entrada primària
        input: Path

        // Arxiu del genoma de referència
        hisat2_index_zip: Path

        // ID de l'informe
        report_id: String
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Entrada primària
        input: Path

        // Arxiu del genoma de referència
        hisat2_index_zip: Path
    }
    ```

Afegiu el valor per defecte de l'ID de l'informe al perfil de prova:

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

A continuació, haurem de preparar les entrades per al procés MultiQC.

#### 2.1.2. Recollir i combinar sortides de QC dels passos anteriors

Hem de donar al procés `MULTIQC` totes les sortides relacionades amb QC dels passos anteriors agrupades.

Per això, utilitzem l'operador `.mix()`, que agrega múltiples canals en un de sol.
Comencem des de `channel.empty()` i barregem tots els canals de sortida que volem combinar.
Això és més net que encadenar `.mix()` directament a un dels canals de sortida, perquè tracta totes les entrades simètricament.

Al nostre workflow, les sortides relacionades amb QC a agregar són:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Les barregem en un sol canal, després utilitzem `.collect()` per agregar els informes de totes les mostres en una sola llista.

Afegiu aquestes línies al cos del workflow després de la crida `HISAT2_ALIGN`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Alineament a un genoma de referència
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Generació d'informes de QC complets
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="38"
        // Alineament a un genoma de referència
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

L'ús de variables intermèdies fa que cada pas sigui clar: `multiqc_files_ch` conté tots els fitxers de QC individuals barrejats en un canal, i `multiqc_files_list` és el paquet recollit llest per passar a MultiQC.

### 2.2. Escriure el procés d'agregació de QC i cridar-lo al workflow

Com abans, hem d'omplir la definició del procés, importar el mòdul i afegir la crida al procés.

#### 2.2.1. Omplir el mòdul per al procés d'agregació de QC

Obriu `modules/multiqc.nf` i examineu l'esquema de la definició del procés.

Aneu endavant i ompliu la definició del procés per vosaltres mateixos utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball contra la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Agrega informes de QC amb MultiQC
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

=== "Després"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Agrega informes de QC amb MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Aquest procés utilitza `#!groovy path '*'` com a qualificador d'entrada per als fitxers de QC.
El comodí `'*'` indica a Nextflow que prepari tots els fitxers recollits al directori de treball sense requerir noms específics.
L'entrada `val output_name` és una cadena que controla el nom del fitxer de l'informe.

La comanda `multiqc .` escaneja el directori actual (on estan tots els fitxers de QC preparats) i genera l'informe.

Un cop hàgiu completat això, el procés està llest per utilitzar.

#### 2.2.2. Incloure el mòdul

Afegiu la declaració d'importació a `rnaseq.nf`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Declaracions INCLUDE de mòduls
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="3"
    // Declaracions INCLUDE de mòduls
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Ara afegiu la crida al procés al workflow.

#### 2.2.3. Afegir la crida al procés

Passeu els fitxers de QC recollits i l'ID de l'informe al procés `MULTIQC`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

El procés MultiQC ara està connectat al workflow.

### 2.3. Actualitzar la gestió de sortides

Hem d'afegir les sortides de MultiQC a la declaració de publicació i configurar on van.

#### 2.3.1. Afegir objectius de publicació per a les sortides de MultiQC

Afegiu les sortides de MultiQC a la secció `publish:`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="52"
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

A continuació, haurem de dir a Nextflow on posar aquestes sortides.

#### 2.3.2. Configurar els nous objectius de sortida

Afegiu entrades per als objectius de MultiQC al bloc `output {}`, publicant-los en un subdirectori `multiqc/`:

=== "Després"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Abans"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

La configuració de sortida està completa.

### 2.4. Executar el workflow

Utilitzem `-resume` perquè els passos de processament anteriors estiguin en memòria cau i només s'executi el nou pas MultiQC.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

S'ha afegit una sola crida a MULTIQC després de les crides de procés en memòria cau.

Podeu trobar les sortides de MultiQC al directori de resultats.

```bash
tree -L 2 results/multiqc
```

```console title="Output"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Aquest últim fitxer `all_single-end.html` és l'informe agregat complet, convenientment empaquetat en un fitxer HTML fàcil de navegar.

### Conclusió

Sabeu com recollir sortides de múltiples canals, agrupar-les amb `.mix()` i `.collect()`, i passar-les a un procés d'agregació.

### Què segueix?

Adaptar el workflow per gestionar dades RNAseq paired-end.

---

## 3. Habilitar el processament de dades RNAseq paired-end

Ara mateix el nostre workflow només pot gestionar dades RNAseq single-end.
És cada vegada més comú veure dades RNAseq paired-end, així que volem poder gestionar-ho.

Fer que el workflow sigui completament agnòstic del tipus de dades requeriria utilitzar característiques del llenguatge Nextflow lleugerament més avançades, així que no ho farem aquí, però podem fer una versió de processament paired-end per demostrar què cal adaptar.

### 3.1. Copiar el workflow i actualitzar les entrades

Comencem copiant el fitxer del workflow single-end i actualitzant-lo per a dades paired-end.

#### 3.1.1. Copiar el fitxer del workflow

Creeu una còpia del fitxer del workflow per utilitzar com a punt de partida per a la versió paired-end.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Ara actualitzeu els paràmetres i la gestió d'entrada al nou fitxer.

#### 3.1.2. Afegir un perfil de prova paired-end

Proporcionem un segon fitxer CSV que conté IDs de mostra i camins de fitxers FASTQ emparellats al directori `data/`.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Afegiu un perfil `test_pe` a `nextflow.config` que apunti a aquest fitxer i utilitzi un ID d'informe paired-end.

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

El perfil de prova per a dades paired-end està llest.

#### 3.1.3. Actualitzar la factory del canal

L'operador `.map()` necessita agafar ambdós camins de fitxer FASTQ i retornar-los com una llista.

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Crea el canal d'entrada des dels continguts d'un fitxer CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Crea el canal d'entrada des dels continguts d'un fitxer CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

La gestió d'entrada està configurada per a dades paired-end.

### 3.2. Adaptar el mòdul FASTQC per a dades paired-end

Copieu el mòdul per crear una versió paired-end:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

L'entrada del procés FASTQC no necessita canviar — quan Nextflow rep una llista de dos fitxers, prepara ambdós i `reads` s'expandeix a ambdós noms de fitxer.
L'únic canvi necessari és al bloc de sortida: com ara obtenim dos informes FastQC per mostra, canviem de patrons basats en `simpleName` a comodins.

=== "Després"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Abans"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Això generalitza el procés d'una manera que el fa capaç de gestionar dades single-end o paired-end.

Actualitzeu la importació a `rnaseq_pe.nf` per utilitzar la versió paired-end:

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

El mòdul FASTQC i la seva importació estan actualitzats per a dades paired-end.

### 3.3. Adaptar el mòdul TRIM_GALORE per a dades paired-end

Copieu el mòdul per crear una versió paired-end:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Aquest mòdul necessita canvis més substancials:

- L'entrada canvia d'un sol camí a una tupla de dos camins
- La comanda afegeix la bandera `--paired` i pren ambdós fitxers de lectura
- La sortida canvia per reflectir les convencions de nomenclatura paired-end de Trim Galore, produint informes FastQC separats per a cada fitxer de lectura

=== "Després"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Abans"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
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
    ```

Actualitzeu la importació a `rnaseq_pe.nf`:

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

El mòdul TRIM_GALORE i la seva importació estan actualitzats per a dades paired-end.

### 3.4. Adaptar el mòdul HISAT2_ALIGN per a dades paired-end

Copieu el mòdul per crear una versió paired-end:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Aquest mòdul necessita canvis similars:

- L'entrada canvia d'un sol camí a una tupla de dos camins
- La comanda HISAT2 canvia de `-U` (no emparellat) a arguments de lectura `-1` i `-2` (emparellats)
- Tots els usos de `reads.simpleName` canvien a `read1.simpleName` ja que ara referenciem un membre específic de la parella

=== "Després"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Abans"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
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
    ```

Actualitzeu la importació a `rnaseq_pe.nf`:

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

El mòdul HISAT2_ALIGN i la seva importació estan actualitzats per a dades paired-end.

### 3.5. Actualitzar l'agregació MultiQC per a sortides paired-end

El procés `TRIM_GALORE` paired-end ara produeix dos canals d'informes FastQC separats (`fastqc_reports_1` i `fastqc_reports_2`) en lloc d'un.
Actualitzeu el bloc `.mix()` a `rnaseq_pe.nf` per incloure ambdós:

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

L'agregació MultiQC ara inclou ambdós conjunts d'informes FastQC paired-end.

### 3.6. Actualitzar la gestió de sortides per a sortides paired-end

La secció `publish:` i el bloc `output {}` també necessiten reflectir els dos canals d'informes FastQC separats del procés `TRIM_GALORE` paired-end.

Actualitzeu la secció `publish:` a `rnaseq_pe.nf`:

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Actualitzeu les entrades corresponents al bloc `output {}`:

=== "Després"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Abans"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

El workflow paired-end ara està completament actualitzat i llest per executar.

### 3.7. Executar el workflow

No utilitzem `-resume` ja que això no es guardaria en memòria cau, i hi ha el doble de dades a processar que abans, però encara hauria de completar-se en menys d'un minut.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Ara tenim dues versions lleugerament divergents del nostre workflow, una per a dades de lectura single-end i una per a dades paired-end.
El següent pas lògic seria fer que el workflow accepti qualsevol tipus de dades sobre la marxa, que està fora de l'abast d'aquest curs, però podem abordar-ho en un seguiment.

---

### Conclusió

Sabeu com adaptar un workflow d'una sola mostra per paral·lelitzar el processament de múltiples mostres, generar un informe de QC complet i adaptar el workflow per utilitzar dades de lectura paired-end.

### Què segueix?

Doneu-vos una gran palmadeta a l'esquena! Heu completat el curs Nextflow per a RNAseq.

Aneu al [resum final del curs](./next_steps.md) per revisar el que heu après i descobrir què ve després.
