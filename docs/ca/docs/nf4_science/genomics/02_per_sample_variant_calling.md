# Part 2: Detecció de variants per mostra

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la Part 1, vau provar les comandes de Samtools i GATK manualment als seus respectius contenidors.
A continuació, encapsularem aquestes mateixes comandes en un workflow de Nextflow.

## Assignació

En aquesta part del curs, desenvoluparem un workflow que fa el següent:

1. Generar un fitxer d'índex per a cada fitxer BAM d'entrada utilitzant [Samtools](https://www.htslib.org/)
2. Executar GATK HaplotypeCaller a cada fitxer BAM d'entrada per generar deteccions de variants per mostra en format VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

Això automatitza els passos de la primera secció de [Part 1: Visió general del mètode](./01_method.md#1-per-sample-variant-calling), on vau executar aquestes comandes manualment als contenidors.

Com a punt de partida, us proporcionem un fitxer de workflow, `genomics.nf`, que esbossa les parts principals del workflow, així com dos fitxers de mòdul, `samtools_index.nf` i `gatk_haplotypecaller.nf`, que esbossen l'estructura de cada procés.

??? full-code "Fitxers base"

    ```groovy title="genomics.nf"
    #!/usr/bin/env nextflow

    // Declaracions INCLUDE de mòduls

    /*
     * Paràmetres del pipeline
     */

    // Entrada principal

    workflow {

        main:
        // Crear canal d'entrada

        // Cridar processos

        publish:
        // Declarar sortides a publicar
    }

    output {
        // Configurar destinacions de publicació
    }
    ```

    ```groovy title="modules/samtools_index.nf"
    #!/usr/bin/env nextflow

    /*
     * Generar fitxer d'índex BAM
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/gatk_haplotypecaller.nf"
    #!/usr/bin/env nextflow

    /*
     * Detectar variants amb GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

Aquests fitxers no són funcionals; el seu propòsit és només servir com a base perquè ompliu les parts interessants del codi.

## Pla de la lliçó

Per fer el procés de desenvolupament més educatiu, ho hem dividit en quatre etapes:

1. **Escriure un workflow d'una sola etapa que executi l'índex de Samtools en un fitxer BAM.**
   Això cobreix la creació d'un mòdul, la seva importació i la seva crida en un workflow.
2. **Afegir un segon procés per executar GATK HaplotypeCaller al fitxer BAM indexat.**
   Això introdueix l'encadenament de sortides de processos a entrades i la gestió de fitxers accessoris.
3. **Adaptar el workflow per executar-se en un lot de mostres.**
   Això cobreix l'execució paral·lela i introdueix les tuples per mantenir els fitxers associats junts.
4. **Fer que el workflow accepti un fitxer de text que contingui un lot de fitxers d'entrada.**
   Això demostra un patró comú per proporcionar entrades en bloc.

Cada pas se centra en un aspecte específic del desenvolupament de workflows.

!!! tip "Consell"

     Assegureu-vos que esteu al directori de treball correcte:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Escriure un workflow d'una sola etapa que executi l'índex de Samtools en un fitxer BAM

Aquest primer pas se centra en els conceptes bàsics: carregar un fitxer BAM i generar un índex per a ell.

Recordeu la comanda `samtools index` de la [Part 1](01_method.md):

```bash
samtools index '<input_bam>'
```

La comanda pren un fitxer BAM com a entrada i produeix un fitxer d'índex `.bai` al seu costat.
L'URI del contenidor era `community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464`.

Encapsularem aquesta informació en Nextflow en tres etapes:

1. Configurar l'entrada
2. Escriure el procés d'indexació i cridar-lo al workflow
3. Configurar la gestió de sortides

### 1.1. Configurar l'entrada

Necessitem declarar un paràmetre d'entrada, crear un perfil de prova per proporcionar un valor per defecte convenient i crear un canal d'entrada.

#### 1.1.1. Afegir una declaració de paràmetre d'entrada

Al fitxer principal del workflow `genomics.nf`, sota la secció `Pipeline parameters`, declareu un paràmetre CLI anomenat `input`.

=== "Després"

    ```groovy title="genomics.nf" linenums="5" hl_lines="4-7"
    /*
     * Paràmetres del pipeline
     */
    params {
        // Entrada principal
        input: Path
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="5"
    /*
     * Paràmetres del pipeline
     */

    // Entrada principal
    ```

Això configura el paràmetre CLI, però no volem escriure la ruta del fitxer cada vegada que executem el workflow durant el desenvolupament.
Hi ha múltiples opcions per proporcionar un valor per defecte; aquí utilitzem un perfil de prova.

#### 1.1.2. Crear un perfil de prova amb un valor per defecte a `nextflow.config`

Un perfil de prova proporciona valors per defecte convenients per provar un workflow sense especificar entrades a la línia de comandes.
Aquesta és una convenció comuna a l'ecosistema Nextflow (vegeu [Hello Config](../../hello_nextflow/06_hello_config.md) per a més detalls).

Afegiu un bloc `profiles` a `nextflow.config` amb un perfil `test` que estableixi el paràmetre `input` a un dels fitxers BAM de prova.

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/bam/reads_mother.bam"
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aquí, estem utilitzant `${projectDir}`, una variable integrada de Nextflow que apunta al directori on es troba l'script del workflow.
Això facilita la referència a fitxers de dades i altres recursos sense codificar rutes absolutes.

#### 1.1.3. Configurar el canal d'entrada

Al bloc workflow, creeu un canal d'entrada a partir del valor del paràmetre utilitzant la factoria de canals `.fromPath` (com s'utilitza a [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Després"

    ```groovy title="genomics.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Crear canal d'entrada (fitxer únic via paràmetre CLI)
        reads_ch = channel.fromPath(params.input)
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="13"
    workflow {

        main:
        // Crear canal d'entrada
    ```

A continuació, necessitarem crear el procés per executar la indexació en aquesta entrada.

### 1.2. Escriure el procés d'indexació i cridar-lo al workflow

Necessitem escriure la definició del procés al fitxer del mòdul, importar-lo al workflow utilitzant una declaració include i cridar-lo a l'entrada.

#### 1.2.1. Omplir el mòdul per al procés d'indexació

Obriu `modules/samtools_index.nf` i examineu l'esquema de la definició del procés.
Hauríeu de reconèixer els elements estructurals principals; si no, considereu llegir [Hello Nextflow](../../hello_nextflow/01_hello_world.md) per refrescar la memòria.

Aneu endavant i ompliu la definició del procés per vosaltres mateixos utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball amb la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/samtools_index.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Generar fitxer d'índex BAM
     */
    process SAMTOOLS_INDEX {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Després"

    ```groovy title="modules/samtools_index.nf" linenums="1" hl_lines="8 11 14 18"
    #!/usr/bin/env nextflow

    /*
     * Generar fitxer d'índex BAM
     */
    process SAMTOOLS_INDEX {

        container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

        input:
        path input_bam

        output:
        path "${input_bam}.bai"

        script:
        """
        samtools index '$input_bam'
        """
    }
    ```

Un cop hàgiu completat això, el procés està complet.
Per utilitzar-lo al workflow, necessitareu importar el mòdul i afegir una crida al procés.

#### 1.2.2. Incloure el mòdul

A `genomics.nf`, afegiu una declaració `include` per fer el procés disponible al workflow:

=== "Després"

    ```groovy title="genomics.nf" linenums="3" hl_lines="2"
    // Declaracions INCLUDE de mòduls
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="3"
    // Declaracions INCLUDE de mòduls
    ```

El procés ara està disponible a l'àmbit del workflow.

#### 1.2.3. Cridar el procés d'indexació a l'entrada

Ara, afegim una crida a `SAMTOOLS_INDEX` al bloc workflow, passant el canal d'entrada com a argument.

=== "Després"

    ```groovy title="genomics.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Crear canal d'entrada (fitxer únic via paràmetre CLI)
        reads_ch = channel.fromPath(params.input)

        // Crear fitxer d'índex per al fitxer BAM d'entrada
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="14"
    workflow {

        main:
        // Crear canal d'entrada (fitxer únic via paràmetre CLI)
        reads_ch = channel.fromPath(params.input)

        // Cridar processos
    ```

El workflow ara carrega l'entrada i executa el procés d'indexació sobre ella.
A continuació, necessitem configurar com es publica la sortida.

### 1.3. Configurar la gestió de sortides

Necessitem declarar quines sortides de processos publicar i especificar on han d'anar.

#### 1.3.1. Declarar una sortida a la secció `publish:`

La secció `publish:` dins del bloc workflow declara quines sortides de processos s'han de publicar.
Assigneu la sortida de `SAMTOOLS_INDEX` a una destinació anomenada `bam_index`.

=== "Després"

    ```groovy title="genomics.nf" linenums="22" hl_lines="2"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="22"
        publish:
        // Declarar sortides a publicar
    }
    ```

A continuació, necessitarem dir a Nextflow on posar la sortida publicada.

#### 1.3.2. Configurar la destinació de sortida al bloc `output {}`

El bloc `output {}` es troba fora del workflow i especifica on es publica cada destinació anomenada.
Afegim una destinació per a `bam_index` que publiqui en un subdirectori `bam/`.

=== "Després"

    ```groovy title="genomics.nf" linenums="26" hl_lines="2-4"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="26"
    output {
        // Configurar destinacions de publicació
    }
    ```

!!! note "Nota"

    Per defecte, Nextflow publica fitxers de sortida com a enllaços simbòlics, cosa que evita duplicacions innecessàries.
    Tot i que els fitxers de dades que estem utilitzant aquí són molt petits, en genòmica poden ser molt grans.
    Els enllaços simbòlics es trenquen quan netegeu el vostre directori `work`, així que per a workflows de producció potser voldreu sobreescriure el mode de publicació per defecte a `'copy'`.

### 1.4. Executar el workflow

En aquest punt, tenim un workflow d'indexació d'un sol pas que hauria de ser completament funcional. Provem que funciona!

Podem executar-lo amb `-profile test` per utilitzar el valor per defecte configurat al perfil de prova i evitar haver d'escriure la ruta a la línia de comandes.

```bash
nextflow run genomics.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Podeu comprovar que el fitxer d'índex s'ha generat correctament mirant al directori de treball o al directori de resultats.

??? abstract "Contingut del directori work"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Contingut del directori de resultats"

    ```console
    results/
    └── bam/
        └── reads_mother.bam.bai -> ...
    ```

Aquí està!

### Conclusió

Sabeu com crear un mòdul que conté un procés, importar-lo a un workflow, cridar-lo amb un canal d'entrada i publicar els resultats.

### Què segueix?

Afegir un segon pas que prengui la sortida del procés d'indexació i l'utilitzi per executar la detecció de variants.

---

## 2. Afegir un segon procés per executar GATK HaplotypeCaller al fitxer BAM indexat

Ara que tenim un índex per al nostre fitxer d'entrada, podem passar a configurar el pas de detecció de variants.

Recordeu la comanda `gatk HaplotypeCaller` de la [Part 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

La comanda pren un fitxer BAM (`-I`), un genoma de referència (`-R`) i un fitxer d'intervals (`-L`), i produeix un fitxer VCF (`-O`) juntament amb el seu índex.
L'eina també espera que l'índex BAM, l'índex de referència i el diccionari de referència estiguin col·locats amb els seus respectius fitxers.
L'URI del contenidor era `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

Seguim les mateixes tres etapes que abans:

1. Configurar les entrades
2. Escriure el procés de detecció de variants i cridar-lo al workflow
3. Configurar la gestió de sortides

### 2.1. Configurar les entrades

El pas de detecció de variants requereix diversos fitxers d'entrada addicionals.
Necessitem declarar paràmetres per a ells, afegir valors per defecte al perfil de prova i crear variables per carregar-los.

#### 2.1.1. Afegir declaracions de paràmetres per a entrades accessòries

Com que el nostre nou procés espera un grapat de fitxers addicionals, afegiu declaracions de paràmetres per a ells a `genomics.nf` sota la secció `Pipeline parameters`:

=== "Després"

    ```groovy title="genomics.nf" linenums="9" hl_lines="5-9"
    params {
        // Entrada principal
        input: Path

        // Fitxers accessoris
        reference: Path
        reference_index: Path
        reference_dict: Path
        intervals: Path
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="9"
    params {
        // Entrada principal
        input: Path
    }
    ```

Com abans, proporcionarem valors per defecte a través del perfil de prova en lloc d'en línia.

#### 2.1.2. Afegir valors per defecte de fitxers accessoris al perfil de prova

Igual que vam fer per a `input` a la secció 1.1.2, afegiu valors per defecte per als fitxers accessoris al perfil de prova a `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="4" hl_lines="3-6"
    test {
        params.input = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.input = "${projectDir}/data/bam/reads_mother.bam"
    }
    ```

A continuació, necessitarem crear variables que carreguin aquestes rutes de fitxers per utilitzar-les al workflow.

#### 2.1.3. Crear variables per als fitxers accessoris

Afegiu variables per a les rutes dels fitxers accessoris dins del bloc workflow:

=== "Després"

    ```groovy title="genomics.nf" linenums="21" hl_lines="7-11"
    workflow {

        main:
        // Crear canal d'entrada (fitxer únic via paràmetre CLI)
        reads_ch = channel.fromPath(params.input)

        // Carregar les rutes de fitxers per als fitxers accessoris (referència i intervals)
        ref_file        = file(params.reference)
        ref_index_file  = file(params.reference_index)
        ref_dict_file   = file(params.reference_dict)
        intervals_file  = file(params.intervals)

        // Crear fitxer d'índex per al fitxer BAM d'entrada
        SAMTOOLS_INDEX(reads_ch)
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="21"
    workflow {

        main:
        // Crear canal d'entrada (fitxer únic via paràmetre CLI)
        reads_ch = channel.fromPath(params.input)

        // Crear fitxer d'índex per al fitxer BAM d'entrada
        SAMTOOLS_INDEX(reads_ch)
    ```

La sintaxi `file()` indica explícitament a Nextflow que gestioni aquestes entrades com a rutes de fitxers.
Podeu aprendre més sobre això a la Side Quest [Working with files](../../side_quests/working_with_files.md).

### 2.2. Escriure el procés de detecció de variants i cridar-lo al workflow

Necessitem escriure la definició del procés al fitxer del mòdul, importar-lo al workflow utilitzant una declaració include i cridar-lo a les lectures d'entrada més la sortida del pas d'indexació i els fitxers accessoris.

#### 2.2.1. Omplir el mòdul per al procés de detecció de variants

Obriu `modules/gatk_haplotypecaller.nf` i examineu l'esquema de la definició del procés.

Aneu endavant i ompliu la definició del procés per vosaltres mateixos utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball amb la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Detectar variants amb GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Després"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="1" hl_lines="8 11-16 19-20 24-28"
    #!/usr/bin/env nextflow

    /*
     * Detectar variants amb GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path input_bam
        path input_bam_index
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

        script:
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    }
    ```

Notareu que aquest procés té més entrades de les que la comanda GATK en si requereix.
GATK sap buscar el fitxer d'índex BAM i els fitxers accessoris del genoma de referència basant-se en convencions de nomenclatura, però Nextflow és agnòstic del domini i no coneix aquestes convencions.
Necessitem llistar-los explícitament perquè Nextflow els prepari al directori de treball en temps d'execució; en cas contrari, GATK llançarà un error sobre fitxers que falten.

De manera similar, llistem explícitament el fitxer d'índex del VCF de sortida (`"${input_bam}.vcf.idx"`) perquè Nextflow en faci un seguiment per a passos posteriors.
Utilitzem la sintaxi `emit:` per assignar un nom a cada canal de sortida, cosa que serà útil quan connectem les sortides al bloc publish.

Un cop hàgiu completat això, el procés està complet.
Per utilitzar-lo al workflow, necessitareu importar el mòdul i afegir una crida al procés.

#### 2.2.2. Importar el nou mòdul

Actualitzeu `genomics.nf` per importar el nou mòdul:

=== "Després"

    ```groovy title="genomics.nf" linenums="3" hl_lines="3"
    // Declaracions INCLUDE de mòduls
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="3"
    // Declaracions INCLUDE de mòduls
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    ```

El procés ara està disponible a l'àmbit del workflow.

#### 2.2.3. Afegir la crida al procés

Afegiu la crida al procés al cos del workflow, sota `main:`:

=== "Després"

    ```groovy title="genomics.nf" linenums="33" hl_lines="4-12"
        // Crear fitxer d'índex per al fitxer BAM d'entrada
        SAMTOOLS_INDEX(reads_ch)

        // Detectar variants del fitxer BAM indexat
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
            ref_file,
            ref_index_file,
            ref_dict_file,
            intervals_file
        )
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="33"
        // Crear fitxer d'índex per al fitxer BAM d'entrada
        SAMTOOLS_INDEX(reads_ch)
    ```

Hauríeu de reconèixer la sintaxi `*.out` de la sèrie de formació Hello Nextflow; estem dient a Nextflow que prengui el canal de sortida de `SAMTOOLS_INDEX` i el connecti a la crida del procés `GATK_HAPLOTYPECALLER`.

!!! note "Nota"

    Noteu que les entrades es proporcionen exactament en el mateix ordre a la crida al procés que estan llistades al bloc d'entrada del procés.
    A Nextflow, les entrades són posicionals, el que significa que _heu_ de seguir el mateix ordre; i per descomptat hi ha d'haver el mateix nombre d'elements.

### 2.3. Configurar la gestió de sortides

Necessitem afegir les noves sortides a la declaració publish i configurar on van.

#### 2.3.1. Afegir destinacions de publicació per a les sortides de detecció de variants

Afegiu les sortides VCF i índex a la secció `publish:`:

=== "Després"

    ```groovy title="genomics.nf" linenums="45" hl_lines="3-4"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="45"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    }
    ```

A continuació, necessitarem dir a Nextflow on posar les noves sortides.

#### 2.3.2. Configurar les noves destinacions de sortida

Afegiu entrades per a les destinacions `vcf` i `vcf_idx` al bloc `output {}`, publicant ambdues en un subdirectori `vcf/`:

=== "Després"

    ```groovy title="genomics.nf" linenums="51" hl_lines="5-10"
    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="49"
    output {
        bam_index {
            path 'bam'
        }
    }
    ```

El VCF i el seu índex es publiquen com a destinacions separades que ambdues van al subdirectori `vcf/`.

### 2.4. Executar el workflow

Executeu el workflow ampliat, afegint `-resume` aquesta vegada perquè no hàgim de tornar a executar el pas d'indexació.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Ara, si mirem la sortida de la consola, veiem els dos processos llistats.

El primer procés es va ometre gràcies a la memòria cau, com s'esperava, mentre que el segon procés es va executar ja que és nou.

Trobareu els fitxers de sortida al directori de resultats (com a enllaços simbòlics al directori de treball).

??? abstract "Contingut del directori"

    ```console
    results/
    ├── bam/
    │   └── reads_mother.bam.bai -> ...
    └── vcf/
        ├── reads_mother.bam.vcf -> ...
        └── reads_mother.bam.vcf.idx -> ...
    ```

Si obriu el fitxer VCF, hauríeu de veure el mateix contingut que al fitxer que vau generar executant la comanda GATK directament al contenidor.

??? abstract "Contingut del fitxer"

    ```console title="reads_mother.bam.vcf" linenums="26"
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
    20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
    20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
    20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
    ```

Aquesta és la sortida que ens interessa generar per a cada mostra del nostre estudi.

### Conclusió

Sabeu com fer un workflow modular de dos passos que fa treball d'anàlisi real i és capaç de gestionar idiosincràsies de formats de fitxers genòmics com els fitxers accessoris.

### Què segueix?

Fer que el workflow gestioni múltiples mostres en bloc.

---

## 3. Adaptar el workflow per executar-se en un lot de mostres

Està molt bé tenir un workflow que pugui automatitzar el processament d'una sola mostra, però què passa si teniu 1000 mostres?
Necessiteu escriure un script bash que faci un bucle per totes les vostres mostres?

No, gràcies a Déu! Només feu un petit ajust al codi i Nextflow també ho gestionarà per vosaltres.

### 3.1. Actualitzar l'entrada per llistar tres mostres

Per executar-se en múltiples mostres, actualitzeu el perfil de prova per proporcionar una matriu de rutes de fitxers en lloc d'una sola.
Aquesta és una manera ràpida de provar l'execució multi-mostra; al següent pas canviarem a un enfocament més escalable utilitzant un fitxer d'entrades.

Primer, comenteu l'anotació de tipus a la declaració del paràmetre, ja que les matrius no poden utilitzar declaracions tipades:

=== "Després"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Entrada principal (matriu de tres mostres)
        input //: Path
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="10"
        // Entrada principal
        input: Path
    ```

Després actualitzeu el perfil de prova per llistar les tres mostres:

=== "Després"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2-6"
    test {
        params.input = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.input = "${projectDir}/data/bam/reads_mother.bam"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

La factoria de canals al cos del workflow (`.fromPath`) accepta múltiples rutes de fitxers igual que una sola, així que no calen altres canvis.

### 3.2. Executar el workflow

Proveu d'executar el workflow ara que la infraestructura està configurada per executar-se en les tres mostres de prova.

```bash
nextflow run genomics.nf -profile test -resume
```

Cosa curiosa: això _pot funcionar_, O _pot fallar_. Per exemple, aquí hi ha una execució que va tenir èxit:

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Si la vostra execució del workflow va tenir èxit, executeu-la de nou fins que obtingueu un error com aquest:

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

    executor >  local (4)
    [01/eea165] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 1 ✔
    [a5/fa9fd0] GATK_HAPLOTYPECALLER (3) | 1 of 3, cached: 1
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Process `GATK_HAPLOTYPECALLER (2)` terminated with an error exit status (2)

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_father.bam         -O reads_father.bam.vcf         -L intervals.bed

    Command exit status:
      2

    Command error:
      ...
      A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
      ...
    ```

Si mireu la sortida d'error de la comanda GATK, hi haurà una línia com aquesta:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Bé, això és estrany, considerant que vam indexar explícitament els fitxers BAM al primer pas del workflow. Podria haver-hi alguna cosa malament amb la infraestructura?

### 3.3. Resoldre el problema

Inspeccionarem els directoris de treball i utilitzarem l'operador `view()` per esbrinar què va anar malament.

#### 3.3.1. Comprovar els directoris de treball per a les crides rellevants

Doneu una ullada dins del directori de treball per a la crida del procés `GATK_HAPLOTYPECALLER` fallida llistada a la sortida de la consola.

??? abstract "Contingut del directori"

    ```console
    work/a5/fa9fd0994b6beede5fb9ea073596c2
    ├── intervals.bed -> /workspaces/training/nf4-science/genomics/data/ref/intervals.bed
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/01/eea16597bd6e810fb4cf89e60f8c2d/reads_father.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    ├── reads_son.bam.vcf
    ├── reads_son.bam.vcf.idx
    ├── ref.dict -> /workspaces/training/nf4-science/genomics/data/ref/ref.dict
    ├── ref.fasta -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta
    └── ref.fasta.fai -> /workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai
    ```

Presteu especial atenció als noms del fitxer BAM i l'índex BAM que estan llistats en aquest directori: `reads_son.bam` i `reads_father.bam.bai`.

Què dimonis? Nextflow ha preparat un fitxer d'índex al directori de treball d'aquesta crida de procés, però és l'incorrecte. Com pot haver passat això?

#### 3.3.2. Utilitzar l'[operador view()](https://www.nextflow.io/docs/latest/reference/operator.html#view) per inspeccionar el contingut dels canals

Afegiu aquestes dues línies al cos del workflow abans de la crida al procés `GATK_HAPLOTYPECALLER` per veure el contingut del canal:

=== "Després"

    ```groovy title="genomics.nf" linenums="36" hl_lines="3-5"
        SAMTOOLS_INDEX(reads_ch)

        // diagnòstics temporals
        reads_ch.view()
        SAMTOOLS_INDEX.out.view()

        // Detectar variants del fitxer BAM indexat
        GATK_HAPLOTYPECALLER(
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="36"
        SAMTOOLS_INDEX(reads_ch)

        // Detectar variants del fitxer BAM indexat
        GATK_HAPLOTYPECALLER(
    ```

Després executeu la comanda del workflow de nou.

```bash
nextflow run genomics.nf -profile test
```

Un cop més, això pot tenir èxit o fallar. Aquí està el que sembla la sortida de les dues crides `.view()` per a una execució fallida:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Les tres primeres línies corresponen al canal d'entrada i la segona, al canal de sortida.
Podeu veure que els fitxers BAM i els fitxers d'índex per a les tres mostres no estan llistats en el mateix ordre!

!!! note "Nota"

    Quan crideu un procés de Nextflow en un canal que conté múltiples elements, Nextflow intentarà paral·lelitzar l'execució tant com sigui possible i recollirà les sortides en l'ordre en què estiguin disponibles.
    La conseqüència és que les sortides corresponents poden ser recollides en un ordre diferent del que es van proporcionar les entrades originals.

Tal com està escrit actualment, el nostre script de workflow assumeix que els fitxers d'índex sortiran del pas d'indexació llistats en el mateix ordre mare/pare/fill que es van donar les entrades.
Però això no està garantit que sigui el cas, per això de vegades (encara que no sempre) els fitxers incorrectes s'aparellen al segon pas.

Per solucionar això, necessitem assegurar-nos que els fitxers BAM i els seus fitxers d'índex viatgin junts pels canals.

!!! tip "Consell"

    Les declaracions `view()` al codi del workflow no fan res, així que no és un problema deixar-les.
    Tanmateix, desordenen la vostra sortida de consola, així que recomanem eliminar-les quan hàgiu acabat de resoldre el problema.

### 3.4. Actualitzar el workflow per gestionar correctament els fitxers d'índex

La solució és agrupar cada fitxer BAM amb el seu índex en una tupla, després actualitzar el procés posterior i la infraestructura del workflow perquè coincideixin.

#### 3.4.1. Canviar la sortida del mòdul SAMTOOLS_INDEX a una tupla

La manera més senzilla d'assegurar que un fitxer BAM i el seu índex es mantinguin estretament associats és empaquetar-los junts en una tupla que surt de la tasca d'índex.

!!! note "Nota"

    Una **tupla** és una llista finita i ordenada d'elements que s'utilitza comunament per retornar múltiples valors d'una funció. Les tuples són particularment útils per passar múltiples entrades o sortides entre processos mentre es preserva la seva associació i ordre.

Actualitzeu la sortida a `modules/samtools_index.nf` per incloure el fitxer BAM:

=== "Després"

    ```groovy title="modules/samtools_index.nf" linenums="14" hl_lines="2"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Abans"

    ```groovy title="modules/samtools_index.nf" linenums="14"
        output:
        path "${input_bam}.bai"
    ```

D'aquesta manera, cada fitxer d'índex estarà estretament acoblat amb el seu fitxer BAM original, i la sortida global del pas d'indexació serà un sol canal que conté parells de fitxers.

#### 3.4.2. Canviar l'entrada del mòdul GATK_HAPLOTYPECALLER per acceptar una tupla

Com que hem canviat la 'forma' de la sortida del primer procés, necessitem actualitzar la definició d'entrada del segon procés perquè coincideixi.

Actualitzeu `modules/gatk_haplotypecaller.nf`:

=== "Després"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11" hl_lines="2"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Abans"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="11"
        input:
        path input_bam
        path input_bam_index
    ```

A continuació, necessitarem actualitzar el workflow per reflectir la nova estructura de tupla a la crida del procés i les destinacions de publicació.

#### 3.4.3. Actualitzar la crida a GATK_HAPLOTYPECALLER al workflow

Ja no necessitem proporcionar el `reads_ch` original al procés `GATK_HAPLOTYPECALLER`, ja que el fitxer BAM ara està agrupat al canal de sortida per `SAMTOOLS_INDEX`.

Actualitzeu la crida a `genomics.nf`:

=== "Després"

    ```groovy title="genomics.nf" linenums="42" hl_lines="2"
        GATK_HAPLOTYPECALLER(
            SAMTOOLS_INDEX.out,
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="42"
        GATK_HAPLOTYPECALLER(
            reads_ch,
            SAMTOOLS_INDEX.out,
    ```

Finalment, necessitem actualitzar les destinacions de publicació per reflectir la nova estructura de sortida.

#### 3.4.4. Actualitzar la destinació de publicació per a la sortida BAM indexada

Com que la sortida de SAMTOOLS_INDEX ara és una tupla que conté tant el fitxer BAM com el seu índex, canvieu el nom de la destinació de publicació de `bam_index` a `indexed_bam` per reflectir millor el seu contingut:

=== "Després"

    ```groovy title="genomics.nf" linenums="46" hl_lines="2 8"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="46"
        publish:
        bam_index = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    }

    output {
        bam_index {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Amb aquests canvis, el BAM i el seu índex estan garantits de viatjar junts, així que l'aparellament sempre serà correcte.

### 3.5. Executar el workflow corregit

Executeu el workflow de nou per assegurar-vos que això funcionarà de manera fiable d'ara endavant.

```bash
nextflow run genomics.nf -profile test
```

Aquesta vegada (i cada vegada) tot hauria de funcionar correctament:

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

El directori de resultats ara conté tant fitxers BAM com BAI per a cada mostra (de la tupla), juntament amb les sortides VCF:

??? abstract "Contingut del directori de resultats"

    ```console
    results/
    ├── bam/
    │   ├── reads_father.bam -> ...
    │   ├── reads_father.bam.bai -> ...
    │   ├── reads_mother.bam -> ...
    │   ├── reads_mother.bam.bai -> ...
    │   ├── reads_son.bam -> ...
    │   └── reads_son.bam.bai -> ...
    └── vcf/
        ├── reads_father.bam.vcf -> ...
        ├── reads_father.bam.vcf.idx -> ...
        ├── reads_mother.bam.vcf -> ...
        ├── reads_mother.bam.vcf.idx -> ...
        ├── reads_son.bam.vcf -> ...
        └── reads_son.bam.vcf.idx -> ...
    ```

Agrupant fitxers associats en tuples, vam assegurar que els fitxers correctes sempre viatgen junts pel workflow.
El workflow ara processa qualsevol nombre de mostres de manera fiable, però llistar-les individualment a la configuració no és molt escalable.
Al següent pas, canviarem a llegir entrades d'un fitxer.

### Conclusió

Sabeu com fer que el vostre workflow s'executi en múltiples mostres (independentment).

### Què segueix?

Fer-ho més fàcil per gestionar mostres en bloc.

---

## 4. Fer que el workflow accepti un fitxer de text que contingui un lot de fitxers d'entrada

Una manera molt comuna de proporcionar múltiples fitxers de dades d'entrada a un workflow és fer-ho amb un fitxer de text que contingui les rutes dels fitxers.
Pot ser tan simple com un fitxer de text que llisti una ruta de fitxer per línia i res més, o el fitxer pot contenir metadades addicionals, en aquest cas sovint s'anomena full de mostres (samplesheet).

Aquí us mostrarem com fer el cas simple.

### 4.1. Examinar el fitxer CSV proporcionat que llista les rutes dels fitxers d'entrada

Ja hem fet un fitxer CSV que llista les rutes dels fitxers d'entrada, anomenat `samplesheet.csv`, que podeu trobar al directori `data/`.

```txt title="samplesheet.csv"
sample_id,reads_bam
NA12878,/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
NA12877,/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
NA12882,/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Aquest fitxer CSV inclou una línia de capçalera que anomena les columnes.

!!! warning "Advertència"

    Les rutes de fitxers al CSV són rutes absolutes que han de coincidir amb el vostre entorn.
    Si no esteu executant això a l'entorn de formació que proporcionem, haureu d'actualitzar les rutes perquè coincideixin amb el vostre sistema.

### 4.2. Actualitzar el paràmetre i el perfil de prova

Canvieu el paràmetre `input` per apuntar al fitxer `samplesheet.csv` en lloc de llistar mostres individuals.

Restaureu l'anotació de tipus al bloc params (ja que és una sola ruta de nou):

=== "Després"

    ```groovy title="genomics.nf" linenums="10" hl_lines="1-2"
        // Entrada principal (fitxer de fitxers d'entrada, un per línia)
        input: Path
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="10"
        // Entrada principal (matriu de tres mostres)
        input
    ```

Després actualitzeu el perfil de prova per apuntar al fitxer de text:

=== "Després"

    ```groovy title="nextflow.config" linenums="4" hl_lines="2"
    test {
        params.input = "${projectDir}/data/samplesheet.csv"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.input = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

La llista de fitxers ja no viu al codi en absolut, cosa que és un gran pas en la direcció correcta.

### 4.3. Actualitzar la factoria de canals per analitzar l'entrada CSV

Actualment, la nostra factoria de canals d'entrada tracta qualsevol fitxer que li donem com les entrades de dades que volem alimentar al procés d'indexació.
Com que ara li estem donant un fitxer que llista rutes de fitxers d'entrada, necessitem canviar el seu comportament per analitzar el fitxer i tractar les rutes de fitxers que conté com les entrades de dades.

Podem fer això utilitzant el mateix patró que vam utilitzar a la [Part 2 de Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): aplicant l'operador [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) per analitzar el fitxer, després una operació `map` per extreure la ruta del fitxer de cada fila.

=== "Després"

    ```groovy title="genomics.nf" linenums="24" hl_lines="1-4"
        // Crear canal d'entrada des d'un fitxer CSV que llista rutes de fitxers d'entrada
        reads_ch = channel.fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> file(row.reads_bam) }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="24"
        // Crear canal d'entrada (fitxer únic via paràmetre CLI)
        reads_ch = channel.fromPath(params.input)
    ```

Una cosa que és nova en comparació amb el que vau trobar al curs Hello Nextflow és que aquest CSV té una línia de capçalera, així que afegim `#!groovy header: true` a la crida `splitCsv()`.
Això ens permet referenciar columnes per nom a l'operació `map`: `#!groovy row.reads_bam` extreu la ruta del fitxer de la columna `reads_bam` de cada fila.

!!! tip "Consell"

    Si no esteu segurs d'entendre què estan fent els operadors aquí, aquesta és una altra gran oportunitat per utilitzar l'operador `.view()` per veure com es veu el contingut del canal abans i després d'aplicar-los.

### 4.4. Executar el workflow

Executeu el workflow una vegada més.

```bash
nextflow run genomics.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    executor >  local (6)
    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

Això hauria de produir el mateix resultat que abans. El nostre workflow simple de detecció de variants ara té totes les característiques bàsiques que volíem.

### Conclusió

Sabeu com fer un workflow modular de múltiples passos per indexar un fitxer BAM i aplicar detecció de variants per mostra utilitzant GATK.

Més generalment, heu après com utilitzar components i lògica essencials de Nextflow per construir un pipeline genòmic simple que fa treball real, tenint en compte les idiosincràsies dels formats de fitxers genòmics i els requisits de les eines.

### Què segueix?

Feu una pausa! Això ha estat molt.

Quan us sentiu refrescats, aneu a la [Part 3](./03_joint_calling.md), on aprendreu com transformar aquest workflow simple de detecció de variants per mostra per aplicar detecció conjunta de variants a les dades.
