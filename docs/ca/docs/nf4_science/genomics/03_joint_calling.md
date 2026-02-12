# Part 3: Crida conjunta sobre una cohort

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Anteriorment, vau construir un pipeline de crida de variants per mostra que processava les dades de cada mostra de manera independent.
Ara l'ampliarem per implementar la crida conjunta de variants, tal com es descriu a la [Part 1](01_method.md) (cas d'ús multi-mostra).

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat la [Part 1: Visió general del mètode](./01_method.md), la [Part 2: Crida de variants per mostra](./02_per_sample_variant_calling.md) i teniu un pipeline `genomics.nf` funcional.

    Si no vau completar la Part 2 o voleu començar de nou per a aquesta part, podeu utilitzar la solució de la Part 2 com a punt de partida.
    Executeu aquestes comandes des del directori `nf4-science/genomics/`:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Això us proporciona un workflow complet de crida de variants per mostra.
    Podeu comprovar que s'executa correctament executant la comanda següent:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Assignació

En aquesta part del curs, ampliarem el workflow per fer el següent:

1. Generar un fitxer d'índex per a cada fitxer BAM d'entrada utilitzant Samtools
2. Executar GATK HaplotypeCaller sobre cada fitxer BAM d'entrada per generar un GVCF de crides de variants genòmiques per mostra
3. Recollir tots els GVCFs i combinar-los en un magatzem de dades GenomicsDB
4. Executar el genotipat conjunt sobre el magatzem de dades GVCF combinat per produir un VCF a nivell de cohort

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

Això implementa el mètode descrit a la [Part 1: Visió general del mètode](./01_method.md) (segona secció que cobreix el cas d'ús multi-mostra) i es construeix directament sobre el workflow produït a la [Part 2: Crida de variants per mostra](./02_per_sample_variant_calling.md).

## Pla de la lliçó

Ho hem dividit en dues etapes:

1. **Modificar el pas de crida de variants per mostra per produir un GVCF.**
   Això cobreix l'actualització de comandes i sortides del procés.
2. **Afegir un pas de genotipat conjunt que combina i genotipa els GVCFs per mostra.**
   Això introdueix l'operador `collect()`, closures de Groovy per a la construcció de línies de comandes, i processos amb múltiples comandes.

Això automatitza els passos de la segona secció de la [Part 1: Visió general del mètode](./01_method.md#2-joint-calling-on-a-cohort), on vau executar aquestes comandes manualment als seus contenidors.

!!! tip "Consell"

     Assegureu-vos que esteu al directori de treball correcte:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modificar el pas de crida de variants per mostra per produir un GVCF

El pipeline de la Part 2 produeix fitxers VCF, però la crida conjunta requereix fitxers GVCF.
Hem d'activar el mode de crida de variants GVCF i actualitzar l'extensió del fitxer de sortida.

Recordeu la comanda de crida de variants GVCF de la [Part 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Comparada amb la comanda base de HaplotypeCaller que vam encapsular a la Part 2, les diferències són el paràmetre `-ERC GVCF` i l'extensió de sortida `.g.vcf`.

### 1.1. Indicar a HaplotypeCaller que emeti un GVCF i actualitzar l'extensió de sortida

Obriu el fitxer de mòdul `modules/gatk_haplotypecaller.nf` per fer dos canvis:

- Afegiu el paràmetre `-ERC GVCF` a la comanda GATK HaplotypeCaller;
- Actualitzeu el camí del fitxer de sortida per utilitzar l'extensió `.g.vcf` corresponent, segons la convenció de GATK.

Assegureu-vos d'afegir una barra invertida (`\`) al final de la línia anterior quan afegiu `-ERC GVCF`.

=== "Després"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Abans"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

També hem d'actualitzar el bloc de sortida per coincidir amb la nova extensió de fitxer.
Com que hem canviat la sortida de la comanda de `.vcf` a `.g.vcf`, el bloc `output:` del procés ha de reflectir el mateix canvi.

### 1.2. Actualitzar l'extensió del fitxer de sortida al bloc de sortides del procés

=== "Després"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Abans"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

També hem d'actualitzar la configuració de publicació i sortida del workflow per reflectir les noves sortides GVCF.

### 1.3. Actualitzar els objectius de publicació per a les noves sortides GVCF

Com que ara estem produint GVCFs en lloc de VCFs, hauríem d'actualitzar la secció `publish:` del workflow per utilitzar noms més descriptius.
També organitzarem els fitxers GVCF en el seu propi subdirectori per claredat.

=== "Després"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Ara actualitzeu el bloc de sortida per coincidir.

### 1.4. Actualitzar el bloc de sortida per a la nova estructura de directoris

També hem d'actualitzar el bloc `output` per posar els fitxers GVCF en un subdirectori `gvcf`.

=== "Després"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="53"
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

Amb el mòdul, els objectius de publicació i el bloc de sortida tots actualitzats, podem provar els canvis.

### 1.5. Executar el pipeline

Executeu el workflow per verificar que els canvis funcionen.

```bash
nextflow run genomics.nf -profile test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

La sortida de Nextflow sembla la mateixa que abans, però els fitxers `.g.vcf` i els seus fitxers d'índex ara estan organitzats en subdirectoris.

??? abstract "Contingut del directori (enllaços simbòlics escurçats)"

    ```console
    results/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Si obriu un dels fitxers GVCF i el desplaceu, podeu verificar que GATK HaplotypeCaller ha produït fitxers GVCF tal com es va sol·licitar.

### Conclusió

Quan canvieu el nom del fitxer de sortida d'una comanda d'eina, el bloc `output:` del procés i la configuració de publicació/sortida s'han d'actualitzar per coincidir.

### Què segueix?

Apreneu a recollir el contingut d'un canal i passar-lo al procés següent com una única entrada.

---

## 2. Afegir un pas de genotipat conjunt

Ara hem de recollir els GVCFs per mostra, combinar-los en un magatzem de dades GenomicsDB i executar el genotipat conjunt per produir un VCF a nivell de cohort.
Tal com es va cobrir a la [Part 1](01_method.md), aquesta és una operació de dues eines: GenomicsDBImport combina els GVCFs, després GenotypeGVCFs produeix les crides de variants finals.
Encapsularem ambdues eines en un únic procés anomenat `GATK_JOINTGENOTYPING`.

Recordeu les dues comandes de la [Part 1](01_method.md):

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

La primera comanda pren els GVCFs per mostra i un fitxer d'intervals, i produeix un magatzem de dades GenomicsDB.
La segona pren aquest magatzem de dades, un genoma de referència, i produeix el VCF final a nivell de cohort.
L'URI del contenidor és el mateix que per a HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Configurar les entrades

El procés de genotipat conjunt necessita dos tipus d'entrades que encara no tenim: un nom de cohort arbitrari, i les sortides GVCF recollides de totes les mostres agrupades juntes.

#### 2.1.1. Afegir un paràmetre `cohort_name`

Hem de proporcionar un nom arbitrari per a la cohort.
Més endavant a la sèrie de formació aprendreu com utilitzar metadades de mostres per a aquest tipus de coses, però per ara simplement declarem un paràmetre CLI utilitzant `params`.

=== "Després"

    ```groovy title="genomics.nf" linenums="18" hl_lines="3-4"
        intervals: Path

        // Nom base per al fitxer de sortida final
        cohort_name: String
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="18"
        intervals: Path
    }
    ```

Utilitzarem aquest paràmetre per nomenar el fitxer de sortida final.

#### 2.1.2. Afegir un valor per defecte per a `cohort_name` al perfil de test

També afegim un valor per defecte per al paràmetre `cohort_name` al perfil de test:

=== "Després"

    ```groovy title="nextflow.config" linenums="4" hl_lines="7"
    test {
        params.input = "${projectDir}/data/samplesheet.csv"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
        params.cohort_name = "family_trio"
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="4"
    test {
        params.input = "${projectDir}/data/samplesheet.csv"
        params.reference = "${projectDir}/data/ref/ref.fasta"
        params.reference_index = "${projectDir}/data/ref/ref.fasta.fai"
        params.reference_dict = "${projectDir}/data/ref/ref.dict"
        params.intervals = "${projectDir}/data/ref/intervals.bed"
    }
    ```

A continuació, haurem de recollir les sortides per mostra perquè puguin ser processades juntes.

#### 2.1.3. Recollir les sortides de HaplotypeCaller entre mostres

Si connectéssim el canal de sortida de `GATK_HAPLOTYPECALLER` directament al nou procés, Nextflow cridaria el procés sobre cada GVCF de mostra per separat.
Volem agrupar tots tres GVCFs (i els seus fitxers d'índex) perquè Nextflow els passi tots junts a una única crida de procés.

Podem fer-ho utilitzant l'operador de canal `collect()`.
Afegiu les línies següents al cos del `workflow`, just després de la crida a GATK_HAPLOTYPECALLER:

=== "Després"

    ```groovy title="genomics.nf" linenums="48" hl_lines="4-6"
            intervals_file
        )

        // Recollir sortides de crida de variants entre mostres
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="48"
            intervals_file
        )
    ```

Desglossant això:

1. Prenem el canal de sortida de `GATK_HAPLOTYPECALLER` utilitzant la propietat `.out`.
2. Com que vam nomenar les sortides utilitzant `emit:` a la secció 1, podem seleccionar els GVCFs amb `.vcf` i els fitxers d'índex amb `.idx`. Sense sortides nomenades, hauríem d'utilitzar `.out[0]` i `.out[1]`.
3. L'operador `collect()` agrupa tots els fitxers en un únic element, així que `all_gvcfs_ch` conté tots tres GVCFs junts, i `all_idxs_ch` conté tots tres fitxers d'índex junts.

Podem recollir els GVCFs i els seus fitxers d'índex per separat (en lloc de mantenir-los junts en tuples) perquè Nextflow prepararà tots els fitxers d'entrada junts per a l'execució, així que els fitxers d'índex estaran presents al costat dels GVCFs.

!!! tip "Consell"

    Podeu utilitzar l'operador `view()` per inspeccionar el contingut dels canals abans i després d'aplicar operadors de canal.

### 2.2. Escriure el procés de genotipat conjunt i cridar-lo al workflow

Seguint el mateix patró que vam utilitzar a la Part 2, escriurem la definició del procés en un fitxer de mòdul, l'importarem al workflow i el cridarem sobre les entrades que acabem de preparar.

#### 2.2.1. Construir una cadena per donar a cada GVCF un argument `-V`

Abans de començar a omplir la definició del procés, hi ha una cosa a resoldre.
La comanda GenomicsDBImport espera un argument `-V` separat per a cada fitxer GVCF, així:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Si escrivíssim `-V ${all_gvcfs_ch}`, Nextflow simplement concatenaria els noms de fitxer i aquesta part de la comanda es veuria així:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Però necessitem que la cadena es vegi així:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Importantment, hem de construir aquesta cadena dinàmicament a partir de qualsevol fitxer que estigui al canal recollit.
Nextflow (via Groovy) proporciona una manera concisa de fer-ho:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Desglossant això:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` itera sobre cada camí de fitxer i afegeix `-V ` al davant, produint `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` els concatena amb espais: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. El resultat s'assigna a una variable local `gvcfs_line` (definida amb `def`), que podem interpolar a la plantilla de comanda.

Aquesta línia va dins del bloc `script:` del procés, abans de la plantilla de comanda.
Podeu col·locar codi Groovy arbitrari entre `script:` i l'obertura `"""` de la plantilla de comanda.

Llavors podreu referir-vos a tota aquesta cadena com `gvcfs_line` al bloc `script:` del procés.

#### 2.2.2. Omplir el mòdul per al procés de genotipat conjunt

A continuació, podem abordar l'escriptura del procés complet.

Obriu `modules/gatk_jointgenotyping.nf` i examineu l'esquema de la definició del procés.

Aneu endavant i ompliu la definició del procés utilitzant la informació proporcionada anteriorment, després comproveu el vostre treball contra la solució a la pestanya "Després" a continuació.

=== "Abans"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Combinar GVCFs en un magatzem de dades GenomicsDB i executar genotipat conjunt per produir crides a nivell de cohort
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Després"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Combinar GVCFs en un magatzem de dades GenomicsDB i executar genotipat conjunt per produir crides a nivell de cohort
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    }
    ```

Hi ha diverses coses que val la pena destacar aquí.

Com anteriorment, diverses entrades estan llistades encara que les comandes no les referencien directament: `all_idxs`, `ref_index` i `ref_dict`.
Llistar-les assegura que Nextflow prepara aquests fitxers al directori de treball al costat dels fitxers que sí apareixen a les comandes, que GATK espera trobar basant-se en convencions de nomenclatura.

La variable `gvcfs_line` utilitza la closure de Groovy descrita anteriorment per construir els arguments `-V` per a GenomicsDBImport.

Aquest procés executa dues comandes en sèrie, tal com ho faríeu al terminal.
GenomicsDBImport combina els GVCFs per mostra en un magatzem de dades, després GenotypeGVCFs llegeix aquest magatzem de dades i produeix el VCF final a nivell de cohort.
El magatzem de dades GenomicsDB (`${cohort_name}_gdb`) és un artefacte intermedi utilitzat només dins del procés; no apareix al bloc de sortida.

Un cop hàgiu completat això, el procés està llest per utilitzar.
Per utilitzar-lo al workflow, haureu d'importar el mòdul i afegir una crida de procés.

#### 2.2.3. Importar el mòdul

Afegiu la declaració d'importació a `genomics.nf`, sota les declaracions d'importació existents:

=== "Després"

    ```groovy title="genomics.nf" linenums="3" hl_lines="4"
    // Declaracions INCLUDE de mòduls
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="3"
    // Declaracions INCLUDE de mòduls
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

El procés ara està disponible a l'àmbit del workflow.

#### 2.2.4. Afegir la crida del procés

Afegiu la crida a `GATK_JOINTGENOTYPING` al cos del workflow, després de les línies `collect()`:

=== "Després"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combinar GVCFs en un magatzem de dades GenomicsDB i aplicar genotipat conjunt
        GATK_JOINTGENOTYPING(
            all_gvcfs_ch,
            all_idxs_ch,
            intervals_file,
            params.cohort_name,
            ref_file,
            ref_index_file,
            ref_dict_file
        )
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="53"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

El procés ara està completament connectat.
A continuació, configurem com es publiquen les sortides.

### 2.3. Configurar la gestió de sortides

Hem de publicar les sortides del VCF conjunt.
Afegiu objectius de publicació i entrades del bloc de sortida per als resultats del genotipat conjunt.

#### 2.3.1. Afegir objectius de publicació per al VCF conjunt

Afegiu el VCF conjunt i el seu índex a la secció `publish:` del workflow:

=== "Després"

    ```groovy title="genomics.nf" linenums="66" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="66"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Ara actualitzeu el bloc de sortida per coincidir.

#### 2.3.2. Afegir entrades del bloc de sortida per al VCF conjunt

Afegiu entrades per als fitxers VCF conjunts.
Els posarem a l'arrel del directori de resultats ja que aquesta és la sortida final.

=== "Després"

    ```groovy title="genomics.nf" linenums="74" hl_lines="11-16"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Abans"

    ```groovy title="genomics.nf" linenums="74"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

Amb el procés, els objectius de publicació i el bloc de sortida tots al seu lloc, podem provar el workflow complet.

### 2.4. Executar el workflow

Executeu el workflow per verificar que tot funciona.

```bash
nextflow run genomics.nf -profile test -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Els dos primers passos estan en memòria cau de l'execució anterior, i el nou pas `GATK_JOINTGENOTYPING` s'executa una vegada sobre les entrades recollides de totes tres mostres.
El fitxer de sortida final, `family_trio.joint.vcf` (i el seu índex), estan al directori de resultats.

??? abstract "Contingut del directori (enllaços simbòlics escurçats)"

    ```console
    results/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Si obriu el fitxer VCF conjunt, podeu verificar que el workflow ha produït les crides de variants esperades.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Ara teniu un workflow de crida conjunta de variants automatitzat i completament reproduïble!

!!! note "Nota"

    Tingueu en compte que els fitxers de dades que us vam proporcionar cobreixen només una petita porció del cromosoma 20.
    La mida real d'un conjunt de crides de variants es comptaria en milions de variants.
    Per això utilitzem només petits subconjunts de dades per a propòsits de formació!

### Conclusió

Sabeu com recollir sortides d'un canal i agrupar-les com una única entrada a un altre procés.
També sabeu com construir una línia de comandes utilitzant closures de Groovy, i com executar múltiples comandes en un únic procés.

### Què segueix?

Doneu-vos una gran palmadeta a l'esquena! Heu completat el curs Nextflow per a Genòmica.

Aneu al [resum final del curs](./next_steps.md) per revisar el que heu après i descobrir què ve després.
