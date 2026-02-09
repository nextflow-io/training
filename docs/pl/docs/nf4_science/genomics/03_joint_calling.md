# Część 3: Wspólne wyszukiwanie wariantów w kohorcie

W Części 2 zbudowałeś/-aś pipeline wyszukiwania wariantów dla pojedynczych próbek, który przetwarzał dane każdej próbki niezależnie.
Teraz rozszerzymy go, aby zaimplementować wspólne wyszukiwanie wariantów, jak omówiono w [Części 1](01_method.md).

## Zadanie

W tej części kursu rozszerzymy workflow'a, aby wykonywał następujące operacje:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Wygeneruj plik indeksu dla każdego pliku BAM wejściowego przy użyciu Samtools
2. Uruchom GATK HaplotypeCaller na każdym pliku BAM wejściowym, aby wygenerować GVCF z wykrytymi wariantami genomowymi dla pojedynczej próbki
3. Zbierz wszystkie pliki GVCF i połącz je w repozytorium danych GenomicsDB
4. Uruchom wspólne genotypowanie na połączonym repozytorium danych GVCF, aby wygenerować plik VCF na poziomie kohorty

Ta część bazuje bezpośrednio na workflow'ie stworzonym w Części 2.

??? info "Jak rozpocząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś/-aś [Część 2: Wyszukiwanie wariantów dla pojedynczych próbek](./02_per_sample_variant_calling.md) i masz działający pipeline `genomics.nf`.

    Jeśli nie ukończyłeś/-aś Części 2 lub chcesz zacząć od początku w tej części, możesz użyć rozwiązania z Części 2 jako punktu wyjścia.
    Uruchom te polecenia z wnętrza katalogu `nf4-science/genomics/`:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    To da Ci kompletny workflow wyszukiwania wariantów dla pojedynczych próbek.
    Możesz sprawdzić, czy działa poprawnie, uruchamiając następujące polecenie:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Plan lekcji

Podzieliliśmy to na dwa kroki:

1. **Zmodyfikuj krok wyszukiwania wariantów dla pojedynczych próbek, aby generował GVCF.**
   Obejmuje to aktualizację poleceń procesu i wyjść.
2. **Dodaj krok wspólnego genotypowania, który łączy i genotypuje pliki GVCF z pojedynczych próbek.**
   To wprowadza operator `collect()`, domknięcia Groovy do konstruowania linii poleceń oraz procesy z wieloma poleceniami.

!!! note "Uwaga"

     Upewnij się, że jesteś we właściwym katalogu roboczym:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Zmodyfikuj krok wyszukiwania wariantów dla pojedynczych próbek, aby generował GVCF

Pipeline z Części 2 generuje pliki VCF, ale wspólne wyszukiwanie wymaga plików GVCF.
Musimy włączyć tryb wyszukiwania wariantów GVCF i zaktualizować rozszerzenie pliku wyjściowego.

Przypomnij sobie polecenie wyszukiwania wariantów GVCF z [Części 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

W porównaniu z podstawowym poleceniem HaplotypeCaller, które opakował(a/i)śmy w Części 2, różnice to parametr `-ERC GVCF` i rozszerzenie wyjściowe `.g.vcf`.

### 1.1. Powiedz HaplotypeCaller, aby emitował GVCF i zaktualizuj rozszerzenie wyjściowe

Otwórz plik modułu `modules/gatk_haplotypecaller.nf`, aby wprowadzić dwie zmiany:

- Dodaj parametr `-ERC GVCF` do polecenia GATK HaplotypeCaller;
- Zaktualizuj ścieżkę pliku wyjściowego, aby używała odpowiedniego rozszerzenia `.g.vcf`, zgodnie z konwencją GATK.

Upewnij się, że dodajesz backslash (`\`) na końcu poprzedniej linii, gdy dodajesz `-ERC GVCF`.

=== "Po"

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

=== "Przed"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Musimy również zaktualizować blok wyjścia, aby odpowiadał nowemu rozszerzeniu pliku.
Ponieważ zmieniliśmy wyjście polecenia z `.vcf` na `.g.vcf`, blok `output:` procesu musi odzwierciedlać tę samą zmianę.

### 1.2. Zaktualizuj rozszerzenie pliku wyjściowego w bloku wyjść procesu

=== "Po"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Przed"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

Musimy również zaktualizować konfigurację publikowania i wyjścia workflow'a, aby odzwierciedlała nowe wyjścia GVCF.

### 1.3. Zaktualizuj cele publikowania dla nowych wyjść GVCF

Ponieważ teraz generujemy pliki GVCF zamiast VCF, powinniśmy zaktualizować sekcję `publish:` workflow'a, aby używać bardziej opisowych nazw.
Zorganizujemy również pliki GVCF w ich własnym podkatalogu dla przejrzystości.

=== "Po"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Teraz zaktualizuj blok wyjścia, aby odpowiadał.

### 1.4. Zaktualizuj blok wyjścia dla nowej struktury katalogów

Musimy również zaktualizować blok `output`, aby umieścić pliki GVCF w podkatalogu `gvcf`.

=== "Po"

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

=== "Przed"

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

Po zaktualizowaniu modułu, celów publikowania i bloku wyjścia możemy przetestować zmiany.

### 1.5. Uruchom pipeline'a

Uruchom workflow'a, aby sprawdzić, czy zmiany działają.

```bash
nextflow run genomics.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Wyjście Nextflow'a wygląda tak samo jak wcześniej, ale pliki `.g.vcf` i ich pliki indeksów są teraz zorganizowane w podkatalogach.

??? abstract "Zawartość katalogu (skrócone dowiązania symboliczne)"

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

Jeśli otworzysz jeden z plików GVCF i przewiniesz go, możesz sprawdzić, że GATK HaplotypeCaller wygenerował pliki GVCF zgodnie z żądaniem.

### Podsumowanie

Gdy zmienisz nazwę pliku wyjściowego polecenia narzędzia, blok `output:` procesu oraz konfiguracja publikowania/wyjścia muszą zostać zaktualizowane, aby odpowiadały.

### Co dalej?

Naucz się zbierać zawartość kanału i przekazywać ją do następnego procesu jako pojedyncze wejście.

---

## 2. Dodaj krok wspólnego genotypowania

Musimy teraz zebrać pliki GVCF z pojedynczych próbek, połączyć je w repozytorium danych GenomicsDB i uruchomić wspólne genotypowanie, aby wygenerować plik VCF na poziomie kohorty.
Jak omówiono w [Części 1](01_method.md), jest to operacja dwunarzędziowa: GenomicsDBImport łączy pliki GVCF, a następnie GenotypeGVCFs generuje ostateczne wykryte warianty.
Opakujemy oba narzędzia w jeden proces o nazwie `GATK_JOINTGENOTYPING`.

Przypomnij sobie dwa polecenia z [Części 1](01_method.md):

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

Pierwsze polecenie przyjmuje pliki GVCF z pojedynczych próbek i plik przedziałów, a generuje repozytorium danych GenomicsDB.
Drugie przyjmuje to repozytorium danych, genom referencyjny i generuje ostateczny plik VCF na poziomie kohorty.
URI kontenera jest taki sam jak dla HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Skonfiguruj wejścia

Proces wspólnego genotypowania potrzebuje dwóch rodzajów wejść, których jeszcze nie mamy: dowolnej nazwy kohorty i zebranych wyjść GVCF ze wszystkich próbek połączonych razem.

#### 2.1.1. Dodaj parametr `cohort_name`

Musimy podać dowolną nazwę dla kohorty.
Później w serii szkoleń nauczysz się, jak używać metadanych próbek do tego rodzaju rzeczy, ale na razie po prostu zadeklarujemy parametr CLI używając `params` i nadamy mu wartość domyślną dla wygody.

=== "Po"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Zbierz wyjścia HaplotypeCaller ze wszystkich próbek

Gdybyśmy podłączyli kanał wyjściowy z `GATK_HAPLOTYPECALLER` bezpośrednio do nowego procesu, Nextflow wywołałby proces na każdym pliku GVCF próbki osobno.
Chcemy połączyć wszystkie trzy pliki GVCF (i ich pliki indeksów) tak, aby Nextflow przekazał je wszystkie razem do pojedynczego wywołania procesu.

Możemy to zrobić używając operatora kanału `collect()`.
Dodaj następujące linie do treści `workflow`, zaraz po wywołaniu GATK_HAPLOTYPECALLER:

=== "Po"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Przed"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Rozbijając to na czynniki:

1. Pobieramy kanał wyjściowy z `GATK_HAPLOTYPECALLER` używając właściwości `.out`.
2. Ponieważ nadaliśmy nazwy wyjściom używając `emit:` w sekcji 1, możemy wybrać pliki GVCF za pomocą `.vcf`, a pliki indeksów za pomocą `.idx`. Bez nazwanych wyjść musielibyśmy użyć `.out[0]` i `.out[1]`.
3. Operator `collect()` łączy wszystkie pliki w jeden element, więc `all_gvcfs_ch` zawiera wszystkie trzy pliki GVCF razem, a `all_idxs_ch` zawiera wszystkie trzy pliki indeksów razem.

Możemy zebrać pliki GVCF i ich pliki indeksów osobno (w przeciwieństwie do trzymania ich razem w krotkach), ponieważ Nextflow umieści wszystkie pliki wejściowe razem do wykonania, więc pliki indeksów będą obecne obok plików GVCF.

!!! tip "Wskazówka"

    Możesz użyć operatora `view()`, aby sprawdzić zawartość kanałów przed i po zastosowaniu operatorów kanału.

### 2.2. Napisz proces wspólnego genotypowania i wywołaj go w workflow'ie

Podążając za tym samym wzorcem, którego użyliśmy w Części 2, napiszemy definicję procesu w pliku modułu, zaimportujemy go do workflow'a i wywołamy go na wejściach, które właśnie przygotowaliśmy.

#### 2.2.1. Skonstruuj ciąg znaków, aby nadać każdemu GVCF argument `-V`

Zanim zaczniemy wypełniać definicję procesu, jest jedna rzecz do rozwiązania.
Polecenie GenomicsDBImport oczekuje osobnego argumentu `-V` dla każdego pliku GVCF, w ten sposób:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Gdybyśmy napisali `-V ${all_gvcfs_ch}`, Nextflow po prostu połączyłby nazwy plików i ta część polecenia wyglądałaby tak:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Ale potrzebujemy, aby ciąg wyglądał tak:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Co ważne, musimy skonstruować ten ciąg dynamicznie z jakichkolwiek plików znajdujących się w zebranym kanale.
Nextflow (poprzez Groovy) zapewnia zwięzły sposób na zrobienie tego:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Rozbijając to na czynniki:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` iteruje po każdej ścieżce pliku i dodaje przed nią `-V `, generując `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` łączy je ze spacjami: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. Wynik jest przypisany do zmiennej lokalnej `gvcfs_line` (zdefiniowanej za pomocą `def`), którą możemy interpolować w szablonie polecenia.

Ta linia trafia wewnątrz bloku `script:` procesu, przed szablonem polecenia.
Możesz umieścić dowolny kod Groovy między `script:` a otwierającym `"""` szablonu polecenia.

Następnie będziesz mógł/-mogła odwoływać się do całego tego ciągu jako `gvcfs_line` w bloku `script:` procesu.

#### 2.2.2. Wypełnij moduł dla procesu wspólnego genotypowania

Teraz możemy zająć się napisaniem pełnego procesu.

Otwórz `modules/gatk_jointgenotyping.nf` i zbadaj zarys definicji procesu.

Wypełnij definicję procesu używając informacji podanych powyżej, a następnie sprawdź swoją pracę z rozwiązaniem w zakładce "Po" poniżej.

=== "Przed"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Połącz pliki GVCF w repozytorium danych GenomicsDB i uruchom wspólne genotypowanie, aby wygenerować wykryte warianty na poziomie kohorty
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

=== "Po"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Połącz pliki GVCF w repozytorium danych GenomicsDB i uruchom wspólne genotypowanie, aby wygenerować wykryte warianty na poziomie kohorty
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

Warto zwrócić uwagę na kilka rzeczy.

Jak wcześniej, kilka wejść jest wymienionych, mimo że polecenia nie odwołują się do nich bezpośrednio: `all_idxs`, `ref_index` i `ref_dict`.
Wymienienie ich zapewnia, że Nextflow umieści te pliki w katalogu roboczym obok plików, które pojawiają się w poleceniach, których GATK oczekuje znaleźć na podstawie konwencji nazewnictwa.

Zmienna `gvcfs_line` używa domknięcia Groovy opisanego powyżej do skonstruowania argumentów `-V` dla GenomicsDBImport.

Ten proces uruchamia dwa polecenia szeregowo, tak jak zrobiłbyś/-aś to w terminalu.
GenomicsDBImport łączy pliki GVCF z pojedynczych próbek w repozytorium danych, a następnie GenotypeGVCFs odczytuje to repozytorium danych i generuje ostateczny plik VCF na poziomie kohorty.
Repozytorium danych GenomicsDB (`${cohort_name}_gdb`) jest artefaktem pośrednim używanym tylko wewnątrz procesu; nie pojawia się w bloku wyjścia.

Po ukończeniu tego proces jest gotowy do użycia.
Aby go użyć w workflow'ie, musisz zaimportować moduł i dodać wywołanie procesu.

#### 2.2.3. Zaimportuj moduł

Dodaj instrukcję importu do `genomics.nf`, poniżej istniejących instrukcji importu:

=== "Po"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Przed"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

Proces jest teraz dostępny w zakresie workflow'a.

#### 2.2.4. Dodaj wywołanie procesu

Dodaj wywołanie `GATK_JOINTGENOTYPING` w treści workflow'a, po liniach `collect()`:

=== "Po"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
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

=== "Przed"

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

Proces jest teraz w pełni podłączony.
Następnie skonfigurujemy sposób publikowania wyjść.

### 2.3. Skonfiguruj obsługę wyjścia

Musimy opublikować wyjścia wspólnego pliku VCF.
Dodaj cele publikowania i wpisy bloku wyjścia dla wyników wspólnego genotypowania.

#### 2.3.1. Dodaj cele publikowania dla wspólnego VCF

Dodaj wspólny VCF i jego indeks do sekcji `publish:` workflow'a:

=== "Po"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Przed"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Teraz zaktualizuj blok wyjścia, aby odpowiadał.

#### 2.3.2. Dodaj wpisy bloku wyjścia dla wspólnego VCF

Dodaj wpisy dla plików wspólnego VCF.
Umieścimy je w katalogu głównym katalogu wyników, ponieważ jest to ostateczne wyjście.

=== "Po"

    ```groovy title="genomics.nf" hl_lines="11-16"
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

=== "Przed"

    ```groovy title="genomics.nf"
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

Z procesem, celami publikowania i blokiem wyjścia wszystko jest na miejscu. Możemy teraz przetestować kompletny workflow.

### 2.4. Uruchom workflow'a

Uruchom workflow'a, aby sprawdzić, czy wszystko działa.

```bash
nextflow run genomics.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Pierwsze dwa kroki są zbuforowane z poprzedniego uruchomienia, a nowy krok `GATK_JOINTGENOTYPING` wykonuje się raz na zebranych wejściach ze wszystkich trzech próbek.
Ostateczny plik wyjściowy, `family_trio.joint.vcf` (i jego indeks), znajdują się w katalogu wyników.

??? abstract "Zawartość katalogu (skrócone dowiązania symboliczne)"

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

Jeśli otworzysz wspólny plik VCF, możesz sprawdzić, że workflow wygenerował oczekiwane wykryte warianty.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Masz teraz zautomatyzowany, całkowicie odtwarzalny workflow wspólnego wyszukiwania wariantów!

!!! note "Uwaga"

    Pamiętaj, że pliki danych, które Ci przekazaliśmy, obejmują tylko malutką część chromosomu 20.
    Rzeczywisty rozmiar zestawu wykrytych wariantów byłby liczony w milionach wariantów.
    Dlatego używamy tylko małych podzbiorów danych do celów szkoleniowych!

### Podsumowanie

Wiesz, jak zbierać wyjścia z kanału i łączyć je jako pojedyncze wejście do innego procesu.
Wiesz również, jak konstruować linię poleceń używając domknięć Groovy i jak uruchamiać wiele poleceń w pojedynczym procesie.

### Co dalej?

Pogratuluj sobie! Ukończyłeś/-aś kurs Nextflow dla Genomiki.

Przejdź do końcowego [podsumowania kursu](./next_steps.md), aby przypomnieć sobie, czego się nauczyłeś/-aś i dowiedzieć się, co dalej.
