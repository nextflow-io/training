# Część 1: Wykrywanie wariantów dla pojedynczych próbek

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W pierwszej części tego kursu pokażemy, jak zbudować prosty pipeline do wykrywania wariantów, który stosuje narzędzie GATK do wykrywania wariantów w indywidualnych próbkach sekwencjonowania.

### Przegląd metody

Wykrywanie wariantów to metoda analizy genomowej, która ma na celu identyfikację zmienności w sekwencji genomu względem genomu referencyjnego.
Tutaj użyjemy narzędzi i metod zaprojektowanych do wykrywania krótkich wariantów, _tzn._ SNP i indeli.

![GATK pipeline](img/gatk-pipeline.png)

Pełny pipeline do wykrywania wariantów zazwyczaj obejmuje wiele kroków, w tym mapowanie do referencji (czasami nazywane wyrównaniem do genomu) oraz filtrowanie i priorytetyzację wariantów.
Dla uproszczenia, w tej części kursu skupimy się tylko na części wykrywania wariantów.

### Zestaw danych

Udostępniamy następujące dane i powiązane zasoby:

- **Genom referencyjny** składający się z niewielkiego regionu ludzkiego chromosomu 20 (z hg19/b37) i jego plików pomocniczych (indeks i słownik sekwencji).
- **Trzy próbki sekwencjonowania całego genomu** odpowiadające trio rodzinnemu (matka, ojciec i syn), które zostały ograniczone do małego fragmentu danych na chromosomie 20, aby zachować małe rozmiary plików.
  To dane z sekwencjonowania Illumina krótkich odczytów, które zostały już zmapowane do genomu referencyjnego, dostarczone w formacie [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Binary Alignment Map, skompresowana wersja SAM, Sequence Alignment Map).
- **Lista przedziałów genomowych**, _tzn._ współrzędnych na genomie, gdzie nasze próbki mają dane odpowiednie do wykrywania wariantów, dostarczona w formacie BED.

### Workflow

W tej części kursu opracujesz workflow, który wykonuje następujące czynności:

1. Generuje plik indeksu dla każdego pliku wejściowego BAM używając [Samtools](https://www.htslib.org/)
2. Uruchamia GATK HaplotypeCaller na każdym pliku wejściowym BAM, aby wygenerować wykrycia wariantów dla pojedynczych próbek w formacie VCF (Variant Call Format)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-1.svg"
</figure>

!!! note "Uwaga"

    Pliki indeksów są powszechną cechą formatów plików bioinformatycznych; zawierają informacje o strukturze pliku głównego, które pozwalają narzędziom takim jak GATK na dostęp do podzbioru danych bez konieczności przeczytania całego pliku.
    Jest to ważne ze względu na to, jak duże mogą być te pliki.

---

## 0. Rozgrzewka: Przetestuj polecenia Samtools i GATK interaktywnie

Najpierw chcemy wypróbować polecenia ręcznie, zanim spróbujemy umieścić je w workflow'ie.
Narzędzia, których potrzebujemy (Samtools i GATK) nie są zainstalowane w środowisku GitHub Codespaces, więc użyjemy ich przez kontenery (zobacz [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Uwaga"

     Upewnij się, że jesteś w katalogu `nf4-science/genomics`, aby ostatnia część ścieżki pokazana po wpisaniu `pwd` była `genomics`.

### 0.1. Indeksuj plik wejściowy BAM używając Samtools

Pobierzemy kontener Samtools, uruchomimy go interaktywnie i wykonamy polecenie `samtools index` na jednym z plików BAM.

#### 0.1.1. Pobierz kontener Samtools

```bash
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

#### 0.1.2. Uruchom kontener Samtools interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

#### 0.1.3. Uruchom polecenie indeksowania

[Dokumentacja Samtools](https://www.htslib.org/doc/samtools-index.html) podaje nam linię poleceń do uruchomienia w celu indeksowania pliku BAM.

Musimy tylko podać plik wejściowy; narzędzie automatycznie wygeneruje nazwę dla pliku wyjściowego, dodając `.bai` do nazwy pliku wejściowego.

```bash
samtools index /data/bam/reads_mother.bam
```

To powinno zakończyć się natychmiast, a teraz powinieneś zobaczyć plik o nazwie `reads_mother.bam.bai` w tym samym katalogu co oryginalny plik wejściowy BAM.

??? abstract "Zawartość katalogu"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    └── reads_son.bam
    ```

#### 0.1.4. Wyjdź z kontenera Samtools

```bash
exit
```

### 0.2. Wykryj warianty używając GATK HaplotypeCaller

Pobierzemy kontener GATK, uruchomimy go interaktywnie i wykonamy polecenie `gatk HaplotypeCaller` na pliku BAM, który właśnie zaindeksowaliśmy.

#### 0.2.1. Pobierz kontener GATK

```bash
docker pull community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

#### 0.2.2. Uruchom kontener GATK interaktywnie

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

#### 0.2.3. Uruchom polecenie wykrywania wariantów

[Dokumentacja GATK](https://gatk.broadinstitute.org/hc/en-us/articles/21905025322523-HaplotypeCaller) podaje nam linię poleceń do uruchomienia w celu wykrycia wariantów w pliku BAM.

Musimy podać plik wejściowy BAM (`-I`) oraz genom referencyjny (`-R`), nazwę dla pliku wyjściowego (`-O`) i listę przedziałów genomowych do analizy (`-L`).

Nie musimy jednak podawać ścieżki do pliku indeksu; narzędzie automatycznie będzie go szukać w tym samym katalogu, w oparciu o ustaloną konwencję nazewnictwa i współlokalizacji.
To samo dotyczy plików pomocniczych genomu referencyjnego (pliki indeksu i słownika sekwencji, `*.fai` i `*.dict`).

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.vcf \
        -L /data/ref/intervals.bed
```

<!--
??? success "Wyjście polecenia"

    ```console

    ```
-->

Plik wyjściowy `reads_mother.vcf` jest tworzony wewnątrz Twojego katalogu roboczego w kontenerze, więc nie zobaczysz go w eksploratorze plików VS Code, chyba że zmienisz ścieżkę pliku wyjściowego.
Jednak jest to mały plik testowy, więc możesz użyć `cat`, aby otworzyć i zobaczyć jego zawartość.
Jeśli przejdziesz na sam początek pliku, znajdziesz nagłówek składający się z wielu linii metadanych, po których następuje lista wykrytych wariantów, po jednym w każdej linii.

```console title="reads_mother.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

Każda linia opisuje możliwy wariant zidentyfikowany w danych sekwencjonowania próbki. Aby uzyskać wskazówki dotyczące interpretacji formatu VCF, zobacz [ten pomocny artykuł](https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/).

Plikowi wyjściowemu VCF towarzyszy plik indeksu o nazwie `reads_mother.vcf.idx`, który został automatycznie utworzony przez GATK.
Ma taką samą funkcję jak plik indeksu BAM, aby umożliwić narzędziom wyszukiwanie i pobieranie podzbiorów danych bez ładowania całego pliku.

#### 0.2.4. Wyjdź z kontenera GATK

```bash
exit
```

### Podsumowanie

Wiesz już, jak przetestować polecenia indeksowania Samtools i wykrywania wariantów GATK w ich odpowiednich kontenerach.

### Co dalej?

Naucz się, jak opakować te same polecenia w dwuetapowy workflow, który używa kontenerów do wykonywania pracy.

---

## 1. Napisz jednoetapowy workflow, który uruchamia Samtools index na pliku BAM

Udostępniamy plik workflow'u, `genomics-1.nf`, który przedstawia główne części workflow'u.
Nie jest funkcjonalny; jego celem jest tylko służenie jako szkielet, którego użyjesz do napisania rzeczywistego workflow'u.

### 1.1. Zdefiniuj proces indeksowania

Zacznijmy od napisania procesu, który nazwiemy `SAMTOOLS_INDEX`, opisującego operację indeksowania.

```groovy title="genomics-1.nf" linenums="9"
/*
 * Generuj plik indeksu BAM
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

Rozpoznasz wszystkie elementy z tego, czego nauczyłeś się w Części 1 i Części 2 tej serii szkoleń.

Ten proces będzie wymagał od nas przekazania ścieżki pliku przez wejście `input_bam`, więc skonfiguruj to następnie.

### 1.2. Dodaj deklarację parametru wejściowego

Na górze pliku, w sekcji `Pipeline parameters`, deklarujemy parametr CLI o nazwie `reads_bam` i nadajemy mu wartość domyślną.
W ten sposób możesz być leniwy i nie podawać wejścia podczas wpisywania polecenia uruchamiającego pipeline (w celach deweloperskich).

```groovy title="genomics-1.nf" linenums="3"
/*
 * Parametry pipeline
 */
params {
    // Główne wejście
    reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
}
```

Teraz masz gotowy proces, a także parametr do podania mu wejścia do przetworzenia, więc połącz te rzeczy razem.

!!! note "Uwaga"

    `${projectDir}` to wbudowana zmienna Nextflow, która wskazuje na katalog, w którym znajduje się obecny skrypt workflow Nextflow (`genomics-1.nf`).

    To ułatwia odwoływanie się do plików, katalogów danych i innych zasobów zawartych w repozytorium workflow'u bez kodowania ścieżek bezwzględnych.

### 1.3. Dodaj blok workflow, aby uruchomić SAMTOOLS_INDEX

W bloku `workflow` musisz skonfigurować **kanał**, aby przekazać wejście do procesu `SAMTOOLS_INDEX`; następnie możesz wywołać sam proces, aby uruchomić go na zawartości tego kanału.

```groovy title="genomics-1.nf" linenums="24"
workflow {

    main:
    // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
    reads_ch = channel.fromPath(params.reads_bam)

    // Utwórz plik indeksu dla wejściowego pliku BAM
    SAMTOOLS_INDEX(reads_ch)

    publish:
    bam_index = SAMTOOLS_INDEX.out
}
```

Blok workflow'u ma dwie sekcje:

- `main:` zawiera operacje na kanałach i wywołania procesów
- `publish:` deklaruje, które wyjścia powinny być opublikowane, przypisując je do nazwanych celów

Zauważ, że używamy tej samej fabryki kanałów `.fromPath`, której używaliśmy w [Hello Channels](../../hello_nextflow/02_hello_channels.md).
Rzeczywiście, robimy coś bardzo podobnego.
Różnica polega na tym, że mówimy Nextflow, aby po prostu załadował samą ścieżkę pliku do kanału jako element wejściowy, zamiast czytać jego zawartość.

### 1.4. Dodaj blok output, aby zdefiniować, gdzie publikowane są wyniki

Po bloku workflow'u dodajemy blok `output`, który określa, gdzie publikować wyjścia workflow'u.

```groovy title="genomics-1.nf" linenums="37"
output {
    bam_index {
        path '.'
    }
}
```

Każdy nazwany cel z sekcji `publish:` (jak `bam_index`) otrzymuje Swój własny blok, w którym można skonfigurować ścieżkę wyjściową względem bazowego katalogu wyjściowego.

!!! note "Uwaga"

    Mimo że pliki danych, których używamy tutaj, są bardzo małe, w genomice mogą być bardzo duże.
    Domyślnie Nextflow tworzy dowiązania symboliczne do plików wyjściowych w katalogu publikacji, co pozwala uniknąć niepotrzebnych kopii plików.
    Możesz zmienić to zachowanie używając opcji `mode` (np. `mode 'copy'`), aby utworzyć rzeczywiste kopie.
    Należy pamiętać, że dowiązania symboliczne przestaną działać po wyczyszczeniu katalogu `work`, więc dla produkcyjnych workflow'ów możesz chcieć użyć `mode 'copy'`.

### 1.5. Skonfiguruj katalog wyjściowy

Bazowy katalog wyjściowy jest ustawiany przez opcję konfiguracyjną `outputDir`. Dodaj go do `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" hl_lines="2"
    docker.enabled = true
    outputDir = 'results_genomics'
    ```

=== "Przed"

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

### 1.6. Uruchom workflow, aby zweryfikować, że krok indeksowania działa

Uruchom workflow! Przypominamy, że nie musisz podawać wejścia w linii poleceń, ponieważ ustawiłeś wartość domyślną dla wejścia podczas deklarowania parametru wejściowego.

```bash
nextflow run genomics-1.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [reverent_sinoussi] DSL2 - revision: 41d43ad7fe

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1) | 1 of 1 ✔
    ```

Możesz sprawdzić, czy plik indeksu został wygenerowany poprawnie, przeglądając katalog roboczy lub katalog wyników.

??? abstract "Zawartość katalogu roboczego"

    ```console
    work/2a/e695367b2f60df09cf826b07192dc3
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    └── reads_mother.bam.bai
    ```

??? abstract "Zawartość katalogu wyników"

    ```console
    results_genomics/
    └── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ```

Oto jest!

### Podsumowanie

Wiesz już, jak opakować narzędzie genomiczne w jednoetapowy workflow Nextflow i uruchomić je używając kontenera.

### Co dalej?

Dodaj drugi krok, który wykorzystuje wyjście pierwszego.

---

## 2. Dodaj drugi proces, aby uruchomić GATK HaplotypeCaller na zaindeksowanym pliku BAM

Teraz, gdy masz indeks dla pliku wejściowego, możesz przejść do skonfigurowania kroku wykrywania wariantów, który jest interesującą częścią workflow'u.

### 2.1. Zdefiniuj proces wykrywania wariantów

Napiszmy proces, który nazwiemy `GATK_HAPLOTYPECALLER`, opisujący operację wykrywania wariantów.

```groovy title="genomics-1.nf" linenums="44"
/*
 * Wykryj warianty używając GATK HaplotypeCaller
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

Zauważ, że wprowadziliśmy tutaj nową składnię (`emit:`), aby jednoznacznie nazwać każdy z naszych kanałów wyjściowych, a powody tego staną się wkrótce jasne.

To polecenie przyjmuje znacznie więcej wejść, ponieważ GATK potrzebuje więcej informacji do wykonania analizy w porównaniu do prostego zadania indeksowania.
Ale zauważ, że jest jeszcze więcej wejść zdefiniowanych w bloku wejść niż jest wymienionych w poleceniu GATK. Dlaczego?

!!! note "Uwaga"

    GATK wie, gdzie szukać pliku indeksu BAM i plików pomocniczych genomu referencyjnego, ponieważ zna konwencje związane z tymi plikami.
    Jednak Nextflow jest zaprojektowany jako niezależny od dziedziny i nie wie nic o wymaganiach formatów plików bioinformatycznych.

Musisz powiedzieć Nextflow wyraźnie, że musi umieścić te pliki w katalogu roboczym w czasie wykonywania; w przeciwnym razie tego nie zrobi, a GATK (prawidłowo) zgłosi błąd dotyczący brakujących plików indeksów.

Podobnie musimy wyraźnie wymienić plik indeksu wyjściowego VCF (plik `"${input_bam}.vcf.idx"`), aby Nextflow wiedział, że ma śledzić ten plik na wypadek, gdyby był potrzebny w kolejnych krokach.

### 2.2. Dodaj definicje dla wejść pomocniczych

Ponieważ Twój nowy proces oczekuje kilku dodatkowych plików, ustaw dla nich parametry CLI w sekcji `Pipeline parameters`, wraz z wartościami domyślnymi (z tych samych powodów co wcześniej).

```groovy title="genomics-1.nf" linenums="8"
    // Pliki pomocnicze
    reference: Path = "${projectDir}/data/ref/ref.fasta"
    reference_index: Path = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict: Path = "${projectDir}/data/ref/ref.dict"
    intervals: Path = "${projectDir}/data/ref/intervals.bed"
```

### 2.3. Utwórz zmienne do przechowywania ścieżek plików pomocniczych

Podczas gdy główne wejścia danych są przesyłane dynamicznie przez kanały, istnieją dwa podejścia do obsługi plików pomocniczych. Zalecanym podejściem jest tworzenie jawnych kanałów, co sprawia, że przepływ danych jest bardziej przejrzysty i spójny. Alternatywnie, funkcja file() może być użyta do tworzenia zmiennych w prostszych przypadkach, szczególnie gdy musisz odwołać się do tego samego pliku w wielu procesach - chociaż pamiętaj, że nadal tworzy to kanały niejawnie. <!-- TODO: Wyjaśnić: czy to jest nadal konieczne z typowanymi wejściami? -->

Dodaj to do bloku workflow'u (po utworzeniu `reads_ch`, wewnątrz sekcji `main:`):

```groovy title="genomics-1.nf" linenums="79"
    // Załaduj ścieżki plików dla plików pomocniczych (referencja i przedziały)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)
```

To sprawi, że ścieżki plików pomocniczych będą dostępne do dostarczenia jako wejście do wszystkich procesów, które ich potrzebują.

### 2.4. Dodaj wywołanie do bloku workflow, aby uruchomić GATK_HAPLOTYPECALLER

Teraz, gdy masz skonfigurowany drugi proces i wszystkie wejścia oraz pliki pomocnicze są gotowe i dostępne, możesz dodać wywołanie procesu `GATK_HAPLOTYPECALLER` w ciele workflow'u.

```groovy title="genomics-1.nf" linenums="88"
    // Wykryj warianty z zaindeksowanego pliku BAM
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
```

Rozpoznasz składnię `*.out` z Części 1 tej serii szkoleń; mówimy Nextflow, aby wziął kanał wyjściowy z `SAMTOOLS_INDEX` i podłączył go do wywołania procesu `GATK_HAPLOTYPECALLER`.

!!! note "Uwaga"

    Zauważ, że wejścia są podawane w dokładnie tej samej kolejności w wywołaniu procesu, jak są wymienione w bloku wejść procesu.
    W Nextflow wejścia są pozycyjne, co oznacza, że _musisz_ zachować tę samą kolejność; i oczywiście musi być taka sama liczba elementów.

### 2.5. Zaktualizuj sekcję publish i blok output

Musisz zaktualizować sekcję `publish:`, aby uwzględnić wyjścia VCF, i dodać odpowiadające cele w bloku `output`.

```groovy title="genomics-1.nf" linenums="99"
    publish:
    bam_index = SAMTOOLS_INDEX.out
    vcf = GATK_HAPLOTYPECALLER.out.vcf
    vcf_idx = GATK_HAPLOTYPECALLER.out.idx
}

output {
    bam_index {
        path '.'
    }
    vcf {
        path '.'
    }
    vcf_idx {
        path '.'
    }
}
```

### 2.6. Uruchom workflow, aby zweryfikować, że krok wykrywania wariantów działa

Uruchom rozszerzony workflow z `-resume`, aby nie musieć ponownie uruchamiać kroku indeksowania.

```bash
nextflow run genomics-1.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [grave_volta] DSL2 - revision: 4790abc96a

    executor >  local (1)
    [2a/e69536] SAMTOOLS_INDEX (1)       | 1 of 1, cached: 1 ✔
    [53/e18e98] GATK_HAPLOTYPECALLER (1) | 1 of 1 ✔
    ```

Teraz, jeśli spojrzysz na wyjście konsoli, zobaczysz wymienione dwa procesy.

Pierwszy proces został pominięty dzięki cache'owaniu, zgodnie z oczekiwaniami, podczas gdy drugi proces został uruchomiony, ponieważ jest zupełnie nowy.

Pliki wyjściowe znajdziesz w katalogu wyników (jako dowiązania symboliczne do katalogu roboczego).

??? abstract "Zawartość katalogu"

    ```console
    results_genomics/
    ├── reads_mother.bam.bai -> */87/908bba*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */cf/36f756*/reads_mother.bam.vcf
    └── reads_mother.bam.vcf.idx -> */cf/36f756*/reads_mother.bam.vcf.idx
    ```

Jeśli otworzysz plik VCF, zobaczysz tę samą zawartość co w pliku wygenerowanym przez uruchomienie polecenia GATK bezpośrednio w kontenerze.

```console title="reads_mother.bam.vcf" linenums="26"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother
20_10037292_10066351	3480	.	C	CT	503.03	.	AC=2;AF=1.00;AN=2;DP=23;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.95;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,18:18:54:517,54,0
20_10037292_10066351	3520	.	AT	A	609.03	.	AC=2;AF=1.00;AN=2;DP=18;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=33.83;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,18:18:54:623,54,0
20_10037292_10066351	3529	.	T	A	155.64	.	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.544;DP=21;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.78;ReadPosRankSum=-1.158;SOR=1.034	GT:AD:DP:GQ:PL	0/1:12,8:20:99:163,0,328
```

To jest wyjście, które zależy Ci na wygenerowaniu dla każdej próbki w Twoim badaniu.

### Podsumowanie

Wiesz już, jak zrobić bardzo prosty dwuetapowy workflow, który wykonuje prawdziwą pracę analityczną i jest w stanie radzić sobie z osobliwościami formatów plików genomicznych, takimi jak pliki pomocnicze.

### Co dalej?

Spraw, aby workflow obsługiwał wiele próbek naraz.

---

## 3. Dostosuj workflow do pracy z zestawem próbek

Dobrze jest mieć workflow, który może zautomatyzować przetwarzanie pojedynczej próbki, ale co jeśli masz 1000 próbek?
Czy musisz napisać skrypt bash, który przechodzi przez wszystkie Twoje próbki w pętli?

Nie! Wystarczy dokonać drobnej zmiany w kodzie, a Nextflow również to dla Ciebie obsłuży.

### 3.1. Przekształć deklarację parametru wejściowego w tablicę wymieniającą trzy próbki

Przekształć tę domyślną ścieżkę pliku w deklaracji wejściowego pliku BAM w tablicę wymieniającą ścieżki plików dla trzech próbek testowych, w górze w sekcji `Pipeline parameters`.

=== "Po"

    ```groovy title="genomics-1.nf" linenums="7"
    // Główne wejście (tablica trzech próbek)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

=== "Przed"

    ```groovy title="genomics-1.nf" linenums="7"
        // Główne wejście
        reads_bam: Path = "${projectDir}/data/bam/reads_mother.bam"
    ```

!!! note "Uwaga"

    Używając typowanych deklaracji parametrów (jak `reads_bam: Path`), nie można przypisać wartości tablicowej.
    Dla tablic pomiń adnotację typu.

I to właściwie wszystko, co musimy zrobić, ponieważ fabryka kanałów, której używamy w ciele workflow'u (`.fromPath`), jest równie chętna do zaakceptowania wielu ścieżek plików do załadowania do kanału wejściowego, jak była do załadowania jednej.

!!! note "Uwaga"

    Normalnie nie chciałbyś umieszczać listy próbek bezpośrednio w pliku workflow'u, ale robimy to tutaj, aby zachować prostotę.
    Przedstawimy bardziej eleganckie sposoby obsługi wejść później w tej serii szkoleń.

### 3.2. Uruchom workflow, aby zweryfikować, że działa na wszystkich trzech próbkach

Spróbuj teraz uruchomić workflow, gdy połączenia są skonfigurowane do pracy na wszystkich trzech próbkach testowych.

```bash
nextflow run genomics-1.nf -resume
```

Ciekawa rzecz: to _może zadziałać_, LUB _może zawieść_. Na przykład, oto uruchomienie, które zakończyło się sukcesem:

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [peaceful_yalow] DSL2 - revision: a256d113ad

    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 1 ✔
    [7a/89bc43] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 1 ✔
    ```

Jeśli Twoje uruchomienie workflow'u zakończyło się sukcesem, uruchom je ponownie, aż otrzymasz błąd taki jak ten:

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [loving_pasteur] DSL2 - revision: d2a8e63076

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

Jeśli spojrzysz na wyjście błędu polecenia GATK, będzie tam linia taka jak ta:

```console
A USER ERROR has occurred: Traversal by intervals was requested but some input files are not indexed.
```

Cóż, to dziwne, biorąc pod uwagę, że wyraźnie zaindeksowałeś pliki BAM w pierwszym kroku workflow'u. Czy może być coś nie tak z połączeniami?

#### 3.2.1. Sprawdź katalogi robocze dla odpowiednich wywołań

Spójrz do środka katalogu roboczego dla nieudanego wywołania procesu `GATK_HAPLOTYPECALLER` wymienionego w wyjściu konsoli.

??? abstract "Zawartość katalogu"

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

Zwróć szczególną uwagę na nazwy pliku BAM i indeksu BAM, które są wymienione w tym katalogu: `reads_son.bam` i `reads_father.bam.bai`.

Co do cholery? Nextflow umieścił plik indeksu w katalogu roboczym tego wywołania procesu, ale jest to zły. Jak to mogło się stać?

#### 3.2.2. Użyj [operatora view()](https://www.nextflow.io/docs/latest/reference/operator.html#view), aby sprawdzić zawartość kanału

Dodaj te dwie linie w ciele workflow'u przed wywołaniem procesu `GATK_HAPLOTYPER`:

```groovy title="genomics-1.nf" linenums="84"
    // tymczasowa diagnostyka
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()
```

Następnie uruchom ponownie polecenie workflow.

```bash
nextflow run genomics-1.nf
```

Ponownie, to może zakończyć się sukcesem lub porażką. Oto udane uruchomienie:

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [fervent_pasteur] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    [a2/dbd8d5] GATK_HAPLOTYPECALLER (3) | 3 of 3 ✔
    ```

A oto nieudane:

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [angry_hamilton] DSL2 - revision: a256d113ad

    /workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
    /workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
    executor >  local (6)
    [4f/7071b0] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    /workspaces/training/nf4-science/genomics/work/4f/7071b082b45dd85b1c9b6b3b32cb69/reads_father.bam.bai
    /workspaces/training/nf4-science/genomics/work/3c/331645a9e20e67edae10da5ba17c7b/reads_son.bam.bai
    /workspaces/training/nf4-science/genomics/work/b4/45a376f0e724be1dc626a6807f73d8/reads_mother.bam.bai
    [a3/cf3a89] GATK_HAPLOTYPECALLER (3) | 1 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (3)'
    ...
    ```

Może być konieczne uruchomienie go kilka razy, aby ponownie wystąpił błąd.
Ten błąd nie będzie się reprodukował konsekwentnie, ponieważ jest zależny od pewnej zmienności w czasach wykonywania poszczególnych wywołań procesów.

Oto jak wygląda wyjście dwóch wywołań `.view()`, które dodaliśmy, dla nieudanego uruchomienia:

```console
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
/workspaces/training/nf4-science/genomics/work/9c/53492e3518447b75363e1cd951be4b/reads_father.bam.bai
/workspaces/training/nf4-science/genomics/work/cc/37894fffdf6cc84c3b0b47f9b536b7/reads_son.bam.bai
/workspaces/training/nf4-science/genomics/work/4d/dff681a3d137ba7d9866e3d9307bd0/reads_mother.bam.bai
```

Pierwsze trzy linie odpowiadają kanałowi wejściowemu, a drugie - kanałowi wyjściowemu.
Widać, że pliki BAM i pliki indeksów dla trzech próbek nie są wymienione w tej samej kolejności!

!!! note "Uwaga"

    Gdy wywołujesz proces na kanale zawierającym wiele elementów, Nextflow będzie próbował sparalelizować wykonywanie tak bardzo, jak to możliwe.
    Wyjścia są zbierane w kolejności, w jakiej staną się dostępne — może to być inna kolejność niż pierwotna kolejność wejść.

Zgodnie z obecnym zapisem, nasz skrypt workflow'u zakłada, że pliki indeksów wyjdą z kroku indeksowania wymienione w tej samej kolejności matka/ojciec/syn, jak podano wejścia.
Ale nie jest to gwarantowane, dlatego czasami (choć nie zawsze) złe pliki zostają sparowane w drugim kroku.

Aby to naprawić, musisz upewnić się, że pliki BAM i ich pliki indeksów przemieszczają się razem przez kanały.

!!! tip "Wskazówka"

    Instrukcje `view()` w kodzie workflow'u nic nie robią, więc nie ma problemu, aby je zostawić.
    Jednak będą zaśmiecać Twoje wyjście konsoli, więc zalecamy ich usunięcie, gdy zakończysz rozwiązywanie problemu.

### 3.3. Zmień wyjście procesu SAMTOOLS_INDEX na krotkę, która trzyma plik wejściowy i jego indeks razem

Najprostszym sposobem, aby upewnić się, że plik BAM i jego indeks pozostają ściśle powiązane, jest spakowanie ich razem w krotkę wychodzącą z zadania indeksowania.

!!! note "Uwaga"

    **Krotka** to skończona, uporządkowana lista elementów, która jest powszechnie używana do zwracania wielu wartości z funkcji. Krotki są szczególnie przydatne do przekazywania wielu wejść lub wyjść między procesami przy zachowaniu ich powiązania i kolejności.

Najpierw zmień wyjście procesu `SAMTOOLS_INDEX`, aby uwzględnić plik BAM w Swojej deklaracji wyjściowej.

=== "Po"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        tuple path(input_bam), path("${input_bam}.bai")
    ```

=== "Przed"

    ```groovy title="genomics-1.nf" linenums="20"
        output:
        path "${input_bam}.bai"
    ```

W ten sposób każdy plik indeksu będzie ściśle połączony z jego oryginalnym plikiem BAM, a ogólne wyjście kroku indeksowania będzie pojedynczym kanałem zawierającym pary plików.

### 3.4. Zmień wejście do procesu GATK_HAPLOTYPECALLER na krotkę

Ponieważ zmieniłeś "kształt" wyjścia pierwszego procesu w workflow'u, musisz zaktualizować definicję wejścia drugiego procesu, aby pasowała.

Konkretnie, tam gdzie wcześniej deklarowaliśmy dwie oddzielne ścieżki wejściowe w bloku wejść procesu `GATK_HAPLOTYPECALLER`, teraz deklarujemy pojedyncze wejście pasujące do struktury krotki emitowanej przez `SAMTOOLS_INDEX`.

=== "Po"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        tuple path(input_bam), path(input_bam_index)
    ```

=== "Przed"

    ```groovy title="genomics-1.nf" linenums="51"
        input:
        path input_bam
        path input_bam_index
    ```

Ponieważ zmieniłeś teraz kształt wejść, których oczekuje `GATK_HAPLOTYPECALLER`, musisz odpowiednio zaktualizować wywołanie procesu w ciele workflow'u.

### 3.5. Zaktualizuj wywołanie GATK_HAPLOTYPECALLER w bloku workflow

Nie musisz już dostarczać oryginalnego `reads_ch` do procesu `GATK_HAPLOTYPECALLER`, ponieważ plik BAM jest teraz spakowany w kanał wyjściowy przez `SAMTOOLS_INDEX`.

W rezultacie możesz po prostu usunąć tę linię.

=== "Po"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
    ```

=== "Przed"

    ```groovy title="genomics-1.nf" linenums="88"
    GATK_HAPLOTYPECALLER(
        reads_ch,
        SAMTOOLS_INDEX.out,
    ```

To wszystkie zmiany w połączeniach, które są konieczne do rozwiązania problemu niezgodności indeksów.

### 3.6. Zaktualizuj sekcję publish i blok output dla krotki

Ponieważ `SAMTOOLS_INDEX.out` jest teraz krotką zawierającą zarówno BAM, jak i jego indeks, oba pliki będą publikowane razem.
Zmień nazwę celu z `bam_index` na `indexed_bam`, aby odzwierciedlić, że teraz zawiera oba pliki.

=== "Po"

    ```groovy title="genomics-1.nf" hl_lines="2"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
    ```

=== "Przed"

    ```groovy title="genomics-1.nf"
        publish:
        bam_index = SAMTOOLS_INDEX.out
    ```

Musisz również zaktualizować blok output, aby używał nowej nazwy celu:

=== "Po"

    ```groovy title="genomics-1.nf" hl_lines="2"
    output {
        indexed_bam {
            path '.'
        }
    ```

=== "Przed"

    ```groovy title="genomics-1.nf"
    output {
        bam_index {
            path '.'
        }
    ```

### 3.7. Uruchom workflow, aby zweryfikować, że działa poprawnie na wszystkich trzech próbkach za każdym razem

Oczywiście, dowód jest w rezultacie, więc uruchom workflow ponownie kilka razy, aby upewnić się, że będzie to działać niezawodnie w przyszłości.

```bash
nextflow run genomics-1.nf
```

Tym razem (i za każdym razem) wszystko powinno działać poprawnie:

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [special_goldstine] DSL2 - revision: 4cbbf6ea3e

    executor >  local (6)
    [d6/10c2c4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [88/1783aa] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Katalog wyników zawiera teraz zarówno pliki BAM, jak i BAI dla każdej próbki (z krotki), wraz z wyjściami VCF:

??? abstract "Zawartość katalogu wyników"

    ```console
    results_genomics/
    ├── reads_father.bam -> */60/e2614c*/reads_father.bam
    ├── reads_father.bam.bai -> */60/e2614c*/reads_father.bam.bai
    ├── reads_father.bam.vcf -> */b8/91b3c8*/reads_father.bam.vcf
    ├── reads_father.bam.vcf.idx -> */b8/91b3c8*/reads_father.bam.vcf.idx
    ├── reads_mother.bam -> */3e/fededc*/reads_mother.bam
    ├── reads_mother.bam.bai -> */3e/fededc*/reads_mother.bam.bai
    ├── reads_mother.bam.vcf -> */32/5ca037*/reads_mother.bam.vcf
    ├── reads_mother.bam.vcf.idx -> */32/5ca037*/reads_mother.bam.vcf.idx
    ├── reads_son.bam -> */3c/36d1c2*/reads_son.bam
    ├── reads_son.bam.bai -> */3c/36d1c2*/reads_son.bam.bai
    ├── reads_son.bam.vcf -> */d7/a6b046*/reads_son.bam.vcf
    └── reads_son.bam.vcf.idx -> */d7/a6b046*/reads_son.bam.vcf.idx
    ```

Jeśli chcesz, możesz użyć `.view()` ponownie, aby zajrzeć do zawartości kanału wyjściowego `SAMTOOLS_INDEX`:

```groovy title="genomics-1.nf" linenums="92"
SAMTOOLS_INDEX.out.view()
```

Zobaczysz, że kanał zawiera trzy oczekiwane krotki (ścieżki plików skrócone dla czytelności).

```console title="Wyjście"
[*/60/e2614c*/reads_father.bam, */60/e2614c*/reads_father.bam.bai]
[*/3e/fededc*/reads_mother.bam, */3e/fededc*/reads_mother.bam.bai]
[*/3c/36d1c2*/reads_son.bam, */3c/36d1c2*/reads_son.bam.bai]
```

To będzie znacznie bezpieczniejsze w przyszłości.

### Podsumowanie

Wiesz już, jak sprawić, aby Twój workflow działał na wielu próbkach (niezależnie).

### Co dalej?

Ułatw obsługę wielu próbek naraz.

---

## 4. Spraw, aby workflow akceptował plik tekstowy zawierający zestaw plików wejściowych

Bardzo powszechnym sposobem dostarczania wielu danych wejściowych do workflow'u jest użycie pliku tekstowego zawierającego ścieżki.
Może to być prosty plik z jedną ścieżką na linię, lub może zawierać dodatkowe metadane — w takim przypadku często nazywany jest arkuszem próbek.

Tutaj pokażemy, jak zrobić prosty przypadek.

### 4.1. Sprawdź dostarczony plik tekstowy wymieniający ścieżki plików wejściowych

Przygotowaliśmy plik tekstowy wymieniający ścieżki plików wejściowych, nazwany `sample_bams.txt`, który możesz znaleźć w katalogu `data/`.

```txt title="sample_bams.txt"
/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_father.bam
/workspaces/training/nf4-science/genomics/data/bam/reads_son.bam
```

Jak widać, wymieniliśmy jedną ścieżkę pliku na linię, i są to ścieżki bezwzględne.

!!! note "Uwaga"

    Pliki, których tutaj używamy, znajdują się po prostu w lokalnym systemie plików Twojego GitHub Codespaces, ale moglibyśmy również wskazywać na pliki w chmurze.

### 4.2. Zaktualizuj wartość domyślną parametru

Zmień wartość domyślną dla parametru wejściowego `reads_bam`, aby wskazywał na plik `sample_bams.txt`.

=== "Po"

    ```groovy title="genomics-1.nf" linenums="7"
        // Główne wejście (plik plików wejściowych, jeden na linię)
        reads_bam: Path = "${projectDir}/data/sample_bams.txt"
    ```

=== "Przed"

    ```groovy title="genomics-1.nf" linenums="7"
    // Główne wejście (tablica trzech próbek)
        reads_bam = [
            "${projectDir}/data/bam/reads_mother.bam",
            "${projectDir}/data/bam/reads_father.bam",
            "${projectDir}/data/bam/reads_son.bam"
        ]
    ```

W ten sposób możesz nadal być leniwy, ale lista plików nie znajduje się już w samym kodzie workflow'u, co jest dużym krokiem we właściwym kierunku.

### 4.3. Zaktualizuj fabrykę kanałów, aby czytała linie z pliku

Obecnie Twoja fabryka kanałów wejściowych traktuje wszystkie pliki, które jej dajesz, jako dane wejściowe, które chcesz przekazać do procesu indeksowania.
Ponieważ teraz podajesz jej plik, który wymienia ścieżki plików wejściowych, musisz zmienić jej zachowanie, aby parsował plik i traktował ścieżki plików, które zawiera, jako dane wejściowe.

Na szczęście możemy to zrobić bardzo prosto, po prostu dodając [operator `.splitText()`](https://www.nextflow.io/docs/latest/reference/operator.html#operator-splittext) do kroku konstrukcji kanału.

=== "Po"

    ```groovy title="genomics-1.nf" linenums="68"
        // Utwórz kanał wejściowy z pliku tekstowego wymieniającego ścieżki plików wejściowych
        reads_ch = channel.fromPath(params.reads_bam).splitText()
    ```

=== "Przed"

    ```groovy title="genomics-1.nf" linenums="68"
        // Utwórz kanał wejściowy (pojedynczy plik przez parametr CLI)
        reads_ch = channel.fromPath(params.reads_bam)
    ```

!!! tip "Wskazówka"

    To kolejna świetna okazja do użycia operatora `.view()`, aby zobaczyć, jak wygląda zawartość kanału przed i po zastosowaniu operatora.

### 4.4. Uruchom workflow, aby zweryfikować, że działa poprawnie

Uruchom workflow jeszcze raz. Powinno to dać ten sam rezultat co wcześniej, prawda?

```bash
nextflow run genomics-1.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-1.nf` [sick_albattani] DSL2 - revision: 46d84642f6

    [18/23b4bb] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [12/f727bb] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    ```

Tak! W rzeczywistości Nextflow poprawnie wykrywa, że wywołania procesów są dokładnie takie same, i nawet nie zadaje sobie trudu ponownego uruchomienia wszystkiego, ponieważ uruchamialiśmy z `-resume`.

I to wszystko! Twój prosty workflow do wykrywania wariantów ma wszystkie podstawowe funkcje, których potrzebowałeś.

### Podsumowanie

Wiesz, jak stworzyć wieloetapowy liniowy workflow do indeksowania pliku BAM i aplikowania wykrywania wariantów dla pojedynczych próbek używając GATK.

Ogólnie rzecz biorąc, nauczyłeś się używać podstawowych komponentów i logiki Nextflow do budowy prostego pipeline'u genomicznego, który wykonuje prawdziwą pracę, biorąc pod uwagę osobliwości formatów plików genomicznych i wymagania narzędzi.

### Co dalej?

Świętuj Swój sukces i weź ekstra długą przerwę!

W następnej części tego kursu nauczysz się używać kilku dodatkowych funkcji Nextflow (w tym więcej operatorów kanałów), aby zastosować wspólne wykrywanie wariantów do danych.
