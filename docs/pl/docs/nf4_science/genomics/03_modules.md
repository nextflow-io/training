# Część 3: Przenoszenie kodu do modułów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W pierwszej części tego kursu zbudowałeś pipeline do wykrywania wariantów, który był całkowicie liniowy i przetwarzał dane każdej próbki niezależnie od innych.

W drugiej części pokazaliśmy, jak używać kanałów i operatorów kanałów do implementacji wspólnego wykrywania wariantów za pomocą GATK, rozbudowując pipeline z Części 1.

W tej części pokażemy, jak przekonwertować kod z tego workflow'u na moduły. Aby przejść tę część szkolenia, musisz ukończyć Część 1 i Część 2, a także [Hello Modules](../../../hello_nextflow/hello_modules.md), które omawia podstawy modułów.

---

## 0. Rozgrzewka

Kiedy zaczęliśmy rozwijać nasz workflow, umieściliśmy wszystko w jednym pliku kodu.
Teraz nadszedł czas, aby zająć się **modularyzacją** naszego kodu, _tzn._ wyodrębnieniem definicji procesów do modułów.

Zaczniemy od tego samego workflow co w Części 2, który udostępniliśmy w pliku `genomics-3.nf`.

!!! note "Uwaga"

     Upewnij się, że jesteś we właściwym katalogu roboczym:
     `cd /workspaces/training/nf4-science/genomics`

Uruchom workflow, aby zweryfikować punkt wyjścia:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Wyjście"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

W katalogu projektu pojawi się teraz katalog `work` oraz katalog `results_genomics`.

### Podsumowanie

Jesteś gotowy, aby rozpocząć modularyzację Swojego workflow'u.

### Co dalej?

Przenieś procesy workflow'u Genomics do modułów.

---

## 1. Przeniesienie procesów do modułów

Jak nauczyłeś się w [Hello Modules](../../../hello_nextflow/hello_modules.md), tworzenie modułu polega na skopiowaniu definicji procesu do osobnego pliku. Plik może znajdować się w wybranym katalogu i mieć wybraną nazwę.

Z powodów, które staną się jasne później (w szczególności przy testowaniu), w tym szkoleniu stosujemy konwencję nazewnictwa pliku `main.nf` i umieszczania go w strukturze katalogów nazwanej według zestawu narzędzi i polecenia.

### 1.1. Utwórz moduł dla procesu `SAMTOOLS_INDEX`

W przypadku procesu `SAMTOOLS_INDEX`, 'samtools' jest zestawem narzędzi, a 'index' jest poleceniem. Utworzymy więc strukturę katalogów `modules/samtools/index` i umieścimy definicję procesu `SAMTOOLS_INDEX` w pliku `main.nf` wewnątrz tego katalogu.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Otwórz plik `main.nf` i skopiuj do niego definicję procesu `SAMTOOLS_INDEX`.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * Generuje plik indeksu BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Następnie usuń definicję procesu `SAMTOOLS_INDEX` z `genomics-3.nf` i dodaj deklarację importu dla modułu przed definicją następnego procesu, w ten sposób:

=== "Po"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Dołącz moduły
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Wywołaj warianty za pomocą GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Przed"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Wywołaj warianty za pomocą GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

Możesz teraz ponownie uruchomić workflow i powinien działać tak samo jak wcześniej. Jeśli podasz flagę `-resume`, żadne nowe zadania nie powinny być nawet uruchamiane:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Utwórz moduły dla procesów `GATK_HAPLOTYPECALLER` i `GATK_JOINTGENOTYPING`

Powtórz te same kroki dla pozostałych procesów.
Dla każdego procesu:

1. Utwórz strukturę katalogów (`modules/gatk/haplotypecaller/` i `modules/gatk/jointgenotyping/`)
2. Utwórz plik `main.nf` zawierający definicję procesu
3. Usuń definicję procesu z `genomics-3.nf`
4. Dodaj deklarację importu dla modułu

Po zakończeniu sprawdź, czy struktura katalogów modułów jest poprawna, uruchamiając:

```bash
tree modules/
```

??? abstract "Zawartość katalogu"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

Będziesz również mieć coś takiego w głównym pliku workflow'u, po sekcji parametrów:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Podsumowanie

Przećwiczyłeś modularyzację workflow'u na przykładzie workflow'u genomiki.

### Co dalej?

Przetestuj zmodularyzowany workflow.

---

## 2. Testowanie zmodularyzowanego workflow

Uruchom zmodularyzowany workflow, aby zweryfikować, że wszystko nadal działa.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Wyjście"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Wszystko nadal działa, w tym możliwość wznowienia pipeline.
Wyniki są nadal publikowane w katalogu `results_genomics`.

```console title="Zawartość katalogu"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Podsumowanie

Zmodularyzowałeś workflow i zweryfikowałeś, że nadal działa tak samo jak wcześniej.

### Co dalej?

Przejrzyj to, czego się nauczyłeś i spójrz na testowanie.

---

## 3. Podsumowanie

Zmodularyzowałeś workflow i nic nie zmieniło się w sposobie działania pipeline'u.
To jest zamierzone: zrestrukturyzowałeś kod bez wpływu na jego funkcjonalność.

Moduły zawierają tylko logikę procesu, co czyni je czystymi i wielokrotnego użytku.
Główny skrypt kontroluje, co jest publikowane i gdzie, podczas gdy moduły pozostają skoncentrowane na Swoim zadaniu obliczeniowym.

Położyłeś fundamenty pod rzeczy, które ułatwią utrzymanie kodu.
Na przykład możesz teraz dodać testy do Swojego pipeline'u używając frameworka nf-test.
To właśnie przyjrzymy się w następnej części tego kursu.
