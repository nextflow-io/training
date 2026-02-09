# Dzielenie i grupowanie

Nextflow zapewnia potężne narzędzia do elastycznej pracy z danymi. Kluczową możliwością jest dzielenie danych na różne strumienie, a następnie grupowanie powiązanych elementów z powrotem. Jest to szczególnie cenne w workflow'ach bioinformatycznych, gdzie musisz przetwarzać różne typy próbek oddzielnie, a następnie łączyć wyniki do analizy.

Pomyśl o tym jak o sortowaniu poczty: oddzielasz listy według miejsca przeznaczenia, przetwarzasz każdą stertę inaczej, a następnie łączysz elementy trafiające do tej samej osoby. Nextflow używa specjalnych operatorów do osiągnięcia tego z danymi naukowymi. To podejście jest również powszechnie znane jako wzorzec **scatter/gather** w obliczeniach rozproszonych i workflow'ach bioinformatycznych.

System kanałów Nextflow'a jest sercem tej elastyczności. Kanały łączą różne części Twojego workflow'a, umożliwiając przepływ danych przez analizę. Możesz tworzyć wiele kanałów z jednego źródła danych, przetwarzać każdy kanał inaczej, a następnie scalać kanały z powrotem, gdy jest to potrzebne. To podejście pozwala projektować workflow'y, które naturalnie odzwierciedlają rozgałęziające się i zbiegające ścieżki złożonych analiz bioinformatycznych.

### Cele szkolenia

W tej misji pobocznej nauczysz się dzielić i grupować dane przy użyciu operatorów kanałów Nextflow'a.
Zaczniemy od pliku CSV zawierającego informacje o próbkach i powiązanych plikach danych, a następnie będziemy manipulować i reorganizować te dane.

Pod koniec tej misji pobocznej będziesz w stanie efektywnie rozdzielać i łączyć strumienie danych, używając następujących technik:

- Odczytywanie danych z plików przy użyciu `splitCsv`
- Filtrowanie i przekształcanie danych za pomocą `filter` i `map`
- Łączenie powiązanych danych przy użyciu `join` i `groupTuple`
- Tworzenie kombinacji danych za pomocą `combine` do przetwarzania równoległego
- Optymalizacja struktury danych przy użyciu `subMap` i strategii deduplikacji
- Budowanie funkcji wielokrotnego użytku z nazwanymi domknięciami, które pomogą Ci manipulować strukturami kanałów

Te umiejętności pomogą Ci budować workflow'y, które mogą efektywnie obsługiwać wiele plików wejściowych i różne typy danych, zachowując jednocześnie czystą, łatwą w utrzymaniu strukturę kodu.

### Wymagania wstępne

Przed podjęciem tej misji pobocznej powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow'a (procesy, kanały, operatory, praca z plikami, metadane)

**Opcjonalnie:** Zalecamy najpierw ukończenie misji pobocznej [Metadane w workflow'ach](./metadata.md).
Obejmuje ona podstawy odczytu plików CSV za pomocą `splitCsv` i tworzenia map metadanych, których będziemy tutaj intensywnie używać.

---

## 0. Rozpocznij pracę

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzyłeś środowisko szkoleniowe zgodnie z opisem w [Konfiguracji środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki do tego samouczka.

```bash
cd side-quests/splitting_and_grouping
```

Możesz ustawić VSCode, aby skupił się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz główny plik workflow'a i katalog `data` zawierający arkusz próbek o nazwie `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Arkusz próbek zawiera informacje o próbkach od różnych pacjentów, w tym identyfikator pacjenta, numer powtórzenia próbki, typ (normalny lub nowotworowy) oraz ścieżki do hipotetycznych plików danych (które w rzeczywistości nie istnieją, ale będziemy udawać, że istnieją).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Ten arkusz próbek zawiera osiem próbek od trzech pacjentów (A, B, C).

Dla każdego pacjenta mamy próbki typu `tumor` (zazwyczaj pochodzące z biopsji nowotworu) lub `normal` (pobrane ze zdrowej tkanki lub krwi).
Jeśli nie znasz analizy nowotworów, wiedz tylko, że odpowiada to modelowi eksperymentalnemu, który wykorzystuje sparowane próbki nowotworowe/normalne do przeprowadzania analiz kontrastowych.

Dla pacjenta A mamy dwa zestawy replik technicznych (powtórzeń).

!!! note

    Nie martw się, jeśli nie znasz tego projektu eksperymentalnego, nie jest to krytyczne dla zrozumienia tego samouczka.

#### Przejrzyj zadanie

Twoim wyzwaniem jest napisanie workflow'a Nextflow'a, który:

1. **Odczyta** dane próbek z pliku CSV i ustrukturyzuje je za pomocą map metadanych
2. **Rozdzieli** próbki na różne kanały na podstawie typu (normalne vs nowotworowe)
3. **Połączy** dopasowane pary nowotworowe/normalne według identyfikatora pacjenta i numeru repliki
4. **Rozdzieli** próbki na interwały genomowe do przetwarzania równoległego
5. **Zgrupuje** powiązane próbki z powrotem do analizy dalszej

Reprezentuje to powszechny wzorzec bioinformatyczny, w którym musisz podzielić dane do niezależnego przetwarzania, a następnie ponownie połączyć powiązane elementy do analizy porównawczej.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moja przestrzeń kodowa działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Odczytaj dane próbek

### 1.1. Odczytaj dane próbek za pomocą `splitCsv` i utwórz mapy metadanych

Zacznijmy od odczytania danych próbek za pomocą `splitCsv` i zorganizowania ich we wzorzec mapy metadanych. W `main.nf` zobaczysz, że już rozpoczęliśmy workflow'a.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    W całym tym samouczku będziemy używać prefiksu `ch_` dla wszystkich zmiennych kanałów, aby wyraźnie wskazać, że są to kanały Nextflow'a.

Jeśli ukończyłeś misję poboczną [Metadane w workflow'ach](./metadata.md), rozpoznasz ten wzorzec. Użyjemy `splitCsv` do odczytu CSV i natychmiastowego ustrukturyzowania danych za pomocą mapy metadanych, aby oddzielić metadane od ścieżek plików.

!!! info

    W tym szkoleniu napotkamy dwa różne koncepty zwane `map`:

    - **Struktura danych**: Mapa Groovy (odpowiednik słowników/haszy w innych językach), która przechowuje pary klucz-wartość
    - **Operator kanału**: Operator `.map()`, który przekształca elementy w kanale

    Będziemy wyjaśniać, który z nich mamy na myśli w kontekście, ale to rozróżnienie jest ważne do zrozumienia podczas pracy z Nextflow'em.

Zastosuj te zmiany do `main.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

To łączy operację `splitCsv` (odczyt CSV z nagłówkami) i operację `map` (strukturyzowanie danych jako krotek `[meta, file]`) w jednym kroku. Zastosuj tę zmianę i uruchom pipeline'a:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Mamy teraz kanał, w którym każdy element jest krotką `[meta, file]` - metadane oddzielone od ścieżek plików. Ta struktura pozwala nam dzielić i grupować nasze obciążenie na podstawie pól metadanych.

---

## 2. Filtruj i przekształcaj dane

### 2.1. Filtruj dane za pomocą `filter`

Możemy użyć [operatora `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) do filtrowania danych na podstawie warunku. Powiedzmy, że chcemy przetwarzać tylko próbki normalne. Możemy to zrobić, filtrując dane na podstawie pola `type`. Wstawmy to przed operatorem `view`.

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Uruchom workflow'a ponownie, aby zobaczyć przefiltrowany wynik:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Pomyślnie przefiltrowaliśmy dane, aby uwzględnić tylko próbki normalne. Podsumujmy, jak to działa.

Operator `filter` przyjmuje domknięcie, które jest stosowane do każdego elementu w kanale. Jeśli domknięcie zwraca `true`, element jest uwzględniany; jeśli zwraca `false`, element jest wykluczany.

W naszym przypadku chcemy zachować tylko próbki, gdzie `meta.type == 'normal'`. Domknięcie używa krotki `meta,file` do odniesienia się do każdej próbki, uzyskuje dostęp do typu próbki za pomocą `meta.type` i sprawdza, czy równa się `'normal'`.

Osiąga się to za pomocą pojedynczego domknięcia, które wprowadziliśmy powyżej:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Utwórz oddzielne przefiltrowane kanały

Obecnie stosujemy filtr do kanału utworzonego bezpośrednio z CSV, ale chcemy filtrować to na więcej niż jeden sposób, więc przepiszmy logikę, aby utworzyć oddzielny przefiltrowany kanał dla próbek normalnych:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Uruchom pipeline'a, aby zobaczyć wyniki:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Pomyślnie przefiltrowaliśmy dane i utworzyliśmy oddzielny kanał dla próbek normalnych.

Utwórzmy również przefiltrowany kanał dla próbek nowotworowych:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Rozdzieliliśmy próbki normalne i nowotworowe na dwa różne kanały i użyliśmy domknięcia dostarczonego do `view()`, aby oznaczyć je inaczej w wyjściu: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Filtrowania danych**: Jak filtrować dane za pomocą `filter`
- **Dzielenia danych**: Jak dzielić dane na różne kanały na podstawie warunku
- **Wyświetlania danych**: Jak używać `view` do wypisywania danych i oznaczania wyjścia z różnych kanałów

Rozdzieliliśmy teraz próbki normalne i nowotworowe na dwa różne kanały. Następnie połączymy próbki normalne i nowotworowe na podstawie pola `id`.

---

## 3. Łączenie kanałów według identyfikatorów

W poprzedniej sekcji rozdzieliliśmy próbki normalne i nowotworowe na dwa różne kanały. Mogłyby być przetwarzane niezależnie przy użyciu określonych procesów lub workflow'ów na podstawie ich typu. Ale co się dzieje, gdy chcemy porównać próbki normalne i nowotworowe od tego samego pacjenta? W tym momencie musimy połączyć je z powrotem, upewniając się, że dopasowujemy próbki na podstawie ich pola `id`.

Nextflow zawiera wiele metod łączenia kanałów, ale w tym przypadku najbardziej odpowiednim operatorem jest [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Jeśli znasz SQL, działa jak operacja `JOIN`, gdzie określamy klucz do łączenia i typ łączenia do wykonania.

### 3.1. Użyj `map` i `join` do łączenia na podstawie identyfikatora pacjenta

Jeśli sprawdzimy dokumentację [`join`](https://www.nextflow.io/docs/latest/operator.html#join), zobaczymy, że domyślnie łączy dwa kanały na podstawie pierwszego elementu w każdej krotce.

#### 3.1.1. Sprawdź strukturę danych

Jeśli nie masz już dostępnego wyjścia konsoli, uruchommy pipeline'a, aby sprawdzić naszą strukturę danych i zobaczyć, jak musimy ją zmodyfikować, aby łączyć na podstawie pola `id`.

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Widzimy, że pole `id` jest pierwszym elementem w każdej mapie metadanych. Aby `join` działał, powinniśmy wyizolować pole `id` w każdej krotce. Następnie możemy po prostu użyć operatora `join` do połączenia dwóch kanałów.

#### 3.1.2. Wyizoluj pole `id`

Aby wyizolować pole `id`, możemy użyć [operatora `map`](https://www.nextflow.io/docs/latest/operator.html#map) do utworzenia nowej krotki z polem `id` jako pierwszym elementem.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Może to być subtelne, ale powinieneś być w stanie zobaczyć, że pierwszy element w każdej krotce to pole `id`.

#### 3.1.3. Połącz dwa kanały

Teraz możemy użyć operatora `join` do połączenia dwóch kanałów na podstawie pola `id`.

Ponownie użyjemy `view` do wypisania połączonych wyjść.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Trudno to stwierdzić, ponieważ jest tak szerokie, ale powinieneś być w stanie zobaczyć, że próbki zostały połączone według pola `id`. Każda krotka ma teraz format:

- `id`: Identyfikator próbki
- `normal_meta_map`: Metadane próbki normalnej, w tym typ, replika i ścieżka do pliku bam
- `normal_sample_file`: Plik próbki normalnej
- `tumor_meta_map`: Metadane próbki nowotworowej, w tym typ, replika i ścieżka do pliku bam
- `tumor_sample`: Próbka nowotworowa, w tym typ, replika i ścieżka do pliku bam

!!! warning

    Operator `join` odrzuci wszelkie niedopasowane krotki. W tym przykładzie upewniliśmy się, że wszystkie próbki były dopasowane dla nowotworu i normalnych, ale jeśli to nie jest prawda, musisz użyć parametru `remainder: true`, aby zachować niedopasowane krotki. Sprawdź [dokumentację](https://www.nextflow.io/docs/latest/operator.html#join), aby uzyskać więcej szczegółów.

Teraz wiesz, jak używać `map` do izolowania pola w krotce i jak używać `join` do łączenia krotek na podstawie pierwszego pola.
Dzięki tej wiedzy możemy pomyślnie łączyć kanały na podstawie wspólnego pola.

Następnie rozważymy sytuację, w której chcesz łączyć na podstawie wielu pól.

### 3.2. Łącz na podstawie wielu pól

Mamy 2 repliki dla próbkiA, ale tylko 1 dla próbkiB i próbkiC. W tym przypadku mogliśmy je skutecznie połączyć, używając pola `id`, ale co by się stało, gdyby były niezsynchronizowane? Moglibyśmy pomylić próbki normalne i nowotworowe z różnych replik!

Aby tego uniknąć, możemy łączyć na podstawie wielu pól. Istnieje wiele sposobów, aby to osiągnąć, ale skupimy się na utworzeniu nowego klucza łączącego, który zawiera zarówno `id` próbki, jak i numer `replicate`.

Zacznijmy od utworzenia nowego klucza łączącego. Możemy to zrobić w ten sam sposób, co wcześniej, używając [operatora `map`](https://www.nextflow.io/docs/latest/operator.html#map) do utworzenia nowej krotki z polami `id` i `repeat` jako pierwszym elementem.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Teraz powinniśmy zobaczyć, że łączenie następuje, ale przy użyciu zarówno pól `id`, jak i `repeat`. Uruchom workflow'a:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Zauważ, jak mamy krotkę dwóch elementów (pola `id` i `repeat`) jako pierwszy element każdego połączonego wyniku. To pokazuje, jak złożone elementy mogą być używane jako klucz łączący, umożliwiając dość skomplikowane dopasowywanie między próbkami z tych samych warunków.

Jeśli chcesz zbadać więcej sposobów łączenia na różnych kluczach, sprawdź [dokumentację operatora join](https://www.nextflow.io/docs/latest/operator.html#join), aby uzyskać dodatkowe opcje i przykłady.

### 3.3. Użyj `subMap` do utworzenia nowego klucza łączącego

Poprzednie podejście traci nazwy pól z naszego klucza łączącego - pola `id` i `repeat` stają się tylko listą wartości. Aby zachować nazwy pól do późniejszego dostępu, możemy użyć [metody `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

Metoda `subMap` wyodrębnia tylko określone pary klucz-wartość z mapy. Tutaj wyodrębnimy tylko pola `id` i `repeat`, aby utworzyć nasz klucz łączący.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Teraz mamy nowy klucz łączący, który nie tylko zawiera pola `id` i `repeat`, ale także zachowuje nazwy pól, dzięki czemu możemy uzyskać do nich dostęp później według nazwy, np. `meta.id` i `meta.repeat`.

### 3.4. Użyj nazwanego domknięcia w map

Aby uniknąć duplikacji i zmniejszyć błędy, możemy użyć nazwanego domknięcia. Nazwane domknięcie pozwala nam utworzyć funkcję wielokrotnego użytku, którą możemy wywołać w wielu miejscach.

Aby to zrobić, najpierw definiujemy domknięcie jako nową zmienną:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Zdefiniowaliśmy transformację map jako nazwaną zmienną, którą możemy ponownie wykorzystać.

Zauważ, że konwertujemy również ścieżkę pliku na obiekt Path przy użyciu `file()`, aby każdy proces otrzymujący ten kanał mógł poprawnie obsłużyć plik (więcej informacji znajdziesz w [Praca z plikami](./working_with_files.md)).

Zaimplementujmy domknięcie w naszym workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Przed"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note

    Operator `map` przełączył się z używania `{ }` na używanie `( )` do przekazywania domknięcia jako argumentu. Dzieje się tak, ponieważ operator `map` oczekuje domknięcia jako argumentu, a `{ }` służy do definiowania anonimowego domknięcia. Podczas wywoływania nazwanego domknięcia użyj składni `( )`.

Uruchom workflow'a jeszcze raz, aby sprawdzić, czy wszystko nadal działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Używanie nazwanego domknięcia pozwala nam ponownie wykorzystać tę samą transformację w wielu miejscach, zmniejszając ryzyko błędów i czyniąc kod bardziej czytelnym i łatwiejszym w utrzymaniu.

### 3.5. Zmniejsz duplikację danych

Mamy dużo zduplikowanych danych w naszym workflow'ie. Każdy element w połączonych próbkach powtarza pola `id` i `repeat`. Ponieważ te informacje są już dostępne w kluczu grupującym, możemy uniknąć tej redundancji. Jako przypomnienie, nasza obecna struktura danych wygląda tak:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Ponieważ pola `id` i `repeat` są dostępne w kluczu grupującym, usuńmy je z reszty każdego elementu kanału, aby uniknąć duplikacji. Możemy to zrobić, używając metody `subMap` do utworzenia nowej mapy tylko z polem `type`. To podejście pozwala nam zachować wszystkie niezbędne informacje, eliminując jednocześnie redundancję w naszej strukturze danych.

=== "Po"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Teraz domknięcie zwraca krotkę, w której pierwszy element zawiera pola `id` i `repeat`, a drugi element zawiera tylko pole `type`. Eliminuje to redundancję poprzez przechowywanie informacji `id` i `repeat` raz w kluczu grupującym, zachowując jednocześnie wszystkie niezbędne informacje.

Uruchom workflow'a, aby zobaczyć, jak to wygląda:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Widzimy, że podajemy pola `id` i `repeat` tylko raz w kluczu grupującym i mamy pole `type` w danych próbki. Nie straciliśmy żadnych informacji, ale udało nam się uczynić zawartość naszego kanału bardziej zwięzłą.

### 3.6. Usuń zbędne informacje

Usunęliśmy zduplikowane informacje powyżej, ale nadal mamy inne zbędne informacje w naszych kanałach.

Na początku rozdzieliliśmy próbki normalne i nowotworowe za pomocą `filter`, a następnie połączyliśmy je na podstawie kluczy `id` i `repeat`. Operator `join` zachowuje kolejność, w jakiej krotki są scalane, więc w naszym przypadku, z próbkami normalnymi po lewej stronie i próbkami nowotworowymi po prawej, wynikowy kanał zachowuje tę strukturę: `id, <elementy normalne>, <elementy nowotworowe>`.

Ponieważ znamy pozycję każdego elementu w naszym kanale, możemy jeszcze bardziej uprościć strukturę, usuwając metadane `[type:normal]` i `[type:tumor]`.

=== "Po"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Uruchom ponownie, aby zobaczyć wynik:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Manipulowania krotkami**: Jak używać `map` do izolowania pola w krotce
- **Łączenia krotek**: Jak używać `join` do łączenia krotek na podstawie pierwszego pola
- **Tworzenia kluczy łączących**: Jak używać `subMap` do tworzenia nowego klucza łączącego
- **Nazwanych domknięć**: Jak używać nazwanego domknięcia w map
- **Łączenia wielu pól**: Jak łączyć na podstawie wielu pól dla bardziej precyzyjnego dopasowania
- **Optymalizacji struktury danych**: Jak usprawnić strukturę kanału poprzez usunięcie zbędnych informacji

Masz teraz workflow'a, który może podzielić arkusz próbek, przefiltrować próbki normalne i nowotworowe, połączyć je według identyfikatora próbki i numeru repliki, a następnie wypisać wyniki.

Jest to powszechny wzorzec w workflow'ach bioinformatycznych, gdzie musisz dopasować próbki lub inne typy danych po niezależnym przetwarzaniu, więc jest to przydatna umiejętność. Następnie przyjrzymy się powtarzaniu próbki wiele razy.

## 4. Rozdziel próbki na interwały

Kluczowym wzorcem w workflow'ach bioinformatycznych jest dystrybucja analizy na regiony genomowe. Na przykład wykrywanie wariantów może być zrównoleglone poprzez podzielenie genomu na interwały (takie jak chromosomy lub mniejsze regiony). Ta strategia paralelizacji znacznie poprawia wydajność pipeline'a poprzez rozłożenie obciążenia obliczeniowego na wiele rdzeni lub węzłów, skracając całkowity czas wykonania.

W następnej sekcji pokażemy, jak rozdzielić nasze dane próbek na wiele interwałów genomowych. Sparujemy każdą próbkę z każdym interwałem, umożliwiając równoległe przetwarzanie różnych regionów genomowych. To pomnoży rozmiar naszego zestawu danych przez liczbę interwałów, tworząc wiele niezależnych jednostek analizy, które można później połączyć z powrotem.

### 4.1. Rozdziel próbki na interwały przy użyciu `combine`

Zacznijmy od utworzenia kanału interwałów. Aby uprościć życie, użyjemy tylko 3 interwałów, które zdefiniujemy ręcznie. W prawdziwym workflow'ie możesz je odczytać z pliku wejściowego lub nawet utworzyć kanał z wieloma plikami interwałów.

=== "Po"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Teraz pamiętaj, chcemy powtórzyć każdą próbkę dla każdego interwału. Jest to czasami określane jako iloczyn kartezjański próbek i interwałów. Możemy to osiągnąć, używając [operatora `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Weźmie on każdy element z kanału 1 i powtórzy go dla każdego elementu w kanale 2. Dodajmy operator combine do naszego workflow'a:

=== "Po"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Teraz uruchommy to i zobaczmy, co się stanie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Sukces! Powtórzyliśmy każdą próbkę dla każdego pojedynczego interwału na naszej liście 3 interwałów. Skutecznie potroiliśmy liczbę elementów w naszym kanale.

Trudno to jednak odczytać, więc w następnej sekcji to uporządkujemy.

### 4.2. Zorganizuj kanał

Możemy użyć operatora `map` do uporządkowania i refaktoryzacji naszych danych próbek, aby były łatwiejsze do zrozumienia. Przenieśmy ciąg interwałów do mapy łączącej w pierwszym elemencie.

=== "Po"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Rozłóżmy na czynniki pierwsze, co robi ta operacja map krok po kroku.

Najpierw używamy nazwanych parametrów, aby uczynić kod bardziej czytelnym. Używając nazw `grouping_key`, `normal`, `tumor` i `interval`, możemy odwoływać się do elementów w krotce według nazwy zamiast według indeksu:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Następnie łączymy `grouping_key` z polem `interval`. `grouping_key` jest mapą zawierającą pola `id` i `repeat`. Tworzymy nową mapę z `interval` i scalamy je za pomocą dodawania map Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Na koniec zwracamy to jako krotkę z trzema elementami: połączoną mapą metadanych, plikiem próbki normalnej i plikiem próbki nowotworowej:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Uruchommy to ponownie i sprawdźmy zawartość kanału:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Używanie `map` do wymuszenia odpowiedniej struktury danych może być trudne, ale jest kluczowe dla efektywnej manipulacji danymi.

Mamy teraz każdą próbkę powtórzoną we wszystkich interwałach genomowych, tworząc wiele niezależnych jednostek analizy, które mogą być przetwarzane równolegle. Ale co, jeśli chcemy połączyć powiązane próbki z powrotem? W następnej sekcji nauczymy się, jak grupować próbki, które mają wspólne atrybuty.

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Rozdzielania próbek na interwały**: Jak używać `combine` do powtarzania próbek na interwałach
- **Tworzenia iloczynów kartezjańskich**: Jak generować wszystkie kombinacje próbek i interwałów
- **Organizowania struktury kanału**: Jak używać `map` do restrukturyzacji danych dla lepszej czytelności
- **Przygotowania do przetwarzania równoległego**: Jak skonfigurować dane do analizy rozproszonej

## 5. Agregowanie próbek przy użyciu `groupTuple`

W poprzednich sekcjach nauczyliśmy się, jak dzielić dane z pliku wejściowego i filtrować według określonych pól (w naszym przypadku próbki normalne i nowotworowe). Ale to obejmuje tylko jeden typ łączenia. Co, jeśli chcemy grupować próbki według określonego atrybutu? Na przykład, zamiast łączyć dopasowane pary normalny-nowotworowy, możemy chcieć przetwarzać wszystkie próbki z "próbkiA" razem, niezależnie od ich typu. Ten wzorzec jest powszechny w workflow'ach bioinformatycznych, gdzie możesz chcieć przetwarzać powiązane próbki oddzielnie ze względu na wydajność przed porównaniem lub połączeniem wyników na końcu.

Nextflow zawiera wbudowane metody do tego, główną, na którą spojrzymy, jest `groupTuple`.

Zacznijmy od zgrupowania wszystkich naszych próbek, które mają te same pola `id` i `interval`, byłoby to typowe dla analizy, w której chcielibyśmy zgrupować repliki techniczne, ale zachować znacząco różne próbki oddzielone.

Aby to zrobić, powinniśmy wyodrębnić nasze zmienne grupujące, abyśmy mogli ich używać w izolacji.

Pierwszy krok jest podobny do tego, co zrobiliśmy w poprzedniej sekcji. Musimy wyizolować naszą zmienną grupującą jako pierwszy element krotki. Pamiętaj, nasz pierwszy element jest obecnie mapą pól `id`, `repeat` i `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Możemy ponownie użyć metody `subMap` z wcześniej, aby wyizolować nasze pola `id` i `interval` z mapy. Jak wcześniej, użyjemy operatora `map`, aby zastosować metodę `subMap` do pierwszego elementu krotki dla każdej próbki.

=== "Po"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Uruchommy to ponownie i sprawdźmy zawartość kanału:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Widzimy, że pomyślnie wyizolowaliśmy pola `id` i `interval`, ale jeszcze nie zgrupowaliśmy próbek.

!!! note

    Odrzucamy tutaj pole `replicate`. Dzieje się tak, ponieważ nie potrzebujemy go do dalszego przetwarzania. Po ukończeniu tego samouczka sprawdź, czy możesz je uwzględnić bez wpływu na późniejsze grupowanie!

Zgrupujmy teraz próbki według tego nowego elementu grupującego, używając [operatora `groupTuple`](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "Po"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

To wszystko! Dodaliśmy tylko jedną linię kodu. Zobaczmy, co się stanie, gdy to uruchomimy:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Zauważ, że nasza struktura danych uległa zmianie i w każdym elemencie kanału pliki są teraz zawarte w krotkach takich jak `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. Dzieje się tak, ponieważ gdy używamy `groupTuple`, Nextflow łączy pojedyncze pliki dla każdej próbki grupy. Jest to ważne do zapamiętania podczas próby obsługi danych dalej.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) jest przeciwieństwem groupTuple. Rozpakowuje elementy w kanale i spłaszcza je. Spróbuj dodać `transpose` i cofnąć grupowanie, które wykonaliśmy powyżej!

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Grupowania powiązanych próbek**: Jak używać `groupTuple` do agregowania próbek według wspólnych atrybutów
- **Izolowania kluczy grupujących**: Jak używać `subMap` do wyodrębniania określonych pól do grupowania
- **Obsługi zgrupowanych struktur danych**: Jak pracować z zagnieżdżoną strukturą utworzoną przez `groupTuple`
- **Obsługi replik technicznych**: Jak grupować próbki, które mają te same warunki eksperymentalne

---

## Podsumowanie

W tej misji pobocznej nauczyłeś się, jak dzielić i grupować dane przy użyciu kanałów.

Modyfikując dane w miarę ich przepływu przez pipeline'a, możesz skonstruować skalowalny pipeline'a bez używania pętli lub instrukcji while, oferując kilka zalet w porównaniu z bardziej tradycyjnymi podejściami:

- Możemy skalować do tylu lub tak niewielu wejść, ile chcemy, bez dodatkowego kodu
- Skupiamy się na obsłudze przepływu danych przez pipeline'a, zamiast iteracji
- Możemy być tak złożeni lub prości, jak to konieczne
- Pipeline'a staje się bardziej deklaratywny, skupiając się na tym, co powinno się wydarzyć, a nie jak powinno się wydarzyć
- Nextflow zoptymalizuje wykonanie dla nas, uruchamiając niezależne operacje równolegle

Opanowanie tych operacji kanałów umożliwi Ci budowanie elastycznych, skalowalnych pipeline'ów, które obsługują złożone relacje danych bez uciekania się do pętli lub programowania iteracyjnego, pozwalając Nextflow'owi optymalizować wykonanie i paralelizować niezależne operacje automatycznie.

### Kluczowe wzorce

1.  **Tworzenie ustrukturyzowanych danych wejściowych:** Rozpoczynając od pliku CSV z mapami metadanych (budując na wzorcach z [Metadane w workflow'ach](./metadata.md))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Dzielenie danych na oddzielne kanały:** Użyliśmy `filter` do podzielenia danych na niezależne strumienie na podstawie pola `type`

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Łączenie dopasowanych próbek:** Użyliśmy `join` do ponownego połączenia powiązanych próbek na podstawie pól `id` i `repeat`

    - Połącz dwa kanały według klucza (pierwszy element krotki)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Wyodrębnij klucz łączący i połącz według tej wartości

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Łącz na podstawie wielu pól przy użyciu subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Dystrybucja na interwały:** Użyliśmy `combine` do tworzenia iloczynów kartezjańskich próbek z interwałami genomowymi do przetwarzania równoległego.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Agregowanie według kluczy grupujących:** Użyliśmy `groupTuple` do grupowania według pierwszego elementu w każdej krotce, zbierając w ten sposób próbki dzielące pola `id` i `interval` oraz scalając repliki techniczne.

    ```groovy
    channel.groupTuple()
    ```

6.  **Optymalizacja struktury danych:** Użyliśmy `subMap` do wyodrębniania określonych pól i utworzyliśmy nazwane domknięcie do tworzenia transformacji wielokrotnego użytku.

    - Wyodrębnij określone pola z mapy

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Użyj nazwanego domknięcia dla transformacji wielokrotnego użytku

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Dodatkowe zasoby

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Co dalej?

Wróć do [menu Misji Pobocznych](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
