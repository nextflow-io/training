# Podział i grupowanie

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow oferuje zaawansowane narzędzia do elastycznej pracy z danymi. Kluczową możliwością jest podział danych na różne strumienie, a następnie ponowne grupowanie powiązanych elementów. Jest to szczególnie przydatne w workflow'ach bioinformatycznych, gdzie trzeba przetwarzać różne typy próbek osobno, zanim połączy się wyniki do analizy.

Wyobraź sobie sortowanie poczty: rozdzielasz listy według miejsca przeznaczenia, przetwarzasz każdy stos inaczej, a następnie łączysz przesyłki kierowane do tej samej osoby. Nextflow używa specjalnych operatorów do realizacji tego wzorca na danych naukowych. Podejście to jest powszechnie znane jako wzorzec **scatter/gather** w obliczeniach rozproszonych i workflow'ach bioinformatycznych.

System kanałów Nextflow'a jest sercem tej elastyczności. Kanały łączą różne części workflow'u, umożliwiając przepływ danych przez analizę. Możesz tworzyć wiele kanałów z jednego źródła danych, przetwarzać każdy z nich inaczej, a następnie łączyć je z powrotem, gdy zajdzie taka potrzeba. Takie podejście pozwala projektować workflow'y, które naturalnie odzwierciedlają rozgałęziające się i zbiegające się ścieżki złożonych analiz bioinformatycznych.

### Cele szkolenia

W tym zadaniu pobocznym nauczysz się dzielić i grupować dane przy użyciu operatorów kanałów Nextflow'a.
Zaczniemy od pliku CSV zawierającego informacje o próbkach i powiązanych plikach danych, a następnie będziemy manipulować tymi danymi i je reorganizować.

Po ukończeniu tego zadania będziesz potrafić efektywnie rozdzielać i łączyć strumienie danych, korzystając z następujących technik:

- Odczytywanie danych z plików przy użyciu `splitCsv`
- Filtrowanie i przekształcanie danych za pomocą `filter` i `map`
- Łączenie powiązanych danych przy użyciu `join` i `groupTuple`
- Tworzenie kombinacji danych za pomocą `combine` do przetwarzania równoległego
- Optymalizacja struktury danych przy użyciu `subMap` i strategii deduplikacji
- Budowanie funkcji wielokrotnego użytku za pomocą nazwanych domknięć do manipulowania strukturami kanałów

Te umiejętności pomogą Ci budować workflow'y, które efektywnie obsługują wiele plików wejściowych i różne typy danych, zachowując przy tym przejrzystą i łatwą w utrzymaniu strukturę kodu.

### Wymagania wstępne

Przed przystąpieniem do tego zadania pobocznego powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow'a (procesy, kanały, operatory, praca z plikami, metadane).

**Opcjonalnie:** Zalecamy wcześniejsze ukończenie zadania pobocznego [Metadata in workflows](../metadata/).
Omawia ono podstawy odczytywania plików CSV za pomocą `splitCsv` i tworzenia map metadanych, z których będziemy tu intensywnie korzystać.

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/splitting_and_grouping
```

Możesz ustawić VSCode tak, aby skupiał się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz tu główny plik workflow'u oraz katalog `data` zawierający arkusz próbek o nazwie `samplesheet.csv`.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Arkusz próbek zawiera informacje o próbkach od różnych pacjentów, w tym identyfikator pacjenta, numer powtórzenia, typ (normalny lub nowotworowy) oraz ścieżki do hipotetycznych plików danych (które w rzeczywistości nie istnieją, ale będziemy udawać, że tak).

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

Arkusz zawiera osiem próbek od trzech pacjentów (A, B, C).

Dla każdego pacjenta mamy próbki typu `tumor` (zazwyczaj pochodzące z biopsji guza) lub `normal` (pobrane ze zdrowej tkanki lub krwi).
Jeśli nie jesteś zaznajomiony z analizą nowotworów, wiedz tylko, że odpowiada to modelowi eksperymentalnemu, który wykorzystuje sparowane próbki nowotworowe i normalne do przeprowadzania analiz kontrastywnych.

Dla pacjenta A mamy dwa zestawy replik technicznych (powtórzeń).

!!! note "Uwaga"

    Nie martw się, jeśli nie znasz tego projektu eksperymentalnego — nie jest to kluczowe dla zrozumienia tego samouczka.

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest napisanie workflow'u Nextflow'a, który:

1. **Odczyta** dane próbek z pliku CSV i ustrukturyzuje je za pomocą map metadanych
2. **Rozdzieli** próbki na różne kanały na podstawie typu (normalny vs nowotworowy)
3. **Połączy** dopasowane pary nowotworowe/normalne według identyfikatora pacjenta i numeru repliki
4. **Rozdzieli** próbki na interwały genomiczne do przetwarzania równoległego
5. **Zgrupuje** powiązane próbki z powrotem do dalszej analizy

Reprezentuje to typowy wzorzec bioinformatyczny, w którym trzeba podzielić dane do niezależnego przetwarzania, a następnie ponownie połączyć powiązane elementy do analizy porównawczej.

#### Lista kontrolna gotowości

Gotowy do działania?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiedni katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Wczytywanie danych próbek

### 1.1. Wczytaj dane próbek za pomocą `splitCsv` i utwórz mapy metadanych

Zacznijmy od wczytania danych próbek za pomocą `splitCsv` i zorganizowania ich zgodnie ze wzorcem mapy metadanych. W pliku `main.nf` zobaczysz, że workflow jest już częściowo napisany.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Uwaga"

    W całym tym samouczku będziemy używać prefiksu `ch_` dla wszystkich zmiennych kanałów, aby wyraźnie wskazać, że są to kanały Nextflow'a.

Jeśli ukończyłeś zadanie poboczne [Metadata in workflows](../metadata/), rozpoznasz ten wzorzec. Użyjemy `splitCsv` do odczytania pliku CSV i natychmiastowego ustrukturyzowania danych za pomocą mapy metadanych, aby oddzielić metadane od ścieżek do plików.

!!! info "Info"

    W tym szkoleniu napotkamy dwa różne koncepty nazywane `map`:

    - **Struktura danych**: Mapa Groovy (odpowiednik słowników/tablic asocjacyjnych w innych językach) przechowująca pary klucz-wartość
    - **Operator kanału**: Operator `.map()`, który przekształca elementy w kanale

    W razie potrzeby wyjaśnimy, o który z nich chodzi, ale to rozróżnienie jest ważne przy pracy z Nextflow'em.

Wprowadź następujące zmiany w `main.nf`:

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

Łączy to operację `splitCsv` (odczyt CSV z nagłówkami) i operację `map` (strukturyzowanie danych jako krotek `[meta, plik]`) w jednym kroku. Wprowadź tę zmianę i uruchom pipeline:

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

Mamy teraz kanał, w którym każdy element jest krotką `[meta, plik]` — metadane oddzielone od ścieżek do plików. Ta struktura pozwala nam dzielić i grupować dane na podstawie pól metadanych.

---

## 2. Filtrowanie i przekształcanie danych

### 2.1. Filtrowanie danych za pomocą `filter`

Możemy użyć [operatora `filter`](https://www.nextflow.io/docs/latest/operator.html#filter), aby filtrować dane na podstawie warunku. Powiedzmy, że chcemy przetwarzać tylko próbki normalne. Możemy to zrobić, filtrując dane na podstawie pola `type`. Wstawmy to przed operatorem `view`.

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

Uruchom workflow ponownie, aby zobaczyć przefiltrowany wynik:

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

Pomyślnie przefiltrowano dane, aby uwzględniały tylko próbki normalne. Przypomnijmy, jak to działa.

Operator `filter` przyjmuje domknięcie, które jest stosowane do każdego elementu w kanale. Jeśli domknięcie zwraca `true`, element jest uwzględniany; jeśli zwraca `false`, element jest wykluczany.

W naszym przypadku chcemy zachować tylko próbki, dla których `meta.type == 'normal'`. Domknięcie używa krotki `meta,file` do odwołania się do każdej próbki, uzyskuje dostęp do jej typu przez `meta.type` i sprawdza, czy jest równy `'normal'`.

Realizuje to pojedyncze domknięcie wprowadzone powyżej:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Tworzenie osobnych przefiltrowanych kanałów

Aktualnie stosujemy filtr do kanału tworzonego bezpośrednio z pliku CSV, ale chcemy filtrować dane na więcej sposobów. Przepiszmy więc logikę, aby utworzyć osobny przefiltrowany kanał dla próbek normalnych:

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

Uruchom pipeline, aby zobaczyć wyniki:

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

Pomyślnie przefiltrowano dane i utworzono osobny kanał dla próbek normalnych.

Utwórzmy teraz przefiltrowany kanał również dla próbek nowotworowych:

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

Rozdzieliliśmy próbki normalne i nowotworowe na dwa różne kanały i użyliśmy domknięcia przekazanego do `view()`, aby oznaczyć je inaczej w wyjściu: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Filtrowania danych**: Jak filtrować dane za pomocą `filter`
- **Podziału danych**: Jak dzielić dane na różne kanały na podstawie warunku
- **Wyświetlania danych**: Jak używać `view` do wydrukowania danych i oznaczania wyjścia z różnych kanałów

Rozdzieliliśmy próbki normalne i nowotworowe na dwa różne kanały. W następnej sekcji połączymy je z powrotem na podstawie pola `id`.

---

## 3. Łączenie kanałów według identyfikatorów

W poprzedniej sekcji rozdzieliliśmy próbki normalne i nowotworowe na dwa różne kanały. Mogłyby być przetwarzane niezależnie przy użyciu specyficznych procesów lub workflow'ów zależnie od ich typu. Co jednak, gdy chcemy porównać próbki normalne i nowotworowe od tego samego pacjenta? W tym momencie musimy je z powrotem połączyć, upewniając się, że próbki są dopasowane na podstawie pola `id`.

Nextflow oferuje wiele metod łączenia kanałów, ale w tym przypadku najbardziej odpowiednim operatorem jest [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Jeśli znasz SQL, działa on jak operacja `JOIN`, gdzie określamy klucz łączenia i typ złączenia.

### 3.1. Użyj `map` i `join` do łączenia na podstawie identyfikatora pacjenta

Sprawdzając [dokumentację `join`](https://www.nextflow.io/docs/latest/operator.html#join), możemy zobaczyć, że domyślnie łączy dwa kanały na podstawie pierwszego elementu każdej krotki.

#### 3.1.1. Sprawdź strukturę danych

Jeśli nie masz już dostępnego wyjścia konsoli, uruchommy pipeline, aby sprawdzić strukturę danych i zobaczyć, jak musimy ją zmodyfikować, aby łączyć na podstawie pola `id`.

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

Widzimy, że pole `id` jest pierwszym elementem każdej mapy metadanych. Aby `join` działał, powinniśmy wyizolować pole `id` w każdej krotce. Następnie możemy po prostu użyć operatora `join` do połączenia dwóch kanałów.

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

Może to być subtelne, ale powinieneś zauważyć, że pierwszym elementem każdej krotki jest pole `id`.

#### 3.1.3. Połącz dwa kanały

Teraz możemy użyć operatora `join` do połączenia dwóch kanałów na podstawie pola `id`.

Ponownie użyjemy `view` do wydrukowania połączonych wyników.

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

Jest to trochę trudne do odczytania ze względu na szerokość, ale powinieneś zauważyć, że próbki zostały połączone na podstawie pola `id`. Każda krotka ma teraz format:

- `id`: Identyfikator próbki
- `normal_meta_map`: Metadane próbki normalnej, w tym typ, numer repliki i ścieżka do pliku BAM
- `normal_sample_file`: Plik próbki normalnej
- `tumor_meta_map`: Metadane próbki nowotworowej, w tym typ, numer repliki i ścieżka do pliku BAM
- `tumor_sample`: Próbka nowotworowa, w tym typ, numer repliki i ścieżka do pliku BAM

!!! warning "Ostrzeżenie"

    Operator `join` odrzuci wszelkie niedopasowane krotki. W tym przykładzie upewniliśmy się, że wszystkie próbki były dopasowane dla nowotworowych i normalnych, ale jeśli tak nie jest, musisz użyć parametru `remainder: true`, aby zachować niedopasowane krotki. Sprawdź [dokumentację](https://www.nextflow.io/docs/latest/operator.html#join), aby uzyskać więcej szczegółów.

Wiesz już, jak używać `map` do izolowania pola w krotce i jak używać `join` do łączenia krotek na podstawie pierwszego pola.
Dzięki tej wiedzy możemy skutecznie łączyć kanały na podstawie wspólnego pola.

Następnie rozważymy sytuację, w której chcemy łączyć na podstawie wielu pól.

### 3.2. Łączenie na podstawie wielu pól

Mamy 2 repliki dla próbki A, ale tylko 1 dla próbek B i C. W tym przypadku byliśmy w stanie skutecznie je połączyć, używając pola `id`, ale co by się stało, gdyby były niesynchronizowane? Moglibyśmy pomylić próbki normalne i nowotworowe z różnych replik!

Aby tego uniknąć, możemy łączyć na podstawie wielu pól. Istnieje kilka sposobów, aby to osiągnąć, ale skupimy się na tworzeniu nowego klucza łączenia, który zawiera zarówno `id` próbki, jak i numer `replicate`.

Zacznijmy od utworzenia nowego klucza łączenia. Możemy to zrobić w ten sam sposób co poprzednio, używając [operatora `map`](https://www.nextflow.io/docs/latest/operator.html#map) do utworzenia nowej krotki z polami `id` i `repeat` jako pierwszym elementem.

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

Teraz powinniśmy zobaczyć, że łączenie odbywa się z użyciem zarówno pól `id`, jak i `repeat`. Uruchom workflow:

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

Zwróć uwagę, że jako pierwszy element każdego połączonego wyniku mamy krotkę dwóch elementów (pola `id` i `repeat`). Pokazuje to, jak złożone elementy mogą być używane jako klucz łączenia, umożliwiając dość skomplikowane dopasowywanie próbek z tych samych warunków.

Jeśli chcesz poznać więcej sposobów łączenia na różnych kluczach, sprawdź [dokumentację operatora join](https://www.nextflow.io/docs/latest/operator.html#join), aby uzyskać dodatkowe opcje i przykłady.

### 3.3. Użyj `subMap` do tworzenia nowego klucza łączenia

Poprzednie podejście traci nazwy pól z klucza łączenia — pola `id` i `repeat` stają się zwykłą listą wartości. Aby zachować nazwy pól do późniejszego dostępu, możemy użyć [metody `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

Metoda `subMap` wyodrębnia tylko określone pary klucz-wartość z mapy. Tutaj wyodrębnimy tylko pola `id` i `repeat`, aby utworzyć nasz klucz łączenia.

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

Mamy teraz nowy klucz łączenia, który nie tylko zawiera pola `id` i `repeat`, ale także zachowuje nazwy pól, dzięki czemu możemy uzyskać do nich dostęp później przez nazwę, np. `meta.id` i `meta.repeat`.

### 3.4. Użyj nazwanego domknięcia w map

Aby uniknąć duplikacji i zmniejszyć ryzyko błędów, możemy użyć nazwanego domknięcia. Nazwane domknięcie pozwala nam tworzyć funkcje wielokrotnego użytku, które możemy wywoływać w wielu miejscach.

W tym celu najpierw definiujemy domknięcie jako nową zmienną:

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

Zwróć uwagę, że konwertujemy też ścieżkę do pliku na obiekt Path za pomocą `file()`, aby każdy proces otrzymujący ten kanał mógł poprawnie obsłużyć plik (więcej informacji znajdziesz w sekcji [Working with files](../working_with_files/)).

Zaimplementujmy domknięcie w naszym workflow'u:

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

!!! note "Uwaga"

    Operator `map` zmienił składnię z `{ }` na `( )` do przekazywania domknięcia jako argumentu. Wynika to z tego, że operator `map` oczekuje domknięcia jako argumentu, a `{ }` służy do definiowania anonimowego domknięcia. Przy wywoływaniu nazwanego domknięcia używaj składni `( )`.

Uruchom workflow jeszcze raz, aby sprawdzić, czy wszystko nadal działa:

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

Użycie nazwanego domknięcia pozwala nam ponownie wykorzystać tę samą transformację w wielu miejscach, zmniejszając ryzyko błędów i czyniąc kod bardziej czytelnym i łatwym w utrzymaniu.

### 3.5. Redukcja duplikacji danych

W naszym workflow'u mamy dużo zduplikowanych danych. Każdy element w połączonych próbkach powtarza pola `id` i `repeat`. Ponieważ informacje te są już dostępne w kluczu grupowania, możemy uniknąć tej redundancji. Dla przypomnienia, nasza aktualna struktura danych wygląda tak:

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

Ponieważ pola `id` i `repeat` są dostępne w kluczu grupowania, usuńmy je z pozostałych elementów każdego kanału, aby uniknąć duplikacji. Możemy to zrobić, używając metody `subMap` do utworzenia nowej mapy zawierającej tylko pole `type`. Takie podejście pozwala nam zachować wszystkie niezbędne informacje, eliminując jednocześnie redundancję w strukturze danych.

=== "Po"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Teraz domknięcie zwraca krotkę, w której pierwszy element zawiera pola `id` i `repeat`, a drugi element zawiera tylko pole `type`. Eliminuje to redundancję poprzez przechowywanie informacji `id` i `repeat` raz w kluczu grupowania, przy jednoczesnym zachowaniu wszystkich niezbędnych danych.

Uruchom workflow, aby zobaczyć, jak to wygląda:

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

Widzimy, że pola `id` i `repeat` są podane tylko raz w kluczu grupowania, a pole `type` jest zawarte w danych próbki. Nie utraciliśmy żadnych informacji, a zawartość kanału stała się bardziej zwięzła.

### 3.6. Usuwanie zbędnych informacji

Powyżej usunęliśmy zduplikowane informacje, ale w naszych kanałach nadal są pewne inne zbędne dane.

Na początku rozdzieliliśmy próbki normalne i nowotworowe za pomocą `filter`, a następnie połączyliśmy je na podstawie kluczy `id` i `repeat`. Operator `join` zachowuje kolejność, w jakiej krotki są łączone, więc w naszym przypadku — z próbkami normalnymi po lewej stronie i nowotworowymi po prawej — wynikowy kanał zachowuje tę strukturę: `id, <elementy normalne>, <elementy nowotworowe>`.

Ponieważ znamy pozycję każdego elementu w naszym kanale, możemy uprościć strukturę, usuwając metadane `[type:normal]` i `[type:tumor]`.

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
- **Tworzenia kluczy łączenia**: Jak używać `subMap` do tworzenia nowego klucza łączenia
- **Nazwanych domknięć**: Jak używać nazwanego domknięcia w map
- **Łączenia na wielu polach**: Jak łączyć na wielu polach dla bardziej precyzyjnego dopasowania
- **Optymalizacji struktury danych**: Jak upraszczać strukturę kanału przez usuwanie zbędnych informacji

Masz teraz workflow, który potrafi podzielić arkusz próbek, przefiltrować próbki normalne i nowotworowe, połączyć je według identyfikatora próbki i numeru repliki, a następnie wydrukować wyniki.

Jest to typowy wzorzec w workflow'ach bioinformatycznych, gdzie trzeba dopasowywać próbki lub inne typy danych po niezależnym przetworzeniu — to przydatna umiejętność. Następnie przyjrzymy się powtarzaniu próbki wiele razy.

## 4. Rozdzielanie próbek na interwały

Kluczowym wzorcem w workflow'ach bioinformatycznych jest rozdzielanie analizy na regiony genomiczne. Na przykład, wywoływanie wariantów można paralelizować, dzieląc genom na interwały (takie jak chromosomy lub mniejsze regiony). Ta strategia paralelizacji znacznie poprawia wydajność pipeline'u, rozdzielając obciążenie obliczeniowe na wiele rdzeni lub węzłów, co skraca całkowity czas wykonania.

W poniższej sekcji pokażemy, jak rozdzielić dane próbek na wiele interwałów genomicznych. Sparujemy każdą próbkę z każdym interwałem, umożliwiając równoległe przetwarzanie różnych regionów genomicznych. Zwiększy to rozmiar naszego zbioru danych o liczbę interwałów, tworząc wiele niezależnych jednostek analizy, które można później ponownie połączyć.

### 4.1. Rozdzielanie próbek na interwały za pomocą `combine`

Zacznijmy od utworzenia kanału interwałów. Dla uproszczenia użyjemy 3 interwałów zdefiniowanych ręcznie. W prawdziwym workflow'u można by je wczytać z pliku wejściowego lub nawet utworzyć kanał z wieloma plikami interwałów.

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

Pamiętaj, że chcemy powtórzyć każdą próbkę dla każdego interwału. Nazywa się to czasem iloczynem kartezjańskim próbek i interwałów. Możemy to osiągnąć, używając [operatora `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Pobiera on każdy element z kanału 1 i powtarza go dla każdego elementu z kanału 2. Dodajmy operator `combine` do naszego workflow'u:

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

Uruchommy to i zobaczmy, co się stanie:

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

Sukces! Powtórzyliśmy każdą próbkę dla każdego interwału z naszej listy 3 interwałów. Efektywnie potroiliśmy liczbę elementów w naszym kanale.

Jest to jednak trochę trudne do odczytania, więc w następnej sekcji to uporządkujemy.

### 4.2. Organizacja kanału

Możemy użyć operatora `map` do uporządkowania i refaktoryzacji danych próbek, aby były łatwiejsze do zrozumienia. Przenieśmy string interwału do mapy grupowania jako pierwszy element.

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

Przeanalizujmy krok po kroku, co robi ta operacja map.

Najpierw używamy nazwanych parametrów, aby kod był bardziej czytelny. Używając nazw `grouping_key`, `normal`, `tumor` i `interval`, możemy odwoływać się do elementów krotki przez nazwę zamiast przez indeks:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Następnie łączymy `grouping_key` z polem `interval`. `grouping_key` to mapa zawierająca pola `id` i `repeat`. Tworzymy nową mapę z `interval` i łączymy je za pomocą dodawania map w Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Na koniec zwracamy to jako krotkę z trzema elementami: połączoną mapę metadanych, plik próbki normalnej i plik próbki nowotworowej:

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

Używanie `map` do przekształcania danych w odpowiednią strukturę może być trudne, ale jest kluczowe dla efektywnej manipulacji danymi.

Mamy teraz każdą próbkę powtórzoną dla wszystkich interwałów genomicznych, tworząc wiele niezależnych jednostek analizy, które można przetwarzać równolegle. Co jednak, jeśli chcemy z powrotem zebrać powiązane próbki? W następnej sekcji nauczymy się grupować próbki, które mają wspólne atrybuty.

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Rozdzielania próbek na interwały**: Jak używać `combine` do powtarzania próbek dla interwałów
- **Tworzenia iloczynów kartezjańskich**: Jak generować wszystkie kombinacje próbek i interwałów
- **Organizowania struktury kanału**: Jak używać `map` do restrukturyzacji danych dla lepszej czytelności
- **Przygotowania do przetwarzania równoległego**: Jak przygotować dane do analizy rozproszonej

## 5. Agregowanie próbek za pomocą `groupTuple`

W poprzednich sekcjach nauczyliśmy się dzielić dane z pliku wejściowego i filtrować według określonych pól (w naszym przypadku próbki normalne i nowotworowe). Obejmuje to jednak tylko jeden typ łączenia. Co jeśli chcemy grupować próbki według określonego atrybutu? Na przykład, zamiast łączyć dopasowane pary normalny-nowotworowy, możemy chcieć przetwarzać wszystkie próbki z „próbki A" razem, niezależnie od ich typu. Ten wzorzec jest powszechny w workflow'ach bioinformatycznych, gdzie możesz chcieć przetwarzać powiązane próbki osobno ze względów wydajnościowych, zanim porównasz lub połączysz wyniki na końcu.

Nextflow zawiera wbudowane metody do tego celu — główną, którą omówimy, jest `groupTuple`.

Zacznijmy od grupowania wszystkich próbek, które mają te same pola `id` i `interval`. Byłoby to typowe dla analizy, w której chcemy grupować repliki techniczne, ale zachować osobno próbki znacząco różniące się od siebie.

W tym celu powinniśmy wyodrębnić nasze zmienne grupowania, abyśmy mogli używać ich w izolacji.

Pierwszy krok jest podobny do tego, co robiliśmy w poprzedniej sekcji. Musimy wyizolować naszą zmienną grupowania jako pierwszy element krotki. Pamiętaj, że nasz pierwszy element to obecnie mapa pól `id`, `repeat` i `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Możemy ponownie użyć metody `subMap` do wyizolowania pól `id` i `interval` z mapy. Podobnie jak poprzednio, użyjemy operatora `map` do zastosowania metody `subMap` do pierwszego elementu krotki dla każdej próbki.

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

Widzimy, że pomyślnie wyizolowaliśmy pola `id` i `interval`, ale próbki nie są jeszcze zgrupowane.

!!! note "Uwaga"

    Odrzucamy tu pole `replicate`. Wynika to z tego, że nie jest ono potrzebne do dalszego przetwarzania. Po ukończeniu tego samouczka sprawdź, czy możesz je uwzględnić bez wpływu na późniejsze grupowanie!

Zgrupujmy teraz próbki według tego nowego elementu grupowania, używając [operatora `groupTuple`](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

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

To wszystko! Dodaliśmy tylko jedną linię kodu. Zobaczmy, co się stanie po uruchomieniu:

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

Zwróć uwagę, że struktura danych uległa zmianie — pliki w każdym elemencie kanału są teraz zawarte w krotkach, takich jak `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. Wynika to z tego, że gdy używamy `groupTuple`, Nextflow łączy pojedyncze pliki dla każdej próbki z grupy. Warto o tym pamiętać przy obsłudze danych w dalszych etapach.

!!! note "Uwaga"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) jest odwrotnością `groupTuple`. Rozpakowuje elementy w kanale i spłaszcza je. Spróbuj dodać `transpose` i cofnąć grupowanie, które wykonaliśmy powyżej!

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Grupowania powiązanych próbek**: Jak używać `groupTuple` do agregowania próbek według wspólnych atrybutów
- **Izolowania kluczy grupowania**: Jak używać `subMap` do wyodrębniania określonych pól do grupowania
- **Obsługi zgrupowanych struktur danych**: Jak pracować z zagnieżdżoną strukturą tworzoną przez `groupTuple`
- **Obsługi replik technicznych**: Jak grupować próbki, które mają te same warunki eksperymentalne

---

## Podsumowanie

W tym zadaniu pobocznym nauczyłeś się dzielić i grupować dane przy użyciu kanałów.

Modyfikując dane w miarę ich przepływu przez pipeline, możesz budować skalowalny pipeline bez używania pętli ani instrukcji while, co oferuje kilka zalet w porównaniu z bardziej tradycyjnymi podejściami:

- Możemy skalować do dowolnej liczby wejść bez dodatkowego kodu
- Skupiamy się na obsłudze przepływu danych przez pipeline, zamiast na iteracji
- Możemy być tak złożeni lub prości, jak wymaga tego sytuacja
- Pipeline staje się bardziej deklaratywny, skupiając się na tym, co powinno się wydarzyć, a nie jak
- Nextflow optymalizuje wykonanie za nas, uruchamiając niezależne operacje równolegle

Opanowanie tych operacji na kanałach pozwoli Ci budować elastyczne, skalowalne pipeline'y obsługujące złożone relacje danych bez uciekania się do pętli lub programowania iteracyjnego, pozwalając Nextflow'owi automatycznie optymalizować wykonanie i paralelizować niezależne operacje.

### Kluczowe wzorce

1.  **Tworzenie ustrukturyzowanych danych wejściowych:** Zaczynając od pliku CSV z mapami metadanych (bazując na wzorcach z [Metadata in workflows](../metadata/))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Podział danych na osobne kanały:** Użyliśmy `filter` do podziału danych na niezależne strumienie na podstawie pola `type`

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Łączenie dopasowanych próbek:** Użyliśmy `join` do ponownego łączenia powiązanych próbek na podstawie pól `id` i `repeat`

    - Łączenie dwóch kanałów według klucza (pierwszego elementu krotki)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Wyodrębnienie klucza łączenia i łączenie według tej wartości

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Łączenie na wielu polach przy użyciu subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Rozdzielanie na interwały:** Użyliśmy `combine` do tworzenia iloczynów kartezjańskich próbek z interwałami genomicznymi do przetwarzania równoległego.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Agregowanie według kluczy grupowania:** Użyliśmy `groupTuple` do grupowania według pierwszego elementu każdej krotki, zbierając tym samym próbki współdzielące pola `id` i `interval` oraz łącząc repliki techniczne.

    ```groovy
    channel.groupTuple()
    ```

6.  **Optymalizacja struktury danych:** Użyliśmy `subMap` do wyodrębniania określonych pól i utworzyliśmy nazwane domknięcie, aby transformacje były wielokrotnego użytku.

    - Wyodrębnianie określonych pól z mapy

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Użycie nazwanego domknięcia do transformacji wielokrotnego użytku

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

Wróć do [menu zadań pobocznych](../) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
