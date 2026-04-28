# Metadane i mapy metadanych

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W każdej analizie naukowej rzadko pracujemy wyłącznie z surowymi plikami danych.
Każdy plik niesie ze sobą dodatkowe informacje: czym jest, skąd pochodzi i co go wyróżnia.
Te dodatkowe informacje nazywamy metadanymi.

Metadane to dane opisujące inne dane.
Śledzą ważne szczegóły dotyczące plików i warunków eksperymentalnych, a także pomagają dostosować analizy do unikalnych cech każdego zestawu danych.

Wyobraź sobie katalog biblioteczny: podczas gdy książki zawierają właściwą treść (surowe dane), karty katalogowe dostarczają kluczowych informacji o każdej książce — kiedy została wydana, kto ją napisał, gdzie ją znaleźć (metadane).
W pipeline'ach Nextflow metadane można wykorzystać do:

- Śledzenia informacji specyficznych dla pliku w trakcie całego workflow'u
- Konfigurowania procesów na podstawie właściwości pliku
- Grupowania powiązanych plików do wspólnej analizy

### Cele szkolenia

W tym zadaniu pobocznym zbadamy, jak obsługiwać metadane w workflow'ach.
Zaczynając od prostego arkusza danych (w bioinformatyce często nazywanego samplesheet), zawierającego podstawowe informacje o plikach, nauczysz się:

- Odczytywać i parsować metadane plików z plików CSV
- Rozumieć, dlaczego interfejs „mapa meta + plik danych" jest powszechnie stosowaną konwencją
- Dodawać nowe pola metadanych podczas wykonywania workflow'u
- Używać metadanych do dostosowywania zachowania procesów i organizowania wyników

Te umiejętności pomogą Ci budować bardziej niezawodne i elastyczne pipeline'y, zdolne do obsługi złożonych relacji między plikami i wymagań przetwarzania.

### Wymagania wstępne

Przed przystąpieniem do tego zadania pobocznego powinieneś/powinnaś:

- Ukończyć samouczek [Hello Nextflow](../../hello_nextflow/index.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow (procesy, kanały, operatory).

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś/zrobiłaś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/metadata
```

Możesz ustawić VSCode tak, aby skupiał się na tym katalogu:

```bash
code .
```

Edytor otworzy się z widokiem na katalog projektu.

#### Przejrzyj materiały

Znajdziesz tu główny plik workflow'u oraz katalog `data` zawierający arkusz danych i kilka plików z danymi.

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

Workflow w pliku `main.nf` to szkielet, który będziesz stopniowo rozbudowywać w pełni działający workflow.

Arkusz danych zawiera ścieżki do plików z danymi oraz powiązane metadane, zorganizowane w 3 kolumnach:

- `id`: oczywiste — identyfikator nadany plikowi
- `character`: nazwa postaci, której użyjemy później do rysowania różnych stworzeń
- `data`: ścieżki do plików `.txt` zawierających pozdrowienia w różnych językach

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Każdy plik z danymi zawiera tekst pozdrowienia w jednym z pięciu języków (fr: francuski, de: niemiecki, es: hiszpański, it: włoski, en: angielski).

Użyjemy narzędzia o nazwie [`COWPY`](https://github.com/jeffbuttars/cowpy) do generowania ASCII art każdej postaci wypowiadającej nagrane pozdrowienie.

??? info "Co robi `COWPY`?"

    `COWPY` to narzędzie wiersza poleceń generujące ASCII art do wyświetlania dowolnych tekstów w zabawny sposób.
    Jest to implementacja w Pythonie klasycznego narzędzia [cowsay](https://en.wikipedia.org/wiki/Cowsay) autorstwa Tony'ego Monroe'a.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Opcjonalnie możesz wybrać postać (lub 'cowacter') zamiast domyślnej krowy.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Dodatkowo użyjemy narzędzia do analizy języka o nazwie `langid`, aby zidentyfikować język każdej postaci i odpowiednio zorganizować wyniki pipeline'u.

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest napisanie workflow'u Nextflow, który:

1. **Wygeneruje ASCII art** każdej postaci
2. **Pogrupuje** wyniki według rodziny językowej (języki germańskie vs romańskie)

To typowy wzorzec workflow'u, w którym metadane specyficzne dla pliku sterują decyzjami przetwarzania — dokładnie ten rodzaj problemu, który mapy metadanych rozwiązują elegancko.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem/ustawiłam odpowiedni katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Podstawowe sposoby wczytywania i używania metadanych

Otwórz plik workflow'u `main.nf`, aby przejrzeć szkielet, który dajemy Ci jako punkt startowy.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

Operator [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) wczytuje każdy wiersz pliku jako element kanału.
To samo podejście stosujemy do wczytywania danych CSV w Hello Nextflow, naszym kursie dla początkujących.
Zajrzyj do [tej sekcji](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file), jeśli potrzebujesz przypomnienia, jak to działa.

Dzięki opcji `header: true` pierwszy wiersz jest traktowany jako nagłówki kolumn, więc każdy element staje się mapą par klucz-wartość z kluczami odpowiadającymi nazwom kolumn.

Ponieważ nie uruchamiamy jeszcze żadnych procesów na danych, bloki `publish` i `output` są na razie tylko szkieletem.

### 1.1. Uruchom workflow

Uruchom workflow, aby zobaczyć, jak jest zorganizowana zawartość kanału po wczytaniu wszystkich danych:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Widać, że operator zbudował mapę par klucz-wartość dla każdego wiersza pliku CSV, używając nagłówków kolumn jako kluczy dla odpowiadających im wartości.

Każdy wpis mapy odpowiada kolumnie w naszym arkuszu danych:

- `id`
- `character`
- `recording`

Dzięki temu łatwo jest uzyskać dostęp do konkretnych pól każdego wiersza.
Na przykład możemy uzyskać dostęp do identyfikatora pliku przez `id` lub do ścieżki pliku txt przez `recording`.

??? info "(Opcjonalnie) Więcej o mapach Groovy"

    W Groovy, języku programowania, na którym zbudowany jest Nextflow, mapa to struktura danych klucz-wartość podobna do słowników w Pythonie, obiektów w JavaScript czy haszy w Ruby.

    Oto uruchamialny skrypt pokazujący, jak w praktyce zdefiniować mapę i uzyskać dostęp do jej zawartości:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Utwórz prostą mapę
    def my_map = [id:'sampleA', character:'squirrel']

    // Wydrukuj całą mapę
    println "map: ${my_map}"

    // Uzyskaj dostęp do poszczególnych wartości za pomocą notacji kropkowej
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Mimo że nie ma właściwego bloku `workflow`, Nextflow może uruchomić ten skrypt jak workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Oto czego możesz się spodziewać w wynikach:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Wybierz konkretne pole za pomocą `map`

Użyjemy operatora `map`, aby iterować po każdym elemencie kanału i wybrać konkretnie pole `character`, do którego możemy uzyskać dostęp po nazwie za pomocą notacji kropkowej.

#### 1.2.1. Dodaj operację map

Aby uzyskać dostęp do kolumny `character`, dodaj operację `map` przed operacją `.view()` w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Ten sposób dostępu do konkretnego pola jest szczegółowo wyjaśniony w [tej sekcji](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) kursu Hello Nextflow, jeśli potrzebujesz przypomnienia.

#### 1.2.2. Uruchom workflow

Uruchom workflow, aby sprawdzić, czy możesz wyświetlić wyodrębnione nazwy postaci.

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Potwierdza to, że możemy uzyskać dostęp do wartości z kolumny `character` dla każdego wiersza.

Teraz zróbmy coś z tymi danymi: użyjmy pól `character` i `recording` razem, aby wygenerować ASCII art za pomocą `COWPY`.

### 1.3. Emituj podkanały za pomocą `multiMap`

Dostarczamy Ci gotowy moduł procesu `COWPY`, więc najpierw musisz zapoznać się z wymaganiami wejściowymi tego procesu.

Możesz otworzyć plik, aby zobaczyć, jak wygląda proces:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Generuj ASCII art za pomocą cowpy
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Jak widać, proces przyjmuje dwa oddzielne wejścia: plik nagrania i nazwę postaci.
Co ważne, mamy wartości dla obu, ale są one obecnie spakowane razem w każdym elemencie kanału.

Jednym ze sposobów wyodrębnienia wielu pól do oddzielnych kanałów jest operator [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap), który dzieli jeden kanał na wiele nazwanych podkanałów w jednej operacji.

#### 1.3.1. Dodaj operację multiMap

Zastąp operację `map` przez `multiMap`:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

Blok `multiMap` definiuje dwa nazwane podkanały (`file` i `character`) dla każdego wiersza, do których możemy uzyskać dostęp jako `ch_datasheet.file` i `ch_datasheet.character`.

#### 1.3.2. Wywołaj COWPY na podkanałach

Teraz dołącz moduł procesu `COWPY` i przekaż każdy podkanał jako oddzielny argument:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Pozwala nam to przekazać oba pola oddzielnie, zgodnie z wymaganiami `COWPY`.

#### 1.3.3. Skonfiguruj publikowanie wyników

Na koniec dodaj wyniki `COWPY` do bloku `publish:`:

=== "Po"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Dzięki temu będziemy mogli łatwo przeglądać wyniki produkowane przez workflow.

#### 1.3.4. Uruchom workflow

Uruchom workflow, aby sprawdzić, czy `COWPY` działa na dostarczonych danych wejściowych:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Jak widać, `COWPY` uruchomił się dla każdego pliku, używając właściwej postaci.

??? abstract "Zawartość katalogu wyników"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Zawartość pliku results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

To podejście działa, ale ma pewne ograniczenie: musieliśmy podzielić kanał na dwa oddzielne podkanały.
Gdybyśmy chcieli przekazać więcej pól do procesu, musielibyśmy wyodrębniać kolejne podkanały — co mogłoby stać się uciążliwe i nieczytelne.

Dobra wiadomość: istnieje prostszy sposób.

### 1.4. Zgrupuj wszystko jako pojedyncze wejście do procesu

Zamiast dzielić pola na oddzielne kanały, możemy zaktualizować proces tak, aby przyjmował wszystkie wejścia jako pojedynczą krotkę, co upraszcza wywołanie procesu.

#### 1.4.1. Zaktualizuj proces COWPY

Zaktualizuj `COWPY`, aby przyjmował krotkę odpowiadającą trzem elementom każdego wiersza:

=== "Po"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Generuj ASCII art za pomocą cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Przed"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Generuj ASCII art za pomocą cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
        """
    }
    ```

Teraz proces przyjmuje jedno wejście zawierające wszystko, co możemy chcieć mu przekazać.

#### 1.4.2. Użyj `map()` do utworzenia krotki wejściowej

Nadal musimy użyć operacji mapowania, aby wyliczyć elementy, które chcemy przekazać w krotce do procesu:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Możesz się zastanawiać, dlaczego nie możemy po prostu przekazać całej mapy Groovy zwracanej przez `splitCsv` bez zmian.
Powodem jest to, że musimy jawnie poinformować Nextflow, że plik nagrania powinien być traktowany jako ścieżka (tzn. musi być odpowiednio przygotowany do użycia).
Dzieje się to na poziomie interfejsu wejściowego `COWPY`, gdzie element `recording` jest jawnie oznaczony jako `path`.

#### 1.4.3. Zaktualizuj wywołanie procesu

Na koniec zastąpmy dwa oddzielne wejścia w wywołaniu procesu pojedynczą krotką, którą właśnie utworzyliśmy:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Nieco upraszcza to wywołanie procesu.

#### 1.4.4. Uruchom workflow

Uruchom workflow, aby sprawdzić, czy `COWPY` nadal poprawnie przetwarza dane:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

Wynikiem jest te same siedem plików `cowpy-*.txt` co poprzednio, teraz produkowanych przy prostszym wywołaniu `COWPY`.

??? abstract "Zawartość katalogu wyników"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Zawartość pliku results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

To nieznaczna poprawa w stosunku do podejścia z `multiMap`.
Nadal jednak musieliśmy rozpakować oryginalną mapę Groovy, aby utworzyć krotkę wejściową, a między procesem a arkuszem danych istnieje ścisłe powiązanie: definicja wejścia `COWPY` bezpośrednio odwołuje się do nazw kolumn `id`, `character` i `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Jeśli współpracownik używa inaczej zorganizowanego arkusza danych — z dodatkowymi kolumnami lub kolumnami w innej kolejności — ten proces nie zadziała bez modyfikacji.
Sprawia to, że jest on kruchy, ponieważ jego struktura wejściowa jest ściśle powiązana z dokładnym składem arkusza danych.

Aby rozwiązać ten problem, potrzebujemy sposobu na przekazanie wszystkich metadanych jako pakietu bez zakodowywania ich dokładnej struktury w interfejsie procesu.

### 1.5. Użyj interfejsu mapa meta + plik

Rozwiązaniem jest rozdzielenie dwóch odrębnych zagadnień w kanale: **metadanych o próbce** i samego **pliku danych**.
Zbierając wszystkie metadane w jedną mapę — „mapę meta" — uzyskujemy spójną krotkę dwuelementową niezależnie od liczby kolumn metadanych w arkuszu danych:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Dodanie lub usunięcie kolumn z arkusza danych zmienia zawartość `meta`, ale kształt krotki `[meta, file]` pozostaje stały.
Procesy przyjmujące tę strukturę nie muszą wiedzieć ani dbać o to, ile pól metadanych istnieje.

#### 1.5.1. Przeorganizuj zawartość krotki w mapę meta

Przestrukturyzujmy operację `map`, aby produkowała krotkę `[meta, file]`:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Zaktualizujemy w następnym kroku

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Zwróć uwagę, że dodaliśmy też instrukcję `view()`, zakomentowaliśmy wywołanie `COWPY` i zastąpiliśmy `COWPY.out` przez `channel.empty()`, ponieważ definicja wejścia procesu nie pasuje jeszcze do nowej struktury.

#### 1.5.2. Uruchom workflow, aby sprawdzić przeorganizowaną zawartość

Uruchom workflow, aby zobaczyć nowy kształt kanału:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Każdy element kanału jest teraz krotką dwuelementową: najpierw mapa meta, potem plik.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Jeśli później dodamy kolumnę `language` do arkusza danych, będzie ona dostępna jako `meta.language` bez konieczności wprowadzania jakichkolwiek zmian w definicji wejścia procesu.

#### 1.5.3. Zaktualizuj proces `COWPY`, aby używał mapy meta

Zaktualizuj `COWPY`, aby przyjmował strukturę krotki `[meta, file]`:

=== "Po"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Generuj ASCII art za pomocą cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Przed"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Generuj ASCII art za pomocą cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Wewnątrz bloku skryptu `meta.character` uzyskuje dostęp do pola `character` z mapy meta.
Każde pole w mapie meta jest dostępne w ten sam sposób.

#### 1.5.4. Zaktualizuj wywołanie procesu

Przywróć wywołanie `COWPY` i podłącz jego wyniki do publikowania:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Zaktualizujemy w następnym kroku

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

Przywróciliśmy też publikowanie wyników.

#### 1.5.5. Uruchom workflow

Uruchom workflow, aby sprawdzić, czy wszystko działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

Katalog wyników zawiera teraz pliki z ASCII art.

??? abstract "Zawartość katalogu"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

??? example "Zawartość pliku results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
       \
        \
            .--.
           |o_o |
           |:_/ |
          //   \ \
         (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Proces otrzymuje teraz wszystkie metadane jako pakiet przez `meta`, używa tego, czego potrzebuje (`meta.character`), i ignoruje resztę.

Jest to standardowy interfejs używany przez wszystkie moduły [nf-core](https://nf-co.re/).
Wzorzec `tuple val(meta), path(file)` pojawia się konsekwentnie w całej bibliotece modułów nf-core, dlatego workflow'y przyjmujące tę konwencję mogą bez trudu korzystać z modułów nf-core.

### Podsumowanie

W tej sekcji nauczyłeś/nauczyłaś się:

- **Jak wczytywać arkusze danych:** Używanie `splitCsv` do parsowania plików CSV z informacjami nagłówkowymi
- **Dlaczego istnieje konwencja mapy meta:** Oddzielanie metadanych od plików danych w krotkach `[meta, file]` utrzymuje stabilną strukturę kanału w miarę ewolucji arkusza danych
- **Jak używać pól mapy meta wewnątrz procesu:** Każde pole w mapie meta jest dostępne za pomocą notacji kropkowej w bloku skryptu

---

## 2. Dodatkowe manipulacje metadanymi

Teraz, gdy interfejs mapy meta jest już na miejscu, możemy wzbogacać ją w miarę przepływu danych przez pipeline.

Użyjemy narzędzia o nazwie [`langid`](https://github.com/saffsd/langid.py) do identyfikacji języka w każdym pliku nagrania.
Dla podanego fragmentu tekstu narzędzie wypisuje na `stdout` przewidywanie języka i wynik prawdopodobieństwa.

### 2.1. Dodaj krok identyfikacji języka

Dostarczamy Ci gotowy moduł procesu o nazwie `IDENTIFY_LANGUAGE`, który opakowuje narzędzie `langid`.

Możesz otworzyć plik modułu, aby przejrzeć jego kod:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
// Użyj langid do przewidzenia języka każdego pliku wejściowego
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Definicja wejścia używa tej samej struktury `tuple val(meta), path(file)`, którą zbudowaliśmy w sekcji 1, więc `ch_datasheet` może bezpośrednio zasilać ten proces bez żadnych adaptacji.

Wyjście dodaje `stdout` jako trzeci element: przechwytuje przewidywanie języka, które `langid` wypisuje na konsolę.
Polecenie `sed` usuwa wynik prawdopodobieństwa i końcowy znak nowej linii, pozostawiając tylko dwuliterowy kod języka.

#### 2.1.1. Dodaj wywołanie `IDENTIFY_LANGUAGE`

Dołącz moduł procesu `IDENTIFY_LANGUAGE` i wywołaj go na kanale arkusza danych:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

Głównym wynikiem tego procesu jest tylko string, więc nie ma plików wyjściowych do opublikowania.
Zamiast tego używamy `IDENTIFY_LANGUAGE.out.view()`, aby podejrzeć wyniki operacji.

#### 2.1.2. Uruchom workflow

Uruchom workflow, aby uzyskać identyfikację języka, używając `-resume`, aby uniknąć ponownego uruchamiania zadań `COWPY`:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Mamy teraz przewidywanie języka dla każdego pliku w zestawie danych.

Zwróć uwagę, że krotka wyjściowa składa się z `[meta, file, lang_id]`, co oznacza, że mapa meta i plik są przenoszone razem z nowym wynikiem.

!!! note "Uwaga"

    Ten wzorzec utrzymywania mapy meta powiązanej z wynikami ułatwia późniejsze kojarzenie wyników między kanałami.
    Nie można polegać na kolejności elementów w kanałach przy prawidłowym kojarzeniu danych.
    Zamiast tego należy używać kluczy.
    Mapy meta zapewniają do tego idealną strukturę.

    Ten przypadek użycia szczegółowo omawiamy w zadaniu pobocznym [Splitting & Grouping](../splitting_and_grouping/index.md).

### 2.2. Wzbogać metadane o wyniki procesu

Przewidywanie języka jest samo w sobie metadaną o danych w pliku.
Zamiast przechowywać je jako oddzielny element, włączmy je z powrotem do mapy meta.

#### 2.2.1. Utwórz nową, rozszerzoną mapę meta

Możemy utworzyć nową mapę meta zastępującą oryginalną za pomocą operatora Groovy `+`:

=== "Po"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Sercem tej operacji jest `#!groovy meta + [lang: lang_id]`.

Ten kod tworzy tymczasową mapę z jedną parą klucz-wartość zawierającą kod języka (`[lang: lang_id]`), a następnie używa operatora Groovy `+`, aby połączyć ją z oryginalną mapą `meta` zawierającą istniejące metadane, tworząc nową, rozszerzoną mapę meta.

Szczegółowe wyjaśnienie znajdziesz w poniższym polu.

??? info "Tworzenie nowej mapy meta za pomocą operatora `+`"

    **Po pierwsze, musisz wiedzieć, że możemy scalić zawartość dwóch map za pomocą operatora Groovy `+`.**

    Powiedzmy, że mamy następujące mapy:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Możemy je scalić w ten sposób:

    ```groovy
    new_map = map1 + map2
    ```

    Zawartość `new_map` będzie następująca:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Świetnie!

    **Ale co, jeśli chcesz dodać pole, które nie jest jeszcze częścią żadnej mapy?**

    Powiedzmy, że zaczynasz ponownie od `map1`, ale przewidywanie języka nie jest w swojej własnej mapie (nie ma `map2`).
    Zamiast tego jest przechowywane w zmiennej o nazwie `lang_id`, a Ty wiesz, że chcesz zapisać jej wartość (`'fr'`) pod kluczem `lang`.

    Możesz zrobić następująco:

    ```groovy
    new_map = map1 + [lang: lang_id]
    ```

    Tutaj `[lang: lang_id]` tworzy nową anonimową mapę w locie, a `map1 + ` scala `map1` z nową anonimową mapą, dając tę samą zawartość `new_map` co poprzednio.

    Eleganckie, prawda?

    **Teraz przełóżmy to na kontekst operacji `channel.map()` w Nextflow.**

    Kod przyjmuje postać:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Robi to następujące rzeczy:

    - `#!groovy map1, lang_id ->` pobiera dwa elementy z krotki
    - `#!groovy map1 + [lang: lang_id]` tworzy nową mapę zgodnie z opisem powyżej

    Wynikiem jest pojedyncza anonimowa mapa o tej samej zawartości co `new_map` w naszym przykładzie.
    Efektywnie przekształciliśmy:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    w:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Mam nadzieję, że widzisz, że jeśli zmienimy `map1` na `meta`, to właściwie wszystko, czego potrzebujemy, aby dodać przewidywanie języka do naszej mapy meta w workflow.

    Jest jednak jedna rzecz do uwzględnienia!

    W przypadku naszego workflow'u **musimy też wziąć pod uwagę obecność obiektu `file` w krotce**, która składa się z `meta, file, lang_id`.

    Kod przyjmuje więc postać:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Jeśli trudno Ci śledzić, dlaczego `file` wydaje się przemieszczać w operacji `map`, wyobraź sobie, że zamiast `#!groovy [meta + [lang: lang_id], file]` ten wiersz brzmi `[new_map, file]`.
    Powinno to wyjaśnić, że po prostu zostawiamy `file` na jego oryginalnym miejscu na drugiej pozycji w krotce. Wzięliśmy wartość `new_info` i wbudowaliśmy ją w mapę na pierwszej pozycji.

    **I to prowadzi nas z powrotem do struktury kanału `tuple val(meta), path(file)`!**

#### 2.2.2. Uruchom workflow

Gdy już rozumiesz, co robi ten kod, uruchom workflow, aby sprawdzić, czy zadziałał:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Tak, to się zgadza!
Starannie przeorganizowaliśmy wynik procesu z `meta, file, lang_id` tak, że `lang_id` jest teraz jednym z kluczy w mapie meta, a krotki kanału ponownie pasują do modelu `meta, file`.

!!! tip "Wskazówka: usuwanie kluczy z mapy meta"

    Możesz usunąć klucz z mapy meta za pomocą metody Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)), która zwraca nową mapę zawierającą tylko wskazane klucze:

    ```groovy
    meta.subMap(['id', 'character'])  // zwraca mapę zawierającą tylko 'id' i 'character'
    ```

    Jest to przydatne, gdy dalszy proces lub moduł nie potrzebuje wszystkich pól, które nagromadziły się w mapie meta.

### 2.3. Przypisz grupę językową za pomocą instrukcji warunkowych

Mając przewidywanie języka w mapie meta, możemy wyprowadzić z niego kolejne metadane.
Języki w naszym zestawie danych należą do dwóch rodzin: germańskiej (angielski, niemiecki) i romańskiej (francuski, hiszpański, włoski).
Dodanie pola `lang_group` udostępni tę klasyfikację w dalszej części pipeline'u.

#### 2.3.1. Dodaj operację `map` z logiką warunkową

Użyjemy drugiej operacji `map` z logiką warunkową, aby przypisać rodzinę językową:

```groovy
.map { meta, file ->

    // tutaj umieszczamy logikę warunkową definiującą lang_group

    [meta + [lang_group: lang_group], file]
}
```

Oto logika, którą chcemy zastosować:

- Zacznij od `lang_group = 'unknown'` jako wartości domyślnej.
- Jeśli `meta.lang` to `'de'` lub `'en'`, ustaw `lang_group` na `'germanic'`.
- W przeciwnym razie, jeśli `meta.lang` należy do `['fr', 'es', 'it']`, ustaw `lang_group` na `'romance'`.

!!! tip "Wskazówka"

    Dostęp do wartości `lang` w operacji map uzyskasz przez `meta.lang`.

Wprowadź następujące zmiany w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        ch_languages.view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Kluczowe punkty:

- `def lang_group = "unknown"` inicjalizuje zmienną z bezpieczną wartością domyślną.
- Struktura `if / else if` obsługuje dwie rodziny językowe; wszystko inne pozostaje jako `'unknown'`.
- `#!groovy .set { ch_languages }` nadaje wynikowemu kanałowi nazwę do użycia w następnym kroku.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Uruchom workflow

Uruchom workflow, aby sprawdzić, czy działa poprawnie:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Mapa meta zawiera teraz cztery pola: `id`, `character`, `lang` i `lang_group`.
Struktura kanału nadal wynosi `[meta, file]`.

### 2.4. Używaj metadanych do nazywania i organizowania wyników

Mając teraz `lang` i `lang_group` w mapie meta, możemy ich użyć do dodania kodu języka do nazw plików wyjściowych i zorganizowania ich w podkatalogi według rodziny językowej.

Wymaga to trzech zmian: zaktualizowania procesu `COWPY`, aby zmieniał nazwy wyników i uwzględniał `meta` w tym, co emituje; zaktualizowania wywołania `COWPY`, aby działał na `ch_languages`; oraz zaktualizowania bloku output, aby określał ścieżkę podkatalogu.

#### 2.4.1. Zaktualizuj proces `COWPY`

Zmień nazwę pliku wyjściowego, używając kodu języka z mapy meta, i dodaj `meta` do wyjścia, aby blok output mógł uzyskać dostęp do `lang_group` przy kierowaniu do podkatalogu:

=== "Po"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Przed"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Pokazuje to, jak możemy korzystać z innych pól metadanych do dostosowywania zachowania procesu bez konieczności modyfikowania definicji wejścia.

#### 2.4.2. Zaktualizuj wywołanie `COWPY`, aby działał na `ch_languages`

Zastąp `COWPY(ch_datasheet)` przez `COWPY(ch_languages)`:

=== "Po"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

Usuwamy też wiersz `ch_languages.view()`, ponieważ nie musimy już podglądać zawartości kanału.

#### 2.4.3. Zaktualizuj blok output

Dodaj domknięcie `path` do bloku `output {}`, aby kierować każdy plik do podkatalogu odpowiadającego jego grupie językowej:

=== "Po"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Pokazuje to, jak możemy używać metadanych do elastycznego organizowania wyników.

#### 2.4.4. Uruchom pełny pipeline

Usuń poprzednie wyniki i uruchom pełny pipeline:

```bash
rm -r results
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

Katalog wyników jest teraz zorganizowany według rodziny językowej, a każdy plik nosi nazwę wykrytego języka:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

Domknięcie `path` w bloku `output {}` otrzymuje każdą krotkę `[meta, file]` i zwraca `meta.lang_group` jako nazwę podkatalogu.
Nazwa samego pliku pochodzi z tego, co proces wypisuje (`#!groovy "${meta.lang}-${input_file}"`).
Oba elementy metadanych (kod języka i grupa językowa) pochodzą z wzbogaconej mapy meta zbudowanej w tej sekcji.

### Podsumowanie

W tej sekcji nauczyłeś/nauczyłaś się:

- **Jak wzbogacać mapę meta o wyniki procesu:** Dodawanie nowych kluczy za pomocą `#!groovy meta + [klucz: wartość]` zachowuje strukturę kanału `[meta, file]` przy jednoczesnym wzbogacaniu metadanych.
- **Jak wyprowadzać metadane z metadanych:** Logika warunkowa wewnątrz operacji `map` może obliczać nowe pola na podstawie istniejących.
- **Jak używać metadanych do organizowania wyników:** Domknięcie `path` w bloku `output {}` może odczytywać z mapy meta, aby kierować pliki do podkatalogów.

---

## 3. Kwestie niezawodności

Gdy wartości metadanych sterują zachowaniem procesu, brakujące lub niekompletne dane mogą powodować trudne do zdiagnozowania problemy.
Oto czego się spodziewać i jak sobie z tym radzić.

### 3.1. Co się dzieje, gdy brakuje wymaganego pola metadanych

Wartość `character` jest wymagana, aby proces `COWPY` produkował poprawne wyniki.
Sposób, w jaki objawia się błąd, zależy od tego, czy kolumna istnieje w arkuszu danych, ale jest pusta, czy też w ogóle jej nie ma.

#### 3.1.1. Kolumna istnieje, ale wartość jest pusta

Powiedzmy, że jeden wpis w arkuszu danych ma puste pole `character`:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Klucz `character` jest tworzony dla wszystkich wpisów podczas parsowania arkusza danych, ale `meta.character` dla `sampleA` będzie pustym stringiem.
Gdy Nextflow podstawia `#!groovy ${meta.character}` do polecenia, narzędzie `COWPY` otrzymuje pusty argument dla `-c` i zgłasza błąd:

??? failure "Wyjście polecenia"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Komunikat błędu (`expected one argument`) wskazuje na pusty argument `-c`.
Sprawdzenie pliku `.command.sh` w katalogu roboczym potwierdza, że polecenie zostało uruchomione z pustą wartością.

#### 3.1.2. Kolumna nie istnieje w arkuszu danych

Jeśli kolumna `character` w ogóle nie istnieje:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Klucz `character` nigdy nie zostanie utworzony w mapie meta.
Gdy skrypt procesu oblicza `#!groovy ${meta.character}`, brakujący klucz zwraca `null`, a Nextflow dosłownie podstawia string `null` do polecenia:

??? failure "Wyjście polecenia"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Wskazówką diagnostyczną jest `cowpy -c null` w wykonanym poleceniu.

### 3.2. Strategie obsługi brakujących metadanych

Istnieją dwa uzupełniające się podejścia, które czynią workflow'y bardziej odpornymi na brakujące metadane.

**1. Walidacja wejść**

Najbardziej niezawodnym rozwiązaniem jest walidacja arkusza danych przed rozpoczęciem jakiegokolwiek przetwarzania, dzięki czemu problemy są wykrywane wcześnie z czytelnym komunikatem błędu, zamiast ujawniać się jako tajemnicze błędy procesu w trakcie uruchomienia.
Szkolenie [Hello nf-core](../../hello_nf-core/05_input_validation.md) omawia, jak dodać walidację wejść za pomocą wtyczki nf-schema. <!-- TODO (future) pending a proper Validation side quest -->

**2. Jawne wejścia procesu dla wymaganych wartości**

Jeśli chcesz, aby interfejs procesu sam komunikował, że dana wartość jest obowiązkowa, rozważ wyodrębnienie jej z mapy meta jako jawnego wejścia:

=== "Definicja procesu"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Wywołanie w workflow'ie"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

To podejście sprawia, że `character` jest widoczną, wymaganą częścią kontraktu procesu.
Każdy, kto czyta moduł, od razu widzi, że wartość postaci musi być dostarczona.
Jeśli pole jest nieobecne, workflow kończy się błędem wyraźnie na poziomie kanału, zanim proces w ogóle zostanie uruchomiony.

Podkreśla to użyteczną zasadę projektowania:

**Używaj mapy meta dla opcjonalnych lub opisowych informacji; wyodrębniaj wymagane wartości jako jawne wejścia.**

Mapa meta utrzymuje czyste i stabilne struktury kanałów, ale dla wartości, które są naprawdę wymagane przez proces, ujawnienie ich jako nazwanych wejść poprawia czytelność i ułatwia prawidłowe używanie modułu w innych kontekstach.

### Podsumowanie

W tej sekcji zobaczyłeś/zobaczyłaś:

- **Jak objawia się brak metadanych:** Puste pole powoduje pusty argument; brakujące pole powoduje dosłowne podstawienie `null` do polecenia.
- **Dwie uzupełniające się strategie:** Walidacja wejść do wczesnego wykrywania problemów oraz jawne wejścia procesu do czytelnego komunikowania wymagań.

---

## Podsumowanie

W tym zadaniu pobocznym zbadałeś/zbadałaś, jak efektywnie pracować z metadanymi w workflow'ach Nextflow.

Wzorzec krotki „mapa meta + plik danych" jest podstawową konwencją w Nextflow, oferującą kilka zalet w porównaniu z przekazywaniem metadanych jako pojedynczych wartości:

- Struktura kanału pozostaje stabilna w miarę ewolucji arkusza danych
- Zachowanie procesu można dostosować dla każdej próbki bez zakodowywania nazw pól na stałe
- Metadane są dostępne przez cały pipeline do nazywania, grupowania i organizowania wyników
- Moduły napisane zgodnie z tym interfejsem są wymienne, w tym moduły nf-core

### Kluczowe wzorce

1.  **Odczytywanie i strukturyzowanie metadanych:** Parsowanie arkusza danych CSV i tworzenie mapy meta.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Rozszerzanie metadanych podczas workflow'u:** Dodawanie nowych kluczy z wyników procesu lub wyprowadzonej logiki.

    ```groovy
    // Z wyników procesu
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // Z logiki warunkowej
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Używanie metadanych wewnątrz procesu:** Dostęp do dowolnego pola za pomocą notacji kropkowej w bloku skryptu.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organizowanie wyników według wartości metadanych:** Użycie domknięcia `path` w bloku `output {}`.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Dodatkowe zasoby

- [operator map](https://www.nextflow.io/docs/latest/operator.html#map)
- [operator multiMap](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [kwalifikator wyjścia stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Co dalej?

Wróć do [menu zadań pobocznych](../index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
