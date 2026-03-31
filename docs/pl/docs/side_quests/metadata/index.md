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
- Tworzyć mapy metadanych i nimi manipulować
- Dodawać nowe pola metadanych podczas wykonywania workflow'u
- Używać metadanych do dostosowywania zachowania procesów

Te umiejętności pomogą Ci budować bardziej niezawodne i elastyczne pipeline'y, zdolne do obsługi złożonych relacji między plikami i wymagań przetwarzania.

### Wymagania wstępne

Przed przystąpieniem do tego zadania pobocznego powinieneś/powinnaś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow (procesy, kanały, operatory).

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś/zrobiłaś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../envsetup/index.md).

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

Udostępniamy Ci również skonteneryzowane narzędzie do analizy języka o nazwie `langid`.

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest napisanie workflow'u Nextflow, który:

1. **Zidentyfikuje** język w każdym pliku automatycznie
2. **Pogrupuje** pliki według rodziny językowej (języki germańskie vs romańskie)
3. **Dostosuje** przetwarzanie każdego pliku na podstawie jego języka i metadanych
4. **Zorganizuje** wyniki według grupy językowej

To typowy wzorzec workflow'u, w którym metadane specyficzne dla pliku sterują decyzjami przetwarzania — dokładnie ten rodzaj problemu, który mapy metadanych rozwiązują elegancko.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem/ustawiłam odpowiedni katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Wczytaj metadane z arkusza danych

Otwórz plik workflow'u `main.nf`, aby przejrzeć szkielet, który dajemy Ci jako punkt startowy.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Widać, że skonfigurowaliśmy podstawową fabrykę kanałów do wczytania przykładowego arkusza danych jako pliku, ale nie odczyta ona jeszcze zawartości tego pliku.
Zacznijmy od dodania tej funkcjonalności.

### 1.1. Wczytaj zawartość za pomocą `splitCsv`

Musimy wybrać operator, który odpowiednio sparsuje zawartość pliku przy minimalnym nakładzie pracy z naszej strony.
Ponieważ nasz arkusz danych jest w formacie CSV, do tego zadania nadaje się operator [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), który wczytuje każdy wiersz pliku jako element kanału.

Wprowadź następujące zmiany, aby dodać operację `splitCsv()` do kodu budującego kanał, a także operację `view()` do sprawdzenia, czy zawartość pliku jest poprawnie wczytywana do kanału.

=== "Po"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Zwróć uwagę, że używamy opcji `header: true`, aby poinformować Nextflow, że pierwszy wiersz pliku CSV ma być traktowany jako wiersz nagłówkowy.

Zobaczmy, co z tego wyjdzie.
Uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

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

Świetnie! Dzięki temu łatwo jest uzyskać dostęp do konkretnych pól każdego pliku.
Na przykład możemy uzyskać dostęp do identyfikatora pliku przez `id` lub do ścieżki pliku txt przez `recording`.

??? info "(Opcjonalnie) Więcej o mapach"

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
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Wybierz konkretne pola za pomocą `map`

Powiedzmy, że chcemy uzyskać dostęp do kolumny `character` z arkusza danych i ją wydrukować.
Możemy użyć operatora `map` Nextflow, aby iterować po każdym elemencie kanału i wybrać konkretnie wpis `character` z obiektu mapy.

Wprowadź następujące zmiany w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Sukces! Wykorzystaliśmy strukturę mapy wywodzącej się z naszego arkusza danych, aby uzyskać dostęp do wartości z poszczególnych kolumn dla każdego wiersza.

Teraz, gdy pomyślnie wczytaliśmy arkusz danych i mamy dostęp do danych w każdym wierszu, możemy zacząć implementować logikę naszego pipeline'u.

### 1.3. Zorganizuj metadane w 'mapę meta'

W obecnym stanie workflow'u pliki wejściowe (pod kluczem `recording`) i powiązane metadane (`id`, `character`) są na równorzędnej pozycji — jakby wszystko było w jednym dużym worku.
Praktyczną konsekwencją jest to, że każdy proces konsumujący ten kanał musiałby być skonfigurowany z uwzględnieniem tej struktury:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

To jest w porządku, dopóki liczba kolumn w arkuszu danych się nie zmienia.
Jednak dodanie nawet jednej kolumny do arkusza spowoduje, że kształt kanału przestanie odpowiadać temu, czego oczekuje proces, a workflow zacznie generować błędy.
Utrudnia to też udostępnianie procesu innym osobom, które mogą mieć nieco inne dane wejściowe, i może prowadzić do konieczności zakodowania na stałe zmiennych w procesie, które nie są potrzebne w bloku skryptu.

Aby uniknąć tego problemu, musimy znaleźć sposób na utrzymanie spójnej struktury kanału niezależnie od liczby kolumn w arkuszu danych.

Możemy to zrobić, zbierając wszystkie metadane w jeden element krotki, który nazwiemy mapą metadanych, lub krócej — 'mapą meta'.

Wprowadź następujące zmiany w operacji `map`:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Przestrukturyzowaliśmy elementy kanału w krotkę składającą się z dwóch elementów: mapy meta i odpowiadającego jej obiektu pliku.

Uruchommy workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Teraz każdy element kanału zawiera najpierw mapę meta, a następnie odpowiadający jej obiekt pliku:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

W rezultacie dodanie kolejnych kolumn do arkusza danych udostępni więcej metadanych w mapie `meta`, ale nie zmieni kształtu kanału.
Pozwala nam to pisać procesy konsumujące kanał bez konieczności zakodowywania na stałe elementów metadanych w specyfikacji wejścia:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Jest to powszechnie stosowana konwencja organizowania metadanych w workflow'ach Nextflow.

### Podsumowanie

W tej sekcji nauczyłeś/nauczyłaś się:

- **Dlaczego metadane są ważne:** Przechowywanie metadanych razem z danymi zachowuje ważne informacje o plikach przez cały czas trwania workflow'u.
- **Jak wczytywać arkusze danych:** Używanie `splitCsv` do odczytu plików CSV z informacjami nagłówkowymi i przekształcania wierszy w ustrukturyzowane dane.
- **Jak tworzyć mapę meta:** Oddzielanie metadanych od danych plikowych przy użyciu struktury krotki `[ [id:wartość, ...], plik ]`.

---

## 2. Manipulowanie metadanymi

Teraz, gdy mamy wczytane metadane, zróbmy z nimi coś pożytecznego!

Użyjemy narzędzia o nazwie [`langid`](https://github.com/saffsd/langid.py) do identyfikacji języka zawartego w pliku nagrania każdej postaci.
Narzędzie jest wstępnie wytrenowane na zestawie języków i dla podanego fragmentu tekstu zwróci przewidywanie języka oraz powiązany wynik prawdopodobieństwa, oba na `stdout`.

### 2.1. Zaimportuj proces i przejrzyj kod

Dostarczamy Ci gotowy moduł procesu o nazwie `IDENTIFY_LANGUAGE`, który opakowuje narzędzie `langid` — wystarczy dodać instrukcję include przed blokiem workflow.

Wprowadź następującą zmianę w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Możesz otworzyć plik modułu, aby przejrzeć jego kod:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

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

Jak widać, definicja wejścia używa tej samej struktury `tuple val(meta), path(file)`, którą właśnie zastosowaliśmy w naszym kanale wejściowym.

Definicja wyjścia jest zbudowana jako krotka o podobnej strukturze do wejścia, z tą różnicą, że zawiera również `stdout` jako trzeci element.
Ten wzorzec `tuple val(meta), path(file), <wyjście>` utrzymuje metadane powiązane zarówno z danymi wejściowymi, jak i wyjściowymi podczas przepływu przez pipeline.

Zwróć uwagę, że używamy tu kwalifikatora wyjścia [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) Nextflow, ponieważ narzędzie drukuje swoje wyniki bezpośrednio na konsolę zamiast zapisywać do pliku. Używamy też `sed` w wierszu poleceń, aby usunąć wynik prawdopodobieństwa, oczyścić string przez usunięcie znaków nowej linii i zwrócić tylko przewidywanie języka.

### 2.2. Dodaj wywołanie `IDENTIFY_LANGUAGE`

Teraz, gdy proces jest dostępny dla workflow'u, możemy dodać wywołanie procesu `IDENTIFY_LANGUAGE`, aby uruchomić go na kanale danych.

Wprowadź następujące zmiany w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Zwróć uwagę, że usunęliśmy oryginalną operację `.view()` z kodu budującego kanał.

Możemy teraz uruchomić workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Doskonale! Mamy teraz przewidywanie języka, którym posługuje się każda postać.

Jak wspomniano wcześniej, w wynikach uwzględniliśmy też plik wejściowy i mapę meta, co oznacza, że oba pozostają powiązane z nowo wygenerowanymi informacjami.
Przyda się to w następnym kroku.

!!! note "Uwaga"

    Ogólnie rzecz biorąc, ten wzorzec utrzymywania mapy meta powiązanej z wynikami ułatwia kojarzenie powiązanych wyników, które mają te same identyfikatory.

    Jak już wiesz, nie można polegać na kolejności elementów w kanałach przy kojarzeniu wyników między nimi.
    Zamiast tego należy używać kluczy do prawidłowego kojarzenia danych, a mapy meta zapewniają do tego idealną strukturę.

    Ten przypadek użycia szczegółowo omawiamy w zadaniu pobocznym [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Wzbogać metadane o wyniki procesu

Ponieważ wyniki, które właśnie uzyskaliśmy, są same w sobie formą metadanych o zawartości plików, warto byłoby dodać je do naszej mapy meta.

Nie chcemy jednak modyfikować istniejącej mapy meta w miejscu.
Z technicznego punktu widzenia jest to _możliwe_, ale niebezpieczne.

Zamiast tego utworzymy nową mapę meta zawierającą zawartość istniejącej mapy meta plus nową parę klucz-wartość `lang: lang_id` przechowującą nowe informacje, używając operatora `+` (funkcja Groovy).
Połączymy to z operacją [`map`](https://www.nextflow.io/docs/latest/operator.html#map), aby zastąpić starą mapę nową.

Oto zmiany, które należy wprowadzić w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Jeśli nie znasz jeszcze operatora `+` lub jeśli coś jest niejasne, poświęć kilka minut na przeczytanie szczegółowego wyjaśnienia poniżej.

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
    new_map = [map1 + [lang: lang_id]]
    ```

    Tutaj `[lang: new_info]` tworzy nową anonimową mapę w locie, a `map1 + ` scala `map1` z nową anonimową mapą, dając tę samą zawartość `new_map` co poprzednio.

    Eleganckie, prawda?

    **Teraz przełóżmy to na kontekst operacji `channel.map()` w Nextflow.**

    Kod przyjmuje postać:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Robi to następujące rzeczy:

    - `map1, lang_id ->` pobiera dwa elementy z krotki
    - `[map1 + [lang: lang_id]]` tworzy nową mapę zgodnie z opisem powyżej

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

    Jeśli trudno Ci śledzić, dlaczego `file` wydaje się przemieszczać w operacji `map`, wyobraź sobie, że zamiast `[meta + [lang: lang_id], file]` ten wiersz brzmi `[new_map, file]`.
    Powinno to wyjaśnić, że po prostu zostawiamy `file` na jego oryginalnym miejscu na drugiej pozycji w krotce. Wzięliśmy wartość `new_info` i wbudowaliśmy ją w mapę na pierwszej pozycji.

    **I to prowadzi nas z powrotem do struktury kanału `tuple val(meta), path(file)`!**

Gdy już rozumiesz, co robi ten kod, uruchom workflow, aby sprawdzić, czy zadziałał:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

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

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Przypisz grupę językową za pomocą instrukcji warunkowych

Teraz, gdy mamy przewidywania języków, użyjmy tych informacji do przypisania nowych grupowań.

W naszych przykładowych danych języki używane przez nasze postacie można podzielić na języki germańskie (angielski, niemiecki) i romańskie (francuski, hiszpański, włoski).
Przydatne byłoby posiadanie tej klasyfikacji w łatwo dostępnym miejscu w dalszej części pipeline'u, więc dodajmy tę informację do mapy meta.

I dobra wiadomość — to kolejny przypadek, który doskonale nadaje się do użycia operatora `map`!

Konkretnie zdefiniujemy zmienną o nazwie `lang_group` i użyjemy prostej logiki warunkowej, aby określić, jaką wartość przypisać do `lang_group` dla każdego elementu danych.

Ogólna składnia będzie wyglądać tak:

```groovy
.map { meta, file ->

    // tutaj umieszczamy logikę warunkową definiującą lang_group

    [meta + [lang_group: lang_group], file]
}
```

Widać, że jest to bardzo podobne do operacji scalania map w locie, której użyliśmy w poprzednim kroku.
Musimy tylko napisać instrukcje warunkowe.

Oto logika warunkowa, którą chcemy zastosować:

- Zdefiniuj zmienną o nazwie `lang_group` z domyślną wartością `'unknown'`.
- Jeśli `lang` to język niemiecki (`'de'`) lub angielski (`'en'`), zmień `lang_group` na `germanic`.
- W przeciwnym razie, jeśli `lang` należy do listy zawierającej francuski (`'fr'`), hiszpański (`'es'`) i włoski (`'it'`), zmień `lang_group` na `romance`.

Spróbuj napisać to samodzielnie, jeśli już wiesz, jak pisać instrukcje warunkowe w Nextflow.

!!! tip "Wskazówka"

    Dostęp do wartości `lang` w operacji map uzyskasz przez `meta.lang`.

Powinieneś/powinnaś wprowadzić następujące zmiany w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
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
            .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Uruchom langid, aby zidentyfikować język każdego pozdrowienia
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Oto kluczowe punkty:

- Używamy `def lang_group = "unknown"`, aby utworzyć zmienną `lang_group` z domyślną wartością `unknown`.
- Używamy struktury `if {} else if {}` dla logiki warunkowej, z alternatywnymi testami `.equals()` dla dwóch języków germańskich i testem przynależności do listy dla trzech języków romańskich.
- Używamy operacji scalania `meta + [lang_group:lang_group]` jak poprzednio, aby wygenerować zaktualizowaną mapę meta.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Gdy wszystko jest jasne, uruchom workflow ponownie, aby zobaczyć wynik:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Jak widać, elementy kanału zachowują strukturę `[meta, file]`, ale mapa meta zawiera teraz tę nową klasyfikację.

### Podsumowanie

W tej sekcji nauczyłeś/nauczyłaś się:

- **Stosować metadane wejściowe do kanałów wyjściowych:** Kopiowanie metadanych w ten sposób pozwala nam później kojarzyć wyniki na podstawie ich zawartości.
- **Tworzyć własne klucze:** Dodałeś/dodałaś dwa nowe klucze do mapy meta, scalając je za pomocą `meta + [nowy_klucz:wartość]` z istniejącą mapą meta — jeden oparty na wartości obliczonej przez proces, drugi na warunku ustawionym w operatorze `map`.

Pozwala to kojarzyć nowe i istniejące metadane z plikami w miarę postępu przez pipeline.
Nawet jeśli nie używasz metadanych bezpośrednio w procesie, utrzymywanie mapy meta powiązanej z danymi ułatwia przechowywanie wszystkich istotnych informacji razem.

---

## 3. Używanie informacji z mapy meta w procesie

Teraz, gdy wiesz, jak tworzyć i aktualizować mapę meta, możemy przejść do naprawdę ciekawej części: faktycznego używania metadanych w procesie.

Konkretnie dodamy drugi krok do naszego workflow'u, który narysuje każde zwierzę jako ASCII art i sprawi, że wypowie nagrany tekst w dymku.
Zrobimy to za pomocą narzędzia o nazwie [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Co robi `cowpy`?"

    `cowpy` to narzędzie wiersza poleceń generujące ASCII art do wyświetlania dowolnych tekstów w zabawny sposób.
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

Jeśli przerabiałeś/przerabiałaś kurs Hello Nextflow, widziałeś/widziałaś już to narzędzie w akcji.
Jeśli nie — nie martw się; omówimy wszystko, co musisz wiedzieć, w miarę postępu.

### 3.1. Zaimportuj proces i przejrzyj kod

Dostarczamy Ci gotowy moduł procesu o nazwie `COWPY`, który opakowuje narzędzie `cowpy` — wystarczy dodać instrukcję include przed blokiem workflow.

Wprowadź następującą zmianę w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Możesz otworzyć plik modułu, aby przejrzeć jego kod:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generuj ASCII art za pomocą cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

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

Jak widać, ten proces jest obecnie zaprojektowany tak, aby przyjmować plik wejściowy (zawierający tekst do wyświetlenia) i wartość określającą postać, która ma być narysowana jako ASCII art — zazwyczaj podawaną na poziomie workflow'u przez parametr wiersza poleceń.

### 3.2. Przekaż pole mapy meta jako wejście

Gdy używaliśmy narzędzia `cowpy` w kursie Hello Nextflow, używaliśmy parametru wiersza poleceń do określenia, jakiej postaci użyć do narysowania końcowego obrazu.
Miało to sens, ponieważ generowaliśmy tylko jeden obraz na uruchomienie pipeline'u.

Jednak w tym samouczku chcemy wygenerować odpowiedni obraz dla każdego przetwarzanego podmiotu, więc użycie parametru wiersza poleceń byłoby zbyt ograniczające.

Dobra wiadomość: mamy kolumnę `character` w naszym arkuszu danych, a co za tym idzie — w naszej mapie meta.
Użyjmy jej do ustawienia postaci, której proces powinien użyć dla każdego wpisu.

W tym celu musimy zrobić trzy rzeczy:

1. Nadać nazwę kanałowi wyjściowemu poprzedniego procesu, aby wygodniej na nim operować.
2. Określić, jak uzyskać dostęp do interesujących nas informacji.
3. Dodać wywołanie drugiego procesu i odpowiednio przekazać informacje.

Zaczynajmy.

#### 3.2.1. Nadaj nazwę poprzedniemu kanałowi wyjściowemu

Poprzednie manipulacje stosowaliśmy bezpośrednio na kanale wyjściowym pierwszego procesu, `IDENTIFY_LANGUAGE.out`.
Aby przekazać zawartość tego kanału do następnego procesu (i zrobić to w sposób czytelny i przejrzysty), chcemy nadać mu własną nazwę: `ch_languages`.

Możemy to zrobić za pomocą operatora [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

W głównym workflow zastąp operator `.view()` przez `.set { ch_languages }` i dodaj wiersz sprawdzający, czy możemy odwoływać się do kanału po nazwie.

=== "Po"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
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

        // Tymczasowo: podejrzyj zawartość ch_languages
        ch_languages.view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
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
            .view()
    ```

Uruchommy to:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Potwierdza to, że możemy teraz odwoływać się do kanału po nazwie.

#### 3.2.2. Uzyskaj dostęp do pliku i metadanych postaci

Z kodu modułu wiemy, że proces `COWPY` oczekuje pliku tekstowego i wartości `character`.
Aby napisać wywołanie procesu `COWPY`, musimy tylko wiedzieć, jak wyodrębnić odpowiedni obiekt pliku i metadane z każdego elementu kanału.

Jak to często bywa, najprostszym sposobem jest użycie operacji `map`.

Nasz kanał zawiera krotki o strukturze `[meta, file]`, więc możemy uzyskać dostęp do obiektu `file` bezpośrednio, a do wartości `character` przechowywanej w mapie meta — przez `meta.character`.

W głównym workflow wprowadź następujące zmiany:

=== "Po"

    ```groovy title="main.nf" linenums="34"
        // Tymczasowo: uzyskaj dostęp do pliku i postaci
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="34"
        // Tymczasowo: podejrzyj zawartość ch_languages
        ch_languages.view()
    ```

Zwróć uwagę, że używamy domknięć (takich jak `{ file -> "File: " + file }`) aby wyniki operacji `.view` były bardziej czytelne.

Uruchommy to:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Ścieżki plików i wartości postaci mogą pojawić się w innej kolejności w Twoich wynikach._

Potwierdza to, że możemy uzyskać dostęp do pliku i postaci dla każdego elementu kanału.

#### 3.2.3. Wywołaj proces `COWPY`

Teraz złóżmy to wszystko razem i faktycznie wywołajmy proces `COWPY` na kanale `ch_languages`.

W głównym workflow wprowadź następujące zmiany:

=== "Po"

    ```groovy title="main.nf" linenums="34"
        // Uruchom cowpy, aby wygenerować ASCII art
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="34"
        // Tymczasowo: uzyskaj dostęp do pliku i postaci
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Widać, że po prostu kopiujemy dwie operacje map (bez instrukcji `.view()`) jako wejścia do wywołania procesu.
Tylko nie zapomnij o przecinku między nimi!

To trochę niezgrabne, ale w następnej sekcji zobaczymy, jak to poprawić.

Uruchommy to:

```bash
nextflow run main.nf -resume
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Jeśli zajrzysz do katalogu wyników, powinieneś/powinnaś zobaczyć poszczególne pliki zawierające ASCII art każdego pozdrowienia wypowiadanego przez odpowiednią postać.

??? abstract "Zawartość katalogu i przykładowego pliku"

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

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Pokazuje to, że byliśmy w stanie użyć informacji z mapy meta do sparametryzowania polecenia w drugim kroku pipeline'u.

Jak jednak wspomniano, część kodu była nieco niezgrabna, ponieważ musieliśmy rozpakowywać metadane jeszcze w kontekście ciała workflow'u.
To podejście sprawdza się dobrze przy używaniu niewielkiej liczby pól z mapy meta, ale słabo skaluje się, gdy chcemy używać ich znacznie więcej.

Istnieje inny operator o nazwie `multiMap()`, który pozwala nieco to usprawnić, ale nawet wtedy nie jest to idealne rozwiązanie.

??? info "(Opcjonalnie) Alternatywna wersja z `multiMap()`"

    Możesz się zastanawiać, dlaczego nie możemy napisać jednej operacji `map()`, która zwraca zarówno `file`, jak i `character` — otóż zwróciłaby je jako krotkę.
    Musieliśmy napisać dwie oddzielne operacje `map()`, aby przekazać elementy `file` i `character` do procesu osobno.

    Technicznie istnieje inny sposób na zrobienie tego przez jedną operację mapowania, używając operatora `multiMap()`, który jest w stanie emitować wiele kanałów.
    Na przykład możesz zastąpić wywołanie `COWPY` powyżej następującym kodem:

    === "Po"

        ```groovy title="main.nf" linenums="34"
            // Uruchom cowpy, aby wygenerować ASCII art
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Przed"

        ```groovy title="main.nf" linenums="34"
            // Uruchom cowpy, aby wygenerować ASCII art
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Daje to dokładnie ten sam wynik.

W obu przypadkach jest to trochę niezręczne, że musimy wykonywać część rozpakowywania na poziomie workflow'u.

Lepiej byłoby móc przekazać całą mapę meta do procesu i wybrać to, czego potrzebujemy, już wewnątrz niego.

### 3.3. Przekaż i użyj całej mapy meta

Celem mapy meta jest przecież przekazywanie wszystkich metadanych razem jako pakietu.
Jedynym powodem, dla którego nie mogliśmy tego zrobić powyżej, jest to, że proces nie jest skonfigurowany do przyjmowania mapy meta.
Ale ponieważ kontrolujemy kod procesu, możemy to zmienić.

Zmodyfikujmy proces `COWPY`, aby przyjmował strukturę krotki `[meta, file]`, której użyliśmy w pierwszym procesie, co pozwoli nam uprościć workflow.

W tym celu musimy zrobić trzy rzeczy:

1. Zmodyfikować definicje wejść modułu procesu `COWPY`
2. Zaktualizować polecenie procesu, aby używało mapy meta
3. Zaktualizować wywołanie procesu w ciele workflow'u

Gotowy/gotowa? Zaczynamy!

#### 3.3.1. Zmodyfikuj wejście modułu `COWPY`

Wprowadź następujące zmiany w pliku modułu `cowpy.nf`:

=== "Po"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Przed"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Pozwala nam to używać struktury krotki `[meta, file]` omówionej wcześniej w samouczku.

Zwróć uwagę, że nie zaktualizowaliśmy definicji wyjścia procesu, aby zwracała mapę meta — w celu uproszczenia samouczka. Możesz jednak zrobić to samodzielnie jako ćwiczenie, wzorując się na procesie `IDENTIFY_LANGUAGE`.

#### 3.3.2. Zaktualizuj polecenie, aby używało pola mapy meta

Cała mapa meta jest teraz dostępna wewnątrz procesu, więc możemy odwoływać się do zawartych w niej informacji bezpośrednio z bloku polecenia.

Wprowadź następujące zmiany w pliku modułu `cowpy.nf`:

=== "Po"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Przed"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Zastąpiliśmy odwołanie do wartości `character` przekazywanej wcześniej jako samodzielne wejście wartością przechowywaną w mapie meta, do której odwołujemy się przez `meta.character`.

Teraz zaktualizujmy odpowiednio wywołanie procesu.

#### 3.3.3. Zaktualizuj wywołanie procesu i uruchom go

Proces oczekuje teraz wejścia w strukturze krotki `[meta, file]`, co jest dokładnie tym, co zwraca poprzedni proces — możemy więc po prostu przekazać cały kanał `ch_languages` do procesu `COWPY`.

Wprowadź następujące zmiany w głównym workflow:

=== "Po"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Uruchom cowpy, aby wygenerować ASCII art
    COWPY(ch_languages)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Uruchom cowpy, aby wygenerować ASCII art
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

To znacznie upraszcza wywołanie!

Usuńmy wyniki poprzedniego uruchomienia i uruchommy workflow:

```bash
rm -r results
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Jeśli zajrzysz do katalogu wyników, powinieneś/powinnaś zobaczyć te same wyniki co poprzednio — czyli poszczególne pliki zawierające ASCII art każdego pozdrowienia wypowiadanego przez odpowiednią postać.

??? abstract "Zawartość katalogu"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Wyniki są więc takie same jak poprzednio, ale kod jest prostszy.

Oczywiście zakłada to, że możesz modyfikować kod procesu.
W niektórych przypadkach możesz być zmuszony/zmuszona do korzystania z istniejących procesów, których nie możesz modyfikować, co ogranicza Twoje możliwości.
Dobra wiadomość dla tych, którzy planują używać modułów z projektu [nf-core](https://nf-co.re/), jest taka, że moduły nf-core są wszystkie skonfigurowane do używania struktury krotki `[meta, file]` jako standardu.

### 3.4. Rozwiązywanie problemów z brakującymi wymaganymi wejściami

Wartość `character` jest wymagana, aby proces `COWPY` działał poprawnie.
Jeśli nie ustawimy dla niej wartości domyślnej w pliku konfiguracyjnym, MUSIMY podać ją w arkuszu danych.

**Co się stanie, jeśli tego nie zrobimy?**
Zależy to od zawartości arkusza danych i wersji workflow'u, której używamy.

#### 3.4.1. Kolumna character istnieje, ale jest pusta

Powiedzmy, że usuwamy wartość character dla jednego z wpisów w naszym arkuszu danych, symulując błąd zbierania danych:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

W obu wersjach workflow'u, których używaliśmy powyżej, klucz `character` zostanie utworzony dla wszystkich wpisów podczas wczytywania arkusza danych, ale dla `sampleA` wartość będzie pustym stringiem.

Spowoduje to błąd.

??? failure "Wyjście polecenia"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

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

Gdy Nextflow uruchamia wiersz poleceń `cowpy` dla tej próbki, `${meta.character}` jest wypełniane pustym stringiem, więc narzędzie `cowpy` zgłasza błąd informujący, że nie podano żadnej wartości dla argumentu `-c`.

#### 3.4.2. Kolumna character nie istnieje w arkuszu danych

Teraz powiedzmy, że całkowicie usuwamy kolumnę `character` z naszego arkusza danych:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

W tym przypadku klucz `character` w ogóle nie zostanie utworzony podczas wczytywania arkusza danych.

##### 3.4.2.1. Wartość dostępna na poziomie workflow'u

Jeśli używamy wersji kodu z sekcji 3.2, Nextflow spróbuje uzyskać dostęp do klucza `character` w mapie meta PRZED wywołaniem procesu `COWPY`.

Nie znajdzie żadnych elementów pasujących do instrukcji, więc w ogóle nie uruchomi `COWPY`.

??? success "Wyjście polecenia"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Z punktu widzenia Nextflow ten workflow zakończył się pomyślnie!
Jednak żadne z oczekiwanych wyników nie zostaną wygenerowane.

##### 3.4.2.2. Wartość dostępna na poziomie procesu

Jeśli używamy wersji z sekcji 3.3, Nextflow przekaże całą mapę meta do procesu `COWPY` i spróbuje uruchomić polecenie.

Spowoduje to błąd, ale inny niż w pierwszym przypadku.

??? failure "Wyjście polecenia"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

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

Dzieje się tak, ponieważ `meta.character` nie istnieje, więc próba uzyskania do niego dostępu zwraca `null`. W rezultacie Nextflow dosłownie wstawia `null` do wiersza poleceń, co oczywiście nie jest rozpoznawane przez narzędzie `cowpy`.

#### 3.4.3. Rozwiązania

Poza ustawieniem wartości domyślnej w konfiguracji workflow'u, możemy zrobić dwie rzeczy, aby obsłużyć to bardziej niezawodnie:

1. Zaimplementować walidację wejść w workflow'u, aby upewnić się, że arkusz danych zawiera wszystkie wymagane informacje. Możesz znaleźć [wprowadzenie do walidacji wejść](../hello_nf-core/05_input_validation.md) w kursie szkoleniowym Hello nf-core. <!-- TODO (future) pending a proper Validation side quest -->

2. Jeśli chcesz mieć pewność, że każdy, kto używa Twojego modułu procesu, może natychmiast zidentyfikować wymagane wejścia, możesz też uczynić wymaganą właściwość metadanych jawnym wejściem.

Oto przykład, jak by to działało.

Po pierwsze, na poziomie procesu zaktualizuj definicję wejścia w następujący sposób:

=== "Po"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Przed"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Następnie na poziomie workflow'u użyj operacji mapowania, aby wyodrębnić właściwość `character` z metadanych i uczynić ją jawnym składnikiem krotki wejściowej:

=== "Po"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

To podejście ma tę zaletę, że jawnie pokazuje, że `character` jest wymagany, i ułatwia ponowne wdrożenie procesu w innych kontekstach.

Podkreśla to ważną zasadę projektowania:

**Używaj mapy meta dla opcjonalnych, opisowych informacji, ale wyodrębniaj wymagane wartości jako jawne wejścia.**

Mapa meta doskonale nadaje się do utrzymywania czystych struktur kanałów i zapobiegania arbitralnym strukturom, ale dla obowiązkowych elementów bezpośrednio przywoływanych w procesie wyodrębnienie ich jako jawnych wejść tworzy bardziej niezawodny i łatwy w utrzymaniu kod.

### Podsumowanie

W tej sekcji nauczyłeś/nauczyłaś się, jak wykorzystywać metadane do dostosowywania wykonania procesu, uzyskując do nich dostęp zarówno na poziomie workflow'u, jak i na poziomie procesu.

---

## Ćwiczenie uzupełniające

Jeśli chcesz poćwiczyć używanie informacji z mapy meta wewnątrz procesu, spróbuj użyć innych informacji z mapy meta, takich jak `lang` i `lang_group`, aby dostosować nazwy i/lub organizację wyników.

Na przykład spróbuj zmodyfikować kod, aby uzyskać następujący wynik:

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

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Podsumowanie

W tym zadaniu pobocznym zbadałeś/zbadałaś, jak efektywnie pracować z metadanymi w workflow'ach Nextflow.

Ten wzorzec jawnego przechowywania metadanych razem z danymi jest podstawową dobrą praktyką w Nextflow, oferującą kilka zalet w porównaniu z zakodowywaniem informacji o plikach na stałe:

- Metadane pliku pozostają powiązane z plikami przez cały czas trwania workflow'u
- Zachowanie procesu można dostosować dla każdego pliku
- Organizacja wyników może odzwierciedlać metadane pliku
- Informacje o plikach można rozszerzać podczas wykonywania pipeline'u

Stosowanie tego wzorca we własnej pracy pozwoli Ci budować niezawodne, łatwe w utrzymaniu workflow'y bioinformatyczne.

### Kluczowe wzorce

1.  **Odczytywanie i strukturyzowanie metadanych:** Odczytywanie plików CSV i tworzenie zorganizowanych map metadanych, które pozostają powiązane z plikami danych.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Rozszerzanie metadanych podczas workflow'u:** Dodawanie nowych informacji do metadanych w miarę postępu pipeline'u przez dodawanie wyników procesów i wyprowadzanie wartości przez logikę warunkową.

    - Dodawanie nowych kluczy na podstawie wyników procesu

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Dodawanie nowych kluczy za pomocą klauzuli warunkowej

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Dostosowywanie zachowania procesu:** Używanie metadanych wewnątrz procesu.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Dodatkowe zasoby

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Co dalej?

Wróć do [menu zadań pobocznych](../) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
