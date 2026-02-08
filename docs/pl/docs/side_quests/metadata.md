# Metadane i mapy meta

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W każdej analizie naukowej rzadko pracujemy jedynie z surowymi plikami danych.
Każdy plik ma Swoje dodatkowe informacje: czym jest, skąd pochodzi i co czyni go wyjątkowym.
Te dodatkowe informacje nazywamy metadanymi.

Metadane to dane opisujące inne dane.
Metadane śledzą ważne szczegóły dotyczące plików i warunków eksperymentalnych oraz pomagają dostosować analizy do unikalnych cech każdego zestawu danych.

Pomyśl o tym jak o katalogu bibliotecznym: podczas gdy książki zawierają faktyczną treść (surowe dane), karty katalogowe dostarczają istotnych informacji o każdej książce — kiedy została opublikowana, kto ją napisał, gdzie ją znaleźć (metadane).
W pipeline'ach Nextflow metadane mogą być wykorzystywane do:

- Śledzenia informacji specyficznych dla pliku w całym workflow
- Konfigurowania procesów na podstawie cech pliku
- Grupowania powiązanych plików do wspólnej analizy

### Cele edukacyjne

W tym side queście zbadamy sposób obsługi metadanych w workflow.
Zaczynając od prostej tabeli danych (często nazywanej samplesheet w bioinformatyce) zawierającej podstawowe informacje o plikach, nauczysz się jak:

- Odczytywać i analizować metadane plików z plików CSV
- Tworzyć i manipulować mapami metadanych
- Dodawać nowe pola metadanych podczas wykonywania workflow
- Wykorzystywać metadane do dostosowywania zachowania procesów

Te umiejętności pomogą Ci tworzyć bardziej solidne i elastyczne pipeline'y, które potrafią obsługiwać złożone relacje między plikami i wymagania przetwarzania.

### Wymagania wstępne

Przed rozpoczęciem tego side questu powinieneś:

- Ukończyć kurs [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory)

---

## 0. Rozpocznij pracę

#### Otwórz środowisko szkoleniowe codespace

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzysz środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki do tego kursu.

```bash
cd side-quests/metadata
```

Możesz ustawić VSCode, aby skupił się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz główny plik workflow oraz katalog `data` zawierający tabelę danych i kilka plików danych.

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

Workflow w pliku `main.nf` to szablon, który stopniowo rozwiniesz w w pełni funkcjonujący workflow.

Tabela danych zawiera ścieżki do plików danych i powiązane metadane, zorganizowane w 3 kolumnach:

- `id`: oczywiste, ID nadane plikowi
- `character`: nazwa postaci, której użyjemy później do rysowania różnych stworzeń
- `data`: ścieżki do plików `.txt`, które zawierają pozdrowienia w różnych językach

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

Każdy plik danych zawiera tekst powitania w jednym z pięciu języków (fr: francuski, de: niemiecki, es: hiszpański, it: włoski, en: angielski).

Dostarczymy Ci również konteneryzowane narzędzie do analizy języka o nazwie `langid`.

#### Przejrzyj zadanie

Twoim wyzwaniem jest napisanie workflow Nextflow, który:

1. **Zidentyfikuje** język w każdym pliku automatycznie
2. **Zgrupuje** pliki według rodziny językowej (języki germańskie kontra języki romańskie)
3. **Dostosuje** przetwarzanie dla każdego pliku na podstawie jego języka i metadanych
4. **Zorganizuje** wyjścia według grupy językowej

To reprezentuje typowy wzorzec workflow, gdzie metadane specyficzne dla pliku sterują decyzjami dotyczącymi przetwarzania; dokładnie taki rodzaj problemu, który mapy metadanych rozwiązują elegancko.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby się zagłębić?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace jest uruchomiony
- [ ] Ustawiłem odpowiednio mój katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Wczytaj metadane z tabeli danych

Otwórz plik workflow `main.nf`, aby przeanalizować szablon workflow, który dajemy Ci jako punkt wyjścia.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Możesz zobaczyć, że skonfigurowaliśmy podstawową fabrykę kanałów, aby wczytać przykładową tabelę danych jako plik, ale to jeszcze nie odczyta zawartości pliku.
Zacznijmy od dodania tego.

### 1.1. Odczytaj zawartość za pomocą `splitCsv`

Musimy wybrać operator, który odpowiednio zanalizuje zawartość pliku przy minimalnym wysiłku z naszej strony.
Ponieważ nasza tabela danych jest w formacie CSV, jest to zadanie dla operatora [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), który wczytuje każdy wiersz w pliku jako element w kanale.

Wprowadź następujące zmiany, aby dodać operację `splitCsv()` do kodu konstrukcji kanału, plus operację `view()`, aby sprawdzić, czy zawartość pliku jest poprawnie wczytywana do kanału.

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

Zauważ, że używamy opcji `header: true`, aby powiedzieć Nextflow, aby odczytał pierwszy wiersz pliku CSV jako wiersz nagłówka.

Zobaczmy, co z tego wychodzi, dobrze?
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

Możemy zobaczyć, że operator skonstruował mapę par klucz-wartość dla każdego wiersza w pliku CSV, z nagłówkami kolumn jako kluczami dla odpowiadających wartości.

Każdy wpis mapy odpowiada kolumnie w naszej tabeli danych:

- `id`
- `character`
- `recording`

To jest świetne! Ułatwia to dostęp do określonych pól z każdego pliku.
Na przykład moglibyśmy uzyskać dostęp do ID pliku za pomocą `id` lub ścieżki do pliku txt za pomocą `recording`.

??? info "(Opcjonalnie) Więcej o mapach"

    W Groovy, języku programowania, na którym zbudowany jest Nextflow, mapa to struktura danych klucz-wartość podobna do słowników w Python, obiektów w JavaScript lub haszy w Ruby.

    Oto skrypt do uruchomienia, który pokazuje, jak można zdefiniować mapę i uzyskać dostęp do jej zawartości w praktyce:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Utwórz prostą mapę
    def my_map = [id:'sampleA', character:'squirrel']

    // Wypisz całą mapę
    println "map: ${my_map}"

    // Uzyskaj dostęp do pojedynczych wartości za pomocą notacji kropkowej
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Mimo że nie ma odpowiedniego bloku `workflow`, Nextflow może to uruchomić tak, jakby to był workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    A oto, co możesz spodziewać się zobaczyć w wyjściu:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Wyodrębnij określone pola za pomocą `map`

Powiedzmy, że chcemy uzyskać dostęp do kolumny `character` z tabeli danych i ją wydrukować.
Możemy użyć operatora Nextflow `map`, aby iterować po każdym elemencie w naszym kanale i specyficznie wybrać wpis `character` z obiektu mapy.

Wprowadź następujące edycje do workflow:

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

Teraz uruchom workflow ponownie:

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

Sukces! Wykorzystaliśmy strukturę mapy pochodzącą z naszej tabeli danych, aby uzyskać dostęp do wartości z poszczególnych kolumn dla każdego wiersza.

Teraz, gdy pomyślnie odczytaliśmy tabelę danych i mamy dostęp do danych w każdym wierszu, możemy zacząć implementować logikę naszego pipeline.

### 1.3. Zorganizuj metadane w 'mapę meta'

W obecnym stanie workflow pliki wejściowe (pod kluczem `recording`) i powiązane metadane (`id`, `character`) są na tym samym poziomie, jak gdyby wszystkie były w jednej dużej torbie.
Praktyczną konsekwencją jest to, że każdy proces, który konsumuje ten kanał, musiałby być skonfigurowany z tą strukturą na myśli:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

To jest w porządku, dopóki liczba kolumn w tabeli danych się nie zmienia.
Jednakże, jeśli dodasz nawet tylko jedną kolumnę do tabeli danych, kształt kanału nie będzie już pasował do tego, czego oczekuje proces, a workflow będzie generował błędy.
Utrudnia to również dzielenie się procesem z innymi, którzy mogą mieć nieco inne dane wejściowe, i możesz skończyć na hardcodowaniu zmiennych do procesu, które nie są potrzebne przez blok script.

Aby uniknąć tego problemu, musimy znaleźć sposób na utrzymanie struktury kanału spójnej niezależnie od tego, ile kolumn zawiera tabela danych.

Możemy to zrobić, zbierając wszystkie metadane w element w krotce, który nazwiemy mapą metadanych, lub prościej 'mapą meta'.

Wprowadź następujące edycje do operacji `map`:

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

Zrestrukturyzowaliśmy nasze elementy kanału w krotkę składającą się z dwóch elementów: mapy meta i odpowiadającego obiektu pliku.

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

Teraz każdy element w kanale zawiera najpierw mapę metadanych, a na drugim miejscu odpowiadający obiekt pliku:

```console title="Przykładowa struktura wyjścia"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

W rezultacie dodanie większej liczby kolumn w tabeli danych udostępni więcej metadanych w mapie `meta`, ale nie zmieni kształtu kanału.
To umożliwia nam pisanie procesów, które konsumują kanał bez konieczności hardcodowania elementów metadanych w specyfikacji wejściowej:

```groovy title="Przykład składni"
    input:
    tuple val(meta), file(recording)
```

Jest to powszechnie stosowana konwencja organizowania metadanych w workflow Nextflow.

### Podsumowanie

W tej sekcji nauczyłeś się:

- **Dlaczego metadane są ważne:** Utrzymywanie metadanych wraz z danymi zachowuje ważne informacje o plikach w całym workflow.
- **Jak odczytywać tabele danych:** Używanie `splitCsv` do odczytu plików CSV z informacjami nagłówkowymi i przekształcania wierszy w ustrukturyzowane dane
- **Jak utworzyć mapę meta:** Oddzielanie metadanych od danych pliku za pomocą struktury krotki `[ [id:value, ...], file ]`

---

## 2. Manipulowanie metadanymi

Teraz, gdy mamy załadowane metadane, zróbmy z nimi coś!

Użyjemy narzędzia o nazwie [`langid`](https://github.com/saffsd/langid.py) do identyfikacji języka zawartego w pliku nagrania każdego stworzenia.
Narzędzie jest wstępnie wytrenowane na zestawie języków i, biorąc fragment tekstu, wyprowadzi przewidywanie języka i powiązany wynik prawdopodobieństwa, oba do `stdout`.

### 2.1. Zaimportuj proces i przeanalizuj kod

Dostarczamy Ci wstępnie napisany moduł procesu o nazwie `IDENTIFY_LANGUAGE`, który opakowuje narzędzie `langid`, więc musisz tylko dodać instrukcję include przed blokiem workflow.

Wprowadź następującą edycję do workflow:

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

Możesz otworzyć plik modułu, aby przeanalizować jego kod:

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

Jak widać, definicja wejściowa używa tej samej struktury `tuple val(meta), path(file)`, którą właśnie zastosowaliśmy do naszego kanału wejściowego.

Definicja wyjściowa jest zbudowana jako krotka o podobnej strukturze do wejścia, z wyjątkiem tego, że zawiera również `stdout` jako trzeci element.
Ten wzorzec `tuple val(meta), path(file), <output>` utrzymuje metadane powiązane zarówno z danymi wejściowymi, jak i wyjściowymi w miarę przepływu przez pipeline.

Zauważ, że używamy tu kwalifikatora wyjściowego [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) Nextflow, ponieważ narzędzie drukuje Swoje wyjście bezpośrednio do konsoli zamiast zapisywać plik; i używamy `sed` w wierszu poleceń, aby usunąć wynik prawdopodobieństwa, oczyścić ciąg znaków poprzez usunięcie znaków nowej linii i zwrócić tylko przewidywanie języka.

### 2.2. Dodaj wywołanie `IDENTIFY_LANGUAGE`

Teraz, gdy proces jest dostępny dla workflow, możemy dodać wywołanie procesu `IDENTIFY_LANGUAGE`, aby uruchomić go na kanale danych.

Wprowadź następujące edycje do workflow:

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

Zauważ, że usunęliśmy oryginalną operację `.view()` w konstrukcji kanału.

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

Doskonale! Mamy teraz przewidywanie, w jakim języku mówi każda postać.

I jak wspomniano wcześniej, uwzględniliśmy również plik wejściowy i mapę meta w wyjściu, co oznacza, że oba pozostają powiązane z nowymi informacjami, które właśnie wygenerowaliśmy.
To okaże się przydatne w następnym kroku.

!!! note

    Bardziej ogólnie, ten wzorzec utrzymywania mapy meta powiązanej z wynikami ułatwia kojarzenie powiązanych wyników, które dzielą te same identyfikatory.

    Jak już się dowiedziałeś, nie możesz polegać na kolejności elementów w kanałach, aby dopasować wyniki między nimi.
    Zamiast tego musisz używać kluczy do prawidłowego kojarzenia danych, a mapy meta zapewniają idealną strukturę do tego celu.

    Badamy ten przypadek użycia szczegółowo w side queście [Splitting & Grouping](./splitting_and_grouping.md).

### 2.3. Rozszerz metadane o wyjścia procesu

Biorąc pod uwagę, że wyniki, które właśnie wygenerowaliśmy, same w sobie są formą metadanych o zawartości plików, przydatne byłoby dodanie ich do naszej mapy meta.

Jednakże nie chcemy modyfikować istniejącej mapy meta w miejscu.
Z technicznego punktu widzenia _możliwe_ jest to zrobienie, ale jest to niebezpieczne.

Zamiast tego stworzymy nową mapę meta zawierającą zawartość istniejącej mapy meta plus nową parę klucz-wartość `lang: lang_id` przechowującą nowe informacje, używając operatora `+` (funkcja Groovy).
I połączymy to z operacją [`map`](https://www.nextflow.io/docs/latest/operator.html#map), aby zastąpić starą mapę nową.

Oto edycje, które musisz wprowadzić do workflow:

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

Jeśli nie jesteś jeszcze zaznajomiony z operatorem `+`, lub jeśli wydaje się to mylące, poświęć kilka minut na przejrzenie szczegółowego wyjaśnienia poniżej.

??? info "Tworzenie nowej mapy meta przy użyciu operatora `+`"

    **Po pierwsze, musisz wiedzieć, że możemy scalić zawartość dwóch map używając operatora Groovy `+`.**

    Powiedzmy, że mamy następujące mapy:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Możemy je scalić w ten sposób:

    ```groovy
    new_map = map1 + map2
    ```

    Zawartość `new_map` będzie:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Świetnie!

    **Ale co jeśli musisz dodać pole, które nie jest jeszcze częścią mapy?**

    Powiedzmy, że zaczynasz ponownie od `map1`, ale przewidywanie języka nie jest w Swojej własnej mapie (nie ma `map2`).
    Zamiast tego jest przechowywane w zmiennej o nazwie `lang_id` i wiesz, że chcesz przechować jej wartość (`'fr'`) z kluczem `lang`.

    Możesz faktycznie zrobić następująco:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Tutaj `[lang: new_info]` tworzy nową nienazwaną mapę w locie, a `map1 + ` scala `map1` z nową nienazwaną mapą, produkując tę samą zawartość `new_map` jak wcześniej.

    Ładne, prawda?

    **Teraz przetransponujmy to w kontekst operacji Nextflow `channel.map()`.**

    Kod staje się:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    To robi następującą rzecz:

    - `map1, lang_id ->` przyjmuje dwa elementy w krotce
    - `[map1 + [lang: lang_id]]` tworzy nową mapę jak szczegółowo opisano powyżej

    Wyjściem jest jedna nienazwana mapa z tą samą zawartością co `new_map` w naszym przykładzie powyżej.
    Więc skutecznie przekształciliśmy:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    w:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Miejmy nadzieję, że możesz zobaczyć, że jeśli zmienimy `map1` na `meta`, to w zasadzie wszystko, czego potrzebujemy, aby dodać przewidywanie języka do naszej mapy meta w naszym workflow.

    Z wyjątkiem jednej rzeczy!

    W przypadku naszego workflow **musimy również uwzględnić obecność obiektu `file` w krotce**, który składa się z `meta, file, lang_id`.

    Więc kod tutaj stałby się:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Jeśli masz trudności ze zrozumieniem, dlaczego `file` wydaje się poruszać w operacji `map`, wyobraź sobie, że zamiast `[meta + [lang: lang_id], file]`, ten wiersz brzmi `[new_map, file]`.
    To powinno wyjaśnić, że po prostu pozostawiamy `file` w jego pierwotnym miejscu na drugiej pozycji w krotce. Po prostu wzięliśmy wartość `new_info` i włożyliśmy ją do mapy, która jest na pierwszej pozycji.

    **I to prowadzi nas z powrotem do struktury kanału `tuple val(meta), path(file)`!**

Gdy będziesz pewien, że rozumiesz, co robi ten kod, uruchom workflow, aby zobaczyć, czy zadziałało:

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
Starannie zreorganizowaliśmy wyjście procesu z `meta, file, lang_id` tak, że `lang_id` jest teraz jednym z kluczy w mapie meta, a krotki kanału ponownie pasują do modelu `meta, file`.

### 2.4. Przypisz grupę językową używając instrukcji warunkowych

Teraz, gdy mamy nasze przewidywania języków, użyjmy tych informacji do przypisania nowych grup.

W naszych przykładowych danych języki używane przez nasze postacie można podzielić na języki germańskie (angielski, niemiecki) i języki romańskie (francuski, hiszpański, włoski).
Przydatne może być posiadanie tej klasyfikacji łatwo dostępnej gdzieś później w pipeline, więc dodajmy te informacje do mapy meta.

I, dobra wiadomość, to kolejny przypadek, który idealnie nadaje się do użycia operatora `map`!

Konkretnie, zdefiniujemy zmienną o nazwie `lang_group`, użyjemy prostej logiki warunkowej do określenia, jaką wartość przypisać `lang_group` dla każdego elementu danych.

Ogólna składnia będzie wyglądać tak:

```groovy
.map { meta, file ->

    // logika warunkowa definiująca lang_group idzie tutaj

    [meta + [lang_group: lang_group], file]
}
```

Możesz zobaczyć, że jest to bardzo podobne do operacji scalania map w locie, którą użyliśmy w poprzednim kroku.
Musimy tylko napisać instrukcje warunkowe.

Oto logika warunkowa, którą chcemy zastosować:

- Zdefiniuj zmienną o nazwie `lang_group` z wartością domyślną `'unknown'`.
- Jeśli `lang` to niemiecki (`'de'`) lub angielski (`'en'`), zmień `lang_group` na `germanic`.
- W przeciwnym razie, jeśli `lang` jest zawarty w liście zawierającej francuski (`'fr'`), hiszpański (`'es'`) i włoski (`'it'`), zmień `lang_group` na `romance`.

Spróbuj napisać to sam, jeśli już wiesz, jak pisać instrukcje warunkowe w Nextflow.

!!! tip

    Możesz uzyskać dostęp do wartości `lang` wewnątrz operacji map za pomocą `meta.lang`.

Powinieneś wprowadzić następujące zmiany do workflow:

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

- Używamy `def lang_group = "unknown"` do utworzenia zmiennej `lang_group` z wartością domyślną ustawioną na `unknown`.
- Używamy struktury `if {} else if {}` dla logiki warunkowej, z alternatywnymi testami `.equals()` dla dwóch języków germańskich oraz testem istnienia w liście dla trzech języków romańskich.
- Używamy operacji scalania `meta + [lang_group:lang_group]` jak wcześniej do wygenerowania zaktualizowanej mapy meta.

Gdy to wszystko ma sens, uruchom workflow ponownie, aby zobaczyć wynik:

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

Jak widać, elementy kanału utrzymują Swoją strukturę `[meta, file]`, ale mapa meta zawiera teraz tę nową klasyfikację.

### Podsumowanie

W tej sekcji nauczyłeś się, jak:

- **Stosować metadane wejściowe do kanałów wyjściowych**: Kopiowanie metadanych w ten sposób pozwala nam później kojarzyć wyniki na podstawie zawartości metadanych.
- **Tworzyć niestandardowe klucze**: Utworzyłeś dwa nowe klucze w Swojej mapie meta, scalając je za pomocą `meta + [new_key:value]` z istniejącą mapą meta. Jeden oparty na obliczonej wartości z procesu, a drugi oparty na warunku ustawionym w operatorze `map`.

Pozwalają one kojarzyć nowe i istniejące metadane z plikami w miarę postępu w pipeline.
Nawet jeśli nie używasz metadanych jako części procesu, utrzymywanie mapy meta powiązanej z danymi w ten sposób ułatwia utrzymanie wszystkich istotnych informacji razem.

---

## 3. Używanie informacji z mapy meta w procesie

Teraz, gdy wiesz, jak tworzyć i aktualizować mapę meta, możemy przejść do naprawdę zabawnej części: faktycznego używania metadanych w procesie.

Dokładniej, dodamy drugi krok do naszego workflow, aby narysować każde zwierzę jako sztukę ASCII i sprawić, że powie nagrany tekst w dymku.
Zrobimy to za pomocą narzędzia o nazwie [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Co robi `cowpy`?"

    `cowpy` to narzędzie wiersza poleceń, które generuje sztukę ASCII do wyświetlania dowolnych wejść tekstowych w zabawny sposób.
    Jest to implementacja w pythonie klasycznego narzędzia [cowsay](https://en.wikipedia.org/wiki/Cowsay) Tony'ego Monroe.

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

    Opcjonalnie możesz wybrać postać (lub 'cowacter') do użycia zamiast domyślnej krowy.

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

Jeśli przeszedłeś kurs Hello Nextflow, już widziałeś to narzędzie w akcji.
Jeśli nie, nie martw się; omówimy wszystko, co musisz wiedzieć w trakcie.

### 3.1. Zaimportuj proces i przeanalizuj kod

Dostarczamy Ci wstępnie napisany moduł procesu o nazwie `COWPY`, który opakowuje narzędzie `cowpy`, więc musisz tylko dodać instrukcję include przed blokiem workflow.

Wprowadź następującą edycję do workflow:

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

Możesz otworzyć plik modułu, aby przeanalizować jego kod:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Wygeneruj grafikę ASCII za pomocą cowpy
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

Jak widać, ten proces jest obecnie zaprojektowany tak, aby przyjmować plik wejściowy (zawierający tekst do wyświetlenia) i wartość określającą postać, która powinna być narysowana w sztuce ASCII, zwykle dostarczaną na poziomie workflow przez parametr wiersza poleceń.

### 3.2. Przekaż pole mapy meta jako wejście

Kiedy używaliśmy narzędzia `cowpy` w kursie Hello Nextflow, użyliśmy parametru wiersza poleceń do określenia, jakiej postaci użyć do narysowania końcowego obrazu.
To miało sens, ponieważ generowaliśmy tylko jeden obraz na uruchomienie pipeline.

Jednakże w tym kursie chcemy wygenerować odpowiedni obraz dla każdego podmiotu, który przetwarzamy, więc użycie parametru wiersza poleceń byłoby zbyt ograniczające.

Dobra wiadomość: mamy kolumnę `character` w naszej tabeli danych i dlatego w naszej mapie meta.
Użyjmy tego do ustawienia postaci, której proces powinien użyć dla każdego wpisu.

W tym celu będziemy musieli zrobić trzy rzeczy:

1. Nadać nazwę kanałowi wyjściowemu wychodzącemu z poprzedniego procesu, abyśmy mogli na nim wygodniej operować.
2. Określić, jak uzyskać dostęp do interesujących nas informacji
3. Dodać wywołanie drugiego procesu i odpowiednio wprowadzić informacje.

Zaczynajmy.

#### 3.2.1. Nazwij poprzedni kanał wyjściowy

Zastosowaliśmy poprzednie manipulacje bezpośrednio na kanale wyjściowym pierwszego procesu, `IDENTIFY_LANGUAGE.out`.
Aby przekazać zawartość tego kanału do następnego procesu (i zrobić to w sposób jasny i łatwy do odczytania) chcemy nadać mu własną nazwę, `ch_languages`.

Możemy to zrobić za pomocą operatora [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

W głównym workflow zastąp operator `.view()` przez `.set { ch_languages }` i dodaj linię testującą, że możemy odwoływać się do kanału po nazwie.

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

        // Tymczasowe: podejrzyj ch_languages
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

To potwierdza, że możemy teraz odwoływać się do kanału po nazwie.

#### 3.2.2. Uzyskaj dostęp do pliku i metadanych postaci

Wiemy z przyjrzenia się kodowi modułu, że proces `COWPY` oczekuje, że zostanie mu podany plik tekstowy i wartość `character`.
Aby napisać wywołanie procesu `COWPY`, musimy tylko wiedzieć, jak wyodrębnić odpowiadający obiekt pliku i metadane z każdego elementu w kanale.

Jak to często bywa, najprostszym sposobem jest użycie operacji `map`.

Nasz kanał zawiera krotki zorganizowane jako `[meta, file]`, więc możemy uzyskać dostęp do obiektu `file` bezpośrednio, a możemy uzyskać dostęp do wartości `character` przechowywanej wewnątrz mapy meta, odwołując się do niej jako `meta.character`.

W głównym workflow wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="main.nf" linenums="34"
        // Tymczasowe: uzyskaj dostęp do pliku i postaci
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="34"
        // Tymczasowe: podejrzyj ch_languages
        ch_languages.view()
    ```

Zauważ, że używamy domknięć (takich jak `{ file -> "File: " + file }`), aby uczynić wyjście operacji `.view` bardziej czytelnym.

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

_Ścieżki plików i wartości postaci mogą pojawić się w innej kolejności w Twoim wyjściu._

To potwierdza, że jesteśmy w stanie uzyskać dostęp do pliku i postaci dla każdego elementu w kanale.

#### 3.2.3. Wywołaj proces `COWPY`

Teraz połóżmy to wszystko razem i faktycznie wywołajmy proces `COWPY` na kanale `ch_languages`.

W głównym workflow wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="main.nf" linenums="34"
        // Uruchom cowpy, aby wygenerować grafikę ASCII
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="34"
        // Tymczasowe: uzyskaj dostęp do pliku i postaci
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Widzisz, że po prostu kopiujemy dwie operacje map (minus instrukcje `.view()`) jako wejścia do wywołania procesu.
Tylko upewnij się, że nie zapomnisz przecinka między nimi!

To trochę niezgrabne, ale zobaczymy, jak to poprawić w następnej sekcji.

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

Jeśli spojrzysz w katalog results, powinieneś zobaczyć poszczególne pliki zawierające sztukę ASCII każdego powitania wymówionego przez odpowiednią postać.

??? abstract "Katalog i przykładowa zawartość pliku"

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

To pokazuje, że byliśmy w stanie użyć informacji w mapie meta do sparametryzowania polecenia w drugim kroku pipeline.

Jednakże, jak wspomniano powyżej, część kodu była trochę niezgrabna, ponieważ musieliśmy rozpakować metadane będąc jeszcze w kontekście ciała workflow.
To podejście działa dobrze dla używania niewielkiej liczby pól z mapy meta, ale słabo by się skalowało, gdybyśmy chcieli użyć znacznie więcej.

Jest inny operator o nazwie `multiMap()`, który pozwala nam to nieco usprawnić, ale nawet wtedy nie jest to idealne.

??? info "(Opcjonalnie) Alternatywna wersja z `multiMap()`"

    W przypadku, gdy się zastanawiasz, nie mogliśmy po prostu napisać pojedynczej operacji `map()`, która wyprowadza zarówno `file`, jak i `character`, ponieważ to zwróciłoby je jako krotkę.
    Musieliśmy napisać dwie oddzielne operacje `map()`, aby przekazać elementy `file` i `character` do procesu osobno.

    Technicznie istnieje inny sposób, aby to zrobić przez pojedynczą operację mapowania, używając operatora `multiMap()`, który jest w stanie emitować wiele kanałów.
    Na przykład możesz zastąpić wywołanie `COWPY` powyżej następującym kodem:

    === "Po"

        ```groovy title="main.nf" linenums="34"
            // Uruchom cowpy, aby wygenerować grafikę ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Przed"

        ```groovy title="main.nf" linenums="34"
            // Uruchom cowpy, aby wygenerować grafikę ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    To produkuje dokładnie taki sam wynik.

W obu przypadkach jest niezręcznie, że musimy robić trochę rozpakowywania na poziomie workflow.

Byłoby lepiej, gdybyśmy mogli przekazać całą mapę meta do procesu i wybrać, czego potrzebujemy, gdy tam jesteśmy.

### 3.3. Przekaż i użyj całej mapy meta

Celem mapy meta jest przecież przekazywanie wszystkich metadanych razem jako pakietu.
Jedynym powodem, dla którego nie mogliśmy tego zrobić powyżej, jest to, że proces nie jest skonfigurowany do przyjmowania mapy meta.
Ale ponieważ kontrolujemy kod procesu, możemy to zmienić.

Zmodyfikujmy proces `COWPY`, aby przyjmował strukturę krotki `[meta, file]`, którą użyliśmy w pierwszym procesie, abyśmy mogli usprawnić workflow.

W tym celu będziemy musieli zrobić trzy rzeczy:

1. Zmodyfikować definicje wejściowe modułu procesu `COWPY`
2. Zaktualizować polecenie procesu, aby używało mapy meta
3. Zaktualizować wywołanie procesu w ciele workflow

Gotowy? Zaczynajmy!

#### 3.3.1. Zmodyfikuj wejście modułu `COWPY`

Wprowadź następujące edycje do pliku modułu `cowpy.nf`:

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

To umożliwia nam użycie struktury krotki `[meta, file]`, którą omówiliśmy wcześniej w kursie.

Zauważ, że nie zaktualizowaliśmy definicji wyjściowej procesu, aby wyprowadzać mapę meta, w celu utrzymania kursu uproszczonego, ale możesz zrobić to sam jako ćwiczenie, śledząc model procesu `IDENTIFY_LANGUAGE`.

#### 3.3.2. Zaktualizuj polecenie, aby używało pola mapy meta

Cała mapa meta jest teraz dostępna wewnątrz procesu, więc możemy odwoływać się do informacji, które zawiera, bezpośrednio z bloku poleceń.

Wprowadź następujące edycje do pliku modułu `cowpy.nf`:

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

Zastąpiliśmy odniesienie do wartości `character` wcześniej przekazywanej jako samodzielne wejście wartością przechowywaną w mapie meta, do której odwołujemy się za pomocą `meta.character`.

Teraz zaktualizujmy odpowiednio wywołanie procesu.

#### 3.3.3. Zaktualizuj wywołanie procesu i uruchom je

Proces teraz oczekuje, że jego wejście użyje struktury krotki `[meta, file]`, która jest tym, co wyprowadza poprzedni proces, więc możemy po prostu przekazać cały kanał `ch_languages` do procesu `COWPY`.

Wprowadź następujące edycje do głównego workflow:

=== "Po"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Uruchom cowpy, aby wygenerować grafikę ASCII
    COWPY(ch_languages)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Uruchom cowpy, aby wygenerować grafikę ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

To znacznie upraszcza wywołanie!

Usuńmy wyniki poprzedniej egzekucji i uruchommy je:

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

Jeśli spojrzysz w katalog results, powinieneś zobaczyć te same wyjścia co wcześniej, _tj._ poszczególne pliki zawierające sztukę ASCII każdego powitania wymówionego przez odpowiednią postać.

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

Więc to produkuje te same wyniki co wcześniej przy prostszym kodzie.

Oczywiście zakłada to, że jesteś w stanie modyfikować kod procesu.
W niektórych przypadkach możesz być zmuszony polegać na istniejących procesach, których nie możesz modyfikować, co ogranicza Twoje opcje.
Dobra wiadomość, jeśli planujesz używać modułów z projektu [nf-core](https://nf-co.re/), jest taka, że wszystkie moduły nf-core są skonfigurowane do używania struktury krotki `[meta, file]` jako standardu.

### 3.4. Rozwiązywanie problemów z brakującymi wymaganymi wejściami

Wartość `character` jest wymagana, aby proces `COWPY` zakończył się pomyślnie.
Jeśli nie ustawimy dla niej wartości domyślnej w pliku konfiguracyjnym, MUSIMY podać dla niej wartość w tabeli danych.

**Co się stanie, jeśli tego nie zrobimy?**
Zależy to od tego, co zawiera tabela danych wejściowych i której wersji workflow używamy.

#### 3.4.1. Kolumna character istnieje, ale jest pusta

Powiedzmy, że usuniemy wartość character dla jednego z wpisów w naszej tabeli danych, aby symulować błąd zbierania danych:

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

Dla obu wersji workflow, których użyliśmy powyżej, klucz `character` zostanie utworzony dla wszystkich wpisów, gdy tabela danych zostanie wczytana, ale dla `sampleA` wartością będzie pusty ciąg znaków.

To spowoduje błąd.

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

Gdy Nextflow uruchamia wiersz poleceń `cowpy` dla tej próbki, `${meta.character}` jest wypełniony pustym ciągiem znaków w wierszu poleceń `cowpy`, więc narzędzie `cowpy` zgłasza błąd mówiący, że nie podano wartości dla argumentu `-c`.

#### 3.4.2. Kolumna character nie istnieje w tabeli danych

Teraz powiedzmy, że usuwamy kolumnę `character` całkowicie z naszej tabeli danych:

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

W tym przypadku klucz `character` w ogóle nie zostanie utworzony, gdy tabela danych zostanie wczytana.

##### 3.4.2.1. Wartość dostępna na poziomie workflow

Jeśli używamy wersji kodu, którą napisaliśmy w sekcji 3.2, Nextflow będzie próbował uzyskać dostęp do klucza `character` w mapie meta PRZED wywołaniem procesu `COWPY`.

Nie znajdzie żadnych elementów, które pasują do instrukcji, więc nie uruchomi `COWPY` w ogóle.

??? success "Wyjście polecenia"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Jeśli chodzi o Nextflow, ten workflow zakończył się pomyślnie!
Jednakże żadne z wyjść, których chcemy, nie zostanie wyprodukowanych.

##### 3.4.2.2. Wartość dostępna na poziomie procesu

Jeśli używamy wersji z sekcji 3.3, Nextflow przekaże całą mapę meta do procesu `COWPY` i spróbuje uruchomić polecenie.

To spowoduje błąd, ale inny w porównaniu z pierwszym przypadkiem.

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

Dzieje się tak, ponieważ `meta.character` nie istnieje, więc nasza próba uzyskania do niego dostępu zwraca `null`. W rezultacie Nextflow dosłownie wstawia `null` do wiersza poleceń, który oczywiście nie jest rozpoznawany przez narzędzie `cowpy`.

#### 3.4.3. Rozwiązania

Poza dostarczeniem wartości domyślnej jako części konfiguracji workflow, możemy zrobić dwie rzeczy, aby obsłużyć to bardziej solidnie:

1. Wdrożyć walidację wejścia do Twojego workflow, aby upewnić się, że tabela danych zawiera wszystkie wymagane informacje. Możesz znaleźć [wprowadzenie do walidacji wejścia](../hello_nf-core/05_input_validation.md) w kursie szkoleniowym Hello nf-core.

2. Jeśli chcesz upewnić się, że każdy, kto używa Twojego modułu procesu, może natychmiast zidentyfikować wymagane wejścia, możesz również uczynić wymaganą właściwość metadanych jawnym wejściem.

Oto przykład, jak by to działało.

Najpierw, na poziomie procesu, zaktualizuj definicję wejściową w następujący sposób:

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

Następnie, na poziomie workflow, użyj operacji mapowania, aby wyodrębnić właściwość `character` z metadanych i uczynić ją jawnym komponentem krotki wejściowej:

=== "Po"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

To podejście ma tę zaletę, że wyraźnie pokazuje, że `character` jest wymagana i ułatwia ponowne wdrożenie procesu w innych kontekstach.

To podkreśla ważną zasadę projektowania:

**Używaj mapy meta dla opcjonalnych, opisowych informacji, ale wyodrębniaj wymagane wartości jako jawne wejścia.**

Mapa meta jest doskonała do utrzymania czystych struktur kanałów i zapobiegania arbitralnym strukturom kanałów, ale dla obowiązkowych elementów, które są bezpośrednio przywoływane w procesie, wyodrębnienie ich jako jawnych wejść tworzy bardziej solidny i łatwiejszy w utrzymaniu kod.

### Podsumowanie

W tej sekcji nauczyłeś się, jak wykorzystywać metadane do dostosowywania wykonywania procesu, uzyskując do nich dostęp albo na poziomie workflow, albo na poziomie procesu.

---

## Ćwiczenie uzupełniające

Jeśli chciałbyś poćwiczyć używanie informacji z mapy meta wewnątrz procesu, spróbuj użyć innych informacji z mapy meta, takich jak `lang` i `lang_group`, aby dostosować sposób nazewnictwa i/lub organizowania wyjść.

Na przykład spróbuj zmodyfikować kod, aby uzyskać ten wynik:

```console title="Zawartość katalogu results"
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

---

## Podsumowanie

W tym side queście zbadałeś sposób skutecznej pracy z metadanymi w workflow Nextflow.

Ten wzorzec utrzymywania metadanych jawnych i powiązanych z danymi jest podstawową najlepszą praktyką w Nextflow, oferującą kilka zalet w porównaniu z hardcodowaniem informacji o plikach:

- Metadane plików pozostają powiązane z plikami w całym workflow
- Zachowanie procesu można dostosować dla każdego pliku
- Organizacja wyjść może odzwierciedlać metadane plików
- Informacje o plikach można rozszerzać podczas wykonywania pipeline

Zastosowanie tego wzorca w Twojej własnej pracy umożliwi Ci budowanie solidnych, łatwych w utrzymaniu workflow bioinformatycznych.

### Kluczowe wzorce

1. **Odczytywanie i strukturyzowanie metadanych:** Odczytywanie plików CSV i tworzenie zorganizowanych map metadanych, które pozostają powiązane z Twoimi plikami danych.

   ```groovy
   channel.fromPath('datasheet.csv')
     .splitCsv(header: true)
     .map { row ->
         [ [id:row.id, character:row.character], row.recording ]
     }
   ```

2. **Rozszerzanie metadanych podczas workflow** Dodawanie nowych informacji do Twoich metadanych w miarę postępu pipeline poprzez dodawanie wyjść procesu i wyprowadzanie wartości przez logikę warunkową.

   - Dodawanie nowych kluczy na podstawie wyjścia procesu

   ```groovy
   .map { meta, file, lang ->
     [ meta + [lang:lang], file ]
   }
   ```

   - Dodawanie nowych kluczy używając klauzuli warunkowej

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

3. **Dostosowywanie zachowania procesu:** Używanie metadanych wewnątrz procesu.

   ```groovy
   cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
   ```

### Dodatkowe zasoby

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Co dalej?

Wróć do [menu Side Quests](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
