# Część 3: Konfiguracja uruchamiania

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ta sekcja zbada, jak zarządzać konfiguracją pipeline'u Nextflow, aby dostosować jego zachowanie, zaadaptować go do różnych środowisk i zoptymalizować wykorzystanie zasobów _bez zmieniania ani jednej linii samego kodu workflow'u_.

Istnieje wiele sposobów, aby to osiągnąć, które można łączyć i są interpretowane zgodnie z kolejnością pierwszeństwa opisaną w dokumentacji [Configuration](https://nextflow.io/docs/latest/config.html).

W tej części kursu pokażemy najprostrzy i najpopularniejszy mechanizm pliku konfiguracyjnego, plik `nextflow.config`, z którym już się spotkałeś w sekcji o kontenerach w Części 2.

Omówimy podstawowe elementy konfiguracji Nextflow, takie jak dyrektywy procesów, executory, profile i pliki parametrów.
Ucząc się efektywnie wykorzystywać te opcje konfiguracji, możesz w pełni wykorzystać elastyczność, skalowalność i wydajność pipeline'ów Nextflow.

Aby przećwiczyć te elementy konfiguracji, będziemy uruchamiać świeżą kopię workflow'u, który ostatnio uruchamialiśmy na końcu Części 2 tego kursu szkoleniowego, przemianowaną na `3-main.nf`.

Jeśli nie znasz pipeline'u Hello lub potrzebujesz przypomnienia, zobacz [tę stronę informacyjną](../info/hello_pipeline.md).

---

## 1. Zarządzaj parametrami wejściowymi workflow

??? example "Scenariusz"

    Pobrałeś pipeline i chcesz go wielokrotnie uruchamiać z tymi samymi plikami wejściowymi i ustawieniami, ale nie chcesz za każdym razem wpisywać wszystkich parametrów.
    Lub może konfigurujesz pipeline dla kolegi, który nie czuje się komfortowo z argumentami wiersza poleceń.

Zaczniemy od aspektu konfiguracji, który jest po prostu rozszerzeniem tego, nad czym pracowaliśmy do tej pory: zarządzania parametrami wejściowymi.

Obecnie nasz workflow jest skonfigurowany do przyjmowania kilku wartości parametrów przez wiersz poleceń, zadeklarowanych w bloku `params` w samym skrypcie workflow'u.
Jeden ma wartość domyślną ustawioną jako część swojej deklaracji.

Jednak możesz chcieć ustawić wartości domyślne dla wszystkich z nich lub nadpisać istniejące ustawienie bez konieczności określania parametrów w wierszu poleceń ani modyfikowania oryginalnego pliku skryptu.

Istnieje wiele sposobów, aby to zrobić; pokażemy trzy podstawowe sposoby, które są bardzo powszechnie używane.

### 1.1. Ustaw wartości w `nextflow.config`

To najprostsze podejście, choć prawdopodobnie najmniej elastyczne, ponieważ główny plik `nextflow.config` nie jest czymś, co chcesz edytować przy każdym uruchomieniu.
Ma jednak zaletę oddzielenia kwestii _deklarowania_ parametrów w workflow (co zdecydowanie tam należy) od dostarczania _wartości domyślnych_, które lepiej pasują do pliku konfiguracyjnego.

Zróbmy to w dwóch krokach.

#### 1.1.1. Utwórz blok `params` w pliku konfiguracyjnym

Wprowadź następujące zmiany w kodzie w pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Zauważ, że nie skopiowaliśmy po prostu bloku `params` z workflow'u do pliku konfiguracyjnego.
Dla parametru `batch`, który miał już zadeklarowaną wartość domyślną, składnia jest nieco inna.
W pliku workflow'u to jest deklaracja typowana.
W konfiguracji to są przypisania wartości.

Technicznie to wystarczy do nadpisania wartości domyślnych nadal określonych w pliku workflow'u.
Możesz zmodyfikować wartość domyślną dla `batch` i uruchomić workflow, aby upewnić się, że wartość ustawiona w pliku konfiguracyjnym nadpisuje tę ustawioną w pliku workflow'u.

Ale w duchu przeniesienia konfiguracji całkowicie do pliku konfiguracyjnego, usuńmy tę wartość domyślną z pliku workflow całkowicie.

#### 1.1.2. Usuń wartość domyślną dla `batch` w pliku workflow

Wprowadź następującą zmianę w kodzie do pliku workflow'u `3-main.nf`:

=== "Po"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Przed"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Teraz sam plik workflow'u nie ustawia żadnych wartości domyślnych dla tych parametrów.

#### 1.1.3. Uruchom pipeline

Przetestujmy, czy działa poprawnie bez określania jakichkolwiek parametrów w wierszu poleceń.

```bash
nextflow run 3-main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

To nadal produkuje to samo wyjście co poprzednio.

Końcowe wyjście grafiki ASCII znajduje się w katalogu `results/3-main/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`, tak samo jak wcześniej.

??? abstract "Zawartość pliku"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Funkcjonalnie ta zmiana niczego nie zmieniła, ale koncepcyjnie jest nieco czystsze mieć wartości domyślne ustawione w pliku konfiguracyjnym.

### 1.2. Użyj pliku konfiguracyjnego specyficznego dla uruchomienia

??? example "Scenariusz"

    Chcesz eksperymentować z różnymi ustawieniami bez modyfikowania głównego pliku konfiguracyjnego.

Możesz to zrobić, tworząc nowy plik `nextflow.config` w podkatalogu, którego użyjesz jako katalog roboczy dla swoich eksperymentów.

#### 1.2.1. Utwórz katalog roboczy z pustą konfiguracją

Zacznijmy od utworzenia nowego katalogu i przejścia do niego:

```bash
mkdir -p tux-run
cd tux-run
```

Następnie utwórz pusty plik konfiguracyjny w tym katalogu:

```bash
touch nextflow.config
```

To tworzy pusty plik.

#### 1.2.2. Skonfiguruj eksperymentalną konfigurację

Teraz otwórz nowy plik i dodaj parametry, które chcesz dostosować:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Zauważ, że ścieżka do pliku wejściowego musi odzwierciedlać strukturę katalogów.

#### 1.2.3. Uruchom pipeline

Możemy teraz uruchomić nasz pipeline z poziomu nowego katalogu roboczego.
Upewnij się, że odpowiednio dostosujesz ścieżkę!

```bash
nextflow run ../3-main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

To utworzy nowy zestaw katalogów w `tux-run/`, w tym `tux-run/work/` i `tux-run/results/`.

W tym uruchomieniu Nextflow łączy `nextflow.config` w naszym bieżącym katalogu z `nextflow.config` w katalogu głównym pipeline'u i tym samym nadpisuje domyślną postać (turkey) postacią tux.

Końcowy plik wyjściowy powinien zawierać postać tux mówiącą powitania.

??? abstract "Zawartość pliku"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

To wszystko; teraz masz przestrzeń do eksperymentowania bez modyfikowania Swojej 'normalnej' konfiguracji.

!!! warning "Ostrzeżenie"

    Upewnij się, że wrócisz do poprzedniego katalogu przed przejściem do następnej sekcji!

    ```bash
    cd ..
    ```

Teraz przyjrzyjmy się innemu przydatnemu sposobowi ustawiania wartości parametrów.

### 1.3. Użyj pliku parametrów

??? example "Scenariusz"

    Musisz udostępnić dokładne parametry uruchomienia współpracownikowi lub zapisać je do publikacji.

Podejście z podkatalogiem działa świetnie do eksperymentowania, ale wymaga trochę konfiguracji i wymaga odpowiedniego dostosowania ścieżek.
Jest prostsze podejście, gdy chcesz uruchomić pipeline z określonym zestawem wartości lub umożliwić komuś innemu zrobienie tego przy minimalnym wysiłku.

Nextflow pozwala nam określić parametry za pomocą [pliku parametrów](https://nextflow.io/docs/latest/config.html#parameter-file) w formacie YAML lub JSON, co czyni bardzo wygodnym zarządzanie i dystrybucję alternatywnych zestawów wartości domyślnych, na przykład, a także wartości parametrów specyficznych dla uruchomienia.

#### 1.3.1. Zbadaj przykładowy plik parametrów

Aby to zademonstrować, dostarczamy przykładowy plik parametrów w bieżącym katalogu, o nazwie `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Ten plik parametrów zawiera parę klucz-wartość dla każdego z danych wejściowych, które chcemy określić.
Zwróć uwagę na użycie dwukropków (`:`) zamiast znaków równości (`=`), jeśli porównujesz składnię z plikiem konfiguracyjnym.
Plik config jest napisany w Groovy, podczas gdy plik parametrów jest napisany w YAML.

!!! info "Informacja"

    Dostarczamy również wersję JSON pliku parametrów jako przykład, ale nie będziemy go tutaj uruchamiać.
    Możesz spróbować tego samodzielnie.

#### 1.3.2. Uruchom pipeline

Aby uruchomić workflow z tym plikiem parametrów, po prostu dodaj `-params-file <nazwa_pliku>` do podstawowego polecenia.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

Końcowy plik wyjściowy powinien zawierać postać stegosaurus mówiącą powitania.

??? abstract "Zawartość pliku"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Użycie pliku parametrów może wydawać się przesadą, gdy masz tylko kilka parametrów do określenia, ale niektóre pipeline'y oczekują dziesiątek parametrów.
W takich przypadkach użycie pliku parametrów pozwoli nam podać wartości parametrów w czasie wykonania bez konieczności wpisywania ogromnych wierszy poleceń i bez modyfikowania skryptu workflow'u.

Ułatwia również dystrybucję zestawów parametrów współpracownikom lub jako materiał pomocniczy do publikacji, na przykład.
To czyni Twoją pracę bardziej odtwarzalną przez innych.

### Podsumowanie

Wiesz, jak wykorzystać kluczowe opcje konfiguracji do zarządzania danymi wejściowymi workflow'u.

### Co dalej?

Dowiedz się, jak zarządzać tym, gdzie i jak Twoje wyjścia workflow'u są publikowane.

---

## 2. Zarządzaj wyjściami workflow'u

??? example "Scenariusz"

    Twój pipeline publikuje wyjścia do zakodowanego na stałe katalogu, ale chcesz organizować wyniki według nazwy projektu lub eksperymentu bez edytowania kodu workflow'u za każdym razem.

Workflow, który odziedzczyliśmy, używa ścieżek dla deklaracji wyjść na poziomie workflow'u, co nie jest szczególnie elastyczne i wymaga dużo powtórzeń.

Przyjrzyjmy się kilku typowym sposobom, w jakie możesz to skonfigurować, aby było bardziej elastyczne.

### 2.1. Dostosuj nazwę katalogu `outputDir`

Każda wersja workflow'u, którą do tej pory uruchomiliśmy, publikowała Swoje wyjścia do innego podkatalogu zakodowanego na stałe w definicjach wyjść.

W Części 1 zmieniliśmy lokalizację tego podkatalogu używając flagi CLI `-output-dir`, ale to nadal jest tylko statyczny ciąg znaków.
Zamiast tego skonfigurujmy to w pliku konfiguracyjnym, gdzie możemy zdefiniować bardziej złożone dynamiczne ścieżki.
Moglibyśmy stworzyć zupełnie nowy parametr do tego celu, ale użyjmy parametru `batch`, ponieważ jest tuż pod ręką.

#### 2.1.1. Ustaw wartość dla `outputDir` w pliku konfiguracyjnym

Ścieżka, której Nextflow używa do publikowania wyjść, jest kontrolowana przez opcję `outputDir`.
Aby zmienić ścieżkę dla wszystkich wyjść, możesz ustawić wartość tej opcji w pliku konfiguracyjnym `nextflow.config`.

Dodaj następujący kod do pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

To zastąpi wbudowaną domyślną ścieżkę, `results/`, na `results_config/` plus wartość parametru `batch` jako podkatalog.

Pamiętaj, że możesz również ustawić tę opcję z wiersza poleceń używając parametru `-output-dir` w poleceniu (w skrócie `-o`), ale wtedy nie możesz użyć wartości parametru `batch`.
Użycie flagi CLI nadpisze `outputDir` w konfiguracji, jeśli jest ustawiony.

#### 2.1.2. Usuń powtarzającą się część zakodowanej na stałe ścieżki

Nadal mamy podkatalog zakodowany na stałe w opcjach wyjść, więc pozbądźmy się go teraz.

Wprowadź następujące zmiany w kodzie w pliku workflow'u:

=== "Po"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

Mogliśmy również po prostu dodać `${params.batch}` do każdej ścieżki zamiast modyfikować domyślną wartość `outputDir`, ale to jest bardziej zwięzłe.

#### 2.1.3. Uruchom pipeline

Przetestujmy, czy działa poprawnie, ustawiając nazwę batch na `outdir` z wiersza poleceń.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

To nadal produkuje to samo wyjście co poprzednio, z wyjątkiem tego, że tym razem znajdujemy nasze wyjścia w `results_config/outdir/`.

??? abstract "Zawartość katalogu"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Możesz połączyć to podejście z niestandardowymi definicjami ścieżek, aby skonstruować dowolną hierarchię katalogów.

### 2.2. Organizuj wyjścia według proces

Jednym z popularnych sposobów dalszej organizacji wyjść jest robienie tego według procesu, _tzn._ tworzenie podkatalogów dla każdego procesu uruchomionego w pipeline'ie.

#### 2.2.1. Zastąp ścieżki wyjść odniesieniem do nazw procesów

Wystarczy odwołać się do nazwy procesu jako `<proces>.name` w deklaracji ścieżki wyjścia.

Wprowadź następujące zmiany w pliku workflow'u:

=== "Po"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

To usuwa pozostałe zakodowane na stałe elementy z konfiguracji ścieżki wyjścia.

#### 2.2.2. Uruchom pipeline

Przetestujmy, czy działa poprawnie, ustawiając nazwę batch na `pnames` z wiersza poleceń.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

To nadal produkuje to samo wyjście co poprzednio, z wyjątkiem tego, że tym razem znajdujemy nasze wyjścia w `results_config/pnames/` i są one pogrupowane według procesu.

??? abstract "Zawartość katalogu"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note "Uwaga"

    Zauważ, że tutaj zatarliśmy rozróżnienie między `intermediates` a końcowymi wyjściami na najwyższym poziomie.
    Możesz mieszać i łączyć te podejścia oraz uwzględniać wiele zmiennych, na przykład ustawiając ścieżkę pierwszego wyjścia jako `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Ustaw tryb publikowania na poziomie workflow

Na koniec, w duchu zmniejszania ilości powtarzającego się kodu, możemy zastąpić deklaracje `mode` dla każdego wyjścia pojedynczą linią w konfiguracji.

#### 2.3.1. Dodaj `workflow.output.mode` do pliku konfiguracyjnego

Dodaj następujący kod do pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

Tak jak opcja `outputDir`, nadanie `workflow.output.mode` wartości w pliku konfiguracyjnym wystarczyłoby do nadpisania tego, co jest ustawione w pliku workflow'u, ale i tak usuńmy niepotrzebny kod.

#### 2.3.2. Usuń tryb wyjścia z pliku workflow'u

Wprowadź następujące zmiany w pliku workflow'u:

=== "Po"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Przed"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

To jest bardziej zwięzłe, prawda?

#### 2.3.3. Uruchom pipeline

Przetestujmy, czy działa poprawnie, ustawiając nazwę batch na `outmode` z wiersza poleceń.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

To nadal produkuje to samo wyjście co poprzednio, z wyjątkiem tego, że tym razem znajdujemy nasze wyjścia w `results_config/outmode/`.
Wszystkie nadal są właściwymi kopiami, nie symlinkami.

??? abstract "Zawartość katalogu"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Głównym powodem, dla którego nadal możesz chcieć użyć sposobu ustawiania trybu dla każdego wyjścia, jest sytuacja, gdy chcesz mieszać i łączyć w tym samym workflow, _tzn._ mieć niektóre wyjścia kopiowane, a niektóre dowiązane symbolicznie.

Jest wiele innych opcji, które możesz dostosować w ten sposób, ale mamy nadzieję, że to daje Ci poczucie zakresu opcji i jak je skutecznie wykorzystać, aby dopasować do swoich preferencji.

### Podsumowanie

Wiesz, jak kontrolować nazewnictwo i strukturę katalogów, w których publikowane są Twoje wyjścia, a także tryb publikowania wyjść workflow'u.

### Co dalej?

Dowiedz się, jak dostosować konfigurację workflow'u do Swojego środowiska obliczeniowego, zaczynając od technologii pakowania oprogramowania.

---

## 3. Wybierz technologię pakowania oprogramowania

Do tej pory przyglądaliśmy się elementom konfiguracji, które kontrolują, jak dane wejściowe wchodzą i skąd dane wyjściowe wychodzą. Teraz czas skupić się bardziej konkretnie na dostosowywaniu konfiguracji workflow'u do środowiska obliczeniowego.

Pierwszym krokiem na tej ścieżce jest określenie, skąd będą pochodzić pakiety oprogramowania, które będą uruchamiane w każdym kroku.
Czy są już zainstalowane w lokalnym środowisku obliczeniowym?
Czy musimy pobrać obrazy i uruchomić je przez system kontenerowy?
Czy musimy pobrać pakiety Conda i zbudować lokalne środowisko Conda?

W pierwszej części tego kursu szkoleniowego (Części 1-4) używaliśmy po prostu lokalnie zainstalowanego oprogramowania w naszym workflow.
Następnie w Części 5 wprowadziliśmy kontenery Docker i plik `nextflow.config`, którego użyliśmy do włączenia korzystania z kontenerów Docker.

Teraz zobaczmy, jak możemy skonfigurować alternatywną opcję pakowania oprogramowania za pomocą pliku `nextflow.config`.

### 3.1. Wyłącz Docker i włącz Conda w pliku config

??? example "Scenariusz"

    Przenosisz Swój pipeline na klaster HPC, gdzie Docker nie jest dozwolony ze względów bezpieczeństwa.
    Klaster obsługuje Singularity i Conda, więc musisz odpowiednio zmienić konfigurację.

Jak wspomniano wcześniej, Nextflow obsługuje wiele technologii kontenerowych, w tym Singularity (który jest szerzej używany na HPC), a także menedżery pakietów oprogramowania, takie jak Conda.

Możemy zmienić nasz plik konfiguracyjny, aby używał Conda zamiast Docker.
Aby to zrobić, zmieńmy wartość `docker.enabled` na `false` i dodajmy dyrektywę włączającą korzystanie z Conda:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

To pozwoli Nextflow tworzyć i wykorzystywać środowiska Conda dla procesów, które mają określone pakiety Conda.
Co oznacza, że teraz musimy dodać jeden z nich do naszego procesu `cowpy`!

### 3.2. Określ pakiet Conda w definicji process

Pobraliśmy już URI dla pakietu Conda zawierającego narzędzie `cowpy`: `conda-forge::cowpy==1.1.5`

Teraz dodajemy URI do definicji procesu `cowpy` używając dyrektywy `conda`:

=== "Po"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Przed"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Dla jasności, nie _zastępujemy_ dyrektywy `docker`, _dodajemy_ alternatywną opcję.

!!! tip "Wskazówka"

    Jest kilka różnych sposobów na uzyskanie URI dla danego pakietu conda.
    Zalecamy korzystanie z zapytania wyszukiwania [Seqera Containers](https://seqera.io/containers/), które da Ci URI, które możesz skopiować i wkleić, nawet jeśli nie planujesz tworzyć z niego kontenera.

### 3.3. Uruchom workflow, aby zweryfikować, że może używać Conda

Wypróbujmy to.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Wyjście polecenia"

    ```console title="Wyjście"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

To powinno działać bez problemu i produkować te same wyjścia co poprzednio w `results_config/conda`.

Za kulisami Nextflow pobrał pakiety Conda i utworzył środowisko, co normalnie wymaga trochę pracy; więc miło, że nie musimy nic z tego robić sami!

!!! info "Informacja"

    To działa szybko, ponieważ pakiet `cowpy` jest dość mały, ale jeśli pracujesz z dużymi pakietami, może to zająć nieco więcej czasu za pierwszym razem i możesz zobaczyć, że wyjście konsoli pozostaje 'zablokowane' przez minutę lub więcej przed zakończeniem.
    To jest normalne i wynika z dodatkowej pracy, którą Nextflow wykonuje za pierwszym razem, gdy używasz nowego pakietu.

Z naszej perspektywy wygląda to tak, jakby działało dokładnie tak samo jak uruchamianie z Docker, mimo że mechanizmy na zapleczu są nieco inne.

To oznacza, że jesteśmy gotowi do uruchamiania ze środowiskami Conda, jeśli zajdzie taka potrzeba.

??? info "Mieszanie i łączenie Docker i Conda"

    Ponieważ te dyrektywy są przypisywane dla każdego procesu, możliwe jest 'mieszanie i łączenie', _tzn._ konfigurowanie niektórych procesów w workflow'ie do uruchamiania z Docker, a innych z Conda, na przykład, jeśli infrastruktura obliczeniowa, której używasz, obsługuje obie.
    W takim przypadku włączyłbyś zarówno Docker, jak i Conda w pliku konfiguracyjnym.
    Jeśli oba są dostępne dla danego procesu, Nextflow będzie priorytetyzować kontenery.

    I jak zauważono wcześniej, Nextflow obsługuje wiele innych technologii pakowania oprogramowania i kontenerów, więc nie jesteś ograniczony tylko do tych dwóch.

### Podsumowanie

Wiesz, jak skonfigurować, jakiego pakietu oprogramowania powinien używać każdy proces i jak przełączać się między technologiami.

### Co dalej?

Dowiedz się, jak zmienić platformę wykonawczą używaną przez Nextflow do faktycznego wykonywania pracy.

---

## 4. Wybierz platformę wykonawczą

??? example "Scenariusz"

    Rozwijałeś i testowałeś Swój pipeline na laptopie, ale teraz musisz uruchomić go na tysiącach próbek.
    Twoja instytucja ma klaster HPC z harmonogramem Slurm, którego chciałbyś użyć.

Do tej pory uruchamialiśmy nasz pipeline z lokalnym executorem.
Ten wykonuje każde zadanie na maszynie, na której działa Nextflow.
Gdy Nextflow się uruchamia, sprawdza dostępne procesory i pamięć.
Jeśli wymagania zadań gotowych do uruchomienia przekraczają dostępne zasoby, Nextflow wstrzyma ostatnie z nich, dopóki jedno lub więcej wcześniejszych nie zakończy się, zwalniając niezbędne moce obliczeniowe.

Lokalny executor jest wygodny i wydajny, ale ograniczony do jednej maszyny.
Dla bardzo dużych obciążeń możesz odkryć, że Twój komputer jest wąskim gardłem, albo dlatego, że masz pojedyncze zadanie wymagające więcej mocy niż masz do dyspozycji, albo dlatego, że masz tak wiele zadań, że oczekiwanie na ich wykonanie przez jedną maszynę zajęłoby zbyt długo.

Nextflow obsługuje [wiele różnych backendów wykonawczych](https://nextflow.io/docs/latest/executor.html), w tym harmonogramy HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor i inne), a także backendy wykonywania w chmurze (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i więcej).

### 4.1. Celowanie w inny backend

Wybór executora jest ustawiany przez dyrektywę procesu o nazwie `executor`.
Domyślnie jest ustawiony na `local`, więc następująca konfiguracja jest domniemana:

```groovy title="Wbudowana konfiguracja"
process {
    executor = 'local'
}
```

Aby ustawić executor do celowania w inny backend, wystarczy określić executor, którego chcesz, używając podobnej składni jak opisano powyżej dla alokacji zasobów (zobacz dokumentację [Executors](https://nextflow.io/docs/latest/executor.html) dla wszystkich opcji).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Ostrzeżenie"

    Nie możemy tego faktycznie przetestować w środowisku szkoleniowym, ponieważ nie jest ono skonfigurowane do łączenia się z HPC.

### 4.2. Radzenie sobie ze składnią specyficzną dla backendu dla parametrów wykonania

Większość platform obliczeniowych o wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry, takie jak żądania alokacji zasobów i ograniczenia (np. liczba procesorów i pamięć) oraz nazwę kolejki zadań do użycia.

Niestety, każdy z tych systemów używa różnych technologii, składni i konfiguracji do określania, jak zadanie powinno być opisane i przesłane do odpowiedniego harmonogramu.

??? abstract "Przykłady"

    Na przykład to samo zadanie wymagające 8 procesorów i 4GB RAM do wykonania w kolejce "my-science-work" musi być wyrażone na różne sposoby w zależności od backendu.

    ```bash title="Config dla SLURM / przesyłanie przez sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config dla PBS / przesyłanie przez qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config dla SGE / przesyłanie przez qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Na szczęście Nextflow to wszystko upraszcza.
Zapewnia znormalizowaną składnię, dzięki której możesz określić odpowiednie właściwości, takie jak `cpus`, `memory` i `queue` tylko raz (zobacz dokumentację [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) dla wszystkich dostępnych opcji).
Następnie, w czasie wykonania, Nextflow użyje tych ustawień do wygenerowania odpowiednich skryptów specyficznych dla backendu na podstawie ustawienia executora.

Omówimy tę znormalizowaną składnię w następnej sekcji.

### Podsumowanie

Teraz wiesz, jak zmienić executor, aby używać różnych rodzajów infrastruktury obliczeniowej.

### Co dalej?

Dowiedz się, jak oceniać i wyrażać alokacje i ograniczenia zasobów w Nextflow.

---

## 5. Kontroluj alokacje zasobów obliczeniowych

??? example "Scenariusz"

    Twój pipeline ciągle zawodzi na klastrze, ponieważ zadania są zabijane za przekroczenie limitów pamięci.
    Lub może jesteś obciążany za zasoby, których nie używasz i chcesz zoptymalizować koszty.

Większość platform obliczeniowych o wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry alokacji zasobów, takie jak liczba procesorów i pamięć.

Domyślnie Nextflow użyje jednego procesora i 2GB pamięci dla każdego procesu.
Odpowiednie dyrektywy procesu nazywają się `cpus` i `memory`, więc następująca konfiguracja jest domniemana:

```groovy title="Wbudowana konfiguracja" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Możesz modyfikować te wartości, zarówno dla wszystkich procesów, jak i dla konkretnych nazwanych procesów, używając dodatkowych dyrektyw process w pliku konfiguracyjnym.
Nextflow przetłumaczy je na odpowiednie instrukcje dla wybranego executora.

Ale skąd wiesz, jakich wartości użyć?

### 5.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów

??? example "Scenariusz"

    Nie wiesz, ile pamięci lub CPU potrzebują Twoje procesy i chcesz uniknąć marnowania zasobów lub zabijania zadań.

Jeśli nie wiesz z góry, ile CPU i pamięci prawdopodobnie będą potrzebować Twoje procesy, możesz przeprowadzić profilowanie, co oznacza uruchomienie workflow z pewnymi domyślnymi alokacjami, zarejestrowanie, ile każdy proces zużył, i na tej podstawie oszacowanie, jak dostosować bazowe alokacje.

Wygodnie, Nextflow zawiera wbudowane narzędzia do tego i chętnie wygeneruje dla Ciebie raport na żądanie.

Aby to zrobić, dodaj `-with-report <nazwa_pliku>.html` do wiersza poleceń.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Raport jest plikiem html, który możesz pobrać i otworzyć w przeglądarce. Możesz również kliknąć go prawym przyciskiem myszy w eksploratorze plików po lewej stronie i kliknąć `Show preview`, aby wyświetlić go w środowisku szkoleniowym.

Poświęć kilka minut na przejrzenie raportu i sprawdź, czy możesz zidentyfikować jakieś możliwości dostosowania zasobów.
Upewnij się, że klikasz na karty pokazujące wyniki wykorzystania jako procent tego, co zostało przydzielone.

Zobacz dokumentację [Reports](https://nextflow.io/docs/latest/reports.html) opisującą wszystkie dostępne funkcje.

### 5.2. Ustaw alokacje zasobów dla wszystkich procesów

Profilowanie pokazuje, że procesy w naszym szkoleniowym workflow'ie są bardzo lekkie, więc zmniejszmy domyślną alokację pamięci do 1GB na proces.

Dodaj następujący kod do pliku `nextflow.config`, przed sekcją parametrów pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

To pomoże zmniejszyć ilość zużywanych zasobów obliczeniowych.

### 5.3. Ustaw alokacje zasobów dla konkretnego procesu

Jednocześnie będziemy udawać, że proces `cowpy` wymaga więcej zasobów niż inne, tylko po to, abyśmy mogli zademonstrować, jak dostosować alokacje dla indywidualnego procesu.

=== "Po"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Z tą konfiguracją wszystkie procesy będą żądać 1GB pamięci i jednego procesora (domniemana wartość domyślna), z wyjątkiem procesu `cowpy`, który będzie żądał 2GB i 2 procesorów.

!!! info "Informacja"

    Jeśli masz maszynę z małą liczbą procesorów i przydzielasz dużą liczbę na proces, możesz zobaczyć, że wywołania procesu są kolejkowane jedno za drugim.
    To dlatego, że Nextflow zapewnia, że nie żądamy więcej procesorów niż jest dostępnych.

### 5.4. Uruchom workflow ze zaktualizowaną konfiguracją

Wypróbujmy to, podając inną nazwę pliku dla raportu profilowania, abyśmy mogli porównać wydajność przed i po zmianach konfiguracji.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Prawdopodobnie nie zauważysz żadnej rzeczywistej różnicy, ponieważ to jest tak małe obciążenie, ale to jest podejście, którego użyjesz do analizy wydajności i wymagań zasobowych rzeczywistego workflow.

Jest to bardzo przydatne, gdy Twoje procesy mają różne wymagania zasobowe. Pozwala Ci odpowiednio dostosować alokacje zasobów, które ustawiasz dla każdego procesu na podstawie rzeczywistych danych, a nie domysłów.

!!! tip "Wskazówka"

    To tylko mały przedsmak tego, co możesz zrobić, aby zoptymalizować wykorzystanie zasobów.
    Sam Nextflow ma wbudowaną naprawdę sprytną [dynamiczną logikę ponownych prób](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) do ponownego uruchamiania zadań, które nie powiodły się z powodu ograniczeń zasobów.
    Dodatkowo Platforma Seqera oferuje narzędzia oparte na sztucznej inteligencji do automatycznej optymalizacji alokacji zasobów.

### 5.5. Dodaj limity zasobów

W zależności od tego, jakiego executora i infrastruktury obliczeniowej używasz, mogą istnieć pewne ograniczenia dotyczące tego, co możesz (lub musisz) przydzielić.
Na przykład Twój klaster może wymagać, abyś pozostał w określonych limitach.

Możesz użyć dyrektywy `resourceLimits`, aby ustawić odpowiednie ograniczenia. Składnia wygląda tak, gdy jest sama w bloku process:

```groovy title="Przykład składni"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow przetłumaczy te wartości na odpowiednie instrukcje w zależności od określonego executora.

Nie będziemy tego uruchamiać, ponieważ nie mamy dostępu do odpowiedniej infrastruktury w środowisku szkoleniowym.
Jednak gdybyś spróbował uruchomić workflow z alokacjami zasobów przekraczającymi te limity, a następnie sprawdził polecenie `sbatch` w pliku skryptu `.command.run`, zobaczyłbyś, że żądania, które faktycznie są wysyłane do executora, są ograniczone do wartości określonych przez `resourceLimits`.

??? info "Instytucjonalne konfiguracje referencyjne"

    Projekt nf-core skompilował [kolekcję plików konfiguracyjnych](https://nf-co.re/configs/) udostępnionych przez różne instytucje na całym świecie, obejmujących szeroką gamę executorów HPC i chmurowych.

    Te udostępnione konfiguracje są wartościowe zarówno dla osób, które tam pracują i mogą po prostu wykorzystać konfigurację swojej instytucji od razu, jak i jako model dla osób, które chcą opracować konfigurację dla własnej infrastruktury.

### Podsumowanie

Wiesz, jak wygenerować raport profilowania do oceny wykorzystania zasobów i jak modyfikować alokacje zasobów dla wszystkich procesów i/lub dla poszczególnych procesów, a także ustawiać ograniczenia zasobów do uruchamiania na HPC.

### Co dalej?

Dowiedz się, jak skonfigurować predefiniowane profile konfiguracji i przełączać się między nimi w czasie wykonania.

---

## 6. Używaj profili do przełączania między predefiniowanymi konfiguracjami

??? example "Scenariusz"

    Regularnie przełączasz się między uruchamianiem pipeline na laptopie do rozwoju i na HPC Swojej instytucji do uruchomień produkcyjnych.
    Masz dość ręcznego zmieniania ustawień konfiguracji za każdym razem, gdy zmieniasz środowiska.

Pokazaliśmy Ci wiele sposobów, w jakie możesz dostosować konfigurację pipeline w zależności od projektu, nad którym pracujesz, lub platformy obliczeniowej, której używasz.

Być może zechcesz przełączać się między alternatywnymi ustawieniami w zależności od tego, jakiej infrastruktury używasz.
Na przykład możesz rozwijać i testować na małą skalę lokalnie na laptopie, a następnie wykonywać pełnoskalowe obciążenia na HPC lub w chmurze.

Nextflow pozwala Ci skonfigurować dowolną liczbę [**profili**](https://nextflow.io/docs/latest/config.html#profiles), które opisują różne konfiguracje, które możesz następnie wybrać w czasie wykonania używając argumentu wiersza poleceń, zamiast modyfikować sam plik konfiguracyjny.

### 6.1. Utwórz profile do przełączania między lokalnym rozwojem a wykonaniem na HPC

Skonfigurujmy dwa alternatywne profile; jeden do uruchamiania małych obciążeń na zwykłym komputerze, gdzie będziemy używać kontenerów Docker, i jeden do uruchamiania na uniwersyteckim HPC z harmonogramem Slurm, gdzie będziemy używać pakietów Conda.

#### 6.1.1. Skonfiguruj profile

Dodaj następujący kod do pliku `nextflow.config`, po sekcji parametrów pipeline, ale przed ustawieniami wyjść:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

Widzisz, że dla uniwersyteckiego HPC określamy również ograniczenia zasobów.

#### 6.1.2. Uruchom workflow z profilem

Aby określić profil w wierszu poleceń Nextflow, używamy argumentu `-profile`.

Spróbujmy uruchomić workflow z konfiguracją `my_laptop`.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Jak widzisz, to pozwala nam bardzo wygodnie przełączać się między konfiguracjami w czasie wykonania.

!!! warning "Ostrzeżenie"

    Profil `univ_hpc` nie będzie działał prawidłowo w środowisku szkoleniowym, ponieważ nie mamy dostępu do harmonogramu Slurm.

Jeśli w przyszłości znajdziemy inne elementy konfiguracji, które zawsze współwystępują z tymi, możemy po prostu dodać je do odpowiedniego profilu(ów).
Możemy również tworzyć dodatkowe profile, jeśli są inne elementy konfiguracji, które chcemy zgrupować.

### 6.2. Utwórz profil parametrów testowych

??? example "Scenariusz"

    Chcesz, aby inni mogli szybko wypróbować Twój pipeline bez zbierania własnych danych wejściowych.

Profile służą nie tylko do konfiguracji infrastruktury.
Możemy ich również używać do ustawiania wartości domyślnych dla parametrów workflow, aby ułatwić innym wypróbowanie workflow bez konieczności samodzielnego zbierania odpowiednich wartości wejściowych.
Możesz rozważyć to jako alternatywę dla używania pliku parametrów.

#### 6.2.1. Skonfiguruj profil

Składnia wyrażania wartości domyślnych w tym kontekście wygląda tak, dla profilu, który nazywamy `test`:

```groovy title="Przykład składni"
    test {
        params.<parametr1>
        params.<parametr2>
        ...
    }
```

Jeśli dodamy profil testowy dla naszego workflow, blok `profiles` staje się:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Tak jak w przypadku profili konfiguracji technicznej, możesz skonfigurować wiele różnych profili określających parametry pod dowolną arbitralną nazwą.

#### 6.2.2. Uruchom workflow lokalnie z profilem testowym

Wygodnie, profile nie wykluczają się wzajemnie, więc możemy określić wiele profili w wierszu poleceń używając następującej składni `-profile <profil1>,<profil2>` (dla dowolnej liczby profili).

Jeśli łączysz profile, które ustawiają wartości dla tych samych elementów konfiguracji i są opisane w tym samym pliku konfiguracyjnym, Nextflow rozwiąże konflikt, używając tej wartości, którą odczytał jako ostatnią (_tzn._ to, co pojawia się później w pliku).
Jeśli konfliktujące ustawienia są ustawione w różnych źródłach konfiguracji, obowiązuje domyślna [kolejność pierwszeństwa](https://www.nextflow.io/docs/latest/config.html).

Spróbujmy dodać profil testowy do naszego poprzedniego polecenia:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

To użyje Docker, gdzie to możliwe, i wyprodukuje wyjścia w `results_config/test`, a tym razem postacią jest komiczny duet `dragonandcow`.

??? abstract "Zawartość pliku"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

To oznacza, że dopóki dystrybuujemy jakiekolwiek pliki danych testowych z kodem workflow'u, każdy może szybko wypróbować workflow bez konieczności dostarczania własnych danych wejściowych przez wiersz poleceń lub plik parametrów.

!!! tip "Wskazówka"

    Możemy wskazać URL-e dla większych plików, które są przechowywane zewnętrznie.
    Nextflow automatycznie je pobierze, o ile istnieje otwarte połączenie.

    Więcej szczegółów znajdziesz w Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Użyj `nextflow config`, aby zobaczyć rozwiązaną konfigurację

Jak zauważono powyżej, czasami ten sam parametr może być ustawiony na różne wartości w profilach, które chcesz połączyć.
I bardziej ogólnie, jest wiele miejsc, w których mogą być przechowywane elementy konfiguracji, a czasami te same właściwości mogą być ustawione na różne wartości w różnych miejscach.

Nextflow stosuje ustaloną [kolejność pierwszeństwa](https://nextflow.io/docs/latest/config.html#configuration-file) do rozwiązywania wszelkich konfliktów, ale może to być trudne do samodzielnego określenia.
I nawet jeśli nic nie jest w konflikcie, może być żmudne sprawdzanie wszystkich możliwych miejsc, gdzie rzeczy mogą być skonfigurowane.

Na szczęście Nextflow zawiera wygodne narzędzie o nazwie `config`, które może zautomatyzować cały ten proces za Ciebie.

Narzędzie `config` zbada całą zawartość w bieżącym katalogu roboczym, zbierze wszystkie pliki konfiguracyjne i wyprodukuje w pełni rozwiązaną konfigurację, której Nextflow użyłby do uruchomienia workflow'u.
To pozwala Ci dowiedzieć się, jakie ustawienia zostaną użyte bez konieczności uruchamiania czegokolwiek.

#### 6.3.1. Rozwiąż domyślną konfigurację

Uruchom to polecenie, aby rozwiązać konfigurację, która zostałaby zastosowana domyślnie.

```bash
nextflow config
```

??? success "Wyjście polecenia"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

To pokazuje Ci bazową konfigurację, którą otrzymujesz, jeśli nie określisz niczego dodatkowego w wierszu poleceń.

#### 6.3.2. Rozwiąż konfigurację z aktywowanymi konkretnymi ustawieniami

Jeśli podasz parametry wiersza poleceń, np. włączając jeden lub więcej profili lub ładując plik parametrów, polecenie dodatkowo je uwzględni.

```bash
nextflow config -profile my_laptop,test
```

??? success "Wyjście polecenia"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

To jest szczególnie przydatne dla złożonych projektów, które obejmują wiele warstw konfiguracji.

### Podsumowanie

Wiesz, jak używać profili do wybierania predefiniowanej konfiguracji w czasie wykonania z minimalnym wysiłkiem.
Bardziej ogólnie, wiesz, jak konfigurować wykonania workflow'u, aby pasowały do różnych platform obliczeniowych i zwiększyć odtwarzalność analiz.

### Co dalej?

Dowiedz się, jak uruchamiać pipeline'y bezpośrednio ze zdalnych repozytoriów, takich jak GitHub.

---

## 7. Uruchamiaj pipeline'y ze zdalnych repozytoriów

??? example "Scenariusz"

    Chcesz uruchomić dobrze ugruntowany pipeline, taki jak te z nf-core, bez konieczności samodzielnego pobierania i zarządzania kodem.

Do tej pory uruchamialiśmy skrypty workflow'u znajdujące się w bieżącym katalogu.
W praktyce często będziesz chciał uruchamiać pipeline'y przechowywane w zdalnych repozytoriach, takich jak GitHub.

Nextflow czyni to prostym: możesz uruchomić dowolny pipeline bezpośrednio z URL repozytorium Git bez wcześniejszego ręcznego pobierania.

### 7.1. Uruchom pipeline z GitHub

Podstawowa składnia do uruchamiania zdalnego pipeline to `nextflow run <repozytorium>`, gdzie `<repozytorium>` może być ścieżką repozytorium GitHub jak `nextflow-io/hello`, pełnym URL-em lub ścieżką do GitLab, Bitbucket lub innych usług hostingu Git.

Spróbuj uruchomić oficjalny demo pipeline Nextflow "hello":

```bash
nextflow run nextflow-io/hello
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Za pierwszym razem, gdy uruchamiasz zdalny pipeline, Nextflow pobiera go i buforuje lokalnie.
Kolejne uruchomienia używają wersji z pamięci podręcznej, chyba że wyraźnie zażądasz aktualizacji.

### 7.2. Określ wersję dla odtwarzalności

Domyślnie Nextflow uruchamia najnowszą wersję z domyślnej gałęzi.
Możesz określić konkretną wersję (tag), gałąź lub commit używając flagi `-r`:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Określanie dokładnych wersji jest niezbędne dla odtwarzalności.

### Podsumowanie

Wiesz, jak uruchamiać pipeline'y bezpośrednio z GitHub i innych zdalnych repozytoriów oraz jak określać wersje dla odtwarzalności.

### Co dalej?

Pochwal się sam!
Wiesz wszystko, co musisz wiedzieć, aby rozpocząć uruchamianie i zarządzanie pipeline'ami Nextflow.

To kończy ten kurs, ale jeśli chcesz kontynuować naukę, mamy dwie główne rekomendacje:

- Jeśli chcesz zagłębić się w tworzenie własnych pipeline'ów, zajrzyj do [Hello Nextflow](../hello_nextflow/index.md), kursu dla początkujących, który obejmuje tę samą ogólną progresję co ten, ale wchodzi w znacznie więcej szczegółów na temat kanałów i operatorów.
- Jeśli chciałbyś kontynuować naukę uruchamiania pipeline'ów Nextflow bez zagłębiania się w kod, zajrzyj do pierwszej części [Hello nf-core](../hello_nf-core/index.md), która wprowadza narzędzia do znajdowania i uruchamiania pipeline'ów z niezwykle popularnego projektu [nf-core](https://nf-co.re/).

Baw się dobrze!

---

## Quiz

<quiz>
Gdy wartości parametrów są ustawione zarówno w pliku workflow'u, jak i w `nextflow.config`, która ma pierwszeństwo?
- [ ] Wartość z pliku workflow'u
- [x] Wartość z pliku konfiguracyjnego
- [ ] Pierwsza napotkana wartość
- [ ] To powoduje błąd

Dowiedz się więcej: [1.1. Ustaw wartości w `nextflow.config`](#11-ustaw-wartosci-w-nextflowconfig)
</quiz>

<quiz>
Jaka jest różnica składniowa między ustawianiem domyślnej wartości parametru w pliku workflow'u a w pliku config?
- [ ] Używają tej samej składni
- [x] Workflow używa deklaracji typowanej (`#!groovy param: Type = value`), config używa przypisania (`#!groovy param = value`)
- [ ] Config używa deklaracji typowanej, workflow używa przypisania
- [ ] Tylko pliki config mogą ustawiać wartości domyślne

Dowiedz się więcej: [1.1. Ustaw wartości w `nextflow.config`](#11-ustaw-wartosci-w-nextflowconfig)
</quiz>

<quiz>
Jak określisz plik parametrów podczas uruchamiania workflow'u?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Dowiedz się więcej: [1.3. Użyj pliku parametrów](#13-uzyj-pliku-parametrow)
</quiz>

<quiz>
Co kontroluje opcja konfiguracji `outputDir`?
- [ ] Lokalizację katalogu roboczego
- [x] Bazową ścieżkę, gdzie publikowane są wyjścia workflow'u
- [ ] Katalog dla plików dziennika
- [ ] Lokalizację plików modułów

Dowiedz się więcej: [2.1. Dostosuj nazwę katalogu outputDir](#21-dostosuj-nazwe-katalogu-outputdir)
</quiz>

<quiz>
Jak odwołujesz się do nazwy procesu dynamicznie w konfiguracji ścieżki wyjścia?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

Dowiedz się więcej: [2.2. Organizuj wyjścia według proces](#22-organizuj-wyjscia-wedlug-process)
</quiz>

<quiz>
Jeśli zarówno Docker, jak i Conda są włączone i proces ma obie dyrektywy, która ma priorytet?
- [x] Docker (kontenery)
- [ ] Conda
- [ ] Pierwsza zdefiniowana w procesie
- [ ] To powoduje błąd

Dowiedz się więcej: [3. Wybierz technologię pakowania oprogramowania](#3-wybierz-technologie-pakowania-oprogramowania)
</quiz>

<quiz>
Jaki jest domyślny executor w Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Dowiedz się więcej: [4. Wybierz platformę wykonawczą](#4-wybierz-platforme-wykonawcza)
</quiz>

<quiz>
Jakie polecenie generuje raport wykorzystania zasobów?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Dowiedz się więcej: [5.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów](#51-uruchom-workflow-aby-wygenerowac-raport-wykorzystania-zasobow)
</quiz>

<quiz>
Jak ustawiasz wymagania zasobowe dla konkretnego procesu o nazwie `cowpy` w pliku config?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Dowiedz się więcej: [5.3. Ustaw alokacje zasobów dla konkretnego proces](#53-ustaw-alokacje-zasobow-dla-konkretnego-process)
</quiz>

<quiz>
Co robi dyrektywa `resourceLimits`?
- [ ] Ustawia minimalne wymagania zasobowe
- [ ] Przydziela zasoby do procesów
- [x] Ogranicza maksymalne zasoby, które można zażądać
- [ ] Monitoruje wykorzystanie zasobów w czasie rzeczywistym

Dowiedz się więcej: [5.5. Dodaj limity zasobów](#55-dodaj-limity-zasobow)
</quiz>

<quiz>
Jak określasz wiele profili w jednym poleceniu?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Dowiedz się więcej: [6. Używaj profili do przełączania między predefiniowanymi konfiguracjami](#6-uzywaj-profili-do-przelaczania-miedzy-predefiniowanymi-konfiguracjami)
</quiz>

<quiz>
Jakie polecenie pokazuje w pełni rozwiązaną konfigurację, której użyłby Nextflow?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Dowiedz się więcej: [6.3. Użyj `nextflow config`, aby zobaczyć rozwiązaną konfigurację](#63-uzyj-nextflow-config-aby-zobaczyc-rozwiazana-konfiguracje)
</quiz>

<quiz>
Do czego mogą być używane profile? (Wybierz wszystkie, które pasują)
- [x] Definiowanie ustawień specyficznych dla infrastruktury (executory, kontenery)
- [x] Ustawianie limitów zasobów dla różnych środowisk
- [x] Dostarczanie parametrów testowych do łatwego testowania workflow
- [ ] Definiowanie nowych procesów

Dowiedz się więcej: [6. Używaj profili do przełączania między predefiniowanymi konfiguracjami](#6-uzywaj-profili-do-przelaczania-miedzy-predefiniowanymi-konfiguracjami)
</quiz>
