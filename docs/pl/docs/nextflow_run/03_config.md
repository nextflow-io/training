# Część 3: Konfiguracja uruchomienia

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej sekcji zbadamy, jak zarządzać konfiguracją pipeline'u Nextflow'a, aby dostosować jego zachowanie, zaadaptować go do różnych środowisk i zoptymalizować wykorzystanie zasobów _bez zmiany ani jednej linii kodu workflow'a_.

Istnieje wiele sposobów, aby to osiągnąć. Można je łączyć, a ich interpretacja następuje zgodnie z kolejnością pierwszeństwa opisaną w dokumentacji [Configuration](https://nextflow.io/docs/latest/config.html).

W tej części szkolenia pokażemy Ci najprostszy i najczęściej używany mechanizm pliku konfiguracyjnego, czyli plik `nextflow.config`, który już poznałeś w sekcji o kontenerach w Części 2.

Omówimy podstawowe elementy konfiguracji Nextflow'a, takie jak dyrektywy procesów, executory, profile i pliki parametrów.
Ucząc się efektywnego wykorzystania tych opcji konfiguracyjnych, będziesz mógł w pełni wykorzystać elastyczność, skalowalność i wydajność pipeline'ów Nextflow'a.

Aby przećwiczyć te elementy konfiguracji, uruchomimy świeżą kopię workflow'a, który ostatnio uruchamialiśmy pod koniec Części 2 tego szkolenia, przemianowaną na `3-main.nf`.

Jeśli nie znasz pipeline'u Hello lub potrzebujesz przypomnienia, zobacz [tę stronę informacyjną](../info/hello_pipeline.md).

---

## 1. Zarządzanie parametrami wejściowymi workflow'a

??? example "Scenariusz"

    Pobrałeś pipeline i chcesz go wielokrotnie uruchamiać z tymi samymi plikami wejściowymi i ustawieniami, ale nie chcesz za każdym razem wpisywać wszystkich parametrów.
    Albo być może konfigurujesz pipeline dla kolegi, który nie czuje się komfortowo z argumentami wiersza poleceń.

Zaczniemy od aspektu konfiguracji, który jest po prostu rozszerzeniem tego, z czym do tej pory pracowaliśmy: zarządzania parametrami wejściowymi.

Obecnie nasz workflow jest skonfigurowany tak, aby przyjmować kilka wartości parametrów przez wiersz poleceń, zadeklarowanych w bloku `params` w samym skrypcie workflow'a.
Jeden z nich ma wartość domyślną ustawioną jako część jego deklaracji.

Możesz jednak chcieć ustawić wartości domyślne dla wszystkich z nich lub nadpisać istniejącą wartość domyślną bez konieczności określania parametrów w wierszu poleceń lub modyfikowania oryginalnego pliku skryptu.

Istnieje wiele sposobów, aby to zrobić; pokażemy Ci trzy podstawowe sposoby, które są bardzo często używane.

### 1.1. Ustawienie wartości w `nextflow.config`

To najprostsze podejście, choć prawdopodobnie najmniej elastyczne, ponieważ główny plik `nextflow.config` nie jest czymś, co chcesz edytować przy każdym uruchomieniu.
Ma jednak tę zaletę, że rozdziela kwestie _deklarowania_ parametrów w workflow'ie (co zdecydowanie tam należy) od dostarczania _wartości domyślnych_, które lepiej pasują do pliku konfiguracyjnego.

Zróbmy to w dwóch krokach.

#### 1.1.1. Utworzenie bloku `params` w pliku konfiguracyjnym

Wprowadź następujące zmiany w kodzie w pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Parametry pipeline'u
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

Zauważ, że nie skopiowaliśmy po prostu bloku `params` z workflow'a do pliku konfiguracyjnego.
Dla parametru `batch`, który miał już zadeklarowaną wartość domyślną, składnia jest nieco inna.
W pliku workflow'a jest to deklaracja z typem.
W konfiguracji są to przypisania wartości.

Technicznie rzecz biorąc, to wystarczy do nadpisania wartości domyślnych nadal określonych w pliku workflow'a.
Możesz zmodyfikować wartość domyślną dla `batch` i uruchomić workflow, aby przekonać się, że wartość ustawiona w pliku konfiguracyjnym nadpisuje tę ustawioną w pliku workflow'a.

Ale w duchu przeniesienia konfiguracji całkowicie do pliku konfiguracyjnego, usuńmy tę wartość domyślną z pliku workflow'a całkowicie.

#### 1.1.2. Usunięcie wartości domyślnej dla `batch` w pliku workflow'a

Wprowadź następującą zmianę w kodzie w pliku workflow'a `3-main.nf`:

=== "Po"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Parametry pipeline'u
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
    * Parametry pipeline'u
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Teraz sam plik workflow'a nie ustawia żadnych wartości domyślnych dla tych parametrów.

#### 1.1.3. Uruchomienie pipeline'u

Przetestujmy, czy działa poprawnie bez określania żadnych parametrów w wierszu poleceń.

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

To nadal produkuje takie samo wyjście jak poprzednio.

Końcowe wyjście ASCII art znajduje się w katalogu `results/3-main/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`, tak jak poprzednio.

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

Funkcjonalnie ta zmiana nic nie zmieniła, ale koncepcyjnie jest nieco czystsza, gdy wartości domyślne są ustawione w pliku konfiguracyjnym.

### 1.2. Użycie pliku konfiguracyjnego specyficznego dla uruchomienia

??? example "Scenariusz"

    Chcesz eksperymentować z różnymi ustawieniami bez modyfikowania głównego pliku konfiguracyjnego.

Możesz to zrobić, tworząc nowy plik `nextflow.config` w podkatalogu, którego będziesz używać jako katalogu roboczego dla swoich eksperymentów.

#### 1.2.1. Utworzenie katalogu roboczego z pustą konfiguracją

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

#### 1.2.2. Skonfigurowanie eksperymentalnej konfiguracji

Teraz otwórz nowy plik i dodaj parametry, które chcesz dostosować:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Zauważ, że ścieżka do pliku wejściowego musi odzwierciedlać strukturę katalogów.

#### 1.2.3. Uruchomienie pipeline'u

Możemy teraz uruchomić nasz pipeline z poziomu naszego nowego katalogu roboczego.
Pamiętaj, aby odpowiednio dostosować ścieżkę!

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

W tym uruchomieniu Nextflow łączy `nextflow.config` z naszego bieżącego katalogu z `nextflow.config` z katalogu głównego pipeline'u i tym samym nadpisuje domyślną postać (indyka) postacią tux.

Końcowy plik wyjściowy powinien zawierać postać tux wypowiadającą pozdrowienia.

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

To wszystko; teraz masz przestrzeń do eksperymentowania bez modyfikowania swojej „normalnej" konfiguracji.

!!! warning

    Pamiętaj, aby wrócić do poprzedniego katalogu przed przejściem do następnej sekcji!

    ```bash
    cd ..
    ```

Teraz przyjrzyjmy się innemu użytecznemu sposobowi ustawiania wartości parametrów.

### 1.3. Użycie pliku parametrów

??? example "Scenariusz"

    Musisz udostępnić dokładne parametry uruchomienia współpracownikowi lub zapisać je do publikacji.

Podejście z podkatalogiem świetnie sprawdza się do eksperymentowania, ale wymaga trochę konfiguracji i wymaga odpowiedniego dostosowania ścieżek.
Istnieje prostsze podejście, gdy chcesz uruchomić swój pipeline z określonym zestawem wartości lub umożliwić komuś innemu zrobienie tego przy minimalnym wysiłku.

Nextflow pozwala nam określić parametry za pomocą [pliku parametrów](https://nextflow.io/docs/latest/config.html#parameter-file) w formacie YAML lub JSON, co czyni bardzo wygodnym zarządzanie i dystrybucję alternatywnych zestawów wartości domyślnych, na przykład, a także wartości parametrów specyficznych dla uruchomienia.

#### 1.3.1. Przejrzenie przykładowego pliku parametrów

Aby to zademonstrować, dostarczamy przykładowy plik parametrów w bieżącym katalogu, o nazwie `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Ten plik parametrów zawiera parę klucz-wartość dla każdego z wejść, które chcemy określić.
Zauważ użycie dwukropków (`:`) zamiast znaków równości (`=`), jeśli porównasz składnię z plikiem konfiguracyjnym.
Plik konfiguracyjny jest napisany w Groovy, podczas gdy plik parametrów jest napisany w YAML.

!!! info

    Dostarczamy również wersję JSON pliku parametrów jako przykład, ale nie będziemy jej tutaj uruchamiać.
    Możesz spróbować tego samodzielnie.

#### 1.3.2. Uruchomienie pipeline'u

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

Końcowy plik wyjściowy powinien zawierać postać stegozaura wypowiadającą pozdrowienia.

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
W takich przypadkach użycie pliku parametrów pozwoli nam dostarczyć wartości parametrów w czasie wykonania bez konieczności wpisywania ogromnych wierszy poleceń i bez modyfikowania skryptu workflow'a.

Ułatwia to również dystrybucję zestawów parametrów współpracownikom lub jako informacje wspierające do publikacji, na przykład.
To sprawia, że Twoja praca jest bardziej odtwarzalna przez innych.

### Podsumowanie

Wiesz, jak wykorzystać kluczowe opcje konfiguracji do zarządzania wejściami workflow'a.

### Co dalej?

Naucz się, jak zarządzać tym, gdzie i jak publikowane są wyjścia Twojego workflow'a.

---

## 2. Zarządzanie wyjściami workflow'a

??? example "Scenariusz"

    Twój pipeline publikuje wyjścia do zakodowanego na stałe katalogu, ale chcesz organizować wyniki według nazwy projektu lub eksperymentu bez edytowania kodu workflow'a za każdym razem.

Workflow, który odziedziczyliśmy, używa ścieżek dla deklaracji wyjść na poziomie workflow'a, co nie jest zbyt elastyczne i wiąże się z dużą ilością powtórzeń.

Przyjrzyjmy się kilku powszechnym sposobom, w jakie możesz to skonfigurować, aby było bardziej elastyczne.

### 2.1. Dostosowanie nazwy katalogu `outputDir`

Każda wersja workflow'a, którą do tej pory uruchamialiśmy, publikowała swoje wyjścia do innego podkatalogu zakodowanego na stałe w definicjach wyjść.

Zmieniliśmy miejsce, w którym znajdował się ten podkatalog w Części 1, używając flagi CLI `-output-dir`, ale to nadal jest tylko statyczny ciąg znaków.
Zamiast tego skonfigurujmy to w pliku konfiguracyjnym, gdzie możemy zdefiniować bardziej złożone dynamiczne ścieżki.
Moglibyśmy utworzyć zupełnie nowy parametr do tego celu, ale użyjmy parametru `batch`, ponieważ jest już dostępny.

#### 2.1.1. Ustawienie wartości dla `outputDir` w pliku konfiguracyjnym

Ścieżka używana przez Nextflow'a do publikowania wyjść jest kontrolowana przez opcję `outputDir`.
Aby zmienić ścieżkę dla wszystkich wyjść, możesz ustawić wartość dla tej opcji w pliku konfiguracyjnym `nextflow.config`.

Dodaj następujący kod do pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Parametry pipeline'u
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Ustawienia wyjścia
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Parametry pipeline'u
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

To zastąpi wbudowaną domyślną ścieżkę `results/` ścieżką `results_config/` plus wartość parametru `batch` jako podkatalog.

Pamiętaj, że możesz również ustawić tę opcję z wiersza poleceń, używając parametru `-output-dir` w swoim poleceniu (`-o` w skrócie), ale wtedy nie mógłbyś użyć wartości parametru `batch`.
Użycie flagi CLI nadpisze `outputDir` w konfiguracji, jeśli jest ustawione.

#### 2.1.2. Usunięcie powtarzającej się części zakodowanej na stałe ścieżki

Nadal mamy podkatalog zakodowany na stałe w opcjach wyjścia, więc pozbądźmy się go teraz.

Wprowadź następujące zmiany w kodzie w pliku workflow'a:

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

Mogliśmy również po prostu dodać `${params.batch}` do każdej ścieżki zamiast modyfikować domyślny `outputDir`, ale to jest bardziej zwięzłe.

#### 2.1.3. Uruchomienie pipeline'u

Przetestujmy, czy działa poprawnie, ustawiając nazwę batcha na `outdir` z wiersza poleceń.

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

To nadal produkuje takie samo wyjście jak poprzednio, z tym że tym razem nasze wyjścia znajdują się w `results_config/outdir/`.

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

Możesz połączyć to podejście z niestandardowymi definicjami ścieżek, aby skonstruować dowolną hierarchię katalogów, jaką chcesz.

### 2.2. Organizacja wyjść według procesu

Jednym z popularnych sposobów dalszej organizacji wyjść jest robienie tego według procesu, _tzn._ tworzenie podkatalogów dla każdego procesu uruchomionego w pipeline'ie.

#### 2.2.1. Zastąpienie ścieżek wyjściowych odwołaniem do nazw procesów

Wszystko, co musisz zrobić, to odwołać się do nazwy procesu jako `<proces>.name` w deklaracji ścieżki wyjściowej.

Wprowadź następujące zmiany w pliku workflow'a:

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

To usuwa pozostałe zakodowane na stałe elementy z konfiguracji ścieżki wyjściowej.

#### 2.2.2. Uruchomienie pipeline'u

Przetestujmy, czy działa poprawnie, ustawiając nazwę batcha na `pnames` z wiersza poleceń.

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

To nadal produkuje takie samo wyjście jak poprzednio, z tym że tym razem nasze wyjścia znajdują się w `results_config/pnames/` i są pogrupowane według procesu.

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

!!! note

    Zauważ, że tutaj wymazaliśmy rozróżnienie między `intermediates` a końcowymi wyjściami znajdującymi się na najwyższym poziomie.
    Możesz mieszać i dopasowywać te podejścia, a nawet uwzględniać wiele zmiennych, na przykład ustawiając ścieżkę pierwszego wyjścia jako `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Ustawienie trybu publikowania na poziomie workflow'a

Na koniec, w duchu zmniejszenia ilości powtarzającego się kodu, możemy zastąpić deklaracje `mode` dla każdego wyjścia pojedynczą linią w konfiguracji.

#### 2.3.1. Dodanie `workflow.output.mode` do pliku konfiguracyjnego

Dodaj następujący kod do pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Ustawienia wyjścia
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Ustawienia wyjścia
    */
    outputDir = "results_config/${params.batch}"
    ```

Podobnie jak opcja `outputDir`, nadanie wartości `workflow.output.mode` w pliku konfiguracyjnym wystarczyłoby do nadpisania tego, co jest ustawione w pliku workflow'a, ale usuńmy niepotrzebny kod mimo wszystko.

#### 2.3.2. Usunięcie trybu wyjścia z pliku workflow'a

Wprowadź następujące zmiany w pliku workflow'a:

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

#### 2.3.3. Uruchomienie pipeline'u

Przetestujmy, czy działa poprawnie, ustawiając nazwę batcha na `outmode` z wiersza poleceń.

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

To nadal produkuje takie samo wyjście jak poprzednio, z tym że tym razem nasze wyjścia znajdują się w `results_config/outmode/`.
Wszystkie są nadal prawidłowymi kopiami, a nie dowiązaniami symbolicznymi.

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

Głównym powodem, dla którego możesz nadal chcieć używać sposobu ustawiania trybu dla każdego wyjścia, jest sytuacja, gdy chcesz mieszać i dopasowywać w ramach tego samego workflow'a, _tzn._ mieć niektóre wyjścia kopiowane, a niektóre jako dowiązania symboliczne.

Istnieje wiele innych opcji, które możesz dostosować w ten sposób, ale mamy nadzieję, że to daje Ci poczucie zakresu opcji i jak je efektywnie wykorzystać, aby dopasować je do Twoich preferencji.

### Podsumowanie

Wiesz, jak kontrolować nazewnictwo i strukturę katalogów, w których publikowane są Twoje wyjścia, a także tryb publikowania wyjść workflow'a.

### Co dalej?

Naucz się, jak dostosować konfigurację workflow'a do swojego środowiska obliczeniowego, zaczynając od technologii pakowania oprogramowania.

---

## 3. Wybór technologii pakowania oprogramowania

Do tej pory przyglądaliśmy się elementom konfiguracji, które kontrolują, jak wchodzą wejścia i gdzie wychodzą wyjścia. Teraz czas skupić się bardziej szczegółowo na dostosowaniu konfiguracji workflow'a do Twojego środowiska obliczeniowego.

Pierwszym krokiem na tej drodze jest określenie, skąd będą pochodzić pakiety oprogramowania, które zostaną uruchomione w każdym kroku.
Czy są już zainstalowane w lokalnym środowisku obliczeniowym?
Czy musimy pobrać obrazy i uruchomić je za pomocą systemu kontenerów?
Czy też musimy pobrać pakiety Conda i zbudować lokalne środowisko Conda?

W samej pierwszej części tego szkolenia (Części 1-4) po prostu używaliśmy lokalnie zainstalowanego oprogramowania w naszym workflow'ie.
Następnie w Części 5 wprowadziliśmy kontenery Docker i plik `nextflow.config`, którego użyliśmy do włączenia używania kontenerów Docker.

Teraz zobaczmy, jak możemy skonfigurować alternatywną opcję pakowania oprogramowania za pomocą pliku `nextflow.config`.

### 3.1. Wyłączenie Dockera i włączenie Condy w pliku konfiguracyjnym

??? example "Scenariusz"

    Przenosisz swój pipeline na klaster HPC, gdzie Docker nie jest dozwolony ze względów bezpieczeństwa.
    Klaster obsługuje Singularity i Condę, więc musisz odpowiednio zmienić swoją konfigurację.

Jak wcześniej zauważono, Nextflow obsługuje wiele technologii kontenerów, w tym Singularity (który jest szerzej używany na HPC), a także menedżery pakietów oprogramowania, takie jak Conda.

Możemy zmienić nasz plik konfiguracyjny, aby używał Condy zamiast Dockera.
Aby to zrobić, zmieńmy wartość `docker.enabled` na `false` i dodajmy dyrektywę włączającą używanie Condy:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

To pozwoli Nextflow'owi tworzyć i wykorzystywać środowiska Conda dla procesów, które mają określone pakiety Conda.
Co oznacza, że teraz musimy dodać jeden z nich do naszego procesu `cowpy`!

### 3.2. Określenie pakietu Conda w definicji procesu

Już pobraliśmy URI dla pakietu Conda zawierającego narzędzie `cowpy`: `conda-forge::cowpy==1.1.5`

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

Żeby było jasne, nie _zastępujemy_ dyrektywy `docker`, _dodajemy_ alternatywną opcję.

!!! tip

    Istnieje kilka różnych sposobów uzyskania URI dla danego pakietu conda.
    Zalecamy użycie zapytania wyszukiwania [Seqera Containers](https://seqera.io/containers/), które da Ci URI, które możesz skopiować i wkleić, nawet jeśli nie planujesz tworzyć z niego kontenera.

### 3.3. Uruchomienie workflow'a w celu weryfikacji, że może używać Condy

Wypróbujmy to.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Wyjście polecenia"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

To powinno działać bez problemu i produkować takie same wyjścia jak poprzednio w `results_config/conda`.

Za kulisami Nextflow pobrał pakiety Conda i utworzył środowisko, co normalnie wymaga trochę pracy; więc miło, że nie musimy tego robić sami!

!!! info

    To działa szybko, ponieważ pakiet `cowpy` jest dość mały, ale jeśli pracujesz z dużymi pakietami, może to zająć trochę dłużej niż zwykle za pierwszym razem i możesz zobaczyć, że wyjście konsoli pozostaje „zablokowane" przez minutę lub dłużej przed zakończeniem.
    To normalne i wynika z dodatkowej pracy, którą Nextflow wykonuje za pierwszym razem, gdy używasz nowego pakietu.

Z naszego punktu widzenia wygląda na to, że działa dokładnie tak samo jak uruchamianie z Dockerem, mimo że na backendzie mechanika jest nieco inna.

To oznacza, że jesteśmy gotowi do uruchamiania ze środowiskami Conda, jeśli zajdzie taka potrzeba.

??? info "Mieszanie i dopasowywanie Dockera i Condy"

    Ponieważ te dyrektywy są przypisywane dla każdego procesu, możliwe jest „mieszanie i dopasowywanie", _tzn._ skonfigurowanie niektórych procesów w workflow'ie do uruchamiania z Dockerem, a innych z Condą, na przykład, jeśli infrastruktura obliczeniowa, której używasz, obsługuje oba.
    W takim przypadku włączyłbyś zarówno Dockera, jak i Condę w swoim pliku konfiguracyjnym.
    Jeśli oba są dostępne dla danego procesu, Nextflow będzie priorytetowo traktować kontenery.

    I jak wcześniej zauważono, Nextflow obsługuje wiele innych technologii pakowania oprogramowania i kontenerów, więc nie jesteś ograniczony tylko do tych dwóch.

### Podsumowanie

Wiesz, jak skonfigurować, którego pakietu oprogramowania powinien używać każdy proces i jak przełączać się między technologiami.

### Co dalej?

Naucz się, jak zmienić platformę wykonawczą używaną przez Nextflow'a do faktycznego wykonywania pracy.

---

## 4. Wybór platformy wykonawczej

??? example "Scenariusz"

    Rozwijałeś i testowałeś swój pipeline na swoim laptopie, ale teraz musisz uruchomić go na tysiącach próbek.
    Twoja instytucja ma klaster HPC z harmonogramem Slurm, którego chciałbyś użyć zamiast tego.

Do tej pory uruchamialiśmy nasz pipeline z executorem lokalnym.
To wykonuje każde zadanie na maszynie, na której działa Nextflow.
Gdy Nextflow się rozpoczyna, sprawdza dostępne procesory i pamięć.
Jeśli zasoby zadań gotowych do uruchomienia przekraczają dostępne zasoby, Nextflow wstrzyma ostatnie zadania od wykonania, dopóki jedno lub więcej wcześniejszych zadań nie zakończy się, zwalniając niezbędne zasoby.

Executor lokalny jest wygodny i wydajny, ale jest ograniczony do tej pojedynczej maszyny. W przypadku bardzo dużych obciążeń możesz odkryć, że Twoja lokalna maszyna jest wąskim gardłem, albo dlatego, że masz pojedyncze zadanie, które wymaga więcej zasobów, niż masz dostępne, albo dlatego, że masz tak wiele zadań, że czekanie, aż pojedyncza maszyna je uruchomi, zajęłoby zbyt długo.

Nextflow obsługuje [wiele różnych backendów wykonawczych](https://nextflow.io/docs/latest/executor.html), w tym harmonogramy HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor i inne), a także backendy wykonawcze w chmurze (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i więcej).

### 4.1. Kierowanie na inny backend

Wybór executora jest ustawiany przez dyrektywę procesu o nazwie `executor`.
Domyślnie jest ustawiony na `local`, więc następująca konfiguracja jest domyślna:

```groovy title="Wbudowana konfiguracja"
process {
    executor = 'local'
}
```

Aby ustawić executor na kierowanie na inny backend, po prostu określiłbyś executor, którego chcesz, używając podobnej składni, jak opisano powyżej dla alokacji zasobów (zobacz [Executors](https://nextflow.io/docs/latest/executor.html) dla wszystkich opcji).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Nie możemy tego faktycznie przetestować w środowisku szkoleniowym, ponieważ nie jest ono skonfigurowane do łączenia się z HPC.

### 4.2. Radzenie sobie ze składnią specyficzną dla backendu dla parametrów wykonania

Większość platform obliczeniowych o wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry, takie jak żądania i ograniczenia alokacji zasobów (np. liczba procesorów i pamięci) oraz nazwę kolejki zadań do użycia.

Niestety, każdy z tych systemów używa różnych technologii, składni i konfiguracji do definiowania, jak zadanie powinno być zdefiniowane i przesłane do odpowiedniego harmonogramu.

??? abstract "Przykłady"

    Na przykład, to samo zadanie wymagające 8 procesorów i 4GB pamięci RAM do wykonania w kolejce "my-science-work" musi być wyrażone w następujący sposób w zależności od backendu.

    ```bash title="Konfiguracja dla SLURM / przesyłanie za pomocą sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Konfiguracja dla PBS / przesyłanie za pomocą qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Konfiguracja dla SGE / przesyłanie za pomocą qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Na szczęście Nextflow to wszystko upraszcza.
Zapewnia ustandaryzowaną składnię, dzięki czemu możesz określić odpowiednie właściwości, takie jak `cpus`, `memory` i `queue` tylko raz (zobacz [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) dla wszystkich dostępnych opcji).
Następnie, w czasie wykonania, Nextflow użyje tych ustawień do wygenerowania odpowiednich skryptów specyficznych dla backendu na podstawie ustawienia executora.

Omówimy tę ustandaryzowaną składnię w następnej sekcji.

### Podsumowanie

Teraz wiesz, jak zmienić executor, aby używać różnych rodzajów infrastruktury obliczeniowej.

### Co dalej?

Naucz się, jak oceniać i wyrażać alokacje i ograniczenia zasobów w Nextflow'ie.

---

## 5. Kontrola alokacji zasobów obliczeniowych

??? example "Scenariusz"

    Twój pipeline ciągle zawodzi na klastrze, ponieważ zadania są zabijane za przekroczenie limitów pamięci.
    Lub być może jesteś obciążany za zasoby, których nie używasz i chcesz zoptymalizować koszty.

Większość platform obliczeniowych o wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry alokacji zasobów, takie jak liczba procesorów i pamięci.

Domyślnie Nextflow użyje pojedynczego procesora i 2GB pamięci dla każdego procesu.
Odpowiednie dyrektywy procesu nazywają się `cpus` i `memory`, więc następująca konfiguracja jest domyślna:

```groovy title="Wbudowana konfiguracja" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Możesz modyfikować te wartości, albo dla wszystkich procesów, albo dla konkretnych nazwanych procesów, używając dodatkowych dyrektyw procesu w swoim pliku konfiguracyjnym.
Nextflow przetłumaczy je na odpowiednie instrukcje dla wybranego executora.

Ale skąd wiesz, jakich wartości użyć?

### 5.1. Uruchomienie workflow'a w celu wygenerowania raportu wykorzystania zasobów

??? example "Scenariusz"

    Nie wiesz, ile pamięci lub procesora potrzebują Twoje procesy i chcesz uniknąć marnowania zasobów lub zabijania zadań.

Jeśli nie wiesz z góry, ile procesora i pamięci Twoje procesy prawdopodobnie będą potrzebować, możesz wykonać profilowanie zasobów, co oznacza uruchomienie workflow'a z pewnymi domyślnymi alokacjami, zapisanie, ile każdy proces użył, a stamtąd oszacowanie, jak dostosować podstawowe alokacje.

Wygodnie, Nextflow zawiera wbudowane narzędzia do tego i chętnie wygeneruje dla Ciebie raport na żądanie.

Aby to zrobić, dodaj `-with-report <nazwa_pliku>.html` do swojego wiersza poleceń.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Raport jest plikiem html, który możesz pobrać i otworzyć w przeglądarce. Możesz również kliknąć go prawym przyciskiem myszy w eksploratorze plików po lewej stronie i kliknąć `Show preview`, aby wyświetlić go w środowisku szkoleniowym.

Poświęć kilka minut na przejrzenie raportu i zobacz, czy możesz zidentyfikować pewne możliwości dostosowania zasobów.
Upewnij się, że kliknąłeś zakładki, które pokazują wyniki wykorzystania jako procent tego, co zostało przydzielone.

Zobacz [Reports](https://nextflow.io/docs/latest/reports.html) dla dokumentacji wszystkich dostępnych funkcji.

### 5.2. Ustawienie alokacji zasobów dla wszystkich procesów

Profilowanie pokazuje, że procesy w naszym workflow'ie szkoleniowym są bardzo lekkie, więc zmniejszmy domyślną alokację pamięci do 1GB na proces.

Dodaj następujący kod do swojego pliku `nextflow.config`, przed sekcją parametrów pipeline'u:

```groovy title="nextflow.config" linenums="4"
/*
* Ustawienia procesu
*/
process {
    memory = 1.GB
}
```

To pomoże zmniejszyć ilość zasobów obliczeniowych, które zużywamy.

### 5.3. Ustawienie alokacji zasobów dla konkretnego procesu

Jednocześnie udamy, że proces `cowpy` wymaga więcej zasobów niż pozostałe, tylko po to, abyśmy mogli zademonstrować, jak dostosować alokacje dla pojedynczego procesu.

=== "Po"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Ustawienia procesu
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
    * Ustawienia procesu
    */
    process {
        memory = 1.GB
    }
    ```

Przy tej konfiguracji wszystkie procesy będą żądać 1GB pamięci i pojedynczego procesora (domyślna wartość domyślna), z wyjątkiem procesu `cowpy`, który będzie żądać 2GB i 2 procesorów.

!!! info

    Jeśli masz maszynę z niewielką liczbą procesorów i przydzielisz dużą liczbę na proces, możesz zobaczyć, że wywołania procesów są kolejkowane jedno za drugim.
    Dzieje się tak, ponieważ Nextflow zapewnia, że nie żądamy więcej procesorów, niż jest dostępnych.

### 5.4. Uruchomienie workflow'a ze zaktualizowaną konfiguracją

Wypróbujmy to, podając inną nazwę pliku dla raportu profilowania, abyśmy mogli porównać wydajność przed i po zmianach konfiguracji.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Prawdopodobnie nie zauważysz żadnej realnej różnicy, ponieważ jest to tak małe obciążenie, ale to jest podejście, którego użyłbyś do analizy wydajności i wymagań zasobowych rzeczywistego workflow'a.

Jest bardzo przydatne, gdy Twoje procesy mają różne wymagania zasobowe. Umożliwia Ci odpowiednie dopasowanie alokacji zasobów, które konfigurujesz dla każdego procesu na podstawie rzeczywistych danych, a nie zgadywania.

!!! tip

    To tylko mały przedsmak tego, co możesz zrobić, aby zoptymalizować wykorzystanie zasobów.
    Sam Nextflow ma naprawdę fajną [dynamiczną logikę ponawiania](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) wbudowaną w celu ponowienia zadań, które nie powiodły się z powodu ograniczeń zasobów.
    Dodatkowo, Seqera Platform oferuje narzędzia oparte na AI do automatycznej optymalizacji alokacji zasobów.

### 5.5. Dodanie limitów zasobów

W zależności od tego, jakiego executora obliczeniowego i infrastruktury obliczeniowej używasz, mogą istnieć pewne ograniczenia dotyczące tego, co możesz (lub musisz) przydzielić.
Na przykład, Twój klaster może wymagać, abyś pozostał w określonych limitach.

Możesz użyć dyrektywy `resourceLimits`, aby ustawić odpowiednie ograniczenia. Składnia wygląda tak, gdy jest sama w bloku procesu:

```groovy title="Przykład składni"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow przetłumaczy te wartości na odpowiednie instrukcje w zależności od executora, który określiłeś.

Nie będziemy tego uruchamiać, ponieważ nie mamy dostępu do odpowiedniej infrastruktury w środowisku szkoleniowym.
Jednak gdybyś spróbował uruchomić workflow z alokacjami zasobów, które przekraczają te limity, a następnie sprawdził polecenie `sbatch` w pliku skryptu `.command.run`, zobaczyłbyś, że żądania, które faktycznie są wysyłane do executora, są ograniczone do wartości określonych przez `resourceLimits`.

??? info "Instytucjonalne konfiguracje referencyjne"

    Projekt nf-core zebrał [kolekcję plików konfiguracyjnych](https://nf-co.re/configs/) udostępnionych przez różne instytucje na całym świecie, obejmujących szeroki zakres executorów HPC i chmurowych.

    Te udostępnione konfiguracje są cenne zarówno dla osób, które tam pracują i mogą zatem po prostu wykorzystać konfigurację swojej instytucji od razu, jak i jako model dla osób, które chcą opracować konfigurację dla własnej infrastruktury.

### Podsumowanie

Wiesz, jak wygenerować raport profilowania w celu oceny wykorzystania zasobów i jak modyfikować alokacje zasobów dla wszystkich procesów i/lub dla poszczególnych procesów, a także ustawiać ograniczenia zasobów do uruchamiania na HPC.

### Co dalej?

Naucz się, jak skonfigurować wstępnie ustawione profile konfiguracyjne i przełączać się między nimi w czasie wykonania.

---

## 6. Użycie profili do przełączania między wstępnie ustawionymi konfiguracjami

??? example "Scenariusz"

    Regularnie przełączasz się między uruchamianiem pipeline'ów na swoim laptopie do rozwoju a na HPC swojej instytucji do uruchomień produkcyjnych.
    Jesteś zmęczony ręczną zmianą ustawień konfiguracji za każdym razem, gdy przełączasz środowiska.

Pokazaliśmy Ci wiele sposobów, w jakie możesz dostosować konfigurację swojego pipeline'u w zależności od projektu, nad którym pracujesz, lub środowiska obliczeniowego, którego używasz.

Możesz chcieć przełączać się między alternatywnymi ustawieniami w zależności od tego, jakiej infrastruktury obliczeniowej używasz. Na przykład, możesz chcieć rozwijać i uruchamiać testy na małą skalę lokalnie na swoim laptopie, a następnie uruchamiać obciążenia na pełną skalę na HPC lub w chmurze.

Nextflow pozwala Ci skonfigurować dowolną liczbę [**profili**](https://nextflow.io/docs/latest/config.html#profiles), które opisują różne konfiguracje, które możesz następnie wybrać w czasie wykonania za pomocą argumentu wiersza poleceń, zamiast konieczności modyfikowania samego pliku konfiguracyjnego.

### 6.1. Utworzenie profili do przełączania między lokalnym rozwojem a wykonaniem na HPC

Skonfigurujmy dwa alternatywne profile; jeden do uruchamiania obciążeń na małą skalę na zwykłym komputerze, gdzie będziemy używać kontenerów Docker, i jeden do uruchamiania na uniwersyteckim HPC z harmonogramem Slurm, gdzie będziemy używać pakietów Conda.

#### 6.1.1. Skonfigurowanie profili

Dodaj następujący kod do swojego pliku `nextflow.config`, po sekcji parametrów pipeline'u, ale przed ustawieniami wyjścia:

```groovy title="nextflow.config" linenums="24"
/*
* Profile
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

#### 6.1.2. Uruchomienie workflow'a z profilem

Aby określić profil w naszym wierszu poleceń Nextflow'a, używamy argumentu `-profile`.

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

Jak widać, pozwala nam to bardzo wygodnie przełączać się między konfiguracjami w czasie wykonania.

!!! warning

    Profil `univ_hpc` nie będzie działał poprawnie w środowisku szkoleniowym, ponieważ nie mamy dostępu do harmonogramu Slurm.

Jeśli w przyszłości znajdziemy inne elementy konfiguracji, które zawsze współwystępują z tymi, możemy po prostu dodać je do odpowiednich profili.
Możemy również tworzyć dodatkowe profile, jeśli istnieją inne elementy konfiguracji, które chcemy zgrupować razem.

### 6.2. Utworzenie profilu parametrów testowych

??? example "Scenariusz"

    Chcesz, aby inni mogli szybko wypróbować Twój pipeline bez zbierania własnych danych wejściowych.

Profile nie służą tylko do konfiguracji infrastruktury.
Możemy również używać ich do ustawiania wartości domyślnych dla parametrów workflow'a, aby ułatwić innym wypróbowanie workflow'a bez konieczności zbierania odpowiednich wartości wejściowych samodzielnie.
Możesz to uznać za alternatywę dla używania pliku parametrów.

#### 6.2.1. Skonfigurowanie profilu

Składnia wyrażania wartości domyślnych w tym kontekście wygląda tak, dla profilu, który nazywamy `test`:

```groovy title="Przykład składni"
    test {
        params.<parametr1>
        params.<parametr2>
        ...
    }
```

Jeśli dodamy profil testowy dla naszego workflow'a, blok `profiles` staje się:

```groovy title="nextflow.config" linenums="24"
/*
* Profile
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

Podobnie jak w przypadku profili konfiguracji technicznej, możesz skonfigurować wiele różnych profili określających parametry pod dowolną arbitralną nazwą, jaką lubisz.

#### 6.2.2. Uruchomienie workflow'a lokalnie z profilem testowym

Wygodnie, profile nie wykluczają się wzajemnie, więc możemy określić wiele profili w naszym wierszu poleceń, używając następującej składni `-profile <profil1>,<profil2>` (dla dowolnej liczby profili).

Jeśli połączysz profile, które ustawiają wartości dla tych samych elementów konfiguracji i są opisane w tym samym pliku konfiguracyjnym, Nextflow rozwiąże konflikt, używając wartości, którą odczytał jako ostatnią (_tzn._ cokolwiek pojawia się później w pliku).
Jeśli sprzeczne ustawienia są ustawione w różnych źródłach konfiguracji, obowiązuje domyślna [kolejność pierwszeństwa](https://www.nextflow.io/docs/latest/config.html).

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

To użyje Dockera tam, gdzie to możliwe i wyprodukuje wyjścia w `results_config/test`, a tym razem postacią jest komediowy duet `dragonandcow`.

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

To oznacza, że dopóki dystrybuujemy jakiekolwiek pliki danych testowych z kodem workflow'a, każdy może szybko wypróbować workflow bez konieczności dostarczania własnych wejść przez wiersz poleceń lub plik parametrów.

!!! tip

    Możemy wskazywać na URL-e dla większych plików, które są przechowywane zewnętrznie.
    Nextflow pobierze je automatycznie, o ile istnieje otwarte połączenie.

    Aby uzyskać więcej szczegółów, zobacz Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Użycie `nextflow config` do zobaczenia rozwiązanej konfiguracji

Jak zauważono powyżej, czasami ten sam parametr może być ustawiony na różne wartości w profilach, które chcesz połączyć.
I bardziej ogólnie, istnieje wiele miejsc, w których elementy konfiguracji mogą być przechowywane, a czasami te same właściwości mogą być ustawione na różne wartości w różnych miejscach.

Nextflow stosuje ustaloną [kolejność pierwszeństwa](https://nextflow.io/docs/latest/config.html#configuration-file) do rozwiązywania wszelkich konfliktów, ale może to być trudne do samodzielnego określenia.
A nawet jeśli nic nie jest sprzeczne, może być żmudne sprawdzanie wszystkich możliwych miejsc, w których rzeczy mogą być skonfigurowane.

Na szczęście Nextflow zawiera wygodne narzędzie użytkowe o nazwie `config`, które może zautomatyzować cały ten proces dla Ciebie.

Narzędzie `config` zbada całą zawartość Twojego bieżącego katalogu roboczego, zbierze wszystkie pliki konfiguracyjne i wyprodukuje w pełni rozwiązaną konfigurację, której Nextflow użyłby do uruchomienia workflow'a.
Pozwala to dowiedzieć się, jakie ustawienia zostaną użyte bez konieczności uruchamiania czegokolwiek.

#### 6.3.1. Rozwiązanie domyślnej konfiguracji

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

To pokazuje Ci podstawową konfigurację, którą otrzymujesz, jeśli nie określisz niczego dodatkowego w wierszu poleceń.

#### 6.3.2. Rozwiązanie konfiguracji z aktywowanymi określonymi ustawieniami

Jeśli podasz parametry wiersza poleceń, np. włączając jeden lub więcej profili lub ładując plik parametrów, polecenie dodatkowo weźmie je pod uwagę.

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

To staje się szczególnie przydatne w przypadku złożonych projektów, które obejmują wiele warstw konfiguracji.

### Podsumowanie

Wiesz, jak używać profili do wybierania wstępnie ustawionej konfiguracji w czasie wykonania przy minimalnym wysiłku.
Bardziej ogólnie, wiesz, jak konfigurować wykonania workflow'a, aby dopasować je do różnych platform obliczeniowych i zwiększyć odtwarzalność Twoich analiz.

### Co dalej?

Naucz się, jak uruchamiać pipeline'y bezpośrednio ze zdalnych repozytoriów, takich jak GitHub.

---

## 7. Uruchamianie pipeline'ów ze zdalnych repozytoriów

??? example "Scenariusz"

    Chcesz uruchomić dobrze ugruntowany pipeline, taki jak te z nf-core, bez konieczności pobierania i zarządzania kodem samodzielnie.

Do tej pory uruchamialiśmy skrypty workflow'ów znajdujące się w bieżącym katalogu.
W praktyce często będziesz chciał uruchamiać pipeline'y przechowywane w zdalnych repozytoriach, takich jak GitHub.

Nextflow sprawia, że jest to proste: możesz uruchomić dowolny pipeline bezpośrednio z URL-a repozytorium Git bez ręcznego pobierania go najpierw.

### 7.1. Uruchomienie pipeline'u z GitHuba

Podstawowa składnia uruchamiania zdalnego pipeline'u to `nextflow run <repozytorium>`, gdzie `<repozytorium>` może być ścieżką repozytorium GitHub, taką jak `nextflow-io/hello`, pełnym URL-em lub ścieżką do GitLab, Bitbucket lub innych usług hostingowych Git.

Spróbuj uruchomić oficjalny demonstracyjny pipeline "hello" Nextflow'a:

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

### 7.2. Określenie wersji dla odtwarzalności

Domyślnie Nextflow uruchamia najnowszą wersję z domyślnej gałęzi.
Możesz określić konkretną wersję (tag), gałąź lub commit, używając flagi `-r`:

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

Wiesz, jak uruchamiać pipeline'y bezpośrednio z GitHuba i innych zdalnych repozytoriów oraz jak określać wersje dla odtwarzalności.

### Co dalej?

Pogratuluj sobie!
Wiesz wszystko, co musisz wiedzieć, aby zacząć uruchamiać i zarządzać pipeline'ami Nextflow'a.

To kończy ten kurs, ale jeśli chcesz kontynuować naukę, mamy dwie główne rekomendacje:

- Jeśli chcesz zagłębić się w rozwijanie własnych pipeline'ów, spójrz na [Hello Nextflow](../hello_nextflow/index.md), kurs dla początkujących, który obejmuje tę samą ogólną progresję co ten, ale wchodzi w znacznie więcej szczegółów na temat kanałów i operatorów.
- Jeśli chciałbyś kontynuować naukę uruchamiania pipeline'ów Nextflow'a bez zagłębiania się w kod, spójrz na pierwszą część [Hello nf-core](../hello_nf-core/index.md), która wprowadza narzędzia do znajdowania i uruchamiania pipeline'ów z niezwykle popularnego projektu [nf-core](https://nf-co.re/).

Baw się dobrze!

---

## Quiz

<quiz>
Gdy wartości parametrów są ustawione zarówno w pliku workflow'a, jak i w `nextflow.config`, która ma pierwszeństwo?
- [ ] Wartość z pliku workflow'a
- [x] Wartość z pliku konfiguracyjnego
- [ ] Pierwsza napotkana wartość
- [ ] Powoduje to błąd

Dowiedz się więcej: [1.1. Ustawienie wartości w `nextflow.config`](#11-ustawienie-wartości-w-nextflowconfig)
</quiz>

<quiz>
Jaka jest różnica składniowa między ustawianiem domyślnego parametru w pliku workflow'a a w pliku konfiguracyjnym?
- [ ] Używają tej samej składni
- [x] Workflow używa deklaracji z typem (`#!groovy param: Type = value`), konfiguracja używa przypisania (`#!groovy param = value`)
- [ ] Konfiguracja używa deklaracji z typem, workflow używa przypisania
- [ ] Tylko pliki konfiguracyjne mogą ustawiać wartości domyślne

Dowiedz się więcej: [1.1. Ustawienie wartości w `nextflow.config`](#11-ustawienie-wartości-w-nextflowconfig)
</quiz>

<quiz>
Jak określić plik parametrów podczas uruchamiania workflow'a?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Dowiedz się więcej: [1.3. Użycie pliku parametrów](#13-użycie-pliku-parametrów)
</quiz>

<quiz>
Co kontroluje opcja konfiguracji `outputDir`?
- [ ] Lokalizację katalogu roboczego
- [x] Podstawową ścieżkę, gdzie publikowane są wyjścia workflow'a
- [ ] Katalog dla plików dziennika
- [ ] Lokalizację plików modułów

Dowiedz się więcej: [2.1. Dostosowanie nazwy katalogu outputDir](#21-dostosowanie-nazwy-katalogu-outputdir)
</quiz>

<quiz>
Jak dynamicznie odwołać się do nazwy procesu w konfiguracji ścieżki wyjściowej?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<proces>.name"`
- [x] `#!groovy path { <proces>.name }`
- [ ] `@processName`

Dowiedz się więcej: [2.2. Organizacja wyjść według procesu](#22-organizacja-wyjść-według-procesu)
</quiz>

<quiz>
Jeśli zarówno Docker, jak i Conda są włączone, a proces ma obie dyrektywy, która jest priorytetowa?
- [x] Docker (kontenery)
- [ ] Conda
- [ ] Pierwsza zdefiniowana w procesie
- [ ] Powoduje to błąd

Dowiedz się więcej: [3. Wybór technologii pakowania oprogramowania](#3-wybór-technologii-pakowania-oprogramowania)
</quiz>

<quiz>
Jaki jest domyślny executor w Nextflow'ie?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Dowiedz się więcej: [4. Wybór platformy wykonawczej](#4-wybór-platformy-wykonawczej)
</quiz>

<quiz>
Jakie polecenie generuje raport wykorzystania zasobów?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Dowiedz się więcej: [5.1. Uruchomienie workflow'a w celu wygenerowania raportu wykorzystania zasobów](#51-uruchomienie-workflow'a-w-celu-wygenerowania-raportu-wykorzystania-zasobów)
</quiz>

<quiz>
Jak ustawić wymagania zasobowe dla konkretnego procesu o nazwie `cowpy` w pliku konfiguracyjnym?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Dowiedz się więcej: [5.3. Ustawienie alokacji zasobów dla konkretnego procesu](#53-ustawienie-alokacji-zasobów-dla-konkretnego-procesu)
</quiz>

<quiz>
Co robi dyrektywa `resourceLimits`?
- [ ] Ustawia minimalne wymagania zasobowe
- [ ] Przydziela zasoby do procesów
- [x] Ogranicza maksymalne zasoby, które mogą być żądane
- [ ] Monitoruje wykorzystanie zasobów w czasie rzeczywistym

Dowiedz się więcej: [5.5. Dodanie limitów zasobów](#55-dodanie-limitów-zasobów)
</quiz>

<quiz>
Jak określić wiele profili w jednym poleceniu?
- [ ] `-profile profil1 -profile profil2`
- [ ] `-profiles profil1,profil2`
- [x] `-profile profil1,profil2`
- [ ] `--profile profil1 --profile profil2`

Dowiedz się więcej: [6. Użycie profili do przełączania między wstępnie ustawionymi konfiguracjami](#6-użycie-profili-do-przełączania-między-wstępnie-ustawionymi-konfiguracjami)
</quiz>

<quiz>
Jakie polecenie pokazuje w pełni rozwiązaną konfigurację, której użyłby Nextflow?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Dowiedz się więcej: [6.3. Użycie `nextflow config` do zobaczenia rozwiązanej konfiguracji](#63-użycie-nextflow-config-do-zobaczenia-rozwiązanej-konfiguracji)
</quiz>

<quiz>
Do czego można używać profili? (Zaznacz wszystkie, które mają zastosowanie)
- [x] Definiowania ustawień specyficznych dla infrastruktury (executory, kontenery)
- [x] Ustawiania limitów zasobów dla różnych środowisk
- [x] Dostarczania parametrów testowych do łatwego testowania workflow'a
- [ ] Definiowania nowych procesów

Dowiedz się więcej: [6. Użycie profili do przełączania między wstępnie ustawionymi konfiguracjami](#6-użycie-profili-do-przełączania-między-wstępnie-ustawionymi-konfiguracjami)
</quiz>
