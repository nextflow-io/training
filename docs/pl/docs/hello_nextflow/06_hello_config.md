# Część 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) na kanale Nextflow w YouTube.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/06_hello_config.md).
///

W tej części poznasz sposoby konfiguracji Twojego pipeline'u w Nextflow, dzięki czemu będziesz mógł dostosować jego zachowanie, zaadaptować go do różnych środowisk i zoptymalizować wykorzystanie zasobów _bez zmiany ani jednej linii kodu samego workflow'a_.

Istnieje wiele sposobów, aby to osiągnąć. Można je łączyć, a ich interpretacja odbywa się zgodnie z [kolejnością pierwszeństwa](https://nextflow.io/docs/latest/config.html) opisaną w dokumentacji konfiguracji.

W tej części kursu pokażemy Ci najprostszy i najczęściej używany mechanizm pliku konfiguracyjnego: [`nextflow.config`](https://nextflow.io/docs/latest/config.html), który już poznałeś w Części 5: Hello Containers.

Omówimy kluczowe elementy konfiguracji Nextflow, takie jak dyrektywy procesów, executory, profile i pliki parametrów.
Ucząc się efektywnego wykorzystania tych opcji konfiguracyjnych, zwiększysz elastyczność, skalowalność i wydajność swoich pipeline'ów.

??? info "Jak zacząć od tej części"

    Ta część kursu zakłada, że ukończyłeś Części 1-5 kursu [Hello Nextflow](./index.md) i masz kompletny, działający pipeline.

    Jeśli zaczynasz kurs od tego momentu, musisz skopiować katalog `modules` oraz plik `nextflow.config` z rozwiązań:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    Plik `nextflow.config` zawiera linię `docker.enabled = true`, która włącza używanie kontenerów Docker.

    Jeśli nie znasz pipeline'u Hello lub potrzebujesz przypomnienia, zobacz [tę stronę informacyjną](../info/hello_pipeline.md).

---

## 0. Rozgrzewka: Uruchom `hello-config.nf`

Jako punkt wyjścia użyjemy skryptu workflow'a `hello-config.nf`.
Jest on równoważny skryptowi powstałemu w wyniku ukończenia Części 5 tego kursu, z tą różnicą, że zmieniliśmy miejsca docelowe wyjść:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Aby upewnić się, że wszystko działa, uruchom skrypt raz przed wprowadzeniem jakichkolwiek zmian:

```bash
nextflow run hello-config.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Jak poprzednio, pliki wyjściowe znajdziesz w katalogu określonym w bloku `output` (`results/hello_config/`).

??? abstract "Zawartość katalogu"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

Końcowe wyjście w postaci grafiki ASCII znajduje się w katalogu `results/hello_config/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Zawartość pliku"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Jeśli to zadziałało, jesteś gotowy, aby nauczyć się konfigurować swoje pipeline'y.

---

## 1. Zarządzaj parametrami wejściowymi workflow'a

Zaczniemy od aspektu konfiguracji, który jest po prostu rozszerzeniem tego, z czym już pracowaliśmy: zarządzania parametrami wejściowymi.

Obecnie nasz workflow jest skonfigurowany tak, aby przyjmować kilka wartości parametrów przez linię poleceń, z wartościami domyślnymi ustawionymi w bloku `params` w samym skrypcie workflow'a.
Możesz jednak chcieć nadpisać te wartości domyślne bez konieczności podawania parametrów w linii poleceń lub modyfikowania oryginalnego pliku skryptu.

Istnieje wiele sposobów, aby to zrobić; pokażemy Ci trzy podstawowe metody, które są bardzo często używane.

### 1.1. Przenieś wartości domyślne do `nextflow.config`

To najprostsza metoda, choć prawdopodobnie najmniej elastyczna, ponieważ główny plik `nextflow.config` nie jest czymś, co chcesz edytować przy każdym uruchomieniu.
Ma jednak tę zaletę, że oddziela kwestie _deklarowania_ parametrów w workflow'ie (co zdecydowanie tam należy) od dostarczania _wartości domyślnych_, które lepiej pasują do pliku konfiguracyjnego.

Zróbmy to w dwóch krokach.

#### 1.1.1. Utwórz blok `params` w pliku konfiguracyjnym

Wprowadź następujące zmiany w pliku `nextflow.config`:

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
Składnia jest nieco inna.
W pliku workflow'a są to deklaracje z typami.
W konfiguracji są to przypisania wartości.

Technicznie rzecz biorąc, to wystarczy do nadpisania wartości domyślnych nadal określonych w pliku workflow'a.
Możesz na przykład zmodyfikować postać i uruchomić workflow, aby przekonać się, że wartość ustawiona w pliku konfiguracyjnym nadpisuje tę ustawioną w pliku workflow'a.

Ale w duchu przenoszenia konfiguracji całkowicie do pliku konfiguracyjnego, usuńmy te wartości z pliku workflow'a.

#### 1.1.2. Usuń wartości z bloku `params` w pliku workflow'a

Wprowadź następujące zmiany w pliku workflow'a `hello-config.nf`:

=== "Po"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
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

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Parametry pipeline'u
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Teraz sam plik workflow'a nie ustawia żadnych wartości domyślnych dla tych parametrów.

#### 1.1.3. Uruchom pipeline

Sprawdźmy, czy działa poprawnie.

```bash
nextflow run hello-config.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

To nadal produkuje takie samo wyjście jak poprzednio.

Końcowe wyjście w postaci grafiki ASCII znajduje się w katalogu `results/hello_config/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`, tak jak poprzednio.

??? abstract "Zawartość pliku"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Funkcjonalnie to przeniesienie nic nie zmieniło, ale koncepcyjnie jest nieco czystsze, gdy wartości domyślne są ustawione w pliku konfiguracyjnym.

### 1.2. Użyj pliku konfiguracyjnego specyficznego dla uruchomienia

To świetnie, ale czasami możesz chcieć przeprowadzić tymczasowe eksperymenty z różnymi wartościami domyślnymi bez ingerencji w główny plik konfiguracyjny.
Możesz to zrobić, tworząc nowy plik `nextflow.config` w podkatalogu, którego użyjesz jako katalogu roboczego dla swoich eksperymentów.

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

Możemy teraz uruchomić nasz pipeline z poziomu naszego nowego katalogu roboczego.
Pamiętaj, aby odpowiednio dostosować ścieżkę!

```bash
nextflow run ../hello-config.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

To utworzy nowy zestaw katalogów w `tux-run/`, w tym `tux-run/work/` i `tux-run/results/`.

W tym uruchomieniu Nextflow łączy `nextflow.config` z naszego bieżącego katalogu z `nextflow.config` z katalogu głównego pipeline'u, a tym samym nadpisuje domyślną postać (turkey) postacią tux.

Końcowy plik wyjściowy powinien zawierać postać tux wypowiadającą pozdrowienia.

??? abstract "Zawartość pliku"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

    Upewnij się, że wróciłeś do poprzedniego katalogu przed przejściem do następnej sekcji!

    ```bash
    cd ..
    ```

Teraz przyjrzyjmy się innemu użytecznemu sposobowi ustawiania wartości parametrów.

### 1.3. Użyj pliku parametrów

Podejście z podkatalogiem świetnie sprawdza się przy eksperymentowaniu, ale wymaga trochę konfiguracji i wymaga odpowiedniego dostosowania ścieżek.
Istnieje prostsze podejście, gdy chcesz uruchomić swój pipeline z określonym zestawem wartości lub umożliwić komuś innemu zrobienie tego przy minimalnym wysiłku.

Nextflow pozwala nam określić parametry za pomocą [pliku parametrów](https://nextflow.io/docs/latest/config.html#params-file) w formacie YAML lub JSON, co ułatwia zarządzanie i dystrybucję alternatywnych zestawów wartości domyślnych, na przykład, a także wartości parametrów specyficznych dla uruchomienia.

#### 1.3.1. Przeanalizuj przykładowy plik parametrów

Aby to zademonstrować, udostępniamy przykładowy plik parametrów w bieżącym katalogu, o nazwie `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Ten plik parametrów zawiera parę klucz-wartość dla każdego z wejść, które chcemy określić.
Zauważ użycie dwukropków (`:`) zamiast znaków równości (`=`), jeśli porównasz składnię z plikiem konfiguracyjnym.
Plik konfiguracyjny jest napisany w Groovy, podczas gdy plik parametrów jest napisany w YAML.

!!! info

    Udostępniamy również wersję JSON pliku parametrów jako przykład, ale nie będziemy jej tutaj uruchamiać.
    Możesz spróbować jej samodzielnie.

#### 1.3.2. Uruchom pipeline

Aby uruchomić workflow z tym plikiem parametrów, po prostu dodaj `-params-file <nazwa_pliku>` do podstawowego polecenia.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Końcowy plik wyjściowy powinien zawierać postać stegozaura wypowiadającą pozdrowienia.

??? abstract "Zawartość pliku"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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

Używanie pliku parametrów może wydawać się przesadą, gdy masz tylko kilka parametrów do określenia, ale niektóre pipeline'y oczekują dziesiątek parametrów.
W takich przypadkach użycie pliku parametrów pozwoli nam dostarczyć wartości parametrów w czasie wykonania bez konieczności wpisywania ogromnych linii poleceń i bez modyfikowania skryptu workflow'a.

Ułatwia to również dystrybucję zestawów parametrów współpracownikom lub jako informacji uzupełniających do publikacji, na przykład.
To sprawia, że Twoja praca jest bardziej odtwarzalna przez innych.

### Podsumowanie

Wiesz, jak wykorzystać kluczowe opcje konfiguracji do zarządzania wejściami workflow'a.

### Co dalej?

Naucz się zarządzać tym, gdzie i jak publikowane są wyjścia Twojego workflow'a.

---

## 2. Zarządzaj wyjściami workflow'a

Do tej pory na sztywno kodowaliśmy wszystkie ścieżki dla deklaracji wyjść na poziomie workflow'a, i jak zauważyliśmy, gdy zaczęliśmy dodawać wiele wyjść, może to wiązać się z pewnym powtarzaniem.

Przyjrzyjmy się kilku powszechnym sposobom konfiguracji tego, aby było bardziej elastyczne.

### 2.1. Dostosuj katalog wyjściowy za pomocą `-output-dir`

Gdy kontrolujemy sposób organizacji naszych „opublikowanych" wyjść, mamy dwa odrębne priorytety:

- Katalog wyjściowy najwyższego poziomu
- Sposób organizacji plików w tym katalogu

Do tej pory używaliśmy domyślnego katalogu najwyższego poziomu: `results`.
Zacznijmy od dostosowania tego, używając opcji CLI `-output-dir`.

#### 2.1.1. Uruchom pipeline z `-output-dir`

Opcja `-output-dir` (skrót: `-o`) nadpisuje domyślny katalog wyjściowy (`results/`) dla wszystkich wyjść workflow'a.
To jest zalecany sposób kontrolowania ścieżki głównej, gdzie publikowane są wyjścia.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

To publikuje wyjścia do `custom-outdir-cli/` zamiast `results/`:

??? abstract "Zawartość katalogu"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Zauważ, że nadal mamy podkatalog `hello_config` z deklaracji `path` w bloku wyjściowym.
Posprzątajmy to.

#### 2.1.2. Usuń zakodowane na sztywno ścieżki z bloku wyjściowego

Prefiks `hello_config/` był zakodowany na sztywno we wcześniejszych rozdziałach, ale ponieważ teraz uczymy się elastycznie konfigurować ścieżki wyjściowe, możemy usunąć to kodowanie.
Dla wyjść, które nie potrzebują podkatalogu, możemy ustawić dyrektywę `path` na pusty ciąg lub całkowicie ją usunąć.

Wprowadź następujące zmiany w pliku workflow'a:

=== "Po"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Uruchom pipeline ponownie:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Teraz wyjścia są publikowane bezpośrednio w `custom-outdir-cli-2/`, bez podkatalogu `hello_config`:

??? abstract "Zawartość katalogu"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip

    Opcja `-output-dir` służy do kontrolowania _gdzie_ trafiają wyjścia, podczas gdy dyrektywa `path` w bloku wyjściowym kontroluje _strukturę podkatalogów_.

### 2.2. Dynamiczne ścieżki wyjściowe

Oprócz zmiany katalogu wyjściowego przez CLI, możemy również ustawić niestandardową wartość domyślną w pliku konfiguracyjnym za pomocą `outputDir`.
To pozwala nam ustawić ścieżkę katalogu dynamicznie - nie tylko używając statycznych ciągów znaków.

#### 2.2.1. Ustaw `outputDir` w pliku konfiguracyjnym

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
    outputDir = "custom-outdir-config/${params.batch}"
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

To ustawia katalog wyjściowy na `custom-outdir-config/` plus wartość parametru `batch` jako podkatalog.
Teraz możesz zmienić lokalizację wyjścia, ustawiając parametr `--batch`:

```bash
nextflow run hello-config.nf --batch my_run
```

To publikuje wyjścia do `custom-outdir-config/my_run/`.

!!! note

    Opcja CLI `-output-dir` ma pierwszeństwo przed ustawieniem konfiguracji `outputDir`.
    Jeśli jest ustawiona, opcja konfiguracyjna zostanie całkowicie zignorowana.

#### 2.2.2. Podkatalogi z nazwami batch i procesów

Możemy również ustawić deklaracje `path` wyjścia podkatalogu dynamicznie, dla każdego wyjścia osobno.

Na przykład możemy zorganizować nasze wyjścia według procesu, odwołując się do `<proces>.name` w deklaracji ścieżki wyjściowej:

=== "Po"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Możemy pójść dalej i tworzyć bardziej złożone ścieżki podkatalogów.

W powyższej edycji wymazaliśmy rozróżnienie między `intermediates` a końcowymi wyjściami na najwyższym poziomie.
Przywróćmy to, a także umieśćmy pliki w podkatalogu `params.batch`.

!!! tip

    Włączenie `params.batch` w `path` bloku wyjściowego, zamiast w `outputDir` konfiguracji, oznacza, że nie zostanie nadpisane przez `-output-dir` w CLI.

Najpierw zaktualizuj plik konfiguracyjny, aby usunąć `${params.batch}` z `outputDir` (ponieważ przenosimy to do deklaracji ścieżek):

=== "Po"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Ustawienia wyjścia
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Ustawienia wyjścia
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Następnie wprowadź następujące zmiany w pliku workflow'a:

=== "Po"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. Uruchom pipeline

Zobaczmy, jak to działa w praktyce, ustawiając zarówno `-output-dir` (lub `-o` w skrócie) na `custom-outdir-config-2`, jak i nazwę batch na `rep2` z linii poleceń:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

To publikuje wyjścia do `custom-outdir-config-2/rep2/`, z określoną ścieżką bazową _i_ podkatalogiem nazwy batch _i_ wynikami pogrupowanymi według procesu:

??? abstract "Zawartość katalogu"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Ustaw tryb publikowania na poziomie workflow'a

Na koniec, w duchu zmniejszania ilości powtarzającego się kodu, możemy zastąpić deklaracje `mode` dla każdego wyjścia pojedynczą linią w konfiguracji.

#### 2.3.1. Dodaj `workflow.output.mode` do pliku konfiguracyjnego

Dodaj następujący kod do pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Ustawienia wyjścia
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Ustawienia wyjścia
    */
    outputDir = "custom-outdir-config/"
    ```

Ustawienie `workflow.output.mode` w pliku konfiguracyjnym wystarczy, aby nadpisać to, co jest ustawione w pliku workflow'a, ale usuńmy niepotrzebny kod mimo wszystko.

#### 2.3.2. Usuń tryb wyjścia z pliku workflow'a

Wprowadź następujące zmiany w pliku workflow'a:

=== "Po"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

To jest bardziej zwięzłe, prawda?

#### 2.3.3. Uruchom pipeline

Sprawdźmy, czy działa poprawnie:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

To publikuje wyjścia do `config-output-mode/`, i nadal wszystkie są prawidłowymi kopiami, a nie dowiązaniami symbolicznymi.

??? abstract "Zawartość katalogu"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

Głównym powodem, dla którego nadal możesz chcieć używać sposobu ustawiania trybu dla każdego wyjścia, jest sytuacja, gdy chcesz mieszać i dopasowywać w ramach tego samego workflow'a, _tzn._ mieć niektóre wyjścia kopiowane, a niektóre jako dowiązania symboliczne.

Istnieje wiele innych opcji, które możesz dostosować w ten sposób, ale mamy nadzieję, że to daje Ci poczucie zakresu opcji i sposobu ich efektywnego wykorzystania zgodnie z Twoimi preferencjami.

### Podsumowanie

Wiesz, jak kontrolować nazewnictwo i strukturę katalogów, w których publikowane są Twoje wyjścia, a także tryb publikowania wyjść workflow'a.

### Co dalej?

Naucz się, jak dostosować konfigurację workflow'a do Twojego środowiska obliczeniowego, zaczynając od technologii pakowania oprogramowania.

---

## 3. Wybierz technologię pakowania oprogramowania

Do tej pory przyglądaliśmy się elementom konfiguracji, które kontrolują sposób wprowadzania danych wejściowych i miejsce, gdzie trafiają dane wyjściowe. Teraz czas skupić się bardziej szczegółowo na dostosowywaniu konfiguracji workflow'a do Twojego środowiska obliczeniowego.

Pierwszym krokiem na tej ścieżce jest określenie, skąd będą pochodzić pakiety oprogramowania, które zostaną uruchomione w każdym kroku.
Czy są już zainstalowane w lokalnym środowisku obliczeniowym?
Czy musimy pobrać obrazy i uruchomić je za pomocą systemu kontenerów?
Czy też musimy pobrać pakiety Conda i zbudować lokalne środowisko Conda?

W samej pierwszej części tego kursu (Części 1-4) po prostu używaliśmy lokalnie zainstalowanego oprogramowania w naszym workflow'ie.
Następnie w Części 5 wprowadziliśmy kontenery Docker i plik `nextflow.config`, którego użyliśmy do włączenia używania kontenerów Docker.

Teraz zobaczmy, jak możemy skonfigurować alternatywną opcję pakowania oprogramowania za pomocą pliku `nextflow.config`.

### 3.1. Wyłącz Docker i włącz Conda w pliku konfiguracyjnym

Wyobraźmy sobie, że pracujemy na klastrze HPC, a administrator nie zezwala na używanie Dockera ze względów bezpieczeństwa.
Na szczęście dla nas Nextflow obsługuje wiele innych technologii kontenerowych, w tym Singularity (która jest szerzej stosowana na HPC), oraz menedżery pakietów oprogramowania, takie jak Conda.

Możemy zmienić nasz plik konfiguracyjny, aby używał [Conda](https://nextflow.io/docs/latest/conda.html) zamiast Dockera.
Aby to zrobić, zmieńmy wartość `docker.enabled` na `false` i dodajmy dyrektywę włączającą używanie Conda:

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

### 3.2. Określ pakiet Conda w definicji procesu

Już pobraliśmy URI dla pakietu Conda zawierającego narzędzie `cowpy`: `conda-forge::cowpy==1.1.5`

Teraz dodajemy URI do definicji procesu `cowpy` za pomocą dyrektywy `conda`:

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

### 3.3. Uruchom workflow, aby sprawdzić, czy może używać Conda

Wypróbujmy to.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Wyjście polecenia"

    ```console title="Wyjście"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

To powinno działać bez problemu i produkować takie same wyjścia jak poprzednio w `custom-outdir-config/conda`.

Za kulisami Nextflow pobrał pakiety Conda i utworzył środowisko, co normalnie wymaga trochę pracy; więc miło, że nie musimy tego robić sami!

!!! note

    To działa szybko, ponieważ pakiet `cowpy` jest dość mały, ale jeśli pracujesz z dużymi pakietami, może to zająć trochę dłużej niż zwykle za pierwszym razem, i możesz zobaczyć, że wyjście konsoli pozostaje „zablokowane" przez minutę lub dłużej przed zakończeniem.
    To normalne i wynika z dodatkowej pracy, którą Nextflow wykonuje przy pierwszym użyciu nowego pakietu.

Z naszego punktu widzenia wygląda na to, że działa dokładnie tak samo jak uruchamianie z Dockerem, mimo że mechanika w tle jest nieco inna.

To oznacza, że jesteśmy gotowi do uruchamiania ze środowiskami Conda, jeśli zajdzie taka potrzeba.

??? info "Mieszanie i dopasowywanie Dockera i Condy"

    Ponieważ te dyrektywy są przypisywane dla każdego procesu, możliwe jest „mieszanie i dopasowywanie", _tzn._ skonfigurowanie niektórych procesów w workflow'ie do uruchamiania z Dockerem, a innych z Condą, na przykład, jeśli infrastruktura obliczeniowa, której używasz, obsługuje oba.
    W takim przypadku włączyłbyś zarówno Dockera, jak i Condę w swoim pliku konfiguracyjnym.
    Jeśli oba są dostępne dla danego procesu, Nextflow będzie priorytetowo traktować kontenery.

    I jak wspomniano wcześniej, Nextflow obsługuje wiele innych technologii pakowania oprogramowania i kontenerów, więc nie jesteś ograniczony tylko do tych dwóch.

### Podsumowanie

Wiesz, jak skonfigurować, który pakiet oprogramowania powinien używać każdy proces, i jak przełączać się między technologiami.

### Co dalej?

Naucz się, jak zmienić platformę wykonawczą używaną przez Nextflow do faktycznego wykonywania pracy.

---

## 4. Wybierz platformę wykonawczą

Do tej pory uruchamialiśmy nasz pipeline z executorem lokalnym.
Wykonuje on każde zadanie na maszynie, na której działa Nextflow.
Gdy Nextflow się uruchamia, sprawdza dostępne procesory i pamięć.
Jeśli zasoby zadań gotowych do uruchomienia przekraczają dostępne zasoby, Nextflow wstrzyma ostatnie zadania od wykonania, dopóki jedno lub więcej wcześniejszych zadań nie zakończy się, zwalniając niezbędne zasoby.

Executor lokalny jest wygodny i wydajny, ale jest ograniczony do tej pojedynczej maszyny. W przypadku bardzo dużych obciążeń możesz odkryć, że Twoja lokalna maszyna jest wąskim gardłem, albo dlatego, że masz pojedyncze zadanie, które wymaga więcej zasobów, niż masz dostępnych, albo dlatego, że masz tak wiele zadań, że czekanie, aż pojedyncza maszyna je uruchomi, zajęłoby zbyt długo.

Nextflow obsługuje [wiele różnych executorów](https://nextflow.io/docs/latest/executor.html), w tym harmonogramy HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor i inne), a także backendy wykonawcze w chmurze (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i więcej).

### 4.1. Kierowanie na inny backend

Wybór executora jest ustawiany przez dyrektywę procesu o nazwie `executor`.
Domyślnie jest ustawiony na `local`, więc następująca konfiguracja jest domyślna:

```groovy title="Wbudowana konfiguracja"
process {
    executor = 'local'
}
```

Aby ustawić executor na inny backend, po prostu określiłbyś executor, którego chcesz, używając podobnej składni, jak opisano powyżej dla alokacji zasobów (zobacz [dokumentację executora](https://nextflow.io/docs/latest/executor.html) dla wszystkich opcji).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Nie możemy tego faktycznie przetestować w środowisku szkoleniowym, ponieważ nie jest ono skonfigurowane do łączenia się z HPC.

### 4.2. Radzenie sobie ze składnią specyficzną dla backendu dla parametrów wykonania

Większość platform obliczeniowych o wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry, takie jak żądania i ograniczenia alokacji zasobów (np. liczba procesorów i pamięci) oraz nazwę kolejki zadań do użycia.

Niestety, każdy z tych systemów używa różnych technologii, składni i konfiguracji do definiowania sposobu, w jaki zadanie powinno być zdefiniowane i przesłane do odpowiedniego harmonogramu.

??? abstract "Przykłady"

    Na przykład to samo zadanie wymagające 8 procesorów i 4 GB pamięci RAM do wykonania w kolejce „my-science-work" musi być wyrażone w następujący sposób w zależności od backendu.

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
Zapewnia ustandaryzowaną składnię, dzięki czemu możesz określić odpowiednie właściwości, takie jak [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) i [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (zobacz [dyrektywy procesu](https://nextflow.io/docs/latest/reference/process.html#process-directives) dla innych właściwości) tylko raz.
Następnie, w czasie wykonania, Nextflow użyje tych ustawień do wygenerowania odpowiednich skryptów specyficznych dla backendu na podstawie ustawienia executora.

Omówimy tę ustandaryzowaną składnię w następnej sekcji.

### Podsumowanie

Teraz wiesz, jak zmienić executor, aby używać różnych rodzajów infrastruktury obliczeniowej.

### Co dalej?

Naucz się, jak oceniać i wyrażać alokacje i ograniczenia zasobów w Nextflow.

---

## 5. Kontroluj alokacje zasobów obliczeniowych

Większość platform obliczeniowych o wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry alokacji zasobów, takie jak liczba procesorów i pamięci.

Domyślnie Nextflow użyje pojedynczego procesora i 2 GB pamięci dla każdego procesu.
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

### 5.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów

Jeśli nie wiesz z góry, ile procesora i pamięci Twoje procesy prawdopodobnie będą potrzebować, możesz wykonać profilowanie zasobów, co oznacza uruchomienie workflow'a z pewnymi domyślnymi alokacjami, zapisanie, ile każdy proces użył, i na tej podstawie oszacowanie, jak dostosować alokacje bazowe.

Wygodnie, Nextflow zawiera wbudowane narzędzia do tego i chętnie wygeneruje dla Ciebie raport na żądanie.

Aby to zrobić, dodaj `-with-report <nazwa_pliku>.html` do swojej linii poleceń.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Raport jest plikiem html, który możesz pobrać i otworzyć w przeglądarce. Możesz również kliknąć go prawym przyciskiem myszy w eksploratorze plików po lewej stronie i kliknąć `Show preview`, aby wyświetlić go w środowisku szkoleniowym.

Poświęć kilka minut na przejrzenie raportu i sprawdź, czy możesz zidentyfikować pewne możliwości dostosowania zasobów.
Upewnij się, że kliknąłeś zakładki, które pokazują wyniki wykorzystania jako procent tego, co zostało przydzielone.

Zobacz [Raporty](https://nextflow.io/docs/latest/reports.html) dla dokumentacji wszystkich dostępnych funkcji.

### 5.2. Ustaw alokacje zasobów dla wszystkich procesów

Profilowanie pokazuje, że procesy w naszym workflow'ie szkoleniowym są bardzo lekkie, więc zmniejszmy domyślną alokację pamięci do 1 GB na proces.

Dodaj następujący kod do swojego pliku `nextflow.config`, przed sekcją parametrów pipeline'u:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Ustawienia procesu
    */
    process {
        memory = 1.GB
    }

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
    docker.enabled = false
    conda.enabled = true

    /*
    * Parametry pipeline'u
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

To pomoże zmniejszyć ilość zasobów obliczeniowych, które zużywamy.

### 5.3. Ustaw alokacje zasobów dla konkretnego procesu

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

Przy tej konfiguracji wszystkie procesy będą żądać 1 GB pamięci i pojedynczego procesora (domyślna wartość), z wyjątkiem procesu `cowpy`, który będzie żądać 2 GB i 2 procesorów.

!!! tip

    Jeśli masz maszynę z niewielką liczbą procesorów i przydzielisz dużą liczbę na proces, możesz zobaczyć, że wywołania procesów są kolejkowane jedno po drugim.
    Dzieje się tak, ponieważ Nextflow zapewnia, że nie żądamy więcej procesorów, niż jest dostępnych.

### 5.4. Uruchom workflow ze zaktualizowaną konfiguracją

Wypróbujmy to, podając inną nazwę pliku dla raportu profilowania, abyśmy mogli porównać wydajność przed i po zmianach konfiguracji.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Prawdopodobnie nie zauważysz żadnej realnej różnicy, ponieważ jest to tak małe obciążenie, ale to jest podejście, którego użyłbyś do analizy wydajności i wymagań zasobowych rzeczywistego workflow'a.

Jest to bardzo przydatne, gdy Twoje procesy mają różne wymagania zasobowe. Umożliwia Ci odpowiednie dopasowanie alokacji zasobów, które konfigurujesz dla każdego procesu, na podstawie rzeczywistych danych, a nie zgadywania.

!!! tip

    To tylko mały przedsmak tego, co możesz zrobić, aby zoptymalizować wykorzystanie zasobów.
    Sam Nextflow ma naprawdę fajną [logikę dynamicznego ponawiania](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) wbudowaną w celu ponowienia zadań, które nie powiodły się z powodu ograniczeń zasobowych.
    Dodatkowo Seqera Platform oferuje narzędzia oparte na AI do automatycznej optymalizacji alokacji zasobów.

### 5.5. Dodaj limity zasobów

W zależności od tego, jakiego executora obliczeniowego i infrastruktury obliczeniowej używasz, mogą istnieć pewne ograniczenia dotyczące tego, co możesz (lub musisz) przydzielić.
Na przykład Twój klaster może wymagać, abyś pozostał w określonych limitach.

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
Jednakże, gdybyś spróbował uruchomić workflow z alokacjami zasobów przekraczającymi te limity, a następnie sprawdził polecenie `sbatch` w pliku skryptu `.command.run`, zobaczyłbyś, że żądania faktycznie wysyłane do executora są ograniczone do wartości określonych przez `resourceLimits`.

??? info "Instytucjonalne konfiguracje referencyjne"

    Projekt nf-core zebrał [kolekcję plików konfiguracyjnych](https://nf-co.re/configs/) udostępnionych przez różne instytucje na całym świecie, obejmujących szeroki zakres executorów HPC i chmurowych.

    Te udostępnione konfiguracje są cenne zarówno dla osób, które tam pracują i mogą zatem po prostu wykorzystać konfigurację swojej instytucji od razu, jak i jako model dla osób, które chcą opracować konfigurację dla własnej infrastruktury.

### Podsumowanie

Wiesz, jak wygenerować raport profilowania, aby ocenić wykorzystanie zasobów i jak modyfikować alokacje zasobów dla wszystkich procesów i/lub dla poszczególnych procesów, a także ustawiać ograniczenia zasobów do uruchamiania na HPC.

### Co dalej?

Naucz się, jak skonfigurować wstępnie ustawione profile konfiguracyjne i przełączać się między nimi w czasie wykonania.

---

## 6. Użyj profili, aby przełączać się między wstępnie ustawionymi konfiguracjami

Pokazaliśmy Ci wiele sposobów dostosowywania konfiguracji pipeline'u w zależności od projektu, nad którym pracujesz, lub środowiska obliczeniowego, którego używasz.

Możesz chcieć przełączać się między alternatywnymi ustawieniami w zależności od tego, jakiej infrastruktury obliczeniowej używasz. Na przykład możesz chcieć rozwijać i uruchamiać testy na małą skalę lokalnie na swoim laptopie, a następnie uruchamiać obciążenia na pełną skalę na HPC lub w chmurze.

Nextflow pozwala skonfigurować dowolną liczbę [profili](https://nextflow.io/docs/latest/config.html#config-profiles), które opisują różne konfiguracje, które możesz następnie wybrać w czasie wykonania za pomocą argumentu linii poleceń, zamiast konieczności modyfikowania samego pliku konfiguracyjnego.

### 6.1. Utwórz profile do przełączania między lokalnym rozwojem a wykonaniem na HPC

Skonfigurujmy dwa alternatywne profile; jeden do uruchamiania obciążeń na małą skalę na zwykłym komputerze, gdzie użyjemy kontenerów Docker, i jeden do uruchamiania na uniwersyteckim HPC z harmonogramem Slurm, gdzie użyjemy pakietów Conda.

#### 6.1.1. Skonfiguruj profile

Dodaj następujący kod do swojego pliku `nextflow.config`, po sekcji parametrów pipeline'u, ale przed ustawieniami wyjścia:

=== "Po"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Parametry pipeline'u
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

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

    /*
    * Ustawienia wyjścia
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="15"
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
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Widzisz, że dla uniwersyteckiego HPC określamy również ograniczenia zasobów.

#### 6.1.2. Uruchom workflow z profilem

Aby określić profil w naszej linii poleceń Nextflow, używamy argumentu `-profile`.

Spróbujmy uruchomić workflow z konfiguracją `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Jak widzisz, pozwala nam to bardzo wygodnie przełączać się między konfiguracjami w czasie wykonania.

!!! warning

    Profil `univ_hpc` nie będzie działał poprawnie w środowisku szkoleniowym, ponieważ nie mamy dostępu do harmonogramu Slurm.

Jeśli w przyszłości znajdziemy inne elementy konfiguracji, które zawsze współwystępują z tymi, możemy po prostu dodać je do odpowiednich profili.
Możemy również tworzyć dodatkowe profile, jeśli istnieją inne elementy konfiguracji, które chcemy zgrupować razem.

### 6.2. Utwórz profil parametrów testowych

Profile nie służą tylko do konfiguracji infrastruktury.
Możemy ich również używać do ustawiania wartości domyślnych dla parametrów workflow'a, aby ułatwić innym wypróbowanie workflow'a bez konieczności samodzielnego zbierania odpowiednich wartości wejściowych.
Możesz to uznać za alternatywę dla używania pliku parametrów.

#### 6.2.1. Skonfiguruj profil

Składnia wyrażania wartości domyślnych w tym kontekście wygląda tak, dla profilu, który nazwiemy `test`:

```groovy title="Przykład składni"
    test {
        params.<parametr1>
        params.<parametr2>
        ...
    }
```

Jeśli dodamy profil testowy dla naszego workflow'a, blok `profiles` stanie się:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
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

Podobnie jak w przypadku profili konfiguracji technicznej, możesz skonfigurować wiele różnych profili określających parametry pod dowolną wybraną nazwą.

#### 6.2.2. Uruchom workflow lokalnie z profilem testowym

Wygodnie, profile nie wykluczają się wzajemnie, więc możemy określić wiele profili w naszej linii poleceń, używając następującej składni `-profile <profil1>,<profil2>` (dla dowolnej liczby profili).

Jeśli połączysz profile, które ustawiają wartości dla tych samych elementów konfiguracji i są opisane w tym samym pliku konfiguracyjnym, Nextflow rozwiąże konflikt, używając wartości, którą odczytał jako ostatnią (_tzn._ cokolwiek pojawia się później w pliku).
Jeśli sprzeczne ustawienia są ustawione w różnych źródłach konfiguracji, obowiązuje domyślna [kolejność pierwszeństwa](https://nextflow.io/docs/latest/config.html).

Spróbujmy dodać profil testowy do naszego poprzedniego polecenia:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

To użyje Dockera tam, gdzie to możliwe, i wyprodukuje wyjścia w `custom-outdir-config/test`, a tym razem postacią jest komediowy duet `dragonandcow`.

??? abstract "Zawartość pliku"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

To oznacza, że dopóki dystrybuujemy jakiekolwiek pliki danych testowych z kodem workflow'a, każdy może szybko wypróbować workflow bez konieczności dostarczania własnych danych wejściowych przez linię poleceń lub plik parametrów.

!!! tip

    Możemy wskazywać na URL-e dla większych plików, które są przechowywane zewnętrznie.
    Nextflow pobierze je automatycznie, o ile istnieje otwarte połączenie.

    Aby uzyskać więcej szczegółów, zobacz Side Quest [Praca z plikami](../side_quests/working_with_files.md)

### 6.3. Użyj `nextflow config`, aby zobaczyć rozwiązaną konfigurację

Jak wspomniano powyżej, czasami ten sam parametr może być ustawiony na różne wartości w profilach, które chcesz połączyć.
I bardziej ogólnie, istnieje wiele miejsc, w których elementy konfiguracji mogą być przechowywane, i czasami te same właściwości mogą być ustawione na różne wartości w różnych miejscach.

Nextflow stosuje ustaloną [kolejność pierwszeństwa](https://nextflow.io/docs/latest/config.html), aby rozwiązać wszelkie konflikty, ale może to być trudne do samodzielnego określenia.
I nawet jeśli nic nie jest sprzeczne, może być żmudne sprawdzanie wszystkich możliwych miejsc, w których rzeczy mogą być skonfigurowane.

Na szczęście Nextflow zawiera wygodne narzędzie użytkowe o nazwie `config`, które może zautomatyzować cały ten proces dla Ciebie.

Narzędzie `config` zbada całą zawartość Twojego bieżącego katalogu roboczego, zbierze wszystkie pliki konfiguracyjne i wyprodukuje w pełni rozwiązaną konfigurację, której Nextflow użyłby do uruchomienia workflow'a.
Pozwala to dowiedzieć się, jakie ustawienia zostaną użyte, bez konieczności uruchamiania czegokolwiek.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

To pokazuje Ci podstawową konfigurację, którą otrzymujesz, jeśli nie określisz niczego dodatkowego w linii poleceń.

#### 6.3.2. Rozwiąż konfigurację z określonymi aktywowanymi ustawieniami

Jeśli podasz parametry linii poleceń, np. włączając jeden lub więcej profili lub ładując plik parametrów, polecenie dodatkowo weźmie je pod uwagę.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

To staje się szczególnie przydatne w przypadku złożonych projektów, które obejmują wiele warstw konfiguracji.

### Podsumowanie

Wiesz, jak używać profili do wybierania wstępnie ustawionej konfiguracji w czasie wykonania przy minimalnym wysiłku.
Bardziej ogólnie, wiesz, jak konfigurować wykonania workflow'a, aby dostosować je do różnych platform obliczeniowych i zwiększyć odtwarzalność Twoich analiz.

### Co dalej?

Świętuj i poklepuj się po plecach! Ukończyłeś swój pierwszy kurs dla deweloperów Nextflow.

Przejdź do końcowego [podsumowania kursu](./next_steps.md), aby przejrzeć to, czego się nauczyłeś, i dowiedzieć się, co będzie dalej.

---

## Quiz

<quiz>
Jak nazywa się plik konfiguracyjny, który Nextflow automatycznie ładuje?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
Co ma pierwszeństwo, gdy ten sam parametr jest ustawiony zarówno w pliku konfiguracyjnym, jak i w linii poleceń?
- [ ] Wartość z pliku konfiguracyjnego
- [x] Wartość z linii poleceń
- [ ] Pierwsza napotkana wartość
- [ ] Żadna; powoduje to błąd

Dowiedz się więcej: [1.1. Przenieś wartości domyślne do `nextflow.config`](#11-przenieś-wartości-domyślne-do-nextflowconfig)
</quiz>

<quiz>
Czy możesz mieć zarówno Dockera, jak i Condę włączone w tej samej konfiguracji?
- [x] Tak, Nextflow może używać obu w zależności od dyrektyw procesu
- [ ] Nie, tylko jeden może być włączony na raz
- [ ] Tak, ale tylko w profilach
- [ ] Nie, wzajemnie się wykluczają
</quiz>

<quiz>
Jeśli zarówno Docker, jak i Conda są włączone, a proces ma obie dyrektywy, który ma priorytet?
- [x] Docker (kontenery)
- [ ] Conda
- [ ] Pierwszy zdefiniowany
- [ ] Powoduje to błąd

Dowiedz się więcej: [3. Wybierz technologię pakowania oprogramowania](#3-wybierz-technologię-pakowania-oprogramowania)
</quiz>

<quiz>
Jaka jest domyślna alokacja pamięci dla procesów Nextflow?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Brak limitu
</quiz>

<quiz>
Jak ustawić wymagania zasobowe dla konkretnego procesu w pliku konfiguracyjnym?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Dowiedz się więcej: [5.3. Ustaw alokacje zasobów dla konkretnego procesu](#53-ustaw-alokacje-zasobów-dla-konkretnego-procesu)
</quiz>

<quiz>
Jaka opcja linii poleceń generuje raport wykorzystania zasobów?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Dowiedz się więcej: [5.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów](#51-uruchom-workflow-aby-wygenerować-raport-wykorzystania-zasobów)
</quiz>

<quiz>
Co robi dyrektywa `resourceLimits`?
- [ ] Ustawia minimalne wymagania zasobowe
- [ ] Przydziela zasoby do procesów
- [x] Ogranicza maksymalne zasoby, które mogą być żądane
- [ ] Monitoruje wykorzystanie zasobów

Dowiedz się więcej: [5.5. Dodaj limity zasobów](#55-dodaj-limity-zasobów)
</quiz>

<quiz>
Jaki jest domyślny executor w Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Dowiedz się więcej: [4. Wybierz platformę wykonawczą](#4-wybierz-platformę-wykonawczą)
</quiz>

<quiz>
Jak określić plik parametrów podczas uruchamiania Nextflow?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Dowiedz się więcej: [1.3. Użyj pliku parametrów](#13-użyj-pliku-parametrów)
</quiz>

<quiz>
Do czego można używać profili? (Zaznacz wszystkie, które mają zastosowanie)
- [x] Definiowania ustawień specyficznych dla infrastruktury
- [x] Ustawiania limitów zasobów dla różnych środowisk
- [x] Dostarczania parametrów testowych
- [ ] Definiowania nowych procesów

Dowiedz się więcej: [6. Użyj profili, aby przełączać się między wstępnie ustawionymi konfiguracjami](#6-użyj-profili-aby-przełączać-się-między-wstępnie-ustawionymi-konfiguracjami)
</quiz>

<quiz>
Jak określić wiele profili w jednym poleceniu?
- [ ] `-profile profil1 -profile profil2`
- [ ] `-profiles profil1,profil2`
- [x] `-profile profil1,profil2`
- [ ] `--profile profil1 --profile profil2`

Dowiedz się więcej: [6. Użyj profili, aby przełączać się między wstępnie ustawionymi konfiguracjami](#6-użyj-profili-aby-przełączać-się-między-wstępnie-ustawionymi-konfiguracjami)
</quiz>
