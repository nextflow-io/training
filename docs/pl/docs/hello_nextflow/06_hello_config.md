# Część 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Obejrzyj [całą playlistę](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/06_hello_config.md).
///
-->

Ta sekcja zbada, jak skonfigurować i zarządzać konfiguracją Twojego pipeline'u Nextflow, abyś mógł dostosować jego zachowanie, zaadaptować go do różnych środowisk i zoptymalizować wykorzystanie zasobów _bez zmiany choćby jednej linii samego kodu workflow_.

Istnieje wiele sposobów, aby to zrobić. Można je używać w kombinacji i są interpretowane zgodnie z [kolejnością pierwszeństwa](https://nextflow.io/docs/latest/config.html) opisaną w dokumentacji konfiguracji.

W tej części kursu pokażemy Ci najprostszy i najczęściej używany mechanizm pliku konfiguracyjnego [`nextflow.config`](https://nextflow.io/docs/latest/config.html), który już spotkałeś w Części 5: Hello Containers.

Omówimy podstawowe komponenty konfiguracji Nextflow, takie jak dyrektywy procesów, executory, profile i pliki parametrów.
Ucząc się efektywnego wykorzystania tych opcji konfiguracyjnych, możesz zwiększyć elastyczność, skalowalność i wydajność Swoich pipeline'ów.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś Części 1-5 kursu [Hello Nextflow](./index.md) i masz kompletny działający pipeline.

    Jeśli zaczynasz kurs od tego miejsca, musisz skopiować katalog `modules` i plik `nextflow.config` z rozwiązań:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    Plik `nextflow.config` zawiera linię `docker.enabled = true`, która włącza użycie kontenerów Docker.

    Jeśli nie znasz pipeline'u Hello lub potrzebujesz przypomnienia, zobacz [tę stronę informacyjną](../info/hello_pipeline.md).

---

## 0. Rozgrzewka: Uruchom `hello-config.nf`

Użyjemy skryptu workflow `hello-config.nf` jako punktu wyjścia.
Jest on równoważny skryptowi utworzonemu podczas pracy nad Częścią 5 tego szkolenia, z tą różnicą, że zmieniliśmy miejsca docelowe wyjść:

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

Aby upewnić się, że wszystko działa, uruchom skrypt raz przed wprowadzeniem zmian:

```bash
nextflow run hello-config.nf
```

??? success "Wynik polecenia"

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

Końcowa grafika ASCII znajduje się w katalogu `results/hello_config/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`.

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

Jeśli to zadziałało, jesteś gotowy do nauki konfigurowania pipeline'ów.

---

## 1. Zarządzanie parametrami wejściowymi workflow

Zaczniemy od aspektu konfiguracji, który jest po prostu rozszerzeniem tego, nad czym pracowaliśmy do tej pory: zarządzania parametrami wejściowymi.

Obecnie nasz workflow jest skonfigurowany do przyjmowania wartości parametrów przez wiersz poleceń, z domyślnymi wartościami ustawionymi w bloku `params` w samym skrypcie workflow.
Możesz jednak chcieć nadpisać te wartości domyślne bez konieczności określania parametrów w wierszu poleceń lub modyfikowania oryginalnego pliku skryptu.

Istnieje wiele sposobów, aby to zrobić; pokażemy Ci trzy podstawowe sposoby, które są bardzo często używane.

### 1.1. Przenieś domyślne wartości do `nextflow.config`

To najprostsze podejście, choć prawdopodobnie najmniej elastyczne, ponieważ główny plik `nextflow.config` nie jest czymś, co chcesz edytować przy każdym uruchomieniu.
Ma jednak tę zaletę, że rozdziela kwestię _deklarowania_ parametrów w workflow (co zdecydowanie tam należy) od dostarczania _domyślnych wartości_, które bardziej pasują do pliku konfiguracyjnego.

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

Zauważ, że nie skopiowaliśmy po prostu bloku `params` z workflow do pliku konfiguracyjnego.
Składnia jest nieco inna.
W pliku workflow są to deklaracje typowane.
W konfiguracji są to przypisania wartości.

Technicznie to wystarczy do nadpisania domyślnych wartości nadal określonych w pliku workflow.
Możesz zmodyfikować postać, na przykład, i uruchomić workflow, aby upewnić się, że wartość ustawiona w pliku konfiguracyjnym nadpisuje tę ustawioną w pliku workflow.

Ale w duchu przeniesienia konfiguracji całkowicie do pliku konfiguracyjnego, usuńmy te wartości z pliku workflow całkowicie.

#### 1.1.2. Usuń wartości z bloku `params` w pliku workflow

Wprowadź następujące zmiany w kodzie w pliku workflow `hello-config.nf`:

=== "Po"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
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

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Teraz sam plik workflow nie ustawia żadnych domyślnych wartości dla tych parametrów.

#### 1.1.3. Uruchom pipeline

Przetestujmy, czy działa poprawnie.

```bash
nextflow run hello-config.nf
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

To nadal produkuje takie same wyjście jak poprzednio.

Końcowa grafika ASCII znajduje się w katalogu `results/hello_config/`, pod nazwą `cowpy-COLLECTED-batch-output.txt`, tak samo jak wcześniej.

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

Funkcjonalnie, to przeniesienie niczego nie zmieniło, ale koncepcyjnie jest nieco czystsze mieć domyślne wartości ustawione w pliku konfiguracyjnym.

### 1.2. Użyj pliku konfiguracyjnego specyficznego dla uruchomienia

To świetnie, ale czasami możesz chcieć przeprowadzić tymczasowe eksperymenty z innymi domyślnymi wartościami bez ingerowania w główny plik konfiguracyjny.
Możesz to zrobić, tworząc nowy plik `nextflow.config` w podkatalogu, którego użyjesz jako katalogu roboczego dla Swoich eksperymentów.

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

Możemy teraz uruchomić nasz pipeline z naszego nowego katalogu roboczego.
Upewnij się, że dostosujesz ścieżkę odpowiednio!

```bash
nextflow run ../hello-config.nf
```

??? success "Wynik polecenia"

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

W tym uruchomieniu Nextflow łączy `nextflow.config` w naszym bieżącym katalogu z `nextflow.config` w katalogu głównym pipeline'u i tym samym nadpisuje domyślną postać (turkey) postacią tux.

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

To wszystko; teraz masz przestrzeń do eksperymentowania bez modyfikowania Swojej 'normalnej' konfiguracji.

!!! warning "Ostrzeżenie"

    Upewnij się, że wrócisz do poprzedniego katalogu przed przejściem do następnej sekcji!

    ```bash
    cd ..
    ```

Teraz spójrzmy na inny użyteczny sposób ustawiania wartości parametrów.

### 1.3. Użyj pliku parametrów

Podejście z podkatalogiem świetnie sprawdza się do eksperymentowania, ale wymaga trochę konfiguracji i wymaga dostosowania ścieżek.
Jest prostsze podejście, gdy chcesz uruchomić pipeline z konkretnym zestawem wartości lub umożliwić komuś innemu zrobienie tego z minimalnym wysiłkiem.

Nextflow pozwala nam określić parametry za pomocą [pliku parametrów](https://nextflow.io/docs/latest/config.html#params-file) w formacie YAML lub JSON, co sprawia, że bardzo wygodne jest zarządzanie i dystrybuowanie alternatywnych zestawów domyślnych wartości, na przykład, a także wartości parametrów specyficznych dla uruchomienia.

#### 1.3.1. Przejrzyj przykładowy plik parametrów

Aby to zademonstrować, dostarczamy przykładowy plik parametrów w bieżącym katalogu o nazwie `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Ten plik parametrów zawiera parę klucz-wartość dla każdego z wejść, które chcemy określić.
Zauważ użycie dwukropków (`:`) zamiast znaków równości (`=`), jeśli porównasz składnię z plikiem konfiguracyjnym.
Plik konfiguracyjny jest napisany w Groovy, podczas gdy plik parametrów jest napisany w YAML.

!!! info "Informacja"

    Dostarczamy również wersję JSON pliku parametrów jako przykład, ale nie będziemy jej tutaj uruchamiać.
    Możesz spróbować samodzielnie.

#### 1.3.2. Uruchom pipeline

Aby uruchomić workflow z tym plikiem parametrów, po prostu dodaj `-params-file <nazwa_pliku>` do podstawowego polecenia.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Końcowy plik wyjściowy powinien zawierać postać stegosaurus wypowiadającą pozdrowienia.

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
W takich przypadkach użycie pliku parametrów pozwoli nam podać wartości parametrów w czasie wykonania bez konieczności wpisywania masywnych poleceń wiersza poleceń i bez modyfikowania skryptu workflow.

Ułatwia to również dystrybucję zestawów parametrów do współpracowników lub jako informacji uzupełniających do publikacji, na przykład.
To sprawia, że Twoja praca jest bardziej powtarzalna przez innych.

### Podsumowanie

Wiesz już, jak wykorzystać kluczowe opcje konfiguracyjne do zarządzania wejściami workflow.

### Co dalej?

Dowiedz się, jak zarządzać tym, gdzie i jak publikowane są wyjścia workflow.

---

## 2. Zarządzanie wyjściami workflow

Do tej pory kodowaliśmy na sztywno wszystkie ścieżki dla deklaracji wyjść na poziomie workflow i, jak zauważyliśmy, gdy zaczęliśmy dodawać wiele wyjść, może to powodować pewne powtórzenia.

Przyjrzyjmy się kilku typowym sposobom konfiguracji, aby było to bardziej elastyczne.

### 2.1. Dostosuj katalog wyjściowy za pomocą `-output-dir`

Gdy kontrolujemy, jak organizowane są nasze 'publikowane' wyjścia, mamy dwa odmienne priorytety:

- Katalog wyjściowy najwyższego poziomu
- Jak pliki są organizowane wewnątrz tego katalogu

Do tej pory używaliśmy domyślnego katalogu najwyższego poziomu: `results`.
Zacznijmy od dostosowania tego, używając opcji wiersza poleceń `-output-dir`.

#### 2.1.1. Uruchom pipeline z `-output-dir`

Opcja `-output-dir` (skrót: `-o`) nadpisuje domyślny katalog wyjściowy (`results/`) dla wszystkich wyjść workflow.
To jest zalecany sposób kontrolowania ścieżki głównej, gdzie publikowane są wyjścia.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Wynik polecenia"

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

Zauważ, że nadal mamy podkatalog `hello_config` z deklaracji `path` w bloku wyjścia.
Uporządkujmy to.

#### 2.1.2. Usuń zakodowane na sztywno ścieżki z bloku wyjścia

Prefiks `hello_config/` był zakodowany na sztywno we wcześniejszych rozdziałach, ale skoro teraz uczymy się elastycznie konfigurować ścieżki wyjściowe, możemy usunąć to zakodowanie na sztywno.
Dla wyjść, które nie potrzebują podkatalogu, możemy ustawić dyrektywę `path` na pusty ciąg lub całkowicie ją usunąć.

Wprowadź następujące zmiany w kodzie w pliku workflow:

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

!!! tip "Wskazówka"

    Opcja `-output-dir` służy do kontrolowania _gdzie_ trafiają wyjścia, podczas gdy dyrektywa `path` w bloku wyjścia kontroluje _strukturę podkatalogów_.

### 2.2. Dynamiczne ścieżki wyjściowe

Oprócz zmiany katalogu wyjściowego przez wiersz poleceń, możemy również ustawić niestandardową wartość domyślną w pliku konfiguracyjnym używając `outputDir`.
To pozwala nam ustawić ścieżkę katalogu dynamicznie - nie tylko używając statycznych ciągów znaków.

#### 2.2.1. Ustaw `outputDir` w pliku konfiguracyjnym

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
    outputDir = "custom-outdir-config/${params.batch}"
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

To ustawia katalog wyjściowy na `custom-outdir-config/` plus wartość parametru `batch` jako podkatalog.
Teraz możesz zmienić lokalizację wyjściową, ustawiając parametr `--batch`:

```bash
nextflow run hello-config.nf --batch my_run
```

To publikuje wyjścia do `custom-outdir-config/my_run/`.

!!! note "Uwaga"

    Opcja wiersza poleceń `-output-dir` ma pierwszeństwo przed ustawieniem konfiguracji `outputDir`.
    Jeśli jest ustawiona, opcja konfiguracyjna będzie całkowicie zignorowana.

#### 2.2.2. Podkatalogi z nazwami partii i procesów

Możemy również ustawiać deklaracje podkatalogów `path` wyjścia dynamicznie, dla każdego wyjścia osobno.

Na przykład możemy organizować nasze wyjścia według procesu, odwołując się do `<proces>.name` w deklaracji ścieżki wyjścia:

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

Możemy pójść dalej i skonstruować bardziej złożone ścieżki podkatalogów.

W powyższej edycji usunęliśmy rozróżnienie między `intermediates` a końcowymi wyjściami na najwyższym poziomie.
Przywróćmy to, a także umieśćmy pliki w podkatalogu `params.batch`.

!!! tip "Wskazówka"

    Umieszczenie `params.batch` w bloku wyjścia `path`, zamiast w konfiguracji `outputDir`, oznacza, że nie zostanie nadpisane przez `-output-dir` w wierszu poleceń.

Najpierw zaktualizuj plik konfiguracyjny, aby usunąć `${params.batch}` z `outputDir` (ponieważ przenosimy to do deklaracji ścieżek):

=== "Po"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Następnie wprowadź następujące zmiany w pliku workflow:

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

Zobaczmy, jak to działa w praktyce, ustawiając zarówno `-output-dir` (lub `-o` w skrócie) na `custom-outdir-config-2`, jak i nazwę partii na `rep2` z wiersza poleceń:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

To publikuje wyjścia do `custom-outdir-config-2/rep2/`, z określoną ścieżką bazową _i_ podkatalogiem nazwy partii _i_ wynikami pogrupowanymi według procesu:

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

### 2.3. Ustaw tryb publikacji na poziomie workflow

Na koniec, w duchu redukcji ilości powtarzającego się kodu, możemy zastąpić deklaracje `mode` dla każdego wyjścia pojedynczą linią w konfiguracji.

#### 2.3.1. Dodaj `workflow.output.mode` do pliku konfiguracyjnego

Dodaj następujący kod do pliku `nextflow.config`:

=== "Po"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    ```

Nadanie `workflow.output.mode` wartości w pliku konfiguracyjnym wystarczy do nadpisania tego, co jest ustawione w pliku workflow, ale i tak usuńmy niepotrzebny kod.

#### 2.3.2. Usuń tryb wyjścia z pliku workflow

Wprowadź następujące zmiany w pliku workflow:

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

To bardziej zwięzłe, prawda?

#### 2.3.3. Uruchom pipeline

Przetestujmy, czy działa poprawnie:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

To publikuje wyjścia do `config-output-mode/`, i nadal są to wszystkie właściwe kopie, a nie dowiązania symboliczne.

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

Głównym powodem, dla którego możesz nadal chcieć używać sposobu ustawiania trybu dla każdego wyjścia osobno, jest sytuacja, gdy chcesz mieszać i dopasowywać w ramach tego samego workflow, tzn. mieć niektóre wyjścia kopiowane, a niektóre dowiązywane symbolicznie.

Jest wiele innych opcji, które możesz dostosować w ten sposób, ale mam nadzieję, że to daje Ci poczucie zakresu opcji i tego, jak je efektywnie wykorzystać zgodnie z Twoimi preferencjami.

### Podsumowanie

Wiesz już, jak kontrolować nazewnictwo i strukturę katalogów, w których publikowane są wyjścia, a także tryb publikacji wyjść workflow.

### Co dalej?

Dowiedz się, jak dostosować konfigurację workflow do środowiska obliczeniowego, zaczynając od technologii pakowania oprogramowania.

---

## 3. Wybierz technologię pakowania oprogramowania

Do tej pory przyglądaliśmy się elementom konfiguracji kontrolującym, jak wejścia trafiają do pipeline'u i gdzie wychodzą wyjścia. Teraz czas skupić się bardziej konkretnie na dostosowaniu konfiguracji workflow do środowiska obliczeniowego.

Pierwszym krokiem na tej ścieżce jest określenie, skąd będą pochodzić pakiety oprogramowania uruchamiane w każdym kroku.
Czy są już zainstalowane w lokalnym środowisku obliczeniowym?
Czy musimy pobrać obrazy i uruchomić je przez system kontenerowy?
A może musimy pobrać pakiety Conda i zbudować lokalne środowisko Conda?

W pierwszej części tego kursu (Części 1-4) używaliśmy po prostu lokalnie zainstalowanego oprogramowania w naszym workflow.
Następnie w Części 5 wprowadziliśmy kontenery Docker i plik `nextflow.config`, którego użyliśmy do włączenia użycia kontenerów Docker.

Teraz zobaczmy, jak możemy skonfigurować alternatywną opcję pakowania oprogramowania przez plik `nextflow.config`.

### 3.1. Wyłącz Docker i włącz Conda w pliku konfiguracyjnym

Udawajmy, że pracujemy na klastrze HPC i administrator nie zezwala na użycie Docker ze względów bezpieczeństwa.
Na szczęście dla nas Nextflow obsługuje wiele innych technologii kontenerowych, w tym Singularity (który jest szerzej używany na HPC), oraz menedżery pakietów oprogramowania takie jak Conda.

Możemy zmienić nasz plik konfiguracyjny, aby używać [Conda](https://nextflow.io/docs/latest/conda.html) zamiast Docker.
W tym celu zmieńmy wartość `docker.enabled` na `false` i dodajmy dyrektywę włączającą użycie Conda:

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

Żeby było jasne, nie _zastępujemy_ dyrektywy `docker`, _dodajemy_ alternatywną opcję.

!!! tip "Wskazówka"

    Jest kilka różnych sposobów na uzyskanie URI dla danego pakietu Conda.
    Zalecamy użycie wyszukiwarki [Seqera Containers](https://seqera.io/containers/), która da Ci URI, które możesz skopiować i wkleić, nawet jeśli nie planujesz tworzyć z niego kontenera.

### 3.3. Uruchom workflow, aby sprawdzić, czy może używać Conda

Spróbujmy.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Wynik polecenia"

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

Za kulisami Nextflow pobrał pakiety Conda i utworzył środowisko, co normalnie wymaga trochę pracy; więc miło jest, że nie musimy niczego robić sami!

!!! note "Uwaga"

    To działa szybko, ponieważ pakiet `cowpy` jest dość mały, ale jeśli pracujesz z dużymi pakietami, może to zająć trochę więcej czasu za pierwszym razem i możesz zobaczyć, że wyjście konsoli pozostaje 'zawieszone' przez minutę lub dłużej przed zakończeniem.
    To normalne i wynika z dodatkowej pracy, którą Nextflow wykonuje przy pierwszym użyciu nowego pakietu.

Z naszej perspektywy wygląda to tak, jakby działało dokładnie tak samo jak uruchamianie z Docker, mimo że w tle mechanika jest nieco inna.

Oznacza to, że jesteśmy gotowi do uruchamiania ze środowiskami Conda w razie potrzeby.

??? info "Mieszanie Docker i Conda"

    Ponieważ te dyrektywy są przypisywane do każdego procesu osobno, możliwe jest 'mieszanie i dopasowywanie', tzn. konfigurowanie niektórych procesów w workflow do uruchamiania z Docker, a innych z Conda, na przykład, jeśli używana infrastruktura obliczeniowa obsługuje oba.
    W takim przypadku włączyłbyś zarówno Docker, jak i Conda w pliku konfiguracyjnym.
    Jeśli oba są dostępne dla danego procesu, Nextflow priorytetyzuje kontenery.

    Jak wspomniano wcześniej, Nextflow obsługuje wiele innych technologii pakowania oprogramowania i kontenerowych, więc nie jesteś ograniczony tylko do tych dwóch.

### Podsumowanie

Wiesz już, jak skonfigurować, jakiego pakietu oprogramowania każdy proces powinien używać, i jak przełączać się między technologiami.

### Co dalej?

Dowiedz się, jak zmienić platformę wykonawczą używaną przez Nextflow do faktycznego wykonywania pracy.

---

## 4. Wybierz platformę wykonawczą

Do tej pory uruchamialiśmy nasz pipeline z lokalnym executorem.
Wykonuje on każde zadanie na komputerze, na którym działa Nextflow.
Gdy Nextflow startuje, sprawdza dostępne procesory i pamięć.
Jeśli zasoby zadań gotowych do uruchomienia przekraczają dostępne zasoby, Nextflow wstrzyma ostatnie zadania przed wykonaniem, dopóki jedno lub więcej wcześniejszych nie zakończy się, zwalniając niezbędne zasoby.

Lokalny executor jest wygodny i wydajny, ale pozostaje ograniczony do pojedynczego komputera. Przy bardzo dużych obciążeniach możesz odkryć, że Twoja lokalna maszyna stanowi wąskie gardło, albo ze względu na pojedyncze zadanie wymagające więcej zasobów niż masz dostępne, albo z powodu tak wielu zadań, że czekanie na pojedynczy komputer zajęłoby zbyt długo.

Nextflow obsługuje [wiele różnych executorów](https://nextflow.io/docs/latest/executor.html), w tym harmonogramy HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor i inne), a także backendy wykonawcze w chmurze (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i więcej).

### 4.1. Kierowanie na inny backend

Wybór executora jest ustawiany przez dyrektywę procesu o nazwie `executor`.
Domyślnie jest ustawiony na `local`, więc następująca konfiguracja jest domyślna:

```groovy title="Wbudowana konfiguracja"
process {
    executor = 'local'
}
```

Aby ustawić executor na inny backend, wystarczy określić żądany executor używając podobnej składni jak opisano powyżej dla alokacji zasobów (zobacz [dokumentację executorów](https://nextflow.io/docs/latest/executor.html) dla wszystkich opcji).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Ostrzeżenie"

    Nie możemy tego przetestować w środowisku szkoleniowym, ponieważ nie jest ono skonfigurowane do łączenia się z HPC.

### 4.2. Radzenie sobie ze składnią specyficzną dla backendu dla parametrów wykonania

Większość platform obliczeniowych wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry, takie jak żądania alokacji zasobów i ograniczenia (np. liczba procesorów i pamięć) oraz nazwę kolejki zadań do użycia.

Niestety, każdy z tych systemów używa różnych technologii, składni i konfiguracji do definiowania, jak zadanie powinno być zdefiniowane i przesłane do odpowiedniego harmonogramu.

??? abstract "Przykłady"

    Na przykład, to samo zadanie wymagające 8 procesorów i 4GB RAM do wykonania w kolejce "my-science-work" musi być wyrażone w różny sposób w zależności od backendu.

    ```bash title="Konfiguracja dla SLURM / wysyłanie przez sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Konfiguracja dla PBS / wysyłanie przez qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Konfiguracja dla SGE / wysyłanie przez qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Na szczęście Nextflow upraszcza to wszystko.
Dostarcza ustandaryzowaną składnię, dzięki której możesz określić odpowiednie właściwości takie jak [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) i [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (zobacz [dyrektywy procesów](https://nextflow.io/docs/latest/reference/process.html#process-directives) dla innych właściwości) tylko raz.
Następnie, w czasie wykonania, Nextflow użyje tych ustawień do wygenerowania odpowiednich skryptów specyficznych dla backendu na podstawie ustawienia executora.

Omówimy tę ustandaryzowaną składnię w następnej sekcji.

### Podsumowanie

Wiesz już, jak zmienić executor, aby używać różnych rodzajów infrastruktury obliczeniowej.

### Co dalej?

Dowiedz się, jak oceniać i wyrażać alokacje zasobów i ograniczenia w Nextflow.

---

## 5. Kontroluj alokacje zasobów obliczeniowych

Większość platform obliczeniowych wysokiej wydajności pozwala (a czasami wymaga), abyś określił pewne parametry alokacji zasobów, takie jak liczba procesorów i pamięć.

Domyślnie Nextflow użyje jednego procesora i 2GB pamięci dla każdego procesu.
Odpowiednie dyrektywy procesu nazywają się `cpus` i `memory`, więc następująca konfiguracja jest domyślna:

```groovy title="Wbudowana konfiguracja" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Możesz modyfikować te wartości, dla wszystkich procesów lub dla konkretnych nazwanych procesów, używając dodatkowych dyrektyw procesu w pliku konfiguracyjnym.
Nextflow przetłumaczy je na odpowiednie instrukcje dla wybranego executora.

Ale skąd wiesz, jakich wartości użyć?

### 5.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów

Jeśli nie wiesz z góry, ile procesora i pamięci Twoje procesy prawdopodobnie będą potrzebować, możesz przeprowadzić profilowanie zasobów, co oznacza, że uruchamiasz workflow z domyślnymi alokacjami, rejestrujesz, ile każdy proces użył, i na tej podstawie szacujesz, jak dostosować bazowe alokacje.

Wygodnie jest to, że Nextflow zawiera wbudowane narzędzia do tego i chętnie wygeneruje dla Ciebie raport na żądanie.

Aby to zrobić, dodaj `-with-report <nazwa_pliku>.html` do wiersza poleceń.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Raport to plik html, który możesz pobrać i otworzyć w przeglądarce. Możesz również kliknąć go prawym przyciskiem myszy w eksploratorze plików po lewej stronie i kliknąć `Show preview`, aby wyświetlić go w środowisku szkoleniowym.

Poświęć kilka minut na przejrzenie raportu i sprawdź, czy możesz zidentyfikować możliwości dostosowania zasobów.
Upewnij się, że klikasz na zakładki pokazujące wyniki wykorzystania jako procent tego, co zostało przydzielone.

Zobacz [Raporty](https://nextflow.io/docs/latest/reports.html) dla dokumentacji opisującej wszystkie dostępne funkcje.

### 5.2. Ustaw alokacje zasobów dla wszystkich procesów

Profilowanie pokazuje, że procesy w naszym szkoleniowym workflow są bardzo lekkie, więc zmniejszmy domyślną alokację pamięci do 1GB na proces.

Dodaj następujący kod do pliku `nextflow.config`, przed sekcją parametrów pipeline'u:

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }

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
    docker.enabled = false
    conda.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

To pomoże zmniejszyć ilość zużywanej mocy obliczeniowej.

### 5.3. Ustaw alokacje zasobów dla konkretnego procesu

Jednocześnie udawajmy, że proces `cowpy` wymaga więcej zasobów niż inne, tylko po to, abyśmy mogli zademonstrować, jak dostosować alokacje dla indywidualnego procesu.

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

Z tą konfiguracją wszystkie procesy będą żądać 1GB pamięci i jednego procesora (domyślna wartość domniemana), z wyjątkiem procesu `cowpy`, który będzie żądał 2GB i 2 procesorów.

!!! tip "Wskazówka"

    Jeśli masz maszynę z małą liczbą procesorów i przydzielasz dużą liczbę na proces, możesz zobaczyć, że wywołania procesów są kolejkowane jedno za drugim.
    To dlatego, że Nextflow zapewnia, że nie żądamy więcej procesorów niż jest dostępnych.

### 5.4. Uruchom workflow ze zaktualizowaną konfiguracją

Spróbujmy, dostarczając inną nazwę pliku dla raportu profilowania, abyśmy mogli porównać wydajność przed i po zmianach konfiguracji.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Prawdopodobnie nie zauważysz żadnej rzeczywistej różnicy, ponieważ to takie małe obciążenie, ale to jest podejście, którego użyłbyś do analizy wydajności i wymagań zasobowych rzeczywistego workflow.

Jest to bardzo przydatne, gdy Twoje procesy mają różne wymagania zasobowe. Daje Ci to możliwość odpowiedniego doboru alokacji zasobów dla każdego procesu na podstawie rzeczywistych danych, a nie zgadywania.

!!! tip "Wskazówka"

    To tylko mały przedsmak tego, co możesz zrobić, aby zoptymalizować wykorzystanie zasobów.
    Sam Nextflow ma wbudowaną naprawdę fajną [dynamiczną logikę ponawiania](https://nextflow.io/docs/latest/process.html#dynamic-task-resources), która ponawia zadania, które nie powiodły się z powodu ograniczeń zasobowych.
    Dodatkowo platforma Seqera oferuje narzędzia oparte na AI do automatycznej optymalizacji alokacji zasobów.

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
Jednak gdybyś spróbował uruchomić workflow z alokacjami zasobów przekraczającymi te limity, a następnie sprawdził polecenie `sbatch` w pliku skryptu `.command.run`, zobaczyłbyś, że żądania faktycznie wysyłane do executora są ograniczone do wartości określonych przez `resourceLimits`.

??? info "Konfiguracje referencyjne instytucji"

    Projekt nf-core skompilował [zbiór plików konfiguracyjnych](https://nf-co.re/configs/) udostępnionych przez różne instytucje na całym świecie, obejmujących szeroki zakres executorów HPC i chmurowych.

    Te współdzielone konfiguracje są wartościowe zarówno dla osób, które tam pracują i mogą zatem po prostu wykorzystać konfigurację Swojej instytucji od razu, jak i jako model dla osób, które chcą opracować konfigurację dla własnej infrastruktury.

### Podsumowanie

Wiesz już, jak generować raport profilowania do oceny wykorzystania zasobów i jak modyfikować alokacje zasobów dla wszystkich procesów i/lub dla indywidualnych procesów, a także ustawiać ograniczenia zasobów dla uruchamiania na HPC.

### Co dalej?

Dowiedz się, jak skonfigurować predefiniowane profile konfiguracji i przełączać się między nimi w czasie wykonania.

---

## 6. Używaj profili do przełączania między predefiniowanymi konfiguracjami

Pokazaliśmy Ci wiele sposobów dostosowywania konfiguracji pipeline'u w zależności od projektu, nad którym pracujesz, lub środowiska, w którym wykonujesz obliczenia.

Możesz chcieć przełączać się między alternatywnymi ustawieniami w zależności od tego, jakiej infrastruktury obliczeniowej używasz. Na przykład możesz chcieć rozwijać i testować małe próbki lokalnie na laptopie, a następnie uruchamiać pełnoskalowe obciążenia na HPC lub w chmurze.

Nextflow pozwala skonfigurować dowolną liczbę [profili](https://nextflow.io/docs/latest/config.html#config-profiles) opisujących różne konfiguracje, które możesz następnie wybrać w czasie wykonania używając argumentu wiersza poleceń, zamiast modyfikować sam plik konfiguracyjny.

### 6.1. Utwórz profile do przełączania między lokalnym rozwojem a wykonaniem na HPC

Skonfigurujmy dwa alternatywne profile; jeden do uruchamiania małych obciążeń na zwykłym komputerze, gdzie będziemy używać kontenerów Docker, i jeden do uruchamiania na uniwersyteckim HPC z harmonogramem Slurm, gdzie będziemy używać pakietów Conda.

#### 6.1.1. Skonfiguruj profile

Dodaj następujący kod do pliku `nextflow.config`, po sekcji parametrów pipeline'u, ale przed ustawieniami wyjść:

=== "Po"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

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

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="15"
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
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Widzisz, że dla uniwersyteckiego HPC określamy również ograniczenia zasobów.

#### 6.1.2. Uruchom workflow z profilem

Aby określić profil w naszym poleceniu Nextflow, używamy argumentu `-profile`.

Spróbujmy uruchomić workflow z konfiguracją `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Jak widać, pozwala nam to bardzo wygodnie przełączać się między konfiguracjami w czasie wykonania.

!!! warning "Ostrzeżenie"

    Profil `univ_hpc` nie będzie działał poprawnie w środowisku szkoleniowym, ponieważ nie mamy dostępu do harmonogramu Slurm.

Jeśli w przyszłości znajdziemy inne elementy konfiguracji, które zawsze współwystępują z tymi, możemy po prostu dodać je do odpowiednich profili.
Możemy również tworzyć dodatkowe profile, jeśli są inne elementy konfiguracji, które chcemy pogrupować razem.

### 6.2. Utwórz profil parametrów testowych

Profile nie są tylko do konfiguracji infrastruktury.
Możemy ich również używać do ustawiania domyślnych wartości dla parametrów workflow, aby ułatwić innym wypróbowanie workflow bez konieczności samodzielnego zbierania odpowiednich wartości wejściowych.
Możesz to uznać za alternatywę dla używania pliku parametrów.

#### 6.2.1. Skonfiguruj profil

Składnia wyrażania domyślnych wartości w tym kontekście wygląda tak, dla profilu, który nazwiemy `test`:

```groovy title="Przykład składni"
    test {
        params.<parametr1>
        params.<parametr2>
        ...
    }
```

Jeśli dodamy profil testowy dla naszego workflow, blok `profiles` staje się:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
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

Podobnie jak w przypadku profili konfiguracji technicznej, możesz skonfigurować wiele różnych profili określających parametry pod dowolną nazwą.

#### 6.2.2. Uruchom workflow lokalnie z profilem testowym

Wygodnie jest to, że profile nie wykluczają się wzajemnie, więc możemy określić wiele profili w naszym poleceniu używając następującej składni `-profile <profil1>,<profil2>` (dla dowolnej liczby profili).

Jeśli łączysz profile, które ustawiają wartości dla tych samych elementów konfiguracji i są opisane w tym samym pliku konfiguracyjnym, Nextflow rozwiąże konflikt, używając wartości, którą wczytał jako ostatnią (tzn. cokolwiek pojawia się później w pliku).
Jeśli konfliktowe ustawienia są ustawione w różnych źródłach konfiguracji, obowiązuje domyślna [kolejność pierwszeństwa](https://nextflow.io/docs/latest/config.html).

Spróbujmy dodać profil testowy do naszego poprzedniego polecenia:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

To użyje Docker tam, gdzie to możliwe, i wyprodukuje wyjścia w `custom-outdir-config/test`, a tym razem postać to komediowy duet `dragonandcow`.

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

To oznacza, że dopóki dystrybuujemy pliki danych testowych wraz z kodem workflow, każdy może szybko wypróbować workflow bez konieczności dostarczania własnych wejść przez wiersz poleceń lub plik parametrów.

!!! tip "Wskazówka"

    Możemy wskazywać na adresy URL dla większych plików przechowywanych zewnętrznie.
    Nextflow pobierze je automatycznie, o ile jest otwarte połączenie.

    Więcej szczegółów znajdziesz w Side Quest [Praca z plikami](../side_quests/working_with_files.md)

### 6.3. Użyj `nextflow config`, aby zobaczyć rozwiązaną konfigurację

Jak wspomniano powyżej, czasami ten sam parametr może być ustawiony na różne wartości w profilach, które chcesz połączyć.
I bardziej ogólnie, jest wiele miejsc, gdzie elementy konfiguracji mogą być przechowywane, a czasami te same właściwości mogą być ustawione na różne wartości w różnych miejscach.

Nextflow stosuje ustaloną [kolejność pierwszeństwa](https://nextflow.io/docs/latest/config.html) do rozwiązywania konfliktów, ale może to być trudne do samodzielnego określenia.
A nawet jeśli nic nie jest w konflikcie, może być nużące przeglądanie wszystkich możliwych miejsc, gdzie rzeczy mogłyby być skonfigurowane.

Na szczęście Nextflow zawiera wygodne narzędzie o nazwie `config`, które może zautomatyzować cały ten proces za Ciebie.

Narzędzie `config` przeszuka całą zawartość w Twoim bieżącym katalogu roboczym, zbierze wszystkie pliki konfiguracyjne i wyprodukuje w pełni rozwiązaną konfigurację, której Nextflow użyłby do uruchomienia workflow.
Pozwala to dowiedzieć się, jakie ustawienia zostaną użyte, bez konieczności uruchamiania czegokolwiek.

#### 6.3.1. Rozwiąż domyślną konfigurację

Uruchom to polecenie, aby rozwiązać konfigurację, która byłaby zastosowana domyślnie.

```bash
nextflow config
```

??? success "Wynik polecenia"

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

To pokazuje bazową konfigurację, którą otrzymujesz, jeśli nie określisz niczego dodatkowego w wierszu poleceń.

#### 6.3.2. Rozwiąż konfigurację z aktywowanymi konkretnymi ustawieniami

Jeśli podasz parametry wiersza poleceń, np. włączając jeden lub więcej profili lub ładując plik parametrów, polecenie dodatkowo je uwzględni.

```bash
nextflow config -profile my_laptop,test
```

??? success "Wynik polecenia"

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

To staje się szczególnie przydatne dla złożonych projektów, które obejmują wiele warstw konfiguracji.

### Podsumowanie

Wiesz już, jak używać profili do wybierania predefiniowanej konfiguracji w czasie wykonania z minimalnym wysiłkiem.
Bardziej ogólnie, wiesz, jak konfigurować wykonania workflow, aby pasowały do różnych platform obliczeniowych i zwiększały powtarzalność Twoich analiz.

### Co dalej?

Świętuj i pogratuluj sobie! Ukończyłeś Swój pierwszy kurs dla deweloperów Nextflow.

Przejdź do końcowego [podsumowania kursu](./next_steps.md), aby przejrzeć, czego się nauczyłeś i dowiedzieć się, co dalej.

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
Co ma pierwszeństwo, gdy ten sam parametr jest ustawiony zarówno w pliku konfiguracyjnym, jak i w wierszu poleceń?
- [ ] Wartość z pliku konfiguracyjnego
- [x] Wartość z wiersza poleceń
- [ ] Pierwsza napotkana wartość
- [ ] Żadna; powoduje to błąd

Dowiedz się więcej: [1.1. Przenieś domyślne wartości do `nextflow.config`](#11-przenies-domyslne-wartosci-do-nextflowconfig)
</quiz>

<quiz>
Czy możesz mieć włączone zarówno Docker, jak i Conda w tej samej konfiguracji?
- [x] Tak, Nextflow może używać obu w zależności od dyrektyw procesu
- [ ] Nie, tylko jedno może być włączone naraz
- [ ] Tak, ale tylko w profilach
- [ ] Nie, wykluczają się wzajemnie
</quiz>

<quiz>
Jeśli zarówno Docker, jak i Conda są włączone, a proces ma obie dyrektywy, które jest priorytetyzowane?
- [x] Docker (kontenery)
- [ ] Conda
- [ ] To, które jest zdefiniowane jako pierwsze
- [ ] Powoduje to błąd

Dowiedz się więcej: [3. Wybierz technologię pakowania oprogramowania](#3-wybierz-technologie-pakowania-oprogramowania)
</quiz>

<quiz>
Jaka jest domyślna alokacja pamięci dla procesów Nextflow?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Bez limitu
</quiz>

<quiz>
Jak ustawiasz wymagania zasobów dla konkretnego procesu w pliku konfiguracyjnym?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Dowiedz się więcej: [5.3. Ustaw alokacje zasobów dla konkretnego procesu](#53-ustaw-alokacje-zasobow-dla-konkretnego-procesu)
</quiz>

<quiz>
Jaka opcja wiersza poleceń generuje raport wykorzystania zasobów?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Dowiedz się więcej: [5.1. Uruchom workflow, aby wygenerować raport wykorzystania zasobów](#51-uruchom-workflow-aby-wygenerowac-raport-wykorzystania-zasobow)
</quiz>

<quiz>
Co robi dyrektywa `resourceLimits`?
- [ ] Ustawia minimalne wymagania zasobów
- [ ] Przydziela zasoby do procesów
- [x] Ogranicza maksymalne zasoby, które mogą być żądane
- [ ] Monitoruje wykorzystanie zasobów

Dowiedz się więcej: [5.5. Dodaj limity zasobów](#55-dodaj-limity-zasobow)
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
Jak określasz plik parametrów podczas uruchamiania Nextflow?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Dowiedz się więcej: [1.3. Użyj pliku parametrów](#13-uzyj-pliku-parametrow)
</quiz>

<quiz>
Do czego mogą być używane profile? (Wybierz wszystkie pasujące)
- [x] Definiowanie ustawień specyficznych dla infrastruktury
- [x] Ustawianie limitów zasobów dla różnych środowisk
- [x] Dostarczanie parametrów testowych
- [ ] Definiowanie nowych procesów

Dowiedz się więcej: [6. Używaj profili do przełączania między predefiniowanymi konfiguracjami](#6-uzywaj-profili-do-przelaczania-miedzy-predefiniowanymi-konfiguracjami)
</quiz>

<quiz>
Jak określasz wiele profili w pojedynczym poleceniu?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Dowiedz się więcej: [6. Używaj profili do przełączania między predefiniowanymi konfiguracjami](#6-uzywaj-profili-do-przelaczania-miedzy-predefiniowanymi-konfiguracjami)
</quiz>
