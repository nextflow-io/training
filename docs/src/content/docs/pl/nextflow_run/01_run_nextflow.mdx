# Część 1: Uruchamianie Nextflow'a

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej części przedstawiamy podstawowe koncepcje uruchamiania pipeline'ów Nextflow.
Zaczynamy od prostego workflow'u Hello World, a następnie przechodzimy do kompletnego, wieloetapowego pipeline'u, który przetwarza wiele danych wejściowych równolegle z użyciem kontenerów.

---

## 1. Hello World

Workflow `1-hello.nf` przyjmuje powitanie jako argument wiersza poleceń i zapisuje je do pliku.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Uruchomienie workflow'u

Uruchom następujące polecenie w terminalu.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

Kluczową linią w wynikach jest linia statusu procesu:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Informuje nas ona, że proces `sayHello` zakończył się pomyślnie i wykonał się raz.
Prefiks `[6d/740edd]` to skrócona ścieżka do katalogu roboczego zadania — więcej na ten temat poniżej.
Blok `Outputs:` wyświetlony po nim zawiera listę wszystkich plików opublikowanych przez pipeline, oznaczonych zgodnie z blokiem `output` omówionym w sekcji [1.4](#14-optional-code-walkthrough) poniżej.

### 1.2. Znajdowanie wyników

Ten workflow jest skonfigurowany tak, aby publikować wyniki w katalogu `results`.
Po uruchomieniu powinieneś/powinnaś znaleźć je w tym miejscu:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Otwórz plik i sprawdź, czy zawiera `Hello World!`.

### 1.3. Eksploracja katalogu `work/`

W tle Nextflow tworzy unikalny katalog zadania dla każdego wywołania procesu, wewnątrz katalogu o nazwie `work/`.
Hash widoczny w wynikach konsoli (`[6d/740edd]`) to ścieżka do tego katalogu.

```bash
ls work/6d/740edd*
```

W środku znajdziesz plik wynikowy wraz z kilkoma ukrytymi plikami dziennika:

- **`.command.sh`**: dokładne polecenie uruchomione przez Nextflow'a
- **`.command.out`** / **`.command.err`**: standardowe wyjście i wyjście błędów procesu
- **`.command.log`**: połączone wyjście dziennika
- **`.exitcode`**: kod wyjścia procesu

Plik `.command.sh` jest szczególnie przydatny podczas debugowania — pokazuje dokładnie to, co zostało wykonane.

### 1.4. Opcjonalnie: Omówienie kodu

Rozumienie kodu nie jest konieczne, jeśli chcesz tylko uruchamiać pipeline'y, ale jeśli jesteś ciekaw/ciekawa, warto rzucić okiem.

??? optional "Kliknij, aby zapoznać się z kodem związanym z tym ćwiczeniem"

    Otwórzmy `1-hello.nf` i przyjrzyjmy się jego głównym komponentom.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Parametry pipeline'u
     */
    params {
        input: String
    }

    workflow {

        main:
        // Wyemituj powitanie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Widzimy następujące elementy:

    - instrukcję `include` wskazującą na moduł z **procesem**
    - blok `params` definiujący parametry pipeline'u
    - blok `workflow` opisujący pracę do wykonania
    - blok `output` opisujący, co zrobić z wynikami

    Przyjrzyjmy się każdemu z nich po kolei.

    ### Moduł z `process`

    Instrukcja `include` mówi Nextflow'owi, aby załadował coś o nazwie `sayHello` z osobnego pliku z kodem.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    W tym pliku znajdziemy definicję procesu o nazwie `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

    **Proces** definiuje pojedynczy krok w pipeline'ie.
    Deklaruje swoje wejścia, wyjścia i skrypt do wykonania.
    Kwalifikator `val` oznacza, że wejście jest zwykłą wartością (string, liczba itp.).
    Kwalifikator `path` oznacza, że wyjście jest ścieżką do pliku.

    Definicję procesu można umieścić w głównym pliku workflow'u, ale przechowywanie ich w osobnych plikach modułów sprawia, że są wielokrotnego użytku: ten sam moduł może być importowany przez wiele skryptów workflow'u.

    ### Blok `params`

    Blok `params` deklaruje parametry wiersza poleceń akceptowane przez workflow:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Każdy parametr zadeklarowany tutaj staje się dostępny w wierszu poleceń z podwójnym myślnikiem (`--input`).
    Obsługiwane typy to `String`, `Integer`, `Float`, `Boolean` i `Path`.

    !!! tip "Wskazówka"

        Parametry workflow'u zawsze używają dwóch myślników (`--input`), aby odróżnić je od własnych flag CLI Nextflow'a, które używają jednego myślnika (np. `-resume`).

    ### Blok `workflow`

    Blok **workflow** definiuje logikę przepływu danych: które procesy uruchomić i w jakiej kolejności.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // Wyemituj powitanie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Tutaj wywoływany jest tylko jeden proces, więc jest to bardzo proste; bardziej realistyczne przykłady omówimy później.

    Sekcja `main:` wywołuje proces `sayHello` z wartością `--input`.
    Sekcja `publish:` zawiera listę wyników, które powinny zostać skopiowane do katalogu z wynikami.

    ### Blok `output`

    Blok `output` na dole pliku określa docelową ścieżkę i tryb kopiowania.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Każdy nazwany wpis odpowiada etykiecie `publish:` w workflow'ie i mapuje ją na podkatalog w `results/`.

### Podsumowanie

Wiesz już, jak uruchomić pipeline Nextflow i znaleźć jego wyniki, a także że praca jest wykonywana w katalogach zadań w `work/`.

### Co dalej?

Dowiedz się, jak Nextflow efektywnie obsługuje wiele danych wejściowych.

---

## 2. Przetwarzanie wielu danych wejściowych

Rzeczywiste pipeline'y zazwyczaj przetwarzają wiele porcji danych, nie tylko jedną.
Workflow `2-inputs.nf` odczytuje dane z pliku CSV i uruchamia `sayHello` raz dla każdego wiersza, równolegle.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Najpierw uruchommy workflow, a następnie przyjrzyjmy się mechanizmowi, którego Nextflow używa do obsługi wielu danych wejściowych.

### 2.1. Uruchomienie workflow'u

Uruchom następujące polecenie w terminalu.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

Wynik `3 of 3` informuje nas, że proces `sayHello` został wywołany trzy razy — raz dla każdego wiersza w pliku CSV.

W katalogu `results` powinny teraz znajdować się trzy pliki wynikowe, po jednym dla każdego powitania:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Otwórz dowolny z plików wynikowych i sprawdź, czy każdy zawiera powitanie.

Skrócony wynik powyżej pokazuje jedną linię podsumowania dla `sayHello`, ale Nextflow faktycznie uruchomił trzy osobne wykonania zadania — po jednym dla każdego wiersza w pliku CSV — i wykonał je równolegle, gdy tylko maszyna miała dostępne zasoby.

Podobnie jak w przypadku pojedynczego zadania, które eksplorowałeś/eksplorowałaś w sekcji [1.3](#13-explore-the-work-directory), każde z tych trzech wykonań otrzymuje własny katalog zadania w `work/`, całkowicie odizolowany od pozostałych:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Każdy plik `.command.sh` zawiera wyłącznie polecenie dla jednego powitania:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Ta izolacja sprawia, że równoległe wykonywanie jest bezpieczne: trzy zadania działające jednocześnie nigdy nie współdzielą katalogu roboczego, więc to, co zapisuje jedno zadanie, nie może kolidować z tym, co zapisuje inne — nawet jeśli produkują pliki o tej samej nazwie.
To również dlatego `-resume` (omówiony dalej) może buforować i ponownie wykorzystywać poszczególne zadania niezależnie: wejścia, wyjścia i dzienniki każdego zadania znajdują się całkowicie w jego własnym katalogu, bez żadnych współdzielonych elementów, które mogłyby się rozsynchornizować.

### 2.2. Ponowne uruchomienie workflow'u z `-ansi-log false`

Domyślnie Nextflow kondensuje wyniki do jednej linii podsumowania na proces.
Aby zobaczyć każde wywołanie procesu wymienione osobno, dodaj `-ansi-log false`:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Widać wszystkie trzy wywołania procesu oraz unikalny podkatalog roboczy utworzony dla każdego z nich.

### 2.3. Użycie `-resume` do pomijania ukończonej pracy

Teraz przejdź do rozszerzonego pliku wejściowego, który dodaje dwa kolejne powitania, i dodaj `-resume` do polecenia:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow uruchomił tylko dwa nowe wejścia.
Trzy powitania przetworzone w poprzednim uruchomieniu zostały zbuforowane i automatycznie ponownie wykorzystane.

Działa to również przy pomijaniu wykonywania procesów dla kroków, które zostały już pomyślnie ukończone w wieloetapowym pipeline'ie.
Na przykład gdy uruchomienie pipeline'u zostało przerwane przez błąd systemowy albo gdy dodałeś/dodałaś nowe kroki do pipeline'u w trakcie jego tworzenia.

Możliwość `-resume` jest szczególnie cenna w długich pipeline'ach, gdzie odzyskiwanie po awarii może zaoszczędzić krytyczny czas i zasoby.

### 2.4. Opcjonalnie: Omówienie kodu

Rozumienie kodu nie jest konieczne, jeśli chcesz tylko uruchamiać pipeline'y, ale jeśli jesteś ciekaw/ciekawa, warto rzucić okiem.

??? optional "Kliknij, aby zapoznać się z kodem związanym z tym ćwiczeniem"

    Kluczowa zmiana w `2-inputs.nf` znajduje się w sekcji `main:` workflow'u:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // Utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Wyemituj powitanie
        sayHello(greeting_ch)
    ```

    To, co tu widzisz, nazywa się **kanałem**: konstruktem kolejki, który obsługuje dane wejściowe w sposób ułatwiający paralelizację operacji.

    - `channel.fromPath(params.input)` tworzy kanał ze ścieżki pliku podanej przez `--input`
    - `.splitCsv()` parsuje plik CSV na wiersze
    - `#!groovy .map { line -> line[0] }` wyodrębnia pierwszą kolumnę z każdego wiersza

    Wynikiem jest kanał zawierający `Hello`, `Bonjour` i `Hola`.
    Po przekazaniu do `sayHello(greeting_ch)` Nextflow automatycznie wywołuje proces raz dla każdego elementu, uruchamiając je równolegle, gdy pozwalają na to zasoby.

### Podsumowanie

Wiesz już, jak przetwarzać wiele danych wejściowych z pliku CSV równolegle i jak używać `-resume`, aby unikać powtarzania ukończonej pracy.

### Co dalej?

Dowiedz się, jak kompletny, wieloetapowy pipeline łączy procesy ze sobą za pomocą kanałów i jak używać kontenerów do zarządzania narzędziami analitycznymi i ich zależnościami.

---

## 3. Uruchamianie wieloetapowego pipeline'u

Do tej pory uruchamiałeś/uruchamiałaś pojedynczy proces, a następnie uruchamiałeś/uruchamiałaś go wielokrotnie równolegle dla zestawu danych wejściowych.
Rzeczywiste pipeline'y zazwyczaj idą dalej: łączą kilka procesów ze sobą, przekazując wyjście jednego do następnego, i często korzystają z więcej niż jednego oprogramowania po drodze.
Workflow `main.nf` łączy oba te elementy w kompletny pipeline.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Każde powitanie przepływa przez wszystkie cztery kroki: `sayHello` zapisuje je do pliku, `convertToUpper` konwertuje tekst na wielkie litery, `collectGreetings` łączy wszystkie wyniki w jeden plik, a `cowpy` generuje ASCII art ze scalonego wyjścia przy użyciu skonteneryzowanego narzędzia.
Nextflow łączy te kroki kanałami: wyjście jednego procesu staje się wejściem następnego, więc cały łańcuch uruchamia się automatycznie w miarę dostępności danych — bez konieczności ręcznego orkiestrowania każdego kroku.

Warto zauważyć, że ten workflow używa modułów: każdy proces jest zdefiniowany we własnym pliku w katalogu `modules/`, a `main.nf` importuje je za pomocą instrukcji `include` zamiast definiować je bezpośrednio.
Dzięki temu każdy proces jest wielokrotnego użytku w wielu workflow'ach bez duplikowania kodu. Aby dowiedzieć się więcej, zapoznaj się z sekcją eksploracji kodu poniżej.

### 3.1. Uruchomienie workflow'u

Uruchom następujące polecenie w terminalu.

```bash
nextflow run main.nf --input data/greetings.csv
```

Parametr `character` domyślnie przyjmuje wartość `turkey` w `nextflow.config`, więc ASCII art przedstawia indyka, chyba że go nadpiszesz (spróbuj dodać `--character tux`).

??? success "Wyjście polecenia"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Uruchomiono cztery procesy, ale nie tę samą liczbę razy.
`sayHello` i `convertToUpper` uruchomiły się po razie dla każdego wejścia (3 of 3): każde powitanie musi zostać zapisane i zamienione na wielkie litery osobno.
`collectGreetings` i `cowpy` uruchomiły się tylko raz (1 of 1): scalanie powitań i generowanie ASCII art ma sens dopiero po zebraniu wszystkich indywidualnych wyników.
Ten kształt „rozejście-a-potem-zejście" — kilka równoległych zadań zasilających mniejszą liczbę zadań dalszych — jest powszechny w rzeczywistych pipeline'ach.

Nextflow nie czeka na zakończenie całego kroku przed rozpoczęciem następnego.
Gdy tylko jedno wyjście `sayHello` jest gotowe, odpowiadające mu zadanie `convertToUpper` może się rozpocząć, więc zadania z różnych procesów działają jednocześnie, a nie w ścisłych partiach.
`collectGreetings` i `cowpy` muszą czekać, ponieważ każdy z nich zależy od dostępności wszystkich wyników z poprzednich kroków.

Katalog `results` odzwierciedla to „zejście", a także to, co autor pipeline'u zdecydował się opublikować i gdzie: przypomnij sobie blok `output` z omówienia kodu w sekcji 1.4, który definiuje tę strukturę.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

Katalog najwyższego poziomu nosi nazwę parametru `batch`, który domyślnie przyjmuje wartość `batch`; zobaczysz, jak się zmienia w późniejszych ćwiczeniach.

Sprawdź plik `cowpy-COLLECTED-batch-output.txt`, aby zobaczyć plik z ASCII art.

??? abstract "Zawartość pliku"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
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

Podobnie jak w sekcji [2.1](#21-run-the-workflow), każde z tych 8 wykonań zadań — we wszystkich czterech procesach — otrzymuje własny katalog w `work/`, całkowicie odizolowany od pozostałych.
`collectGreetings` dobrze ilustruje, dlaczego to ma znaczenie: zależy od wyników wszystkich trzech zadań `convertToUpper`, które znajdują się w trzech różnych katalogach zadań. Dlatego Nextflow tworzy dowiązania symboliczne do tych plików wewnątrz własnego katalogu `collectGreetings`, zamiast odczytywać je bezpośrednio z katalogów zadań poprzedzających:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Każde zadanie widzi wyłącznie konkretne pliki, których potrzebuje — niezależnie od ich pochodzenia — i nigdy wewnętrznej zawartości katalogu innego zadania.
W całym pipeline'ie ta sama izolacja, którą widziałeś/widziałaś przy pojedynczym procesie w sekcji [2.1](#21-run-the-workflow), pozwala Nextflow'owi bezpiecznie uruchamiać każde zadanie z każdego procesu jednocześnie.

!!! note "Uwaga"

    Krok `cowpy` działa wewnątrz kontenera Docker, a nie opiera się na oprogramowaniu zainstalowanym lokalnie.
    Kontener pakuje aplikację wraz ze wszystkim, czego potrzebuje do działania, więc nie musisz samodzielnie instalować zależności i nimi zarządzać, a pipeline zachowuje się tak samo na każdej maszynie, która może uruchomić kontener.
    Nextflow obsługuje również Conda jako alternatywę dla kontenerów; zapoznaj się z [Częścią 2](./02_configure_pipeline.md), aby dowiedzieć się, jak między nimi przełączać.

### 3.2. Opcjonalnie: Omówienie kodu

Rozumienie kodu nie jest konieczne, jeśli chcesz tylko uruchamiać pipeline'y, ale jeśli jesteś ciekaw/ciekawa, warto rzucić okiem.

??? optional "Kliknij, aby zapoznać się z kodem związanym z tym ćwiczeniem"

    ### Jak dane przepływają z jednego kroku do następnego

    Każdy proces przekazuje swój kanał wyjściowy do następnego:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // Utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    Wzorzec `processName.out` odnosi się do kanału wyjściowego procesu.

    Operator `.collect()` zbiera wszystkie indywidualne wyjścia z `convertToUpper` w jeden element kanału przed przekazaniem ich do `collectGreetings`.

    ### Używanie modułów z procesami

    `main.nf` nie definiuje żadnego kodu procesu bezpośrednio.
    Zamiast tego importuje każdy proces z jego własnego pliku w katalogu `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Każdy plik modułu zawiera jedną definicję procesu, zbudowaną tak samo jak moduł `sayHello` w sekcji [1.4](#14-optional-code-walkthrough).
    Przechowywanie procesów w osobnych plikach sprawia, że są wielokrotnego użytku w wielu workflow'ach bez duplikowania kodu.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Używanie skonteneryzowanego oprogramowania

    Proces `cowpy` działa wewnątrz kontenera Docker określonego w jego pliku modułu:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

    Nextflow automatycznie pobiera obraz, uruchamia skrypt wewnątrz kontenera i czyści po sobie.
    Docker jest włączony dla tego projektu w `nextflow.config`:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Ta jedna linia włącza Docker dla każdego procesu w pipeline'ie, który ma określony kontener.

### Podsumowanie

Uruchomiłeś/uruchomiłaś kompletny, wieloetapowy pipeline, który przetwarza wiele danych wejściowych równolegle przy użyciu skonteneryzowanego narzędzia.

### Co dalej?

Przejdź do [Części 2](./02_configure_pipeline.md), gdzie dowiesz się, jak konfigurować zachowanie pipeline'u za pomocą `nextflow.config`.

---

## Podsumowanie

W tej części nauczyłeś/nauczyłaś się:

- Uruchamiać workflow Nextflow i znajdować jego wyniki
- Eksplorować katalog `work/` i jego pliki dziennika
- Przetwarzać wiele danych wejściowych z pliku CSV równolegle
- Używać `-resume` do pomijania ukończonej pracy przy dodawaniu nowych danych wejściowych
- Uruchamiać wieloetapowy pipeline korzystający ze skonteneryzowanego narzędzia
