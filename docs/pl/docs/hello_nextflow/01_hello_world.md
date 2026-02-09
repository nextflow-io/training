# Część 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) na kanale YouTube Nextflow'a.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/01_hello_world.md).
///

W tej pierwszej części szkolenia Hello Nextflow zaczynamy od bardzo prostego, niezależnego od dziedziny przykładu Hello World, który będziemy stopniowo rozbudowywać, aby zademonstrować użycie podstawowej logiki i komponentów Nextflow'a.

??? info "Czym jest przykład Hello World?"

    „Hello World!" to minimalistyczny przykład, który ma na celu zademonstrowanie podstawowej składni i struktury języka programowania lub frameworka programistycznego.
    Przykład zazwyczaj polega na wypisaniu frazy „Hello, World!" do urządzenia wyjściowego, takiego jak konsola lub terminal, lub zapisaniu jej do pliku.

---

## 0. Rozgrzewka: Uruchom przykład Hello World bezpośrednio

Zademonstrujmy to prostym poleceniem, które uruchomimy bezpośrednio w terminalu, aby pokazać, co robi, zanim opakujemy je w Nextflow'a.

!!! tip

    Pamiętaj, że powinieneś teraz znajdować się w katalogu `hello-nextflow/`, jak opisano na stronie [Pierwsze kroki](00_orientation.md).

### 0.1. Spraw, aby terminal przywitał się z Tobą

Uruchom następujące polecenie w swoim terminalu.

```bash
echo 'Hello World!'
```

??? success "Wyjście polecenia"

    ```console
    Hello World!
    ```

Wypisuje to tekst 'Hello World' bezpośrednio w terminalu.

### 0.2. Zapisz wyjście do pliku

Uruchamianie pipeline'ów polega głównie na odczytywaniu danych z plików i zapisywaniu wyników do innych plików, więc zmodyfikujmy polecenie, aby zapisywało wyjście tekstowe do pliku, co uczyni przykład nieco bardziej praktycznym.

```bash
echo 'Hello World!' > output.txt
```

??? success "Wyjście polecenia"

    ```console

    ```

Nie wypisuje to niczego do terminala.

### 0.3. Znajdź wyjście

Tekst 'Hello World' powinien teraz znajdować się w pliku wyjściowym, który określiliśmy, o nazwie `output.txt`.
Możesz otworzyć go w eksploratorze plików lub z wiersza poleceń, używając na przykład narzędzia `cat`.

??? abstract "Zawartość pliku"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

To właśnie spróbujemy odtworzyć za pomocą naszego pierwszego workflow'a Nextflow'a.

### Podsumowanie

Wiesz już, jak uruchomić proste polecenie w terminalu, które wypisuje tekst, i opcjonalnie, jak sprawić, aby zapisywało wyjście do pliku.

### Co dalej?

Dowiedz się, jak wyglądałoby to samo zapisane jako workflow Nextflow'a.

---

## 1. Przeanalizuj skrypt i uruchom go

Udostępniamy Ci w pełni funkcjonalny, choć minimalistyczny skrypt workflow'a o nazwie `hello-world.nf`, który robi to samo co wcześniej (zapisuje 'Hello World!'), ale z użyciem Nextflow'a.

Aby Cię wprowadzić, otwórzmy skrypt workflow'a, abyś mógł poczuć, jak jest zbudowany.
Następnie uruchomimy go i poszukamy jego wyjść.

### 1.1. Przeanalizuj kod

Znajdziesz skrypt `hello-world.nf` w swoim bieżącym katalogu, którym powinien być `hello-nextflow`. Otwórz go w panelu edytora.

??? full-code "Pełny plik kodu"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo, aby wypisać 'Hello World!' do pliku
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // Wyemituj powitanie
        sayHello()
    }
    ```

Skrypt workflow'a Nextflow'a zazwyczaj zawiera jedną lub więcej definicji [**process**](https://nextflow.io/docs/latest/process.html) oraz sam [**workflow**](https://nextflow.io/docs/latest/workflow.html), a także kilka opcjonalnych bloków (nieobecnych tutaj), które przedstawimy później.

Każdy **process** opisuje, jakie operacje powinien wykonać odpowiedni krok w pipeline'ie, podczas gdy **workflow** opisuje logikę przepływu danych łączącą poszczególne kroki.

Najpierw przyjrzymy się bliżej blokowi **process**, a następnie blokowi **workflow**.

#### 1.1.1. Definicja `process`

Pierwszy blok kodu opisuje **process**.

Definicja procesu zaczyna się od słowa kluczowego `process`, po którym następuje nazwa procesu, a na końcu ciało procesu ograniczone nawiasami klamrowymi.
Ciało procesu musi zawierać blok skryptu, który określa polecenie do uruchomienia, które może być wszystkim, co można uruchomić w terminalu wiersza poleceń.

```groovy title="hello-world.nf" linenums="3"
/*
* Użyj echo, aby wypisać 'Hello World!' do pliku
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Tutaj mamy **process** o nazwie `sayHello`, który zapisuje swoje **wyjście** do pliku o nazwie `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

To bardzo minimalna definicja procesu, która zawiera tylko definicję `output` i `script` do wykonania.

Definicja `output` zawiera kwalifikator `path`, który mówi Nextflow'owi, że powinno to być traktowane jako ścieżka (obejmuje zarówno ścieżki katalogów, jak i pliki).
Innym powszechnym kwalifikatorem jest `val`.

Co ważne, definicja wyjścia nie _określa_, jakie wyjście zostanie utworzone.
Po prostu _deklaruje_, jakie jest oczekiwane wyjście, aby Nextflow mógł go poszukać po zakończeniu wykonania.
Jest to konieczne do weryfikacji, czy polecenie zostało wykonane pomyślnie, oraz do przekazania wyjścia do procesów dalszych, jeśli jest to potrzebne. Wyjście wygenerowane, które nie pasuje do tego, co zadeklarowano w bloku wyjścia, nie zostanie przekazane do procesów dalszych.

!!! warning

    Ten przykład jest kruchy, ponieważ zakodowaliśmy na stałe nazwę pliku wyjściowego w dwóch oddzielnych miejscach (w skrypcie i w blokach wyjścia).
    Jeśli zmienimy jedno, ale nie drugie, skrypt się zepsuje.
    Później nauczysz się sposobów używania zmiennych, aby złagodzić ten problem.

W rzeczywistym pipeline'ie proces zazwyczaj zawiera dodatkowe bloki, takie jak dyrektywy i wejścia, które przedstawimy za chwilę.

#### 1.1.2. Definicja `workflow`

Drugi blok kodu opisuje sam **workflow**.
Definicja workflow'a zaczyna się od słowa kluczowego `workflow`, po którym następuje opcjonalna nazwa, a następnie ciało workflow'a ograniczone nawiasami klamrowymi.

Tutaj mamy **workflow**, który składa się z bloku `main:` (który mówi 'to jest główne ciało workflow'a') zawierającego wywołanie procesu `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // Wyemituj powitanie
    sayHello()
}
```

To bardzo minimalna definicja **workflow'a**.
W rzeczywistym pipeline'ie workflow zazwyczaj zawiera wiele wywołań **procesów** połączonych **kanałami**, a procesy oczekują jednego lub więcej zmiennych **wejść**.

Nauczysz się, jak dodawać zmienne wejścia później w tym module szkoleniowym; a nauczysz się, jak dodawać więcej procesów i łączyć je kanałami w Części 3 tego kursu.

!!! tip

    Technicznie linia `main:` nie jest wymagana dla prostych workflow'ów takich jak ten, więc możesz napotkać workflow'y, które jej nie mają.
    Ale będziemy jej potrzebować, aby skorzystać z wyjść na poziomie workflow'a, więc równie dobrze możemy ją uwzględnić od początku.

### 1.2. Uruchom workflow'a

Patrzenie na kod nie jest tak zabawne jak jego uruchamianie, więc wypróbujmy to w praktyce.

#### 1.2.1. Uruchom workflow'a i monitoruj wykonanie

W terminalu uruchom następujące polecenie:

```bash
nextflow run hello-world.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Jeśli Twoje wyjście konsoli wygląda mniej więcej tak, to gratulacje, właśnie uruchomiłeś swój pierwszy workflow Nextflow'a!

Najważniejszym wyjściem jest tutaj ostatnia linia, która jest podświetlona w powyższym wyjściu:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Mówi nam to, że proces `sayHello` został pomyślnie wykonany raz (`1 of 1 ✔`).

Co ważne, ta linia mówi Ci również, gdzie znaleźć wyjście wywołania procesu `sayHello`.
Spójrzmy na to teraz.

#### 1.2.2. Znajdź wyjście i logi w katalogu `work`

Gdy uruchamiasz Nextflow'a po raz pierwszy w danym katalogu, tworzy on katalog o nazwie `work`, w którym zapisze wszystkie pliki (i wszelkie dowiązania symboliczne) wygenerowane w trakcie wykonania.

W katalogu `work` Nextflow organizuje wyjścia i logi dla każdego wywołania procesu.
Dla każdego wywołania procesu Nextflow tworzy zagnieżdżony podkatalog, nazwany hashem w celu uczynienia go unikalnym, gdzie przygotuje wszystkie niezbędne wejścia (domyślnie używając dowiązań symbolicznych), zapisze pliki pomocnicze oraz zapisze logi i wszelkie wyjścia procesu.

Ścieżka do tego podkatalogu jest pokazana w skróconej formie w nawiasach kwadratowych w wyjściu konsoli.
Patrząc na to, co otrzymaliśmy dla uruchomienia pokazanego powyżej, linia logu konsoli dla procesu sayHello zaczyna się od `[65/7be2fa]`. Odpowiada to następującej ścieżce katalogu: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Spójrzmy, co tam jest.

??? abstract "Zawartość katalogu"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Nie widzisz tego samego?"

    Dokładne nazwy podkatalogów będą różne w Twoim systemie.

    Jeśli przeglądasz zawartość podkatalogu zadania w eksploratorze plików VSCode, zobaczysz wszystkie pliki od razu.
    Jednak pliki logów są ustawione jako niewidoczne w terminalu, więc jeśli chcesz użyć `ls` lub `tree` do ich wyświetlenia, musisz ustawić odpowiednią opcję wyświetlania niewidocznych plików.

    ```bash
    tree -a work
    ```

Pierwszą rzeczą, na którą chcesz spojrzeć, jest rzeczywiste wyjście workflow'a, czyli plik `output.txt` wyprodukowany przez proces `sayHello`.
Otwórz go, a znajdziesz powitanie `Hello World!`, które było celem naszego minimalistycznego workflow'a.

??? abstract "Zawartość pliku"

    ```console title="output.txt"
    Hello World!
    ```

Zadziałało!

Oczywiście może się wydawać, że to dużo kodu opakowującego dla tak małego wyniku, ale wartość całego tego kodu opakowującego stanie się bardziej oczywista, gdy zaczniemy czytać pliki wejściowe i łączyć wiele kroków.

Mimo to spójrzmy również na pozostałe pliki w tym katalogu. Są to pliki pomocnicze i logów wyprodukowane przez Nextflow'a w ramach wykonania zadania.

- **`.command.begin`**: Metadane związane z początkiem wykonania wywołania procesu
- **`.command.err`**: Komunikaty błędów (`stderr`) wyemitowane przez wywołanie procesu
- **`.command.log`**: Pełne wyjście logu wyemitowane przez wywołanie procesu
- **`.command.out`**: Regularne wyjście (`stdout`) wywołania procesu
- **`.command.run`**: Pełny skrypt uruchomiony przez Nextflow'a w celu wykonania wywołania procesu
- **`.command.sh`**: Polecenie, które faktycznie zostało uruchomione przez wywołanie procesu
- **`.exitcode`**: Kod wyjścia wynikający z polecenia

Plik `.command.sh` jest szczególnie przydatny, ponieważ mówi Ci, jakie główne polecenie wykonał Nextflow, nie uwzględniając całej księgowości i konfiguracji zadania/środowiska.

??? abstract "Zawartość pliku"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

To pasuje do tego, co uruchomiliśmy wcześniej ręcznie.

W tym przypadku jest to bardzo proste, ponieważ polecenie procesu było zakodowane na stałe, ale później w kursie zobaczysz polecenia procesów, które obejmują interpolację zmiennych.
To sprawia, że szczególnie cenne jest możliwość zobaczenia dokładnie, jak Nextflow zinterpretował kod i jakie polecenie zostało wyprodukowane, gdy rozwiązujesz problemy z nieudanym uruchomieniem.

### 1.3. Uruchom workflow'a ponownie

Spróbuj ponownie uruchomić workflow'a kilka razy, a następnie spójrz na katalogi zadań w `work/`.

??? abstract "Zawartość katalogu"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Widzisz, że dla każdego uruchomienia został utworzony nowy podkatalog z pełnym zestawem plików wyjściowych i logów.
To pokazuje, że uruchomienie tego samego workflow'a kilka razy nie nadpisze wyników poprzednich uruchomień.

### Podsumowanie

Wiesz, jak rozszyfrować prosty skrypt Nextflow'a, uruchomić go i znaleźć wyjście oraz odpowiednie pliki logów w katalogu work.

### Co dalej?

Dowiedz się, jak opublikować wyjścia workflow'a w wygodniejszej lokalizacji.

---

## 2. Publikuj wyjścia

Jak właśnie się dowiedziałeś, wyjście wyprodukowane przez nasz pipeline jest zakopane w katalogu roboczym kilka warstw w głąb.
Jest to zrobione celowo; Nextflow kontroluje ten katalog i nie powinniśmy z nim wchodzić w interakcje.
Jednak to sprawia, że niewygodnie jest pobierać wyjścia, na których nam zależy.

Na szczęście Nextflow zapewnia sposób publikowania wyjść do wyznaczonego katalogu za pomocą [definicji wyjść workflow'a](https://nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Podstawowe użycie

Będzie to obejmować dwa nowe fragmenty kodu:

1. Blok `publish:` wewnątrz ciała `workflow`, deklarujący wyjścia procesów.
2. Blok `output` do skryptu określający opcje wyjścia, takie jak tryb i lokalizacja.

#### 2.1.1. Zadeklaruj wyjście procesu `sayHello`

Musimy dodać blok `publish:` do ciała workflow'a (ten sam rodzaj elementu kodu co blok `main:`) i wylistować wyjście procesu `sayHello()`.

W pliku skryptu workflow'a `hello-world.nf` dodaj następujące linie kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // Wyemituj powitanie
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // Wyemituj powitanie
        sayHello()
    }
    ```

Widzisz, że możemy odwołać się do wyjścia procesu po prostu robiąc `sayHello().out` i przypisać mu dowolną nazwę, `first_output`.

#### 2.1.2. Dodaj blok `output:` do skryptu

Teraz musimy tylko dodać blok `output:`, w którym zostanie określona ścieżka katalogu wyjściowego. Zauważ, że ten nowy blok znajduje się **poza** i **poniżej** bloku `workflow` w skrypcie.

W pliku skryptu workflow'a `hello-world.nf` dodaj następujące linie kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // Wyemituj powitanie
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // Wyemituj powitanie
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Możemy użyć tego do przypisania określonych ścieżek do dowolnych wyjść procesów zadeklarowanych w bloku `workflow`.
Później nauczysz się sposobów generowania wyrafinowanych struktur katalogów wyjściowych, ale na razie po prostu kodujemy na stałe minimalną ścieżkę dla uproszczenia.

#### 2.1.3. Uruchom workflow'a

Teraz uruchom zmodyfikowany skrypt workflow'a:

```bash
nextflow run hello-world.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

Wyjście terminala powinno wyglądać znajomo. Zewnętrznie nic się nie zmieniło.

Jednak sprawdź swój eksplorator plików: tym razem Nextflow utworzył nowy katalog o nazwie `results/`.

??? abstract "Zawartość katalogu"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

Wewnątrz katalogu `results` znajdujemy dowiązanie symboliczne do `output.txt` wyprodukowanego w katalogu work przez polecenie, które właśnie uruchomiliśmy.

To pozwala nam łatwo pobrać pliki wyjściowe bez konieczności przeszukiwania podkatalogu work.

### 2.2. Ustaw niestandardową lokalizację

Posiadanie domyślnej lokalizacji jest świetne, ale możesz chcieć dostosować, gdzie wyniki są zapisywane i jak są zorganizowane.

Na przykład możesz chcieć zorganizować swoje wyjścia w podkatalogi.
Najprostszym sposobem, aby to zrobić, jest przypisanie określonej ścieżki wyjściowej dla każdego wyjścia.

#### 2.2.1. Zmodyfikuj ścieżkę wyjściową

Po raz kolejny modyfikowanie zachowania publikowania dla określonego wyjścia jest naprawdę proste.
Aby ustawić niestandardową lokalizację, po prostu edytuj `path` odpowiednio:

=== "Po"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Ponieważ jest to ustawione na poziomie indywidualnego wyjścia, możesz określić różne lokalizacje i podkatalogi, aby dostosować je do swoich potrzeb.

#### 2.2.2. Uruchom workflow'a ponownie

Wypróbujmy to.

```bash
nextflow run hello-world.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Tym razem wynik zostaje zapisany w określonym podkatalogu.

??? abstract "Zawartość katalogu"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Widzisz, że wynik z poprzedniego wykonania nadal tam jest.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

Możesz użyć tylu poziomów zagnieżdżenia, ile chcesz.
Możliwe jest również użycie nazwy procesu lub innych zmiennych do nazwania katalogów używanych do organizowania wyników, i możliwe jest zmienienie domyślnej nazwy katalogu wyjściowego najwyższego poziomu (który jest kontrolowany przez flagę CLI `-o` lub zmienną konfiguracyjną `outputDir`).
Omówimy te opcje później w szkoleniu.

### 2.3. Ustaw tryb publikowania na kopiowanie

Domyślnie wyjścia są publikowane jako dowiązania symboliczne z katalogu `work`.
To oznacza, że istnieje tylko jeden plik w systemie plików.

Jest to świetne, gdy masz do czynienia z bardzo dużymi plikami, dla których nie chcesz przechowywać wielu kopii.
Jednak jeśli w którymś momencie usuniesz katalog work (wkrótce omówimy operacje czyszczenia), stracisz dostęp do pliku.
Więc musisz mieć plan zapisywania kopii wszelkich ważnych plików w bezpiecznym miejscu.

Jedną łatwą opcją jest przełączenie trybu publikowania na kopiowanie dla wyjść, na których Ci zależy.

#### 2.3.1. Dodaj dyrektywę trybu

Ten fragment jest naprawdę prosty.
Po prostu dodaj `mode 'copy'` do odpowiedniej definicji wyjścia na poziomie workflow'a:

=== "Po"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

To ustawia tryb publikowania dla tego konkretnego wyjścia.

#### 2.3.2. Uruchom workflow'a ponownie

Wypróbujmy to.

```bash
nextflow run hello-world.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Tym razem, jeśli spojrzysz na wyniki, plik jest właściwą kopią zamiast tylko dowiązania symbolicznego.

??? abstract "Zawartość katalogu"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Ponieważ to również jest ustawione na poziomie indywidualnego wyjścia, pozwala Ci ustawić tryb publikowania w sposób szczegółowy.
Będzie to szczególnie przydatne później, gdy przejdziemy do pipeline'ów wieloetapowych, gdzie możesz chcieć kopiować tylko końcowe wyjścia i pozostawić wyjścia pośrednie jako dowiązania symboliczne, na przykład.

Jak wspomniano wcześniej, istnieją inne, bardziej wyrafinowane opcje kontrolowania sposobu publikowania wyjść.
Pokażemy Ci, jak z nich korzystać w odpowiednim czasie na Twojej drodze z Nextflow'em.

### 2.4. Uwaga o dyrektywach `publishDir` na poziomie procesu

Do niedawna ustalonym sposobem publikowania wyjść było robienie tego na poziomie każdego indywidualnego procesu za pomocą dyrektywy `publishDir`.

Aby osiągnąć to, co właśnie zrobiliśmy dla wyjść procesu `sayHello`, zamiast tego dodalibyśmy następującą linię do definicji procesu:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Nadal znajdziesz ten wzorzec kodu wszędzie w starszych pipeline'ach Nextflow'a i modułach procesów, więc ważne jest, aby być tego świadomym.
Jednak nie zalecamy używania go w żadnej nowej pracy, ponieważ w końcu będzie niedozwolony w przyszłych wersjach języka Nextflow.

### Podsumowanie

Wiesz, jak opublikować wyjścia workflow'a w wygodniejszej lokalizacji.

### Co dalej?

Dowiedz się, jak dostarczyć zmienne wejście za pomocą parametru wiersza poleceń i efektywnie wykorzystywać wartości domyślne.

---

## 3. Użyj zmiennego wejścia przekazanego w wierszu poleceń

W swoim obecnym stanie nasz workflow używa powitania zakodowanego na stałe w poleceniu procesu.
Chcemy dodać trochę elastyczności, używając zmiennej wejściowej, abyśmy mogli łatwiej zmieniać powitanie w czasie wykonania.

Wymaga to od nas wprowadzenia trzech zestawów zmian w naszym skrypcie:

1. Zmienić proces, aby oczekiwał zmiennego wejścia
2. Skonfigurować parametr wiersza poleceń do przechwytywania danych wejściowych użytkownika
3. Przekazać wejście do procesu w ciele workflow'a

Wprowadźmy te zmiany po kolei.

### 3.1. Zmień proces `sayHello`, aby oczekiwał zmiennego wejścia

Musimy edytować definicję procesu, aby (1) akceptowała zmienną wejściową i (2) używała tej zmiennej w wierszu poleceń.

#### 3.1.1. Dodaj blok wejścia do definicji procesu

Najpierw dostosujmy definicję procesu, aby akceptowała wejście o nazwie `greeting`.

W bloku procesu wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

Zmienna `greeting` jest poprzedzona `val`, aby powiedzieć Nextflow'owi, że jest to wartość (a nie ścieżka).

#### 3.1.2. Edytuj polecenie procesu, aby używało zmiennej wejściowej

Teraz zamieniamy oryginalną wartość zakodowaną na stałe na wartość zmiennej wejściowej, którą oczekujemy otrzymać.

W bloku procesu wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

Symbol `$` i nawiasy klamrowe (`{ }`) mówią Nextflow'owi, że jest to nazwa zmiennej, która musi zostać zastąpiona rzeczywistą wartością wejściową (=interpolowana).

!!! tip

    Nawiasy klamrowe (`{ }`) były technicznie opcjonalne w poprzednich wersjach Nextflow'a, więc możesz zobaczyć starsze workflow'y, gdzie jest to zapisane jako `echo '$greeting' > output.txt`.

Teraz, gdy proces `sayHello()` jest gotowy do przyjęcia zmiennego wejścia, potrzebujemy sposobu na dostarczenie wartości wejściowej do wywołania procesu na poziomie workflow'a.

### 3.2. Skonfiguruj parametr wiersza poleceń do przechwytywania danych wejściowych użytkownika

Moglibyśmy po prostu zakodować wejście bezpośrednio, robiąc wywołanie procesu `sayHello('Hello World!')`.
Jednak gdy wykonujemy prawdziwą pracę z naszym workflow'em, będziemy chcieli móc kontrolować jego wejścia z wiersza poleceń, więc możemy zrobić coś takiego:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Na szczęście Nextflow ma wbudowany system parametrów workflow'a zwany [`params`](https://nextflow.io/docs/latest/config.html#params), który ułatwia deklarowanie i używanie parametrów CLI.

Ogólna składnia polega na zadeklarowaniu `params.<nazwa_parametru>`, aby powiedzieć Nextflow'owi, że ma oczekiwać parametru `--<nazwa_parametru>` w wierszu poleceń.

Tutaj chcemy utworzyć parametr o nazwie `--input`, więc musimy zadeklarować `params.input` gdzieś w workflow'ie.
W zasadzie możemy to napisać gdziekolwiek; ale ponieważ będziemy chcieli przekazać to do wywołania procesu `sayHello()`, możemy podłączyć to tam bezpośrednio, pisząc `sayHello(params.input)`.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // Wyemituj powitanie
    sayHello(params.input)
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // Wyemituj powitanie
    sayHello()
    ```

To mówi Nextflow'owi, aby uruchomił proces `sayHello` na wartości dostarczonej przez parametr `--input`.

W efekcie wykonaliśmy kroki (2) i (3) nakreślone na początku sekcji za jednym zamachem.

### 3.3. Uruchom polecenie workflow'a

Uruchommy to!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Jeśli dokonałeś wszystkich tych edycji poprawnie, powinieneś otrzymać kolejne pomyślne wykonanie.

Upewnij się, że otworzysz plik wyjściowy, aby sprawdzić, czy masz teraz nową wersję powitania.

??? abstract "Zawartość pliku"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Et voilà!

Zauważ, jak nowe wykonanie nadpisało plik wyjściowy opublikowany w katalogu `results`.
Jednak wyniki poprzednich uruchomień są nadal zachowane w katalogach zadań w `work`.

!!! tip

    Możesz łatwo odróżnić parametry na poziomie Nextflow'a od parametrów na poziomie pipeline'a.

    - Parametry, które dotyczą pipeline'a, zawsze przyjmują podwójny myślnik (`--`).
    - Parametry, które modyfikują ustawienie Nextflow'a, _np._ funkcja `-resume`, której użyliśmy wcześniej, przyjmują pojedynczy myślnik (`-`).

### 3.4. Użyj wartości domyślnych dla parametrów wiersza poleceń

Ok, to było wygodne, ale w wielu przypadkach ma sens dostarczenie wartości domyślnej dla danego parametru, aby nie trzeba było go określać dla każdego uruchomienia.

#### 3.4.1. Ustaw wartość domyślną dla parametru CLI

Nadajmy parametrowi `input` wartość domyślną, deklarując ją przed definicją workflow'a.

```groovy title="hello-world.nf" linenums="20"
/*
 * Parametry pipeline'a
 */
params {
    input: String = 'Holà mundo!'
}
```

Jak widzisz, możemy określić typ wejścia, którego oczekuje workflow (Nextflow 25.10.2 i nowsze).
Składnia to `nazwa: Typ = wartość_domyślna`.
Obsługiwane typy obejmują `String`, `Integer`, `Float`, `Boolean` i `Path`.

!!! info

    W starszych workflow'ach możesz zobaczyć, że cały blok `params` jest zapisany po prostu jako `input = 'Holà mundo!'`.

W miarę dodawania kolejnych parametrów do swojego pipeline'a powinieneś dodać je wszystkie do tego bloku, niezależnie od tego, czy musisz nadać im wartość domyślną, czy nie.
To ułatwi znalezienie wszystkich konfigurowalnych parametrów na pierwszy rzut oka.

#### 3.4.2. Uruchom workflow'a ponownie bez określania parametru

Teraz, gdy masz ustawioną wartość domyślną, możesz uruchomić workflow'a ponownie bez konieczności określania wartości w wierszu poleceń.

```bash
nextflow run hello-world.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

Wyjście będzie w tym samym miejscu co poprzednio, ale zawartość powinna zostać zaktualizowana nowym tekstem.

??? abstract "Zawartość pliku"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow użył wartości domyślnej parametru powitania do utworzenia wyjścia.

#### 3.4.3. Nadpisz wartość domyślną

Jeśli podasz parametr w wierszu poleceń, wartość CLI nadpisze wartość domyślną.

Wypróbuj to:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Po raz kolejny powinieneś znaleźć odpowiednie zaktualizowane wyjście w swoim katalogu wyników.

??? abstract "Zawartość pliku"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note

    W Nextflow'ie istnieje wiele miejsc, w których możesz określić wartości parametrów.
    Jeśli ten sam parametr jest ustawiony na różne wartości w wielu miejscach, Nextflow określi, jakiej wartości użyć, na podstawie kolejności pierwszeństwa, która jest opisana [tutaj](https://www.nextflow.io/docs/latest/config.html).

    Omówimy to bardziej szczegółowo w Części 6 (Konfiguracja).

### Podsumowanie

Wiesz, jak używać prostego zmiennego wejścia dostarczonego w czasie wykonania za pomocą parametru wiersza poleceń, a także jak konfigurować, używać i nadpisywać wartości domyślne.

### Co dalej?

Dowiedz się, jak wygodniej zarządzać wykonaniami.

---

## 4. Zarządzaj wykonaniami workflow'a

Wiedza o tym, jak uruchamiać workflow'y i pobierać wyjścia jest świetna, ale szybko odkryjesz, że jest kilka innych aspektów zarządzania workflow'em, które ułatwią Ci życie, szczególnie jeśli rozwijasz własne workflow'y.

Tutaj pokazujemy Ci, jak używać funkcji [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) na wypadek, gdybyś musiał ponownie uruchomić ten sam workflow, jak sprawdzać log poprzednich wykonań za pomocą [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) i jak usuwać starsze katalogi robocze za pomocą [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean).

### 4.1. Ponownie uruchom workflow'a z `-resume`

Czasami będziesz chciał ponownie uruchomić pipeline, który już wcześniej uruchomiłeś, bez powtarzania kroków, które już zostały pomyślnie ukończone.

Nextflow ma opcję zwaną [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html), która pozwala Ci to zrobić.
Konkretnie, w tym trybie wszelkie procesy, które zostały już uruchomione z dokładnie tym samym kodem, ustawieniami i wejściami, zostaną pominięte.
To oznacza, że Nextflow uruchomi tylko procesy, które dodałeś lub zmodyfikowałeś od ostatniego uruchomienia, lub którym dostarczasz nowe ustawienia lub wejścia.

Istnieją dwie kluczowe zalety tego podejścia:

- Jeśli jesteś w trakcie rozwijania swojego pipeline'a, możesz iterować szybciej, ponieważ musisz uruchomić tylko proces(y), nad którymi aktywnie pracujesz, aby przetestować swoje zmiany.
- Jeśli uruchamiasz pipeline w produkcji i coś pójdzie nie tak, w wielu przypadkach możesz naprawić problem i ponownie uruchomić pipeline, a on wznowi działanie od punktu awarii, co może zaoszczędzić Ci dużo czasu i mocy obliczeniowej.

Aby z tego skorzystać, po prostu dodaj `-resume` do swojego polecenia i uruchom je:

```bash
nextflow run hello-world.nf -resume
```

??? success "Wyjście polecenia"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

Wyjście konsoli powinno wyglądać znajomo, ale jest jedna rzecz, która jest trochę inna niż wcześniej.

Poszukaj fragmentu `cached:`, który został dodany w linii statusu procesu (linia 5), co oznacza, że Nextflow rozpoznał, że już wykonał tę pracę i po prostu ponownie użył wyniku z poprzedniego pomyślnego uruchomienia.

Możesz również zobaczyć, że hash podkatalogu roboczego jest taki sam jak w poprzednim uruchomieniu.
Nextflow dosłownie wskazuje Ci poprzednie wykonanie i mówi „Już to zrobiłem tam".

!!! tip

    Gdy ponownie uruchamiasz pipeline z `resume`, Nextflow nie nadpisuje żadnych plików opublikowanych poza katalogiem roboczym przez żadne wykonania, które zostały uruchomione pomyślnie wcześniej.

### 4.2. Sprawdź log poprzednich wykonań

Niezależnie od tego, czy rozwijasz nowy pipeline, czy uruchamiasz pipeline'y w produkcji, w pewnym momencie prawdopodobnie będziesz musiał wyszukać informacje o poprzednich uruchomieniach.
Oto jak to zrobić.

Za każdym razem, gdy uruchamiasz workflow Nextflow'a, linia zostaje zapisana do pliku logu o nazwie `history`, w ukrytym katalogu o nazwie `.nextflow` w bieżącym katalogu roboczym.

??? abstract "Zawartość pliku"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Ten plik podaje Ci znacznik czasu, nazwę uruchomienia, status, ID rewizji, ID sesji i pełny wiersz poleceń dla każdego uruchomienia Nextflow'a, które zostało uruchomione z bieżącego katalogu roboczego.

Wygodniejszym sposobem dostępu do tych informacji jest użycie polecenia `nextflow log`.

```bash
nextflow log
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

To wyprowadzi zawartość pliku logu do terminala, wzbogaconą o linię nagłówka.

Zauważysz, że ID sesji zmienia się za każdym razem, gdy uruchamiasz nowe polecenie `nextflow run`, Z WYJĄTKIEM sytuacji, gdy używasz opcji `-resume`.
W takim przypadku ID sesji pozostaje takie samo.

Nextflow używa ID sesji do grupowania informacji o cache'owaniu uruchomień w katalogu `cache`, również znajdującym się w `.nextflow`.

### 4.3. Usuń starsze katalogi robocze

Podczas procesu rozwoju zazwyczaj uruchomisz swój szkicowy pipeline wiele razy, co może prowadzić do gromadzenia wielu plików w wielu podkatalogach.

Na szczęście Nextflow zawiera pomocne podpolecenie `clean`, które może automatycznie usuwać podkatalogi robocze dla poprzednich uruchomień, na których Ci już nie zależy.

#### 4.3.1. Określ kryteria usuwania

Istnieje wiele [opcji](https://www.nextflow.io/docs/latest/reference/cli.html#clean) określania, co usunąć.

Tutaj pokazujemy Ci przykład, który usuwa wszystkie podkatalogi z uruchomień przed danym uruchomieniem, określonym za pomocą jego nazwy uruchomienia.

Wyszukaj najnowsze pomyślne uruchomienie, w którym nie użyłeś `-resume`; w naszym przypadku nazwa uruchomienia to `golden_cantor`.

Nazwa uruchomienia to wygenerowany maszynowo dwuczęściowy ciąg pokazany w nawiasach kwadratowych w linii wyjścia konsoli `Launching (...)`.
Możesz również użyć logu Nextflow'a, aby wyszukać uruchomienie na podstawie jego znacznika czasu i/lub wiersza poleceń.

#### 4.3.2. Wykonaj próbne uruchomienie

Najpierw używamy flagi próbnego uruchomienia `-n`, aby sprawdzić, co zostanie usunięte przy danym poleceniu:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Wyjście polecenia"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Twoje wyjście będzie miało różne nazwy katalogów zadań i może mieć inną liczbę linii, ale powinno wyglądać podobnie do przykładu.

Jeśli nie widzisz żadnych linii wyjściowych, albo nie podałeś prawidłowej nazwy uruchomienia, albo nie ma poprzednich uruchomień do usunięcia. Upewnij się, że zmieniłeś `golden_cantor` w przykładowym poleceniu na odpowiednią najnowszą nazwę uruchomienia w swoim logu.

#### 4.3.3. Przystąp do usuwania

Jeśli wyjście wygląda zgodnie z oczekiwaniami i chcesz przystąpić do usuwania, uruchom ponownie polecenie z flagą `-f` zamiast `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Wyjście polecenia"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Wyjście powinno być podobne do poprzedniego, ale teraz mówiące 'Removed' zamiast 'Would remove'.
Zauważ, że to nie usuwa dwuznakowych podkatalogów (jak `a3/` powyżej), ale opróżnia ich zawartość.

!!! Warning

    Usuwanie podkatalogów roboczych z poprzednich uruchomień usuwa je z cache'a Nextflow'a i usuwa wszelkie wyjścia, które były przechowywane w tych katalogach.
    To oznacza, że łamie zdolność Nextflow'a do wznowienia wykonania bez ponownego uruchamiania odpowiednich procesów.

    Jesteś odpowiedzialny za zapisanie wszelkich wyjść, na których Ci zależy lub na których planujesz polegać! To jest główny powód, dla którego wolimy używać trybu `copy` zamiast trybu `symlink` dla dyrektywy `publish`.

### Podsumowanie

Wiesz, jak publikować wyjścia do określonego katalogu, ponownie uruchamiać pipeline bez powtarzania kroków, które zostały już uruchomione w identyczny sposób, i używać polecenia `nextflow clean` do czyszczenia starych katalogów roboczych.

Bardziej ogólnie, wiesz, jak interpretować prosty workflow Nextflow'a, zarządzać jego wykonaniem i pobierać wyjścia.

### Co dalej?

Zrób sobie małą przerwę, zasłużyłeś na to!

Gdy będziesz gotowy, przejdź do [**Części 2: Hello Channels**](./02_hello_channels.md), aby dowiedzieć się, jak używać kanałów do wprowadzania wejść do swojego workflow'a, co pozwoli Ci skorzystać z wbudowanego paralelizmu przepływu danych Nextflow'a i innych potężnych funkcji.

---

## Quiz

<quiz>
Jakie są minimalne wymagane komponenty procesu Nextflow'a?
- [ ] Tylko bloki wejścia i wyjścia
- [x] Bloki wyjścia i skryptu
- [ ] Bloki wejścia, wyjścia i skryptu
- [ ] Tylko blok skryptu

Dowiedz się więcej: [1.1.1. Definicja procesu](#111-definicja-process)
</quiz>

<quiz>
Jaki jest cel bloku wyjścia w procesie?
- [ ] Wypisywanie wyników do konsoli
- [ ] Zapisywanie plików do katalogu roboczego
- [x] Deklarowanie oczekiwanych wyjść z procesu
- [ ] Definiowanie zmiennych środowiskowych

Dowiedz się więcej: [1.1.1. Definicja procesu](#111-definicja-process)
</quiz>

<quiz>
Jakie polecenie jest używane do uruchomienia workflow'a Nextflow'a?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Patrząc na katalog roboczy zadania, który plik zawiera faktyczne polecenie, które zostało wykonane?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Dowiedz się więcej: [1.2.2. Znajdź wyjście i logi w katalogu `work`](#122-znajdź-wyjście-i-logi-w-katalogu-work)
</quiz>

<quiz>
Co robi flaga `-resume`?
- [ ] Restartuje workflow'a od początku
- [ ] Wstrzymuje workflow'a
- [x] Pomija procesy, które zostały już pomyślnie ukończone
- [ ] Tworzy kopię zapasową workflow'a

Dowiedz się więcej: [4.1. Ponownie uruchom workflow'a z `-resume`](#41-ponownie-uruchom-workflowa-z--resume)
</quiz>

<quiz>
Jaki jest domyślny tryb publikowania wyjść workflow'a?
- [ ] Kopiowanie plików do katalogu wyjściowego
- [x] Tworzenie dowiązań symbolicznych w katalogu wyjściowym
- [ ] Przenoszenie plików do katalogu wyjściowego
- [ ] Kompresowanie plików w katalogu wyjściowym

Dowiedz się więcej: [2.3. Ustaw tryb publikowania na kopiowanie](#23-ustaw-tryb-publikowania-na-kopiowanie)
</quiz>

<quiz>
Jak przekazać wartość parametru do workflow'a Nextflow'a z wiersza poleceń?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Dowiedz się więcej: [3.2. Skonfiguruj parametr wiersza poleceń do przechwytywania danych wejściowych użytkownika](#32-skonfiguruj-parametr-wiersza-poleceń-do-przechwytywania-danych-wejściowych-użytkownika)
</quiz>

<quiz>
Jak odwołać się do zmiennej wewnątrz bloku skryptu Nextflow'a?
- [ ] Użyj składni `%variable%`
- [x] Użyj składni `#!groovy ${variable}`
- [ ] Użyj składni `{{variable}}`
- [ ] Użyj składni `[variable]`
</quiz>
