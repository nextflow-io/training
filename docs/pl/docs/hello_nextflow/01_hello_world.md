# Część 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/01_hello_world.md).
///
-->

W tej pierwszej części kursu szkoleniowego Hello Nextflow łagodnie wprowadzamy temat prostym, niezależnym od dziedziny przykładem Hello World.
Będziemy go stopniowo rozbudowywać, aby zademonstrować kluczowe elementy logiki i komponentów Nextflow.

??? info "Czym jest przykład Hello World?"

    "Hello World!" to minimalistyczny przykład, który ma na celu zademonstrowanie podstawowej składni i struktury języka programowania lub frameworku oprogramowania.
    Przykład zazwyczaj polega na wypisaniu frazy "Hello, World!" na urządzeniu wyjściowym, takim jak konsola lub terminal, lub zapisaniu jej do pliku.

---

## 0. Rozgrzewka: Uruchom przykład Hello World bezpośrednio

Zademonstrujmy to prostym poleceniem, które uruchomimy bezpośrednio w terminalu, aby pokazać, co robi, zanim opakujemy je w Nextflow.

!!! tip "Wskazówka"

    Pamiętaj, że powinieneś teraz znajdować się w katalogu `hello-nextflow/`, jak opisano na stronie [Rozpoczęcie pracy](00_orientation.md).

### 0.1. Spraw, aby terminal powiedział "hello"

Uruchom następujące polecenie w terminalu.

```bash
echo 'Hello World!'
```

??? success "Wyjście polecenia"

    ```console
    Hello World!
    ```

To wyświetla tekst 'Hello World' bezpośrednio w terminalu.

### 0.2. Zapisz wyjście do pliku

Pipeline'y zazwyczaj odczytują dane z plików wejściowych i zapisują wyniki do dokumentów wyjściowych.
Zmodyfikujmy więc polecenie, aby zapisać tekst do pliku, czyniąc przykład bardziej odpowiednim.

```bash
echo 'Hello World!' > output.txt
```

??? success "Wyjście polecenia"

    ```console

    ```

To nie wyświetla niczego w terminalu.

### 0.3. Znajdź wyjście

Tekst 'Hello World' powinien teraz znajdować się w pliku wyjściowym, który określiliśmy, o nazwie `output.txt`.
Możesz go otworzyć w eksploratorze plików lub z wiersza poleceń używając na przykład narzędzia `cat`.

??? abstract "Zawartość pliku"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

To jest to, co spróbujemy odtworzyć w naszym pierwszym workflow'ie Nextflow.

### Podsumowanie

Wiesz już, jak uruchomić proste polecenie w terminalu, które wyświetla tekst, i opcjonalnie, jak sprawić, by zapisało wyjście do pliku.

### Co dalej?

Dowiedz się, jak wyglądałoby to napisane jako workflow Nextflow.

---

## 1. Przeanalizuj skrypt i uruchom go

Dostarczamy Ci w pełni funkcjonalny, choć minimalistyczny skrypt workflow'u o nazwie `hello-world.nf`, który robi to samo co wcześniej (wypisuje 'Hello World!'), ale z Nextflow.

Na początek otwórzmy skrypt workflow'u, abyś mógł zorientować się w jego strukturze.
Następnie uruchomimy go i poszukamy jego wyjść.

### 1.1. Przeanalizuj kod

Skrypt `hello-world.nf` znajdziesz w Swoim bieżącym katalogu, którym powinien być `hello-nextflow`. Otwórz go w panelu edytora.

??? full-code "Pełny plik kodu"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo do wypisania 'Hello World!' do pliku
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
        // wyemituj pozdrowienie
        sayHello()
    }
    ```

Skrypt workflow'u Nextflow zazwyczaj zawiera jedną lub więcej definicji **process** oraz sam **workflow**, plus kilka opcjonalnych bloków (nieobecnych tutaj), które wprowadzimy później.

Każdy **process** opisuje, jakie operacje powinien wykonać odpowiedni krok w pipeline'ie, podczas gdy **workflow** opisuje logikę przepływu danych, która łączy poszczególne kroki.

Najpierw przyjrzymy się bliżej blokowi **process**, a następnie blokowi **workflow**.

#### 1.1.1. Definicja `process`

Pierwszy blok kodu opisuje **process**.

Definicja procesu zaczyna się od słowa kluczowego `process`, po którym następuje nazwa i ciało ograniczone nawiasami klamrowymi.
Ciało musi zawierać blok skryptu określający polecenie do wykonania - może to być cokolwiek, co uruchomiłbyś w terminalu.

```groovy title="hello-world.nf" linenums="3"
/*
* Użyj echo do wypisania 'Hello World!' do pliku
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

To jest bardzo minimalna definicja procesu, która zawiera tylko definicję `output` i blok `script` do wykonania.

Definicja `output` zawiera kwalifikator `path`, który mówi Nextflow, że powinien to traktować jako ścieżkę (obejmuje zarówno ścieżki katalogów, jak i pliki).
Innym popularnym kwalifikatorem jest `val`.

Co ważne, definicja wyjścia nie _określa_, jakie wyjście zostanie utworzone.
Po prostu _deklaruje_, jakie jest oczekiwane wyjście, aby Nextflow mógł go szukać po zakończeniu wykonywania.
Jest to niezbędne do weryfikacji pomyślnego wykonania polecenia i przekazania wyniku do procesów downstream. Wyjście utworzone, które nie pasuje do tego, co jest zadeklarowane w bloku output, nie zostanie przekazane do procesów downstream.

!!! warning "Ostrzeżenie"

    Ten przykład jest kruchy, ponieważ zakodowaliśmy na sztywno nazwę pliku wyjściowego w dwóch oddzielnych miejscach (bloki `script` i `output`).
    Jeśli zmienimy jedno, a nie drugie, skrypt się zepsuje.
    Później nauczysz się sposobów na użycie zmiennych, aby złagodzić ten problem.

W rzeczywistym pipeline'ie proces zazwyczaj zawiera dodatkowe bloki, takie jak dyrektywy i wejścia, które wprowadzimy za chwilę.

#### 1.1.2. Definicja `workflow`

Drugi blok kodu opisuje sam **workflow**.
Definicja workflow'u zaczyna się od słowa kluczowego `workflow`, po którym następuje opcjonalna nazwa, a następnie ciało workflow'u ograniczone nawiasami klamrowymi.

Tutaj mamy **workflow**, który składa się z bloku `main:` (który mówi 'to jest główne ciało workflow'u') zawierającego wywołanie procesu `sayHello`.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // wyemituj pozdrowienie
    sayHello()
}
```

To jest bardzo minimalna definicja **workflow**.
W rzeczywistym pipeline'ie workflow zazwyczaj zawiera wiele wywołań **procesów** połączonych **kanałami**, a procesy oczekują jednego lub więcej zmiennych **wejść**.

Nauczysz się, jak dodawać zmienne wejścia później w tym module szkoleniowym; a w Części 3 tego kursu nauczysz się, jak dodawać więcej procesów i łączyć je kanałami.

!!! tip "Wskazówka"

    Technicznie linia `main:` nie jest wymagana dla prostych workflow'ów takich jak ten, więc możesz napotkać workflow'y, które jej nie mają.
    Ale będziemy jej potrzebować, aby skorzystać z wyjść na poziomie workflow'u, więc równie dobrze możemy ją uwzględnić od początku.

### 1.2. Uruchom workflow

Patrzenie na kod nie jest tak zabawne jak uruchamianie go, więc wypróbujmy to w praktyce.

#### 1.2.1. Uruchom workflow i monitoruj wykonywanie

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

Jeśli Twoje wyjście konsoli wygląda mniej więcej tak, to gratulacje, właśnie uruchomiłeś Swój pierwszy workflow Nextflow!

Najważniejszym wyjściem tutaj jest ostatnia linia, która jest podświetlona w wyjściu powyżej:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

To mówi nam, że proces `sayHello` został pomyślnie wykonany raz (`1 of 1 ✔`).

Co ważne, ta linia mówi również, gdzie znaleźć wyjście wywołania procesu `sayHello`.
Przyjrzyjmy się temu teraz.

#### 1.2.2. Znajdź wyjście i dzienniki w katalogu `work`

Kiedy uruchamiasz Nextflow po raz pierwszy w danym katalogu, tworzy on folder o nazwie `work`, w którym będzie zapisywał wszystkie pliki (i wszelkie dowiązania symboliczne) wygenerowane w trakcie wykonywania.

W katalogu `work` Nextflow organizuje wyjścia i dzienniki dla każdego wywołania procesu.
Dla każdego wywołania procesu Nextflow tworzy zagnieżdżony podkatalog, nazwany hashem, aby był unikalny, gdzie przygotuje wszystkie niezbędne wejścia (domyślnie używając dowiązań symbolicznych), zapisze pliki pomocnicze i zapisze dzienniki oraz wszelkie wyjścia procesu.

Ścieżka do tego podkatalogu jest pokazana w skróconej formie w nawiasach kwadratowych w wyjściu konsoli.
Patrząc na to, co otrzymaliśmy dla uruchomienia pokazanego powyżej, linia dziennika konsoli dla procesu sayHello zaczyna się od `[65/7be2fa]`. To odpowiada następującej ścieżce katalogu: `work/65/7be2fad5e71e5f49998f795677fd68`

Zobaczmy, co tam jest.

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
    Jednak pliki dziennika są ustawione jako niewidoczne w terminalu, więc jeśli chcesz użyć `ls` lub `tree` do ich przeglądania, będziesz musiał ustawić odpowiednią opcję wyświetlania niewidocznych plików.

    ```bash
    tree -a work
    ```

Pierwszą rzeczą, na którą chcesz spojrzeć, jest faktyczne wyjście workflow'u, czyli plik `output.txt` utworzony przez proces `sayHello`.
Otwórz go, a znajdziesz pozdrowienie `Hello World!`, które było celem naszego minimalistycznego workflow'u.

??? abstract "Zawartość pliku"

    ```console title="output.txt"
    Hello World!
    ```

Zadziałało!

Trzeba przyznać, że może to wydawać się dużo kodu opakowującego dla tak małego wyniku, ale wartość całego tego kodu opakowującego stanie się bardziej oczywista, gdy zaczniemy wczytywać pliki wejściowe i łączyć wiele kroków.

To powiedziawszy, przyjrzyjmy się również innym plikom w tym katalogu. To są pliki pomocnicze i dziennika tworzone przez Nextflow jako część wykonywania zadania.

- **`.command.begin`**: Metadane związane z początkiem wykonywania wywołania procesu
- **`.command.err`**: Komunikaty o błędach (`stderr`) emitowane przez wywołanie procesu
- **`.command.log`**: Pełne wyjście dziennika emitowane przez wywołanie procesu
- **`.command.out`**: Zwykłe wyjście (`stdout`) wywołania procesu
- **`.command.run`**: Pełny skrypt uruchomiony przez Nextflow do wykonania wywołania procesu
- **`.command.sh`**: Polecenie, które zostało faktycznie uruchomione przez wywołanie procesu
- **`.exitcode`**: Kod wyjścia wynikający z polecenia

Plik `.command.sh` jest szczególnie przydatny, ponieważ mówi Ci główne polecenie, które Nextflow wykonał, nie włączając w to całej księgowości i konfiguracji zadania/środowiska.

??? abstract "Zawartość pliku"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

To pasuje do tego, co wcześniej uruchomiliśmy ręcznie.

W tym przypadku jest to bardzo proste, ponieważ polecenie procesu było zakodowane na sztywno, ale później w kursie zobaczysz polecenia procesów, które zawierają interpolację zmiennych.
To sprawia, że szczególnie cenne jest móc zobaczyć dokładnie, jak Nextflow zinterpretował kod i jakie polecenie zostało wygenerowane, gdy rozwiązujesz problemy z nieudanym uruchomieniem.

### 1.3. Uruchom workflow ponownie

Spróbuj uruchomić workflow jeszcze kilka razy, a następnie spójrz na katalogi zadań w `work/`.

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

Widzisz, że dla każdego uruchomienia został utworzony nowy podkatalog z kompletnym zestawem plików wyjściowych i dziennika.
To pokazuje, że uruchomienie tego samego workflow'u kilka razy nie nadpisze wyników poprzednich uruchomień.

### Podsumowanie

Wiesz, jak odczytać prosty skrypt Nextflow, uruchomić go i znaleźć wyjście oraz odpowiednie pliki dziennika w katalogu work.

### Co dalej?

Naucz się publikować wyjścia workflow'u do bardziej wygodnej lokalizacji.

---

## 2. Publikuj wyjścia

Jak właśnie się dowiedziałeś, wyjście utworzone przez nasz pipeline jest zakopane w katalogu roboczym kilka poziomów w głąb.
Jest to zrobione celowo; Nextflow kontroluje ten katalog i nie powinniśmy z nim wchodzić w interakcję.
Jednak to sprawia, że jest niewygodne pobieranie wyjść, na których nam zależy.

Na szczęście Nextflow zapewnia sposób publikowania wyjść do wyznaczonego katalogu za pomocą [definicji wyjść na poziomie workflow'u](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs).

### 2.1. Podstawowe użycie

To będzie wymagało dwóch nowych fragmentów kodu:

1. Bloku `publish:` wewnątrz ciała `workflow`, deklarującego wyjścia procesu.
2. Bloku `output` w skrypcie określającego opcje wyjścia, takie jak tryb i lokalizacja.

#### 2.1.1. Zadeklaruj wyjście procesu `sayHello`

Musimy dodać blok `publish:` do ciała workflow'u (ten sam rodzaj elementu kodu co blok `main:`) i wymienić wyjście procesu `sayHello()`.

W pliku skryptu workflow'u `hello-world.nf` dodaj następujące linie kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // wyemituj pozdrowienie
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // wyemituj pozdrowienie
        sayHello()
    }
    ```

Widzisz, że możemy odwołać się do wyjścia procesu po prostu robiąc `sayHello().out` i przypisać mu dowolną nazwę, `first_output`.

#### 2.1.2. Dodaj blok `output:` do skryptu

Teraz musimy tylko dodać blok `output:`, w którym zostanie określona ścieżka katalogu wyjściowego. Zauważ, że ten nowy blok znajduje się **poza** i **poniżej** bloku `workflow` w skrypcie.

W pliku skryptu workflow'u `hello-world.nf` dodaj następujące linie kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // wyemituj pozdrowienie
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
        // wyemituj pozdrowienie
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Możemy użyć tego do przypisania określonych ścieżek do dowolnych wyjść procesów zadeklarowanych w bloku `workflow`.
Później nauczysz się sposobów generowania zaawansowanych struktur katalogów wyjściowych, ale na razie po prostu zakodowujemy minimalną ścieżkę dla prostoty.

#### 2.1.3. Uruchom workflow

Teraz uruchom zmodyfikowany skrypt workflow'u:

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

Jednak sprawdź eksplorator plików: tym razem Nextflow utworzył nowy katalog o nazwie `results/`.

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

Wewnątrz katalogu `results` znajdziemy dowiązanie symboliczne do pliku `output.txt` utworzonego w katalogu work przez polecenie, które właśnie uruchomiliśmy.

To pozwala nam łatwo pobierać pliki wyjściowe bez konieczności przekopywania się przez podkatalogi work.

### 2.2. Ustaw niestandardową lokalizację

Posiadanie domyślnej lokalizacji jest świetne, ale możesz chcieć dostosować, gdzie wyniki są zapisywane i jak są zorganizowane.

Na przykład możesz chcieć zorganizować swoje wyjścia w podkatalogi.
Najprostszym sposobem na to jest przypisanie określonej ścieżki wyjściowej dla każdego wyjścia.

#### 2.2.1. Zmodyfikuj ścieżkę wyjściową

Po raz kolejny modyfikacja zachowania publikowania dla określonego wyjścia jest naprawdę prosta.
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

Ponieważ jest to ustawiane na poziomie pojedynczego wyjścia, możesz określić różne lokalizacje i podkatalogi, aby dopasować się do Swoich potrzeb.

#### 2.2.2. Uruchom workflow ponownie

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

Tym razem wynik jest zapisywany w określonym podkatalogu.

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
Możliwe jest również użycie nazwy procesu lub innych zmiennych do nazywania katalogów używanych do organizowania wyników, a także możliwe jest zmienienie domyślnej nazwy katalogu wyjściowego najwyższego poziomu (która jest kontrolowana przez specjalną zmienną `outputDir`).
Omówimy te opcje w późniejszych szkoleniach.

### 2.3. Ustaw tryb publikowania na kopiowanie

Domyślnie wyjścia są publikowane jako dowiązania symboliczne z katalogu `work`.
To oznacza, że na systemie plików znajduje się tylko jeden plik.

To świetne, gdy masz do czynienia z bardzo dużymi plikami, dla których nie chcesz przechowywać wielu kopii.
Jednak jeśli w pewnym momencie usuniesz katalog work (wkrótce omówimy operacje czyszczenia), stracisz dostęp do pliku.
Więc musisz mieć plan zapisywania kopii wszystkich ważnych plików w bezpiecznym miejscu.

Jedną z łatwych opcji jest przełączenie trybu publikowania na kopiowanie dla wyjść, na których Ci zależy.

#### 2.3.1. Dodaj dyrektywę mode

Ta część jest naprawdę prosta.
Po prostu dodaj `mode 'copy'` do odpowiedniej definicji wyjścia na poziomie workflow'u:

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

#### 2.3.2. Uruchom workflow ponownie

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

Tym razem, jeśli spojrzysz na wyniki, plik jest właściwą kopią, a nie tylko dowiązaniem symbolicznym.

??? abstract "Zawartość katalogu"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Ponieważ to również jest ustawiane na poziomie pojedynczego wyjścia, pozwala to ustawić tryb publikowania w sposób szczegółowy.
To będzie szczególnie przydatne później, gdy przejdziemy do wieloetapowych pipeline'ów, gdzie możesz chcieć kopiować tylko końcowe wyjścia i zostawiać pośrednie wyjścia jako dowiązania symboliczne, na przykład.

Jak zauważono wcześniej, istnieją inne, bardziej zaawansowane opcje kontrolowania sposobu publikowania wyjść.
Pokażemy Ci, jak ich używać w odpowiednim czasie na Twojej drodze z Nextflow.

### 2.4. Uwaga o dyrektywach `publishDir` na poziomie procesu

Do niedawna ustalonym sposobem publikowania wyjść było robienie tego na poziomie każdego pojedynczego procesu za pomocą dyrektywy `publishDir`.

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

Nadal znajdziesz ten wzorzec kodu wszędzie w starszych pipeline'ach Nextflow i modułach procesów, więc ważne jest, aby być tego świadomym.
Jednak nie zalecamy używania go w żadnej nowej pracy, ponieważ ostatecznie zostanie niedozwolony w przyszłych wersjach języka Nextflow.

### Podsumowanie

Wiesz, jak publikować wyjścia workflow'u do bardziej wygodnej lokalizacji.

### Co dalej?

Naucz się dostarczać zmienne wejście poprzez parametr wiersza poleceń i efektywnie wykorzystywać wartości domyślne.

---

## 3. Użyj zmiennego wejścia przekazywanego z wiersza poleceń

W obecnym stanie nasz workflow używa pozdrowienia zakodowanego na sztywno w poleceniu procesu.
Chcemy dodać elastyczność za pomocą zmiennej wejściowej, co ułatwi zmianę tej wartości w czasie wykonywania.

Wymaga to wprowadzenia trzech zestawów zmian w naszym skrypcie:

1. Zmiana procesu, aby oczekiwał zmiennego wejścia
2. Skonfigurowanie parametru wiersza poleceń do przechwytywania danych wejściowych użytkownika
3. Przekazanie wejścia do procesu w ciele workflow'u

Wprowadźmy te zmiany po kolei.

### 3.1. Zmień proces `sayHello`, aby oczekiwał zmiennego wejścia

Musimy edytować definicję procesu, aby (1) akceptowała zmienną wejściową i (2) używała tej zmiennej w wierszu poleceń.

#### 3.1.1. Dodaj blok wejściowy do definicji procesu

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

Zmienna `greeting` jest poprzedzona `val`, aby powiedzieć Nextflow, że to wartość (nie ścieżka).

#### 3.1.2. Edytuj polecenie procesu, aby używać zmiennej wejściowej

Teraz zamieniamy oryginalną wartość zakodowaną na sztywno na wartość zmiennej wejściowej, którą oczekujemy otrzymać.

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

Symbol `$` i nawiasy klamrowe (`{ }`) mówią Nextflow, że to jest nazwa zmiennej, która musi być zastąpiona faktyczną wartością wejściową (=interpolowana).

!!! tip "Wskazówka"

    Nawiasy klamrowe (`{ }`) były technicznie opcjonalne w poprzednich wersjach Nextflow, więc możesz zobaczyć starsze workflow'y, gdzie to jest zapisane jako `echo '$greeting' > output.txt`.

Teraz, gdy proces `sayHello()` jest gotowy do przyjęcia zmiennego wejścia, potrzebujemy sposobu na dostarczenie wartości wejściowej do wywołania procesu na poziomie workflow'u.

### 3.2. Skonfiguruj parametr wiersza poleceń do przechwytywania danych wejściowych użytkownika

Moglibyśmy po prostu zakodować wejście na sztywno, robiąc wywołanie procesu `sayHello('Hello World!')`.
Jednak gdy wykonujemy prawdziwą pracę z naszym workflow'em, będziemy chcieli móc kontrolować jego wejścia z wiersza poleceń.

Dobra wiadomość: Nextflow ma wbudowany system parametrów workflow'u o nazwie `params`, który ułatwia deklarowanie i używanie parametrów CLI.

Ogólna składnia to zadeklarowanie `params.<nazwa_parametru>`, aby powiedzieć Nextflow, że ma oczekiwać parametru `--<nazwa_parametru>` w wierszu poleceń.

Tutaj chcemy utworzyć parametr o nazwie `--input`, więc musimy zadeklarować `params.input` gdzieś w workflow'ie.
Zasadniczo możemy to napisać gdziekolwiek; ale ponieważ zamierzamy przekazać go do wywołania procesu `sayHello()`, możemy wstawić go tam bezpośrednio, pisząc `sayHello(params.input)`.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // wyemituj pozdrowienie
    sayHello(params.input)
    ```

=== "Przed"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // wyemituj pozdrowienie
    sayHello()
    ```

To mówi Nextflow, aby uruchomił proces `sayHello` na wartości dostarczonej przez parametr `--input`.

W efekcie osiągnęliśmy kroki (2) i (3) opisane na początku sekcji za jednym razem.

### 3.3. Uruchom polecenie workflow'u

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

Jeśli wprowadzono wszystkie te edycje poprawnie, powinieneś uzyskać kolejne pomyślne wykonanie.

Koniecznie otwórz plik wyjściowy, aby sprawdzić, czy masz teraz nową wersję pozdrowienia.

??? abstract "Zawartość pliku"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà!

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Zauważ, że nowe wykonanie nadpisało plik wyjściowy opublikowany w katalogu `results`.
Jednak wyniki poprzednich uruchomień są nadal zachowane w katalogach zadań w `work`.

!!! tip "Wskazówka"

    Możesz łatwo odróżnić parametry na poziomie Nextflow od parametrów na poziomie pipeline'u.

    - Parametry, które mają zastosowanie do pipeline'u, zawsze mają podwójny myślnik (`--`).
    - Parametry, które modyfikują ustawienie Nextflow, _np._ funkcja `-resume`, której użyliśmy wcześniej, mają pojedynczy myślnik (`-`).

### 3.4. Użyj wartości domyślnych dla parametrów wiersza poleceń

Ok, to było wygodne, ale w wielu przypadkach ma sens dostarczenie wartości domyślnej dla danego parametru, abyś nie musiał go określać przy każdym uruchomieniu.

#### 3.4.1. Ustaw wartość domyślną dla parametru CLI

Nadajmy parametrowi `input` wartość domyślną, deklarując go przed definicją workflow'u.

```groovy title="hello-world.nf" linenums="20"
/*
 * Parametry pipeline'u
 */
params {
    input: String = 'Holà mundo!'
}
```

Jak widzisz, możemy określić typ wejścia, którego oczekuje workflow (Nextflow 25.10.2 i nowsze).
Składnia to `nazwa: Typ = wartość_domyślna`.
Obsługiwane typy to `String`, `Integer`, `Float`, `Boolean` i `Path`.

!!! info "Informacja"

    W starszych workflow'ach możesz zobaczyć, że cały blok `params` jest zapisany po prostu jako `input = 'Holà mundo!'`.

Gdy dodajesz więcej parametrów do Swojego pipeline'u, powinieneś dodawać je wszystkie do tego bloku, niezależnie od tego, czy musisz nadać im wartość domyślną.
To ułatwi znalezienie wszystkich konfigurowalnych parametrów na pierwszy rzut oka.

#### 3.4.2. Uruchom workflow ponownie bez określania parametru

Teraz, gdy masz ustawioną wartość domyślną, możesz uruchomić workflow ponownie bez konieczności określania wartości w wierszu poleceń.

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

Wyjście będzie w tym samym miejscu co poprzednio, ale zawartość powinna być zaktualizowana o nowy tekst.

??? abstract "Zawartość pliku"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow użył domyślnej wartości parametru greeting do utworzenia wyjścia.

#### 3.4.3. Nadpisz wartość domyślną

Jeśli podasz parametr w wierszu poleceń, wartość CLI nadpisze wartość domyślną.

Wypróbuj:

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

Po raz kolejny powinieneś znaleźć odpowiednie zaktualizowane wyjście w Swoim katalogu wyników.

??? abstract "Zawartość pliku"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Uwaga"

    W Nextflow istnieje wiele miejsc, gdzie możesz określić wartości dla parametrów.
    Jeśli ten sam parametr jest ustawiony na różne wartości w wielu miejscach, Nextflow określi, jakiej wartości użyć, na podstawie kolejności pierwszeństwa opisanej [tutaj](https://www.nextflow.io/docs/latest/config.html).

    Omówimy to bardziej szczegółowo w Części 6 (Configuration).

### Podsumowanie

Wiesz, jak używać prostego zmiennego wejścia dostarczanego w czasie wykonywania przez parametr wiersza poleceń, a także jak konfigurować, używać i nadpisywać wartości domyślne.

### Co dalej?

Naucz się wygodniej zarządzać wykonywaniem.

---

## 4. Zarządzaj wykonywaniem workflow'u

Wiedza o tym, jak uruchamiać workflow'y i pobierać wyjścia, jest świetna, ale szybko przekonasz się, że jest kilka innych aspektów zarządzania workflow'em, które ułatwią Ci życie, szczególnie jeśli tworzysz własne workflow'y.

Tutaj pokażemy Ci, jak używać funkcji `resume`, gdy musisz ponownie uruchomić ten sam workflow, jak przeglądać dziennik poprzednich wykonań za pomocą `nextflow log` i jak usuwać starsze katalogi work za pomocą `nextflow clean`.

<!-- Any other cool options we should include? Added log -->

### 4.1. Uruchom ponownie workflow z `-resume`

Czasami zechcesz ponownie uruchomić pipeline bez powtarzania kroków zakończonych pomyślnie.

Nextflow ma opcję o nazwie `-resume`, która to umożliwia.
W tym trybie wszelkie procesy wykonane wcześniej z dokładnie tym samym kodem, ustawieniami i wejściami zostaną pominięte.
Oznacza to, że Nextflow uruchomi tylko procesy dodane lub zmodyfikowane od ostatniego uruchomienia, lub te dla których dostarczasz nowe ustawienia.

Są dwie kluczowe zalety takiego postępowania:

- W fazie rozwoju możesz szybciej iterować, uruchamiając tylko te procesy, nad którymi aktywnie pracujesz.
- W produkcji, gdy coś pójdzie nie tak, często wystarczy naprawić problem i wznowić wykonanie od punktu awarii, oszczędzając czas i zasoby obliczeniowe.

Aby jej użyć, po prostu dodaj `-resume` do Swojego polecenia i uruchom go:

```bash
nextflow run hello-world.nf -resume
```

??? success "Wyjście polecenia"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

Wyjście konsoli powinno wyglądać znajomo, ale jest jedna rzecz, która jest nieco inna w porównaniu do poprzedniego.

Szukaj części `cached:`, która została dodana w linii statusu procesu (linia 5), co oznacza, że Nextflow rozpoznał, że już wykonał tę pracę i po prostu ponownie użył wyniku z poprzedniego pomyślnego uruchomienia.

Możesz również zobaczyć, że hash podkatalogu work jest taki sam jak w poprzednim uruchomieniu.
Nextflow dosłownie wskazuje Ci poprzednie wykonanie i mówi "Już to zrobiłem tam."

!!! tip "Wskazówka"

    Gdy ponownie uruchamiasz pipeline z `resume`, Nextflow nie nadpisuje żadnych plików opublikowanych poza katalogiem work przez żadne wykonania, które zostały pomyślnie uruchomione wcześniej.

### 4.2. Przeglądaj dziennik poprzednich wykonań

Niezależnie od tego, czy tworzysz nowy pipeline, czy uruchamiasz pipeline'y w produkcji, w pewnym momencie prawdopodobnie będziesz musiał wyszukać informacje o poprzednich uruchomieniach.
Oto jak to zrobić.

Za każdym razem, gdy uruchamiasz workflow nextflow, linia jest zapisywana do pliku dziennika o nazwie `history`, w ukrytym katalogu o nazwie `.nextflow` w bieżącym katalogu roboczym.

??? abstract "Zawartość pliku"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Ten plik zawiera znacznik czasu, nazwę uruchomienia, status, ID rewizji, ID sesji i pełny wiersz poleceń dla każdego uruchomienia Nextflow, które zostało uruchomione z bieżącego katalogu roboczego.

Wygodniejszym sposobem uzyskania dostępu do tych informacji jest użycie polecenia `nextflow log`.

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

To wyświetli zawartość pliku dziennika w terminalu, wzbogaconą o linię nagłówka.

Zauważysz, że ID sesji zmienia się za każdym razem, gdy uruchamiasz nowe polecenie `nextflow run`, Z WYJĄTKIEM sytuacji, gdy używasz opcji `-resume`.
W takim przypadku ID sesji pozostaje takie samo.

Nextflow używa ID sesji do grupowania informacji o cache'owaniu uruchomień w katalogu `cache`, również znajdującym się w `.nextflow`.

### 4.3. Usuń starsze katalogi work

Podczas procesu tworzenia zazwyczaj uruchomisz Swój szkic pipeline'u wiele razy, co może prowadzić do nagromadzenia wielu plików w wielu podkatalogach.

Na szczęście Nextflow zawiera pomocne podpolecenie `clean`, które może automatycznie usunąć podkatalogi work dla poprzednich uruchomień, które już Cię nie interesują.

#### 4.3.1. Określ kryteria usuwania

Istnieje wiele [opcji](https://www.nextflow.io/docs/latest/reference/cli.html#clean) do określenia, co usunąć.

Tutaj pokazujemy przykład, który usuwa wszystkie podkatalogi z uruchomień przed danym uruchomieniem, określonym za pomocą jego nazwy uruchomienia.

Wyszukaj najnowsze pomyślne uruchomienie, w którym nie używałeś `-resume`; w naszym przypadku nazwa uruchomienia to `golden_cantor`.

Nazwa uruchomienia to generowany maszynowo dwuczęściowy ciąg pokazany w nawiasach kwadratowych w linii wyjścia konsoli `Launching (...)`.
Możesz również użyć dziennika Nextflow, aby wyszukać uruchomienie na podstawie jego znacznika czasu i/lub wiersza poleceń.

#### 4.3.2. Wykonaj próbne uruchomienie

Najpierw używamy flagi próbnego uruchomienia `-n`, aby sprawdzić, co zostanie usunięte przy danym poleceniu:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Wyjście polecenia"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Twoje wyjście będzie miało inne nazwy katalogów zadań i może mieć inną liczbę linii, ale powinno wyglądać podobnie do przykładu.

Jeśli nie widzisz żadnych wyświetlonych linii, albo nie podałeś prawidłowej nazwy uruchomienia, albo nie ma poprzednich uruchomień do usunięcia. Upewnij się, że zmienisz `golden_cantor` w przykładowym poleceniu na odpowiednią najnowszą nazwę uruchomienia w Twoim dzienniku.

#### 4.3.3. Kontynuuj z usuwaniem

Jeśli wyjście wygląda zgodnie z oczekiwaniami i chcesz kontynuować usuwanie, uruchom polecenie ponownie z flagą `-f` zamiast `-n`:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Wyjście polecenia"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Wyjście powinno być podobne do poprzedniego, ale teraz mówi 'Removed' zamiast 'Would remove'.
Zauważ, że to nie usuwa dwuznakowych podkatalogów (jak `a3/` powyżej), ale opróżnia ich zawartość.

!!! warning "Ostrzeżenie"

    Usunięcie podkatalogów work z poprzednich uruchomień usuwa je z cache'u Nextflow i usuwa wszelkie wyjścia, które były przechowywane w tych katalogach.
    Oznacza to, że psuje to zdolność Nextflow do wznawiania wykonywania bez ponownego uruchamiania odpowiednich procesów.

    Jesteś odpowiedzialny za zapisywanie wszelkich wyjść, na których Ci zależy lub na których planujesz polegać! To jest główny powód, dla którego wolimy używać trybu `copy` zamiast trybu `symlink` dla dyrektywy `publish`.

### Podsumowanie

Wiesz, jak publikować wyjścia do określonego katalogu, ponownie uruchamiać pipeline bez powtarzania kroków, które już zostały uruchomione w identyczny sposób, i używać polecenia `nextflow clean` do czyszczenia starych katalogów work.

Ogólniej rzecz biorąc, wiesz, jak interpretować prosty workflow Nextflow, zarządzać jego wykonywaniem i pobierać wyjścia.

### Co dalej?

Zrób sobie małą przerwę, zasłużyłeś na to!

Gdy będziesz gotowy, przejdź do [**Części 2: Hello Channels**](./02_hello_channels.md), aby nauczyć się, jak używać kanałów do zasilania wejść do workflow'u, co pozwoli Ci skorzystać z wbudowanego równoległości przepływu danych Nextflow i innych potężnych funkcji.

---

## Quiz

<quiz>
Jakie są minimalne wymagane komponenty procesu Nextflow?
- [ ] Tylko bloki wejścia i wyjścia
- [x] Bloki wyjścia i skryptu
- [ ] Bloki wejścia, wyjścia i skryptu
- [ ] Tylko blok skryptu

Dowiedz się więcej: [1.1.1. Definicja process](#111-definicja-process)
</quiz>

<quiz>
Jaki jest cel bloku output w procesie?
- [ ] Wypisywanie wyników do konsoli
- [ ] Zapisywanie plików do katalogu work
- [x] Deklarowanie oczekiwanych wyjść z procesu
- [ ] Definiowanie zmiennych środowiskowych

Dowiedz się więcej: [1.1.1. Definicja process](#111-definicja-process)
</quiz>

<quiz>
Jakiego polecenia używa się do uruchomienia workflow'u Nextflow?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Patrząc na katalog work zadania, który plik zawiera faktyczne polecenie, które zostało wykonane?

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

Dowiedz się więcej: [1.2.2. Znajdź wyjście i dzienniki w katalogu `work`](#122-znajdz-wyjscie-i-dzienniki-w-katalogu-work)
</quiz>

<quiz>
Co robi flaga `-resume`?
- [ ] Restartuje workflow od początku
- [ ] Wstrzymuje workflow
- [x] Pomija procesy, które już zakończyły się pomyślnie
- [ ] Tworzy kopię zapasową workflow'u

Dowiedz się więcej: [4.1. Uruchom ponownie workflow z `-resume`](#41-uruchom-ponownie-workflow-z--resume)
</quiz>

<quiz>
Jaki jest domyślny tryb publikowania wyjść workflow'u?
- [ ] Kopiowanie plików do katalogu wyjściowego
- [x] Tworzenie dowiązań symbolicznych w katalogu wyjściowym
- [ ] Przenoszenie plików do katalogu wyjściowego
- [ ] Kompresowanie plików w katalogu wyjściowym

Dowiedz się więcej: [2.3. Ustaw tryb publikowania na kopiowanie](#23-ustaw-tryb-publikowania-na-kopiowanie)
</quiz>

<quiz>
Jak przekazujesz wartość parametru do workflow'u Nextflow z wiersza poleceń?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Dowiedz się więcej: [3.2. Skonfiguruj parametr wiersza poleceń do przechwytywania danych wejściowych użytkownika](#32-skonfiguruj-parametr-wiersza-polecen-do-przechwytywania-danych-wejsciowych-uzytkownika)
</quiz>

<quiz>
Jak odwołujesz się do zmiennej wewnątrz bloku skryptu w Nextflow?
- [ ] Używa składni `%variable%`
- [x] Używa składni `#!groovy ${variable}`
- [ ] Używa składni `{{variable}}`
- [ ] Używa składni `[variable]`
</quiz>
