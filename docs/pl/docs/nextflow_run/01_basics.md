# Część 1: Uruchamianie podstawowych operacji

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej pierwszej części kursu szkoleniowego Nextflow Run zaczynamy od bardzo prostego, niezależnego od dziedziny przykładu Hello World, który wykorzystamy do zademonstrowania podstawowych operacji i wskazania odpowiadających im komponentów kodu Nextflow'a.

??? info "Czym jest przykład Hello World?"

    „Hello World!" to minimalistyczny przykład, który ma na celu zademonstrowanie podstawowej składni i struktury języka programowania lub frameworka.
    Przykład zazwyczaj polega na wypisaniu frazy „Hello, World!" do urządzenia wyjściowego, takiego jak konsola lub terminal, lub zapisaniu jej do pliku.

---

## 1. Uruchom Hello World bezpośrednio

Zademonstrujmy tę koncepcję prostym poleceniem, które uruchomimy bezpośrednio w terminalu, aby pokazać, co robi, zanim opakujemy je w Nextflow'a.

!!! tip "Wskazówka"

    Pamiętaj, że powinieneś teraz znajdować się w katalogu `nextflow-run/`, jak opisano na stronie [Pierwsze kroki](00_orientation.md).

### 1.1. Spraw, aby terminal powiedział hello

Uruchom następujące polecenie w swoim terminalu.

```bash
echo 'Hello World!'
```

??? success "Wyjście polecenia"

    ```console
    Hello World!
    ```

To wypisuje tekst 'Hello World' bezpośrednio w terminalu.

### 1.2. Zapisz wyjście do pliku

Uruchamianie pipeline'ów polega głównie na odczytywaniu danych z plików i zapisywaniu wyników do innych plików, więc zmodyfikujmy polecenie, aby zapisywało wyjście tekstowe do pliku, czyniąc przykład nieco bardziej praktycznym.

```bash
echo 'Hello World!' > output.txt
```

??? success "Wyjście polecenia"

    ```console

    ```

To nie wypisuje niczego do terminala.

### 1.3. Znajdź wyjście

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

Dowiedz się, czego potrzeba, aby uruchomić workflow Nextflow'a, który osiąga ten sam rezultat.

---

## 2. Uruchom workflow'a

Udostępniamy Ci skrypt workflow'a o nazwie `1-hello.nf`, który przyjmuje powitanie wejściowe przez argument wiersza poleceń o nazwie `--input` i tworzy plik tekstowy zawierający to powitanie.

Nie będziemy jeszcze przeglądać kodu; najpierw zobaczmy, jak wygląda jego uruchomienie.

### 2.1. Uruchom workflow'a i monitoruj wykonanie

W terminalu uruchom następujące polecenie.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Jeśli Twoje wyjście konsoli wygląda mniej więcej tak, to gratulacje, właśnie uruchomiłeś swój pierwszy workflow Nextflow'a!

Najważniejszym wyjściem jest tutaj ostatnia linia, która jest podświetlona w powyższym wyjściu:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

To mówi nam, że proces `sayHello` został pomyślnie wykonany raz (`1 of 1 ✔`).

To świetnie, ale możesz się zastanawiać: gdzie jest wyjście?

### 2.2. Znajdź plik wyjściowy w katalogu `results`

Ten workflow jest skonfigurowany tak, aby publikować swoje wyjście do katalogu wyników.
Jeśli spojrzysz na swój bieżący katalog, zobaczysz, że gdy uruchomiłeś workflow'a, Nextflow utworzył nowy katalog o nazwie `results`, a także podkatalog o nazwie `1-hello` pod nim, zawierający plik o nazwie `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Otwórz plik; jego zawartość powinna odpowiadać ciągowi znaków, który podałeś w wierszu poleceń.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Świetnie, nasz workflow zrobił to, co powinien!

### 2.3. Zapisz wyniki do innego katalogu

Domyślnie Nextflow zapisze wyjścia pipeline'u do katalogu o nazwie `results` w Twojej bieżącej ścieżce.
Aby zmienić miejsce, do którego publikowane są Twoje pliki, użyj flagi CLI `-output-dir` (lub `-o` w skrócie)

!!! danger "Uwaga"

    Zauważ, że `--input` ma dwa myślniki, a `-output-dir` ma jeden!
    Dzieje się tak, ponieważ `--input` to _parametr_ pipeline'u, a `-output-dir` to podstawowa flaga CLI Nextflow'a.
    Więcej o tym później.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Powinieneś zobaczyć, że Twoje wyjścia są teraz publikowane do katalogu o nazwie `hello_results` zamiast `results`:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Pliki w tym katalogu są dokładnie takie same jak wcześniej, różni się tylko katalog najwyższego poziomu.
Należy jednak pamiętać, że w obu przypadkach „opublikowany" wynik jest kopią (lub w niektórych przypadkach dowiązaniem symbolicznym) rzeczywistego wyjścia wytworzonego przez Nextflow'a podczas wykonywania workflow'a.

Teraz zajrzymy pod maskę, aby zobaczyć, gdzie Nextflow faktycznie wykonał pracę.

!!! Warning "Ostrzeżenie"

    Nie wszystkie workflow'y będą skonfigurowane do publikowania wyjść do katalogu wyników, a nazwy katalogów i struktura mogą być różne.
    Nieco dalej w tej sekcji pokażemy Ci, jak dowiedzieć się, gdzie określone jest to zachowanie.

### 2.4. Znajdź oryginalne wyjście i logi w katalogu `work/`

Gdy uruchamiasz workflow'a, Nextflow tworzy odrębny „katalog zadania" dla każdego wywołania każdego procesu w workflow'ie (= każdego kroku w pipeline'ie).
Dla każdego z nich przygotuje niezbędne wejścia, wykona odpowiednie instrukcje i zapisze wyjścia oraz pliki logów w tym jednym katalogu, który jest automatycznie nazywany przy użyciu skrótu, aby uczynić go unikalnym.

Wszystkie te katalogi zadań będą znajdować się w katalogu o nazwie `work` w Twoim bieżącym katalogu (gdzie uruchamiasz polecenie).

To może brzmieć zagmatwanie, więc zobaczmy, jak to wygląda w praktyce.

Wracając do wyjścia konsoli dla workflow'a, który uruchomiliśmy wcześniej, mieliśmy tę linię:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Widzisz, jak linia zaczyna się od `[a3/1e1535]`?
To jest skrócona forma ścieżki katalogu zadania dla tego jednego wywołania procesu i mówi Ci, gdzie znaleźć wyjście wywołania procesu `sayHello` w ścieżce katalogu `work/`.

Możesz znaleźć pełną ścieżkę, wpisując następujące polecenie (zastępując `a3/1e1535` tym, co widzisz w swoim własnym terminalu) i naciskając klawisz tab, aby automatycznie uzupełnić ścieżkę, lub dodając gwiazdkę:

```bash
ls work/a3/1e1535*
```

To powinno dać pełną ścieżkę katalogu: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Spójrzmy, co tam jest.

??? abstract "Zawartość katalogu"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
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

W `work/` znajdują się dwa zestawy katalogów, z dwóch różnych uruchomień pipeline'u, które wykonaliśmy.
Każde wykonanie zadania otrzymuje swój własny, izolowany katalog do pracy.
W tym przypadku pipeline zrobił to samo za każdym razem, więc zawartość każdego katalogu zadania jest identyczna.

Powinieneś natychmiast rozpoznać plik `output.txt`, który w rzeczywistości jest oryginalnym wyjściem procesu `sayHello`, które zostało opublikowane w katalogu `results`.
Jeśli go otworzysz, znajdziesz ponownie powitanie `Hello World!`.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

A co z tymi wszystkimi innymi plikami?

To pliki pomocnicze i logi, które Nextflow zapisał w ramach wykonania zadania:

- **`.command.begin`**: Plik wartowniczy utworzony natychmiast po uruchomieniu zadania.
- **`.command.err`**: Komunikaty o błędach (`stderr`) wyemitowane przez wywołanie procesu
- **`.command.log`**: Kompletne wyjście logu wyemitowane przez wywołanie procesu
- **`.command.out`**: Regularne wyjście (`stdout`) wywołania procesu
- **`.command.run`**: Pełny skrypt uruchomiony przez Nextflow'a w celu wykonania wywołania procesu
- **`.command.sh`**: Polecenie, które faktycznie zostało uruchomione przez wywołanie procesu
- **`.exitcode`**: Kod wyjścia wynikający z polecenia

Plik `.command.sh` jest szczególnie przydatny, ponieważ pokazuje główne polecenie wykonane przez Nextflow'a, nie uwzględniając całej księgowości i konfiguracji zadania/środowiska.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

To potwierdza, że workflow skomponował to samo polecenie, które uruchomiliśmy bezpośrednio w wierszu poleceń wcześniej.

Gdy coś pójdzie nie tak i musisz rozwiązać problem, może być przydatne spojrzenie na skrypt `command.sh`, aby sprawdzić dokładnie, jakie polecenie skomponował Nextflow na podstawie instrukcji workflow'a, interpolacji zmiennych i tak dalej.

### 2.5. Uruchom ponownie workflow'a z różnymi powitaniami

Spróbuj uruchomić ponownie workflow'a kilka razy z różnymi wartościami dla argumentu `--input`, a następnie spójrz na katalogi zadań.

??? abstract "Zawartość katalogu"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Widzisz, że dla każdego uruchomienia został utworzony nowy podkatalog z kompletnym zestawem plików wyjściowych i logów.

W przeciwieństwie do tego, jeśli spojrzysz na katalog `results`, nadal jest tam tylko jeden zestaw wyników, a zawartość pliku wyjściowego odpowiada temu, co uruchomiłeś ostatnio.

??? abstract "Zawartość katalogu"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

To pokazuje, że opublikowane wyniki zostaną nadpisane przez kolejne wykonania, podczas gdy katalogi zadań w `work/` są zachowywane.

### Podsumowanie

Wiesz, jak uruchomić prosty skrypt Nextflow'a, monitorować jego wykonanie i znaleźć jego wyjścia.

### Co dalej?

Naucz się czytać podstawowy skrypt Nextflow'a i identyfikować, jak jego komponenty odnoszą się do jego funkcjonalności.

---

## 3. Zbadaj startowy skrypt workflow'a Hello World

To, co tam zrobiliśmy, to w zasadzie potraktowanie skryptu workflow'a jak czarnej skrzynki.
Teraz, gdy zobaczyliśmy, co robi, otwórzmy skrzynkę i zajrzyjmy do środka.

Naszym celem nie jest zapamiętanie składni kodu Nextflow'a, ale wyrobienie podstawowej intuicji, jakie są główne komponenty i jak są zorganizowane.

### 3.1. Zbadaj ogólną strukturę kodu

Znajdziesz skrypt `1-hello.nf` w swoim bieżącym katalogu, którym powinien być `nextflow-run`. Otwórz go w panelu edytora.

??? full-code "Pełny plik kodu"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo, aby wypisać 'Hello World!' do pliku
    */
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }

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

Skrypt workflow'a Nextflow'a zazwyczaj zawiera jedną lub więcej definicji **process**, sam **workflow** oraz kilka opcjonalnych bloków, takich jak **params** i **output**.

Każdy **process** opisuje, jakie operacje powinien wykonać odpowiadający mu krok w pipeline'ie, podczas gdy **workflow** opisuje logikę przepływu danych, która łączy różne kroki.

Przyjrzyjmy się najpierw bliżej blokowi **process**, a następnie spojrzymy na blok **workflow**.

### 3.2. Definicja `process`

Pierwszy blok kodu opisuje [**process**](https://nextflow.io/docs/latest/process.html).
Definicja procesu zaczyna się od słowa kluczowego `process`, po którym następuje nazwa procesu, a na końcu ciało procesu ograniczone nawiasami klamrowymi.
Ciało procesu musi zawierać blok skryptu, który określa polecenie do uruchomienia, które może być wszystkim, co byłbyś w stanie uruchomić w terminalu wiersza poleceń.

```groovy title="1-hello.nf" linenums="3"
/*
* Użyj echo, aby wypisać powitanie do pliku
*/
process sayHello {

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '${greeting}' > output.txt
    """
}
```

Tutaj mamy **process** o nazwie `sayHello`, który przyjmuje zmienną **input** o nazwie `greeting` i zapisuje swoje **output** do pliku o nazwie `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

To bardzo minimalna definicja procesu, która zawiera tylko definicję `input`, definicję `output` i `script` do wykonania.

Definicja `input` zawiera kwalifikator `val`, który mówi Nextflow'owi, aby oczekiwał wartości jakiegoś rodzaju (może to być ciąg znaków, liczba, cokolwiek).

Definicja `output` zawiera kwalifikator `path`, który mówi Nextflow'owi, że powinno to być traktowane jako ścieżka (obejmuje zarówno ścieżki katalogów, jak i pliki).

### 3.3. Definicja `workflow`

Drugi blok kodu opisuje sam [**workflow**](https://nextflow.io/docs/latest/workflow.html).
Definicja workflow'a zaczyna się od słowa kluczowego `workflow`, po którym następuje opcjonalna nazwa, a następnie ciało workflow'a ograniczone nawiasami klamrowymi.

Tutaj mamy **workflow**, który składa się z bloku `main:` i bloku `publish:`.
Blok `main:` jest głównym ciałem workflow'a, a blok `publish:` wymienia wyjścia, które powinny zostać opublikowane w katalogu `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // Wyemituj powitanie
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

W tym przypadku blok `main:` zawiera wywołanie procesu `sayHello` i podaje mu wejście o nazwie `params.input` do użycia jako powitanie.

Jak omówimy bardziej szczegółowo za chwilę, `params.input` przechowuje wartość, którą podaliśmy parametrowi `--input` w naszym wierszu poleceń.

Blok `publish:` wymienia wyjście wywołania procesu `sayHello()`, do którego odnosi się jako `sayHello.out` i nadaje nazwę `first_output` (może to być cokolwiek, czego chce autor workflow'a).

To bardzo minimalna definicja **workflow**.
W rzeczywistym pipeline'ie workflow zazwyczaj zawiera wiele wywołań **procesów** połączonych **kanałami** i mogą być ustawione wartości domyślne dla zmiennych wejściowych.

Zajmiemy się tym w części 2 kursu.
Na razie przyjrzyjmy się bliżej, jak nasz workflow obsługuje wejścia i wyjścia.

### 3.4. System `params` parametrów wiersza poleceń

`params.input`, które dostarczamy do wywołania procesu `sayHello()`, to zgrabny fragment kodu Nextflow'a i warto poświęcić mu dodatkową minutę.

Jak wspomniano powyżej, w ten sposób przekazujemy wartość parametru wiersza poleceń `--input` do wywołania procesu `sayHello()`.
W rzeczywistości samo zadeklarowanie `params.someParameterName` wystarczy, aby nadać workflow'owi parametr o nazwie `--someParameterName` z wiersza poleceń.

Tutaj sformalizowaliśmy tę deklarację parametru, konfigurując blok `params`, który określa typ wejścia oczekiwanego przez workflow (Nextflow 25.10.2 i nowsze).

```groovy title="1-hello.nf" linenums="20"
/*
 * Parametry pipeline'u
 */
params {
    input: String
}
```

Obsługiwane typy obejmują `String`, `Integer`, `Float`, `Boolean` i `Path`.
Aby dowiedzieć się więcej, zobacz [Parametry workflow'a](https://nextflow.io/docs/latest/config.html#workflow-parameters) w dokumentacji referencyjnej Nextflow'a.

!!! tip "Wskazówka"

    Pamiętaj, że parametry _workflow'a_ zadeklarowane przy użyciu systemu `params` zawsze przyjmują dwa myślniki w wierszu poleceń (`--`).
    To odróżnia je od flag CLI _poziomu Nextflow'a_, które przyjmują tylko jeden myślnik (`-`).

### 3.5. Dyrektywa `publish`

Na drugim końcu workflow'a już rzuciliśmy okiem na blok `publish:`.
To jedna połowa systemu obsługi wyjścia; druga połowa to blok `output` znajdujący się poniżej.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

To określa, że wyjście `first_output` wymienione w bloku `publish:` powinno zostać skopiowane do podkatalogu o nazwie `1-hello` w domyślnym katalogu wyjściowym `results`.

Linia `mode 'copy'` nadpisuje domyślne zachowanie systemu, którym jest utworzenie dowiązania symbolicznego (lub symlinku) do oryginalnego pliku w katalogu `work/` zamiast właściwej kopii.

Istnieje więcej opcji niż wyświetlone tutaj do kontrolowania zachowania publikowania; omówimy kilka z nich później.
Zobaczysz również, że gdy workflow generuje wiele wyjść, każde z nich jest wymienione w ten sposób w bloku `output`.

Aby dowiedzieć się więcej, zobacz [Publikowanie wyjść](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) w dokumentacji referencyjnej Nextflow'a.

??? info "Starsza składnia publikowania wyjść przy użyciu `publishDir`"

    Do niedawna ustalonym sposobem publikowania wyjść było robienie tego na poziomie każdego indywidualnego procesu przy użyciu dyrektywy `publishDir`.

    Nadal znajdziesz ten wzorzec kodu wszędzie w starszych pipeline'ach Nextflow'a i modułach procesów, więc ważne jest, aby być tego świadomym.

    Zamiast mieć blok `publish:` w workflow'ie i blok `output` na najwyższym poziomie, zobaczyłbyś linię `publishDir` w definicji procesu `sayHello`:

    ```groovy title="Przykład składni" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Jednak nie zalecamy używania tego w żadnej nowej pracy, ponieważ w przyszłych wersjach języka Nextflow będzie to ostatecznie niedozwolone.

### Podsumowanie

Wiesz już, jak zbudowany jest prosty workflow Nextflow'a i jak podstawowe komponenty odnoszą się do jego funkcjonalności.

### Co dalej?

Naucz się wygodnie zarządzać wykonaniami workflow'a.

---

## 4. Zarządzaj wykonaniami workflow'a

Wiedza o tym, jak uruchamiać workflow'y i pobierać wyjścia, jest świetna, ale szybko odkryjesz, że jest kilka innych aspektów zarządzania workflow'ami, które ułatwią Ci życie.

Tutaj pokażemy Ci, jak wykorzystać funkcję `resume` na wypadek, gdybyś musiał ponownie uruchomić ten sam workflow, jak sprawdzić logi wykonania za pomocą `nextflow log` i jak usunąć starsze katalogi robocze za pomocą `nextflow clean`.

### 4.1. Uruchom ponownie workflow'a z `-resume`

Czasami będziesz chciał ponownie uruchomić pipeline, który już wcześniej uruchomiłeś, bez powtarzania jakiejkolwiek pracy, która została już pomyślnie ukończona.

Nextflow ma opcję o nazwie `-resume`, która pozwala Ci to zrobić.
W szczególności w tym trybie wszystkie procesy, które zostały już uruchomione z dokładnie tym samym kodem, ustawieniami i wejściami, zostaną pominięte.
Oznacza to, że Nextflow uruchomi tylko procesy, które dodałeś lub zmodyfikowałeś od ostatniego uruchomienia, lub którym dostarczasz nowe ustawienia lub wejścia.

Istnieją dwie kluczowe zalety tego podejścia:

- Jeśli jesteś w trakcie tworzenia pipeline'u, możesz iterować szybciej, ponieważ musisz uruchomić tylko proces(y), nad którymi aktywnie pracujesz, aby przetestować swoje zmiany.
- Jeśli uruchamiasz pipeline w produkcji i coś pójdzie nie tak, w wielu przypadkach możesz naprawić problem i ponownie uruchomić pipeline, a on wznowi działanie od punktu awarii, co może zaoszczędzić Ci dużo czasu i mocy obliczeniowej.

Aby z tego skorzystać, po prostu dodaj `-resume` do swojego polecenia i uruchom je:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

Wyjście konsoli powinno wyglądać znajomo, ale jest jedna rzecz, która jest nieco inna niż wcześniej.

Poszukaj fragmentu `cached:`, który został dodany w linii statusu procesu (linia 5), co oznacza, że Nextflow rozpoznał, że już wykonał tę pracę i po prostu ponownie wykorzystał wynik z poprzedniego pomyślnego uruchomienia.

Możesz również zobaczyć, że skrót podkatalogu roboczego jest taki sam jak w poprzednim uruchomieniu.
Nextflow dosłownie wskazuje Ci poprzednie wykonanie i mówi: „Już to zrobiłem tam".

!!! tip "Wskazówka"

    Gdy ponownie uruchamiasz pipeline z `resume`, Nextflow nie nadpisuje żadnych plików opublikowanych poza katalogiem roboczym przez żadne wykonania, które zostały wcześniej pomyślnie uruchomione.

    Aby dowiedzieć się więcej, zobacz [Cache i resume](https://nextflow.io/docs/latest/cache-and-resume.html) w dokumentacji referencyjnej Nextflow'a.

### 4.2. Sprawdź log przeszłych wykonań

Za każdym razem, gdy uruchamiasz workflow Nextflow'a, linia zostaje zapisana do pliku logu o nazwie `history`, w ukrytym katalogu o nazwie `.nextflow` w bieżącym katalogu roboczym.

??? abstract "Zawartość pliku"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Ten plik podaje znacznik czasu, nazwę uruchomienia, status, ID rewizji, ID sesji i pełny wiersz poleceń dla każdego uruchomienia Nextflow'a, które zostało uruchomione z bieżącego katalogu roboczego.

Wygodniejszym sposobem dostępu do tych informacji jest użycie polecenia [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

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

Nextflow używa ID sesji do grupowania informacji o buforowaniu uruchomień w katalogu `cache`, również znajdującym się w `.nextflow`.

### 4.3. Usuń starsze katalogi robocze

Jeśli uruchamiasz wiele pipeline'ów, możesz w końcu zgromadzić bardzo wiele plików w wielu podkatalogach.
Ponieważ podkatalogi są nazywane losowo, trudno powiedzieć po ich nazwach, które są starsze, a które nowsze.

Na szczęście Nextflow zawiera pomocne polecenie o nazwie [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean), które może automatycznie usunąć podkatalogi robocze z przeszłych uruchomień, na których Ci już nie zależy.

#### 4.3.1. Określ kryteria usuwania

Istnieje wiele opcji określania, co usunąć, które możesz zbadać w dokumentacji linkowanej powyżej.
Tutaj pokazujemy przykład, który usuwa wszystkie podkatalogi z uruchomień przed danym uruchomieniem, określonym przy użyciu jego nazwy uruchomienia.

Wyszukaj najnowsze pomyślne uruchomienie, w którym nie użyłeś `-resume`; w naszym przypadku nazwa uruchomienia to `backstabbing_swartz`.

Nazwa uruchomienia to wygenerowany maszynowo dwuczęściowy ciąg pokazany w nawiasach kwadratowych w linii wyjścia konsoli `Launching (...)`.
Możesz również użyć logu Nextflow'a, aby wyszukać uruchomienie na podstawie jego znacznika czasu i/lub wiersza poleceń.

#### 4.3.2. Wykonaj próbne uruchomienie

Najpierw używamy flagi próbnego uruchomienia `-n`, aby sprawdzić, co zostanie usunięte przy danym poleceniu:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Wyjście polecenia"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Twoje wyjście będzie miało różne nazwy katalogów zadań i może mieć inną liczbę linii, ale powinno wyglądać podobnie do przykładu.

Jeśli nie widzisz żadnych linii wyjściowych, albo nie podałeś prawidłowej nazwy uruchomienia, albo nie ma przeszłych uruchomień do usunięcia. Upewnij się, że zmieniłeś `backstabbing_swartz` w przykładowym poleceniu na odpowiadającą najnowszą nazwę uruchomienia w Twoim logu.

#### 4.3.3. Przystąp do usuwania

Jeśli wyjście wygląda zgodnie z oczekiwaniami i chcesz przystąpić do usuwania, uruchom ponownie polecenie z flagą `-f` zamiast `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Wyjście polecenia"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Wyjście powinno być podobne do poprzedniego, ale teraz mówiące 'Removed' zamiast 'Would remove'.
Zauważ, że to nie usuwa dwuznakowych podkatalogów (jak `eb/` powyżej), ale opróżnia ich zawartość.

!!! Warning "Ostrzeżenie"

    Usuwanie podkatalogów roboczych z przeszłych uruchomień usuwa je z bufora Nextflow'a i usuwa wszelkie wyjścia, które były przechowywane w tych katalogach.
    Oznacza to, że przerywa zdolność Nextflow'a do wznowienia wykonania bez ponownego uruchamiania odpowiednich procesów.

    Jesteś odpowiedzialny za zapisanie wszelkich wyjść, na których Ci zależy! To główny powód, dla którego wolimy używać trybu `copy` zamiast trybu `symlink` dla dyrektywy `publish`.

### Podsumowanie

Wiesz, jak ponownie uruchomić pipeline bez powtarzania kroków, które zostały już uruchomione w identyczny sposób, sprawdzić log wykonania i użyć polecenia `nextflow clean` do czyszczenia starych katalogów roboczych.

### Co dalej?

Zrób sobie małą przerwę! Właśnie przyswaiłeś elementy składowe składni Nextflow'a i podstawowe instrukcje użytkowania.

W następnej sekcji tego szkolenia przyjrzymy się czterem kolejno bardziej realistycznym wersjom pipeline'u Hello World, które pokażą, jak Nextflow pozwala efektywnie przetwarzać wiele wejść, uruchamiać workflow'y składające się z wielu połączonych ze sobą kroków, wykorzystywać modułowe komponenty kodu oraz używać kontenerów dla większej powtarzalności i przenośności.

---

## Quiz

<quiz>
W linii wyjścia konsoli `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, co reprezentuje `[a3/7be2fa]`?
- [ ] Numer wersji procesu
- [ ] Unikalny identyfikator uruchomienia
- [x] Skróconą ścieżkę do katalogu roboczego zadania
- [ ] Sumę kontrolną pliku wyjściowego

Dowiedz się więcej: [2.3. Znajdź oryginalne wyjście i logi w katalogu `work/`](#23-znajdz-oryginalne-wyjscie-i-logi-w-katalogu-work)
</quiz>

<quiz>
Jaki jest cel pliku `.command.sh` w katalogu zadania?
- [ ] Przechowuje ustawienia konfiguracyjne zadania
- [x] Pokazuje faktyczne polecenie, które zostało wykonane przez proces
- [ ] Zawiera komunikaty o błędach z nieudanych zadań
- [ ] Wymienia pliki wejściowe przygotowane dla zadania

Dowiedz się więcej: [2.3. Znajdź oryginalne wyjście i logi w katalogu `work/`](#23-znajdz-oryginalne-wyjscie-i-logi-w-katalogu-work)
</quiz>

<quiz>
Co dzieje się z opublikowanymi wynikami, gdy ponownie uruchamiasz workflow'a bez `-resume`?
- [ ] Są zachowywane w oddzielnych katalogach z znacznikami czasu
- [x] Zostają nadpisane przez nowe wykonanie
- [ ] Nextflow zapobiega nadpisywaniu i kończy się niepowodzeniem
- [ ] Są automatycznie archiwizowane

Dowiedz się więcej: [2.4. Uruchom ponownie workflow'a z różnymi powitaniami](#24-uruchom-ponownie-workflow-a-z-roznymi-powitaniami)
</quiz>

<quiz>
Co wskazuje to wyjście konsoli?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Zadanie nie powiodło się i zostało pominięte
- [ ] Zadanie czeka w kolejce
- [x] Nextflow ponownie wykorzystał wyniki z poprzedniego identycznego wykonania
- [ ] Zadanie zostało ręcznie anulowane

Dowiedz się więcej: [4.1. Uruchom ponownie workflow'a z `-resume`](#41-uruchom-ponownie-workflow-a-z--resume)
</quiz>

<quiz>
Gdzie Nextflow przechowuje historię wykonań, którą wyświetla polecenie `nextflow log`?
- [ ] W katalogu wyników
- [ ] W katalogu roboczym
- [x] W pliku `.nextflow/history`
- [ ] W `nextflow.config`

Dowiedz się więcej: [4.2. Sprawdź log przeszłych wykonań](#42-sprawdz-log-przeszlych-wykonan)
</quiz>

<quiz>
Jaki jest cel bloku `params` w pliku workflow'a?
- [ ] Definiowanie wymagań zasobowych procesu
- [ ] Konfigurowanie executora
- [x] Deklarowanie i typowanie parametrów wejściowych workflow'a
- [ ] Określanie opcji publikowania wyjść

Dowiedz się więcej: [3.4. System `params` parametrów wiersza poleceń](#34-system-params-parametrow-wiersza-polecen)
</quiz>

<quiz>
W bloku `output` workflow'a, co robi `mode 'copy'`?
- [ ] Tworzy kopię zapasową katalogu roboczego
- [x] Tworzy pełną kopię plików zamiast dowiązań symbolicznych
- [ ] Kopiuje skrypt workflow'a do wyników
- [ ] Włącza przyrostowe kopiowanie plików

Dowiedz się więcej: [3.5. Dyrektywa `publish`](#35-dyrektywa-publish)
</quiz>

<quiz>
Jaka jest zalecana flaga do użycia z poleceniem `nextflow clean` przed faktycznym usunięciem plików?
- [x] `-n` (próbne uruchomienie) do podglądu tego, co zostanie usunięte
- [ ] `-v` (szczegółowe) do wyświetlenia szczegółowego wyjścia
- [ ] `-a` (wszystko) do wybrania wszystkich katalogów
- [ ] `-q` (ciche) do pominięcia ostrzeżeń

Dowiedz się więcej: [4.3. Usuń starsze katalogi robocze](#43-usun-starsze-katalogi-robocze)
</quiz>
