# Część 1: Podstawowe operacje

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W tej pierwszej części kursu szkoleniowego Nextflow Run, łagodnie wprowadzamy temat za pomocą bardzo podstawowego, niezależnego od domeny przykładu Hello World, którego użyjemy do zademonstrowania podstawowych operacji i wskazania odpowiednich komponentów kodu Nextflow.

??? info "Czym jest przykład Hello World?"

    "Hello World!" to minimalistyczny przykład, który ma na celu zademonstrowanie podstawowej składni i struktury języka programowania lub frameworka oprogramowania.
    Przykład zazwyczaj polega na wypisaniu frazy "Hello, World!" na urządzenie wyjściowe, takie jak konsola lub terminal, lub zapisaniu jej do pliku.

---

## 1. Uruchom Hello World bezpośrednio

Zademonstrujmy tę koncepcję prostym poleceniem, które uruchamiamy bezpośrednio w terminalu, aby pokazać, co robi, zanim opakujemy je w Nextflow.

!!! tip "Wskazówka"

    Pamiętaj, że powinieneś teraz znajdować się w katalogu `nextflow-run/`, jak opisano na stronie [Rozpoczęcie pracy](00_orientation.md).

### 1.1. Spraw, aby terminal powiedział cześć

Uruchom następujące polecenie w terminalu.

```bash
echo 'Hello World!'
```

??? success "Wyjście polecenia"

    ```console
    Hello World!
    ```

To wypisuje tekst 'Hello World' bezpośrednio w terminalu.

### 1.2. Zapisz wyjście do pliku

Uruchamianie pipeline głównie polega na odczytywaniu danych z plików i zapisywaniu wyników do innych plików, więc zmodyfikujmy polecenie, aby zapisać wyjście tekstowe do pliku, co uczyni przykład bardziej odpowiednim.

```bash
echo 'Hello World!' > output.txt
```

??? success "Wyjście polecenia"

    ```console

    ```

To nie wypisuje niczego do terminala.

### 1.3. Znajdź wyjście

Tekst 'Hello World' powinien teraz znajdować się w pliku wyjściowym, który określiliśmy, o nazwie `output.txt`.
Możesz go otworzyć w eksploratorze plików lub z wiersza poleceń, używając na przykład narzędzia `cat`.

??? abstract "Zawartość pliku"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

To jest to, co będziemy próbowali odtworzyć za pomocą naszego pierwszego workflow Nextflow.

### Podsumowanie

Teraz wiesz, jak uruchomić proste polecenie w terminalu, które wypisuje tekst, i opcjonalnie, jak sprawić, aby zapisało wyjście do pliku.

### Co dalej?

Dowiedz się, co trzeba zrobić, aby uruchomić workflow Nextflow, który osiąga ten sam wynik.

---

## 2. Uruchom workflow

Dostarczamy skrypt workflow o nazwie `1-hello.nf`, który przyjmuje powitanie wejściowe przez argument wiersza poleceń o nazwie `--input` i produkuje plik tekstowy zawierający to powitanie.

Nie będziemy jeszcze patrzeć na kod; najpierw zobaczmy, jak wygląda jego uruchomienie.

### 2.1. Uruchom workflow i monitoruj wykonanie

W terminalu uruchom następujące polecenie:

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

Jeśli wyjście konsoli wygląda mniej więcej tak, to gratulacje, właśnie uruchomiłeś swój pierwszy workflow Nextflow!

Najważniejsze wyjście tutaj to ostatnia linia, która jest podświetlona w powyższym wyjściu:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

To mówi nam, że process `sayHello` został pomyślnie wykonany raz (`1 of 1 ✔`).

Świetnie, ale możesz się zastanawiać: gdzie jest wyjście?

### 2.2. Znajdź plik wyjściowy w katalogu `results`

Ten workflow jest skonfigurowany do publikowania swojego wyjścia w katalogu results.
Jeśli spojrzysz na swój bieżący katalog, zobaczysz, że gdy uruchomiłeś workflow, Nextflow utworzył nowy katalog o nazwie `results`, a także podkatalog o nazwie `1-hello` pod nim, zawierający plik o nazwie `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Otwórz plik; zawartość powinna odpowiadać ciągowi znaków, który podałeś w wierszu poleceń.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Świetnie, nasz workflow zrobił to, co miał zrobić!

Jednak pamiętaj, że 'opublikowany' wynik jest kopią (lub w niektórych przypadkach dowiązaniem symbolicznym) rzeczywistego wyjścia wyprodukowanego przez Nextflow podczas wykonywania workflow.

Więc teraz zajrzymy pod maskę, aby zobaczyć, gdzie Nextflow faktycznie wykonał pracę.

!!! warning "Ostrzeżenie"

    Nie wszystkie workflow będą skonfigurowane do publikowania wyjść do katalogu results, i/lub nazwy katalogów i struktura mogą być inne.
    Trochę dalej w tej sekcji pokażemy, jak dowiedzieć się, gdzie to zachowanie jest określone.

### 2.3. Znajdź oryginalne wyjście i logi w katalogu `work/`

Gdy uruchamiasz workflow, Nextflow tworzy odrębny 'katalog zadania' dla każdego pojedynczego wywołania każdego process w workflow (=każdego kroku w pipeline).
Dla każdego z nich przygotowuje niezbędne dane wejściowe, wykonuje odpowiednie instrukcje i zapisuje wyjścia i pliki dziennika w tym jednym katalogu, który jest automatycznie nazywany przy użyciu skrótu hash, aby uczynić go unikalnym.

Wszystkie te katalogi zadań będą znajdować się w katalogu o nazwie `work` w bieżącym katalogu (skąd uruchamiasz polecenie).

To może brzmieć myląco, więc zobaczmy, jak to wygląda w praktyce.

Wracając do wyjścia konsoli dla workflow, który uruchomiliśmy wcześniej, mieliśmy tę linię:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Widzisz, jak linia zaczyna się od `[a3/7be2fa]`?
To jest skrócona forma ścieżki katalogu zadania dla tego jednego wywołania process i mówi, gdzie znaleźć wyjście wywołania `sayHello` process w ścieżce katalogu `work/`.

Możesz znaleźć pełną ścieżkę, wpisując następujące polecenie (zastępując `a3/7be2fa` tym, co widzisz w swoim terminalu) i naciskając klawisz tab, aby automatycznie uzupełnić ścieżkę, lub dodając gwiazdkę:

```bash
ls work/a3/7be2fa*
```

To powinno dać pełną ścieżkę katalogu: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Zobaczmy, co jest w środku.

??? abstract "Zawartość katalogu"

    ```console
    work
    └── a3
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

    Dokładne nazwy podkatalogów będą inne w twoim systemie.

    Jeśli przeglądasz zawartość podkatalogu zadania w eksploratorze plików VSCode, zobaczysz wszystkie pliki od razu.
    Jednak pliki dziennika są ustawione jako niewidoczne w terminalu, więc jeśli chcesz użyć `ls` lub `tree` do ich wyświetlenia, musisz ustawić odpowiednią opcję wyświetlania niewidocznych plików.

    ```bash
    tree -a work
    ```

Powinieneś natychmiast rozpoznać plik `output.txt`, który jest w rzeczywistości oryginalnym wyjściem process `sayHello`, które zostało opublikowane do katalogu `results`.
Jeśli go otworzysz, znajdziesz ponownie powitanie `Hello World!`.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt"
Hello World!
```

A co z wszystkimi tymi innymi plikami?

To są pliki pomocnicze i dziennika, które Nextflow zapisał jako część wykonania zadania:

- **`.command.begin`**: Plik wartowniczy tworzony zaraz po uruchomieniu zadania.
- **`.command.err`**: Komunikaty o błędach (`stderr`) emitowane przez wywołanie process
- **`.command.log`**: Kompletne wyjście dziennika emitowane przez wywołanie process
- **`.command.out`**: Zwykłe wyjście (`stdout`) przez wywołanie process
- **`.command.run`**: Pełny skrypt uruchomiony przez Nextflow do wykonania wywołania process
- **`.command.sh`**: Polecenie, które zostało faktycznie uruchomione przez wywołanie process
- **`.exitcode`**: Kod wyjścia wynikający z polecenia

Plik `.command.sh` jest szczególnie przydatny, ponieważ pokazuje główne polecenie, które wykonał Nextflow, nie włączając całej księgowości i konfiguracji zadania/środowiska.

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

To potwierdza, że workflow złożył to samo polecenie, które uruchomiliśmy bezpośrednio w wierszu poleceń wcześniej.

Gdy coś pójdzie nie tak i musisz rozwiązać problem, przydatne może być spojrzenie na skrypt `command.sh`, aby sprawdzić dokładnie, jakie polecenie Nextflow złożył na podstawie instrukcji workflow, interpolacji zmiennych itd.

### 2.4. Uruchom ponownie workflow z różnymi powitaniami

Spróbuj uruchomić workflow kilka razy z różnymi wartościami argumentu `--input`, a następnie spójrz na katalogi zadań.

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
    ├── 67
    │   ├── 134e6317f90726c6c17ad53234a32b
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

W przeciwieństwie do tego, jeśli spojrzysz na katalog `results`, nadal jest tylko jeden zestaw wyników, a zawartość pliku wyjściowego odpowiada temu, co uruchomiłeś ostatnio.

??? abstract "Zawartość katalogu"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

To pokazuje, że opublikowane wyniki zostaną nadpisane przez kolejne wykonania, podczas gdy katalogi zadań w `work/` są zachowane.

### Podsumowanie

Wiesz, jak uruchomić prosty skrypt Nextflow, monitorować jego wykonanie i znajdować jego wyjścia.

### Co dalej?

Naucz się czytać podstawowy skrypt Nextflow i identyfikować, jak jego komponenty odnoszą się do jego funkcjonalności.

---

## 3. Zbadaj skrypt startowy workflow Hello World

To, co tam zrobiliśmy, to w zasadzie traktowanie skryptu workflow jak czarnej skrzynki.
Teraz, gdy widzieliśmy, co robi, otwórzmy skrzynkę i zajrzyjmy do środka.

Naszym celem tutaj nie jest zapamiętywanie składni kodu Nextflow, ale wyrobienie sobie podstawowej intuicji na temat tego, jakie są główne komponenty i jak są zorganizowane.

### 3.1. Zbadaj ogólną strukturę kodu

Znajdziesz skrypt `1-hello.nf` w bieżącym katalogu, którym powinien być `nextflow-run`. Otwórz go w panelu edytora.

??? full-code "Pełny plik kodu"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo do wypisania 'Hello World!' do pliku
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
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // wyemituj pozdrowienie
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

Skrypt workflow Nextflow zazwyczaj zawiera jedną lub więcej definicji **process**, sam **workflow** i kilka opcjonalnych bloków, takich jak **params** i **output**.

Każdy **process** opisuje, jakie operacje powinien wykonać odpowiedni krok w pipeline, podczas gdy **workflow** opisuje logikę przepływu danych, która łączy różne kroki.

Przyjrzyjmy się bliżej najpierw blokowi **process**, a potem spojrzymy na blok **workflow**.

### 3.2. Definicja `process`

Pierwszy blok kodu opisuje **process**.
Definicja process zaczyna się od słowa kluczowego `process`, po którym następuje nazwa process, a na końcu ciało process ograniczone nawiasami klamrowymi.
Ciało process musi zawierać blok script, który określa polecenie do uruchomienia, które może być czymkolwiek, co można by uruchomić w terminalu wiersza poleceń.

```groovy title="1-hello.nf" linenums="3"
/*
* Użyj echo do wypisania pozdrowienia do pliku
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
--8<-- "docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

To jest bardzo minimalna definicja process, która zawiera tylko definicję `input`, definicję `output` i `script` do wykonania.

Definicja `input` zawiera kwalifikator `val`, który mówi Nextflow, że oczekuje wartości jakiegoś rodzaju (może to być ciąg znaków, liczba, cokolwiek).

Definicja `output` zawiera kwalifikator `path`, który mówi Nextflow, że powinien to być traktowany jako ścieżka (obejmuje zarówno ścieżki katalogów, jak i pliki).

### 3.3. Definicja `workflow`

Drugi blok kodu opisuje sam **workflow**.
Definicja workflow zaczyna się od słowa kluczowego `workflow`, po którym następuje opcjonalna nazwa, a następnie ciało workflow ograniczone nawiasami klamrowymi.

Tutaj mamy **workflow**, który składa się z bloku `main:` i bloku `publish:`.
Blok `main:` jest głównym ciałem workflow, a blok `publish:` wymienia wyjścia, które powinny być opublikowane do katalogu `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // wyemituj pozdrowienie
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

W tym przypadku blok `main:` zawiera wywołanie process `sayHello` i przekazuje mu dane wejściowe o nazwie `params.input` do użycia jako powitanie.

Jak omówimy bardziej szczegółowo za chwilę, `params.input` przechowuje wartość, którą podaliśmy parametrowi `--input` w naszym wierszu poleceń.

Blok `publish:` wymienia wyjście wywołania process `sayHello()`, które odnosi się jako `sayHello.out` i nadaje mu nazwę `first_output` (może to być cokolwiek, co chce autor workflow).

To jest bardzo minimalna definicja **workflow**.
W rzeczywistym pipeline workflow zazwyczaj zawiera wiele wywołań **process** połączonych przez **channel** i mogą być ustawione wartości domyślne dla zmiennych wejściowych.

Zajmiemy się tym w Części 2 kursu.
Na razie przyjrzyjmy się bliżej temu, jak nasz workflow obsługuje dane wejściowe i wyjściowe.

### 3.4. System parametrów wiersza poleceń `params`

`params.input`, który przekazujemy do wywołania process `sayHello()`, to zgrabny fragment kodu Nextflow i warto poświęcić mu dodatkową minutę.

Jak wspomniano powyżej, w ten sposób przekazujemy wartość parametru wiersza poleceń `--input` do wywołania process `sayHello()`.
W rzeczywistości samo zadeklarowanie `params.someParameterName` wystarczy, aby dać workflow parametr o nazwie `--someParameterName` z wiersza poleceń.

Tutaj sformalizowaliśmy tę deklarację parametru, konfigurując blok `params`, który określa typ danych wejściowych, jakich oczekuje workflow (Nextflow 25.10.2 i nowsze).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

Obsługiwane typy obejmują `String`, `Integer`, `Float`, `Boolean` i `Path`.

!!! tip "Wskazówka"

    Parametry workflow zadeklarowane przy użyciu systemu `params` zawsze przyjmują dwa myślniki w wierszu poleceń (`--`).
    To odróżnia je od parametrów poziomu Nextflow, które przyjmują tylko jeden myślnik (`-`).

### 3.5. Dyrektywa `publish`

Po drugiej stronie workflow, już rzuciliśmy okiem na blok `publish:`.
To jedna połowa systemu obsługi wyjść; druga połowa to blok `output` znajdujący się poniżej.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

To określa, że wyjście `first_output` wymienione w bloku `publish:` powinno być skopiowane do podkatalogu o nazwie `1-hello` w domyślnym katalogu wyjściowym `results`.

Linia `mode 'copy'` nadpisuje domyślne zachowanie systemu, które polega na tworzeniu dowiązania symbolicznego (lub symlinku) do oryginalnego pliku w katalogu `work/` zamiast właściwej kopii.

Jest więcej opcji niż pokazano tutaj do kontrolowania zachowania publikowania; omówimy kilka później.
Zobaczysz również, że gdy workflow generuje wiele wyjść, każde z nich jest wymieniane w ten sposób w bloku `output`.

??? info "Starsza składnia publikowania wyjść przy użyciu `publishDir`"

    Aż do niedawna ustaloną metodą publikowania wyjść było robienie tego na poziomie każdego indywidualnego process przy użyciu dyrektywy `publishDir`.

    Nadal znajdziesz ten wzorzec kodu wszędzie w starszych pipeline Nextflow i modułach process, więc ważne jest, aby o nim wiedzieć.

    Zamiast mieć blok `publish:` w workflow i blok `output` na najwyższym poziomie, zobaczyłbyś linię `publishDir` w definicji process `sayHello`:

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

    Jednak nie zalecamy używania tego w żadnej nowej pracy, ponieważ zostanie to ostatecznie zabronione w przyszłych wersjach języka Nextflow.

### Podsumowanie

Teraz wiesz, jak prosty workflow Nextflow jest zbudowany i jak podstawowe komponenty odnoszą się do jego funkcjonalności.

### Co dalej?

Naucz się wygodnie zarządzać wykonaniami workflow.

---

## 4. Zarządzaj wykonaniami workflow

Wiedza o tym, jak uruchamiać workflow i pobierać wyjścia, jest świetna, ale szybko odkryjesz, że jest kilka innych aspektów zarządzania workflow, które ułatwią ci życie.

Tutaj pokażemy, jak wykorzystać funkcję `resume`, gdy musisz ponownie uruchomić ten sam workflow, jak sprawdzić logi wykonania za pomocą `nextflow log` i jak usunąć starsze katalogi robocze za pomocą `nextflow clean`.

### 4.1. Ponownie uruchom workflow z `-resume`

Czasami będziesz chciał ponownie uruchomić pipeline, który już wcześniej uruchomiłeś, bez powtarzania pracy, która została już pomyślnie ukończona.

Nextflow ma opcję o nazwie `-resume`, która pozwala to zrobić.
Konkretnie, w tym trybie wszelkie procesy, które zostały już uruchomione z dokładnie tym samym kodem, ustawieniami i danymi wejściowymi, zostaną pominięte.
To oznacza, że Nextflow uruchomi tylko procesy, które dodałeś lub zmodyfikowałeś od ostatniego uruchomienia, lub do których przekazujesz nowe ustawienia lub dane wejściowe.

Są dwie kluczowe zalety robienia tego:

- Jeśli jesteś w trakcie rozwijania pipeline, możesz iterować szybciej, ponieważ musisz uruchomić tylko process(y), nad którymi aktywnie pracujesz, aby przetestować swoje zmiany.
- Jeśli uruchamiasz pipeline w produkcji i coś pójdzie nie tak, w wielu przypadkach możesz naprawić problem i ponownie uruchomić pipeline, a on wznowi działanie od punktu awarii, co może zaoszczędzić dużo czasu i zasobów obliczeniowych.

Aby go użyć, po prostu dodaj `-resume` do polecenia i uruchom je:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

Wyjście konsoli powinno wyglądać znajomo, ale jest jedna rzecz, która jest trochę inna niż wcześniej.

Poszukaj fragmentu `cached:`, który został dodany w linii statusu process (linia 5), co oznacza, że Nextflow rozpoznał, że już wykonał tę pracę i po prostu ponownie użył wyniku z poprzedniego pomyślnego uruchomienia.

Możesz również zobaczyć, że hash podkatalogu roboczego jest taki sam jak w poprzednim uruchomieniu.
Nextflow dosłownie wskazuje ci poprzednie wykonanie i mówi "Już to zrobiłem tam".

!!! tip "Wskazówka"

    Gdy ponownie uruchamiasz pipeline z `resume`, Nextflow nie nadpisuje żadnych plików opublikowanych poza katalogiem roboczym przez jakiekolwiek wykonania, które zostały wcześniej pomyślnie uruchomione.

### 4.2. Sprawdź dziennik poprzednich wykonań

Za każdym razem, gdy uruchamiasz workflow nextflow, linia jest zapisywana do pliku dziennika o nazwie `history`, w ukrytym katalogu o nazwie `.nextflow` w bieżącym katalogu roboczym.

??? abstract "Zawartość pliku"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Ten plik zawiera znacznik czasu, nazwę uruchomienia, status, ID rewizji, ID sesji i pełne polecenie wiersza poleceń dla każdego uruchomienia Nextflow, które zostało uruchomione z bieżącego katalogu roboczego.

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

To wypisze zawartość pliku dziennika do terminala, wzbogaconą o linię nagłówka.

Zauważysz, że ID sesji zmienia się za każdym razem, gdy uruchamiasz nowe polecenie `nextflow run`, Z WYJĄTKIEM gdy używasz opcji `-resume`.
W takim przypadku ID sesji pozostaje takie samo.

Nextflow używa ID sesji do grupowania informacji o pamięci podręcznej uruchomienia w katalogu `cache`, również znajdującym się w `.nextflow`.

### 4.3. Usuń starsze katalogi robocze

Jeśli uruchamiasz dużo pipeline, możesz zgromadzić bardzo wiele plików w wielu podkatalogach.
Ponieważ podkatalogi są nazywane losowo, trudno jest powiedzieć po ich nazwach, które są starsze, a które nowsze.

Na szczęście Nextflow zawiera pomocne podpolecenie `clean`, które może automatycznie usuwać podkatalogi robocze dla poprzednich uruchomień, na których ci już nie zależy.

#### 4.3.1. Określ kryteria usuwania

Istnieje wiele [opcji](https://www.nextflow.io/docs/latest/reference/cli.html#clean) do określenia, co usunąć.

Tutaj pokazujemy przykład, który usuwa wszystkie podkatalogi z uruchomień przed danym uruchomieniem, określonym przy użyciu jego nazwy uruchomienia.

Znajdź najnowsze pomyślne uruchomienie, w którym nie użyłeś `-resume`; w naszym przypadku nazwa uruchomienia to `backstabbing_swartz`.

Nazwa uruchomienia to automatycznie generowany dwuczęściowy ciąg znaków pokazany w nawiasach kwadratowych w linii wyjścia konsoli `Launching (...)`.
Możesz również użyć dziennika Nextflow, aby wyszukać uruchomienie na podstawie jego znacznika czasu i/lub wiersza poleceń.

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

Twoje wyjście będzie miało inne nazwy katalogów zadań i może mieć inną liczbę linii, ale powinno wyglądać podobnie do przykładu.

Jeśli nie widzisz żadnych wypisanych linii, albo nie podałeś prawidłowej nazwy uruchomienia, albo nie ma poprzednich uruchomień do usunięcia. Upewnij się, że zmienisz `backstabbing_swartz` w przykładowym poleceniu na odpowiednią najnowszą nazwę uruchomienia w twoim dzienniku.

#### 4.3.3. Kontynuuj usuwanie

Jeśli wyjście wygląda zgodnie z oczekiwaniami i chcesz kontynuować usuwanie, uruchom ponownie polecenie z flagą `-f` zamiast `-n`:

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

!!! warning "Ostrzeżenie"

    Usuwanie podkatalogów roboczych z poprzednich uruchomień usuwa je z pamięci podręcznej Nextflow i kasuje wszelkie wyjścia, które były przechowywane w tych katalogach.
    To oznacza, że psuje zdolność Nextflow do wznowienia wykonania bez ponownego uruchamiania odpowiednich procesów.

    Jesteś odpowiedzialny za zapisanie wszelkich wyjść, na których ci zależy! To jest główny powód, dla którego wolimy używać trybu `copy` zamiast trybu `symlink` dla dyrektywy `publish`.

### Podsumowanie

Wiesz, jak ponownie uruchomić pipeline bez powtarzania kroków, które były już uruchomione w identyczny sposób, sprawdzać dziennik wykonania i używać polecenia `nextflow clean` do czyszczenia starych katalogów roboczych.

### Co dalej?

Zrób sobie małą przerwę! Właśnie przyswoiłeś budulce składni Nextflow i podstawowe instrukcje użycia.

W następnej sekcji tego szkolenia przyjrzymy się czterem kolejnym, coraz bardziej realistycznym wersjom pipeline Hello World, które pokażą, jak Nextflow pozwala efektywnie przetwarzać wiele danych wejściowych, uruchamiać workflow złożone z wielu połączonych ze sobą kroków, wykorzystywać modułowe komponenty kodu i używać kontenerów dla większej odtwarzalności i przenośności.

---

## Quiz

<quiz>
W linii wyjścia konsoli `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, co reprezentuje `[a3/7be2fa]`?
- [ ] Numer wersji process
- [ ] Unikalny identyfikator uruchomienia
- [x] Skrócona ścieżka do katalogu roboczego zadania
- [ ] Suma kontrolna pliku wyjściowego

Dowiedz się więcej: [2.3. Znajdź oryginalne wyjście i logi w katalogu `work/`](#23-znajdz-oryginalne-wyjscie-i-logi-w-katalogu-work)
</quiz>

<quiz>
Jaki jest cel pliku `.command.sh` w katalogu zadania?
- [ ] Przechowuje ustawienia konfiguracyjne zadania
- [x] Pokazuje rzeczywiste polecenie, które zostało wykonane przez process
- [ ] Zawiera komunikaty o błędach z nieudanych zadań
- [ ] Wymienia pliki wejściowe przygotowane dla zadania

Dowiedz się więcej: [2.3. Znajdź oryginalne wyjście i logi w katalogu `work/`](#23-znajdz-oryginalne-wyjscie-i-logi-w-katalogu-work)
</quiz>

<quiz>
Co dzieje się z opublikowanymi wynikami, gdy ponownie uruchamiasz workflow bez `-resume`?
- [ ] Są zachowane w oddzielnych katalogach ze znacznikami czasu
- [x] Są nadpisywane przez nowe wykonanie
- [ ] Nextflow zapobiega nadpisywaniu i kończy się niepowodzeniem
- [ ] Są automatycznie archiwizowane

Dowiedz się więcej: [2.4. Uruchom ponownie workflow z różnymi powitaniami](#24-uruchom-ponownie-workflow-z-roznymi-powitaniami)
</quiz>

<quiz>
Co wskazuje to wyjście konsoli?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Zadanie nie powiodło się i zostało pominięte
- [ ] Zadanie czeka w kolejce
- [x] Nextflow ponownie użył wyników z poprzedniego identycznego wykonania
- [ ] Zadanie zostało ręcznie anulowane

Dowiedz się więcej: [4.1. Ponownie uruchom workflow z `-resume`](#41-ponownie-uruchom-workflow-z--resume)
</quiz>

<quiz>
Gdzie Nextflow przechowuje historię wykonania, którą wyświetla polecenie `nextflow log`?
- [ ] W katalogu results
- [ ] W katalogu work
- [x] W pliku `.nextflow/history`
- [ ] W `nextflow.config`

Dowiedz się więcej: [4.2. Sprawdź dziennik poprzednich wykonań](#42-sprawdz-dziennik-poprzednich-wykonan)
</quiz>

<quiz>
Jaki jest cel bloku `params` w pliku workflow?
- [ ] Definiowanie wymagań zasobów process
- [ ] Konfigurowanie executora
- [x] Deklarowanie i typowanie parametrów wejściowych workflow
- [ ] Określanie opcji publikowania wyjść

Dowiedz się więcej: [3.4. System parametrów wiersza poleceń `params`](#34-system-parametrow-wiersza-polecen-params)
</quiz>

<quiz>
Co robi `mode 'copy'` w bloku `output` workflow?
- [ ] Tworzy kopię zapasową katalogu roboczego
- [x] Tworzy pełną kopię plików zamiast dowiązań symbolicznych
- [ ] Kopiuje skrypt workflow do results
- [ ] Włącza przyrostowe kopiowanie plików

Dowiedz się więcej: [3.5. Dyrektywa `publish`](#35-dyrektywa-publish)
</quiz>

<quiz>
Jaka jest zalecana flaga do użycia z poleceniem `nextflow clean` przed faktycznym usunięciem plików?
- [x] `-n` (próbne uruchomienie) aby zobaczyć, co zostałoby usunięte
- [ ] `-v` (szczegółowe) aby zobaczyć szczegółowe wyjście
- [ ] `-a` (wszystko) aby wybrać wszystkie katalogi
- [ ] `-q` (cicho) aby pominąć ostrzeżenia

Dowiedz się więcej: [4.3. Usuń starsze katalogi robocze](#43-usun-starsze-katalogi-robocze)
</quiz>
