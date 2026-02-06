# Część 2: Uruchamianie prawdziwych pipeline'ów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W Części 1 tego kursu (Uruchamianie podstawowych operacji) zaczęliśmy od przykładowego workflow'u, który miał tylko minimalne funkcje, aby utrzymać niską złożoność kodu.
Na przykład `1-hello.nf` używał parametru wiersza poleceń (`--input`) do przekazywania pojedynczej wartości na raz.

Jednak większość rzeczywistych pipeline'ów wykorzystuje bardziej zaawansowane funkcje, aby umożliwić efektywne przetwarzanie dużych ilości danych na dużą skalę i stosowanie wielu kroków przetwarzania połączonych czasami złożoną logiką.

W tej części szkolenia demonstrujemy kluczowe funkcje rzeczywistych pipeline'ów, wypróbowując rozszerzone wersje oryginalnego pipeline'u Hello World.

## 1. Przetwarzanie danych wejściowych z pliku

W rzeczywistym pipeline'ie zazwyczaj chcemy przetwarzać wiele punktów danych (lub serii danych) zawartych w jednym lub więcej plikach wejściowych.
I gdziekolwiek to możliwe, chcemy uruchamiać przetwarzanie niezależnych danych równolegle, aby skrócić czas oczekiwania na analizę.

Aby zademonstrować, jak Nextflow to robi, przygotowaliśmy plik CSV o nazwie `greetings.csv`, który zawiera kilka powitań wejściowych, naśladując rodzaj danych kolumnowych, które możesz chcieć przetwarzać w prawdziwej analizie danych.
Zauważ, że liczby nie mają znaczenia, są tam tylko w celach ilustracyjnych.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Napisaliśmy również ulepszoną wersję oryginalnego workflow'u, teraz o nazwie `2a-inputs.nf`, która odczyta plik CSV, wyodrębni powitania i zapisze każde z nich do oddzielnego pliku.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Uruchommy najpierw workflow, a potem przyjrzymy się odpowiedniemu kodowi Nextflow.

### 1.1. Uruchom workflow

Uruchom następujące polecenie w terminalu.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Ekscytująco, wydaje się to wskazywać, że wykonano '3 of 3' wywołań procesu, co jest zachęcające, ponieważ w dostarczonym pliku CSV były trzy wiersze danych.
To sugeruje, że proces `sayHello()` został wywołany trzy razy, raz dla każdego wiersza wejściowego.

### 1.2. Znajdź opublikowane wyjścia w katalogu `results`

Przyjrzyjmy się katalogowi 'results', aby zobaczyć, czy nasz workflow nadal zapisuje kopię naszych wyjść tam.

??? abstract "Zawartość katalogu"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Tak! Widzimy nowy katalog o nazwie `2a-inputs` z trzema plikami wyjściowymi o różnych nazwach, co jest wygodne.

Możesz otworzyć każdy z nich, aby upewnić się, że zawierają odpowiedni ciąg powitania.

??? abstract "Zawartość plików"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

To potwierdza, że każde powitanie w pliku wejściowym zostało odpowiednio przetworzone.

### 1.3. Znajdź oryginalne wyjścia i logi

Być może zauważyłeś, że powyższe wyjście konsoli odnosiło się tylko do jednego katalogu zadania.
Czy to oznacza, że wszystkie trzy wywołania `sayHello()` zostały wykonane w tym jednym katalogu zadania?

#### 1.3.1. Zbadaj katalog zadania podany w terminalu

Zajrzyjmy do tego katalogu zadania `8e/0eb066`.

??? abstract "Zawartość katalogu"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Znajdujemy tylko wyjście odpowiadające jednemu z powitań (oraz pliki pomocnicze, jeśli włączymy wyświetlanie ukrytych plików).

Więc co się dzieje?

Domyślnie system logowania ANSI zapisuje informacje o statusie dla wszystkich wywołań tego samego procesu w tej samej linii.
W rezultacie pokazał nam tylko jedną z trzech ścieżek katalogów zadań (`8e/0eb066`) w wyjściu konsoli.
Są dwie inne, które tam nie są wymienione.

#### 1.3.2. Spraw, aby terminal pokazywał więcej szczegółów

Możemy zmodyfikować zachowanie logowania, aby zobaczyć pełną listę wywołań procesu, dodając `-ansi-log false` do polecenia w następujący sposób:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Tym razem widzimy wszystkie trzy uruchomienia procesu i ich powiązane podkatalogi robocze wymienione w wyjściu.
Wyłączenie logowania ANSI również uniemożliwiło Nextflow używanie kolorów w wyjściu terminala.

Zauważ, że sposób raportowania statusu jest nieco inny między dwoma trybami logowania.
W trybie skondensowanym Nextflow raportuje, czy wywołania zostały pomyślnie ukończone, czy nie.
W tym rozszerzonym trybie raportuje tylko, że zostały przesłane.

To potwierdza, że proces `sayHello()` jest wywoływany trzy razy i dla każdego wywołania tworzony jest oddzielny katalog zadania.

Jeśli zajrzymy do każdego z katalogów zadań wymienionych tam, możemy zweryfikować, że każdy odpowiada jednemu z powitań.

??? abstract "Zawartość katalogu"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

To potwierdza, że każde wywołanie procesu jest wykonywane w izolacji od wszystkich innych.
Ma to wiele zalet, w tym unikanie kolizji, jeśli proces produkuje jakieś pliki pośrednie o nieunikatowych nazwach.

!!! tip "Wskazówka"

    Dla złożonego workflow'u lub dużej liczby danych wejściowych wyświetlanie pełnej listy do terminala może być nieco przytłaczające, więc ludzie normalnie nie używają `-ansi-log false` w rutynowym użyciu.

### 1.4. Zbadaj kod workflow'u

Więc ta wersja workflow'u jest w stanie odczytać plik CSV z danymi wejściowymi, przetwarzać dane wejściowe osobno i nazywać wyjścia unikalnie.

Przyjrzyjmy się, co to umożliwia w kodzie workflow'u.

??? full-code "Pełny plik kodu"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo do wypisania 'Hello World!' do pliku
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Ponownie, nie musisz zapamiętywać składni kodu, ale dobrze jest nauczyć się rozpoznawać kluczowe komponenty workflow'u, które zapewniają ważną funkcjonalność.

#### 1.4.1. Ładowanie danych wejściowych z CSV

To jest najciekawsza część: jak przeszliśmy od pobierania pojedynczej wartości z wiersza poleceń do wczytywania pliku CSV, parsowania go i obsługi zawartych w nim pojedynczych powitań?

W Nextflow robimy to za pomocą [**kanału**](https://nextflow.io/docs/latest/channel.html): konstrukcji kolejki zaprojektowanej do efektywnego zarządzania danymi wejściowymi i przekazywania ich z jednego kroku do drugiego w wieloetapowych workflow'ach, zapewniając jednocześnie wbudowaną paralelizację i wiele dodatkowych korzyści.

Rozłóżmy to na czynniki.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // utwórz kanał dla danych wejściowych z pliku CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // wyemituj powitanie
    sayHello(greeting_ch)
```

Ten kod tworzy kanał o nazwie `greeting_ch`, który odczytuje plik CSV, parsuje go i wyodrębnia pierwszą kolumnę z każdego wiersza.
Wynikiem jest kanał zawierający `Hello`, `Bonjour` i `Holà`.

??? tip "Jak to działa?"

    Oto co ta linia oznacza prostym językiem:

    - `channel.fromPath` to **fabryka kanału**, która tworzy kanał ze ścieżki(ek) pliku
    - `(params.input)` określa, że ścieżka pliku jest podana przez `--input` w wierszu poleceń

    Innymi słowy, ta linia mówi Nextflow: weź ścieżkę pliku podaną z `--input` i przygotuj się do traktowania jej zawartości jako danych wejściowych.

    Następnie kolejne dwie linie stosują **operatory**, które wykonują faktyczne parsowanie pliku i ładowanie danych do odpowiedniej struktury danych:

    - `.splitCsv()` mówi Nextflow, aby sparsował plik CSV na tablicę reprezentującą wiersze i kolumny
    - `.map { line -> line[0] }` mówi Nextflow, aby wziął tylko element z pierwszej kolumny z każdego wiersza

    Więc w praktyce, zaczynając od następującego pliku CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Przekształciliśmy to w tablicę, która wygląda tak:

    ```txt title="Zawartość tablicy"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    A następnie wzięliśmy pierwszy element z każdego z trzech wierszy i załadowaliśmy je do kanału Nextflow, który teraz zawiera: `Hello`, `Bonjour` i `Holà`.

    Jeśli chcesz zrozumieć kanały i operatory dogłębnie, w tym jak je samodzielnie pisać, zobacz [Hello Nextflow Część 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Wywołaj proces na każdym powitaniu

Następnie, w ostatniej linii bloku `main:` workflow'u, przekazujemy załadowany kanał `greeting_ch` jako dane wejściowe do procesu `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // utwórz kanał dla danych wejściowych z pliku CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // wyemituj powitanie
    sayHello(greeting_ch)
```

To mówi Nextflow, aby uruchomił proces indywidualnie na każdym elemencie w kanale, _tzn._ na każdym powitaniu.
A ponieważ Nextflow jest tak inteligentny, uruchomi te wywołania procesu równolegle, jeśli to możliwe, w zależności od dostępnej infrastruktury obliczeniowej.

To umożliwia efektywne i skalowalne przetwarzanie dużej ilości danych (wielu próbek lub punktów danych, cokolwiek jest Twoją jednostką badawczą) przy stosunkowo niewielkiej ilości kodu.

#### 1.4.3. Jak nazywane są wyjścia

Na koniec warto szybko spojrzeć na kod procesu, aby zobaczyć, jak uzyskujemy unikalne nazwy plików wyjściowych.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

Widzisz, że w porównaniu z wersją tego procesu w `1-hello.nf`, deklaracja wyjścia i odpowiednia część polecenia zmieniły się, aby uwzględnić wartość powitania w nazwie pliku wyjściowego.

To jest jeden ze sposobów na upewnienie się, że nazwy plików wyjściowych nie będą kolidować, gdy zostaną opublikowane do wspólnego katalogu results.

I to jest jedyna zmiana, którą musieliśmy wprowadzić wewnątrz deklaracji procesu!

### Podsumowanie

Rozumiesz na podstawowym poziomie, jak kanały i operatory umożliwiają nam efektywne przetwarzanie wielu danych wejściowych.

### Co dalej?

Odkryj, jak konstruowane są wieloetapowe workflow'y i jak działają.

---

## 2. Uruchamianie wieloetapowych workflow'ów

Większość rzeczywistych workflow'ów obejmuje więcej niż jeden krok.
Opierając się na tym, czego właśnie nauczyliśmy się o kanałach, przyjrzyjmy się, jak Nextflow używa kanałów i operatorów do łączenia procesów w wieloetapowym workflow'ie.

W tym celu dostarczamy przykładowy workflow, który łączy trzy oddzielne kroki i demonstruje:

1. Przepływ danych z jednego procesu do następnego
2. Zbieranie wyjść z wielu wywołań procesu do jednego wywołania procesu

Konkretnie, stworzyliśmy rozszerzoną wersję workflow'u o nazwie `2b-multistep.nf`, która przyjmuje każde powitanie wejściowe, konwertuje je na wielkie litery, a następnie zbiera wszystkie powitania wielkimi literami do jednego pliku wyjściowego.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Jak poprzednio, najpierw uruchomimy workflow, a potem przyjrzymy się kodowi, aby zobaczyć, co jest nowe.

### 2.1. Uruchom workflow

Uruchom następujące polecenie w terminalu:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Widzisz, że zgodnie z obietnicą, wiele kroków zostało uruchomionych jako część workflow'u; pierwsze dwa (`sayHello` i `convertToUpper`) były prawdopodobnie uruchomione na każdym indywidualnym powitaniu, a trzeci (`collectGreetings`) został uruchomiony tylko raz, na wyjściach wszystkich trzech wywołań `convertToUpper`.

### 2.2. Znajdź wyjścia

Zweryfikujmy, że to faktycznie się stało, patrząc na katalog `results`.

??? abstract "Zawartość katalogu"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Jak widzisz, mamy nowy katalog o nazwie `2b-multistep` i zawiera znacznie więcej plików niż wcześniej.
Niektóre pliki zostały zgrupowane w podkatalogu o nazwie `intermediates`, podczas gdy dwa pliki znajdują się na najwyższym poziomie.

Te dwa to końcowe wyniki wieloetapowego workflow'u.
Poświęć chwilę, aby spojrzeć na nazwy plików i sprawdzić ich zawartość, aby potwierdzić, że są takie, jakich oczekujesz.

??? abstract "Zawartość plików"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

Pierwszy zawiera nasze trzy powitania, wielkimi literami i zebrane z powrotem do jednego pliku, zgodnie z obietnicą.
Drugi to plik raportu, który podsumowuje pewne informacje o uruchomieniu.

### 2.3. Zbadaj kod

Przyjrzyjmy się kodowi i zidentyfikujmy kluczowe wzorce dla wieloetapowych workflow'ów.

??? full-code "Pełny plik kodu"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo do wypisania 'Hello World!' do pliku
    */
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

    /*
    * Użyj narzędzia zamiany tekstu do przekształcenia powitania na wielkie litery
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Zbierz powitania pisane wielkimi literami do jednego pliku wyjściowego
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // wyemituj powitanie
        sayHello(greeting_ch)
        // przekształć powitanie na wielkie litery
        convertToUpper(sayHello.out)
        // zbierz wszystkie powitania do jednego pliku
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Dużo się tam dzieje, ale najbardziej oczywistą różnicą w porównaniu z poprzednią wersją workflow'u jest to, że teraz jest wiele definicji procesów i odpowiednio kilka wywołań procesów w bloku workflow'u.

Przyjrzyjmy się bliżej i zobaczmy, czy możemy zidentyfikować najciekawsze elementy.

#### 2.3.1. Wizualizacja struktury workflow'u

Jeśli używasz VSCode z rozszerzeniem Nextflow, możesz uzyskać pomocny diagram pokazujący, jak procesy są połączone, klikając mały link `DAG preview` wyświetlany tuż nad blokiem workflow w dowolnym skrypcie Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

To daje ładny przegląd tego, jak procesy są połączone i co produkują.

Widzisz, że oprócz oryginalnego procesu `sayHello`, mamy teraz również `convertToUpper` i `collectGreetings`, które odpowiadają nazwom procesów, które widzieliśmy w wyjściu konsoli.
Dwie nowe definicje procesów są zbudowane w ten sam sposób co proces `sayHello`, z wyjątkiem tego, że `collectGreetings` przyjmuje dodatkowy parametr wejściowy o nazwie `batch` i produkuje dwa wyjścia.

Nie będziemy wchodzić w szczegóły kodu dla każdego z nich, ale jeśli jesteś ciekawy, możesz sprawdzić szczegóły w [Części 2 Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Na razie przyjrzyjmy się, jak procesy są ze sobą połączone.

#### 2.3.2. Jak procesy są połączone

Naprawdę interesującą rzeczą do przyjrzenia się jest to, jak wywołania procesów są połączone ze sobą w bloku `main:` workflow'u.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // utwórz kanał dla danych wejściowych z pliku CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // wyemituj powitanie
    sayHello(greeting_ch)
    // przekształć powitanie na wielkie litery
    convertToUpper(sayHello.out)
    // zbierz wszystkie powitania do jednego pliku
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Widzisz, że pierwsze wywołanie procesu, `sayHello(greeting_ch)`, jest niezmienione.
Następnie kolejne wywołanie procesu, do `convertToUpper`, odnosi się do wyjścia `sayHello` jako `sayHello.out`.

Wzorzec jest prosty: `processName.out` odnosi się do kanału wyjściowego procesu, który może być przekazany bezpośrednio do następnego procesu.
W ten sposób przekazujemy dane z jednego kroku do następnego w Nextflow.

#### 2.3.3. Proces może przyjmować wiele danych wejściowych

Trzecie wywołanie procesu, do `collectGreetings`, jest nieco inne.

```groovy title="2b-multistep.nf" linenums="77"
    // zbierz wszystkie powitania do jednego pliku
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Widzisz, że to wywołanie otrzymuje dwa dane wejściowe, `convertToUpper.out.collect()` i `params.batch`.
Ignorując na razie część `.collect()`, możemy to uogólnić jako `collectGreetings(input1, input2)`.

To odpowiada dwóm deklaracjom wejściowym w module procesu:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Gdy Nextflow parsuje to, przypisze pierwsze dane wejściowe w wywołaniu do `path input_files`, a drugie do `val batch_name`.

Więc teraz wiesz, że proces może przyjmować wiele danych wejściowych i jak wygląda wywołanie w bloku workflow'u.

Teraz przyjrzyjmy się bliżej temu pierwszemu wejściu, `convertToUpper.out.collect()`.

#### 2.3.4. Co robi `collect()` w wywołaniu `collectGreetings`

Aby przekazać wyjście `sayHello` do `convertToUpper`, po prostu odwołaliśmy się do kanału wyjściowego `sayHello` jako `sayHello.out`. Ale dla następnego kroku widzimy odniesienie do `convertToUpper.out.collect()`.

Co to jest ten fragment `collect()` i co robi?

To oczywiście operator. Tak jak operatory `splitCsv` i `map`, które napotkaliśmy wcześniej.
Tym razem operator nazywa się `collect` i jest stosowany do kanału wyjściowego produkowanego przez `convertToUpper`.

Operator `collect` służy do agregowania wyjść z wielu wywołań tego samego procesu i pakowania ich w pojedynczy element kanału.

W kontekście tego workflow'u pobiera trzy powitania wielkimi literami z `convertToUpper.out` --które są trzema oddzielnymi elementami i normalnie byłyby obsługiwane w oddzielnych wykonaniach przez następny proces-- i łączy je w jeden element.
To umożliwia zebranie wszystkich powitań z powrotem do tego samego pliku.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

W przeciwieństwie do tego, gdybyśmy nie zastosowali `collect()` do wyjścia `convertToUpper()` przed przekazaniem go do `collectGreetings()`, Nextflow po prostu uruchomiłby `collectGreetings()` niezależnie na każdym powitaniu, co nie osiągnęłoby naszego celu.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Jest wiele innych [operatorów](https://nextflow.io/docs/latest/reference/operator.html) dostępnych do stosowania transformacji na zawartości kanału między wywołaniami procesów.

To daje twórcom pipeline'ów dużo elastyczności w dostosowywaniu logiki przepływu ich pipeline'u.
Wadą jest to, że czasami może to utrudnić rozszyfrowanie tego, co pipeline robi.

#### 2.3.5. Parametr wejściowy może mieć wartość domyślną

Być może zauważyłeś, że `collectGreetings` przyjmuje drugie wejście, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // zbierz wszystkie powitania do jednego pliku
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

To przekazuje parametr CLI o nazwie `--batch` do workflow'u.
Jednak gdy uruchomiliśmy workflow wcześniej, nie określiliśmy parametru `--batch`.

Co się dzieje?
Przyjrzyj się blokowi `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

W workflow'ie jest skonfigurowana wartość domyślna, więc nie musimy jej podawać.
Ale jeśli podamy jedną w wierszu poleceń, wartość, którą określimy, zostanie użyta zamiast domyślnej.

Spróbuj:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Wyjście polecenia"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Powinieneś zobaczyć nowe końcowe wyjścia nazwane Twoją własną nazwą batch.

??? abstract "Zawartość katalogu"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

To jest aspekt konfiguracji danych wejściowych, który omówimy bardziej szczegółowo w Części 3, ale na razie ważne jest, aby wiedzieć, że parametry wejściowe mogą mieć wartości domyślne.

#### 2.3.6. Proces może produkować wiele wyjść

W definicji procesu `collectGreetings` widzimy następujące deklaracje wyjść:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Które są następnie przywoływane przez nazwę podaną z `emit:` w bloku `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

To ułatwia przekazywanie konkretnych wyjść indywidualnie do innych procesów w workflow'ie, w połączeniu z różnymi operatorami.

#### 2.3.7. Opublikowane wyjścia można organizować

W bloku `output` użyliśmy własnych ścieżek do grupowania wyników pośrednich, aby ułatwić wybranie tylko końcowych wyjść workflow'u.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Istnieją bardziej zaawansowane sposoby organizowania opublikowanych wyjść; omówimy kilka w części o konfiguracji.

!!! tip "Chcesz dowiedzieć się więcej o budowaniu workflow'ów?"

    Szczegółowe omówienie budowania wieloetapowych workflow'ów znajdziesz w [Hello Nextflow Część 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Podsumowanie

Rozumiesz na podstawowym poziomie, jak wieloetapowe workflow'y są konstruowane przy użyciu kanałów i operatorów oraz jak działają.
Widziałeś również, że procesy mogą przyjmować wiele danych wejściowych i produkować wiele wyjść, oraz że mogą być publikowane w uporządkowany sposób.

### Co dalej?

Dowiedz się, jak pipeline'y Nextflow mogą być modułowe, aby promować ponowne wykorzystanie kodu i łatwość konserwacji.

---

## 3. Uruchamianie zmodularyzowanych pipeline'ów

Jak dotąd wszystkie workflow'y, które oglądaliśmy, składały się z jednego pojedynczego pliku workflow'u zawierającego cały odpowiedni kod.

Jednak rzeczywiste pipeline'y zazwyczaj korzystają z _modularyzacji_, co oznacza, że kod jest podzielony na różne pliki.
To może uczynić ich rozwój i konserwację bardziej efektywnymi i zrównoważonymi.

Tutaj zademonstrujemy najpopularniejszą formę modułowości w Nextflow, czyli użycie **modułów**.

W Nextflow [**moduł**](https://nextflow.io/docs/latest/module.html) to pojedyncza definicja procesu zamknięta w samodzielnym pliku kodu.
Aby go użyć w workflow'ie, wystarczy dodać jednoliniową instrukcję importu do pliku kodu workflow'u; następnie można zintegrować proces w normalny sposób.
To umożliwia ponowne wykorzystanie definicji procesów w wielu workflow'ach bez tworzenia wielu kopii kodu.

Do tej pory uruchamialiśmy workflow'y, które miały wszystkie swoje procesy zawarte w monolitycznym pliku kodu.
Teraz zobaczymy, jak to wygląda, gdy procesy są przechowywane w indywidualnych modułach.

Oczywiście znowu przygotowaliśmy odpowiedni workflow do celów demonstracyjnych, o nazwie `2c-modules.nf`, wraz z zestawem modułów znajdujących się w katalogu `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Zawartość katalogu"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Widzisz, że są cztery pliki Nextflow, każdy nazwany po jednym z procesów.
Na razie możesz zignorować plik `cowpy.nf`; zajmiemy się nim później.

### 3.1. Zbadaj kod

Tym razem najpierw przyjrzymy się kodowi.
Zacznij od otwarcia pliku workflow `2c-modules.nf`.

??? full-code "Pełny plik kodu"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // wyemituj powitanie
        sayHello(greeting_ch)
        // przekształć powitanie na wielkie litery
        convertToUpper(sayHello.out)
        // zbierz wszystkie powitania do jednego pliku
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Widzisz, że logika workflow'u jest dokładnie taka sama jak w poprzedniej wersji workflow'u.
Jednak kod procesu zniknął z pliku workflow'u, a zamiast tego są instrukcje `include` wskazujące na oddzielne pliki w katalogu `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Dołącz moduły
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Otwórz jeden z tych plików, a znajdziesz kod dla odpowiedniego procesu.

??? full-code "Pełny plik kodu"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Użyj echo do wypisania 'Hello World!' do pliku
    */
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

Jak widzisz, kod procesu nie zmienił się; został po prostu skopiowany do indywidualnego pliku modułu zamiast być w głównym pliku workflow'u.
To samo dotyczy pozostałych dwóch procesów.

Zobaczmy więc, jak wygląda uruchomienie tej nowej wersji.

### 3.2. Uruchom workflow

Uruchom to polecenie w terminalu, z flagą `-resume`:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79484f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/4358gc] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Zauważysz, że wykonania procesu wszystkie pomyślnie użyły pamięci podręcznej, co oznacza, że Nextflow rozpoznał, że już wykonał żądaną pracę, mimo że kod został podzielony, a główny plik workflow'u został przemianowany.

Nic z tego nie ma znaczenia dla Nextflow; liczy się skrypt zadania, który jest generowany po zebraniu i ocenie całego kodu.

!!! tip "Wskazówka"

    Możliwe jest również zamknięcie sekcji workflow'u jako 'subworkflow', który można zaimportować do większego pipeline'u, ale to wykracza poza zakres tego kursu.

    Więcej o rozwijaniu komponowalnych workflow'ów możesz dowiedzieć się w Side Quest o [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Podsumowanie

Wiesz, jak procesy mogą być przechowywane w samodzielnych modułach, aby promować ponowne wykorzystanie kodu i poprawić łatwość konserwacji.

### Co dalej?

Naucz się używać kontenerów do zarządzania zależnościami oprogramowania.

---

## 4. Używanie konteneryzowanego oprogramowania

Jak dotąd workflow'y, których używaliśmy jako przykładów, musiały tylko uruchamiać bardzo podstawowe operacje przetwarzania tekstu przy użyciu narzędzi UNIX dostępnych w naszym środowisku.

Jednak rzeczywiste pipeline'y zazwyczaj wymagają specjalistycznych narzędzi i pakietów, które nie są domyślnie zawarte w większości środowisk.
Zazwyczaj musiałbyś zainstalować te narzędzia, zarządzać ich zależnościami i rozwiązywać wszelkie konflikty.

To wszystko jest bardzo żmudne i denerwujące.
Znacznie lepszym sposobem rozwiązania tego problemu jest użycie **kontenerów**.

**Kontener** to lekka, samodzielna, wykonywalna jednostka oprogramowania utworzona z **obrazu** kontenera, która zawiera wszystko, co potrzebne do uruchomienia aplikacji, w tym kod, biblioteki systemowe i ustawienia.

!!! tip "Wskazówka"

    Uczymy tego przy użyciu technologii [Docker](https://www.docker.com/get-started/), ale Nextflow obsługuje również kilka innych technologii kontenerowych.
    Więcej o wsparciu Nextflow dla kontenerów możesz dowiedzieć się [tutaj](https://nextflow.io/docs/latest/container.html).

### 4.1. Użyj kontenera bezpośrednio

Najpierw spróbujmy bezpośrednio wejść w interakcję z kontenerem.
To pomoże utrwalić zrozumienie tego, czym są kontenery, zanim zaczniemy ich używać w Nextflow.

#### 4.1.1. Pobierz obraz kontenera

Aby użyć kontenera, zazwyczaj pobierasz lub "ściągasz" obraz kontenera z rejestru kontenerów, a następnie uruchamiasz obraz kontenera, aby utworzyć instancję kontenera.

Ogólna składnia jest następująca:

```bash title="Składnia"
docker pull '<container>'
```

- `docker pull` to instrukcja dla systemu kontenerowego, aby pobrać obraz kontenera z repozytorium.
- `'<container>'` to adres URI obrazu kontenera.

Jako przykład pobierzmy obraz kontenera zawierający [cowpy](https://github.com/jeffbuttars/cowpy), pythonową implementację narzędzia o nazwie `cowsay`, które generuje grafikę ASCII do wyświetlania dowolnych tekstowych danych wejściowych w zabawny sposób.

Istnieją różne repozytoria, w których można znaleźć opublikowane kontenery.
Użyliśmy usługi [Seqera Containers](https://seqera.io/containers/), aby wygenerować ten obraz kontenera Docker z pakietu Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Uruchom pełne polecenie pull:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Wyjście polecenia"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

To mówi systemowi, aby pobrał określony obraz.
Po zakończeniu pobierania masz lokalną kopię obrazu kontenera.

#### 4.1.2. Uruchom kontener

Kontenery można uruchamiać jako jednorazowe polecenie, ale można ich również używać interaktywnie, co daje Ci wiersz poleceń wewnątrz kontenera i pozwala eksperymentować z poleceniem.

Ogólna składnia jest następująca:

```bash title="Składnia"
docker run --rm '<container>' [polecenie narzędzia]
```

- `docker run --rm '<container>'` to instrukcja dla systemu kontenerowego, aby uruchomić instancję kontenera z obrazu kontenera i wykonać w niej polecenie.
- `--rm` mówi systemowi, aby wyłączył instancję kontenera po zakończeniu polecenia.

W pełni złożone polecenie wykonania kontenera wygląda tak:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Uruchom to polecenie, a powinieneś zobaczyć, że Twój znak zachęty zmienia się na coś takiego jak `(base) root@b645838b3314:/tmp#`, co wskazuje, że jesteś teraz wewnątrz kontenera.

Możesz to zweryfikować, uruchamiając `ls`, aby wylistować zawartość katalogu:

```bash
ls /
```

??? success "Wyjście polecenia"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Widzisz, że system plików wewnątrz kontenera jest inny niż system plików na Twoim systemie hosta.

!!! tip "Wskazówka"

    Gdy uruchamiasz kontener, jest on domyślnie odizolowany od systemu hosta.
    To oznacza, że kontener nie może uzyskać dostępu do żadnych plików w systemie hosta, chyba że wyraźnie na to pozwolisz, określając, że chcesz zamontować wolumin jako część polecenia `docker run` używając następującej składni:

    ```bash title="Składnia"
    -v <ścieżka_zewnętrzna>:<ścieżka_wewnętrzna>
    ```

    To skutecznie ustanawia tunel przez ścianę kontenera, którego możesz użyć, aby uzyskać dostęp do tej części Twojego systemu plików.

    Jest to omówione bardziej szczegółowo w [Części 5 Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Uruchom narzędzie `cowpy`

Z wnętrza kontenera możesz uruchomić polecenie `cowpy` bezpośrednio.

```bash
cowpy "Hello Containers"
```

??? success "Wyjście polecenia"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

To produkuje grafikę ASCII domyślnej postaci krowy (lub 'cowacter') z dymkiem mowy zawierającym określony przez nas tekst.

Teraz, gdy przetestowałeś podstawowe użycie, możesz spróbować podać mu jakieś parametry.
Na przykład dokumentacja narzędzia mówi, że możemy ustawić postać za pomocą `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Wyjście polecenia"

    ```console
    __________________
    < Hello Containers >
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

Tym razem wyjście grafiki ASCII pokazuje pingwina Linuxa, Tux, ponieważ określiliśmy parametr `-c tux`.

Ponieważ jesteś wewnątrz kontenera, możesz uruchamiać polecenie cowpy tyle razy, ile chcesz, zmieniając parametry wejściowe, bez martwienia się o instalację jakichkolwiek bibliotek w samym systemie.

??? tip "Inne dostępne postacie"

    Użyj flagi '-c', aby wybrać inną postać, w tym:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Możesz się tym pobawić.
Gdy skończysz, wyjdź z kontenera używając polecenia `exit`:

```bash
exit
```

Znajdziesz się z powrotem w normalnej powłoce.

### 4.2. Użyj kontenera w workflow'ie

Gdy uruchamiamy pipeline, chcemy wskazać Nextflow, jakiego obrazu użyć w każdym kroku, a co ważne, oczekujemy też, że automatycznie obsłuży całą tę pracę: pobierze obraz, uruchomi instancję, wykona polecenie i usunie ją po zakończeniu.

Dobra wiadomość: to dokładnie to, co Nextflow zrobi za nas.
Musimy tylko określić kontener dla każdego procesu.

Aby zademonstrować, jak to działa, stworzyliśmy kolejną wersję naszego workflow'u, która uruchamia `cowpy` na pliku zebranych powitań wyprodukowanych w trzecim kroku.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

To powinno wyprodukować plik zawierający grafikę ASCII z trzema powitaniami w dymku mowy.

#### 4.2.1. Zbadaj kod

Workflow jest bardzo podobny do poprzedniego, plus dodatkowy krok do uruchomienia `cowpy`.

??? full-code "Pełny plik kodu"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // wyemituj powitanie
        sayHello(greeting_ch)
        // przekształć powitanie na wielkie litery
        convertToUpper(sayHello.out)
        // zbierz wszystkie powitania do jednego pliku
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // wygeneruj grafikę ASCII powitań za pomocą cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Widzisz, że ten workflow importuje proces `cowpy` z pliku modułu i wywołuje go na wyjściu wywołania `collectGreetings()`, plus parametr wejściowy o nazwie `params.character`.

```groovy title="2d-container.nf" linenums="31"
// wygeneruj grafikę ASCII powitań za pomocą cowpy
cowpy(collectGreetings.out.outfile, params.character)
```

Proces `cowpy`, który opakowuje polecenie cowpy do generowania grafiki ASCII, jest zdefiniowany w module `cowpy.nf`.

??? full-code "Pełny plik kodu"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Wygeneruj grafikę ASCII za pomocą cowpy (https://github.com/jeffbuttars/cowpy)
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

Proces `cowpy` wymaga dwóch danych wejściowych: ścieżki do pliku wejściowego zawierającego tekst do umieszczenia w dymku mowy (`input_file`) oraz wartości zmiennej character.

Co ważne, zawiera również linię `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, która wskazuje na URI kontenera, którego użyliśmy wcześniej.

#### 4.2.2. Sprawdź, czy Docker jest włączony w konfiguracji

Trochę wyprzedzimy Część 3 tego kursu szkoleniowego, wprowadzając plik konfiguracyjny `nextflow.config`, który jest jednym z głównych sposobów, jakie Nextflow oferuje do konfigurowania wykonywania workflow'u.
Gdy plik o nazwie `nextflow.config` jest obecny w bieżącym katalogu, Nextflow automatycznie go załaduje i zastosuje zawartą w nim konfigurację.

W tym celu dołączyliśmy plik `nextflow.config` z pojedynczą linią kodu, która włącza Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Ta konfiguracja mówi Nextflow, aby używał Docker dla każdego procesu, który określa kompatybilny kontener.

!!! tip "Wskazówka"

    Technicznie możliwe jest włączenie wykonywania Docker z wiersza poleceń, na podstawie pojedynczego uruchomienia, używając parametru `-with-docker <container>`.
    Jednak to pozwala nam określić tylko jeden kontener dla całego workflow'u, podczas gdy podejście, które właśnie pokazaliśmy, pozwala nam określić inny kontener dla każdego procesu.
    To drugie jest znacznie lepsze dla modułowości, konserwacji kodu i odtwarzalności.

#### 4.2.3. Uruchom workflow

Podsumowując, oto co zamierzamy uruchomić:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Myślisz, że zadziała?

Uruchommy workflow z flagą `-resume` i określmy, że chcemy, aby postacią był indyk.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

Pierwsze trzy kroki użyły pamięci podręcznej, ponieważ już je wcześniej uruchomiliśmy, ale proces `cowpy` jest nowy, więc faktycznie zostaje uruchomiony.

Możesz znaleźć wyjście kroku `cowpy` w katalogu `results`.

??? abstract "Zawartość pliku"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
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

Widzisz, że postać mówi wszystkie powitania, ponieważ uruchomiła się na pliku zebranych powitań wielkimi literami.

Co ważniejsze, mogliśmy to uruchomić jako część naszego pipeline'u bez konieczności prawidłowej instalacji cowpy i wszystkich jego zależności.
I teraz możemy udostępnić pipeline współpracownikom i sprawić, że uruchomią go na swojej infrastrukturze bez konieczności instalowania czegokolwiek, poza Docker lub jedną z jego alternatyw (taką jak Singularity/Apptainer) jak wspomniano powyżej.

#### 4.2.4. Sprawdź, jak Nextflow uruchomił konteneryzowane zadanie

Na zakończenie tej sekcji przyjrzyjmy się podkatalogowi roboczemu dla jednego z wywołań procesu `cowpy`, aby uzyskać nieco więcej wglądu w to, jak Nextflow pracuje z kontenerami pod maską.

Sprawdź wyjście z polecenia `nextflow run`, aby znaleźć ścieżkę do podkatalogu roboczego dla procesu `cowpy`.
Patrząc na to, co otrzymaliśmy dla uruchomienia pokazanego powyżej, linia dziennika konsoli dla procesu `cowpy` zaczyna się od `[7f/caf718]`.
To odpowiada następującej skróconej ścieżce katalogu: `work/7f/caf718`.

W tym katalogu znajdziesz plik `.command.run`, który zawiera wszystkie polecenia, które Nextflow uruchomił w Twoim imieniu w trakcie wykonywania pipeline'u.

??? abstract "Zawartość pliku"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Jeśli poszukasz `nxf_launch` w tym pliku, powinieneś zobaczyć coś takiego:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

To polecenie uruchomienia pokazuje, że Nextflow używa bardzo podobnego polecenia `docker run` do uruchomienia wywołania procesu, jak to zrobiliśmy, gdy uruchomiliśmy je ręcznie.
Montuje również odpowiedni podkatalog roboczy do kontenera, ustawia katalog roboczy wewnątrz kontenera odpowiednio i uruchamia nasz szablonowany skrypt bash w pliku `.command.sh`.

To potwierdza, że cała ciężka praca, którą musieliśmy wykonać ręcznie w poprzedniej sekcji, jest teraz wykonywana za nas przez Nextflow!

### Podsumowanie

Rozumiesz, jaką rolę odgrywają kontenery w zarządzaniu wersjami narzędzi oprogramowania i zapewnianiu odtwarzalności.

Ogólniej, masz podstawowe zrozumienie tego, czym są główne komponenty rzeczywistych pipeline'ów Nextflow i jak są zorganizowane.
Znasz podstawy tego, jak Nextflow może efektywnie przetwarzać wiele danych wejściowych, uruchamiać workflow'y składające się z wielu kroków połączonych ze sobą, wykorzystywać modułowe komponenty kodu i wykorzystywać kontenery dla większej odtwarzalności i przenośności.

### Co dalej?

Zrób kolejną przerwę! To była duża porcja informacji o tym, jak działają pipeline'y Nextflow.

W ostatniej sekcji tego szkolenia zagłębimy się głębiej w temat konfiguracji.
Dowiesz się, jak skonfigurować wykonanie Twojego pipeline'a, aby dopasować się do Twojej infrastruktury, a także zarządzać konfiguracją danych wejściowych i parametrów.

---

## Quiz

<quiz>
Dlaczego Nextflow tworzy oddzielny katalog zadania dla każdego wywołania procesu?
- [ ] Aby poprawić szybkość wykonywania
- [ ] Aby zmniejszyć użycie pamięci
- [x] Aby izolować wykonania i uniknąć kolizji między wyjściami
- [ ] Aby włączyć równoległą kompresję plików

Dowiedz się więcej: [1.3. Znajdź oryginalne wyjścia i logi](#13-znajdź-oryginalne-wyjścia-i-logi)
</quiz>

<quiz>
Co robi opcja `-ansi-log false` podczas uruchamiania workflow'u?
- [ ] Wyłącza całe wyjście konsoli
- [x] Usuwa kolor z wyjścia
- [x] Pokazuje wszystkie ścieżki katalogów zadań zamiast kondensowania ich w jednej linii
- [ ] Włącza tryb verbose debugowania

Dowiedz się więcej: [1.3.2. Spraw, aby terminal pokazywał więcej szczegółów](#132-spraw-aby-terminal-pokazywał-więcej-szczegółów)

Możesz również użyć jednej z następujących zmiennych środowiskowych, jeśli wolisz ten styl:

```bash
export NXF_ANSI_LOG=0
# lub
export NO_COLOR=1
```

</quiz>

<quiz>
W kodzie `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, co robi `#!groovy .map { line -> line[0] }`?
- [ ] Filtruje puste linie
- [ ] Sortuje linie alfabetycznie
- [x] Wyodrębnia pierwszą kolumnę z każdego wiersza CSV
- [ ] Liczy liczbę linii

Dowiedz się więcej: [1.4.1. Ładowanie danych wejściowych z CSV](#141-ładowanie-danych-wejściowych-z-csv)
</quiz>

<quiz>
Dlaczego ważne jest uwzględnienie wartości wejściowej w nazwach plików wyjściowych (np. `#!groovy "${greeting}-output.txt"`)?
- [ ] Aby poprawić szybkość przetwarzania
- [ ] Aby włączyć funkcjonalność resume
- [x] Aby zapobiec nadpisywaniu plików wyjściowych podczas przetwarzania wielu danych wejściowych
- [ ] Aby ułatwić kompresję plików

Dowiedz się więcej: [1.4.3. Jak nazywane są wyjścia](#143-jak-nazywane-są-wyjścia)
</quiz>

<quiz>
Jaki jest cel instrukcji `include` w zmodularizowanym workflow'ie?
- [ ] Aby skopiować kod procesu do pliku workflow'u
- [x] Aby zaimportować definicję procesu z zewnętrznego pliku modułu
- [ ] Aby dołączyć ustawienia konfiguracji
- [ ] Aby dodać komentarze dokumentacji

Dowiedz się więcej: [3. Uruchamianie zmodularyzowanych pipeline'ów](#3-uruchamianie-zmodularyzowanych-pipelineów)
</quiz>

<quiz>
Co się dzieje, gdy zmodularyzujesz workflow i uruchomisz go z `-resume`?
- [ ] Pamięć podręczna jest wyłączona dla procesów modułowych
- [ ] Wszystkie zadania muszą być ponownie wykonane
- [x] Pamięć podręczna działa normalnie w oparciu o wygenerowane skrypty zadań
- [ ] Tylko główny plik workflow'u jest buforowany

Dowiedz się więcej: [3.2. Uruchom workflow](#32-uruchom-workflow)
</quiz>

<quiz>
Co określa dyrektywa `container` w definicji procesu?
- [ ] Katalog roboczy procesu
- [ ] Maksymalną alokację pamięci
- [x] URI obrazu kontenera do użycia podczas uruchamiania procesu
- [ ] Format pliku wyjściowego

Dowiedz się więcej: [4.2. Użyj kontenera w workflow'ie](#42-użyj-kontenera-w-workflow)
</quiz>

<quiz>
Co zawiera funkcja `nxf_launch` w pliku `.command.run`?
- [ ] Informacje o wersji Nextflow
- [ ] Parametry workflow'u
- [x] Polecenie `docker run` z montowaniem woluminów i ustawieniami kontenera
- [ ] Deklaracje wejściowe procesu

Dowiedz się więcej: [4.2.4. Sprawdź, jak Nextflow uruchomił konteneryzowane zadanie](#424-sprawdź-jak-nextflow-uruchomił-konteneryzowane-zadanie)
</quiz>

<quiz>
Co Nextflow automatycznie obsługuje podczas uruchamiania konteneryzowanego procesu? (Zaznacz wszystkie, które mają zastosowanie)
- [x] Pobieranie obrazu kontenera, jeśli to konieczne
- [x] Montowanie katalogu roboczego do kontenera
- [x] Uruchamianie skryptu procesu wewnątrz kontenera
- [x] Czyszczenie instancji kontenera po wykonaniu

Dowiedz się więcej: [4. Używanie konteneryzowanego oprogramowania](#4-używanie-konteneryzowanego-oprogramowania)
</quiz>
