# Część 3: Hello Workflow

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=pl" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Obejrzyj [całą playlistę](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/03_hello_workflow.md).
///

Większość rzeczywistych workflow składa się z więcej niż jednego kroku.
W tym module szkoleniowym nauczysz się łączyć procesy w wieloetapowy workflow.

Poznasz następujące techniki Nextflow'a:

1. Przepływ danych z jednego procesu do następnego
2. Zbieranie wyjść z wielu wywołań procesu do pojedynczego wywołania
3. Przekazywanie dodatkowych parametrów do procesu
4. Obsługę wielu wyjść pochodzących z procesu

Dla demonstracji będziemy kontynuować rozbudowę domenowo-agnostycznego przykładu Hello World z Części 1 i 2.
Tym razem wprowadzimy następujące zmiany w naszym workflow, aby lepiej odzwierciedlić sposób budowania rzeczywistych workflow:

1. Dodamy drugi krok konwertujący pozdrowienie na wielkie litery.
2. Dodamy trzeci krok zbierający wszystkie przekształcone pozdrowienia i zapisujący je do pojedynczego pliku.
3. Dodamy parametr do nazwania końcowego pliku wyjściowego i przekażemy go jako dodatkowe wejście do kroku zbierania.
4. Sprawimy, że krok zbierania będzie również raportował prostą statystykę o tym, co zostało przetworzone.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś Części 1-2 kursu [Hello Nextflow](./index.md), ale jeśli czujesz się komfortowo z podstawami omówionymi w tych sekcjach, możesz zacząć od tego miejsca bez dodatkowych przygotowań.

---

## 0. Rozgrzewka: Uruchom `hello-workflow.nf`

Użyjemy skryptu workflow'a `hello-workflow.nf` jako punktu wyjścia.
Jest on równoważny skryptowi utworzonemu podczas pracy nad Częścią 2 tego szkolenia, z tą różnicą, że usunęliśmy instrukcje `view()` i zmieniliśmy miejsce docelowe wyjścia:

```groovy title="hello-workflow.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
```

Ten diagram podsumowuje obecne działanie workflow'a.
Powinien wyglądać znajomo, z tą różnicą, że teraz wyraźnie pokazujemy, że wyjścia procesu są pakowane w kanał, podobnie jak wejścia.
Za chwilę wykorzystamy ten kanał wyjściowy.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-workflow-channels.svg"
</figure>

Aby upewnić się, że wszystko działa, uruchom skrypt raz przed wprowadzeniem zmian:

```bash
nextflow run hello-workflow.nf
```

??? success "Wynik polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [admiring_lamarr] DSL2 - revision: 4d4053520d

    executor >  local (3)
    [b1/5826b5] process > sayHello (2) [100%] 3 of 3 ✔
    ```

Jak poprzednio, pliki wyjściowe znajdziesz w lokalizacji określonej w bloku `output`.
W tym rozdziale jest to katalog `results/hello_workflow/`.

??? abstract "Zawartość katalogu"

    ```console
    results/hello_workflow
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    └── Holà-output.txt
    ```

Jeśli to zadziałało, jesteś gotowy do nauki składania wieloetapowego workflow'a.

---

## 1. Dodaj drugi krok do workflow'a

Dodamy krok konwertujący każde pozdrowienie na wielkie litery.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep.svg"
</figure>

W tym celu musimy wykonać trzy rzeczy:

- Zdefiniować polecenie, którego użyjemy do konwersji na wielkie litery.
- Napisać nowy proces opakowujący polecenie konwersji.
- Wywołać nowy proces w bloku workflow i skonfigurować go tak, aby przyjmował wyjście procesu `sayHello()` jako wejście.

### 1.1. Zdefiniuj polecenie konwersji na wielkie litery i przetestuj je w terminalu

Do konwersji pozdrowień na wielkie litery użyjemy klasycznego narzędzia UNIX o nazwie `tr` (od 'text replacement'), z następującą składnią:

```bash title="Syntax"
tr '[a-z]' '[A-Z]'
```

To bardzo naiwny one-liner do zamiany tekstu, który nie uwzględnia liter z akcentami, więc na przykład 'Holà' zostanie zamienione na 'HOLà', ale wystarczająco dobrze posłuży do demonstracji koncepcji Nextflow'a i to jest najważniejsze.

Aby to przetestować, możemy uruchomić polecenie `echo 'Hello World'` i przekierować jego wyjście do polecenia `tr`:

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Wyjściem jest plik tekstowy o nazwie `UPPER-output.txt`, który zawiera wersję ciągu `Hello World` napisaną wielkimi literami.

??? abstract "Zawartość pliku"

    ```console title="UPPER-output.txt"
    HELLO WORLD
    ```

To właśnie zamierzamy osiągnąć za pomocą naszego workflow'a.

### 1.2. Napisz krok konwersji jako proces Nextflow'a

Możemy wzorować nasz nowy proces na pierwszym, ponieważ chcemy użyć wszystkich tych samych komponentów.

Dodaj następującą definicję procesu do skryptu workflow, tuż pod pierwszym procesem:

```groovy title="hello-workflow.nf" linenums="20"
/*
 * Użyj narzędzia zamiany tekstu do przekształcenia pozdrowienia na wielkie litery
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
```

W tym procesie komponujemy nazwę drugiego pliku wyjściowego na podstawie nazwy pliku wejściowego, podobnie jak zrobiliśmy to oryginalnie dla wyjścia pierwszego procesu.

### 1.3. Dodaj wywołanie nowego procesu w bloku workflow

Teraz musimy powiedzieć Nextflow'owi, aby faktycznie wywołał proces, który właśnie zdefiniowaliśmy.

W bloku workflow wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="44" hl_lines="10-11"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)
        // convert the greeting to uppercase
        convertToUpper()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="44"
    workflow {

        main:
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emit a greeting
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

To jeszcze nie jest funkcjonalne, ponieważ nie określiliśmy, co powinno być wejściem do procesu `convertToUpper()`.

### 1.4. Przekaż wyjście pierwszego procesu do drugiego procesu

Teraz musimy sprawić, aby wyjście procesu `sayHello()` trafiało do procesu `convertToUpper()`.

Wygodnie jest to, że Nextflow automatycznie pakuje wyjście procesu do kanału, jak pokazano na diagramie w sekcji rozgrzewki.
Możemy odnieść się do kanału wyjściowego procesu jako `<process>.out`.

Więc wyjście procesu `sayHello` to kanał o nazwie `sayHello.out`, który możemy bezpośrednio podłączyć do wywołania `convertToUpper()`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-multistep-connector.svg"
</figure>

W bloku workflow wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="53" hl_lines="2"
        // convert the greeting to uppercase
        convertToUpper()
    ```

W prostym przypadku takim jak ten (jedno wyjście do jednego wejścia), to wszystko, co musimy zrobić, aby połączyć dwa procesy!

### 1.5. Skonfiguruj publikowanie wyjść workflow'a

Na koniec zaktualizujmy wyjścia workflow'a, aby publikować również wyniki z drugiego procesu.

#### 1.5.1. Zaktualizuj sekcję `publish:` bloku `workflow`

W bloku `workflow` wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="56" hl_lines="3"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="56"
        publish:
        first_output = sayHello.out
    }
    ```

Logika jest taka sama jak poprzednio.

#### 1.5.2. Zaktualizuj blok `output`

W bloku `output` wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="61" hl_lines="6-9"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="61"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Ponownie logika jest taka sama jak wcześniej.

To pokazuje, że możesz kontrolować ustawienia wyjścia na bardzo szczegółowym poziomie, dla każdego indywidualnego wyjścia.
Spróbuj zmienić ścieżki lub tryb publikacji dla jednego z procesów, aby zobaczyć, co się stanie.

Oczywiście oznacza to, że powtarzamy tutaj pewne informacje, co może stać się niewygodne, jeśli chcielibyśmy zaktualizować lokalizację dla wszystkich wyjść w ten sam sposób.
Później w kursie dowiesz się, jak konfigurować te ustawienia dla wielu wyjść w uporządkowany sposób.

### 1.6. Uruchom workflow z `-resume`

Przetestujmy to używając flagi `-resume`, ponieważ już pomyślnie uruchomiliśmy pierwszy krok workflow'a.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Wynik polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [high_cantor] DSL2 - revision: d746983511

    executor >  local (3)
    [ab/816321] process > sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [e0/ecf81b] process > convertToUpper (3) [100%] 3 of 3 ✔
    ```

W wyniku konsoli pojawia się teraz dodatkowa linia odpowiadająca nowemu procesowi, który właśnie dodaliśmy.

Wyjścia znajdziesz w katalogu `results/hello_workflow` zgodnie z ustawieniami w bloku `output`.

??? abstract "Zawartość katalogu"

    ```console
    results/hello_workflow/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

To wygodne! Ale warto też zajrzeć do katalogu roboczego jednego z wywołań drugiego procesu.

??? abstract "Zawartość katalogu"

    ```console
    work/e0/ecf81b4cacc648b9b994218d5b29d7/
    ├── Holà-output.txt -> /workspaces/training/hello-nextflow/work/ab/81632178cd37e9e815959278808819/Holà-output.txt
    └── UPPER-Holà-output.txt
    ```

Zauważ, że są tam dwa pliki `*-output`: wyjście pierwszego procesu oraz wyjście drugiego.

Wyjście pierwszego procesu znajduje się tam, ponieważ Nextflow **przygotował** je tam, aby mieć wszystko potrzebne do wykonania w tym samym podkatalogu.

Jednak w rzeczywistości jest to dowiązanie symboliczne wskazujące na oryginalny plik w podkatalogu wywołania pierwszego procesu.
Domyślnie, podczas uruchamiania na pojedynczej maszynie, jak to robimy tutaj, Nextflow używa dowiązań symbolicznych zamiast kopii do przygotowywania plików wejściowych i pośrednich.

Teraz, zanim przejdziemy dalej, pomyśl o tym, że wszystko, co zrobiliśmy, to połączenie wyjścia `sayHello` z wejściem `convertToUpper`, dzięki czemu dwa procesy mogły być uruchomione szeregowo.
Nextflow wykonał za nas ciężką pracę związaną z obsługą poszczególnych plików wejściowych i wyjściowych oraz przekazywaniem ich między dwoma poleceniami.

To jeden z powodów, dla których kanały w Nextflow są tak potężne: zajmują się rutynową pracą związaną z łączeniem kroków workflow'a.

### Podsumowanie

Wiesz już, jak łączyć procesy w łańcuch, przekazując wyjście jednego kroku jako wejście do następnego.

### Co dalej?

Dowiedz się, jak zbierać wyjścia z przetwarzanych partiami wywołań procesu i przekazywać je do pojedynczego procesu.

---

## 2. Dodaj trzeci krok do zbierania wszystkich pozdrowień

Gdy używamy procesu do zastosowania transformacji na każdym z elementów w kanale, jak tutaj z wieloma pozdrowieniami, czasami zachodzi potrzeba zebrania elementów z kanału wyjściowego tego procesu i przekazania ich do innego procesu wykonującego analizę lub podsumowanie.

Dla demonstracji dodamy nowy krok do naszego pipeline'u, który zbiera wszystkie pozdrowienia zapisane wielkimi literami, wyprodukowane przez proces `convertToUpper`, i zapisuje je do pojedynczego pliku.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

Nie zdradzając niespodzianki, ale będzie to wymagało użycia bardzo przydatnego operatora.

### 2.1. Zdefiniuj polecenie zbierania i przetestuj je w terminalu

Krok zbierania, który chcemy dodać do naszego workflow'a, użyje polecenia `cat` do połączenia wielu pozdrowień zapisanych wielkimi literami w pojedynczy plik.

Uruchommy polecenie samodzielnie w terminalu, aby sprawdzić, czy działa zgodnie z oczekiwaniami, tak jak robiliśmy to wcześniej.

Uruchom w terminalu:

```bash
echo 'Hello' | tr '[a-z]' '[A-Z]' > UPPER-Hello-output.txt
echo 'Bonjour' | tr '[a-z]' '[A-Z]' > UPPER-Bonjour-output.txt
echo 'Holà' | tr '[a-z]' '[A-Z]' > UPPER-Holà-output.txt
cat UPPER-Hello-output.txt UPPER-Bonjour-output.txt UPPER-Holà-output.txt > COLLECTED-output.txt
```

Wyjściem jest plik tekstowy o nazwie `COLLECTED-output.txt`, który zawiera wersje oryginalnych pozdrowień zapisane wielkimi literami.

??? abstract "Zawartość pliku"

    ```console title="COLLECTED-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

To jest wynik, który chcemy osiągnąć za pomocą naszego workflow'a.

### 2.2. Utwórz nowy proces do kroku zbierania

Utwórzmy nowy proces i nazwijmy go `collectGreetings()`.
Możemy zacząć go pisać na podstawie tego, co widzieliśmy wcześniej.

#### 2.2.1. Napisz 'oczywiste' części procesu

Dodaj następującą definicję procesu do skryptu workflow:

```groovy title="hello-workflow.nf" linenums="37"
/*
 * Zbierz pozdrowienia pisane wielkimi literami do jednego pliku wyjściowego
 */
process collectGreetings {

    input:
    ???

    output:
    path "COLLECTED-output.txt"

    script:
    """
    cat ??? > 'COLLECTED-output.txt'
    """
}
```

To jest to, co możemy napisać z pewnością na podstawie tego, czego się dotąd nauczyłeś.
Ale to nie jest funkcjonalne!
Brakuje definicji wejść i pierwszej połowy polecenia skryptu, ponieważ musimy wymyślić, jak to napisać.

#### 2.2.2. Zdefiniuj wejścia do `collectGreetings()`

Musimy zebrać pozdrowienia ze wszystkich wywołań procesu `convertToUpper()`.
Co wiemy, że możemy uzyskać z poprzedniego kroku workflow'a?

Kanał wyprodukowany przez `convertToUpper()` będzie zawierał ścieżki do poszczególnych plików zawierających pozdrowienia zapisane wielkimi literami.
To odpowiada jednemu slotowi wejściowemu; nazwijmy go `input_files` dla uproszczenia.

W bloku procesu wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          path input_files
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="2"
          input:
          ???
    ```

Zauważ, że używamy przedrostka `path`, mimo że oczekujemy, że będzie to zawierało wiele plików.

#### 2.2.3. Skomponuj polecenie konkatenacji

Tu sprawa może być nieco trudniejsza, ponieważ musimy być w stanie obsłużyć dowolną liczbę plików wejściowych.
Konkretnie, nie możemy napisać polecenia z góry, więc musimy powiedzieć Nextflow'owi, jak je skomponować w czasie wykonania na podstawie tego, jakie wejścia trafią do procesu.

Innymi słowy, jeśli mamy kanał wejściowy zawierający element `[file1.txt, file2.txt, file3.txt]`, potrzebujemy, aby Nextflow zamienił to na `cat file1.txt file2.txt file3.txt`.

Na szczęście Nextflow z przyjemnością zrobi to za nas, jeśli po prostu napiszemy `cat ${input_files}` w poleceniu skryptu.

W bloku procesu wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="54" hl_lines="3"
        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="54"
        script:
        """
        cat ??? > 'COLLECTED-output.txt'
        """
    ```

W teorii powinno to obsłużyć dowolną liczbę plików wejściowych.

!!! tip "Wskazówka"

    Niektóre narzędzia wiersza poleceń wymagają podania argumentu (np. `-input`) dla każdego pliku wejściowego.
    W takim przypadku musielibyśmy wykonać trochę dodatkowej pracy, aby skomponować polecenie.
    Przykład tego możesz zobaczyć w kursie [Nextflow dla genomiki](../../nf4_science/genomics/).

### 2.3. Dodaj krok zbierania do workflow'a

Teraz powinniśmy tylko wywołać proces zbierania na wyjściu kroku konwersji na wielkie litery.
To również jest kanał, o nazwie `convertToUpper.out`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-connector.svg"
</figure>

#### 2.3.1. Połącz wywołania procesów

W bloku workflow wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="75" hl_lines="4 5"
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out)
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="75"
        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

To łączy wyjście `convertToUpper()` z wejściem `collectGreetings()`.

#### 2.3.2. Uruchom workflow z `-resume`

Spróbujmy.

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Wynik polecenia"

    ```console hl_lines="8"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [mad_gilbert] DSL2 - revision: 6acfd5e28d

    executor >  local (3)
    [79/33b2f0] sayHello (2)         | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (3)   | 3 of 3, cached: 3 ✔
    [47/50fe4a] collectGreetings (1) | 3 of 3 ✔
    ```

Uruchomienie zakończyło się pomyślnie, włącznie z trzecim krokiem.

Jednak spójrz na liczbę wywołań `collectGreetings()` w ostatniej linii.
Oczekiwaliśmy tylko jednego, ale są trzy.

Teraz spójrz na zawartość końcowego pliku wyjściowego.

??? abstract "Zawartość pliku"

    ```console title="results/COLLECTED-output.txt"
    Holà
    ```

O nie. Krok zbierania został uruchomiony indywidualnie dla każdego pozdrowienia, co NIE jest tym, czego chcieliśmy.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-no-operator.svg"
</figure>

Musimy coś zrobić, aby wyraźnie powiedzieć Nextflow'owi, że chcemy, aby ten trzeci krok działał na wszystkich elementach w kanale wyprodukowanym przez `convertToUpper()`.

### 2.4. Użyj operatora do zebrania pozdrowień w pojedyncze wejście

Tak, ponownie odpowiedzią na nasz problem jest operator.

Konkretnie użyjemy trafnie nazwanego operatora [`collect()`](https://nextflow.io/docs/latest/reference/operator.html#collect).

#### 2.4.1. Dodaj operator `collect()`

Tym razem będzie to wyglądać nieco inaczej, ponieważ nie dodajemy operatora w kontekście fabryki kanałów; dodajemy go do kanału wyjściowego.

Bierzemy `convertToUpper.out` i dołączamy operator `collect()`, co daje nam `convertToUpper.out.collect()`.
Możemy to bezpośrednio podłączyć do wywołania procesu `collectGreetings()`.

W bloku workflow wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out)
    }
    ```

#### 2.4.2. Dodaj instrukcje `view()`

Dodajmy również kilka instrukcji `view()`, aby zwizualizować stany kanału przed i po zastosowaniu operatora.

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="73"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    }
    ```

Instrukcje `view()` mogą być umieszczone gdziekolwiek chcesz; umieściliśmy je zaraz po wywołaniu dla czytelności.

#### 2.4.3. Uruchom workflow ponownie z `-resume`

Spróbujmy:

```bash
nextflow run hello-workflow.nf -resume
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    Before collect: /workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt
    Before collect: /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt
    After collect: [/workspaces/training/hello-nextflow/work/b3/d52708edba8b864024589285cb3445/UPPER-Bonjour-output.txt, /workspaces/training/hello-nextflow/work/99/79394f549e3040dfc2440f69ede1fc/UPPER-Hello-output.txt, /workspaces/training/hello-nextflow/work/aa/56bfe7cf00239dc5badc1d04b60ac4/UPPER-Holà-output.txt]
    ```

Uruchomienie zakończyło się pomyślnie, chociaż wynik w logu może wyglądać nieco bardziej chaotycznie (uporządkowaliśmy go dla czytelności).

Tym razem trzeci krok został wywołany tylko raz!
Patrząc na wynik instrukcji `view()`, widzimy:

- Trzy instrukcje `Before collect:`, po jednej dla każdego pozdrowienia: w tym momencie ścieżki plików są indywidualnymi elementami w kanale.
- Pojedyncza instrukcja `After collect:`: trzy ścieżki plików są teraz spakowane w pojedynczy element.

Możemy to podsumować następującym diagramem:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-WITH-operator.svg"
</figure>

Na koniec możesz zajrzeć do zawartości pliku wyjściowego, aby upewnić się, że wszystko zadziałało poprawnie.

??? abstract "Zawartość pliku"

    ```console title="results/COLLECTED-output.txt"
    BONJOUR
    HELLO
    HOLà
    ```

Tym razem mamy wszystkie trzy pozdrowienia w końcowym pliku wyjściowym. Sukces!

!!! note "Uwaga"

    Jeśli uruchomisz to kilka razy bez `-resume`, zobaczysz, że kolejność pozdrowień zmienia się od uruchomienia do uruchomienia.
    To pokazuje, że kolejność, w jakiej elementy przepływają przez wywołania procesów, nie jest gwarantowana jako stała.

#### 2.4.4. Usuń instrukcje `view()` dla czytelności

Zanim przejdziesz do następnej sekcji, zalecamy usunięcie instrukcji `view()`, aby nie zaśmiecać wyjścia konsoli.

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="73"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="73" hl_lines="4-6"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())

        // optional view statements
        convertToUpper.out.view { contents -> "Before collect: $contents" }
        convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    ```

To jest zasadniczo odwrotna operacja od punktu 2.4.2.

### Podsumowanie

Wiesz już, jak zbierać wyjścia z partii wywołań procesu i przekazywać je do wspólnego kroku analizy lub podsumowania.

Podsumowując, oto co zbudowałeś do tej pory:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect.svg"
</figure>

### Co dalej?

Dowiedz się, jak przekazać więcej niż jedno wejście do procesu.

---

## 3. Przekaż dodatkowe parametry do procesu

Chcemy mieć możliwość nazwania końcowego pliku wyjściowego konkretną nazwą, aby móc przetwarzać kolejne partie pozdrowień bez nadpisywania poprzednich rezultatów.

W tym celu wprowadzimy następujące udoskonalenia do workflow'a:

- Zmodyfikujemy proces zbierający, aby akceptował zdefiniowaną przez użytkownika nazwę pliku wyjściowego (`batch_name`)
- Dodamy parametr wiersza poleceń do workflow'a (`--batch`) i przekażemy go do procesu zbierającego

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-batch.svg"
</figure>

### 3.1. Zmodyfikuj proces zbierający

Musimy zadeklarować dodatkowe wejście i zintegrować je z nazwą pliku wyjściowego.

#### 3.1.1. Zadeklaruj dodatkowe wejście

Dobra wiadomość: możemy zadeklarować tyle zmiennych wejściowych, ile chcemy w definicji procesu.
Nazwijmy tę zmienną `batch_name`.

W bloku procesu wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="42" hl_lines="3"
        input:
        path input_files
        val batch_name
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="42"
        input:
        path input_files
    ```

Możesz skonfigurować procesy tak, aby oczekiwały tylu wejść, ile chcesz.
W tej chwili wszystkie są ustawione jako wymagane wejścia; _musisz_ podać wartość, aby workflow działał.

Dowiesz się, jak zarządzać wymaganymi i opcjonalnymi wejściami później na swojej drodze z Nextflow.

#### 3.1.2. Użyj zmiennej `batch_name` w nazwie pliku wyjściowego

Możemy wstawić zmienną do nazwy pliku wyjściowego w ten sam sposób, w jaki komponowaliśmy dynamiczne nazwy plików wcześniej.

W bloku procesu wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-${batch_name}-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 6"
        output:
        path "COLLECTED-output.txt"

        script:
        """
        cat ${input_files} > 'COLLECTED-output.txt'
        """
    ```

To konfiguruje proces do używania wartości `batch_name` do generowania konkretnej nazwy pliku dla końcowego wyjścia workflow'a.

### 3.2. Dodaj parametr wiersza poleceń `batch`

Teraz potrzebujemy sposobu na dostarczenie wartości dla `batch_name` i przekazanie jej do wywołania procesu.

#### 3.2.1. Użyj `params` do skonfigurowania parametru

Wiesz już, jak używać systemu `params` do deklarowania parametrów CLI.
Użyjmy tego do zadeklarowania parametru `batch` (z domyślną wartością, bo jesteśmy leniwi).

W sekcji parametrów pipeline'u wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="6"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="55"
    /*
     * Pipeline parameters
     */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

Tak jak demonstrowaliśmy dla `--input`, możesz nadpisać tę domyślną wartość, określając wartość z `--batch` w wierszu poleceń.

#### 3.2.2. Przekaż parametr `batch` do procesu

Aby dostarczyć wartość parametru do procesu, musimy dodać ją w wywołaniu procesu.

W bloku workflow wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="74" hl_lines="2"
        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect())
    ```

Widzisz, że aby dostarczyć wiele wejść do procesu, po prostu wymieniasz je w nawiasach wywołania, oddzielone przecinkami.

!!! warning "Ostrzeżenie"

    MUSISZ dostarczyć wejścia do procesu w DOKŁADNIE TEJ SAMEJ KOLEJNOŚCI, w jakiej są wymienione w bloku definicji wejść procesu.

### 3.3. Uruchom workflow

Spróbujmy uruchomić to z nazwą partii w wierszu poleceń.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Wynik polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [confident_rutherford] DSL2 - revision: bc58af409c

    executor >  local (1)
    [79/33b2f0] sayHello (2)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [b5/f19efe] collectGreetings   | 1 of 1 ✔
    ```

Uruchomienie zakończyło się pomyślnie i wygenerowało oczekiwane wyjście:

??? abstract "Zawartość pliku"

    ```console title="results/COLLECTED-trio-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

Teraz, o ile odpowiednio określimy parametr, kolejne uruchomienia na innych partiach wejść nie nadpiszą poprzednich wyników.

### Podsumowanie

Wiesz już, jak przekazać więcej niż jedno wejście do procesu.

### Co dalej?

Dowiedz się, jak emitować wiele wyjść i wygodnie je obsługiwać.

---

## 4. Dodaj wyjście do kroku zbierającego

Do tej pory używaliśmy procesów, które produkowały tylko jedno wyjście każdy.
Mogliśmy uzyskać dostęp do ich wyjść bardzo wygodnie za pomocą składni `<process>.out`, której używaliśmy zarówno w kontekście przekazywania wyjścia do następnego procesu (np. `convertToUpper(sayHello.out)`), jak i w kontekście sekcji `publish:` (np. `first_output = sayHello.out`).

Co się dzieje, gdy proces produkuje więcej niż jedno wyjście?
Jak obsługujemy wiele wyjść?
Czy możemy wybrać i użyć konkretnego wyjścia?

Wszystko to doskonałe pytania, a krótka odpowiedź brzmi: tak, możemy!

Wiele wyjść zostanie spakowanych w oddzielne kanały.
Tym strumieniom wyjściowym można nadać nazwy, co ułatwia późniejsze odwoływanie się do nich indywidualnie, albo uzyskać do nich dostęp przez indeks.

Zagłębmy się w przykład.

Dla celów demonstracyjnych załóżmy, że chcemy policzyć liczbę pozdrowień zbieranych dla danej partii wejść i zaraportować to w pliku.

### 4.1. Zmodyfikuj proces, aby liczył i wyprowadzał liczbę pozdrowień

To będzie wymagało dwóch kluczowych zmian w definicji procesu: potrzebujemy sposobu na policzenie pozdrowień i zapisanie pliku raportu, a następnie musimy dodać ten plik raportu do bloku `output` procesu.

#### 4.1.1. Policz liczbę zebranych pozdrowień

Wygodnie jest to, że Nextflow pozwala nam dodawać dowolny kod w bloku `script:` definicji procesu, co jest bardzo przydatne do robienia takich rzeczy jak ta.

Oznacza to, że możemy użyć wbudowanej funkcji Nextflow'a `size()`, aby uzyskać liczbę plików w tablicy `input_files`, i zapisać wynik do pliku za pomocą polecenia `echo`.

W bloku procesu `collectGreetings` wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="55" hl_lines="2 5"
        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="55"
        script:
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        """
    ```

Zmienna `count_greetings` zostanie obliczona w czasie wykonania.

#### 4.1.2. Wyemituj plik raportu i nazwij wyjścia

W zasadzie wszystko, co musimy zrobić, to dodać plik raportu do bloku `output:`.

Jednak przy okazji dodamy również tagi `emit:` do naszych deklaracji wyjść. Pozwoli nam to wybierać wyjścia po nazwie zamiast używać indeksów pozycyjnych.

W bloku procesu wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="46" hl_lines="2 3"
        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="46"
        output:
        path "COLLECTED-${batch_name}-output.txt"
    ```

Tagi `emit:` są opcjonalne i mogliśmy dodać tag tylko do jednego z wyjść.
Ale czemu nie obu?

!!! tip "Wskazówka"

    Jeśli nie nazwiesz wyjść procesu za pomocą `emit:`, nadal możesz uzyskać do nich dostęp indywidualnie, używając ich odpowiedniego indeksu (liczonego od zera).
    Na przykład, użyłbyś `<process>.out[0]`, aby uzyskać pierwsze wyjście, `<process>.out[1]`, aby uzyskać drugie wyjście, i tak dalej.

    Preferujemy nazywanie wyjść, ponieważ w przeciwnym razie zbyt łatwo jest przez pomyłkę użyć złego indeksu, szczególnie gdy proces produkuje wiele wyjść.

### 4.2. Zaktualizuj wyjścia workflow'a

Teraz, gdy mamy dwa wyjścia z procesu `collectGreetings`, wyjście `collectGreetings.out` zawiera dwa kanały:

- `collectGreetings.out.outfile` zawiera końcowy plik wyjściowy
- `collectGreetings.out.report` zawiera plik raportu

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-report.svg"
</figure>

Musimy odpowiednio zaktualizować wyjścia workflow'a.

#### 4.2.1. Zaktualizuj sekcję `publish:`

W bloku `workflow` wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4 5"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="80" hl_lines="4"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out
    ```

Jak widać, odnoszenie się do konkretnych wyjść procesu jest teraz trywialne.
Gdy dodamy jeszcze jeden krok do naszego pipeline'u w Części 5 (Kontenery), będziemy mogli łatwo odnieść się do `collectGreetings.out.outfile` i przekazać go do nowego procesu (spoiler: nowy proces nazywa się `cowpy`).

Ale na razie dokończmy aktualizację wyjść na poziomie workflow'a.

#### 4.2.2. Zaktualizuj blok `output`

W bloku `output` wprowadź następującą zmianę w kodzie:

=== "Po"

    ```groovy title="hello-workflow.nf" linenums="86" hl_lines="14-17"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
        batch_report {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

=== "Przed"

    ```groovy title="hello-workflow.nf" linenums="80"
    output {
        first_output {
            path 'hello_workflow'
            mode 'copy'
        }
        uppercased {
            path 'hello_workflow'
            mode 'copy'
        }
        collected {
            path 'hello_workflow'
            mode 'copy'
        }
    }
    ```

Nie musimy aktualizować definicji wyjścia `collected`, ponieważ ta nazwa się nie zmieniła.
Musimy tylko dodać nowe wyjście.

### 4.3. Uruchom workflow

Spróbujmy uruchomić to z bieżącą partią pozdrowień.

```bash
nextflow run hello-workflow.nf -resume --batch trio
```

??? success "Wynik polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-workflow.nf` [ecstatic_wilson] DSL2 - revision: c80285f8c8

    executor >  local (1)
    [c5/4c6ca9] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [0e/6cbc59] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [02/61ead2] collectGreetings   [100%] 1 of 1 ✔
    ```

Jeśli zajrzysz do katalogu `results/hello_workflow/`, znajdziesz nowy plik raportu `trio-report.txt`.
Otwórz go, aby sprawdzić, czy workflow poprawnie zaraportował liczbę przetworzonych pozdrowień.

??? abstract "Zawartość pliku"

    ```txt title="trio-report.txt"
    There were 3 greetings in this batch.
    ```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-collect-4-way.svg"
</figure>

Możesz dodać więcej pozdrowień do pliku CSV i przetestować, co się stanie.

### Podsumowanie

Wiesz już, jak sprawić, aby proces emitował wiele nazwanych wyjść i jak odpowiednio obsługiwać je na poziomie workflow'a.

Ogólnie rzecz biorąc, rozumiesz kluczowe zasady związane z łączeniem procesów na typowe sposoby.

### Co dalej?

Weź dodatkową długą przerwę, zasłużyłeś na to.

Gdy będziesz gotowy, przejdź do [**Część 4: Hello Modules**](./04_hello_modules.md), aby dowiedzieć się, jak modularyzować Swój kod dla lepszej utrzymywalności i efektywności.

---

## Quiz

<quiz>
Jak uzyskujesz dostęp do wyjścia procesu w bloku workflow?
- [ ] `process.output`
- [ ] `output.processName`
- [x] `processName.out`
- [ ] `get(processName)`

Dowiedz się więcej: [1.4. Przekaż wyjście pierwszego procesu do drugiego procesu](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Co determinuje kolejność wykonywania procesów w Nextflow?
- [ ] Kolejność, w jakiej procesy są zapisane w bloku workflow
- [ ] Kolejność alfabetyczna według nazwy procesu
- [x] Zależności danych między procesami
- [ ] Losowa kolejność dla równoległego wykonania

Dowiedz się więcej: [1.4. Przekaż wyjście pierwszego procesu do drugiego procesu](#14-pass-the-output-of-the-first-process-to-the-second-process)
</quiz>

<quiz>
Który operator powinien zastąpić `???`, aby zebrać wszystkie wyjścia w pojedynczą listę dla procesu następczego?

```groovy hl_lines="4"
workflow {
    greetings_ch = Channel.of('Hello', 'Bonjour', 'Hola')
    SAYHELLO(greetings_ch)
    GATHER_ALL(SAYHELLO.out.???)
}
```

- [ ] `flatten()`
- [x] `collect()`
- [ ] `mix()`
- [ ] `join()`

Dowiedz się więcej: [2.4. Użyj operatora do zebrania pozdrowień w pojedyncze wejście](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
Kiedy należy używać operatora `collect()`?
- [ ] Gdy chcesz przetwarzać elementy równolegle
- [ ] Gdy potrzebujesz filtrować zawartość kanału
- [x] Gdy proces następczy potrzebuje wszystkich elementów z procesu poprzedzającego
- [ ] Gdy chcesz rozdzielić dane między wiele procesów

Dowiedz się więcej: [2.4. Użyj operatora do zebrania pozdrowień w pojedyncze wejście](#24-use-an-operator-to-collect-the-greetings-into-a-single-input)
</quiz>

<quiz>
Jak uzyskujesz dostęp do nazwanego wyjścia z procesu?
- [ ] `processName.outputName`
- [ ] `processName.get(outputName)`
- [x] `processName.out.outputName`
- [ ] `output.processName.outputName`

Dowiedz się więcej: [4.1.2. Wyemituj plik raportu i nazwij wyjścia](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Jaka jest poprawna składnia do nazywania wyjścia w procesie?
- [ ] `name: outputName`
- [ ] `output: outputName`
- [x] `emit: outputName`
- [ ] `label: outputName`

Dowiedz się więcej: [4.1.2. Wyemituj plik raportu i nazwij wyjścia](#412-emit-the-report-file-and-name-outputs)
</quiz>

<quiz>
Przy dostarczaniu wielu wejść do procesu, co musi być prawdą?
- [ ] Wszystkie wejścia muszą być tego samego typu
- [ ] Wejścia muszą być podane w kolejności alfabetycznej
- [x] Kolejność wejść musi odpowiadać kolejności zdefiniowanej w bloku input
- [ ] Tylko dwa wejścia mogą być podane jednocześnie

Dowiedz się więcej: [3. Przekaż dodatkowe parametry do procesu](#3-pass-more-than-one-input-to-a-process)
</quiz>
