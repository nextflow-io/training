# Część 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/02_hello_channels.md).
///
-->

W Części 1 tego kursu (Hello World) pokazaliśmy Ci, jak dostarczyć zmienne wejście do procesu, podając je bezpośrednio w wywołaniu procesu: `sayHello(params.input)`.
To było celowo uproszczone podejście.
W praktyce takie podejście ma poważne ograniczenia; mianowicie działa tylko dla bardzo prostych sytuacji, gdy potrzebujemy uruchomić proces tylko raz, na pojedynczej wartości.
W większości realistycznych scenariuszy użycia workflow'u przetwarzamy wiele elementów (dane eksperymentalne dla wielu próbek, na przykład), więc potrzebujemy bardziej wyrafinowanego sposobu obsługi wejść.

Do tego służą **kanały** Nextflow.
Kanały to kolejki zaprojektowane do efektywnej obsługi wejść i przekazywania ich z jednego kroku do drugiego w wieloetapowych workflow'ach, zapewniając jednocześnie wbudowaną równoległość i wiele dodatkowych korzyści.

W tej części kursu nauczysz się, jak używać kanału do obsługi wielu wejść z różnych źródeł.
Nauczysz się również używać **operatorów** do transformowania zawartości kanału w razie potrzeby.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś Część 1 kursu [Hello Nextflow](./index.md), ale jeśli czujesz się komfortowo z podstawami omówionymi w tej sekcji, możesz zacząć od tego miejsca bez robienia czegokolwiek specjalnego.

---

## 0. Rozgrzewka: Uruchom `hello-channels.nf`

Użyjemy skryptu workflow'u `hello-channels.nf` jako punktu wyjścia.
Jest on równoważny ze skryptem powstałym w wyniku pracy przez Część 1 tego kursu szkoleniowego, z wyjątkiem tego, że zmieniliśmy miejsce docelowe wyjścia:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Aby upewnić się, że wszystko działa, uruchom skrypt raz przed wprowadzeniem jakichkolwiek zmian:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Jak poprzednio, plik wyjściowy o nazwie `output.txt` znajdziesz w katalogu `results/hello_channels` (jak określono w bloku `output` skryptu workflow'u, pokazanym powyżej).

??? abstract "Zawartość katalogu"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Zawartość pliku"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Jeśli to zadziałało, jesteś gotowy do nauki o kanałach.

---

## 1. Dostarczaj zmienne wejścia przez kanał jawnie

Zamierzamy utworzyć **kanał** do przekazania zmiennego wejścia do procesu `sayHello()` zamiast polegać na niejawnej obsłudze, która ma pewne ograniczenia.

### 1.1. Utwórz kanał wejściowy

Istnieje wiele **fabryk kanałów**, których możemy użyć do skonfigurowania kanału.
Aby na razie zachować prostotę, użyjemy najbardziej podstawowej fabryki kanałów o nazwie `channel.of`, która utworzy kanał zawierający pojedynczą wartość.
Funkcjonalnie będzie podobne do poprzedniej konfiguracji, ale zamiast polegać na niejawnym tworzeniu kanału przez Nextflow, robimy to teraz jawnie.

Oto linia kodu, której użyjemy:

```console title="Składnia"
greeting_ch = channel.of('Hello Channels!')
```

To tworzy kanał o nazwie `greeting_ch` używając fabryki kanałów `channel.of()`, która konfiguruje prosty kanał kolejki, i ładuje ciąg `'Hello Channels!'` jako wartość pozdrowienia.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Uwaga"

    Tymczasowo wracamy do zakodowanych na sztywno ciągów zamiast używać parametru CLI ze względu na czytelność. Wrócimy do używania parametrów CLI, gdy omówimy, co dzieje się na poziomie kanału.

W bloku workflow dodaj kod fabryki kanałów:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
        // wyemituj pozdrowienie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // wyemituj pozdrowienie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

To jeszcze nie jest funkcjonalne, ponieważ nie przełączyliśmy jeszcze wejścia do wywołania procesu.

### 1.2. Dodaj kanał jako wejście do wywołania procesu

Teraz musimy faktycznie podłączyć nasz nowo utworzony kanał do wywołania procesu `sayHello()`, zastępując parametr CLI, który podawaliśmy bezpośrednio wcześniej.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
        // wyemituj pozdrowienie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

To mówi Nextflow, aby uruchomił proces `sayHello` na zawartości kanału `greeting_ch`.

Teraz nasz workflow jest właściwie funkcjonalny; jest jawnym odpowiednikiem napisania `sayHello('Hello Channels!')`.

### 1.3. Uruchom workflow

Uruchommy to!

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Jeśli wprowadziłeś obie edycje poprawnie, powinieneś uzyskać pomyślne wykonanie.
Możesz sprawdzić katalog wyników, aby upewnić się, że wynik jest nadal taki sam jak poprzednio.

??? abstract "Zawartość pliku"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Zwiększyliśmy więc elastyczność naszego workflow'u, osiągając ten sam końcowy wynik.
Może się wydawać, że piszemy więcej kodu bez wymiernej korzyści, ale wartość stanie się jasna, gdy tylko zaczniemy obsługiwać więcej wejść.

Jako podgląd tego, spójrzmy na jeszcze jedną rzecz, zanim przejdziemy dalej: jedną małą, ale wygodną korzyść z używania jawnego kanału do zarządzania wejściem danych.

### 1.4. Użyj `view()` do inspekcji zawartości kanału

Kanały Nextflow są zbudowane w sposób, który pozwala nam operować na ich zawartości za pomocą operatorów, które omówimy szczegółowo później w tym rozdziale.

Na razie pokażemy Ci tylko, jak używać super prostego operatora o nazwie [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) do inspekcji zawartości kanału.
Możesz myśleć o `view()` jako o narzędziu do debugowania, jak instrukcja `print()` w Pythonie lub jej odpowiednik w innych językach.

Dodaj tę małą linię do bloku workflow:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Dokładna liczba spacji nie ma znaczenia, o ile jest wielokrotnością 4; dążymy tylko do wyrównania początku instrukcji `.view()` do części `.of()` konstrukcji kanału.

Teraz uruchom workflow ponownie:

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Jak widzisz, to wyświetla zawartość kanału do konsoli.
Tutaj mamy tylko jeden element, ale gdy zaczniemy ładować wiele wartości do kanału w następnej sekcji, zobaczysz, że jest ustawione na wyświetlanie jednego elementu na linię.

### Podsumowanie

Wiesz, jak używać podstawowej fabryki kanałów do dostarczenia wejścia do procesu.

### Co dalej?

Naucz się używać kanałów, aby workflow iterował po wielu wartościach wejściowych.

---

## 2. Zmodyfikuj workflow, aby działał na wielu wartościach wejściowych

Workflow'y zazwyczaj działają na partiach wejść, które mają być przetwarzane masowo, więc chcemy ulepszyć workflow, aby akceptował wiele wartości wejściowych.

### 2.1. Załaduj wiele pozdrowień do kanału wejściowego

Dogodnie, fabryka kanałów `channel.of()`, której używamy, chętnie przyjmuje więcej niż jedną wartość, więc nie musimy jej w ogóle modyfikować.
Możemy po prostu załadować wiele wartości do kanału.

Niech to będą `'Hello'`, `'Bonjour'` i `'Holà'`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Na diagramie kanał jest reprezentowany na zielono, a kolejność elementów jest reprezentowana jak kulki w rurze: pierwsza załadowana jest po prawej, potem druga w środku, potem trzecia po lewej._

#### 2.1.1. Dodaj więcej pozdrowień

Przed blokiem workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // utwórz kanał dla danych wejściowych
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // utwórz kanał dla danych wejściowych
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

Dokumentacja mówi, że to powinno działać. Czy to naprawdę może być takie proste?

#### 2.1.2. Uruchom polecenie i spójrz na wyjście dziennika

Sprawdźmy.

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Z pewnością wygląda na to, że uruchomił się bez problemów.
Monitor wykonywania pokazuje, że `3 of 3` wywołania zostały wykonane dla procesu `sayHello`, i widzimy trzy pozdrowienia wyliczone przez instrukcję `view()`, jedno na linię, jak obiecano.

Jednak w katalogu wyników nadal jest tylko jedno wyjście:

??? abstract "Zawartość katalogu"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Zawartość pliku"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Powinieneś zobaczyć tam jedno z trzech pozdrowień, ale to, które otrzymałeś, może być inne niż pokazane tutaj.
Czy możesz się domyślić, dlaczego tak może być?

Patrząc wstecz na monitor wykonywania, dał nam tylko jedną ścieżkę podkatalogu (`f4/c9962c`).
Zajrzyjmy tam.

??? abstract "Zawartość katalogu"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Zawartość pliku"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

To nawet nie jest to samo pozdrowienie, które mamy w katalogu wyników! Co się dzieje?

W tym momencie musimy Ci powiedzieć, że domyślnie system logowania ANSI zapisuje logowanie z wielu wywołań do tego samego procesu w tej samej linii.
Więc status ze wszystkich trzech wywołań procesu sayHello() ląduje w tym samym miejscu.

Na szczęście możemy wyłączyć to zachowanie, aby zobaczyć pełną listę wywołań procesów.

#### 2.1.3. Uruchom polecenie ponownie z opcją `-ansi-log false`

Aby rozwinąć logowanie do wyświetlania jednej linii na wywołanie procesu, dodaj `-ansi-log false` do polecenia.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Tym razem widzimy wszystkie trzy uruchomienia procesów i ich powiązane podkatalogi work wymienione w wyjściu.

To znacznie lepiej, przynajmniej dla prostego workflow'u.
Dla złożonego workflow'u lub dużej liczby wejść, posiadanie pełnej listy wyświetlanej w terminalu byłoby nieco przytłaczające.
Dlatego `-ansi-log false` nie jest domyślnym zachowaniem.

!!! tip "Wskazówka"

    Sposób raportowania statusu jest nieco inny między dwoma trybami logowania.
    W trybie skondensowanym Nextflow raportuje, czy wywołania zakończyły się pomyślnie, czy nie.
    W tym rozwiniętym trybie raportuje tylko, że zostały przesłane.

W każdym razie, teraz mamy podkatalogi każdego wywołania procesu, możemy szukać ich logów i wyjść.

??? abstract "Zawartość katalogu"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Zawartość pliku"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

To pokazuje, że wszystkie trzy procesy uruchomiły się pomyślnie (juhu).

To powiedziawszy, nadal mamy problem, że w katalogu wyników jest tylko jeden plik wyjściowy.

Możesz pamiętać, że zakodowaliśmy na sztywno nazwę wyjściową dla procesu `sayHello`, więc wszystkie trzy wywołania utworzyły dokument o nazwie `output.txt`.

Dopóki pliki wyjściowe pozostają w podkatalogach work, odizolowane od innych procesów, jest to w porządku.
Ale gdy są publikowane do tego samego katalogu wyników, którykolwiek został tam skopiowany jako pierwszy, jest nadpisywany przez następny, i tak dalej.

### 2.2. Upewnij się, że nazwy plików wyjściowych będą unikalne

Możemy kontynuować publikowanie wszystkich wyjść do tego samego katalogu wyników, ale musimy upewnić się, że będą miały unikalne nazwy.
Konkretnie, musimy zmodyfikować pierwszy proces, aby generował nazwę pliku dynamicznie, tak aby końcowe nazwy plików były unikalne.

Więc jak sprawić, żeby nazwy plików były unikalne?
Powszechnym sposobem na to jest użycie jakiegoś unikalnego fragmentu metadanych z wejść (otrzymanych z kanału wejściowego) jako części nazwy pliku wyjściowego.
Tutaj, dla wygody, użyjemy po prostu samego pozdrowienia, ponieważ to tylko krótki ciąg, i dodamy go na początku podstawowej nazwy pliku wyjściowego.

#### 2.2.1. Skonstruuj dynamiczną nazwę pliku wyjściowego

W bloku procesu wprowadź następujące zmiany kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
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

Upewnij się, że zastąpisz `output.txt` zarówno w definicji `output`, jak i w bloku polecenia `script:`.

!!! tip "Wskazówka"

    W definicji output MUSISZ użyć podwójnych cudzysłowów wokół wyrażenia nazwy pliku wyjściowego (NIE pojedynczych cudzysłowów), w przeciwnym razie to się nie powiedzie.

To powinno generować unikalną nazwę pliku wyjściowego za każdym razem, gdy proces jest wywoływany, tak aby można go było odróżnić od wyjść z innych wywołań tego samego procesu w katalogu wyjściowym.

#### 2.2.2. Uruchom workflow

Uruchommy to. Zauważ, że wracamy do uruchamiania z domyślnymi ustawieniami logu ANSI.

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Wracając do widoku podsumowania, wyjście jest ponownie podsumowane w jednej linii.
Spójrz na katalog `results`, aby zobaczyć, czy wszystkie wyjściowe pozdrowienia tam są.

??? abstract "Zawartość katalogu"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Tak! I każdy ma oczekiwaną zawartość.

??? abstract "Zawartość pliku"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Sukces! Teraz możemy dodawać tyle pozdrowień, ile chcemy, bez martwienia się o nadpisywanie plików wyjściowych.

!!! tip "Wskazówka"

    W praktyce nazywanie plików na podstawie samych danych wejściowych jest prawie zawsze niepraktyczne.
    Lepszym sposobem na generowanie dynamicznych nazw plików jest przekazywanie metadanych do procesu wraz z plikami wejściowymi.
    Metadane są zazwyczaj dostarczane przez 'arkusz próbek' lub odpowiedniki.
    Nauczysz się, jak to zrobić później w Swoim szkoleniu Nextflow (zobacz [Side quest dotyczący metadanych](../side_quests/metadata.md)).

### Podsumowanie

Wiesz, jak przeprowadzać wiele elementów wejściowych przez kanał.

### Co dalej?

Naucz się używać operatora do transformowania zawartości kanału.

---

## 3. Dostarcz wiele wejść przez tablicę

Właśnie pokazaliśmy Ci, jak obsługiwać wiele elementów wejściowych, które były zakodowane na sztywno bezpośrednio w fabryce kanałów.
Co jeśli chcielibyśmy dostarczyć te wiele wejść w inny sposób?

Na przykład, wyobraź sobie, że skonfigurowaliśmy zmienną wejściową zawierającą tablicę elementów w ten sposób:

`greetings_array = ['Hello','Bonjour','Holà']`

Czy możemy załadować to do naszego kanału wyjściowego i oczekiwać, że zadziała?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Dowiedzmy się.

### 3.1. Dostarcz tablicę wartości jako wejście do kanału

Zdrowy rozsądek sugeruje, że powinniśmy móc po prostu przekazać tablicę wartości zamiast pojedynczej wartości.
Spróbujmy; będziemy musieli skonfigurować zmienną wejściową i załadować ją do fabryki kanałów.

#### 3.1.1. Skonfiguruj zmienną wejściową

Weźmy zmienną `greetings_array`, którą właśnie sobie wyobraziłeś, i uczyńmy ją rzeczywistością, dodając ją do bloku workflow:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

To jeszcze nie jest funkcjonalne, dodaliśmy tylko deklarację tablicy.

#### 3.1.2. Ustaw tablicę pozdrowień jako wejście do fabryki kanałów

Teraz zamierzamy zastąpić wartości `'Hello','Bonjour','Holà'` aktualnie zakodowane na sztywno w fabryce kanałów tablicą `greetings_array`, którą właśnie utworzyliśmy.

W bloku workflow wprowadź następującą zmianę:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

To powinno być teraz funkcjonalne.

#### 3.1.3. Uruchom workflow

Spróbujmy go uruchomić:

```bash
nextflow run hello-channels.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

O nie! Jest błąd!

Spójrz na wyjście `view()` i komunikaty o błędach.

Wygląda na to, że Nextflow próbował uruchomić pojedyncze wywołanie procesu, używając `[Hello, Bonjour, Holà]` jako wartości ciągu, zamiast używać trzech ciągów w tablicy jako osobnych wartości.

Więc to 'opakowanie' powoduje problem.
Jak sprawić, żeby Nextflow rozpakował tablicę i załadował poszczególne ciągi do kanału?

### 3.2. Użyj operatora do transformowania zawartości kanału

Tutaj wchodzą do gry **[operatory](https://www.nextflow.io/docs/latest/reference/operator.html)**.
Już używałeś operatora `.view()`, który po prostu patrzy na to, co tam jest.
Teraz przyjrzymy się operatorom, które pozwalają nam działać na zawartości kanału.

Jeśli przejrzysz [listę operatorów](https://www.nextflow.io/docs/latest/reference/operator.html) w dokumentacji Nextflow, znajdziesz [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten), który robi dokładnie to, czego potrzebujemy: rozpakowuje zawartość tablicy i emituje je jako pojedyncze elementy.

#### 3.2.1. Dodaj operator `flatten()`

Aby zastosować operator `flatten()` do naszego kanału wejściowego, dodajemy go do deklaracji fabryki kanałów.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Tutaj dodaliśmy operator w następnej linii dla czytelności, ale możesz dodać operatory w tej samej linii co fabryka kanałów, jeśli wolisz, w ten sposób:
`greeting_ch = channel.of(greetings_array).view().flatten()`

#### 3.2.2. Doprecyzuj instrukcje `view()`

Moglibyśmy uruchomić to od razu, aby sprawdzić, czy działa, ale przy okazji doprecyzujemy sposób, w jaki sprawdzamy zawartość kanału.

Chcemy móc porównać, jak zawartość wygląda przed i po zastosowaniu operatora `flatten()`, więc dodamy drugą instrukcję oraz trochę kodu, aby były wyraźniej oznaczone w wyjściu.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Widzisz, że dodaliśmy drugą instrukcję `.view`, a dla każdej z nich zastąpiliśmy puste nawiasy (`()`) nawiasami klamrowymi zawierającymi kod, taki jak `{ greeting -> "Before flatten: $greeting" }`.

Są to tak zwane _closures_. Kod, który zawierają, będzie wykonywany dla każdego elementu w kanale.
Definiujemy tymczasową zmienną dla wewnętrznej wartości, tutaj nazwaną `greeting` (ale mogłaby mieć dowolną nazwę), która jest używana tylko w zakresie tej closure.

W tym przykładzie `$greeting` reprezentuje każdy pojedynczy element załadowany do kanału.
To spowoduje ładnie oznaczone wyjście konsoli.

!!! info "Informacja"

    W niektórych pipeline'ach możesz zobaczyć specjalną zmienną o nazwie `$it` używaną wewnątrz closures operatorów.
    Jest to _niejawna_ zmienna, która pozwala na skrócony dostęp do wewnętrznej zmiennej,
    bez potrzeby definiowania jej za pomocą `->`.

    Preferujemy być jawni, aby pomóc w czytelności kodu, więc składnia `$it` jest odradzana i będzie stopniowo wycofywana z języka Nextflow.

#### 3.2.3. Uruchom workflow

W końcu możesz spróbować uruchomić workflow ponownie!

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

Tym razem działa I daje nam dodatkowy wgląd w to, jak zawartość kanału wygląda przed i po uruchomieniu operatora `flatten()`.

- Widzisz, że otrzymujemy pojedynczą instrukcję `Before flatten:`, ponieważ w tym momencie kanał zawiera jeden element, oryginalną tablicę.
  Następnie otrzymujemy trzy oddzielne instrukcje `After flatten:`, jedną dla każdego pozdrowienia, które są teraz pojedynczymi elementami w kanale.

Co ważne, oznacza to, że każdy element może być teraz przetwarzany osobno przez workflow.

!!! tip "Wskazówka"

    Technicznie możliwe jest osiągnięcie tych samych wyników przez użycie innej fabryki kanałów, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), która zawiera niejawny krok mapowania w Swojej operacji.
    Tutaj zdecydowaliśmy się tego nie używać, aby zademonstrować użycie operatora na prostym przypadku użycia.

### Podsumowanie

Wiesz, jak używać operatora takiego jak `flatten()` do transformowania zawartości kanału i jak używać operatora `view()` do inspekcji zawartości kanału przed i po zastosowaniu operatora.

### Co dalej?

Naucz się, jak sprawić, żeby workflow przyjmował plik jako źródło wartości wejściowych.

---

## 4. Odczytaj wartości wejściowe z pliku CSV

Realistycznie rzecz biorąc, rzadko, jeśli w ogóle, będziemy zaczynać od tablicy wartości.
Najprawdopodobniej będziemy mieć jeden lub więcej plików zawierających dane, które muszą być przetworzone, w jakimś rodzaju strukturyzowanego formatu.

Przygotowaliśmy plik CSV o nazwie `greetings.csv`, który zawiera kilka pozdrowień wejściowych, naśladując rodzaj danych kolumnowych, które możesz chcieć przetworzyć w rzeczywistej analizie danych, przechowywany w `data/`.
(Liczby nie mają znaczenia, są tam tylko w celach ilustracyjnych.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Naszym następnym zadaniem jest dostosowanie naszego workflow'u do odczytania wartości z tego pliku.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Zobaczmy, jak możemy to zrobić.

### 4.1. Zmodyfikuj skrypt, aby oczekiwał pliku CSV jako źródła pozdrowień

Aby zacząć, musimy wprowadzić dwie kluczowe zmiany w skrypcie:

- Przełączyć parametr wejściowy, aby wskazywał na plik CSV
- Przełączyć fabrykę kanałów na taką, która jest zaprojektowana do obsługi pliku

#### 4.1.1. Przełącz parametr wejściowy, aby wskazywał na plik CSV

Pamiętasz parametr `params.input`, który skonfigurowaliśmy w Części 1?
Zaktualizujemy go, aby wskazywał na plik CSV zawierający nasze pozdrowienia.

Wprowadź następującą edycję deklaracji parametru:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    input: String = 'Holà mundo!'
    ```

To zakłada, że plik jest współlokalizowany z kodem workflow'u.
Nauczysz się, jak radzić sobie z innymi lokalizacjami danych później na swojej drodze z Nextflow.

#### 4.1.2. Przełącz na fabrykę kanałów zaprojektowaną do obsługi pliku

Ponieważ teraz chcemy użyć pliku zamiast prostych ciągów jako wejścia, nie możemy użyć fabryki kanałów `channel.of()` z poprzedniego.
Musimy przełączyć się na użycie nowej fabryki kanałów, [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#channel-path), która ma wbudowaną funkcjonalność do obsługi ścieżek plików.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych from a CSV file
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // zadeklaruj tablicę pozdrowień wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Zauważysz, że przełączyliśmy wejście kanału z powrotem na `param.input` i usunęliśmy deklarację `greetings_array`, ponieważ nie będziemy już jej potrzebować.
Zakomentowaliśmy też `flatten()` i drugą instrukcję `view()`.

#### 4.1.3. Uruchom workflow

Spróbujmy uruchomić workflow z nową fabryką kanałów i plikiem wejściowym.

```bash
nextflow run hello-channels.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

O nie, to nie działa. Spójrz na początek wyjścia konsoli i komunikat o błędzie.
Część `Command executed:` jest tu szczególnie pomocna.

To może wyglądać trochę znajomo.
Wygląda na to, że Nextflow próbował uruchomić pojedyncze wywołanie procesu, używając samej ścieżki pliku jako wartości ciągu.
Więc poprawnie rozpoznał ścieżkę pliku, ale faktycznie nie sparsował jego zawartości, czego chcieliśmy.

Jak sprawić, żeby Nextflow otworzył plik i załadował jego zawartość do kanału?

Brzmi jakbyśmy potrzebowali kolejnego [operatora](https://www.nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Użyj operatora `splitCsv()` do parsowania pliku

Przeglądając ponownie listę operatorów, znajdujemy [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitCsv), który jest zaprojektowany do parsowania i dzielenia tekstu w formacie CSV.

#### 4.2.1. Zastosuj `splitCsv()` do kanału

Aby zastosować operator, dodajemy go do linii fabryki kanałów, jak poprzednio.

W bloku workflow wprowadź następującą zmianę kodu, aby zastąpić `flatten()` (zakomentowane) przez `splitCsv()`:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych from a CSV file
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych from a CSV file
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Jak widzisz, zaktualizowaliśmy też instrukcje before/after `view()`.
Technicznie mogliśmy użyć tej samej nazwy zmiennej (`greeting`), ale zaktualizowaliśmy ją na coś bardziej odpowiedniego (`csv`), aby kod był bardziej czytelny dla innych.

#### 4.2.2. Uruchom workflow ponownie

Spróbujmy uruchomić workflow z dodaną logiką parsowania CSV.

```bash
nextflow run hello-channels.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Co ciekawe, to też się nie udaje, ale z innym błędem.
Tym razem Nextflow sparsował zawartość pliku (juhu!), ale załadował każdy wiersz jako tablicę, a każda tablica jest elementem w kanale.

Musimy mu powiedzieć, żeby wziął tylko pierwszą kolumnę w każdym wierszu.
Więc jak to rozpakować?

Wcześniej użyliśmy `flatten()` do rozpakowania zawartości kanału, ale to by tu nie zadziałało, ponieważ flatten rozpakowuje _wszystko_ (możesz spróbować sam, jeśli chcesz zobaczyć na własne oczy).

Zamiast tego użyjemy innego operatora o nazwie `map()`, który jest naprawdę użyteczny i pojawia się często w pipeline'ach Nextflow.

### 4.3. Użyj operatora `map()` do wyodrębnienia pozdrowień

Operator [`map()`](https://www.nextflow.io/docs/latest/reference/operator.html#map) to bardzo poręczne małe narzędzie, które pozwala nam robić różne mapowania na zawartości kanału.

W tym przypadku użyjemy go do wyodrębnienia tego jednego elementu, który chcemy z każdego wiersza w naszym pliku danych.
Oto jak wygląda składnia:

```groovy title="Składnia"
.map { row -> row[0] }
```

Oznacza to 'dla każdego wiersza w kanale, weź 0-ty (pierwszy) element, który zawiera'.

Więc zastosujmy to do naszego parsowania CSV.

#### 4.3.1. Zastosuj `map()` do kanału

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych from a CSV file
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych from a CSV file
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Widzisz, że dodaliśmy kolejne wywołanie `view()`, aby potwierdzić, że operator robi to, czego oczekujemy.

#### 4.3.2. Uruchom workflow

Uruchommy to jeszcze raz:

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

Tym razem powinno się uruchomić bez błędu.

Patrząc na wyjście instrukcji `view()`, widzisz następujące rzeczy:

- Pojedynczą instrukcję `Before splitCsv:`: w tym momencie kanał zawiera jeden element, oryginalną ścieżkę pliku.
- Trzy oddzielne instrukcje `After splitCsv:`: jedną dla każdego pozdrowienia, ale każde jest zawarte w tablicy odpowiadającej tej linii w pliku.
- Trzy oddzielne instrukcje `After map:`: jedną dla każdego pozdrowienia, które są teraz pojedynczymi elementami w kanale.

Zauważ, że linie mogą pojawiać się w innej kolejności w Twoim wyjściu.

Możesz też spojrzeć na pliki wyjściowe, aby zweryfikować, że każde pozdrowienie zostało poprawnie wyodrębnione i przetworzone przez workflow.

Osiągnęliśmy ten sam wynik co poprzednio, ale teraz mamy znacznie większą elastyczność, aby dodawać więcej elementów do kanału pozdrowień, które chcemy przetworzyć, modyfikując plik wejściowy, bez modyfikowania żadnego kodu.
Nauczysz się bardziej zaawansowanych podejść do obsługi złożonych wejść w późniejszym szkoleniu.

### Podsumowanie

Wiesz, jak używać konstruktora kanału `.fromPath()` i operatorów `splitCsv()` oraz `map()` do wczytania pliku wartości wejściowych i odpowiedniego ich obsłużenia.

Ogólniej rzecz biorąc, masz podstawowe zrozumienie tego, jak Nextflow używa **kanałów** do zarządzania wejściami do procesów i **operatorów** do transformowania ich zawartości.

### Co dalej?

Zrób sobie dużą przerwę, ciężko pracowałeś w tej części!

Gdy będziesz gotowy, przejdź do [**Części 3: Hello Workflow**](./03_hello_workflow.md), aby nauczyć się, jak dodawać więcej kroków i łączyć je w prawdziwy workflow.

---

## Quiz

<quiz>
Czym jest kanał w Nextflow?
- [ ] Specyfikacją ścieżki pliku
- [ ] Definicją procesu
- [x] Strukturą podobną do kolejki do przekazywania danych między procesami
- [ ] Ustawieniem konfiguracyjnym

Dowiedz się więcej: [1.1. Utwórz kanał wejściowy](#11-utworz-kanal-wejsciowy)
</quiz>

<quiz>
Co wyświetli ten kod?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (pojedyncza lista)
- [x] Każdy element w osobnej linii: `Hello`, `Bonjour`, `Hola`
- [ ] Nic (kanały domyślnie nie wypisują)
- [ ] Błąd (nieprawidłowa składnia)

Dowiedz się więcej: [1.1. Utwórz kanał wejściowy](#11-utworz-kanal-wejsciowy)
</quiz>

<quiz>
Gdy kanał zawiera wiele wartości, jak Nextflow obsługuje wykonywanie procesu?
- [ ] Proces uruchamia się raz ze wszystkimi wartościami
- [x] Proces uruchamia się raz dla każdej wartości w kanale
- [ ] Proces uruchamia się tylko z pierwszą wartością
- [ ] Proces uruchamia się tylko z ostatnią wartością

Dowiedz się więcej: [2. Zmodyfikuj workflow, aby działał na wielu wartościach wejściowych](#2-zmodyfikuj-workflow-aby-dzialal-na-wielu-wartosciach-wejsciowych)
</quiz>

<quiz>
Co robi operator `flatten()`?
- [ ] Łączy wiele kanałów w jeden
- [ ] Sortuje elementy kanału
- [x] Rozpakowuje tablice na pojedyncze elementy
- [ ] Usuwa duplikaty elementów

Dowiedz się więcej: [3.2.1. Dodaj operator `flatten()`](#321-dodaj-operator-flatten)
</quiz>

<quiz>
Jaki jest cel operatora `view()`?
- [ ] Filtrowanie zawartości kanału
- [ ] Transformowanie elementów kanału
- [x] Inspekcja i debugowanie zawartości kanału
- [ ] Zapisywanie zawartości kanału do pliku

Dowiedz się więcej: [1.4. Użyj `view()` do inspekcji zawartości kanału](#14-uzyj-view-do-inspekcji-zawartosci-kanalu)
</quiz>

<quiz>
Co robi `splitCsv()`?
- [ ] Tworzy plik CSV z zawartości kanału
- [ ] Dzieli ciąg przez przecinki
- [x] Parsuje plik CSV na tablice reprezentujące każdy wiersz
- [ ] Łączy wiele plików CSV

Dowiedz się więcej: [4.2. Użyj operatora `splitCsv()` do parsowania pliku](#42-uzyj-operatora-splitcsv-do-parsowania-pliku)
</quiz>

<quiz>
Jaki jest cel operatora `map()`?
- [ ] Filtrowanie elementów z kanału
- [ ] Łączenie wielu kanałów
- [x] Transformowanie każdego elementu w kanale
- [ ] Liczenie elementów w kanale

Dowiedz się więcej: [4.3. Użyj operatora `map()` do wyodrębnienia pozdrowień](#43-uzyj-operatora-map-do-wyodrebnienia-pozdrowien)
</quiz>

<quiz>
Dlaczego ważne jest używanie dynamicznych nazw plików wyjściowych przy przetwarzaniu wielu wejść?
- [ ] Aby poprawić wydajność
- [ ] Aby zmniejszyć zużycie miejsca na dysku
- [x] Aby zapobiec nadpisywaniu plików wyjściowych przez siebie nawzajem
- [ ] Aby włączyć funkcjonalność resume

Dowiedz się więcej: [2.2. Upewnij się, że nazwy plików wyjściowych będą unikalne](#22-upewnij-sie-ze-nazwy-plikow-wyjsciowych-beda-unikalne)
</quiz>
