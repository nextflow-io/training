# Część 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/02_hello_channels.md).
///

W Części 1 tego kursu (Hello World) pokazaliśmy Ci, jak przekazać zmienną wejściową do procesu, podając ją bezpośrednio w wywołaniu procesu: `sayHello(params.input)`.
To było celowo uproszczone podejście.
W praktyce takie podejście ma poważne ograniczenia; mianowicie działa tylko w bardzo prostych przypadkach, gdy chcemy uruchomić proces tylko raz, na pojedynczej wartości.
W większości realistycznych przypadków użycia workflow'ów chcemy przetwarzać wiele wartości (na przykład dane eksperymentalne dla wielu próbek), więc potrzebujemy bardziej zaawansowanego sposobu obsługi danych wejściowych.

Do tego służą [**kanały**](https://nextflow.io/docs/latest/channel.html) Nextflow'a.
Kanały to kolejki zaprojektowane do efektywnej obsługi danych wejściowych i przekazywania ich z jednego kroku do drugiego w wieloetapowych workflow'ach, zapewniając jednocześnie wbudowaną paralelizację i wiele dodatkowych korzyści.

W tej części kursu nauczysz się, jak używać kanału do obsługi wielu danych wejściowych z różnych źródeł.
Nauczysz się również używać [**operatorów**](https://nextflow.io/docs/latest/reference/operator.html) do przekształcania zawartości kanałów według potrzeb.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś Część 1 kursu [Hello Nextflow](./index.md), ale jeśli czujesz się komfortowo z podstawami omówionymi w tamtej sekcji, możesz zacząć od tego miejsca bez żadnych specjalnych przygotowań.

---

## 0. Rozgrzewka: Uruchom `hello-channels.nf`

Będziemy używać skryptu workflow'a `hello-channels.nf` jako punktu wyjścia.
Jest on równoważny skryptowi powstałemu w wyniku przejścia przez Część 1 tego kursu szkoleniowego, z tą różnicą, że zmieniliśmy miejsce docelowe wyjścia:

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

Jak poprzednio, znajdziesz plik wyjściowy o nazwie `output.txt` w katalogu `results/hello_channels` (zgodnie z blokiem `output` w skrypcie workflow'a, pokazanym powyżej).

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

Jeśli to zadziałało, jesteś gotowy, aby poznać kanały.

---

## 1. Przekaż zmienne wejściowe przez kanał jawnie

Zamierzamy utworzyć **kanał**, aby przekazać zmienną wejściową do procesu `sayHello()` zamiast polegać na niejawnej obsłudze, która ma pewne ograniczenia.

### 1.1. Utwórz kanał wejściowy

Istnieje wiele [**fabryk kanałów**](https://nextflow.io/docs/latest/reference/channel.html), których możemy użyć do skonfigurowania kanału.
Aby na razie zachować prostotę, użyjemy najbardziej podstawowej fabryki kanałów, zwanej [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of), która utworzy kanał zawierający pojedynczą wartość.
Funkcjonalnie będzie to podobne do tego, jak mieliśmy to skonfigurowane wcześniej, ale zamiast pozwalać Nextflow'owi tworzyć kanał niejawnie, robimy to teraz jawnie.

Oto linia kodu, której użyjemy:

```console title="Składnia"
greeting_ch = channel.of('Hello Channels!')
```

Tworzy to kanał o nazwie `greeting_ch` przy użyciu fabryki kanałów `channel.of()`, która konfiguruje prosty kanał kolejki i ładuje ciąg znaków `'Hello Channels!'` do użycia jako wartość powitania.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Uwaga"

    Tymczasowo wracamy do zakodowanych na stałe ciągów znaków zamiast używać parametru CLI dla czytelności. Wrócimy do używania parametrów CLI, gdy omówimy, co dzieje się na poziomie kanału.

W bloku workflow dodaj kod fabryki kanałów:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
        // wyemituj powitanie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // wyemituj powitanie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

To jeszcze nie jest funkcjonalne, ponieważ nie przełączyliśmy jeszcze wejścia do wywołania procesu.

### 1.2. Dodaj kanał jako wejście do wywołania procesu

Teraz musimy faktycznie podłączyć nasz nowo utworzony kanał do wywołania procesu `sayHello()`, zastępując parametr CLI, który wcześniej podawaliśmy bezpośrednio.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
        // wyemituj powitanie
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
        // wyemituj powitanie
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

To mówi Nextflow'owi, aby uruchomił proces `sayHello` na zawartości kanału `greeting_ch`.

Teraz nasz workflow jest w pełni funkcjonalny; jest to jawny odpowiednik zapisu `sayHello('Hello Channels!')`.

### 1.3. Uruchom workflow

Uruchommy go!

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

Jeśli wykonałeś obie edycje poprawnie, powinieneś uzyskać pomyślne wykonanie.
Możesz sprawdzić katalog wyników, aby upewnić się, że rezultat jest nadal taki sam jak wcześniej.

??? abstract "Zawartość pliku"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Zwiększyliśmy więc elastyczność naszego workflow'a, osiągając ten sam rezultat końcowy.
Może się wydawać, że piszemy więcej kodu bez wymiernych korzyści, ale wartość stanie się jasna, gdy zaczniemy obsługiwać więcej danych wejściowych.

Jako zapowiedź tego, spójrzmy na jeszcze jedną rzecz, zanim przejdziemy dalej: jedną małą, ale wygodną korzyść z używania jawnego kanału do zarządzania danymi wejściowymi.

### 1.4. Użyj `view()` do sprawdzenia zawartości kanału

Kanały Nextflow'a są zbudowane w sposób, który pozwala nam operować na ich zawartości za pomocą operatorów, które omówimy szczegółowo później w tym rozdziale.

Na razie pokażemy Ci tylko, jak używać bardzo prostego operatora o nazwie [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) do sprawdzenia zawartości kanału.
Możesz myśleć o `view()` jako o narzędziu debugowania, podobnym do instrukcji `print()` w Pythonie lub jej odpowiednika w innych językach.

Dodaj tę małą linię do bloku workflow:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // wyemituj powitanie
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
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Dokładna liczba spacji nie ma znaczenia, o ile jest wielokrotnością 4; po prostu staramy się wyrównać początek instrukcji `.view()` do części `.of()` konstrukcji kanału.

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

Jak widzisz, wyświetla to zawartość kanału w konsoli.
Tutaj mamy tylko jeden element, ale gdy zaczniemy ładować wiele wartości do kanału w następnej sekcji, zobaczysz, że jest to ustawione tak, aby wyświetlać jeden element na linię.

### Podsumowanie

Wiesz, jak używać podstawowej fabryki kanałów do przekazywania danych wejściowych do procesu.

### Co dalej?

Naucz się, jak używać kanałów, aby workflow iterował po wielu wartościach wejściowych.

---

## 2. Zmodyfikuj workflow, aby działał na wielu wartościach wejściowych

Workflow'y zazwyczaj działają na partiach danych wejściowych, które mają być przetwarzane zbiorczo, więc chcemy ulepszyć workflow, aby akceptował wiele wartości wejściowych.

### 2.1. Załaduj wiele powitań do kanału wejściowego

Wygodnie, fabryka kanałów `channel.of()`, której używaliśmy, z przyjemnością akceptuje więcej niż jedną wartość, więc nie musimy jej w ogóle modyfikować.
Możemy po prostu załadować wiele wartości do kanału.

Niech będą to `'Hello'`, `'Bonjour'` i `'Holà'`.

#### 2.1.1. Dodaj więcej powitań

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

Dokumentacja mówi nam, że to powinno zadziałać. Czy naprawdę może być tak proste?

#### 2.1.2. Uruchom polecenie i spójrz na wyjście logu

Spróbujmy.

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

Z pewnością wydaje się, że uruchomił się bez problemu.
Monitor wykonania pokazuje, że wykonano `3 z 3` wywołań procesu `sayHello`, a widzimy trzy powitania wyliczone przez instrukcję `view()`, jedno na linię, jak obiecano.

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

Powinieneś zobaczyć jedno z trzech powitań, chociaż to, które otrzymałeś, może różnić się od tego, co pokazano tutaj.
Czy możesz pomyśleć, dlaczego tak może być?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Na diagramie kanał jest reprezentowany na zielono, a kolejność elementów jest przedstawiona jak kulki w rurze: pierwszy załadowany jest po prawej, następnie drugi w środku, a trzeci po lewej._

Patrząc wstecz na monitor wykonania, podał nam tylko jedną ścieżkę podkatalogu (`f4/c9962c`).
Spójrzmy tam.

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

To nawet nie jest to samo powitanie, które otrzymaliśmy w katalogu wyników! Co się dzieje?

W tym momencie musimy Ci powiedzieć, że domyślnie system logowania ANSI zapisuje logi z wielu wywołań tego samego procesu w tym samym miejscu.
Więc status ze wszystkich trzech wywołań procesu `sayHello()` trafia w to samo miejsce.

Na szczęście możemy wyłączyć to zachowanie, aby zobaczyć pełną listę wywołań procesów.

#### 2.1.3. Uruchom polecenie ponownie z opcją `-ansi-log false`

Aby rozwinąć logowanie i wyświetlić jedną linię na wywołanie procesu, dodaj `-ansi-log false` do polecenia.

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

Tym razem widzimy wszystkie trzy uruchomienia procesów i ich powiązane podkatalogi robocze wymienione w wyjściu.

To znacznie lepiej, przynajmniej dla prostego workflow'a.
W przypadku złożonego workflow'a lub dużej liczby danych wejściowych, wyświetlanie pełnej listy w terminalu byłoby nieco przytłaczające.
Dlatego `-ansi-log false` nie jest domyślnym zachowaniem.

!!! tip "Wskazówka"

    Sposób raportowania statusu różni się nieco między dwoma trybami logowania.
    W trybie skondensowanym Nextflow raportuje, czy wywołania zostały zakończone pomyślnie, czy nie.
    W tym rozszerzonym trybie raportuje tylko, że zostały przesłane.

W każdym razie, teraz gdy mamy podkatalogi każdego wywołania procesu, możemy poszukać ich logów i wyjść.

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

??? abstract "Zawartość plików"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

To pokazuje, że wszystkie trzy procesy zostały uruchomione pomyślnie (hurra).

Mimo to nadal mamy problem, że w katalogu wyników jest tylko jeden plik wyjściowy.

Możesz pamiętać, że zakodowaliśmy na stałe nazwę pliku wyjściowego dla procesu `sayHello`, więc wszystkie trzy wywołania wygenerowały plik o nazwie `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Dopóki pliki wyjściowe pozostają w podkatalogach roboczych, odizolowane od innych procesów, jest to w porządku.
Ale gdy są publikowane do tego samego katalogu wyników, ten, który został tam skopiowany jako pierwszy, zostaje nadpisany przez następny i tak dalej.

### 2.2. Upewnij się, że nazwy plików wyjściowych będą unikalne

Możemy nadal publikować wszystkie wyjścia do tego samego katalogu wyników, ale musimy upewnić się, że będą miały unikalne nazwy.
Konkretnie, musimy zmodyfikować pierwszy proces, aby dynamicznie generował nazwę pliku, tak aby końcowe nazwy plików były unikalne.

Jak więc sprawić, aby nazwy plików były unikalne?
Powszechnym sposobem jest użycie unikalnego fragmentu metadanych z danych wejściowych (otrzymanych z kanału wejściowego) jako części nazwy pliku wyjściowego.
Tutaj, dla wygody, użyjemy po prostu samego powitania, ponieważ jest to tylko krótki ciąg znaków, i dodamy go przed podstawową nazwą pliku wyjściowego.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

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

Upewnij się, że zastąpiłeś `output.txt` zarówno w definicji wyjścia, jak i w bloku polecenia `script:`.

!!! tip "Wskazówka"

    W definicji wyjścia MUSISZ użyć podwójnych cudzysłowów wokół wyrażenia nazwy pliku wyjściowego (NIE pojedynczych cudzysłowów), w przeciwnym razie nie zadziała.

To powinno generować unikalną nazwę pliku wyjściowego za każdym razem, gdy proces jest wywoływany, tak aby można było ją odróżnić od wyjść z innych wywołań tego samego procesu w katalogu wyjściowym.

#### 2.2.2. Uruchom workflow

Uruchommy go. Zauważ, że wracamy do uruchamiania z domyślnymi ustawieniami logu ANSI.

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
Spójrz na katalog `results`, aby sprawdzić, czy wszystkie powitania wyjściowe są tam.

??? abstract "Zawartość katalogu"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Tak! I każdy ma oczekiwaną zawartość.

??? abstract "Zawartość plików"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Sukces! Teraz możemy dodać tyle powitań, ile chcemy, nie martwiąc się o nadpisywanie plików wyjściowych.

!!! tip "Wskazówka"

    W praktyce nazywanie plików na podstawie samych danych wejściowych jest prawie zawsze niepraktyczne.
    Lepszym sposobem generowania dynamicznych nazw plików jest przekazywanie metadanych do procesu wraz z plikami wejściowymi.
    Metadane są zazwyczaj dostarczane za pomocą 'arkusza próbek' lub odpowiedników.
    Nauczysz się, jak to zrobić, później w szkoleniu Nextflow (zobacz [Metadata side quest](../side_quests/metadata.md)).

### Podsumowanie

Wiesz, jak przekazać wiele elementów wejściowych przez kanał.

### Co dalej?

Naucz się używać operatora do przekształcania zawartości kanału.

---

## 3. Przekaż wiele danych wejściowych przez tablicę

Właśnie pokazaliśmy Ci, jak obsługiwać wiele elementów wejściowych, które były zakodowane bezpośrednio w fabryce kanałów.
Co jeśli chcielibyśmy przekazać te wiele danych wejściowych w inny sposób?

Na przykład wyobraź sobie, że skonfigurowaliśmy zmienną wejściową zawierającą tablicę elementów w ten sposób:

`greetings_array = ['Hello','Bonjour','Holà']`

Czy możemy załadować to do naszego kanału wyjściowego i oczekiwać, że to zadziała?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Przekonajmy się.

### 3.1. Przekaż tablicę wartości jako wejście do kanału

Zdrowy rozsądek sugeruje, że powinniśmy móc po prostu przekazać tablicę wartości zamiast pojedynczej wartości.
Spróbujmy; będziemy musieli skonfigurować zmienną wejściową i załadować ją do fabryki kanałów.

#### 3.1.1. Skonfiguruj zmienną wejściową

Weźmy zmienną `greetings_array`, którą właśnie sobie wyobraziliśmy, i uczyńmy ją rzeczywistością, dodając ją do bloku workflow:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // wyemituj powitanie
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
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

To jeszcze nie jest funkcjonalne, po prostu dodaliśmy deklarację tablicy.

#### 3.1.2. Ustaw tablicę powitań jako wejście do fabryki kanałów

Teraz zamierzamy zastąpić wartości `'Hello','Bonjour','Holà'` obecnie zakodowane w fabryce kanałów zmienną `greetings_array`, którą właśnie utworzyliśmy.

W bloku workflow wprowadź następującą zmianę:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

To powinno być teraz funkcjonalne.

#### 3.1.3. Uruchom workflow

Spróbujmy to uruchomić:

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

Wygląda na to, że Nextflow próbował uruchomić pojedyncze wywołanie procesu, używając `[Hello, Bonjour, Holà]` jako pojedynczej wartości ciągu, zamiast używać trzech ciągów w tablicy jako oddzielnych wartości.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Więc to 'opakowanie' powoduje problem.
Jak sprawić, aby Nextflow rozpakował tablicę i załadował poszczególne ciągi do kanału?

### 3.2. Użyj operatora do przekształcenia zawartości kanału

Tu wchodzą w grę [**operatory**](https://nextflow.io/docs/latest/reference/operator.html).
Już używałeś operatora `.view()`, który po prostu patrzy, co tam jest.
Teraz przyjrzymy się operatorom, które pozwalają nam działać na zawartości kanału.

Jeśli przejrzysz [listę operatorów](https://nextflow.io/docs/latest/reference/operator.html) w dokumentacji Nextflow'a, znajdziesz [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten), który robi dokładnie to, czego potrzebujemy: rozpakować zawartość tablicy i wyemitować je jako pojedyncze elementy.

#### 3.2.1. Dodaj operator `flatten()`

Aby zastosować operator `flatten()` do naszego kanału wejściowego, dołączamy go do deklaracji fabryki kanałów.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Tutaj dodaliśmy operator w następnej linii dla czytelności, ale możesz dodać operatory w tej samej linii co fabryka kanałów, jeśli wolisz, w ten sposób:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. Udoskonal instrukcje `view()`

Moglibyśmy to uruchomić od razu, aby sprawdzić, czy działa, ale przy okazji zamierzamy udoskonalić sposób, w jaki sprawdzamy zawartość kanału.

Chcemy móc skontrastować, jak wygląda zawartość przed i po zastosowaniu operatora `flatten()`, więc dodamy drugi, I dodamy trochę kodu, aby uzyskać wyraźniejsze etykiety w wyjściu.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Przed flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Po flatten: $greeting" }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Widzisz, że dodaliśmy drugą instrukcję `.view`, a dla każdej z nich zastąpiliśmy puste nawiasy (`()`) nawiasami klamrowymi zawierającymi kod, takimi jak `{ greeting -> "Przed flatten: $greeting" }`.

Nazywa się to _domknięciami_. Kod, który zawierają, będzie wykonywany dla każdego elementu w kanale.
Definiujemy zmienną tymczasową dla wartości wewnętrznej, tutaj nazwaną `greeting` (ale może to być dowolna nazwa), która jest używana tylko w zakresie tego domknięcia.

W tym przykładzie `$greeting` reprezentuje każdy pojedynczy element załadowany w kanale.
To spowoduje ładnie oznaczone wyjście konsoli.

!!! info

    W niektórych pipeline'ach możesz zobaczyć specjalną zmienną o nazwie `$it` używaną wewnątrz domknięć operatorów.
    Jest to zmienna _niejawna_, która pozwala na skrócony dostęp do zmiennej wewnętrznej,
    bez konieczności definiowania jej za pomocą `->`.

    Wolimy być jawni, aby pomóc w przejrzystości kodu, dlatego składnia `$it` jest odradzana i będzie stopniowo wycofywana z języka Nextflow.

#### 3.2.3. Uruchom workflow

Na koniec możesz spróbować uruchomić workflow ponownie!

```bash
nextflow run hello-channels.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Przed flatten: [Hello, Bonjour, Holà]
    Po flatten: Hello
    Po flatten: Bonjour
    Po flatten: Holà
    ```

Tym razem działa I daje nam dodatkowy wgląd w to, jak wygląda zawartość kanału przed i po uruchomieniu operatora `flatten()`.

- Pojedyncza instrukcja `Przed flatten:`, ponieważ w tym momencie kanał zawiera jeden element, oryginalną tablicę.
- Trzy oddzielne instrukcje `Po flatten:`, po jednej dla każdego powitania, które są teraz pojedynczymi elementami w kanale.

Co ważne, oznacza to, że każdy element może być teraz przetwarzany oddzielnie przez workflow.

!!! tip "Wskazówka"

    Technicznie możliwe jest osiągnięcie tych samych rezultatów przy użyciu innej fabryki kanałów, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), która zawiera niejawny krok mapowania w swojej operacji.
    Tutaj zdecydowaliśmy się tego nie używać, aby zademonstrować użycie operatora w prostym przypadku użycia.

### Podsumowanie

Wiesz, jak używać operatora takiego jak `flatten()` do przekształcania zawartości kanału i jak używać operatora `view()` do sprawdzania zawartości kanału przed i po zastosowaniu operatora.

### Co dalej?

Naucz się, jak sprawić, aby workflow przyjmował plik jako źródło wartości wejściowych.

---

## 4. Odczytaj wartości wejściowe z pliku CSV

Realistycznie, rzadko, jeśli w ogóle, będziemy zaczynać od tablicy wartości.
Najprawdopodobniej będziemy mieli jeden lub więcej plików zawierających dane, które muszą być przetworzone, w jakimś ustrukturyzowanym formacie.

Przygotowaliśmy plik CSV o nazwie `greetings.csv`, który zawiera kilka powitań wejściowych, naśladując rodzaj danych kolumnowych, które możesz chcieć przetworzyć w rzeczywistej analizie danych, przechowywanych w `data/`.
(Liczby nie mają znaczenia, są tam tylko w celach ilustracyjnych.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Naszym następnym zadaniem jest dostosowanie naszego workflow'a do odczytu wartości z tego pliku.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Zobaczmy, jak możemy to osiągnąć.

### 4.1. Zmodyfikuj skrypt, aby oczekiwał pliku CSV jako źródła powitań

Na początek będziemy musieli wprowadzić dwie kluczowe zmiany w skrypcie:

- Przełączyć parametr wejściowy, aby wskazywał na plik CSV
- Przełączyć fabrykę kanałów na taką, która jest zaprojektowana do obsługi pliku

#### 4.1.1. Przełącz parametr wejściowy, aby wskazywał na plik CSV

Pamiętasz parametr `params.input`, który skonfigurowaliśmy w Części 1?
Zaktualizujemy go, aby wskazywał na plik CSV zawierający nasze powitania.

Wprowadź następującą edycję do deklaracji parametru:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Parametry pipeline'u
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Parametry pipeline'u
     */
    input: String = 'Holà mundo!'
    ```

Zakłada to, że plik znajduje się w tym samym miejscu co kod workflow'a.
Nauczysz się, jak radzić sobie z innymi lokalizacjami danych później w Twojej drodze z Nextflow.

#### 4.1.2. Przełącz na fabrykę kanałów zaprojektowaną do obsługi pliku

Ponieważ teraz chcemy użyć pliku zamiast prostych ciągów jako wejścia, nie możemy użyć fabryki kanałów `channel.of()` z wcześniej.
Musimy przełączyć się na używanie nowej fabryki kanałów, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath), która ma wbudowaną funkcjonalność do obsługi ścieżek plików.

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Przed flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Po flatten: $greeting" }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // zadeklaruj tablicę powitań wejściowych
        greetings_array = ['Hello','Bonjour','Holà']
        // utwórz kanał dla danych wejściowych
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Przed flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "Po flatten: $greeting" }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Zauważysz, że przełączyliśmy wejście kanału z powrotem na `param.input` i usunęliśmy deklarację `greetings_array`, ponieważ nie będziemy jej już potrzebować.
Zakomentowaliśmy również `flatten()` i drugą instrukcję `view()`.

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
    Przed flatten: /workspaces/training/hello-nextflow/data/greetings.csv
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

O nie, nie działa. Spójrz na początek wyjścia konsoli i komunikat o błędzie.
Część `Command executed:` jest szczególnie pomocna tutaj.

To może wyglądać trochę znajomo.
Wygląda na to, że Nextflow próbował uruchomić pojedyncze wywołanie procesu, używając samej ścieżki pliku jako wartości ciągu.
Więc poprawnie rozwiązał ścieżkę pliku, ale faktycznie nie przeanalizował jego zawartości, czego chcieliśmy.

Jak sprawić, aby Nextflow otworzył plik i załadował jego zawartość do kanału?

Wygląda na to, że potrzebujemy kolejnego [operatora](https://nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Użyj operatora `splitCsv()` do przeanalizowania pliku

Przeglądając ponownie listę operatorów, znajdujemy [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv), który jest zaprojektowany do analizowania i dzielenia tekstu w formacie CSV.

#### 4.2.1. Zastosuj `splitCsv()` do kanału

Aby zastosować operator, dołączamy go do linii fabryki kanałów jak wcześniej.

W bloku workflow wprowadź następującą zmianę kodu, aby zastąpić `flatten()` przez `splitcsv()` (odkomentowane):

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Przed splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Po splitCsv: $csv" }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Przed flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "Po flatten: $greeting" }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Jak widzisz, zaktualizowaliśmy również instrukcje `view()` przed/po.
Technicznie mogliśmy użyć tej samej nazwy zmiennej (`greeting`), ale zaktualizowaliśmy ją do czegoś bardziej odpowiedniego (`csv`), aby kod był bardziej czytelny dla innych.

#### 4.2.2. Uruchom workflow ponownie

Spróbujmy uruchomić workflow z dodaną logiką analizowania CSV.

```bash
nextflow run hello-channels.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Przed splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Po splitCsv: [Hello, English, 123]
    Po splitCsv: [Bonjour, French, 456]
    Po splitCsv: [Holà, Spanish, 789]
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

Co ciekawe, to również się nie powiodło, ale z innym błędem.
Tym razem Nextflow przeanalizował zawartość pliku (hurra!), ale załadował każdy wiersz jako tablicę, a każda tablica jest elementem w kanale.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Musimy powiedzieć mu, aby wziął tylko pierwszą kolumnę w każdym wierszu.
Jak więc to rozpakować?

Wcześniej używaliśmy `flatten()` do rozpakowania zawartości kanału, ale to by tutaj nie zadziałało, ponieważ flatten rozpakuje _wszystko_ (możesz spróbować, jeśli chcesz zobaczyć sam).

Zamiast tego użyjemy innego operatora o nazwie `map()`, który jest naprawdę przydatny i pojawia się często w pipeline'ach Nextflow.

### 4.3. Użyj operatora `map()` do wyodrębnienia powitań

Operator [`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) to bardzo przydatne małe narzędzie, które pozwala nam wykonywać wszelkiego rodzaju mapowania na zawartości kanału.

W tym przypadku użyjemy go do wyodrębnienia tego jednego elementu, który chcemy z każdego wiersza w naszym pliku danych.
Oto jak wygląda składnia:

```groovy title="Składnia"
.map { row -> row[0] }
```

To oznacza 'dla każdego wiersza w kanale, weź 0-ty (pierwszy) element, który zawiera'.

Zastosujmy to więc do naszego analizowania CSV.

#### 4.3.1. Zastosuj `map()` do kanału

W bloku workflow wprowadź następującą zmianę kodu:

=== "Po"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Przed splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Po splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "Po map: $csv" }
        // wyemituj powitanie
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Przed"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Przed splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "Po splitCsv: $csv" }
        // wyemituj powitanie
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
    Przed splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    Po splitCsv: [Hello, English, 123]
    Po splitCsv: [Bonjour, French, 456]
    Po splitCsv: [Holà, Spanish, 789]
    Po map: Hello
    Po map: Bonjour
    Po map: Holà
    ```

Tym razem powinno się uruchomić bez błędu.

Patrząc na wyjście instrukcji `view()`, widzisz następujące:

- Pojedyncza instrukcja `Przed splitCsv:`: w tym momencie kanał zawiera jeden element, oryginalną ścieżkę pliku.
- Trzy oddzielne instrukcje `Po splitCsv:`: po jednej dla każdego powitania, ale każde jest zawarte w tablicy, która odpowiada temu wierszowi w pliku.
- Trzy oddzielne instrukcje `Po map:`: po jednej dla każdego powitania, które są teraz pojedynczymi elementami w kanale.

_Zauważ, że linie mogą pojawić się w innej kolejności w Twoim wyjściu._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

Możesz również spojrzeć na pliki wyjściowe, aby sprawdzić, czy każde powitanie zostało poprawnie wyodrębnione i przetworzone przez workflow.

Osiągnęliśmy ten sam rezultat co wcześniej, ale teraz mamy znacznie większą elastyczność, aby dodać więcej elementów do kanału powitań, które chcemy przetworzyć, modyfikując plik wejściowy, bez modyfikowania jakiegokolwiek kodu.
Nauczysz się bardziej zaawansowanych podejść do obsługi złożonych danych wejściowych w późniejszym szkoleniu.

### Podsumowanie

Wiesz, jak używać konstruktora kanałów `.fromPath()` oraz operatorów `splitCsv()` i `map()` do odczytu pliku wartości wejściowych i odpowiedniej ich obsługi.

Bardziej ogólnie, masz podstawowe zrozumienie, jak Nextflow używa **kanałów** do zarządzania danymi wejściowymi do procesów i **operatorów** do przekształcania ich zawartości.
Widziałeś również, jak kanały obsługują równoległe wykonywanie niejawnie.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### Co dalej?

Zrób dużą przerwę, ciężko pracowałeś w tej części!

Kiedy będziesz gotowy, przejdź do [**Części 3: Hello Workflow**](./03_hello_workflow.md), aby nauczyć się, jak dodać więcej kroków i połączyć je razem w odpowiedni workflow.

---

## Quiz

<quiz>
Czym jest kanał w Nextflow?
- [ ] Specyfikacją ścieżki pliku
- [ ] Definicją procesu
- [x] Strukturą podobną do kolejki do przekazywania danych między procesami
- [ ] Ustawieniem konfiguracji

Dowiedz się więcej: [1.1. Utwórz kanał wejściowy](#11-utwórz-kanał-wejściowy)
</quiz>

<quiz>
Co wyświetli ten kod?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (pojedyncza lista)
- [x] Każdy element w osobnej linii: `Hello`, `Bonjour`, `Hola`
- [ ] Nic (kanały domyślnie nie drukują)
- [ ] Błąd (nieprawidłowa składnia)

Dowiedz się więcej: [1.1. Utwórz kanał wejściowy](#11-utwórz-kanał-wejściowy)
</quiz>

<quiz>
Gdy kanał zawiera wiele wartości, jak Nextflow obsługuje wykonywanie procesu?
- [ ] Proces uruchamia się raz ze wszystkimi wartościami
- [x] Proces uruchamia się raz dla każdej wartości w kanale
- [ ] Proces uruchamia się tylko z pierwszą wartością
- [ ] Proces uruchamia się tylko z ostatnią wartością

Dowiedz się więcej: [2. Zmodyfikuj workflow, aby działał na wielu wartościach wejściowych](#2-zmodyfikuj-workflow-aby-działał-na-wielu-wartościach-wejściowych)
</quiz>

<quiz>
Co robi operator `flatten()`?
- [ ] Łączy wiele kanałów w jeden
- [ ] Sortuje elementy kanału
- [x] Rozpakuje tablice na pojedyncze elementy
- [ ] Usuwa zduplikowane elementy

Dowiedz się więcej: [3.2.1. Dodaj operator `flatten()`](#321-dodaj-operator-flatten)
</quiz>

<quiz>
Jaki jest cel operatora `view()`?
- [ ] Filtrowanie zawartości kanału
- [ ] Przekształcanie elementów kanału
- [x] Sprawdzanie i debugowanie zawartości kanału
- [ ] Zapisywanie zawartości kanału do pliku

Dowiedz się więcej: [1.4. Użyj `view()` do sprawdzenia zawartości kanału](#14-użyj-view-do-sprawdzenia-zawartości-kanału)
</quiz>

<quiz>
Co robi `splitCsv()`?
- [ ] Tworzy plik CSV z zawartości kanału
- [ ] Dzieli ciąg znaków po przecinkach
- [x] Analizuje plik CSV na tablice reprezentujące każdy wiersz
- [ ] Łączy wiele plików CSV

Dowiedz się więcej: [4.2. Użyj operatora `splitCsv()` do przeanalizowania pliku](#42-użyj-operatora-splitcsv-do-przeanalizowania-pliku)
</quiz>

<quiz>
Jaki jest cel operatora `map()`?
- [ ] Filtrowanie elementów z kanału
- [ ] Łączenie wielu kanałów
- [x] Przekształcanie każdego elementu w kanale
- [ ] Liczenie elementów w kanale

Dowiedz się więcej: [4.3. Użyj operatora `map()` do wyodrębnienia powitań](#43-użyj-operatora-map-do-wyodrębnienia-powitań)
</quiz>

<quiz>
Dlaczego ważne jest używanie dynamicznych nazw plików wyjściowych podczas przetwarzania wielu danych wejściowych?
- [ ] Aby poprawić wydajność
- [ ] Aby zmniejszyć przestrzeń dyskową
- [x] Aby zapobiec nadpisywaniu się plików wyjściowych
- [ ] Aby włączyć funkcjonalność wznowienia

Dowiedz się więcej: [2.2. Upewnij się, że nazwy plików wyjściowych będą unikalne](#22-upewnij-się-że-nazwy-plików-wyjściowych-będą-unikalne)
</quiz>
