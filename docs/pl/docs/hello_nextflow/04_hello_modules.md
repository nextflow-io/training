# Część 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zaproponuj poprawki](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Zobacz [całą playlistę](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) na kanale YouTube Nextflow.

:green_book: Transkrypcja wideo jest dostępna [tutaj](./transcripts/04_hello_modules.md).
///

Ta sekcja omawia, jak zorganizować kod workflow'a, aby uczynić rozwój i utrzymanie Twojego pipeline'u bardziej efektywnym i trwałym.
W szczególności pokażemy, jak używać [**modułów**](https://nextflow.io/docs/latest/module.html).

W Nextflow'ie **moduł** to samodzielny plik kodu, często zawierający definicję pojedynczego procesu.
Aby użyć modułu w workflow'ie, wystarczy dodać jedną linię z instrukcją `include` do pliku z kodem workflow'a; następnie możesz zintegrować proces z workflow'em w taki sam sposób, jak zwykle.
Dzięki temu możliwe jest ponowne wykorzystanie definicji procesów w wielu workflow'ach bez tworzenia wielu kopii kodu.

Kiedy zaczęliśmy rozwijać nasz workflow, napisaliśmy wszystko w jednym pliku kodu.
Teraz przeniesiemy procesy do osobnych modułów.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

To sprawi, że nasz kod będzie bardziej współdzielny, elastyczny i łatwiejszy w utrzymaniu.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś Części 1-3 kursu [Hello Nextflow](./index.md), ale jeśli czujesz się komfortowo z podstawami omówionymi w tych sekcjach, możesz zacząć od tego miejsca bez żadnych specjalnych przygotowań.

---

## 0. Rozgrzewka: Uruchom `hello-modules.nf`

Jako punkt wyjścia użyjemy skryptu workflow'a `hello-modules.nf`.
Jest on równoważny skryptowi powstałemu w wyniku przejścia przez Część 3 tego kursu szkoleniowego, z tą różnicą, że zmieniliśmy miejsca docelowe wyjścia:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Aby upewnić się, że wszystko działa, uruchom skrypt raz przed wprowadzeniem jakichkolwiek zmian:

```bash
nextflow run hello-modules.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Jak poprzednio, pliki wyjściowe znajdziesz w katalogu określonym w bloku `output` (tutaj: `results/hello_modules/`).

??? abstract "Zawartość katalogu"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Jeśli to zadziałało, jesteś gotowy, aby nauczyć się modularyzować kod workflow'a.

---

## 1. Utwórz katalog do przechowywania modułów

Najlepszą praktyką jest przechowywanie modułów w dedykowanym katalogu.
Możesz nazwać ten katalog dowolnie, ale konwencją jest nazywanie go `modules/`.

```bash
mkdir modules
```

---

## 2. Utwórz moduł dla `sayHello()`

W najprostszej formie przekształcenie istniejącego procesu w moduł to niewiele więcej niż operacja kopiuj-wklej.
Utworzymy szkielet pliku dla modułu, skopiujemy odpowiedni kod, a następnie usuniemy go z głównego pliku workflow'a.

Następnie wystarczy dodać instrukcję `include`, aby Nextflow wiedział, że ma pobrać odpowiedni kod w czasie wykonania.

### 2.1. Utwórz szkielet pliku dla nowego modułu

Utwórzmy pusty plik dla modułu o nazwie `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

To daje nam miejsce na umieszczenie kodu procesu.

### 2.2. Przenieś kod procesu `sayHello` do pliku modułu

Skopiuj całą definicję procesu z pliku workflow'a do pliku modułu.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Użyj echo, aby wypisać 'Hello World!' do pliku
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

Po wykonaniu tej czynności usuń definicję procesu z pliku workflow'a.

### 2.3. Dodaj deklarację include przed blokiem workflow

Składnia dołączania procesu z modułu jest dość prosta:

```groovy title="Składnia: deklaracja include"
include { <NAZWA_PROCESU> } from '<ścieżka_do_modułu>'
```

Wstawmy to powyżej bloku `params` i wypełnijmy odpowiednio.

=== "Po"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Parametry pipeline'u
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Przed"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Parametry pipeline'u
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Jak widzisz, wypełniliśmy nazwę procesu, `sayHello`, oraz ścieżkę do pliku zawierającego kod modułu, `./modules/sayHello.nf`.

### 2.4. Uruchom workflow

Uruchamiamy workflow z zasadniczo tym samym kodem i danymi wejściowymi co wcześniej, więc uruchommy z flagą `-resume` i zobaczmy, co się stanie.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Powinno to wykonać się bardzo szybko, ponieważ wszystko jest w pamięci podręcznej.
Możesz sprawdzić opublikowane wyniki.

Nextflow rozpoznał, że nadal jest to ta sama praca do wykonania, nawet jeśli kod jest podzielony na wiele plików.

### Podsumowanie

Wiesz, jak wyodrębnić proces do lokalnego modułu i wiesz, że nie wpływa to na możliwość wznowienia workflow'a.

### Co dalej?

Poćwicz tworzenie kolejnych modułów.
Kiedy zrobisz jeden, możesz zrobić milion więcej...
Ale na razie zróbmy jeszcze dwa.

---

## 3. Zmodularyzuj proces `convertToUpper()`

### 3.1. Utwórz szkielet pliku dla nowego modułu

Utwórz pusty plik dla modułu o nazwie `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Przenieś kod procesu `convertToUpper` do pliku modułu

Skopiuj całą definicję procesu z pliku workflow'a do pliku modułu.

```groovy title="modules/convertToUpper.nf" linenums="1"
/*
 * Użyj narzędzia do zamiany tekstu, aby przekonwertować powitanie na wielkie litery
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

Po wykonaniu tej czynności usuń definicję procesu z pliku workflow'a.

### 3.3. Dodaj deklarację include przed blokiem `params`

Wstaw deklarację include powyżej bloku `params` i wypełnij ją odpowiednio.

=== "Po"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Parametry pipeline'u
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Przed"

    ```groovy title="hello-modules.nf" linenums="23"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Parametry pipeline'u
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

To powinno zacząć wyglądać bardzo znajomo.

### 3.4. Uruchom workflow ponownie

Uruchom to z flagą `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Powinno to nadal generować takie samo wyjście jak wcześniej.

Dwa gotowe, jeszcze jeden!

---

## 4. Zmodularyzuj proces `collectGreetings()`

### 4.1. Utwórz szkielet pliku dla nowego modułu

Utwórz pusty plik dla modułu o nazwie `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Przenieś kod procesu `collectGreetings` do pliku modułu

Skopiuj całą definicję procesu z pliku workflow'a do pliku modułu.

```groovy title="modules/collectGreetings.nf" linenums="1"
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
```

Po wykonaniu tej czynności usuń definicję procesu z pliku workflow'a.

### 4.3. Dodaj deklarację include przed blokiem `params`

Wstaw deklarację include powyżej bloku `params` i wypełnij ją odpowiednio.

=== "Po"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Parametry pipeline'u
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Przed"

    ```groovy title="hello-modules.nf" linenums="3"
    // Dołącz moduły
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Parametry pipeline'u
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Ostatni!

### 4.4. Uruchom workflow

Uruchom to z flagą `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Powinno to nadal generować takie samo wyjście jak wcześniej.

### Podsumowanie

Wiesz, jak zmodularyzować wiele procesów w workflow'ie.

Gratulacje, wykonałeś całą tę pracę i absolutnie nic się nie zmieniło w sposobie działania pipeline'u!

Żarty na bok, teraz Twój kod jest bardziej modularny, a jeśli zdecydujesz się napisać kolejny pipeline, który wywołuje jeden z tych procesów, wystarczy wpisać jedną krótką instrukcję `include`, aby użyć odpowiedniego modułu.
To lepsze niż kopiowanie i wklejanie kodu, ponieważ jeśli później zdecydujesz się ulepszyć moduł, wszystkie Twoje pipeline'y odziedziczą te ulepszenia.

### Co dalej?

Zrób sobie krótką przerwę, jeśli masz na to ochotę.

Kiedy będziesz gotowy, przejdź do [**Części 5: Hello Containers**](./05_hello_containers.md), aby dowiedzieć się, jak używać kontenerów do wygodniejszego i bardziej odtwarzalnego zarządzania zależnościami oprogramowania.

---

## Quiz

<quiz>
Czym jest moduł w Nextflow'ie?
- [ ] Plikiem konfiguracyjnym
- [x] Samodzielnym plikiem, który może zawierać definicje procesów
- [ ] Definicją workflow'a
- [ ] Operatorem kanału

Dowiedz się więcej: [2. Utwórz moduł dla `sayHello()`](#2-utworz-modul-dla-sayhello)
</quiz>

<quiz>
Jaka konwencja jest zazwyczaj stosowana do przechowywania plików modułów?
- [ ] W tym samym katalogu co workflow
- [ ] W katalogu `bin/`
- [x] W katalogu `modules/`
- [ ] W katalogu `lib/`

Dowiedz się więcej: [1. Utwórz katalog do przechowywania modułów](#1-utworz-katalog-do-przechowywania-modulow)
</quiz>

<quiz>
Jaka jest poprawna składnia użycia modułu?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Dowiedz się więcej: [2.3. Dodaj deklarację include](#23-dodaj-deklaracje-include-przed-blokiem-workflow)
</quiz>

<quiz>
Co dzieje się z funkcjonalnością `-resume` podczas używania modułów?
- [ ] Przestaje działać
- [ ] Wymaga dodatkowej konfiguracji
- [x] Działa tak samo jak wcześniej
- [ ] Działa tylko dla modułów lokalnych
</quiz>

<quiz>
Jakie są korzyści z używania modułów? (Zaznacz wszystkie poprawne odpowiedzi)
- [x] Możliwość ponownego wykorzystania kodu w różnych workflow'ach
- [x] Łatwiejsze utrzymanie
- [x] Lepsza organizacja kodu workflow'a
- [ ] Szybsze wykonanie
</quiz>
