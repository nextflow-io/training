# Przetwarzanie plików wejściowych

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Workflow'y do analiz naukowych często wymagają przetwarzania dużej liczby plików.
Nextflow dostarcza potężnych narzędzi do wydajnej obsługi plików, pomagając organizować i przetwarzać dane przy minimalnej ilości kodu.

### Cele szkolenia

W tym side queście przyjrzymy się, jak Nextflow obsługuje pliki — od podstawowych operacji po bardziej zaawansowane techniki pracy z kolekcjami plików.
Nauczysz się wyodrębniać metadane z nazw plików, co jest powszechnym wymogiem w pipeline'ach do analiz naukowych.

Po ukończeniu tego side questa będziesz potrafić:

- Tworzyć obiekty Path z ciągów znaków reprezentujących ścieżki, używając metody `file()` Nextflow'a
- Uzyskiwać dostęp do atrybutów pliku, takich jak nazwa, rozszerzenie i katalog nadrzędny
- Obsługiwać zarówno pliki lokalne, jak i zdalne w sposób przezroczysty, używając URI
- Używać kanałów do automatyzacji obsługi plików za pomocą `channel.fromPath()` i `channel.fromFilePairs()`
- Wyodrębniać i strukturyzować metadane z nazw plików przy użyciu manipulacji ciągami znaków
- Grupować powiązane pliki za pomocą dopasowywania wzorców i wyrażeń glob
- Integrować operacje na plikach z procesami Nextflow'a z właściwą obsługą wejścia
- Organizować wyjścia procesów przy użyciu struktur katalogów opartych na metadanych

Te umiejętności pomogą Ci budować workflow'y zdolne do obsługi różnych rodzajów plików wejściowych z dużą elastycznością.

### Wymagania wstępne

Przed przystąpieniem do tego side questa powinieneś:

- Ukończyć samouczek [Hello Nextflow](../../hello_nextflow/) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow'a (procesy, kanały, operatory).

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/working_with_files
```

Możesz ustawić VSCode tak, aby skupiał się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz tu prosty plik workflow'u o nazwie `main.nf`, katalog `modules` zawierający dwa pliki modułów oraz katalog `data` z przykładowymi plikami danych.

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Ten katalog zawiera dane sekwencjonowania paired-end od trzech pacjentów (A, B, C).

Dla każdego pacjenta mamy próbki typu `tumor` (zazwyczaj pochodzące z biopsji guza) lub `normal` (pobrane ze zdrowej tkanki lub krwi).
Jeśli nie jesteś zaznajomiony z analizą nowotworów, wiedz tylko, że odpowiada to modelowi eksperymentalnemu, który wykorzystuje sparowane próbki guz/normalny do przeprowadzania analiz kontrastywnych.

Dla pacjenta A mamy dwa zestawy replik technicznych (powtórzeń).

Pliki z danymi sekwencjonowania są nazwane zgodnie z typową konwencją `_R1_` i `_R2_` dla tak zwanych „odczytów w przód" i „odczytów wstecz".

_Nie martw się, jeśli nie jesteś zaznajomiony z tym projektem eksperymentalnym — nie jest to kluczowe dla zrozumienia tego samouczka._

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest napisanie workflow'u Nextflow'a, który:

1. **Wczyta** pliki wejściowe przy użyciu metod obsługi plików Nextflow'a
2. **Wyodrębni** metadane (ID pacjenta, replikę, typ próbki) ze struktury nazwy pliku
3. **Pogrupuje** sparowane pliki (R1/R2) przy użyciu `channel.fromFilePairs()`
4. **Przetworzy** pliki za pomocą dostarczonego modułu analitycznego
5. **Zorganizuje** wyjścia w strukturę katalogów opartą na wyodrębnionych metadanych

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiedni katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Podstawowe operacje na plikach

### 1.1. Identyfikacja typu obiektu za pomocą `.class`

Przyjrzyj się plikowi workflow'u `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

To mini-workflow (bez żadnych procesów), który odwołuje się do ścieżki pojedynczego pliku, a następnie wypisuje ją na konsolę wraz z jej klasą.

??? info "Czym jest `.class`?"

    W Nextflow'ie `.class` mówi nam, z jakim typem obiektu mamy do czynienia. To jak pytanie „czym to jest?" — żeby sprawdzić, czy to ciąg znaków, liczba, plik czy coś innego.
    Pomoże nam to zilustrować różnicę między zwykłym ciągiem znaków a obiektem Path w kolejnych sekcjach.

Uruchommy workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Jak widać, Nextflow wypisał ścieżkę dokładnie tak, jak ją napisaliśmy.

To tylko wyjście tekstowe — Nextflow nie zrobił z nim jeszcze nic specjalnego.
Potwierdziliśmy też, że z punktu widzenia Nextflow'a jest to jedynie ciąg znaków (klasy `java.lang.String`).
Ma to sens, ponieważ nie powiedzieliśmy jeszcze Nextflow'owi, że odpowiada on plikowi.

### 1.2. Tworzenie obiektu Path za pomocą file()

Możemy poinformować Nextflow'a, jak obsługiwać pliki, tworząc [obiekty Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) z ciągów znaków reprezentujących ścieżki.

W naszym workflow'ie możemy przekonwertować ciąg znaków `data/patientA_rep1_normal_R1_001.fastq.gz` na obiekt Path przy użyciu metody `file()`, która zapewnia dostęp do właściwości i operacji na plikach.

Edytuj `main.nf`, opakowując ciąg znaków w `file()` w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Teraz uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Tym razem widzisz pełną ścieżkę bezwzględną zamiast ścieżki względnej, którą podaliśmy jako wejście.

Nextflow przekonwertował nasz ciąg znaków na obiekt Path i rozwiązał go do rzeczywistej lokalizacji pliku w systemie.
Ścieżka do pliku będzie teraz bezwzględna, jak w `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Zauważ też, że klasa obiektu Path to `sun.nio.fs.UnixPath` — tak Nextflow reprezentuje pliki lokalne.
Jak zobaczymy później, pliki zdalne będą miały inne nazwy klas (np. `nextflow.file.http.XPath` dla plików HTTP), ale wszystkie działają dokładnie tak samo i można ich używać identycznie w workflow'ach.

!!! tip "Wskazówka"

    **Kluczowa różnica:**

    - **Ciąg znaków ścieżki**: Zwykły tekst, który Nextflow traktuje jako znaki
    - **Obiekt Path**: Inteligentne odwołanie do pliku, z którym Nextflow może pracować

    Pomyśl o tym tak: ciąg znaków ścieżki to jak zapisanie adresu na kartce papieru, podczas gdy obiekt Path to jak wczytanie adresu do urządzenia GPS, które wie, jak tam dotrzeć i może powiedzieć Ci szczegóły dotyczące trasy.

### 1.3. Dostęp do atrybutów pliku

Dlaczego to jest pomocne? Teraz, gdy Nextflow rozumie, że `myFile` jest obiektem Path, a nie tylko ciągiem znaków, możemy uzyskać dostęp do różnych atrybutów tego obiektu.

Zaktualizujmy nasz workflow, aby wypisywał wbudowane atrybuty pliku:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Powyżej widzisz różne atrybuty pliku wypisane na konsolę.

### 1.4. Przekazanie pliku do procesu

Różnica między ciągami znaków a obiektami Path staje się krytyczna, gdy zaczynasz budować rzeczywiste workflow'y z procesami.
Do tej pory zweryfikowaliśmy, że Nextflow traktuje teraz nasz plik wejściowy jako plik — sprawdźmy jednak, czy możemy faktycznie uruchomić na nim coś w procesie.

#### 1.4.1. Zaimportuj proces i przejrzyj kod

Dostarczamy Ci gotowy moduł procesu o nazwie `COUNT_LINES`, który przyjmuje plik na wejściu i liczy, ile wierszy się w nim znajduje.

Aby użyć procesu w workflow'ie, wystarczy dodać instrukcję include przed blokiem workflow:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Możesz otworzyć plik modułu, aby przejrzeć jego kod:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Jak widać, to dość prosty skrypt, który rozpakowuje plik i liczy, ile wierszy zawiera.

??? info "Co robi `debug true`?"

    Dyrektywa `debug true` w definicji procesu powoduje, że Nextflow wypisuje wyjście ze skryptu (np. liczbę wierszy „40") bezpośrednio w dzienniku wykonania.
    Bez niej widziałbyś tylko status wykonania procesu, ale nie rzeczywiste wyjście ze skryptu.

    Więcej informacji na temat debugowania procesów Nextflow'a znajdziesz w side queście [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Dodaj wywołanie `COUNT_LINES`

Teraz, gdy proces jest dostępny dla workflow'u, możemy dodać wywołanie procesu `COUNT_LINES`, aby uruchomić go na pliku wejściowym.

Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Policz wiersze w pliku
        COUNT_LINES(myFile)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Teraz uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Widać, że jesteśmy w stanie prawidłowo operować na pliku wewnątrz procesu.

Konkretnie, Nextflow z powodzeniem wykonał następujące operacje:

- Umieścił plik w katalogu roboczym (staging)
- Zdekompresował plik .gz
- Policzył wiersze (40 wierszy w tym przypadku)
- Zakończył bez błędu

Kluczem do tej sprawnej operacji jest to, że jawnie informujemy Nextflow'a, że nasze wejście jest plikiem i powinno być tak traktowane.

### 1.5. Rozwiązywanie podstawowych błędów obsługi plików wejściowych

To często sprawia trudności osobom zaczynającym pracę z Nextflow'em, więc poświęćmy chwilę na przyjrzenie się temu, co się dzieje, gdy zrobimy to źle.

Są dwa główne miejsca, w których można popełnić błąd w obsłudze plików: na poziomie workflow'u i na poziomie procesu.

#### 1.5.1. Błąd na poziomie workflow'u

Zobaczmy, co się stanie, jeśli wrócimy do traktowania pliku jako ciągu znaków przy podawaniu wejścia w bloku workflow.

Wprowadź następujące zmiany w workflow'ie, pamiętając o zakomentowaniu instrukcji wypisujących atrybuty specyficzne dla ścieżki:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Policz wiersze w pliku
        COUNT_LINES(myFile)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Policz wiersze w pliku
        COUNT_LINES(myFile)
    ```

Teraz uruchom workflow:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Oto najważniejszy fragment:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Gdy podajesz wejście typu `path`, Nextflow sprawdza, czy przekazujesz rzeczywiste odwołania do plików, a nie tylko ciągi znaków.
Ten błąd informuje Cię, że `'data/patientA_rep1_normal_R1_001.fastq.gz'` nie jest prawidłową wartością ścieżki, ponieważ jest ciągiem znaków, a nie obiektem Path.

Nextflow natychmiast wykrył problem i zatrzymał się, zanim w ogóle uruchomił proces.

#### 1.5.2. Błąd na poziomie procesu

Innym miejscem, w którym możemy zapomnieć o poinformowaniu Nextflow'a, że chcemy, aby traktował wejście jako plik, jest definicja procesu.

!!! warning "Ostrzeżenie"

    Aby ten test działał poprawnie, pozostaw workflow w stanie z błędem (używając zwykłego ciągu znaków zamiast `file()`).
    W połączeniu z `val` w procesie spowoduje to błąd pokazany poniżej.

Wprowadź następującą zmianę w module:

=== "Po"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Przed"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

Teraz uruchom workflow ponownie:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Wynik zawiera wiele szczegółów dotyczących błędu, ponieważ proces jest skonfigurowany do wypisywania informacji debugowania, jak wspomniano powyżej.

Oto najistotniejsze fragmenty:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

System nie mógł znaleźć pliku, choć jeśli sprawdzisz ścieżkę, plik o tej nazwie istnieje w tej lokalizacji.

Gdy uruchomiliśmy workflow, Nextflow przekazał wartość ciągu znaków do skryptu, ale nie umieścił (_staging_) rzeczywistego pliku w katalogu roboczym.
Proces próbował więc użyć ścieżki względnej `data/patientA_rep1_normal_R1_001.fastq.gz`, ale ten plik nie istnieje w katalogu roboczym procesu.

Oba przykłady razem pokazują, jak ważne jest informowanie Nextflow'a, że wejście powinno być traktowane jako plik.

!!! note "Uwaga"

    Pamiętaj, aby wrócić i naprawić oba celowe błędy przed przejściem do następnej sekcji.

### Podsumowanie

- Ciągi znaków ścieżki vs obiekty Path: ciągi znaków to zwykły tekst, obiekty Path to inteligentne odwołania do plików
- Metoda `file()` konwertuje ciąg znaków ścieżki na obiekt Path, z którym Nextflow może pracować
- Możesz uzyskać dostęp do właściwości pliku, takich jak `name`, `simpleName`, `extension` i `parent`, [używając atrybutów pliku](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Używanie obiektów Path zamiast ciągów znaków pozwala Nextflow'owi prawidłowo zarządzać plikami w Twoim workflow'ie
- Wyniki obsługi wejścia procesu: właściwa obsługa plików wymaga obiektów Path, a nie ciągów znaków, aby zapewnić prawidłowe umieszczenie plików i dostępność dla procesów

---

## 2. Używanie plików zdalnych

Jedną z kluczowych cech Nextflow'a jest możliwość płynnego przełączania się między plikami lokalnymi (na tej samej maszynie) a plikami zdalnymi dostępnymi przez internet.

Jeśli robisz to prawidłowo, nigdy nie powinieneś musieć zmieniać logiki workflow'u, aby obsługiwać pliki pochodzące z różnych lokalizacji.
Aby użyć pliku zdalnego, wystarczy podać odpowiedni prefiks w ścieżce pliku podczas dostarczania go do workflow'u.

Na przykład `/path/to/data` nie ma prefiksu, co wskazuje, że jest to „zwykła" lokalna ścieżka pliku, natomiast `s3://path/to/data` zawiera prefiks `s3://`, wskazując, że plik znajduje się w magazynie obiektów Amazon S3.

Obsługiwanych jest wiele różnych protokołów:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Aby użyć któregokolwiek z nich, wystarczy podać odpowiedni prefiks w ciągu znaków, który technicznie nazywa się wtedy Uniform Resource Identifier (URI), a nie ścieżką pliku.
Nextflow zajmie się uwierzytelnianiem i umieszczaniem plików we właściwym miejscu, pobieraniem lub przesyłaniem oraz wszystkimi innymi oczekiwanymi operacjami na plikach.

Kluczową zaletą tego systemu jest to, że umożliwia przełączanie się między środowiskami bez zmiany logiki pipeline'u.
Na przykład możesz rozwijać workflow z małym, lokalnym zestawem testowym, a następnie przełączyć się na pełnowymiarowy zestaw testowy w zdalnym magazynie, po prostu zmieniając URI.

### 2.1. Użycie pliku z internetu

Przetestujmy to, zastępując lokalną ścieżkę w naszym workflow'ie ścieżką HTTPS wskazującą na kopię tych samych danych przechowywaną w Github.

!!! warning "Ostrzeżenie"

    Zadziała to tylko wtedy, gdy masz aktywne połączenie z internetem.

Otwórz ponownie `main.nf` i zmień ścieżkę wejściową w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Użyj pliku zdalnego z internetu
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Uruchommy workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Działa! Widać, że niewiele się zmieniło.

Jedyna różnica w wyjściu konsoli polega na tym, że klasa obiektu Path to teraz `nextflow.file.http.XPath`, podczas gdy dla ścieżki lokalnej była to `sun.nio.fs.UnixPath`.
Nie musisz pamiętać tych klas — wspominamy o tym tylko po to, aby pokazać, że Nextflow identyfikuje i obsługuje różne lokalizacje w odpowiedni sposób.

Za kulisami Nextflow pobrał plik do katalogu staging znajdującego się w katalogu roboczym.
Ten umieszczony plik może być następnie traktowany jako plik lokalny i dowiązany symbolicznie do odpowiedniego katalogu procesu.

Możesz to zweryfikować, przeglądając zawartość katalogu roboczego znajdującego się pod wartością skrótu procesu.

??? abstract "Zawartość katalogu roboczego"

    Jeśli skrót procesu to `8a/2ab7ca`, możesz przejrzeć katalog roboczy:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Dowiązanie symboliczne wskazuje na umieszczoną kopię pliku zdalnego, którą Nextflow pobrał automatycznie.

Pamiętaj, że w przypadku większych plików krok pobierania zajmie więcej czasu niż uruchamianie na plikach lokalnych.
Nextflow sprawdza jednak, czy ma już umieszczoną kopię, aby uniknąć niepotrzebnych pobrań.
Jeśli więc uruchomisz workflow ponownie na tym samym pliku i nie usunąłeś umieszczonej kopii, Nextflow użyje jej.

To pokazuje, jak łatwo jest przełączać się między danymi lokalnymi i zdalnymi w Nextflow'ie, co jest kluczową cechą tego narzędzia.

!!! note "Uwaga"

    Jednym ważnym wyjątkiem od tej zasady jest to, że nie można używać wzorców glob ani ścieżek katalogów z HTTPS, ponieważ HTTPS nie może wylistować wielu plików — musisz podać dokładne adresy URL plików.
    Inne protokoły przechowywania, takie jak blob storage (`s3://`, `az://`, `gs://`), obsługują zarówno wzorce glob, jak i ścieżki katalogów.

    Oto jak można używać wzorców glob z magazynem w chmurze:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 ze wzorcami glob - dopasuje wiele plików
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage ze wzorcami glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage ze wzorcami glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Pokażemy Ci, jak pracować z wzorcami glob w praktyce w następnej sekcji.

### 2.2. Powrót do pliku lokalnego

Wrócimy teraz do używania lokalnych przykładowych plików w pozostałej części tego side questa, więc przełączmy wejście workflow'u z powrotem na oryginalny plik:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Podsumowanie

- Dane zdalne są dostępne za pomocą URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow automatycznie pobierze i umieści dane we właściwym miejscu, o ile ścieżki są przekazywane do procesów
- Nie pisz logiki do pobierania ani przesyłania plików zdalnych!
- Pliki lokalne i zdalne tworzą różne typy obiektów, ale działają identycznie
- **Ważne**: HTTP/HTTPS działa tylko z pojedynczymi plikami (bez wzorców glob)
- Magazyn w chmurze (S3, Azure, GCS) obsługuje zarówno pojedyncze pliki, jak i wzorce glob
- Możesz płynnie przełączać się między lokalnymi i zdalnymi źródłami danych bez zmiany logiki kodu (o ile protokół obsługuje wymagane operacje)

---

## 3. Używanie fabryki kanałów `fromPath()`

Do tej pory pracowaliśmy z jednym plikiem na raz, ale w Nextflow'ie zazwyczaj będziemy chcieli tworzyć kanał wejściowy z wieloma plikami do przetworzenia.

Naiwnym sposobem na to byłoby połączenie metody `file()` z [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) w taki sposób:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

To działa, ale jest nieporęczne.

!!! tip "Wskazówka: kiedy używać `file()` vs `channel.fromPath()`"

    - Używaj `file()`, gdy potrzebujesz pojedynczego obiektu Path do bezpośredniej manipulacji (sprawdzanie, czy plik istnieje, odczytywanie jego atrybutów lub przekazywanie do pojedynczego wywołania procesu)
    - Używaj `channel.fromPath()`, gdy potrzebujesz kanału mogącego przechowywać wiele plików, szczególnie ze wzorcami glob, lub gdy pliki będą przepływać przez wiele procesów

Tu właśnie przydaje się [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): wygodna fabryka kanałów, która łączy całą potrzebną funkcjonalność do generowania kanału z jednego lub więcej statycznych ciągów znaków pliku, a także wzorców glob.

### 3.1. Dodanie fabryki kanałów

Zaktualizujmy nasz workflow, aby używał `channel.fromPath`.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Wypisz atrybuty pliku
        /* Na razie zakomentuj, wrócimy do tego!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Policz wiersze w pliku
        // COUNT_LINES(myFile)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Utwórz obiekt Path z ciągu znaków reprezentującego ścieżkę
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Policz wiersze w pliku
        COUNT_LINES(myFile)
    ```

Zakomentowaliśmy też kod wypisujący atrybuty i dodaliśmy instrukcję `.view`, aby wypisywać tylko nazwę pliku.

Uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Jak widać, ścieżka pliku jest wczytywana jako obiekt typu `Path` w kanale.
Jest to podobne do tego, co zrobiłaby metoda `file()`, z tą różnicą, że mamy teraz kanał, do którego możemy wczytać więcej plików.

Użycie `channel.fromPath()` to wygodny sposób tworzenia nowego kanału wypełnionego listą plików.

### 3.2. Wyświetlanie atrybutów plików w kanale

W pierwszym podejściu do używania fabryki kanałów uprościliśmy kod i wypisaliśmy tylko nazwę pliku.

Wróćmy do wypisywania pełnych atrybutów pliku:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Policz wiersze w pliku
        COUNT_LINES(ch_files)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Policz wiersze w pliku
        // COUNT_LINES(ch_files)
    ```

Ponownie włączamy też wywołanie procesu `COUNT_LINES`, aby sprawdzić, czy przetwarzanie plików nadal działa poprawnie z naszym podejściem opartym na kanałach.

Uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Takie same wyniki jak poprzednio, ale teraz plik jest w kanale, więc możemy dodać więcej.

### 3.3. Używanie wzorca glob do dopasowania wielu plików

Istnieje kilka sposobów wczytania większej liczby plików do kanału.
Pokażemy Ci, jak używać wzorców glob, które są wygodnym sposobem dopasowywania i pobierania nazw plików i katalogów na podstawie znaków wieloznacznych.
Proces dopasowywania tych wzorców nazywa się „globbingiem" lub „rozwijaniem nazw plików".

!!! note "Uwaga"

    Jak wspomniano wcześniej, Nextflow obsługuje globbing do zarządzania plikami wejściowymi i wyjściowymi w większości przypadków, z wyjątkiem ścieżek HTTPS, ponieważ HTTPS nie może wylistować wielu plików.

Powiedzmy, że chcemy pobrać oba pliki z pary plików powiązanych z danym pacjentem, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Ponieważ jedyną różnicą między nazwami plików jest numer repliki, tj. liczba po `R`, możemy użyć znaku wieloznacznego `*` jako zastępnika dla tej liczby:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

To jest wzorzec glob, którego potrzebujemy.

Teraz wystarczy zaktualizować ścieżkę pliku w fabryce kanałów, aby używała tego wzorca:

=== "Po"

    ```groovy title="main.nf" linenums="7"
      // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7"
      // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow automatycznie rozpozna, że jest to wzorzec glob i odpowiednio go obsłuży.

Uruchom workflow, aby to przetestować:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Jak widać, mamy teraz dwa obiekty Path w naszym kanale, co pokazuje, że Nextflow poprawnie rozwinął nazwy plików i wczytał oraz przetworzył oba pliki zgodnie z oczekiwaniami.

Używając tej metody, możemy pobierać dowolną liczbę plików, po prostu zmieniając wzorzec glob. Gdybyśmy uczynili go bardziej ogólnym, na przykład zastępując wszystkie zmienne części nazw plików przez `*` (np. `data/patient*_rep*_*_R*_001.fastq.gz`), moglibyśmy pobrać wszystkie przykładowe pliki z katalogu `data`.

### Podsumowanie

- `channel.fromPath()` tworzy kanał z plikami pasującymi do wzorca
- Każdy plik jest emitowany jako osobny element w kanale
- Możemy używać wzorca glob do dopasowania wielu plików
- Pliki są automatycznie konwertowane na obiekty Path z pełnymi atrybutami
- Metoda `.view()` umożliwia inspekcję zawartości kanału

---

## 4. Wyodrębnianie podstawowych metadanych z nazw plików

W większości dziedzin naukowych bardzo powszechne jest kodowanie metadanych w nazwach plików zawierających dane.
Na przykład w bioinformatyce pliki zawierające dane sekwencjonowania są często nazywane w sposób kodujący informacje o próbce, warunku, replice i numerze odczytu.

Jeśli nazwy plików są skonstruowane zgodnie ze spójną konwencją, możesz wyodrębniać te metadane w ustandaryzowany sposób i używać ich w toku analizy.
To duże „jeśli", oczywiście, i powinieneś być bardzo ostrożny, gdy polegasz na strukturze nazw plików; ale w rzeczywistości podejście to jest bardzo szeroko stosowane, więc przyjrzyjmy się, jak to się robi w Nextflow'ie.

W przypadku naszych przykładowych danych wiemy, że nazwy plików zawierają spójnie ustrukturyzowane metadane.
Na przykład nazwa pliku `patientA_rep1_normal_R2_001` koduje następujące informacje:

- ID pacjenta: `patientA`
- ID repliki: `rep1`
- typ próbki: `normal` (w przeciwieństwie do `tumor`)
- zestaw odczytów: `R1` (w przeciwieństwie do `R2`)

Zmodyfikujemy nasz workflow, aby pobierał te informacje w trzech krokach:

1. Pobranie `simpleName` pliku, który zawiera metadane
2. Rozdzielenie metadanych za pomocą metody `tokenize()`
3. Użycie mapy do organizacji metadanych

!!! warning "Ostrzeżenie"

    Nigdy nie powinieneś kodować wrażliwych informacji w nazwach plików, takich jak imiona pacjentów lub inne cechy identyfikujące, ponieważ może to naruszyć prywatność pacjentów lub inne stosowne ograniczenia bezpieczeństwa.

### 4.1. Pobranie `simpleName`

`simpleName` to atrybut pliku odpowiadający nazwie pliku pozbawionej ścieżki i rozszerzenia.

Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Pobiera to `simpleName` i kojarzy go z pełnym obiektem pliku za pomocą operacji `map()`.

Uruchom workflow, aby sprawdzić, czy działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Każdy element w kanale jest teraz krotką zawierającą `simpleName` i oryginalny obiekt pliku.

### 4.2. Wyodrębnienie metadanych z `simpleName`

W tym momencie metadane, których szukamy, są osadzone w `simpleName`, ale nie możemy bezpośrednio uzyskać dostępu do poszczególnych elementów.
Musimy więc podzielić `simpleName` na jego składniki.
Na szczęście te składniki są po prostu oddzielone podkreśleniami w oryginalnej nazwie pliku, więc możemy zastosować powszechną metodę Nextflow'a o nazwie `tokenize()`, która jest idealna do tego zadania.

Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Metoda `tokenize()` podzieli ciąg znaków `simpleName` wszędzie tam, gdzie znajdzie podkreślenia, i zwróci listę zawierającą podciągi.

Uruchom workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Teraz krotka dla każdego elementu w naszym kanale zawiera listę metadanych (np. `[patientA, rep1, normal, R1, 001]`) i oryginalny obiekt pliku.

Świetnie!
Rozłożyliśmy informacje o pacjencie z pojedynczego ciągu znaków na listę ciągów.
Możemy teraz obsługiwać każdą część informacji o pacjencie osobno.

### 4.3. Użycie mapy do organizacji metadanych

Nasze metadane to na razie zwykła lista.
Jest wystarczająco łatwa w użyciu, ale trudna do odczytania.

```console
[patientA, rep1, normal, R1, 001]
```

Czym jest element o indeksie 3? Czy możesz to stwierdzić bez odwoływania się do oryginalnego opisu struktury metadanych?

To doskonała okazja do użycia magazynu klucz-wartość, gdzie każdy element ma zestaw kluczy i powiązanych z nimi wartości, dzięki czemu można łatwo odwoływać się do każdego klucza, aby uzyskać odpowiednią wartość.

W naszym przykładzie oznacza to przejście od tej organizacji:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Do tej:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

W Nextflow'ie nazywa się to [mapą](https://nextflow.io/docs/latest/script.html#maps).

Przekonwertujmy teraz naszą płaską listę na mapę.
Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Wczytaj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Kluczowe zmiany to:

- **Przypisanie destrukturyzujące**: `def (patient, replicate, type, readNum) = ...` wyodrębnia stokenizowane wartości do nazwanych zmiennych w jednej linii
- **Składnia literału mapy**: `[id: patient, replicate: ...]` tworzy mapę, w której każdy klucz (np. `id`) jest powiązany z wartością (np. `patient`)
- **Zagnieżdżona struktura**: Zewnętrzna lista `[..., myFile]` łączy mapę metadanych z oryginalnym obiektem pliku

Uprościliśmy też kilka ciągów metadanych za pomocą metody zastępowania ciągów `replace()`, aby usunąć zbędne znaki (np. `replicate.replace('rep', '')`, aby zachować tylko numer z ID replik).

Uruchommy workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Teraz metadane są starannie oznaczone (np. `[id:patientA, replicate:1, type:normal, readNum:2]`), więc znacznie łatwiej jest stwierdzić, co jest czym.

Będzie też znacznie łatwiej faktycznie korzystać z elementów metadanych w workflow'ie, a nasz kod stanie się łatwiejszy do odczytania i utrzymania.

### Podsumowanie

- Możemy obsługiwać nazwy plików w Nextflow'ie z mocą pełnego języka programowania
- Możemy traktować nazwy plików jako ciągi znaków, aby wyodrębniać istotne informacje
- Użycie metod takich jak `tokenize()` i `replace()` pozwala nam manipulować ciągami znaków w nazwie pliku
- Operacja `.map()` przekształca elementy kanału, zachowując strukturę
- Ustrukturyzowane metadane (mapy) sprawiają, że kod jest bardziej czytelny i łatwiejszy w utrzymaniu niż listy pozycyjne

Następnie przyjrzymy się, jak obsługiwać sparowane pliki danych.

---

## 5. Obsługa sparowanych plików danych

Wiele projektów eksperymentalnych produkuje sparowane pliki danych, które korzystają z jawnie sparowanej obsługi.
Na przykład w bioinformatyce dane sekwencjonowania są często generowane w postaci sparowanych odczytów, czyli ciągów sekwencji pochodzących z tego samego fragmentu DNA (często nazywanych „w przód" i „wstecz", ponieważ są odczytywane z przeciwnych końców).

Tak jest w przypadku naszych przykładowych danych, gdzie R1 i R2 odnoszą się do dwóch zestawów odczytów.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow dostarcza wyspecjalizowaną fabrykę kanałów do pracy z takimi sparowanymi plikami o nazwie `channel.fromFilePairs()`, która automatycznie grupuje pliki na podstawie wspólnego wzorca nazewnictwa. Pozwala to na ściślejsze powiązanie sparowanych plików przy mniejszym wysiłku.

Zmodyfikujemy nasz workflow, aby to wykorzystać.
Zajmie to dwa kroki:

1. Przełączenie fabryki kanałów na `channel.fromFilePairs()`
2. Wyodrębnienie i mapowanie metadanych

### 5.1. Przełączenie fabryki kanałów na `channel.fromFilePairs()`

Aby użyć `channel.fromFilePairs`, musimy podać wzorzec, którego Nextflow powinien używać do identyfikowania dwóch członków pary.

Wracając do naszych przykładowych danych, możemy sformalizować wzorzec nazewnictwa w następujący sposób:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Jest to podobne do wzorca glob, którego używaliśmy wcześniej, z tą różnicą, że ten konkretnie wylicza podciągi (albo `1`, albo `2` następujące bezpośrednio po R), które identyfikują dwóch członków pary.

Zaktualizujmy odpowiednio workflow `main.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Na razie zakomentuj mapowanie, wrócimy do tego!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Przełączyliśmy fabrykę kanałów i dostosowaliśmy wzorzec dopasowywania plików, a przy okazji zakomentowaliśmy operację map.
Dodamy ją z powrotem później, z kilkoma modyfikacjami.

Uruchom workflow, aby go przetestować:

```bash
nextflow run main.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Ups, tym razem uruchomienie się nie powiodło!

Istotna część komunikatu o błędzie to:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Dzieje się tak, ponieważ zmieniliśmy fabrykę kanałów.
Do tej pory oryginalny kanał wejściowy zawierał tylko ścieżki plików.
Wszystkie manipulacje metadanymi, które wykonywaliśmy, w rzeczywistości nie wpływały na zawartość kanału.

Teraz, gdy używamy fabryki kanałów `.fromFilePairs`, zawartość wynikowego kanału jest inna.
Widzimy tylko jeden element kanału, złożony z krotki zawierającej dwa elementy: część `simpleName` wspólną dla obu plików, która służy jako identyfikator, oraz krotkę zawierającą dwa obiekty pliku, w formacie `id, [ file1, file2 ]`.

To świetnie, ponieważ Nextflow wykonał ciężką pracę wyodrębnienia nazwy pacjenta, badając wspólny prefiks i używając go jako identyfikatora pacjenta.

Jednak psuje to nasz obecny workflow.
Gdybyśmy chcieli nadal uruchamiać `COUNT_LINES` w ten sam sposób bez zmiany procesu, musielibyśmy zastosować operację mapowania, aby wyodrębnić ścieżki plików.
Nie zamierzamy tego robić, ponieważ naszym ostatecznym celem jest użycie innego procesu, `ANALYZE_READS`, który odpowiednio obsługuje pary plików.

Zakomentujmy więc (lub usuńmy) wywołanie `COUNT_LINES` i przejdźmy dalej.

=== "Po"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Policz wiersze w pliku
        // COUNT_LINES(ch_files)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Policz wiersze w pliku
        COUNT_LINES(ch_files)
    ```

Możesz też zakomentować lub usunąć instrukcję include dla `COUNT_LINES`, ale nie będzie to miało żadnego efektu funkcjonalnego.

Teraz uruchommy workflow ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Tym razem workflow zakończył się sukcesem!

Musimy jednak jeszcze wyodrębnić pozostałe metadane z pola `id`.

### 5.2. Wyodrębnienie i organizacja metadanych z par plików

Nasza poprzednia operacja `map` nie zadziała, ponieważ nie pasuje do struktury danych, ale możemy ją zmodyfikować.

Mamy już dostęp do rzeczywistego identyfikatora pacjenta w ciągu znaków, którego `fromFilePairs()` użył jako identyfikatora, więc możemy go użyć do wyodrębnienia metadanych bez pobierania `simpleName` z obiektu Path jak poprzednio.

Odkomentuj operację map w workflow'ie i wprowadź następujące zmiany:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Na razie zakomentuj mapowanie, wrócimy do tego!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Tym razem mapa zaczyna się od `id, files` zamiast tylko `myFile`, a `tokenize()` jest stosowane do `id` zamiast do `myFile.simpleName`.

Zauważ też, że usunęliśmy `readNum` z linii `tokenize()`; wszelkie podciągi, których konkretnie nie nazwiemy (zaczynając od lewej), zostaną po cichu pominięte.
Możemy to zrobić, ponieważ sparowane pliki są teraz ściśle powiązane, więc nie potrzebujemy już `readNum` w mapie metadanych.

Uruchommy workflow:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

I mamy to: mapa metadanych (`[id:patientA, replicate:1, type:normal]`) na pierwszej pozycji krotki wyjściowej, a po niej krotka sparowanych plików, zgodnie z zamierzeniem.

Oczywiście, to pobierze i przetworzy tylko tę konkretną parę plików.
Jeśli chcesz poeksperymentować z przetwarzaniem wielu par, możesz spróbować dodać znaki wieloznaczne do wzorca wejściowego i zobaczyć, co się stanie.
Na przykład spróbuj użyć `data/patientA_rep1_*_R{1,2}_001.fastq.gz`.

### Podsumowanie

- [`channel.fromFilePairs()` automatycznie znajduje i paruje powiązane pliki](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Upraszcza to obsługę odczytów paired-end w Twoim pipeline'ie
- Sparowane pliki mogą być grupowane jako krotki `[id, [file1, file2]]`
- Wyodrębnianie metadanych można wykonać z ID sparowanego pliku, a nie z poszczególnych plików

---

## 6. Używanie operacji na plikach w procesach

Teraz złóżmy to wszystko razem w prostym procesie, aby utrwalić sposób używania operacji na plikach wewnątrz procesu Nextflow'a.

Dostarczamy Ci gotowy moduł procesu o nazwie `ANALYZE_READS`, który przyjmuje krotkę metadanych i parę plików wejściowych i je analizuje.
Możemy sobie wyobrazić, że wykonuje dopasowanie sekwencji, wykrywanie wariantów lub dowolny inny krok, który ma sens dla tego typu danych.

Zaczynajmy.

### 6.1. Zaimportuj proces i przejrzyj kod

Aby użyć tego procesu w workflow'ie, wystarczy dodać instrukcję include modułu przed blokiem workflow.

Wprowadź następującą zmianę w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Możesz otworzyć plik modułu, aby przejrzeć jego kod:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note "Uwaga"

    Dyrektywy `tag` i `publishDir` używają składni closure (`{ ... }`) zamiast interpolacji ciągów znaków (`"${...}"`).
    Dzieje się tak, ponieważ te dyrektywy odwołują się do zmiennych wejściowych (`meta`), które nie są dostępne aż do czasu wykonania.
    Składnia closure odracza ewaluację do momentu faktycznego uruchomienia procesu.

!!! note "Uwaga"

    Naszą mapę metadanych nazywamy `meta` zgodnie z konwencją.
    Aby dowiedzieć się więcej o mapach meta, zapoznaj się z side questem [Metadata and meta maps](../metadata/).

### 6.2. Wywołanie procesu w workflow'ie

Teraz, gdy proces jest dostępny dla workflow'u, możemy dodać wywołanie procesu `ANALYZE_READS`, aby go uruchomić.

Aby uruchomić go na naszych przykładowych danych, musimy zrobić dwie rzeczy:

1. Nadać nazwę przemapowanemu kanałowi
2. Dodać wywołanie procesu

#### 6.2.1. Nadanie nazwy przemapowanemu kanałowi wejściowemu

Wcześniej stosowaliśmy manipulacje mapowaniem bezpośrednio do kanału wejściowego.
Aby przekazać przemapowaną zawartość do procesu `ANALYZE_READS` (i zrobić to w sposób przejrzysty i łatwy do odczytania), chcemy utworzyć nowy kanał o nazwie `ch_samples`.

Możemy to zrobić za pomocą operatora [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

W głównym workflow'ie zastąp operator `.view()` przez `.set { ch_samples }` i dodaj linię testującą, czy możemy odwoływać się do kanału po nazwie.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Tymczasowo: podejrzyj zawartość ch_samples
        ch_samples.view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Uruchommy to:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Potwierdza to, że możemy teraz odwoływać się do kanału po nazwie.

#### 6.2.2. Wywołanie procesu na danych

Teraz faktycznie wywołajmy proces `ANALYZE_READS` na kanale `ch_samples`.

W głównym workflow'ie wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="main.nf" linenums="23"
        // Uruchom analizę
        ANALYZE_READS(ch_samples)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="23"
        // Tymczasowo: podejrzyj zawartość ch_samples
        ch_samples.view()
    ```

Uruchommy to:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Ten proces jest skonfigurowany do publikowania wyjść w katalogu `results`, więc zajrzyj tam.

??? abstract "Zawartość katalogu i pliku"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

Proces przyjął nasze wejścia i zgodnie z zamierzeniem utworzył nowy plik zawierający metadane pacjenta.
Znakomicie!

### 6.3. Przetwarzanie wielu pacjentów

Oczywiście, to przetwarza tylko jedną parę plików dla jednego pacjenta, co nie jest dokładnie tym rodzajem wysokiej przepustowości, na którą liczysz w Nextflow'ie.
Prawdopodobnie będziesz chciał przetwarzać znacznie więcej danych na raz.

Pamiętaj, że `channel.fromPath()` przyjmuje _glob_ jako wejście, co oznacza, że może przyjąć dowolną liczbę plików pasujących do wzorca.
Jeśli więc chcemy uwzględnić wszystkich pacjentów, możemy po prostu zmodyfikować ciąg wejściowy, aby obejmował więcej pacjentów, jak wspomniano wcześniej.

Udawajmy, że chcemy być jak najbardziej zachłanni.
Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Wczytaj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Uruchom pipeline ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Katalog wyników powinien teraz zawierać wyniki dla wszystkich dostępnych danych.

??? abstract "Zawartość katalogu"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Sukces! Przeanalizowaliśmy wszystkich pacjentów za jednym razem! Prawda?

Może nie do końca.
Jeśli przyjrzysz się bliżej, mamy problem: mamy dwie repliki dla pacjenta A, ale tylko jeden plik wyjściowy!
Za każdym razem nadpisujemy plik wyjściowy.

### 6.4. Zapewnienie unikalności publikowanych plików

Ponieważ mamy dostęp do metadanych pacjenta, możemy ich użyć, aby zapewnić unikalność publikowanych plików, włączając różnicujące metadane — albo w strukturze katalogów, albo w samych nazwach plików.

Wprowadź następującą zmianę w workflow'ie:

=== "Po"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Przed"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Pokazujemy tu opcję używania dodatkowych poziomów katalogów, aby uwzględnić typy próbek i repliki, ale możesz też poeksperymentować z robieniem tego na poziomie nazwy pliku.

Teraz uruchom pipeline jeszcze raz, ale najpierw usuń katalog wyników, aby mieć czyste środowisko pracy:

```bash
rm -r results
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Sprawdź teraz katalog wyników:

??? abstract "Zawartość katalogu"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

I mamy to — wszystkie nasze metadane, starannie zorganizowane. To sukces!

Jest wiele więcej rzeczy, które możesz zrobić, gdy masz metadane wczytane do mapy:

1. Tworzyć zorganizowane katalogi wyjściowe na podstawie atrybutów pacjenta
2. Podejmować decyzje w procesach na podstawie właściwości pacjenta
3. Dzielić, łączyć i ponownie łączyć dane na podstawie wartości metadanych

Ten wzorzec jawnego przechowywania metadanych i dołączania ich do danych (zamiast kodowania w nazwach plików) jest podstawową najlepszą praktyką w Nextflow'ie, która umożliwia budowanie solidnych i łatwych w utrzymaniu workflow'ów analitycznych.
Więcej na ten temat możesz dowiedzieć się w side queście [Metadata and meta maps](../metadata/).

### Podsumowanie

- Dyrektywa `publishDir` może organizować wyjścia na podstawie wartości metadanych
- Metadane w krotkach umożliwiają ustrukturyzowaną organizację wyników
- To podejście tworzy łatwe w utrzymaniu workflow'y z przejrzystą proweniencją danych
- Procesy mogą przyjmować krotki metadanych i plików jako wejście
- Dyrektywa `tag` zapewnia identyfikację procesu w dziennikach wykonania
- Struktura workflow'u oddziela tworzenie kanałów od wykonywania procesów

---

## Podsumowanie

W tym side queście nauczyłeś się pracować z plikami w Nextflow'ie — od podstawowych operacji po bardziej zaawansowane techniki obsługi kolekcji plików.

Stosowanie tych technik we własnej pracy pozwoli Ci budować wydajniejsze i łatwiejsze w utrzymaniu workflow'y, szczególnie przy pracy z dużą liczbą plików o złożonych konwencjach nazewnictwa.

### Kluczowe wzorce

1.  **Podstawowe operacje na plikach:** Tworzyliśmy obiekty Path za pomocą `file()` i uzyskiwaliśmy dostęp do atrybutów pliku, takich jak nazwa, rozszerzenie i katalog nadrzędny, ucząc się różnicy między ciągami znaków a obiektami Path.

    - Utwórz obiekt Path za pomocą `file()`

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Pobierz atrybuty pliku

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Używanie plików zdalnych**: Nauczyliśmy się, jak przezroczyście przełączać się między plikami lokalnymi i zdalnymi za pomocą URI, demonstrując zdolność Nextflow'a do obsługi plików z różnych źródeł bez zmiany logiki workflow'u.

    - Plik lokalny

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **Wczytywanie plików za pomocą fabryki kanałów `fromPath()`:** Tworzyliśmy kanały ze wzorców plików za pomocą `channel.fromPath()` i przeglądaliśmy ich atrybuty, w tym typy obiektów.

    - Utwórz kanał ze wzorca pliku

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Pobierz atrybuty pliku

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Wyodrębnianie metadanych pacjenta z nazw plików:** Używaliśmy `tokenize()` i `replace()` do wyodrębniania i strukturyzowania metadanych z nazw plików, konwertując je na zorganizowane mapy.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Uproszczenie za pomocą channel.fromFilePairs:** Używaliśmy `channel.fromFilePairs()` do automatycznego parowania powiązanych plików i wyodrębniania metadanych z ID sparowanych plików.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Używanie operacji na plikach w procesach:** Integrowaliśmy operacje na plikach z procesami Nextflow'a z właściwą obsługą wejścia, używając `publishDir` do organizowania wyjść na podstawie metadanych.

    - Powiąż mapę meta z wejściami procesu

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Organizuj wyjścia na podstawie metadanych

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Dodatkowe zasoby

- [Dokumentacja Nextflow'a: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Co dalej?

Wróć do [menu Side Quests](../) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
