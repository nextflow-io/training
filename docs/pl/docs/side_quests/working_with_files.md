# Przetwarzanie plików wejściowych

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Workflow'y analizy naukowej często wymagają przetwarzania dużej liczby plików.
Nextflow dostarcza potężne narzędzia do efektywnej obsługi plików, pomagając Ci organizować i przetwarzać dane przy użyciu minimalnej ilości kodu.

### Cele szkolenia

W tym zadaniu pobocznym zbadamy, jak Nextflow obsługuje pliki, od podstawowych operacji na plikach po bardziej zaawansowane techniki pracy z kolekcjami plików.
Nauczysz się, jak wydobywać metadane z nazw plików, co jest powszechnym wymogiem w pipeline'ach analizy naukowej.

Pod koniec tego zadania pobocznego będziesz potrafić:

- Tworzyć obiekty Path z ciągów znaków ścieżek plików przy użyciu metody `file()` w Nextflow
- Uzyskiwać dostęp do atrybutów plików, takich jak nazwa, rozszerzenie i katalog nadrzędny
- Obsługiwać zarówno pliki lokalne, jak i zdalne w sposób przezroczysty przy użyciu URI
- Używać kanałów do automatyzacji obsługi plików za pomocą `channel.fromPath()` i `channel.fromFilePairs()`
- Wydobywać i strukturyzować metadane z nazw plików przy użyciu manipulacji ciągami znaków
- Grupować powiązane pliki przy użyciu dopasowywania wzorców i wyrażeń glob
- Integrować operacje na plikach z procesami Nextflow przy odpowiedniej obsłudze wejść
- Organizować wyjścia procesów przy użyciu struktur katalogów opartych na metadanych

Te umiejętności pomogą Ci budować workflow'y, które mogą obsługiwać różne rodzaje wejść plikowych z dużą elastycznością.

### Wymagania wstępne

Przed podjęciem tego zadania pobocznego powinieneś:

- Ukończyć samouczek [Hello Nextflow](../../hello_nextflow/) lub równoważny kurs dla początkujących
- Czuć się komfortowo z podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory)

---

## 0. Rozpocznij pracę

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzyłeś środowisko szkoleniowe zgodnie z opisem w [Konfiguracji środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki dla tego samouczka.

```bash
cd side-quests/working_with_files
```

Możesz ustawić VSCode, aby skupił się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz prosty plik workflow'a o nazwie `main.nf`, katalog `modules` zawierający dwa pliki modułów oraz katalog `data` zawierający przykładowe pliki danych.

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

Ten katalog zawiera dane sekwencjonowania sparowanego z trzech pacjentów (A, B, C).

Dla każdego pacjenta mamy próbki typu `tumor` (zazwyczaj pochodzące z biopsji guza) lub `normal` (pobrane ze zdrowej tkanki lub krwi).
Jeśli nie znasz analizy nowotworów, wiedz tylko, że odpowiada to modelowi eksperymentalnemu, który wykorzystuje sparowane próbki guz/normalny do przeprowadzania analiz kontrastowych.

Dla pacjenta A mamy dwa zestawy replik technicznych (powtórzeń).

Pliki danych sekwencjonowania są nazwane zgodnie z typową konwencją `_R1_` i `_R2_` dla tzw. 'odczytów do przodu' i 'odczytów do tyłu'.

_Nie martw się, jeśli nie znasz tego projektu eksperymentalnego, nie jest to krytyczne dla zrozumienia tego samouczka._

#### Przejrzyj zadanie

Twoim wyzwaniem jest napisanie workflow'a Nextflow, który:

1. **Załaduje** pliki wejściowe przy użyciu metod obsługi plików Nextflow
2. **Wydobędzie** metadane (ID pacjenta, replika, typ próbki) ze struktury nazwy pliku
3. **Zgrupuje** sparowane pliki (R1/R2) razem przy użyciu `channel.fromFilePairs()`
4. **Przetworzy** pliki za pomocą dostarczonego modułu analizy
5. **Zorganizuje** wyjścia w strukturę katalogów opartą na wydobytych metadanych

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moja przestrzeń kodowa działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Podstawowe operacje na plikach

### 1.1. Zidentyfikuj typ obiektu za pomocą `.class`

Spójrz na plik workflow'a `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Utwórz obiekt Path z ciągu znaków ścieżki
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

To mini-workflow (bez żadnych procesów), który odwołuje się do pojedynczej ścieżki pliku w swoim workflow'ie, a następnie wypisuje ją do konsoli wraz z jej klasą.

??? info "Co to jest `.class`?"

    W Nextflow `.class` mówi nam, z jakim typem obiektu pracujemy. To jak pytanie "co to za rzecz?", aby dowiedzieć się, czy jest to ciąg znaków, liczba, plik czy coś innego.
    Pomoże nam to zilustrować różnicę między zwykłym ciągiem znaków a obiektem Path w następnych sekcjach.

Uruchommy workflow'a:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Jak widzisz, Nextflow wypisał ścieżkę jako ciąg znaków dokładnie tak, jak ją napisaliśmy.

To tylko wyjście tekstowe; Nextflow nie zrobił jeszcze z tym nic specjalnego.
Potwierdziliśmy również, że dla Nextflow'a jest to tylko ciąg znaków (klasy `java.lang.String`).
Ma to sens, ponieważ nie powiedzieliśmy jeszcze Nextflow'owi, że odpowiada to plikowi.

### 1.2. Utwórz obiekt Path za pomocą file()

Możemy powiedzieć Nextflow'owi, jak obsługiwać pliki, tworząc [obiekty Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) z ciągów znaków ścieżek.

W naszym workflow'ie możemy przekonwertować ciąg znaków ścieżki `data/patientA_rep1_normal_R1_001.fastq.gz` na obiekt Path przy użyciu metody `file()`, która zapewnia dostęp do właściwości i operacji na plikach.

Edytuj `main.nf`, aby opakować ciąg znaków za pomocą `file()` w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Teraz uruchom workflow'a ponownie:

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
Ścieżka pliku będzie teraz bezwzględna, jak w `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Zauważ również, że klasa obiektu Path to `sun.nio.fs.UnixPath`: to sposób Nextflow'a na reprezentowanie plików lokalnych.
Jak zobaczymy później, pliki zdalne będą miały różne nazwy klas (takie jak `nextflow.file.http.XPath` dla plików HTTP), ale wszystkie działają dokładnie w ten sam sposób i mogą być używane identycznie w Twoich workflow'ach.

!!! tip "Wskazówka"

    **Kluczowa różnica:**

    - **Ciąg znaków ścieżki**: Tylko tekst, który Nextflow traktuje jako znaki
    - **Obiekt Path**: Inteligentne odniesienie do pliku, z którym Nextflow może pracować

    Pomyśl o tym w ten sposób: ciąg znaków ścieżki jest jak napisanie adresu na papierze, podczas gdy obiekt Path jest jak załadowanie adresu w urządzeniu GPS, które wie, jak tam nawigować i może powiedzieć Ci szczegóły o podróży.

### 1.3. Uzyskaj dostęp do atrybutów pliku

Dlaczego jest to pomocne? Cóż, teraz gdy Nextflow rozumie, że `myFile` jest obiektem Path, a nie tylko ciągiem znaków, możemy uzyskać dostęp do różnych atrybutów obiektu Path.

Zaktualizujmy nasz workflow'a, aby wypisać wbudowane atrybuty pliku:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Utwórz obiekt Path z ciągu znaków ścieżki
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
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Uruchom workflow'a:

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

Widzisz różne atrybuty pliku wypisane do konsoli powyżej.

### 1.4. Przekaż plik do procesu

Różnica między ciągami znaków a obiektami Path staje się krytyczna, gdy zaczynasz budować rzeczywiste workflow'y z procesami.
Do tej pory zweryfikowaliśmy, że Nextflow traktuje teraz nasz plik wejściowy jako plik, ale zobaczmy, czy możemy faktycznie uruchomić coś na tym pliku w procesie.

#### 1.4.1. Zaimportuj proces i zbadaj kod

Dostarczamy Ci wstępnie napisany moduł procesu o nazwie `COUNT_LINES`, który przyjmuje plik wejściowy i liczy, ile linii zawiera.

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

Możesz otworzyć plik modułu, aby zbadać jego kod:

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

Jak widzisz, to dość prosty mały skrypt, który rozpakowuje plik i liczy, ile linii zawiera.

??? info "Co robi `debug true`?"

    Dyrektywa `debug true` w definicji procesu powoduje, że Nextflow wypisuje wyjście ze skryptu (jak liczba linii "40") bezpośrednio w logu wykonania.
    Bez tego widziałbyś tylko status wykonania procesu, ale nie rzeczywiste wyjście ze skryptu.

    Więcej informacji o debugowaniu procesów Nextflow znajdziesz w zadaniu pobocznym [Debugowanie workflow'ów Nextflow](debugging.md).

#### 1.4.2. Dodaj wywołanie `COUNT_LINES`

Teraz, gdy proces jest dostępny dla workflow'a, możemy dodać wywołanie procesu `COUNT_LINES`, aby uruchomić go na pliku wejściowym.

Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Policz linie w pliku
        COUNT_LINES(myFile)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

A teraz uruchom workflow'a:

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

To pokazuje, że jesteśmy w stanie odpowiednio operować na pliku wewnątrz procesu.

Konkretnie, Nextflow pomyślnie wykonał następujące operacje:

- Umieścił plik w katalogu roboczym
- Rozpakował plik .gz
- Policzył linie (40 linii w tym przypadku)
- Zakończył bez błędu

Kluczem do tego płynnego działania jest to, że wyraźnie mówimy Nextflow'owi, że nasze wejście jest plikiem i powinno być traktowane jako takie.

### 1.5. Rozwiązuj podstawowe błędy wejścia pliku

To często sprawia problemy nowicjuszom w Nextflow, więc poświęćmy kilka minut na przyjrzenie się, co się dzieje, gdy robisz to źle.

Są dwa główne miejsca, w których możesz źle obsłużyć plik: na poziomie workflow'a i na poziomie procesu.

#### 1.5.1. Błąd na poziomie workflow'a

Zobaczmy, co się stanie, jeśli wrócimy do traktowania pliku jako ciągu znaków, gdy określamy wejście w bloku workflow.

Wprowadź następujące zmiany w workflow'ie, upewniając się, że zakomentowałeś instrukcje wypisywania specyficzne dla ścieżki:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Policz linie w pliku
        COUNT_LINES(myFile)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Policz linie w pliku
        COUNT_LINES(myFile)
    ```

A teraz uruchom workflow'a:

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

To jest ważny fragment:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Gdy określasz wejście `path`, Nextflow sprawdza, czy przekazujesz rzeczywiste odniesienia do plików, a nie tylko ciągi znaków.
Ten błąd mówi Ci, że `'data/patientA_rep1_normal_R1_001.fastq.gz'` nie jest prawidłową wartością ścieżki, ponieważ jest to ciąg znaków, a nie obiekt Path.

Nextflow natychmiast wykrył problem i zatrzymał się przed rozpoczęciem procesu.

#### 1.5.2. Błąd na poziomie procesu

Drugim miejscem, w którym możemy zapomnieć określić, że chcemy, aby Nextflow traktował wejście jako plik, jest definicja procesu.

!!! warning "Zachowaj błąd workflow'a z 1.5.1"

    Aby ten test działał poprawnie, zachowaj workflow'a w jego uszkodzonym stanie (używając zwykłego ciągu znaków zamiast `file()`).
    W połączeniu z `val` w procesie daje to błąd pokazany poniżej.

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

A teraz uruchom workflow'a ponownie:

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

To pokazuje wiele szczegółów o błędzie, ponieważ proces jest ustawiony na wypisywanie informacji debugowania, jak wspomniano powyżej.

To są najbardziej istotne sekcje:

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

To mówi, że system nie mógł znaleźć pliku; jednak jeśli sprawdzisz ścieżkę, istnieje plik o tej nazwie w tej lokalizacji.

Gdy to uruchomiliśmy, Nextflow przekazał wartość ciągu znaków do skryptu, ale nie _umieścił_ rzeczywistego pliku w katalogu roboczym.
Więc proces próbował użyć względnego ciągu znaków, `data/patientA_rep1_normal_R1_001.fastq.gz`, ale ten plik nie istnieje w katalogu roboczym procesu.

Razem te dwa przykłady pokazują, jak ważne jest powiedzenie Nextflow'owi, czy wejście powinno być obsługiwane jako plik.

!!! note "Uwaga"

    Upewnij się, że cofniesz oba celowe błędy przed przejściem do następnej sekcji.

### Podsumowanie

- Ciągi znaków ścieżek vs obiekty Path: Ciągi znaków to tylko tekst, obiekty Path to inteligentne odniesienia do plików
- Metoda `file()` konwertuje ciąg znaków ścieżki na obiekt Path, z którym Nextflow może pracować
- Możesz uzyskać dostęp do właściwości pliku, takich jak `name`, `simpleName`, `extension` i `parent` [używając atrybutów pliku](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Używanie obiektów Path zamiast ciągów znaków pozwala Nextflow'owi prawidłowo zarządzać plikami w Twoim workflow'ie
- Wyniki wejścia procesu: Prawidłowa obsługa plików wymaga obiektów Path, a nie ciągów znaków, aby zapewnić, że pliki są poprawnie umieszczane i dostępne do użycia przez procesy

---

## 2. Używanie plików zdalnych

Jedną z kluczowych funkcji Nextflow'a jest możliwość płynnego przełączania się między plikami lokalnymi (na tej samej maszynie) a plikami zdalnymi dostępnymi przez internet.

Jeśli robisz to dobrze, nigdy nie powinieneś musieć zmieniać logiki swojego workflow'a, aby obsłużyć pliki pochodzące z różnych lokalizacji.
Wszystko, co musisz zrobić, aby użyć pliku zdalnego, to określić odpowiedni prefiks w ścieżce pliku, gdy dostarczasz go do workflow'a.

Na przykład `/path/to/data` nie ma prefiksu, co wskazuje, że jest to 'normalna' lokalna ścieżka pliku, podczas gdy `s3://path/to/data` zawiera prefiks `s3://`, wskazując, że znajduje się w magazynie obiektów S3 Amazon.

Obsługiwanych jest wiele różnych protokołów:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Aby użyć któregokolwiek z nich, po prostu określ odpowiedni prefiks w ciągu znaków, który jest wtedy technicznie nazywany Uniform Resource Identifier (URI) zamiast ścieżki pliku.
Nextflow obsłuży uwierzytelnianie i umieszczanie plików we właściwym miejscu, pobieranie lub przesyłanie oraz wszystkie inne operacje na plikach, których można się spodziewać.

Kluczową zaletą tego systemu jest to, że umożliwia nam przełączanie się między środowiskami bez zmiany jakiejkolwiek logiki pipeline'a.
Na przykład możesz rozwijać z małym, lokalnym zestawem testowym przed przełączeniem się na pełnoskalowy zestaw testowy znajdujący się w zdalnym magazynie, po prostu zmieniając URI.

### 2.1. Użyj pliku z internetu

Przetestujmy to, zamieniając lokalną ścieżkę, którą dostarczamy do naszego workflow'a, na ścieżkę HTTPS wskazującą na kopię tych samych danych, która jest przechowywana w Github.

!!! warning "Ostrzeżenie"

    To zadziała tylko, jeśli masz aktywne połączenie internetowe.

Otwórz ponownie `main.nf` i zmień ścieżkę wejściową w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Używanie pliku zdalnego z internetu
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
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Uruchommy workflow'a:

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

Działa! Widzisz, że bardzo niewiele się zmieniło.

Jedyną różnicą w wyjściu konsoli jest to, że klasa obiektu ścieżki to teraz `nextflow.file.http.XPath`, podczas gdy dla ścieżki lokalnej klasa była `sun.nio.fs.UnixPath`.
Nie musisz pamiętać tych klas; wspominamy o tym tylko po to, aby pokazać, że Nextflow identyfikuje i obsługuje różne lokalizacje odpowiednio.

Za kulisami Nextflow pobrał plik do katalogu tymczasowego znajdującego się w katalogu roboczym.
Ten umieszczony plik może być następnie traktowany jako plik lokalny i dowiązany symbolicznie do odpowiedniego katalogu procesu.

Możesz to zweryfikować, przeglądając zawartość katalogu roboczego znajdującego się pod wartością hash procesu.

??? abstract "Zawartość katalogu roboczego"

    Jeśli hash procesu był `8a/2ab7ca`, możesz zbadać katalog roboczy:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Dowiązanie symboliczne wskazuje na umieszczoną kopię pliku zdalnego, którą Nextflow automatycznie pobrał.

Zauważ, że w przypadku większych plików krok pobierania zajmie dodatkowy czas w porównaniu z uruchamianiem na plikach lokalnych.
Jednak Nextflow sprawdza, czy ma już umieszczoną kopię, aby uniknąć niepotrzebnych pobrań.
Więc jeśli uruchomisz ponownie na tym samym pliku i nie usunąłeś umieszczonego pliku, Nextflow użyje umieszczonej kopii.

To pokazuje, jak łatwo jest przełączać się między danymi lokalnymi i zdalnymi przy użyciu Nextflow'a, co jest kluczową funkcją Nextflow'a.

!!! note "Uwaga"

    Jedynym ważnym wyjątkiem od tej zasady jest to, że nie możesz używać wzorców glob ani ścieżek katalogów z HTTPS, ponieważ HTTPS nie może wyświetlić wielu plików, więc musisz określić dokładne adresy URL plików.
    Jednak inne protokoły magazynowania, takie jak magazyn obiektów (`s3://`, `az://`, `gs://`), mogą używać zarówno globów, jak i ścieżek katalogów.

    Oto jak możesz używać wzorców glob z magazynem w chmurze:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // S3 ze wzorcami glob - dopasowałoby wiele plików
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage ze wzorcami glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage ze wzorcami glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Pokażemy Ci, jak pracować z globami w praktyce w następnej sekcji.

### 2.2. Przełącz się z powrotem na plik lokalny

Wrócimy do używania naszych lokalnych przykładowych plików przez resztę tego zadania pobocznego, więc przełączmy wejście workflow'a z powrotem na oryginalny plik:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utwórz obiekt Path z ciągu znaków ścieżki
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
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Podsumowanie

- Dostęp do danych zdalnych odbywa się przy użyciu URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow automatycznie pobierze i umieści dane we właściwym miejscu, o ile te ścieżki są przekazywane do procesów
- Nie pisz logiki do pobierania lub przesyłania plików zdalnych!
- Pliki lokalne i zdalne tworzą różne typy obiektów, ale działają identycznie
- **Ważne**: HTTP/HTTPS działa tylko z pojedynczymi plikami (bez wzorców glob)
- Magazyn w chmurze (S3, Azure, GCS) obsługuje zarówno pojedyncze pliki, jak i wzorce glob
- Możesz płynnie przełączać się między lokalnymi i zdalnymi źródłami danych bez zmiany logiki kodu (o ile protokół obsługuje wymagane operacje)

---

## 3. Używanie fabryki kanałów `fromPath()`

Do tej pory pracowaliśmy z pojedynczym plikiem na raz, ale w Nextflow zazwyczaj będziemy chcieli utworzyć kanał wejściowy z wieloma plikami wejściowymi do przetworzenia.

Naiwnym sposobem na to byłoby połączenie metody `file()` z [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) w ten sposób:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

To działa, ale jest niezgrabne.

!!! tip "Kiedy używać `file()` vs `channel.fromPath()`"

    - Używaj `file()`, gdy potrzebujesz pojedynczego obiektu Path do bezpośredniej manipulacji (sprawdzanie, czy plik istnieje, odczytywanie jego atrybutów lub przekazywanie do pojedynczego wywołania procesu)
    - Używaj `channel.fromPath()`, gdy potrzebujesz kanału, który może zawierać wiele plików, szczególnie ze wzorcami glob, lub gdy pliki będą przepływać przez wiele procesów

Tu wkracza [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): wygodna fabryka kanałów, która łączy całą funkcjonalność potrzebną do wygenerowania kanału z jednego lub więcej statycznych ciągów znaków plików, a także wzorców glob.

### 3.1. Dodaj fabrykę kanałów

Zaktualizujmy nasz workflow'a, aby używał `channel.fromPath`.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Wypisz atrybuty pliku
        /* Zakomentuj je na razie, wrócimy do nich!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Policz linie w pliku
        // COUNT_LINES(myFile)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Utwórz obiekt Path z ciągu znaków ścieżki
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Policz linie w pliku
        COUNT_LINES(myFile)
    ```

Zakomentowaliśmy również kod, który wypisuje atrybuty na razie, i dodaliśmy instrukcję `.view`, aby wypisać tylko nazwę pliku.

Uruchom workflow'a:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Jak widzisz, ścieżka pliku jest ładowana jako obiekt typu `Path` w kanale.
Jest to podobne do tego, co zrobiłby `file()`, z tym że teraz mamy kanał, do którego możemy załadować więcej plików, jeśli chcemy.

Używanie `channel.fromPath()` to wygodny sposób tworzenia nowego kanału wypełnionego listą plików.

### 3.2. Wyświetl atrybuty plików w kanale

W naszym pierwszym podejściu do używania fabryki kanałów uprościliśmy kod i po prostu wypisaliśmy nazwę pliku.

Wróćmy do wypisywania pełnych atrybutów pliku:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Policz linie w pliku
        COUNT_LINES(ch_files)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Policz linie w pliku
        // COUNT_LINES(ch_files)
    ```

Ponownie włączamy również wywołanie procesu `COUNT_LINES`, aby zweryfikować, że przetwarzanie plików nadal działa poprawnie z naszym podejściem opartym na kanałach.

Uruchom workflow'a:

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

I oto jesteś, te same wyniki co wcześniej, ale teraz mamy plik w kanale, więc możemy dodać więcej.

### 3.3. Używanie glob do dopasowania wielu plików

Istnieje kilka sposobów, w jakie moglibyśmy załadować więcej plików do kanału.
Tutaj pokażemy Ci, jak używać wzorców glob, które są wygodnym sposobem dopasowywania i pobierania nazw plików i katalogów na podstawie znaków wieloznacznych.
Proces dopasowywania tych wzorców nazywa się "globowaniem" lub "rozszerzaniem nazw plików".

!!! note "Uwaga"

    Jak wspomniano wcześniej, Nextflow obsługuje globowanie do zarządzania plikami wejściowymi i wyjściowymi w większości przypadków, z wyjątkiem ścieżek plików HTTPS, ponieważ HTTPS nie może wyświetlić wielu plików.

Powiedzmy, że chcemy pobrać oba pliki w parze plików powiązanych z danym pacjentem, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Ponieważ jedyną różnicą między nazwami plików jest numer repliki, _tj._ liczba po `R`, możemy użyć znaku wieloznacznego `*`, aby zastąpić liczbę w następujący sposób:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

To jest wzorzec glob, którego potrzebujemy.

Teraz wszystko, co musimy zrobić, to zaktualizować ścieżkę pliku w fabryce kanałów, aby używała tego wzorca glob w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="7"
      // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7"
      // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow automatycznie rozpozna, że jest to wzorzec glob i obsłuży go odpowiednio.

Uruchom workflow'a, aby to przetestować:

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

Jak widzisz, mamy teraz dwa obiekty Path w naszym kanale, co pokazuje, że Nextflow poprawnie wykonał rozszerzenie nazw plików i załadował oraz przetworzył oba pliki zgodnie z oczekiwaniami.

Używając tej metody, możemy pobrać tyle lub tak mało plików, ile chcemy, po prostu zmieniając wzorzec glob. Gdybyśmy uczynili go bardziej hojnym, na przykład zastępując wszystkie zmienne części nazw plików przez `*` (_np._ `data/patient*_rep*_*_R*_001.fastq.gz`), moglibyśmy pobrać wszystkie przykładowe pliki w katalogu `data`.

### Podsumowanie

- `channel.fromPath()` tworzy kanał z plikami pasującymi do wzorca
- Każdy plik jest emitowany jako oddzielny element w kanale
- Możemy użyć wzorca glob do dopasowania wielu plików
- Pliki są automatycznie konwertowane na obiekty Path z pełnymi atrybutami
- Metoda `.view()` pozwala na inspekcję zawartości kanału

---

## 4. Wydobywanie podstawowych metadanych z nazw plików

W większości dziedzin naukowych bardzo powszechne jest kodowanie metadanych w nazwach plików zawierających dane.
Na przykład w bioinformatyce pliki zawierające dane sekwencjonowania są często nazywane w sposób, który koduje informacje o próbce, warunku, replice i numerze odczytu.

Jeśli nazwy plików są konstruowane zgodnie ze spójną konwencją, możesz wydobyć te metadane w ustandaryzowany sposób i użyć ich w trakcie analizy.
To duże "jeśli", oczywiście, i powinieneś być bardzo ostrożny, gdy polegasz na strukturze nazw plików; ale rzeczywistość jest taka, że to podejście jest bardzo szeroko stosowane, więc przyjrzyjmy się, jak to się robi w Nextflow.

W przypadku naszych przykładowych danych wiemy, że nazwy plików zawierają spójnie ustrukturyzowane metadane.
Na przykład nazwa pliku `patientA_rep1_normal_R2_001` koduje następujące informacje:

- ID pacjenta: `patientA`
- ID repliki: `rep1`
- typ próbki: `normal` (w przeciwieństwie do `tumor`)
- zestaw odczytów: `R1` (w przeciwieństwie do `R2`)

Zmodyfikujemy nasz workflow'a, aby pobrać te informacje w trzech krokach:

1. Pobrać `simpleName` pliku, który zawiera metadane
2. Rozdzielić metadane przy użyciu metody zwanej `tokenize()`
3. Użyć mapy do zorganizowania metadanych

!!! warning "Ostrzeżenie"

    Nigdy nie powinieneś kodować wrażliwych informacji w nazwach plików, takich jak nazwiska pacjentów lub inne cechy identyfikujące, ponieważ może to naruszyć prywatność pacjentów lub inne odpowiednie ograniczenia bezpieczeństwa.

### 4.1. Pobierz `simpleName`

`simpleName` to atrybut pliku, który odpowiada nazwie pliku pozbawionej ścieżki i rozszerzenia.

Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

To pobiera `simpleName` i kojarzy go z pełnym obiektem pliku przy użyciu operacji `map()`.

Uruchom workflow'a, aby przetestować, czy działa:

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

### 4.2. Wydobądź metadane z `simplename`

W tym momencie metadane, których chcemy, są osadzone w `simplename`, ale nie możemy bezpośrednio uzyskać dostępu do poszczególnych elementów.
Więc musimy podzielić `simplename` na jego składniki.
Na szczęście te składniki są po prostu oddzielone podkreśleniami w oryginalnej nazwie pliku, więc możemy zastosować powszechną metodę Nextflow zwaną `tokenize()`, która jest idealna do tego zadania.

Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Metoda `tokenize()` podzieli ciąg znaków `simpleName` wszędzie tam, gdzie znajdzie podkreślenia, i zwróci listę zawierającą podciągi.

Uruchom workflow'a:

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

Teraz krotka dla każdego elementu w naszym kanale zawiera listę metadanych (_np._ `[patientA, rep1, normal, R1, 001]`) i oryginalny obiekt pliku.

To świetnie!
Rozbiliśmy informacje o naszym pacjencie z pojedynczego ciągu znaków na listę ciągów znaków.
Możemy teraz obsługiwać każdą część informacji o pacjencie osobno.

### 4.3. Użyj mapy do zorganizowania metadanych

Nasze metadane są obecnie tylko płaską listą.
Łatwo się z nimi pracuje, ale trudno je czytać.

```console
[patientA, rep1, normal, R1, 001]
```

Co to jest element pod indeksem 3? Czy możesz powiedzieć bez odniesienia się do oryginalnego wyjaśnienia struktury metadanych?

To świetna okazja do użycia magazynu klucz-wartość, gdzie każdy element ma zestaw kluczy i powiązanych z nimi wartości, więc możesz łatwo odwołać się do każdego klucza, aby uzyskać odpowiednią wartość.

W naszym przykładzie oznacza to przejście z tej organizacji:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Do tej:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

W Nextflow nazywa się to [mapą](https://nextflow.io/docs/latest/script.html#maps).

Przekonwertujmy teraz naszą płaską listę na mapę.
Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Załaduj pliki za pomocą channel.fromPath
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
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Kluczowe zmiany to:

- **Przypisanie destrukturyzujące**: `def (patient, replicate, type, readNum) = ...` wydobywa tokenizowane wartości do nazwanych zmiennych w jednej linii
- **Składnia literału mapy**: `[id: patient, replicate: ...]` tworzy mapę, w której każdy klucz (jak `id`) jest powiązany z wartością (jak `patient`)
- **Struktura zagnieżdżona**: Zewnętrzna lista `[..., myFile]` paruje mapę metadanych z oryginalnym obiektem pliku

Uprościliśmy również kilka ciągów znaków metadanych przy użyciu metody zastępowania ciągów znaków zwanej `replace()`, aby usunąć niektóre znaki, które są niepotrzebne (_np._ `replicate.replace('rep', '')`, aby zachować tylko liczbę z ID replik).

Uruchommy workflow'a ponownie:

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

Teraz metadane są starannie oznaczone (_np._ `[id:patientA, replicate:1, type:normal, readNum:2]`), więc o wiele łatwiej jest powiedzieć, co jest czym.

Będzie również o wiele łatwiej faktycznie wykorzystać elementy metadanych w workflow'ie i sprawi, że nasz kod będzie łatwiejszy do odczytania i bardziej łatwy w utrzymaniu.

### Podsumowanie

- Możemy obsługiwać nazwy plików w Nextflow z mocą pełnego języka programowania
- Możemy traktować nazwy plików jako ciągi znaków, aby wydobyć odpowiednie informacje
- Użycie metod takich jak `tokenize()` i `replace()` pozwala nam manipulować ciągami znaków w nazwie pliku
- Operacja `.map()` przekształca elementy kanału, zachowując strukturę
- Ustrukturyzowane metadane (mapy) sprawiają, że kod jest bardziej czytelny i łatwiejszy w utrzymaniu niż listy pozycyjne

Następnie przyjrzymy się, jak obsługiwać sparowane pliki danych.

---

## 5. Obsługa sparowanych plików danych

Wiele projektów eksperymentalnych tworzy sparowane pliki danych, które korzystają z obsługi w sposób wyraźnie sparowany.
Na przykład w bioinformatyce dane sekwencjonowania są często generowane w postaci sparowanych odczytów, co oznacza ciągi sekwencji pochodzące z tego samego fragmentu DNA (często nazywane 'do przodu' i 'do tyłu', ponieważ są odczytywane z przeciwnych końców).

Tak jest w przypadku naszych przykładowych danych, gdzie R1 i R2 odnoszą się do dwóch zestawów odczytów.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow zapewnia wyspecjalizowaną fabrykę kanałów do pracy ze sparowanymi plikami o nazwie `channel.fromFilePairs()`, która automatycznie grupuje pliki na podstawie wspólnego wzorca nazewnictwa. Pozwala to na ściślejsze powiązanie sparowanych plików przy mniejszym wysiłku.

Zmodyfikujemy nasz workflow'a, aby to wykorzystać.
Zajmie to dwa kroki:

1. Przełączyć fabrykę kanałów na `channel.fromFilePairs()`
2. Wydobyć i zmapować metadane

### 5.1. Przełącz fabrykę kanałów na `channel.fromFilePairs()`

Aby użyć `channel.fromFilePairs`, musimy określić wzorzec, którego Nextflow powinien użyć do identyfikacji dwóch członków w parze.

Wracając do naszych przykładowych danych, możemy sformalizować wzorzec nazewnictwa w następujący sposób:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Jest to podobne do wzorca glob, którego użyliśmy wcześniej, z tym że to konkretnie wylicza podciągi (albo `1`, albo `2` pojawiające się zaraz po R), które identyfikują dwóch członków pary.

Zaktualizujmy odpowiednio workflow'a `main.nf`:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Załaduj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Zakomentuj mapowanie na razie, wrócimy do niego!
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
        // Załaduj pliki za pomocą channel.fromFilePairs
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
Dodamy to z powrotem później, z kilkoma modyfikacjami.

Uruchom workflow'a, aby go przetestować:

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

Ojej, tym razem uruchomienie nie powiodło się!

Odpowiedni fragment komunikatu o błędzie jest tutaj:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Dzieje się tak, ponieważ zmieniliśmy fabrykę kanałów.
Do tej pory oryginalny kanał wejściowy zawierał tylko ścieżki plików.
Cała manipulacja metadanymi, którą wykonywaliśmy, nie wpłynęła faktycznie na zawartość kanału.

Teraz, gdy używamy fabryki kanałów `.fromFilePairs`, zawartość wynikowego kanału jest inna.
Widzimy tylko jeden element kanału, złożony z krotki zawierającej dwa elementy: część `simpleName` wspólną dla dwóch plików, która służy jako identyfikator, oraz krotkę zawierającą dwa obiekty plików, w formacie `id, [ file1, file2 ]`.

To świetnie, ponieważ Nextflow wykonał ciężką pracę wydobycia nazwy pacjenta poprzez zbadanie wspólnego prefiksu i użycie go jako identyfikatora pacjenta.

Jednak to psuje nasz obecny workflow'a.
Gdybyśmy chcieli nadal uruchamiać `COUNT_LINES` w ten sam sposób bez zmiany procesu, musielibyśmy zastosować operację mapowania, aby wydobyć ścieżki plików.
Ale nie zamierzamy tego robić, ponieważ naszym ostatecznym celem jest użycie innego procesu, `ANALYZE_READS`, który odpowiednio obsługuje pary plików.

Więc po prostu zakomentujmy (lub usuńmy) wywołanie `COUNT_LINES` i przejdźmy dalej.

=== "Po"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Policz linie w pliku
        // COUNT_LINES(ch_files)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Policz linie w pliku
        COUNT_LINES(ch_files)
    ```

Możesz również zakomentować lub usunąć instrukcję include `COUNT_LINES`, ale nie będzie to miało żadnego efektu funkcjonalnego.

Teraz uruchommy workflow'a ponownie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Hurra, tym razem workflow'a się udaje!

Jednak nadal musimy wydobyć resztę metadanych z pola `id`.

### 5.2. Wydobądź i zorganizuj metadane z par plików

Nasza operacja `map` z wcześniej nie zadziała, ponieważ nie pasuje do struktury danych, ale możemy ją zmodyfikować, aby działała.

Mamy już dostęp do rzeczywistego identyfikatora pacjenta w ciągu znaków, którego `fromFilePairs()` użył jako identyfikatora, więc możemy go użyć do wydobycia metadanych bez pobierania `simpleName` z obiektu Path, jak robiliśmy wcześniej.

Odkomentuj operację map w workflow'ie i wprowadź następujące zmiany:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // Załaduj pliki za pomocą channel.fromFilePairs
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
        // Załaduj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Zakomentuj mapowanie na razie, wrócimy do niego!
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

Zauważ również, że usunęliśmy `readNum` z linii `tokenize()`; wszelkie podciągi, których nie nazwiemy konkretnie (zaczynając od lewej), zostaną po cichu porzucone.
Możemy to zrobić, ponieważ sparowane pliki są teraz ściśle powiązane, więc nie potrzebujemy już `readNum` w mapie metadanych.

Uruchommy workflow'a:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

I oto jest: mamy mapę metadanych (`[id:patientA, replicate:1, type:normal]`) w pierwszej pozycji krotki wyjściowej, po której następuje krotka sparowanych plików, zgodnie z zamierzeniem.

Oczywiście to pobierze i przetworzy tylko tę konkretną parę plików.
Jeśli chcesz poeksperymentować z przetwarzaniem wielu par, możesz spróbować dodać symbole wieloznaczne do wzorca wejściowego i zobaczyć, co się stanie.
Na przykład spróbuj użyć `data/patientA_rep1_*_R{1,2}_001.fastq.gz`

### Podsumowanie

- [`channel.fromFilePairs()` automatycznie znajduje i paruje powiązane pliki](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- To upraszcza obsługę odczytów sparowanych w Twoim pipeline'ie
- Sparowane pliki mogą być grupowane jako krotki `[id, [file1, file2]]`
- Wydobywanie metadanych można wykonać z ID sparowanego pliku, a nie z poszczególnych plików

---

## 6. Używanie operacji na plikach w procesach

Teraz połączmy to wszystko w prostym procesie, aby wzmocnić, jak używać operacji na plikach wewnątrz procesu Nextflow.

Dostarczamy Ci wstępnie napisany moduł procesu o nazwie `ANALYZE_READS`, który przyjmuje krotkę metadanych i parę plików wejściowych i je analizuje.
Moglibyśmy sobie wyobrazić, że to robi dopasowanie sekwencji lub wykrywanie wariantów lub jakikolwiek inny krok, który ma sens dla tego typu danych.

Zaczynajmy.

### 6.1. Zaimportuj proces i zbadaj kod

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

Możesz otworzyć plik modułu, aby zbadać jego kod:

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

    Dyrektywy `tag` i `publishDir` używają składni domknięcia (`{ ... }`) zamiast interpolacji ciągów znaków (`"${...}"`).
    Dzieje się tak, ponieważ te dyrektywy odwołują się do zmiennych wejściowych (`meta`), które nie są dostępne do czasu wykonania.
    Składnia domknięcia odkłada ewaluację do momentu, gdy proces faktycznie się uruchomi.

!!! note "Uwaga"

    Nazywamy naszą mapę metadanych `meta` zgodnie z konwencją.
    Aby głębiej zanurzyć się w mapy meta, zobacz zadanie poboczne [Metadane i mapy meta](./metadata.md).

### 6.2. Wywołaj proces w workflow'ie

Teraz, gdy proces jest dostępny dla workflow'a, możemy dodać wywołanie procesu `ANALYZE_READS`, aby go uruchomić.

Aby uruchomić go na naszych przykładowych danych, będziemy musieli zrobić dwie rzeczy:

1. Nadać nazwę przemapowanemu kanałowi
2. Dodać wywołanie procesu

#### 6.2.1. Nazwij przemapowany kanał wejściowy

Wcześniej zastosowaliśmy manipulacje mapowania bezpośrednio do kanału wejściowego.
Aby przekazać przemapowaną zawartość do procesu `ANALYZE_READS` (i zrobić to w sposób jasny i łatwy do odczytania), chcemy utworzyć nowy kanał o nazwie `ch_samples`.

Możemy to zrobić przy użyciu operatora [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

W głównym workflow'ie zastąp operator `.view()` przez `.set { ch_samples }` i dodaj linię testującą, że możemy odwołać się do kanału po nazwie.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // Załaduj pliki za pomocą channel.fromFilePairs
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

        // Tymczasowo: zajrzyj do ch_samples
        ch_samples.view()
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // Załaduj pliki za pomocą channel.fromFilePairs
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

To potwierdza, że możemy teraz odwołać się do kanału po nazwie.

#### 6.2.2. Wywołaj proces na danych

Teraz faktycznie wywołajmy proces `ANALYZE_READS` na kanale `ch_samples`.

W głównym workflow'ie wprowadź następujące zmiany w kodzie:

=== "Po"

    ```groovy title="main.nf" linenums="23"
        // Uruchom analizę
        ANALYZE_READS(ch_samples)
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="23"
        // Tymczasowo: zajrzyj do ch_samples
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

Ten proces jest skonfigurowany do publikowania swoich wyjść do katalogu `results`, więc zajrzyj tam.

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

Proces wziął nasze wejścia i utworzył nowy plik zawierający metadane pacjenta, zgodnie z projektem.
Wspaniale!

### 6.3. Uwzględnij znacznie więcej pacjentów

Oczywiście to przetwarza tylko pojedynczą parę plików dla pojedynczego pacjenta, co nie jest dokładnie tym rodzajem wysokiej przepustowości, na którą liczysz z Nextflow.
Prawdopodobnie będziesz chciał przetwarzać znacznie więcej danych na raz.

Pamiętaj, że `channel.fromPath()` akceptuje _glob_ jako wejście, co oznacza, że może akceptować dowolną liczbę plików pasujących do wzorca.
Dlatego jeśli chcemy uwzględnić wszystkich pacjentów, możemy po prostu zmodyfikować ciąg wejściowy, aby uwzględnić więcej pacjentów, jak wspomniano mimochodem wcześniej.

Udajmy, że chcemy być tak chciwi, jak to możliwe.
Wprowadź następujące zmiany w workflow'ie:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Załaduj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // Załaduj pliki za pomocą channel.fromFilePairs
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Uruchom pipeline'a ponownie:

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

Może nie.
Jeśli przyjrzysz się bliżej, mamy problem: mamy dwie repliki dla patientA, ale tylko jeden plik wyjściowy!
Nadpisujemy plik wyjściowy za każdym razem.

### 6.4. Uczyń opublikowane pliki unikalnymi

Ponieważ mamy dostęp do metadanych pacjenta, możemy ich użyć, aby uczynić opublikowane pliki unikalnymi, włączając różnicujące metadane, albo w strukturze katalogów, albo w samych nazwach plików.

Wprowadź następującą zmianę w workflow'ie:

=== "Po"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Przed"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Tutaj pokazujemy opcję użycia dodatkowych poziomów katalogów, aby uwzględnić typy próbek i repliki, ale możesz poeksperymentować z robieniem tego na poziomie nazwy pliku również.

Teraz uruchom pipeline'a jeszcze raz, ale upewnij się, że najpierw usuniesz katalog wyników, aby dać sobie czystą przestrzeń roboczą:

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

I oto jest, wszystkie nasze metadane, starannie zorganizowane. To sukces!

Jest o wiele więcej, co możesz zrobić, gdy masz swoje metadane załadowane do mapy w ten sposób:

1. Tworzyć zorganizowane katalogi wyjściowe na podstawie atrybutów pacjenta
2. Podejmować decyzje w procesach na podstawie właściwości pacjenta
3. Dzielić, łączyć i ponownie łączyć dane na podstawie wartości metadanych

Ten wzorzec utrzymywania metadanych w sposób jawny i dołączony do danych (zamiast kodowania w nazwach plików) jest podstawową najlepszą praktyką w Nextflow, która umożliwia budowanie solidnych, łatwych w utrzymaniu workflow'ów analizy.
Możesz dowiedzieć się więcej o tym w zadaniu pobocznym [Metadane i mapy meta](./metadata.md).

### Podsumowanie

- Dyrektywa `publishDir` może organizować wyjścia na podstawie wartości metadanych
- Metadane w krotkach umożliwiają ustrukturyzowaną organizację wyników
- To podejście tworzy łatwe w utrzymaniu workflow'y z jasnym pochodzeniem danych
- Procesy mogą przyjmować krotki metadanych i plików jako wejście
- Dyrektywa `tag` zapewnia identyfikację procesu w logach wykonania
- Struktura workflow'a oddziela tworzenie kanałów od wykonywania procesów

---

## Podsumowanie

W tym zadaniu pobocznym nauczyłeś się, jak pracować z plikami w Nextflow, od podstawowych operacji po bardziej zaawansowane techniki obsługi kolekcji plików.

Zastosowanie tych technik we własnej pracy umożliwi Ci budowanie bardziej wydajnych i łatwych w utrzymaniu workflow'ów, szczególnie podczas pracy z dużą liczbą plików o złożonych konwencjach nazewnictwa.

### Kluczowe wzorce

1.  **Podstawowe operacje na plikach:** Utworzyliśmy obiekty Path za pomocą `file()` i uzyskaliśmy dostęp do atrybutów plików, takich jak nazwa, rozszerzenie i katalog nadrzędny, ucząc się różnicy między ciągami znaków a obiektami Path.

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

2.  **Używanie plików zdalnych**: Nauczyliśmy się, jak przezroczyście przełączać się między plikami lokalnymi i zdalnymi przy użyciu URI, demonstrując zdolność Nextflow'a do obsługi plików z różnych źródeł bez zmiany logiki workflow'a.

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

3.  **Ładowanie plików przy użyciu fabryki kanałów `fromPath()`:** Utworzyliśmy kanały ze wzorców plików za pomocą `channel.fromPath()` i wyświetliliśmy ich atrybuty plików, w tym typy obiektów.

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

4.  **Wydobywanie metadanych pacjenta z nazw plików:** Użyliśmy `tokenize()` i `replace()` do wydobycia i strukturyzacji metadanych z nazw plików, konwertując je na zorganizowane mapy.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **Uproszczenie za pomocą channel.fromFilePairs:** Użyliśmy `channel.fromFilePairs()` do automatycznego parowania powiązanych plików i wydobywania metadanych z ID sparowanych plików.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Używanie operacji na plikach w procesach:** Zintegrowaliśmy operacje na plikach z procesami Nextflow przy odpowiedniej obsłudze wejść, używając `publishDir` do organizowania wyjść na podstawie metadanych.

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

- [Dokumentacja Nextflow: Praca z plikami](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Co dalej?

Wróć do [menu zadań pobocznych](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
