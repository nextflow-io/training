# Przetwarzanie plików wejściowych

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Przepływy pracy analizy naukowej często obejmują przetwarzanie dużej liczby plików.
Nextflow zapewnia potężne narzędzia do efektywnej obsługi plików, pomagając organizować i przetwarzać dane przy minimalnym nakładzie kodu.

### Cele nauki

W tym side queście zbadamy, jak Nextflow obsługuje pliki, od podstawowych operacji na plikach po bardziej zaawansowane techniki pracy z kolekcjami plików.
Nauczysz się wyodrębniać metadane z nazw plików, co jest powszechnym wymaganiem w pipeline'ach analizy naukowej.

Pod koniec tego side questa będziesz potrafić:

- Tworzyć obiekty Path z ciągów znaków ścieżek plików za pomocą metody `file()` w Nextflow
- Uzyskiwać dostęp do atrybutów pliku, takich jak nazwa, rozszerzenie i katalog nadrzędny
- Obsługiwać zarówno pliki lokalne, jak i zdalne w sposób przejrzysty przy użyciu URI
- Używać kanałów do automatyzacji obsługi plików za pomocą `channel.fromPath()` i `channel.fromFilePairs()`
- Wyodrębniać i strukturyzować metadane z nazw plików za pomocą manipulacji ciągami znaków
- Grupować powiązane pliki przy użyciu dopasowywania wzorców i wyrażeń glob
- Integrować operacje na plikach z procesami Nextflow z prawidłową obsługą wejścia
- Organizować wyjścia procesów przy użyciu struktur katalogów opartych na metadanych

Te umiejętności pomogą Ci budować workflow'y, które mogą obsługiwać różne rodzaje plików wejściowych z dużą elastycznością.

### Wymagania wstępne

Przed podjęciem tego side questa powinieneś:

- Ukończyć tutorial [Hello Nextflow](../../hello_nextflow/) lub równoważny kurs dla początkujących.
- Czuć się komfortowo z podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory)

---

## 0. Rozpocznij

#### Otwórz środowisko szkoleniowe codespace

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzyłeś środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki do tego tutoriala.

```bash
cd side-quests/working_with_files
```

Możesz ustawić VSCode, aby skupić się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz prosty plik workflow o nazwie `main.nf`, katalog `modules` zawierający dwa pliki modułów oraz katalog `data` zawierający przykładowe pliki danych.

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
Jeśli nie znasz analizy nowotworów, wiedz tylko, że odpowiada to modelowi eksperymentalnemu wykorzystującemu sparowane próbki guz/normalne do wykonywania analiz kontrastowych.

Konkretnie dla pacjenta A mamy dwa zestawy replikatów technicznych (powtórzeń).

Pliki danych sekwencjonowania są nazwane zgodnie z typową konwencją `_R1_` i `_R2_` dla tzw. 'odczytów w przód' i 'odczytów w tył'.

_Nie martw się, jeśli nie znasz tego projektu eksperymentalnego, nie jest to krytyczne dla zrozumienia tego tutoriala._

#### Przejrzyj zadanie

Twoim wyzwaniem jest napisanie workflow Nextflow, który będzie:

1. **Wczytywać** pliki wejściowe przy użyciu metod obsługi plików Nextflow
2. **Wyodrębniać** metadane (ID pacjenta, replikat, typ próbki) ze struktury nazwy pliku
3. **Grupować** sparowane pliki (R1/R2) razem przy użyciu `channel.fromFilePairs()`
4. **Przetwarzać** pliki za pomocą dostarczonego modułu analizy
5. **Organizować** wyjścia w strukturę katalogów opartą na wyodrębnionych metadanych

#### Lista gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace działa
- [ ] Ustawiłem katalog roboczy odpowiednio
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Podstawowe operacje na plikach

### 1.1. Zidentyfikuj typ obiektu za pomocą `.class`

Spójrz na plik workflow `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Utwórz obiekt Path ze ścieżki tekstowej
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

To mini-workflow (bez żadnych procesów), który odnosi się do pojedynczej ścieżki pliku w Swoim workflow, następnie wypisuje ją do konsoli wraz z jej klasą.

??? info "Co to jest `.class`?"

    W Nextflow, `.class` mówi nam, z jakim typem obiektu pracujemy. To jak pytanie "co to za rzecz?" aby dowiedzieć się, czy to ciąg znaków, liczba, plik czy coś innego.
    To pomoże nam zilustrować różnicę między zwykłym ciągiem znaków a obiektem Path w następnych sekcjach.

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

Jak widać, Nextflow wypisał ścieżkę ciągu dokładnie tak, jak ją napisaliśmy.

To jest tylko wyjście tekstowe; Nextflow nie zrobił z tym jeszcze nic specjalnego.
Potwierdziliśmy również, że dla Nextflow jest to tylko ciąg znaków (klasy `java.lang.String`).
To ma sens, ponieważ nie powiedzieliśmy jeszcze Nextflow, że odpowiada to plikowi.

### 1.2. Utwórz obiekt Path za pomocą file()

Możemy powiedzieć Nextflow, jak obsługiwać pliki, tworząc [obiekty Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) z ciągów znaków ścieżek.

W naszym workflow możemy przekonwertować ciąg ścieżki `data/patientA_rep1_normal_R1_001.fastq.gz` na obiekt Path używając metody `file()`, która zapewnia dostęp do właściwości i operacji na plikach.

Edytuj `main.nf`, aby opakować ciąg za pomocą `file()` w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Utwórz obiekt Path ze ścieżki tekstowej
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Utwórz obiekt Path ze ścieżki tekstowej
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

Nextflow przekonwertował nasz ciąg na obiekt Path i rozwiązał go do rzeczywistej lokalizacji pliku w systemie.
Ścieżka pliku będzie teraz bezwzględna, jak w `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Zauważ również, że klasa obiektu Path to `sun.nio.fs.UnixPath`: to sposób Nextflow na reprezentowanie plików lokalnych.
Jak zobaczymy później, pliki zdalne będą miały inne nazwy klas (takie jak `nextflow.file.http.XPath` dla plików HTTP), ale wszystkie działają dokładnie w ten sam sposób i mogą być używane identycznie w Twoich workflow.

!!! tip

    **Kluczowa różnica:**

    - **Ciąg ścieżki**: Tylko tekst, który Nextflow traktuje jako znaki
    - **Obiekt Path**: Inteligentne odniesienie do pliku, z którym Nextflow może pracować

    Pomyśl o tym w ten sposób: ciąg ścieżki jest jak napisanie adresu na papierze, podczas gdy obiekt Path jest jak załadowanie adresu w urządzeniu GPS, które wie, jak tam dotrzeć i może podać szczegóły dotyczące trasy.

### 1.3. Dostęp do atrybutów pliku

Dlaczego to jest pomocne? Teraz, gdy Nextflow rozumie, że `myFile` jest obiektem Path, a nie tylko ciągiem znaków, możemy uzyskać dostęp do różnych atrybutów obiektu Path.

Zaktualizujmy nasz workflow, aby wypisał wbudowane atrybuty pliku:

=== "Po"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Utwórz obiekt Path ze ścieżki tekstowej
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
        // Utwórz obiekt Path ze ścieżki tekstowej
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

Widzisz różne atrybuty pliku wypisane w konsoli powyżej.

### 1.4. Przekaż plik do procesu

Różnica między ciągami znaków a obiektami Path staje się krytyczna, gdy zaczniesz budować rzeczywiste workflow z procesami.
Do tej pory zweryfikowaliśmy, że Nextflow traktuje teraz nasz plik wejściowy jako plik, ale zobaczmy, czy możemy faktycznie uruchomić coś na tym pliku w procesie.

#### 1.4.1. Zaimportuj proces i zbadaj kod

Udostępniamy wstępnie napisany moduł procesu o nazwie `COUNT_LINES`, który przyjmuje plik wejściowy i liczy, ile ma linii.

Aby użyć procesu w workflow, wystarczy dodać instrukcję include przed blokiem workflow:

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

Jak widać, to dość prosty mały skrypt, który rozpakowuje plik i liczy, ile zawiera linii.

??? info "Co robi `debug true`?"

    Dyrektywa `debug true` w definicji procesu powoduje, że Nextflow wypisuje wyjście ze skryptu (jak liczba linii "40") bezpośrednio w logu wykonania.
    Bez tego zobaczyłbyś tylko status wykonania procesu, ale nie rzeczywiste wyjście ze skryptu.

    Więcej informacji na temat debugowania workflow Nextflow znajdziesz w side queście [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Dodaj wywołanie `COUNT_LINES`

Teraz, gdy proces jest dostępny dla workflow, możemy dodać wywołanie procesu `COUNT_LINES`, aby uruchomić go na pliku wejściowym.

Wprowadź następujące edycje w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Utwórz obiekt Path ze ścieżki tekstowej
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
        // Utwórz obiekt Path ze ścieżki tekstowej
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

A teraz uruchom workflow:

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

Konkretnie, Nextflow wykonał następujące operacje pomyślnie:

- Przeniósł plik do katalogu roboczego
- Zdekompresował plik .gz
- Policzył linie (40 linii w tym przypadku)
- Zakończył bez błędu

Kluczem do tej płynnej operacji jest to, że wyraźnie mówimy Nextflow, że nasze wejście jest plikiem i powinno być traktowane jako takie.

### 1.5. Rozwiązywanie problemów z podstawowymi błędami wejścia pliku

To często wprowadza w błąd osoby nowe w Nextflow, więc poświęćmy kilka minut na przyjrzenie się, co się dzieje, gdy robisz to źle.

Są dwa główne miejsca, w których możesz źle obsłużyć plik: na poziomie workflow i na poziomie procesu.

#### 1.5.1. Błąd na poziomie workflow

Zobaczmy, co się stanie, jeśli wrócimy do traktowania pliku jako ciągu znaków podczas określania wejścia w bloku workflow.

Wprowadź następujące edycje w workflow, upewniając się, że zakomentowałeś instrukcje wypisywania specyficzne dla ścieżki:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Utwórz obiekt Path ze ścieżki tekstowej
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
        // Utwórz obiekt Path ze ścieżki tekstowej
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

A teraz uruchom workflow:

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

To jest ważna część:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Kiedy określasz wejście `path`, Nextflow waliduje, że przekazujesz rzeczywiste odniesienia do plików, a nie tylko ciągi znaków.
Ten błąd mówi Ci, że `'data/patientA_rep1_normal_R1_001.fastq.gz'` nie jest prawidłową wartością ścieżki, ponieważ jest to ciąg znaków, a nie obiekt Path.

Nextflow natychmiast wykrył problem i zatrzymał się przed rozpoczęciem procesu.

#### 1.5.2. Błąd na poziomie procesu

Drugim miejscem, w którym możemy zapomnieć określić, że chcemy, aby Nextflow traktował wejście jako plik, jest definicja procesu.

!!! warning "Zachowaj błąd workflow z 1.5.1"

    Aby ten test działał poprawnie, zachowaj workflow w jego uszkodzonym stanie (używając zwykłego ciągu zamiast `file()`).
    W połączeniu z `val` w procesie, to powoduje błąd pokazany poniżej.

Wprowadź następującą edycję w module:

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

A teraz uruchom workflow ponownie:

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

To pokazuje wiele szczegółów dotyczących błędu, ponieważ proces jest ustawiony na wyświetlanie informacji debugowania, jak wspomniano powyżej.

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

Mówi to, że system nie mógł znaleźć pliku; jednak jeśli sprawdzisz ścieżkę, istnieje plik o tej nazwie w tej lokalizacji.

Kiedy to uruchomiliśmy, Nextflow przekazał wartość ciągu do skryptu, ale nie _przeniósł_ rzeczywistego pliku do katalogu roboczego.
Więc proces próbował użyć względnego ciągu, `data/patientA_rep1_normal_R1_001.fastq.gz`, ale ten plik nie istnieje w katalogu roboczym procesu.

Razem wzięte, te dwa przykłady pokazują, jak ważne jest poinformowanie Nextflow, czy wejście powinno być obsługiwane jako plik.

!!! note

    Upewnij się, że cofniesz i naprawisz oba celowe błędy przed kontynuowaniem następnej sekcji.

### Wnioski

- Ciągi ścieżek vs obiekty Path: Ciągi to tylko tekst, obiekty Path to inteligentne odniesienia do plików
- Metoda `file()` konwertuje ciąg ścieżki na obiekt Path, z którym Nextflow może pracować
- Możesz uzyskać dostęp do właściwości pliku, takich jak `name`, `simpleName`, `extension` i `parent` [używając atrybutów pliku](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Używanie obiektów Path zamiast ciągów pozwala Nextflow prawidłowo zarządzać plikami w Twoim workflow
- Wyniki wejścia procesu: Prawidłowa obsługa plików wymaga obiektów Path, a nie ciągów, aby zapewnić, że pliki są prawidłowo przenoszone i dostępne do użycia przez procesy.

---

## 2. Używanie plików zdalnych

Jedną z kluczowych funkcji Nextflow jest możliwość płynnego przełączania między plikami lokalnymi (na tej samej maszynie) a plikami zdalnymi dostępnymi przez internet.

Jeśli robisz to dobrze, nigdy nie powinieneś potrzebować zmieniać logiki Swojego workflow, aby obsługiwać pliki pochodzące z różnych lokalizacji.
Wszystko, co musisz zrobić, aby użyć pliku zdalnego, to określić odpowiedni prefiks w ścieżce pliku podczas dostarczania go do workflow.

Na przykład, `/path/to/data` nie ma prefiksu, co wskazuje, że jest to 'normalna' lokalna ścieżka pliku, podczas gdy `s3://path/to/data` zawiera prefiks `s3://`, wskazując, że znajduje się w magazynie obiektów S3 Amazon.

Obsługiwanych jest wiele różnych protokołów:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Aby użyć któregokolwiek z nich, po prostu określ odpowiedni prefiks w ciągu, który jest następnie technicznie nazywany identyfikatorem URI (Uniform Resource Identifier) zamiast ścieżką pliku.
Nextflow zajmie się uwierzytelnianiem i przeniesieniem plików we właściwe miejsce, pobieraniem lub przesyłaniem oraz wszystkimi innymi operacjami na plikach, których można oczekiwać.

Kluczową zaletą tego systemu jest to, że umożliwia nam przełączanie między środowiskami bez zmiany jakiejkolwiek logiki pipeline.
Na przykład możesz rozwijać z małym, lokalnym zestawem testowym przed przełączeniem na pełnoskalowy zestaw testowy znajdujący się w zdalnym magazynie, po prostu zmieniając URI.

### 2.1. Użyj pliku z internetu

Przetestujmy to, zamieniając lokalną ścieżkę, którą dostarczamy do naszego workflow, na ścieżkę HTTPS wskazującą na kopię tych samych danych przechowywanych w Github.

!!! warning

    To będzie działać tylko wtedy, gdy masz aktywne połączenie internetowe.

Otwórz `main.nf` ponownie i zmień ścieżkę wejściową w następujący sposób:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Użycie zdalnego pliku z internetu
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
        // Utwórz obiekt Path ze ścieżki tekstowej
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

Działa! Możesz zobaczyć, że niewiele się zmieniło.

Jedną różnicą w wyjściu konsoli jest to, że klasa obiektu ścieżki to teraz `nextflow.file.http.XPath`, podczas gdy dla ścieżki lokalnej klasa była `sun.nio.fs.UnixPath`.
Nie musisz pamiętać tych klas; wspominamy to tylko po to, aby pokazać, że Nextflow identyfikuje i obsługuje różne lokalizacje odpowiednio.

Za kulisami Nextflow pobrał plik do katalogu tymczasowego znajdującego się w katalogu roboczym.
Ten przygotowany plik może być następnie traktowany jako plik lokalny i dowiązany symbolicznie do odpowiedniego katalogu procesu.

Możesz zweryfikować, że tak się stało, przeglądając zawartość katalogu roboczego znajdującego się pod wartością hash procesu.

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

    Dowiązanie symboliczne wskazuje na kopię przygotowaną pliku zdalnego, którą Nextflow automatycznie pobrał.

Zauważ, że dla większych plików krok pobierania zajmie trochę więcej czasu w porównaniu z uruchomieniem na plikach lokalnych.
Jednak Nextflow sprawdza, czy już ma przygotowaną kopię, aby uniknąć niepotrzebnych pobrań.
Więc jeśli uruchomisz ponownie na tym samym pliku i nie usunąłeś przygotowanego pliku, Nextflow użyje przygotowanej kopii.

To pokazuje, jak łatwo jest przełączać między danymi lokalnymi i zdalnymi za pomocą Nextflow, co jest kluczową funkcją Nextflow.

!!! note

    Jednym ważnym wyjątkiem od tej zasady jest to, że nie możesz używać wzorców glob ani ścieżek katalogów z HTTPS, ponieważ HTTPS nie może wymieniać wielu plików, więc musisz określić dokładne adresy URL plików.
    Jednak inne protokoły magazynowania, takie jak magazyn obiektów blob (`s3://`, `az://`, `gs://`) mogą używać zarówno globów, jak i ścieżek katalogów.

    Oto jak możesz używać wzorców glob z magazynem w chmurze:

    ```groovy title="Przykłady magazynu w chmurze (nie do uruchomienia w tym środowisku)"
    // S3 with glob patterns - would match multiple files
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage ze wzorcami glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage ze wzorcami glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Pokażemy Ci, jak pracować z globami w praktyce w następnej sekcji.

### 2.2. Wróć do pliku lokalnego

Wrócimy do używania naszych lokalnych przykładowych plików przez resztę tego side questa, więc przełączmy wejście workflow z powrotem na oryginalny plik:

=== "Po"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Utwórz obiekt Path ze ścieżki tekstowej
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
        // Utwórz obiekt Path ze ścieżki tekstowej
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Wypisz atrybuty pliku
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Wnioski

- Dostęp do danych zdalnych odbywa się przy użyciu URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow automatycznie pobierze i przeniesie dane we właściwe miejsce, o ile te ścieżki są przekazywane do procesów
- Nie pisz logiki do pobierania lub przesyłania plików zdalnych!
- Pliki lokalne i zdalne produkują różne typy obiektów, ale działają identycznie
- **Ważne**: HTTP/HTTPS działają tylko z pojedynczymi plikami (bez wzorców glob)
- Magazyn w chmurze (S3, Azure, GCS) obsługuje zarówno pojedyncze pliki, jak i wzorce glob
- Możesz bezproblemowo przełączać między lokalnymi i zdalnymi źródłami danych bez zmiany logiki kodu (o ile protokół obsługuje wymagane operacje)

---

## 3. Używanie fabryki kanałów `fromPath()`

Do tej pory pracowaliśmy z jednym plikiem na raz, ale w Nextflow zazwyczaj będziemy chcieli utworzyć kanał wejściowy z wieloma plikami wejściowymi do przetworzenia.

Naiwnym sposobem byłoby połączenie metody `file()` z [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) w ten sposób:

```groovy title="Przykład składni"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

To działa, ale jest nieporęczne.

!!! tip "Kiedy używać `file()` vs `channel.fromPath()`"

    - Użyj `file()`, gdy potrzebujesz pojedynczego obiektu Path do bezpośredniej manipulacji (sprawdzanie, czy plik istnieje, odczytywanie jego atrybutów lub przekazywanie do pojedynczego wywołania procesu)
    - Użyj `channel.fromPath()`, gdy potrzebujesz kanału, który może przechowywać wiele plików, szczególnie ze wzorcami glob, lub gdy pliki będą przepływać przez wiele procesów

Tu wchodzi [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): wygodna fabryka kanałów, która łączy wszystkie funkcje potrzebne do wygenerowania kanału z jednego lub więcej statycznych ciągów plików oraz wzorców glob.

### 3.1. Dodaj fabrykę kanałów

Zaktualizujmy nasz workflow, aby używał `channel.fromPath`.

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Wypisz atrybuty pliku
        /* Comment these out for now, we'll come back to them!
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
        // Utwórz obiekt Path ze ścieżki tekstowej
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

Zakomentowaliśmy również kod wypisujący atrybuty na razie i dodaliśmy instrukcję `.view` do wypisania tylko nazwy pliku.

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

Jak widać, ścieżka pliku jest ładowana jako obiekt typu `Path` w kanale.
Jest to podobne do tego, co zrobiłby `file()`, z wyjątkiem tego, że teraz mamy kanał, do którego możemy załadować więcej plików, jeśli chcemy.

Używanie `channel.fromPath()` to wygodny sposób tworzenia nowego kanału wypełnionego listą plików.

### 3.2. Wyświetl atrybuty plików w kanale

W naszej pierwszej próbie użycia fabryki kanałów uprościliśmy kod i po prostu wypisaliśmy nazwę pliku.

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

Ponownie włączamy również wywołanie procesu `COUNT_LINES`, aby sprawdzić, czy przetwarzanie plików nadal działa poprawnie z naszym podejściem opartym na kanałach.

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

I masz to, te same wyniki co wcześniej, ale teraz mamy plik w kanale, więc możemy dodać więcej.

### 3.3. Używanie globu do dopasowywania wielu plików

Jest kilka sposobów, w jakie moglibyśmy załadować więcej plików do kanału.
Tutaj pokażemy Ci, jak używać wzorców glob, które są wygodnym sposobem dopasowywania i pobierania nazw plików i katalogów na podstawie znaków wieloznacznych.
Proces dopasowywania tych wzorców nazywa się "globbingiem" lub "rozszerzaniem nazw plików".

!!! note

    Jak wspomniano wcześniej, Nextflow obsługuje globbing do zarządzania plikami wejściowymi i wyjściowymi w większości przypadków, z wyjątkiem ścieżek plików HTTPS, ponieważ HTTPS nie może wymieniać wielu plików.

Powiedzmy, że chcemy pobrać oba pliki w parze plików związanych z danym pacjentem, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Ponieważ jedyną różnicą między nazwami plików jest numer replikatu, _czyli_ liczba po `R`, możemy użyć znaku wieloznacznego `*` do zastąpienia liczby w następujący sposób:

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

Jak widać, mamy teraz dwa obiekty Path w naszym kanale, co pokazuje, że Nextflow poprawnie wykonał rozszerzenie nazw plików i załadował i przetworzył oba pliki zgodnie z oczekiwaniami.

Używając tej metody, możemy pobrać tyle lub tak niewiele plików, ile chcemy, po prostu zmieniając wzorzec glob. Gdybyśmy uczynili go bardziej ogólnym, na przykład zastępując wszystkie zmienne części nazw plików przez `*` (_np._ `data/patient*_rep*_*_R*_001.fastq.gz`) moglibyśmy pobrać wszystkie przykładowe pliki w katalogu `data`.

### Wnioski

- `channel.fromPath()` tworzy kanał z plikami pasującymi do wzorca
- Każdy plik jest emitowany jako oddzielny element w kanale
- Możemy użyć wzorca glob do dopasowania wielu plików
- Pliki są automatycznie konwertowane na obiekty Path z pełnymi atrybutami
- Metoda `.view()` umożliwia inspekcję zawartości kanału

---

## 4. Wyodrębnianie podstawowych metadanych z nazw plików

W większości dziedzin naukowych bardzo powszechne jest kodowanie metadanych w nazwach plików zawierających dane.
Na przykład w bioinformatyce pliki zawierające dane sekwencjonowania są często nazywane w sposób kodujący informacje o próbce, warunku, replikacie i numerze odczytu.

Jeśli nazwy plików są konstruowane zgodnie z spójną konwencją, możesz wyodrębnić te metadane w znormalizowany sposób i użyć ich w trakcie analizy.
To jest duże „jeśli", oczywiście, i powinieneś być bardzo ostrożny, ilekroć polegasz na strukturze nazw plików; ale rzeczywistość jest taka, że to podejście jest bardzo szeroko stosowane, więc przyjrzyjmy się, jak to się robi w Nextflow.

W przypadku naszych przykładowych danych wiemy, że nazwy plików zawierają konsekwentnie ustrukturyzowane metadane.
Na przykład nazwa pliku `patientA_rep1_normal_R2_001` koduje następujące:

- ID pacjenta: `patientA`
- ID replikatu: `rep1`
- typ próbki: `normal` (w przeciwieństwie do `tumor`)
- zestaw odczytów: `R1` (w przeciwieństwie do `R2`)

Zmodyfikujemy nasz workflow, aby pobrać te informacje w trzech krokach:

1. Pobierz `simpleName` pliku, który zawiera metadane
2. Rozdziel metadane za pomocą metody zwanej `tokenize()`
3. Użyj mapy do zorganizowania metadanych

!!! warning

    Nigdy nie powinieneś kodować wrażliwych informacji w nazwach plików, takich jak nazwiska pacjentów lub inne cechy identyfikujące, ponieważ może to zagrozić prywatności pacjentów lub innym odpowiednim ograniczeniom bezpieczeństwa.

### 4.1. Pobierz `simpleName`

`simpleName` jest atrybutem pliku, który odpowiada nazwie pliku pozbawionej ścieżki i rozszerzenia.

Wprowadź następujące edycje w workflow:

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

To pobiera `simpleName` i kojarzy go z pełnym obiektem pliku za pomocą operacji `map()`.

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

### 4.2. Wyodrębnij metadane z `simplename`

W tym momencie metadane, które chcemy, są osadzone w `simplename`, ale nie możemy bezpośrednio uzyskać dostępu do poszczególnych elementów.
Więc musimy podzielić `simplename` na jego komponenty.
Na szczęście te komponenty są po prostu oddzielone podkreśleniami w oryginalnej nazwie pliku, więc możemy zastosować powszechną metodę Nextflow zwaną `tokenize()`, która jest idealna do tego zadania.

Wprowadź następujące edycje w workflow:

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

Metoda `tokenize()` podzieli ciąg `simpleName` wszędzie tam, gdzie znajdzie podkreślenia, i zwróci listę zawierającą podciągi.

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

Teraz krotka dla każdego elementu w naszym kanale zawiera listę metadanych (_np._ `[patientA, rep1, normal, R1, 001]`) i oryginalny obiekt pliku.

To świetnie!
Rozbiliśmy informacje o pacjencie z pojedynczego ciągu na listę ciągów.
Możemy teraz obsługiwać każdą część informacji o pacjencie osobno.

### 4.3. Użyj mapy do zorganizowania metadanych

Nasze metadane to obecnie po prostu płaska lista.
Łatwo ją użyć, ale trudno odczytać.

```console
[patientA, rep1, normal, R1, 001]
```

Co to jest element o indeksie 3? Czy możesz powiedzieć bez odniesienia się do oryginalnego wyjaśnienia struktury metadanych?

To świetna okazja do użycia magazynu klucz-wartość, gdzie każdy element ma zestaw kluczy i powiązanych z nimi wartości, więc możesz łatwo odwoływać się do każdego klucza, aby uzyskać odpowiednią wartość.

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

W Nextflow nazywa się to [map](https://nextflow.io/docs/latest/script.html#maps).

Przekonwertujmy teraz naszą płaską listę na mapę.
Wprowadź następujące edycje w workflow:

=== "Po"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Załaduj pliki za pomocą channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, read
