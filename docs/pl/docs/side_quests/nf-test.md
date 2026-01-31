# Testowanie z nf-test

Możliwość systematycznego testowania, czy każda część przepływu pracy działa zgodnie z założeniami, jest kluczowa dla odtwarzalności i długoterminowej konserwacji, a także może być ogromną pomocą podczas procesu twórczego.

Poświęćmy chwilę na omówienie, dlaczego testowanie jest tak ważne. Jeśli tworzysz przepływ pracy, jedną z pierwszych rzeczy, które zrobisz, jest pobranie danych testowych, o których wiesz, że są poprawne i powinny dać wynik. Dodajesz pierwszy proces do potoku i podłączasz go do swoich danych wejściowych, aby działał. Następnie, aby sprawdzić, czy wszystko działa, uruchamiasz go na danych testowych. Zakładając, że działa, przechodzisz do następnego procesu i ponownie uruchamiasz dane testowe. Powtarzasz ten proces, aż otrzymasz potok, z którego jesteś zadowolony.

Następnie, być może dodajesz prosty parametr prawda lub fałsz, taki jak `--skip_process`. Teraz musisz uruchomić potok dwa razy, raz z każdym parametrem, aby upewnić się, że działa zgodnie z oczekiwaniami. Ale czekaj, jak sprawdzić, czy `--skip_process` faktycznie pomija proces? Musimy przejrzeć wyniki lub sprawdzić pliki logów! To jest uciążliwe i podatne na błędy.

W miarę rozwoju potoku szybko stanie się on tak złożony, że ręczne testowanie każdej iteracji będzie wolne i podatne na błędy. Co więcej, jeśli znajdziesz błąd, bardzo trudno będzie dokładnie określić, skąd w potoku pochodzi błąd. Właśnie tutaj przydaje się testowanie.

Testowanie pozwala systematycznie sprawdzać, czy każda część potoku działa zgodnie z oczekiwaniami. Korzyści dla programisty z dobrze napisanych testów są ogromne:

- **Pewność**: Ponieważ testy obejmują cały potok, możesz być pewny, że zmiana czegoś nie wpłynie na nic innego
- **Zaufanie**: Kiedy wielu programistów pracuje nad potokiem, wiedzą, że inni programiści nie zepsuli potoku ani żadnego komponentu.
- **Przejrzystość**: Testy pokazują, gdzie potok zawodzi i ułatwiają wyśledzenie problemu. Działają również jako forma dokumentacji, pokazując, jak uruchomić proces lub przepływ pracy.
- **Szybkość**: Ponieważ testy są zautomatyzowane, można je uruchamiać bardzo szybko i wielokrotnie. Możesz szybko iterować z mniejszą obawy o wprowadzenie nowych błędów.

Możemy napisać wiele różnych rodzajów testów:

1. **Testy na poziomie modułu**: Dla pojedynczych procesów
2. **Testy na poziomie przepływu pracy**: Dla pojedynczego przepływu pracy
3. **Testy na poziomie potoku**: Dla potoku jako całości
4. **Testy wydajnościowe**: Dla szybkości i wydajności potoku
5. **Testy obciążeniowe**: Ocena wydajności potoku w ekstremalnych warunkach w celu określenia jego granic

Testowanie pojedynczych procesów jest analogiczne do testów jednostkowych w innych językach. Testowanie przepływu pracy lub całego potoku jest analogiczne do tego, co nazywa się testami integracyjnymi w innych językach, gdzie testujemy interakcje komponentów.

[**nf-test**](https://www.nf-test.com/) to narzędzie, które pozwala pisać testy na poziomie modułu, przepływu pracy i potoku. Krótko mówiąc, pozwala systematycznie sprawdzać, czy każda pojedyncza część potoku działa zgodnie z oczekiwaniami, _w izolacji_.

### Cele nauki

W tym side queście nauczysz się używać nf-test do pisania testów na poziomie przepływu pracy dla potoku, a także testów na poziomie modułu dla trzech procesów, które wywołuje.

Pod koniec tego side questa będziesz w stanie efektywnie używać następujących technik:

- Inicjalizowanie nf-test w projekcie
- Generowanie testów na poziomie modułu i przepływu pracy
- Dodawanie typowych rodzajów asercji
- Rozumienie, kiedy używać migawek a asercji zawartości
- Uruchamianie testów dla całego projektu

Te umiejętności pomogą Ci wdrożyć kompleksową strategię testowania w projektach potoków, zapewniając, że są bardziej solidne i łatwe w utrzymaniu.

### Wymagania wstępne

Przed podjęciem tego side questa powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory, praca z plikami, metadane)

---

## 0. Rozpoczęcie pracy

#### Otwórz codespace szkoleniowy

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzyłeś środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki dla tego samouczka.

```bash
cd side-quests/nf-test
```

Możesz ustawić VSCode tak, aby skupił się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz główny plik przepływu pracy i plik CSV o nazwie `greetings.csv`, który zawiera dane wejściowe do potoku.

```console title="Zawartość katalogu"
.
├── greetings.csv
└── main.nf
```

Aby uzyskać szczegółowy opis plików, zobacz [rozgrzewkę z Hello Nextflow](../hello_nextflow/00_orientation.md).

Przepływ pracy, który będziemy testować, jest podzbiorem przepływu pracy Hello zbudowanego w [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "Co robi przepływ pracy Hello Nextflow?"

    Jeśli nie ukończyłeś szkolenia [Hello Nextflow](../hello_nextflow/index.md), oto krótki przegląd tego, co robi ten prosty przepływ pracy.

    Przepływ pracy pobiera plik CSV zawierający pozdrowienia, przeprowadza na nich cztery kolejne kroki transformacji i wyprowadza pojedynczy plik tekstowy zawierający obraz ASCII zabawnej postaci wypowiadającej pozdrowienia.

    Cztery kroki są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w oddzielnych plikach modułów.

    1. **`sayHello`:** Zapisuje każde pozdrowienie do własnego pliku wyjściowego (np. "Hello-output.txt")
    2. **`convertToUpper`:** Konwertuje każde pozdrowienie na wielkie litery (np. "HELLO")
    3. **`collectGreetings`:** Zbiera wszystkie pozdrowienia z wielkimi literami do jednego pliku wsadowego
    4. **`cowpy`:** Generuje grafikę ASCII za pomocą narzędzia `cowpy`

    Wyniki są publikowane w katalogu o nazwie `results/`, a końcowe wyjście potoku (po uruchomieniu z domyślnymi parametrami) to zwykły plik tekstowy zawierający grafikę ASCII postaci wypowiadającej pozdrowienia z wielkimi literami.

    W tym side queście używamy pośredniej formy przepływu pracy Hello, która zawiera tylko dwa pierwsze procesy.

Podzbiór, z którym będziemy pracować, składa się z dwóch procesów: `sayHello` i `convertToUpper`.
Pełny kod przepływu pracy możesz zobaczyć poniżej.

??? example "Kod przepływu pracy"

    ```groovy title="main.nf"
    /*
    * Pipeline parameters
    */
    params.input_file = "greetings.csv"

    /*
    * Use echo to print 'Hello World!' to standard out
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Use a text replace utility to convert the greeting to uppercase
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### Uruchom przepływ pracy

Uruchommy przepływ pracy, aby upewnić się, że działa zgodnie z oczekiwaniami.

```bash
nextflow run main.nf
```

```console title="Wynik uruchomienia przepływu pracy"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

GRATULACJE! Właśnie uruchomiłeś test!

"Czekaj, co? Po prostu uruchomiłem przepływ pracy i zadziałał! Jak to jest test?"

Dobre pytanie!

Rozłóżmy to na czynniki pierwsze.

Uruchomiłeś przepływ pracy z domyślnymi parametrami, potwierdziłeś, że działa i jesteś zadowolony z wyników. To jest istota testowania. Jeśli pracowałeś przez kurs szkoleniowy Hello Nextflow, zauważysz, że zawsze zaczynaliśmy każdą sekcję od uruchomienia przepływu pracy, którego używaliśmy jako punktu wyjścia, aby potwierdzić, że wszystko jest poprawnie skonfigurowane.

Testowanie oprogramowania zasadniczo wykonuje ten proces za nas.

#### Przejrzyj zadanie

Twoim wyzwaniem jest dodanie standardowych testów do tego przepływu pracy za pomocą nf-test, aby ułatwić weryfikację, że każda część nadal działa zgodnie z oczekiwaniami w przypadku wprowadzenia jakichkolwiek dalszych zmian.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby się zanurzyć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace jest uruchomiony
- [ ] Odpowiednio ustawiłem katalog roboczy
- [ ] Pomyślnie uruchomiłem przepływ pracy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz rozpocząć.

---

## 1. Inicjalizacja `nf-test`

Pakiet `nf-test` zapewnia polecenie inicjalizacji, które konfiguruje kilka rzeczy, aby mogliśmy zacząć tworzyć testy dla naszego projektu.

```bash
nf-test init
```

Powinno to dać następujące wyjście:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Tworzy również katalog `tests` zawierający szkielet pliku konfiguracyjnego.

### 1.1. Wygeneruj nf-test

`nf-test` zawiera zestaw narzędzi do budowania plików nf-test, oszczędzając nam większość pracy. Są one dostępne pod podpoleceniem `generate`. Wygenerujmy test dla potoku:

```bash
nf-test generate pipeline main.nf
```

```console title="Wyjście"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

To utworzy plik `main.nf.test` w katalogu `tests`. To jest nasz plik testowy na poziomie potoku. Jeśli uruchomisz `tree tests/`, powinieneś zobaczyć coś takiego:

```console title="Zawartość katalogu testów"
tests/
├── main.nf.test
└── nextflow.config
```

Plik `main.nf.test` to nasz plik testowy na poziomie potoku. Otwórzmy go i przyjrzyjmy się zawartości.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Poświęćmy chwilę na zrozumienie struktury pliku testowego.

Blok `nextflow_pipeline` jest punktem wejścia dla wszystkich testów na poziomie potoku. Zawiera następujące elementy:

- `name`: Nazwa testu.
- `script`: Ścieżka do skryptu potoku.

Blok `test` to rzeczywisty test. Zawiera następujące elementy:

- `when`: Warunki, w których test powinien być uruchomiony. Obejmuje to parametry, które będą używane do uruchomienia potoku.
- `then`: Asercje, które powinny być wykonane. Obejmuje to oczekiwane wyniki potoku.

Mówiąc prostym językiem, logika testu brzmi następująco:
"**Gdy** te _parametry_ są dostarczane do tego _potoku_, **wtedy** oczekujemy zobaczyć te wyniki."

To nie jest test funkcjonalny, pokażemy, jak go przekształcić w taki w następnej sekcji.

### Uwaga o nazewnictwie testów

W powyższym przykładzie użyliśmy domyślnej nazwy "Should run without failures", która jest odpowiednia dla podstawowego testu, który tylko sprawdza, czy potok działa pomyślnie. Jednak gdy dodajemy bardziej szczegółowe przypadki testowe, powinniśmy używać bardziej opisowych nazw, które wskazują, co faktycznie testujemy. Na przykład:

- "Should convert input to uppercase" - podczas testowania konkretnej funkcjonalności
- "Should handle empty input gracefully" - podczas testowania przypadków brzegowych
- "Should respect max memory parameter" - podczas testowania ograniczeń zasobów
- "Should create expected output files" - podczas testowania generowania plików

Dobre nazwy testów powinny:

1. Zaczynać się od "Should", aby było jasne, jakie jest oczekiwane zachowanie
2. Opisywać konkretną funkcjonalność lub scenariusz, który jest testowany
3. Być wystarczająco jasne, aby w przypadku niepowodzenia testu wiedzieć, jaka funkcjonalność jest uszkodzona

Gdy później dodamy więcej asercji i konkretnych przypadków testowych, użyjemy tych bardziej opisowych nazw, aby było jasne, co weryfikuje każdy test.

### 1.2. Uruchom test

Uruchommy test, aby zobaczyć, co się stanie.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

Test nie powiódł się! Co się stało?

1. nf-test próbował uruchomić potok w obecnej postaci, używając ustawień w bloku `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test sprawdził status potoku i porównał go z blokiem `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Zauważ, jak nf-test zgłosił, że potok nie powiódł się i dostarczył komunikat o błędzie z Nextflow:

```console title="Błąd"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Jaki był więc problem? Pamiętaj, że potok ma plik greetings.csv w katalogu projektu. Gdy nf-test uruchamia potok, będzie szukał tego pliku, ale nie może go znaleźć. Plik jest tam, co się dzieje? Cóż, jeśli spojrzymy na ścieżkę, możemy zobaczyć, że test odbywa się w ścieżce `./nf-test/tests/longHashString/`. Podobnie jak Nextflow, nf-test tworzy nowy katalog dla każdego testu, aby wszystko było izolowane. Plik danych nie znajduje się tam, więc musimy poprawić ścieżkę do pliku w oryginalnym teście.

Wróćmy do pliku testowego i zmieńmy ścieżkę do pliku w bloku `when`.

Możesz się zastanawiać, jak będziemy wskazywać na katalog główny potoku w teście. Ponieważ jest to powszechna sytuacja, nf-test ma szereg zmiennych globalnych, których możemy użyć, aby ułatwić sobie życie. Pełną listę znajdziesz [tutaj](https://www.nf-test.com/docs/testcases/global_variables/), ale tymczasem użyjemy zmiennej `projectDir`, co oznacza katalog główny projektu potoku.

_Przed:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Po:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Uruchommy test ponownie, aby zobaczyć, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Potok przechodzi"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Sukces! Potok działa pomyślnie i test przechodzi. Uruchom go tyle razy, ile chcesz, a zawsze otrzymasz ten sam wynik!

Domyślnie wyjście Nextflow jest ukryte, ale aby przekonać się, że nf-test na pewno uruchamia przepływ pracy, możesz użyć flagi `--verbose`:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Potok uruchamia wszystkie procesy"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Dodaj asercje

Prostym sprawdzeniem jest upewnienie się, że nasz potok uruchamia wszystkie oczekiwane procesy i nie pomija żadnego w sposób cichy. Pamiętaj, że nasz potok uruchamia 6 procesów, jeden o nazwie `sayHello` i jeden o nazwie `convertToUpper` dla każdego z 3 pozdrowień.

Dodajmy asercję do naszego testu, aby sprawdzić, czy potok uruchamia oczekiwaną liczbę procesów. Zaktualizujemy również nazwę naszego testu, aby lepiej odzwierciedlała to, co testujemy.

**Przed:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**Po:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

Nazwa testu teraz lepiej odzwierciedla to, co faktycznie weryfikujemy - nie tylko to, że potok działa bez awarii, ale że uruchamia oczekiwaną liczbę procesów.

Uruchommy test ponownie, aby zobaczyć, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Potok przechodzi z asercjami"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Sukces! Potok działa pomyślnie i test przechodzi. Teraz zaczęliśmy testować szczegóły potoku, a także ogólny status.

### 1.4. Testowanie wyjścia

Dodajmy asercję do naszego testu, aby sprawdzić, czy plik wyjściowy został utworzony. Dodamy ją jako oddzielny test z informacyjną nazwą, aby ułatwić interpretację wyników.

**Przed:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**Po:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

Uruchom test ponownie, aby zobaczyć, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Potok przechodzi z asercjami plików"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Sukces! Testy przechodzą, ponieważ potok zakończył się pomyślnie, uruchomiona została prawidłowa liczba procesów, a pliki wyjściowe zostały utworzone. Powinno to również pokazać, jak przydatne jest dostarczanie tych informacyjnych nazw dla testów.

To tylko powierzchnia, możemy dalej pisać asercje, aby sprawdzić szczegóły potoku, ale na razie przejdźmy do testowania wewnętrznych elementów potoku.

### Wnioski

Wiesz, jak napisać nf-test dla potoku.

### Co dalej?

Naucz się testować proces Nextflow.

---

## 2. Testowanie procesu Nextflow

Nie musimy pisać testów dla każdej części potoku, ale im więcej testów mamy, tym bardziej kompleksowo możemy ocenić potok i tym bardziej możemy być pewni, że działa zgodnie z oczekiwaniami. W tej sekcji będziemy testować oba procesy w potoku jako pojedyncze jednostki.

### 2.1. Testowanie procesu `sayHello`

Zacznijmy od procesu `sayHello`.

Użyjmy ponownie polecenia `nf-test generate`, aby wygenerować testy dla procesu.

```bash
nf-test generate process main.nf
```

```console title="Wyjście"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Skupmy się teraz na procesie `sayhello` w pliku `main.sayhello.nf.test`.

Otwórzmy plik i przyjrzyjmy się zawartości.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Jak poprzednio, zaczynamy od szczegółów testu, po których następują bloki `when` i `then`. Jednak mamy również dodatkowy blok `process`, który pozwala nam zdefiniować wejścia do procesu.

Uruchommy test, aby zobaczyć, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Test procesu nie powiódł się"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

Test nie powiódł się, ponieważ proces `sayHello` deklaruje 1 wejście, ale został wywołany z 0 argumentami. Naprawmy to, dodając wejście do procesu. Pamiętaj z [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (i sekcji rozgrzewkowej powyżej), że nasz proces `sayHello` przyjmuje pojedyncze wejście wartości, które musimy podać. Powinniśmy również poprawić nazwę testu, aby lepiej odzwierciedlała to, co testujemy.

**Przed:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Po:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Uruchommy test ponownie, aby zobaczyć, czy działa.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Sukces! Test przechodzi, ponieważ proces `sayHello` działał pomyślnie i wyjście zostało utworzone.

### 2.2. Sprawdź migawkę utworzoną przez test

Jeśli spojrzymy na plik `tests/main.sayhello.nf.test`, zobaczymy, że używa metody `snapshot()` w bloku asercji:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

To mówi nf-test, aby utworzyć migawkę wyjścia procesu `sayHello`. Przyjrzyjmy się zawartości pliku migawki.

```console title="Zawartość pliku migawki"
code tests/main.sayhello.nf.test.snap
```

Nie wydrukujemy tego tutaj, ale powinieneś zobaczyć plik JSON zawierający szczegóły procesu i wyjść procesu. W szczególności możemy zobaczyć linię, która wygląda następująco:

```json title="Zawartość pliku migawki"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

To reprezentuje wyjścia utworzone przez proces `sayHello`, które testujemy jawnie. Jeśli ponownie uruchomimy test, program sprawdzi, czy nowe wyjście pasuje do wyjścia, które zostało pierwotnie zarejestrowane. To jest szybki, prosty sposób testowania, że wyjścia procesu się nie zmieniają, dlatego nf-test zapewnia to jako wartość domyślną.

!!!warning

    To oznacza, że musimy być pewni, że wyjście, które rejestrujemy w oryginalnym uruchomieniu, jest poprawne!

Jeśli w trakcie przyszłego rozwoju coś w kodzie się zmieni, co spowoduje, że wyjście będzie inne, test nie powiedzie się i będziemy musieli określić, czy zmiana jest oczekiwana, czy nie.

- Jeśli okaże się, że coś w kodzie się zepsuło, będziemy musieli to naprawić, oczekując, że naprawiony kod przejdzie test.
- Jeśli jest to oczekiwana zmiana (np. narzędzie zostało ulepszone, a wyniki są lepsze), będziemy musieli zaktualizować migawkę, aby zaakceptować nowe wyjście jako odniesienie do dopasowania. nf-test ma parametr `--update-snapshot` w tym celu.

Możemy uruchomić test ponownie i zobaczyć, że test powinien przejść:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Sukces! Test przechodzi, ponieważ proces `sayHello` działał pomyślnie, a wyjście pasowało do migawki.

### 2.3. Alternatywa dla migawek: bezpośrednie asercje zawartości

Chociaż migawki są świetne do wychwytywania wszelkich zmian w wyjściu, czasami chcesz zweryfikować konkretną zawartość bez bycia tak restrykcyjnym co do całego dopasowania pliku. Na przykład:

- Gdy części wyjścia mogą się zmieniać (znaczniki czasu, losowe identyfikatory itp.), ale musi być obecna określona kluczowa zawartość
- Gdy chcesz sprawdzić konkretne wzorce lub wartości w wyjściu
- Gdy chcesz uczynić test bardziej jawnym co do tego, co stanowi sukces

Oto jak moglibyśmy zmodyfikować nasz test, aby sprawdzić konkretną zawartość:

**Przed:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Po:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

Zauważ, że nf-test widzi wyjścia procesu jako listę list, więc `process.out[0][0]` pobiera pierwszą część pierwszego elementu kanału (lub 'emisji') z tego procesu.

To podejście:

- Wyraźnie pokazuje, czego dokładnie oczekujemy w wyjściu
- Jest bardziej odporne na nieistotne zmiany w wyjściu
- Zapewnia lepsze komunikaty o błędach, gdy testy nie przechodzą
- Pozwala na bardziej złożone walidacje (wzorce wyrażeń regularnych, porównania numeryczne itp.)

Uruchommy test, aby zobaczyć, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Test procesu nie powiódł się"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Testowanie procesu `convertToUpper`

Otwórzmy plik `tests/main.converttoupper.nf.test` i przyjrzyjmy się zawartości:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

To jest podobny test do procesu `sayHello`, ale testuje proces `convertToUpper`. Wiemy, że ten nie powiedzie się, ponieważ podobnie jak z `sayHello`, proces `convertToUpper` przyjmuje pojedyncze wejście ścieżki, ale go nie określiliśmy.

Teraz musimy dostarczyć pojedynczy plik wejściowy do procesu convertToUpper, który zawiera tekst, który chcemy przekonwertować na wielkie litery. Jest wiele sposobów, w jakie moglibyśmy to zrobić:

- Moglibyśmy utworzyć dedykowany plik do testowania
- Moglibyśmy ponownie użyć istniejącego pliku data/greetings.csv
- Moglibyśmy utworzyć go w locie w teście

Na razie ponownie użyjmy istniejącego pliku data/greetings.csv, używając przykładu, którego użyliśmy w teście na poziomie potoku. Jak poprzednio, możemy nazwać test, aby lepiej odzwierciedlał to, co testujemy, ale tym razem pozostawmy "migawkę" zawartości zamiast sprawdzania konkretnych ciągów (jak zrobiliśmy w innym procesie).

**Przed:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Po:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

I uruchom test!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Zauważ, że utworzyliśmy plik migawki dla procesu `convertToUpper` w `tests/main.converttoupper.nf.test.snap`. Jeśli ponownie uruchomimy test, powinniśmy zobaczyć, że nf-test ponownie przechodzi.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Wnioski

Wiesz, jak pisać testy dla procesu Nextflow i je uruchamiać.

### Co dalej?

Naucz się uruchamiać testy dla wszystkiego na raz!

## 3. Uruchomienie testów dla całego repozytorium

Uruchamianie nf-test dla każdego komponentu jest w porządku, ale pracochłonne i podatne na błędy. Czy nie możemy po prostu przetestować wszystkiego na raz?

Tak, możemy!

Uruchommy nf-test na całym repozytorium.

### 3.1. Uruchomienie nf-test na całym repozytorium

Możemy uruchomić nf-test na całym repozytorium, wykonując polecenie `nf-test test`.

```bash
nf-test test .
```

Zauważ, że używamy tylko `.`, aby uruchomić wszystko z naszego bieżącego katalogu. To będzie zawierać każdy test!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

Spójrz na to! Uruchomiliśmy 4 testy, 1 dla każdego procesu i 2 dla całego potoku za pomocą jednego polecenia. Wyobraź sobie, jak potężne jest to w dużej bazie kodu!

---

## Podsumowanie

W tym side queście nauczyłeś się wykorzystywać funkcje nf-test do tworzenia i uruchamiania testów dla pojedynczych procesów, a także testów end-to-end dla całego potoku.
Jesteś teraz świadomy dwóch głównych podejść do walidacji wyjścia, migawek i bezpośrednich asercji zawartości, oraz kiedy używać któregokolwiek z nich.
Wiesz również, jak uruchamiać testy pojedynczo lub dla całego projektu.

Zastosowanie tych technik we własnej pracy pozwoli Ci zapewnić, że:

- Twój kod działa zgodnie z oczekiwaniami
- Zmiany nie psują istniejącej funkcjonalności
- Inni programiści mogą wnosić wkład z pewnością
- Problemy mogą być szybko identyfikowane i naprawiane
- Zawartość wyjścia odpowiada oczekiwaniom

### Kluczowe wzorce

1. Testy na poziomie potoku:
   - Podstawowe testowanie sukcesu
   - Weryfikacja liczby procesów
   - Sprawdzanie istnienia plików wyjściowych
2. Testy na poziomie procesu
3. Dwa podejścia do walidacji wyjścia:
   - Używanie migawek do pełnej weryfikacji wyjścia
   - Używanie bezpośrednich asercji zawartości do sprawdzania konkretnej zawartości
4. Uruchamianie wszystkich testów w repozytorium za pomocą jednego polecenia

### Dodatkowe zasoby

Sprawdź [dokumentację nf-test](https://www.nf-test.com/), aby poznać bardziej zaawansowane funkcje testowania i najlepsze praktyki. Możesz chcieć:

- Dodać bardziej kompleksowe asercje do swoich testów
- Napisać testy dla przypadków brzegowych i warunków błędów
- Skonfigurować ciągłą integrację do automatycznego uruchamiania testów
- Dowiedzieć się więcej o innych typach testów, takich jak testy przepływów pracy i modułów
- Zbadać bardziej zaawansowane techniki walidacji zawartości

**Pamiętaj:** Testy są żywą dokumentacją tego, jak powinien zachowywać się Twój kod. Im więcej testów napiszesz i im bardziej szczegółowe będą Twoje asercje, tym bardziej możesz być pewny niezawodności swojego potoku.

---

## Co dalej?

Wróć do [menu Side Questów](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
