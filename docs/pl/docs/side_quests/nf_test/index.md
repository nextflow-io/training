# Testowanie z nf-test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Możliwość systematycznego sprawdzania, czy każda część workflow'u działa zgodnie z oczekiwaniami, jest kluczowa dla reprodukowalności i długoterminowego utrzymania kodu — a podczas samego procesu tworzenia może być ogromnym ułatwieniem.

Poświęćmy chwilę na omówienie, dlaczego testowanie jest tak ważne. Tworząc workflow, jedną z pierwszych rzeczy, które robisz, jest zebranie danych testowych, o których wiesz, że są poprawne i powinny dać określony wynik. Dodajesz pierwszy proces do pipeline'u i podłączasz go do wejść, żeby działał. Następnie, żeby sprawdzić, czy wszystko gra, uruchamiasz go na danych testowych. Zakładając, że to działa, przechodzisz do kolejnego procesu i znowu uruchamiasz dane testowe. Powtarzasz ten cykl, aż uzyskasz pipeline, z którego jesteś zadowolony.

Potem może dodajesz prosty parametr logiczny, np. `--skip_process`. Teraz musisz uruchomić pipeline dwukrotnie — raz z każdym parametrem — żeby upewnić się, że działa zgodnie z oczekiwaniami. Ale chwila, jak sprawdzić, czy `--skip_process` rzeczywiście pomija dany proces? Trzeba zagłębić się w wyjścia albo przejrzeć pliki logów! To uciążliwe i podatne na błędy.

W miarę rozwijania pipeline'u szybko stanie się on na tyle złożony, że ręczne testowanie każdej iteracji będzie powolne i zawodne. Co więcej, jeśli znajdziesz błąd, bardzo trudno będzie ustalić, w którym dokładnie miejscu pipeline'u on powstaje. Tu właśnie wkracza testowanie.

Testowanie pozwala systematycznie sprawdzać, czy każda część pipeline'u działa zgodnie z oczekiwaniami. Korzyści z dobrze napisanych testów są dla programisty ogromne:

- **Pewność**: Ponieważ testy obejmują cały pipeline, możesz być pewny, że zmiana jednej rzeczy nie wpłynie na pozostałe.
- **Zaufanie**: Gdy nad pipeline'em pracuje wielu programistów, każdy wie, że inni nie zepsuli pipeline'u ani żadnego z jego komponentów.
- **Przejrzystość**: Testy wskazują, gdzie pipeline zawodzi, i ułatwiają namierzenie problemu. Pełnią też rolę dokumentacji, pokazując, jak uruchomić dany proces lub workflow.
- **Szybkość**: Ponieważ testy są zautomatyzowane, można je uruchamiać bardzo szybko i wielokrotnie. Możesz iterować sprawnie, z mniejszą obawą o wprowadzanie nowych błędów.

Istnieje wiele różnych rodzajów testów, które możemy pisać:

1. **Testy na poziomie modułu**: Dla poszczególnych procesów
2. **Testy na poziomie workflow'u**: Dla pojedynczego workflow'u
3. **Testy na poziomie pipeline'u**: Dla całego pipeline'u
4. **Testy wydajnościowe**: Dla szybkości i efektywności pipeline'u
5. **Testy obciążeniowe**: Ocena działania pipeline'u w ekstremalnych warunkach w celu określenia jego limitów

Testowanie poszczególnych procesów jest analogiczne do testów jednostkowych w innych językach. Testowanie workflow'u lub całego pipeline'u odpowiada temu, co w innych językach nazywa się testami integracyjnymi — sprawdzamy w nich interakcje między komponentami.

[**nf-test**](https://www.nf-test.com/) to narzędzie umożliwiające pisanie testów na poziomie modułu, workflow'u i pipeline'u. Krótko mówiąc, pozwala systematycznie sprawdzać, czy każda indywidualna część pipeline'u działa zgodnie z oczekiwaniami — _w izolacji_.

### Cele szkolenia

W tym zadaniu dodatkowym nauczysz się używać nf-test do pisania testu na poziomie workflow'u dla pipeline'u, a także testów na poziomie modułu dla trzech wywoływanych przez niego procesów.

Po ukończeniu tego zadania będziesz potrafić efektywnie stosować następujące techniki:

- Inicjalizować nf-test w swoim projekcie
- Generować testy na poziomie modułu i workflow'u
- Dodawać typowe rodzaje asercji
- Rozumieć, kiedy używać snapshotów, a kiedy asercji treści
- Uruchamiać testy dla całego projektu

Te umiejętności pomogą Ci wdrożyć kompleksową strategię testowania w projektach pipeline'owych, czyniąc je bardziej niezawodnymi i łatwymi w utrzymaniu.

### Wymagania wstępne

Przed przystąpieniem do tego zadania dodatkowego powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow (procesy, kanały, operatory, praca z plikami, metadane).

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/nf-test
```

Możesz ustawić VSCode tak, żeby skupiał się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz tu główny plik workflow'u oraz plik CSV o nazwie `greetings.csv`, zawierający dane wejściowe do pipeline'u.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Szczegółowy opis plików znajdziesz w [rozgrzewce z Hello Nextflow](../hello_nextflow/00_orientation.md).

Workflow, który będziemy testować, jest podzbiorem workflow'u Hello zbudowanego w [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "Co robi workflow Hello Nextflow?"

    Jeśli nie przechodziłeś szkolenia [Hello Nextflow](../hello_nextflow/index.md), oto krótki przegląd tego, co robi ten prosty workflow.

    Workflow przyjmuje plik CSV zawierający pozdrowienia, przepuszcza je przez cztery kolejne kroki transformacji i zwraca pojedynczy plik tekstowy zawierający obrazek ASCII z zabawną postacią wypowiadającą te pozdrowienia.

    Cztery kroki są zaimplementowane jako procesy Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) przechowywane w osobnych plikach modułów.

    1. **`sayHello`:** Zapisuje każde pozdrowienie do własnego pliku wyjściowego (np. "Hello-output.txt")
    2. **`convertToUpper`:** Konwertuje każde pozdrowienie na wielkie litery (np. "HELLO")
    3. **`collectGreetings`:** Zbiera wszystkie pozdrowienia pisane wielkimi literami do jednego pliku zbiorczego
    4. **`cowpy`:** Generuje grafikę ASCII przy użyciu narzędzia `cowpy`

    Wyniki są publikowane w katalogu `results/`, a końcowe wyjście pipeline'u (przy uruchomieniu z domyślnymi parametrami) to zwykły plik tekstowy zawierający grafikę ASCII postaci wypowiadającej pozdrowienia zapisane wielkimi literami.

    W tym zadaniu dodatkowym używamy pośredniej formy workflow'u Hello, która zawiera tylko dwa pierwsze procesy. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

Podzbiór, z którym będziemy pracować, składa się z dwóch procesów: `sayHello` i `convertToUpper`.
Pełny kod workflow'u możesz zobaczyć poniżej.

??? example "Kod workflow'u"

    ```groovy title="main.nf"
    /*
    * Parametry pipeline'u
    */
    params.input_file = "greetings.csv"

    /*
    * Użyj echo, żeby wypisać 'Hello World!' na standardowe wyjście
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
    * Użyj narzędzia do zamiany tekstu, żeby przekonwertować pozdrowienie na wielkie litery
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

        // utwórz kanał dla danych wejściowych z pliku CSV
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // wyemituj pozdrowienie
        sayHello(greeting_ch)

        // przekonwertuj pozdrowienie na wielkie litery
        convertToUpper(sayHello.out)
    }
    ```

#### Uruchom workflow

Uruchommy workflow, żeby upewnić się, że działa zgodnie z oczekiwaniami.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

GRATULACJE! Właśnie uruchomiłeś test!

„Chwila, co? Po prostu uruchomiłem workflow i zadziałał! Jak to jest test?"

Dobre pytanie!

Przeanalizujmy, co właśnie się stało.

Uruchomiłeś workflow z domyślnymi parametrami, potwierdziłeś, że działa, i jesteś zadowolony z wyników. To jest właśnie istota testowania. Jeśli przechodziłeś kurs Hello Nextflow, zauważyłeś zapewne, że zawsze zaczynaliśmy każdą sekcję od uruchomienia workflow'u, którego używaliśmy jako punktu startowego — żeby potwierdzić, że wszystko jest poprawnie skonfigurowane.

Testowanie oprogramowania zasadniczo wykonuje ten proces za nas.

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest dodanie standaryzowanych testów do tego workflow'u przy użyciu nf-test, aby łatwo było weryfikować, czy każda część nadal działa zgodnie z oczekiwaniami w przypadku wprowadzenia dalszych zmian.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, żeby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem odpowiednio swój katalog roboczy
- [ ] Uruchomiłem workflow pomyślnie
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Inicjalizacja `nf-test`

Pakiet `nf-test` udostępnia polecenie inicjalizacji, które konfiguruje kilka rzeczy potrzebnych do rozpoczęcia tworzenia testów dla naszego projektu.

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

Polecenie tworzy też katalog `tests` zawierający szkielet pliku konfiguracyjnego.

### 1.1. Generowanie testu nf-test

`nf-test` zawiera zestaw narzędzi do budowania plików nf-test, oszczędzając nam większości pracy. Dostępne są one jako podpolecenie `generate`. Wygenerujmy test dla pipeline'u:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Spowoduje to utworzenie pliku `main.nf.test` w katalogu `tests`. To jest nasz plik testowy na poziomie pipeline'u. Po uruchomieniu `tree tests/` powinieneś zobaczyć coś takiego:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

Plik `main.nf.test` to nasz plik testowy na poziomie pipeline'u. Otwórzmy go i przyjrzyjmy się jego zawartości.

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

Blok `nextflow_pipeline` jest punktem wejścia dla wszystkich testów na poziomie pipeline'u. Zawiera:

- `name`: Nazwa testu.
- `script`: Ścieżka do skryptu pipeline'u.

Blok `test` to właściwy test. Zawiera:

- `when`: Warunki, w których test powinien być uruchomiony. Obejmuje parametry, które zostaną użyte do uruchomienia pipeline'u.
- `then`: Asercje, które powinny zostać sprawdzone. Zawiera oczekiwane wyniki działania pipeline'u.

Mówiąc po ludzku, logika testu brzmi następująco:
„**Gdy** do tego _pipeline'u_ zostaną przekazane te _parametry_, **wtedy** oczekujemy tych wyników."

To nie jest jeszcze funkcjonalny test — w następnej sekcji pokażemy, jak go takim uczynić.

### Uwaga o nazwach testów

W powyższym przykładzie użyliśmy domyślnej nazwy „Should run without failures", która jest odpowiednia dla podstawowego testu sprawdzającego jedynie, czy pipeline uruchamia się pomyślnie. Jednak w miarę dodawania bardziej szczegółowych przypadków testowych powinniśmy używać bardziej opisowych nazw, wskazujących, co faktycznie testujemy. Na przykład:

- „Should convert input to uppercase" — przy testowaniu konkretnej funkcjonalności
- „Should handle empty input gracefully" — przy testowaniu przypadków brzegowych
- „Should respect max memory parameter" — przy testowaniu ograniczeń zasobów
- „Should create expected output files" — przy testowaniu generowania plików

Dobre nazwy testów powinny:

1. Zaczynać się od „Should", żeby jasno określić oczekiwane zachowanie
2. Opisywać konkretną funkcjonalność lub scenariusz, który jest testowany
3. Być na tyle jasne, że gdy test nie przejdzie, od razu wiadomo, jaka funkcjonalność jest zepsuta

W miarę dodawania kolejnych asercji i szczegółowych przypadków testowych będziemy używać tych bardziej opisowych nazw, żeby jasno określić, co każdy test weryfikuje.

### 1.2. Uruchomienie testu

Uruchommy test, żeby zobaczyć, co się stanie.

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

Test nie przeszedł! Co się stało?

1. nf-test próbował uruchomić pipeline w obecnym stanie, używając ustawień z bloku `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test sprawdził status pipeline'u i porównał go z blokiem `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Zwróć uwagę, jak nf-test poinformował o niepowodzeniu pipeline'u i podał komunikat błędu z Nextflow:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Na czym polegał problem? Pamiętaj, że pipeline ma plik `greetings.csv` w katalogu projektu. Gdy nf-test uruchamia pipeline, szuka tego pliku, ale nie może go znaleźć. Plik jest tam, więc co się dzieje? Jeśli spojrzymy na ścieżkę, zobaczymy, że test odbywa się w ścieżce `./nf-test/tests/longHashString/`. Podobnie jak Nextflow, nf-test tworzy nowy katalog dla każdego testu, żeby wszystko było izolowane. Plik z danymi nie znajduje się tam, więc musimy poprawić ścieżkę do pliku w oryginalnym teście.

Wróćmy do pliku testowego i zmieńmy ścieżkę do pliku w bloku `when`.

Możesz się zastanawiać, jak wskazać korzeń pipeline'u w teście. Ponieważ jest to częsta sytuacja, nf-test udostępnia szereg zmiennych globalnych, które ułatwiają nam życie. Pełną listę znajdziesz [tutaj](https://www.nf-test.com/docs/testcases/global_variables/), ale na razie użyjemy zmiennej `projectDir`, która oznacza korzeń projektu pipeline'u.

=== "Po"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
    when {
        params {
            input_file = "${projectDir}/greetings.csv"
        }
    }
    ```

=== "Przed"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
    when {
        params {
            // define parameters here. Example:
            // outdir = "tests/results"
        }
    }
    ```

Uruchommy test ponownie, żeby sprawdzić, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Sukces! Pipeline uruchamia się pomyślnie i test przechodzi. Uruchamiaj go tyle razy, ile chcesz — zawsze otrzymasz ten sam wynik!

Domyślnie wyjście Nextflow jest ukryte, ale żeby przekonać się, że nf-test rzeczywiście uruchamia workflow, możesz użyć flagi `--verbose`:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
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

### 1.3. Dodawanie asercji

Prostym sprawdzeniem jest upewnienie się, że nasz pipeline uruchamia wszystkie oczekiwane procesy i żadnego nie pomija po cichu. Pamiętaj, że nasz pipeline uruchamia 6 procesów — jeden `sayHello` i jeden `convertToUpper` dla każdego z 3 pozdrowień.

Dodajmy asercję do naszego testu, żeby sprawdzić, czy pipeline uruchamia oczekiwaną liczbę procesów. Zaktualizujemy też nazwę testu, żeby lepiej odzwierciedlała to, co testujemy.

=== "Po"

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

=== "Przed"

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

Nazwa testu teraz lepiej odzwierciedla to, co faktycznie weryfikujemy — nie tylko to, że pipeline uruchamia się bez błędów, ale że uruchamia oczekiwaną liczbę procesów.

Uruchommy test ponownie, żeby sprawdzić, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Sukces! Pipeline uruchamia się pomyślnie i test przechodzi. Zaczęliśmy teraz testować szczegóły pipeline'u, a nie tylko jego ogólny status.

### 1.4. Testowanie wyjścia

Dodajmy asercję do naszego testu, żeby sprawdzić, czy plik wyjściowy został utworzony. Dodamy go jako osobny test z informacyjną nazwą, żeby wyniki były łatwiejsze do interpretacji.

=== "Po"

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

=== "Przed"

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

Uruchom test ponownie, żeby sprawdzić, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Sukces! Testy przechodzą, ponieważ pipeline zakończył się pomyślnie, uruchomiono właściwą liczbę procesów i pliki wyjściowe zostały utworzone. Powinieneś też teraz zobaczyć, jak przydatne jest nadawanie testom informacyjnych nazw.

To tylko wierzchołek góry lodowej — możemy pisać kolejne asercje sprawdzające szczegóły pipeline'u, ale na razie przejdźmy do testowania jego wewnętrznych komponentów.

### Podsumowanie

Wiesz już, jak pisać testy nf-test dla pipeline'u.

### Co dalej?

Naucz się testować proces Nextflow.

---

## 2. Testowanie procesu Nextflow

Nie musimy pisać testów dla każdej części pipeline'u, ale im więcej testów mamy, tym bardziej kompleksowe jest nasze pokrycie i tym większą pewność możemy mieć, że pipeline działa zgodnie z oczekiwaniami. W tej sekcji przetestujemy oba procesy pipeline'u jako indywidualne jednostki.

### 2.1. Testowanie procesu `sayHello`

Zacznijmy od procesu `sayHello`.

Użyjmy ponownie polecenia `nf-test generate`, żeby wygenerować testy dla procesu.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Skupmy się teraz na procesie `sayhello` w pliku `main.sayhello.nf.test`.

Otwórzmy plik i przyjrzyjmy się jego zawartości.

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

Podobnie jak wcześniej, zaczynamy od szczegółów testu, po których następują bloki `when` i `then`. Mamy jednak dodatkowy blok `process`, który pozwala nam zdefiniować dane wejściowe do procesu.

Uruchommy test, żeby sprawdzić, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
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

Test nie przechodzi, ponieważ proces `sayHello` deklaruje 1 wejście, ale został wywołany z 0 argumentami. Naprawmy to, dodając wejście do procesu. Pamiętaj z [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (i sekcji rozgrzewki powyżej), że nasz proces `sayHello` przyjmuje pojedyncze wejście wartości, które musimy dostarczyć. Powinniśmy też poprawić nazwę testu, żeby lepiej odzwierciedlała to, co testujemy.

=== "Po"

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

=== "Przed"

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

Uruchommy test ponownie, żeby sprawdzić, czy działa.

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

Sukces! Test przechodzi, ponieważ proces `sayHello` uruchomił się pomyślnie i wyjście zostało utworzone.

### 2.2. Sprawdzenie snapshotu utworzonego przez test

Jeśli spojrzymy na plik `tests/main.sayhello.nf.test`, zobaczymy, że używa metody `snapshot()` w bloku asercji:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Mówi to nf-test, żeby utworzył snapshot wyjścia procesu `sayHello`. Przyjrzyjmy się zawartości pliku snapshotu.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

Nie będziemy go tu drukować, ale powinieneś zobaczyć plik JSON zawierający szczegóły procesu i jego wyjść. W szczególności możemy zobaczyć linię wyglądającą tak:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Reprezentuje to wyjścia utworzone przez proces `sayHello`, które testujemy wprost. Jeśli ponownie uruchomimy test, program sprawdzi, czy nowe wyjście pasuje do wyjścia pierwotnie zarejestrowanego. To szybki, prosty sposób testowania, czy wyjścia procesów się nie zmieniają — dlatego nf-test udostępnia go jako domyślny.

!!!warning "Ostrzeżenie"

    Oznacza to, że musimy być pewni, że wyjście zarejestrowane w pierwotnym uruchomieniu jest poprawne!

Jeśli w toku dalszego rozwoju coś w kodzie zmieni się tak, że wyjście będzie inne, test nie przejdzie i będziemy musieli ustalić, czy zmiana jest oczekiwana, czy nie.

- Jeśli okaże się, że coś w kodzie się zepsuło, będziemy musieli to naprawić, oczekując, że poprawiony kod przejdzie test.
- Jeśli jest to oczekiwana zmiana (np. narzędzie zostało ulepszone i wyniki są lepsze), będziemy musieli zaktualizować snapshot, żeby zaakceptować nowe wyjście jako referencję do porównania. nf-test ma do tego celu parametr `--update-snapshot`.

Możemy uruchomić test ponownie i zobaczyć, że powinien przejść:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Sukces! Test przechodzi, ponieważ proces `sayHello` uruchomił się pomyślnie i wyjście pasuje do snapshotu.

### 2.3. Alternatywa dla snapshotów: bezpośrednie asercje treści

Snapshoty świetnie nadają się do wykrywania wszelkich zmian w wyjściu, ale czasem chcemy zweryfikować konkretną treść bez tak rygorystycznego wymagania, żeby cały plik pasował. Na przykład:

- Gdy części wyjścia mogą się zmieniać (znaczniki czasu, losowe identyfikatory itp.), ale pewna kluczowa treść musi być obecna
- Gdy chcemy sprawdzić konkretne wzorce lub wartości w wyjściu
- Gdy chcemy uczynić test bardziej precyzyjnym co do tego, co stanowi sukces

Oto jak możemy zmodyfikować nasz test, żeby sprawdzał konkretną treść:

=== "Po"

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

=== "Przed"

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

Zwróć uwagę, że nf-test widzi wyjścia procesu jako listę list, więc `process.out[0][0]` pobiera pierwszą część pierwszego elementu kanału (czyli pierwszej „emisji") z tego procesu.

To podejście:

- Jasno określa, czego oczekujemy w wyjściu
- Jest bardziej odporne na nieistotne zmiany w wyjściu
- Dostarcza lepszych komunikatów błędów, gdy testy nie przechodzą
- Umożliwia bardziej złożone walidacje (wzorce regex, porównania numeryczne itp.)

Uruchommy test, żeby sprawdzić, czy działa.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Testowanie procesu `convertToUpper`

Otwórzmy plik `tests/main.converttoupper.nf.test` i przyjrzyjmy się jego zawartości:

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

To podobny test do tego dla procesu `sayHello`, ale testuje proces `convertToUpper`. Wiemy, że ten test nie przejdzie, ponieważ — podobnie jak w przypadku `sayHello` — proces `convertToUpper` przyjmuje pojedyncze wejście ścieżki, którego nie podaliśmy.

Musimy teraz dostarczyć pojedynczy plik wejściowy do procesu `convertToUpper`, zawierający tekst, który chcemy przekonwertować na wielkie litery. Możemy to zrobić na wiele sposobów:

- Możemy utworzyć dedykowany plik do testów
- Możemy ponownie użyć istniejącego pliku `data/greetings.csv`
- Możemy go utworzyć w locie w ramach testu

Na razie ponownie użyjmy istniejącego pliku `data/greetings.csv`, korzystając z przykładu użytego w teście na poziomie pipeline'u. Jak poprzednio, możemy nazwać test tak, żeby lepiej odzwierciedlał to, co testujemy, ale tym razem zostawmy „snapshot" treści zamiast sprawdzać konkretne ciągi znaków (jak zrobiliśmy to w przypadku drugiego procesu).

=== "Po"

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

=== "Przed"

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

I uruchommy test!

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

Zwróć uwagę, że utworzyliśmy plik snapshotu dla procesu `convertToUpper` w `tests/main.converttoupper.nf.test.snap`. Jeśli uruchomimy test ponownie, powinniśmy zobaczyć, że nf-test znowu przechodzi.

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

### Podsumowanie

Wiesz już, jak pisać testy dla procesu Nextflow i je uruchamiać.

### Co dalej?

Naucz się uruchamiać wszystkie testy naraz!

## 3. Uruchamianie testów dla całego repozytorium

Uruchamianie nf-test dla każdego komponentu osobno jest możliwe, ale żmudne i podatne na błędy. Czy nie możemy przetestować wszystkiego naraz?

Możemy!

Uruchommy nf-test na całym repozytorium.

### 3.1. Uruchamianie nf-test na całym repozytorium

Możemy uruchomić nf-test na całym repozytorium, używając polecenia `nf-test test`.

```bash
nf-test test .
```

Zwróć uwagę, że używamy po prostu `.`, żeby uruchomić wszystko z naszego bieżącego katalogu. Obejmie to każdy test!

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

Spójrz na to! Uruchomiliśmy 4 testy — 1 dla każdego procesu i 2 dla całego pipeline'u — jednym poleceniem. Wyobraź sobie, jak potężne jest to w przypadku dużej bazy kodu!

---

## Podsumowanie

W tym zadaniu dodatkowym nauczyłeś się wykorzystywać funkcje nf-test do tworzenia i uruchamiania testów dla poszczególnych procesów, a także testów end-to-end dla całego pipeline'u.
Znasz teraz dwa główne podejścia do walidacji wyjść — snapshoty i bezpośrednie asercje treści — oraz wiesz, kiedy używać każdego z nich.
Wiesz też, jak uruchamiać testy jeden po drugim lub dla całego projektu.

Stosowanie tych technik we własnej pracy pozwoli Ci zapewnić, że:

- Twój kod działa zgodnie z oczekiwaniami
- Zmiany nie psują istniejącej funkcjonalności
- Inni programiści mogą wnosić wkład z pewnością siebie
- Problemy można szybko zidentyfikować i naprawić
- Treść wyjść spełnia oczekiwania

### Kluczowe wzorce

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Testy na poziomie pipeline'u:
   - Podstawowe testowanie sukcesu
   - Weryfikacja liczby procesów
   - Sprawdzanie istnienia plików wyjściowych
2. Testy na poziomie procesu
3. Dwa podejścia do walidacji wyjść:
   - Używanie snapshotów do pełnej weryfikacji wyjścia
   - Używanie bezpośrednich asercji treści do sprawdzania konkretnej zawartości
4. Uruchamianie wszystkich testów w repozytorium jednym poleceniem

### Dodatkowe zasoby

Zajrzyj do [dokumentacji nf-test](https://www.nf-test.com/), żeby poznać bardziej zaawansowane funkcje testowania i najlepsze praktyki. Możesz chcieć:

- Dodać bardziej kompleksowe asercje do swoich testów
- Pisać testy dla przypadków brzegowych i warunków błędów
- Skonfigurować ciągłą integrację, żeby testy uruchamiały się automatycznie
- Poznać inne rodzaje testów, takie jak testy workflow'u i modułów
- Zgłębić bardziej zaawansowane techniki walidacji treści

**Pamiętaj:** Testy to żywa dokumentacja tego, jak Twój kod powinien się zachowywać. Im więcej testów piszesz i im bardziej szczegółowe są Twoje asercje, tym większą pewność możesz mieć co do niezawodności swojego pipeline'u.

---

## Co dalej?

Wróć do [menu zadań dodatkowych](../) lub kliknij przycisk w prawym dolnym rogu strony, żeby przejść do następnego tematu na liście.
