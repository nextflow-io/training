# Workflow'y Workflow'ów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Podczas tworzenia pipeline'u często zdarza się, że piszemy podobne sekwencje procesów dla różnych typów danych lub etapów analizy. Można wtedy skończyć na kopiowaniu i wklejaniu tych sekwencji, co prowadzi do zduplikowanego kodu, trudnego w utrzymaniu — albo stworzyć jeden ogromny workflow, który jest trudny do zrozumienia i modyfikacji.

Jedną z najpotężniejszych funkcji Nextflow'a jest możliwość komponowania złożonych pipeline'ów z mniejszych, wielokrotnego użytku modułów workflow'ów. Takie podejście modularne sprawia, że pipeline'y są łatwiejsze do tworzenia, testowania i utrzymania.

### Cele szkolenia

W tym side queście zbadamy, jak tworzyć moduły workflow'ów, które można testować i używać osobno, składać je w większy pipeline oraz zarządzać przepływem danych między modułami.

Po ukończeniu tego side questu będziesz potrafić:

- Rozkładać złożone pipeline'y na logiczne, wielokrotnego użytku jednostki
- Testować każdy moduł workflow'u niezależnie
- Łączyć workflow'y w celu tworzenia nowych pipeline'ów
- Współdzielić wspólne moduły workflow'ów między różnymi pipeline'ami
- Pisać kod bardziej czytelny i łatwiejszy w utrzymaniu

Te umiejętności pomogą Ci budować złożone pipeline'y przy zachowaniu przejrzystej, łatwej w utrzymaniu struktury kodu.

### Wymagania wstępne

Przed przystąpieniem do tego side questu powinieneś/powinnaś:

- Ukończyć samouczek [Hello Nextflow](../../hello_nextflow/index.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow'a (procesy, kanały, operatory, moduły).

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś/zrobiłaś, otwórz środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/workflows_of_workflows
```

Możesz ustawić VSCode tak, aby skupiał się na tym katalogu:

```bash
code .
```

Edytor otworzy się z widokiem na katalog projektu.

#### Przejrzyj materiały

Znajdziesz katalog `modules` z definicjami procesów, katalog `workflows` z dwoma gotowymi skryptami workflow'ów oraz plik `main.nf`, który będziesz stopniowo aktualizować:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

Katalog `modules/` zawiera definicje poszczególnych procesów, a katalog `workflows/` — dwa gotowe skrypty workflow'ów, z którymi będziesz pracować w tym side queście.

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest złożenie tych modułów w dwa osobne workflow'y, które następnie skomponujemy w główny workflow:

- `GREETING_WORKFLOW` — waliduje nazwy, tworzy powitania i dodaje znaczniki czasu
- `TRANSFORM_WORKFLOW` — konwertuje tekst na wielkie litery i odwraca go

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem/ustawiłam odpowiedni katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Dodanie greeting workflow do pipeline'u

Greeting workflow waliduje nazwy i generuje powitania ze znacznikami czasu.

### 1.1. Przejrzyj i uruchom greeting workflow

Otwórz `workflows/greeting.nf` i przyjrzyj się kodowi:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Łańcuch procesów: walidacja -> tworzenie powitania -> dodanie znacznika czasu
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Jest to kompletny, samodzielny workflow o strukturze podobnej do tych, które widziałeś/widziałaś w samouczku 'Hello Nextflow'.
Nazwy wejściowe są w nim zakodowane na stałe, trzy procesy są połączone w łańcuch, a dwa wyjścia są publikowane.

Uruchom go, aby sprawdzić, czy wszystko działa:

```bash
nextflow run workflows/greeting.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Aby uczynić go kompozytowalnym z innymi workflow'ami, konieczne jest wprowadzenie kilku zmian.

### 1.2. Uczyń workflow kompozytowalnym

Aby workflow był kompozytowalny, należy wprowadzić cztery zmiany:
workflow otrzymuje nazwę, wejścia przenosi się do bloku `take:`, wyjścia przenosi się do bloku `emit:`,
a samodzielne bloki `publish:`/`output {}` są usuwane (należą do entry workflow).

Omówmy te zmiany po kolei.

#### 1.2.1. Nadaj workflow'owi nazwę

Nadaj workflow'owi nazwę, aby można go było importować z nadrzędnego workflow'u.

=== "Po"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Przed"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Dzięki nazwie workflow można go importować do innych skryptów.

#### 1.2.2. Zadeklaruj wejścia przy użyciu `take:`

Zastąp zakodowaną na stałe deklarację kanału blokiem `take:`, który deklaruje oczekiwane wejścia workflow'u.
Blok `take:` umieszcza się przed `main:`, a wiersz `names_ch = channel.of(...)` jest usuwany.

=== "Po"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // Kanał wejściowy z nazwami

        main:
        // Łańcuch procesów: walidacja -> tworzenie powitania -> dodanie znacznika czasu
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Przed"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Łańcuch procesów: walidacja -> tworzenie powitania -> dodanie znacznika czasu
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

Blok `take:` deklaruje kanał wyłącznie przez nazwę — szczegóły dotyczące jego zawartości zostaną zdefiniowane przez nadrzędny workflow.

#### 1.2.3. Zadeklaruj wyjścia przy użyciu `emit:`

Zastąp sekcję `publish:` i usuń blok `output {}`, zastępując je blokiem `emit:`, który nadaje wyjściom nazwy.

=== "Po"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Oryginalne powitania
        timestamped = timestamped_ch // Powitania ze znacznikami czasu
    }
    ```

=== "Przed"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

Blok `emit:` udostępnia nazwane wyjścia, do których nadrzędne workflow'y mogą uzyskiwać dostęp za pomocą `GREETING_WORKFLOW.out.greetings` i `GREETING_WORKFLOW.out.timestamped`.

#### 1.2.4. Zweryfikuj wynik i przetestuj

Po wprowadzeniu wszystkich trzech zmian kompletny plik powinien wyglądać następująco:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // Kanał wejściowy z nazwami

    main:
    // Łańcuch procesów: walidacja -> tworzenie powitania -> dodanie znacznika czasu
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Oryginalne powitania
    timestamped = timestamped_ch // Powitania ze znacznikami czasu
}
```

Spróbuj teraz uruchomić go bezpośrednio:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Pojawia się tu kluczowy koncept: **entry workflow**.
Nextflow używa nienazwanego bloku `workflow {}` jako punktu wejścia podczas bezpośredniego uruchamiania skryptu.
`GREETING_WORKFLOW` ma nazwę, więc Nextflow nie wie, jak go samodzielnie uruchomić.

Jest to zamierzone — kompozytowalne workflow'y są zaprojektowane do wywoływania z entry workflow, a nie do bezpośredniego uruchamiania.
Rozwiązaniem jest stworzenie entry workflow w `main.nf`, który importuje i wywołuje `GREETING_WORKFLOW`.

### 1.3. Zaktualizuj i przetestuj główny workflow

Zaktualizujmy teraz główny workflow, aby wywoływał greeting workflow.

#### 1.3.1. Dołącz greeting workflow i wywołaj go

Dodaj instrukcję `include`, zaktualizuj ciało workflow'u tak, aby wywoływało `GREETING_WORKFLOW`, i zastąp placeholder `channel.empty()` w sekcji `publish:`:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Uruchom greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

Entry workflow pozostaje nienazwany, dzięki czemu Nextflow użyje go jako punktu wejścia pipeline'u.

#### 1.3.2. Zaktualizuj blok output

Dodaj dyrektywę `path`, aby kierować publikowane powitania do podkatalogu `greetings/`:

=== "Po"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. Uruchom workflow

Uruchom workflow, aby sprawdzić, czy działa poprawnie:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    ```

??? abstract "Zawartość katalogu"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "Zawartość pliku"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Pliki z powitaniami są publikowane do `results/greetings/`.
Główny workflow wywołuje `GREETING_WORKFLOW` i przekazuje jego wyjście bezpośrednio do sekcji `publish:`.

### Podsumowanie

W tej sekcji poznałeś/poznałaś kilka ważnych konceptów:

- **Nazwane workflow'y**: Tworzenie nazwanego workflow'u (`GREETING_WORKFLOW`), który można importować i ponownie używać
- **Interfejsy workflow'ów**: Definiowanie wyraźnych wejść za pomocą `take:` i wyjść za pomocą `emit:` w celu stworzenia kompozytowalnego workflow'u
- **Entry points**: Zrozumienie, że Nextflow potrzebuje nienazwanego entry workflow, aby uruchomić skrypt
- **Kompozycja workflow'ów**: Importowanie i używanie nazwanego workflow'u wewnątrz innego workflow'u
- **Przestrzenie nazw workflow'ów**: Dostęp do wyjść workflow'u za pomocą notacji `.out` (`GREETING_WORKFLOW.out.greetings`)

Masz teraz działający greeting workflow, który:

- Przyjmuje kanał nazw jako wejście
- Waliduje każdą nazwę
- Tworzy powitanie dla każdej prawidłowej nazwy
- Dodaje znaczniki czasu do powitań
- Udostępnia zarówno oryginalne, jak i opatrzone znacznikami czasu powitania jako wyjścia

Takie podejście modularne pozwala testować greeting workflow niezależnie lub używać go jako komponentu w większych pipeline'ach.

---

## 2. Dodanie transform workflow do pipeline'u

Transform workflow stosuje transformacje tekstu do powitań ze znacznikami czasu.

### 2.1. Przejrzyj i uruchom workflow

Otwórz `workflows/transform.nf` i przyjrzyj się kodowi:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Zastosuj transformacje sekwencyjnie
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Ten samodzielny workflow odczytuje pliki z powitaniami ze znacznikami czasu z katalogu `results/` wyprodukowanego przez `greeting.nf`, konwertuje je na wielkie litery, a następnie odwraca tekst.

Uruchom go, aby sprawdzić, czy działa poprawnie z wynikami z sekcji 1.1:

```bash
nextflow run workflows/transform.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

Aby uczynić go kompozytowalnym z `GREETING_WORKFLOW`, należy zastosować te same trzy zmiany co w sekcji 1.2.

### 2.2. Uczyń go kompozytowalnym

Zastosuj te same trzy zmiany co w sekcji 1.2: nadaj workflow'owi nazwę, zastąp zakodowane na stałe wejście blokiem `take:` i zastąp `publish:`/`output {}` blokiem `emit:`.

Gotowy plik powinien wyglądać następująco:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Kanał wejściowy z wiadomościami

    main:
    // Zastosuj transformacje sekwencyjnie
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Powitania wielkimi literami
    reversed = reversed_ch // Odwrócone powitania wielkimi literami
}
```

Transform workflow jest teraz kompozytowalny i gotowy do zaimportowania do głównego workflow'u.

### 2.3. Zaktualizuj i przetestuj główny workflow

Zaktualizujmy teraz główny workflow, aby wywoływał transform workflow.

#### 2.3.1. Dołącz transform workflow i wywołaj go

Dodaj instrukcję include, wywołanie `TRANSFORM_WORKFLOW` połączone z powitaniami ze znacznikami czasu oraz dwa nowe wpisy w sekcji `publish:`:

=== "Po"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Uruchom greeting workflow
        GREETING_WORKFLOW(names)

        // Uruchom transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Uruchom greeting workflow
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Spowoduje to uruchomienie transform workflow na powitaniach ze znacznikami czasu.

#### 2.3.2. Zaktualizuj blok output

Dodaj wpisy `upper` i `reversed` do bloku `output {}`, każdy z dyrektywą `path` wskazującą odpowiedni podkatalog:

=== "Po"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Przed"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Spowoduje to opublikowanie końcowych wyników do odpowiednich katalogów.

#### 2.3.3. Uruchom kompletny pipeline

Uruchom pipeline, aby sprawdzić, czy wszystko działa:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Zawartość katalogu"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "Zawartość pliku"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

Pipeline działa od początku do końca: powitanie zostało zamienione na wielkie litery i odwrócone.

### Podsumowanie

Powinieneś/powinnaś mieć teraz kompletny pipeline, który:

- Przetwarza nazwy przez greeting workflow
- Przekazuje powitania ze znacznikami czasu do transform workflow
- Produkuje zarówno wersje wielkimi literami, jak i odwrócone wersje powitań

---

## Podsumowanie

W tym side queście zbadaliśmy potężny koncept kompozycji workflow'ów w Nextflow, który pozwala budować złożone pipeline'y z mniejszych, wielokrotnego użytku komponentów.

Takie podejście modularne oferuje kilka zalet w porównaniu z monolitycznymi pipeline'ami:

- Każdy workflow można rozwijać, testować i debugować niezależnie
- Workflow'y można ponownie używać w różnych pipeline'ach
- Ogólna struktura pipeline'u staje się bardziej czytelna i łatwiejsza w utrzymaniu
- Zmiany w jednym workflow'ie niekoniecznie wpływają na inne, jeśli interfejsy pozostają spójne
- Entry points można konfigurować tak, aby uruchamiały różne części pipeline'u w zależności od potrzeb

_Warto jednak pamiętać, że choć wywoływanie workflow'ów jest nieco podobne do wywoływania procesów, nie jest tym samym. Nie można na przykład uruchomić workflow'u N razy, wywołując go z kanałem o rozmiarze N — należy przekazać kanał o rozmiarze N do workflow'u i iterować wewnętrznie._

Stosowanie tych technik we własnej pracy pozwoli Ci budować bardziej zaawansowane pipeline'y Nextflow, które mogą obsługiwać złożone zadania przetwarzania danych, pozostając przy tym łatwymi w utrzymaniu i skalowalnymi.

### Kluczowe wzorce

1.  **Struktura workflow'u**: Zdefiniowaliśmy wyraźne wejścia i wyjścia dla każdego workflow'u przy użyciu składni `take:` i `emit:`, tworząc dobrze zdefiniowane interfejsy między komponentami, a logikę workflow'u umieściliśmy wewnątrz bloku `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Kanały wejściowe są deklarowane tutaj
            input_ch

        main:
            // Logika workflow'u jest tutaj
            // Tu wywołuje się procesy i manipuluje kanałami
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Kanały wyjściowe są deklarowane tutaj
            output_ch = result_ch
    }
    ```

2.  **Importowanie workflow'ów:** Zbudowaliśmy dwa niezależne moduły workflow'ów i zaimportowaliśmy je do głównego pipeline'u za pomocą instrukcji include.

    - Importowanie pojedynczego workflow'u

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Importowanie wielu workflow'ów

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Importowanie z aliasem, aby uniknąć konfliktów nazw

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry points**: Nextflow wymaga nienazwanego entry workflow, aby wiedzieć, od czego zacząć wykonanie. Ten entry workflow wywołuje Twoje nazwane workflow'y.

    - Nienazwany workflow (entry point)

    ```groovy
    workflow {
        // To jest punkt wejścia podczas uruchamiania skryptu
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Nazwany workflow (wywoływany z entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Musi być wywoływany z entry workflow
    }
    ```

4.  **Zarządzanie przepływem danych:** Nauczyliśmy się, jak uzyskiwać dostęp do wyjść workflow'u za pomocą notacji przestrzeni nazw (`WORKFLOW_NAME.out.channel_name`) i przekazywać je do innych workflow'ów.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Dodatkowe zasoby

- [Dokumentacja Nextflow Workflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Dokumentacja operatorów kanałów](https://www.nextflow.io/docs/latest/operator.html)
- [Dokumentacja strategii obsługi błędów](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Co dalej?

Wróć do [menu Side Quests](../index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
