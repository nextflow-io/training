# Workflow'y workflow'ów

Podczas tworzenia pipeline'u często zdarza się, że tworzysz podobne sekwencje procesów dla różnych typów danych lub etapów analizy. Możesz w końcu kopiować i wklejać te sekwencje procesów, co prowadzi do zduplikowanego kodu, który jest trudny w utrzymaniu; albo możesz stworzyć jeden masywny workflow, który jest trudny do zrozumienia i modyfikacji.

Jedną z najpotężniejszych funkcji Nextflow'a jest możliwość komponowania złożonych pipeline'ów z mniejszych, wielokrotnego użytku modułów workflow'ów. To modularne podejście sprawia, że pipeline'y są łatwiejsze do rozwijania, testowania i utrzymania.

### Cele szkoleniowe

W tej misji pobocznej zbadamy, jak tworzyć moduły workflow'ów, które można testować i używać osobno, komponować te moduły w większy pipeline oraz zarządzać przepływem danych między modułami.

Pod koniec tej misji pobocznej będziesz w stanie:

- Rozbijać złożone pipeline'y na logiczne, wielokrotnego użytku jednostki
- Testować każdy moduł workflow'a niezależnie
- Mieszać i dopasowywać workflow'y, aby tworzyć nowe pipeline'y
- Dzielić się wspólnymi modułami workflow'ów między różnymi pipeline'ami
- Sprawić, by Twój kod był bardziej łatwy w utrzymaniu i zrozumieniu

Te umiejętności pomogą Ci budować złożone pipeline'y, zachowując czystą, łatwą w utrzymaniu strukturę kodu.

### Wymagania wstępne

Przed podjęciem tej misji pobocznej powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow'a (procesy, kanały, operatory, moduły)

---

## 0. Rozpocznij pracę

#### Otwórz przestrzeń kodową szkolenia

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzysz środowisko szkoleniowe zgodnie z opisem w [Konfiguracji środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki do tego samouczka.

```bash
cd side-quests/workflows_of_workflows
```

Możesz ustawić VSCode, aby skupił się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz katalog `modules` zawierający kilka definicji procesów, które rozwijają to, czego nauczyłeś się w 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### Przejrzyj zadanie

Twoim wyzwaniem jest złożenie tych modułów w dwa oddzielne workflow'y, które następnie skomponujemy w główny workflow:

- `GREETING_WORKFLOW`, który waliduje nazwy, tworzy powitania i dodaje znaczniki czasu
- `TRANSFORM_WORKFLOW`, który konwertuje tekst na wielkie litery i odwraca go

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moja przestrzeń kodowa działa
- [ ] Ustawiłem odpowiednio mój katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Utwórz workflow powitań

Zacznijmy od stworzenia workflow'a, który waliduje nazwy i generuje powitania ze znacznikami czasu.

### 1.1. Utwórz strukturę workflow'a

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Dodaj kod pierwszego (pod)workflow'a

Dodaj ten kod do `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

To jest kompletny workflow o strukturze podobnej do tych, które widziałeś w samouczku 'Hello Nextflow', który możemy testować niezależnie. Spróbujmy tego teraz:

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

Działa zgodnie z oczekiwaniami, ale aby uczynić go komponowalnym, musimy zmienić kilka rzeczy.

### 1.3. Uczyń workflow komponowalnym

Komponowalne workflow'y mają kilka różnic w porównaniu z tymi, które widziałeś w samouczku 'Hello Nextflow':

- Blok `workflow` musi być nazwany
- Wejścia są deklarowane przy użyciu słowa kluczowego `take:`
- Zawartość workflow'a jest umieszczana wewnątrz bloku `main:`
- Wyjścia są deklarowane przy użyciu słowa kluczowego `emit:`

Zaktualizujmy workflow powitań, aby pasował do tej struktury. Zmień kod na następujący:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Input channel with names

    main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
}
```

Widzisz, że workflow jest teraz nazwany i ma bloki `take:` oraz `emit:`, a to są połączenia, których użyjemy do skomponowania workflow'a wyższego poziomu.
Zawartość workflow'a jest również umieszczona wewnątrz bloku `main:`. Zauważ również, że usunęliśmy deklarację kanału wejściowego `names_ch`, ponieważ jest on teraz przekazywany jako argument do workflow'a.

Przetestujmy workflow ponownie, aby sprawdzić, czy działa zgodnie z oczekiwaniami:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

To informuje Cię o kolejnej nowej koncepcji, 'workflow'ie wejściowym'. Workflow wejściowy to workflow, który jest wywoływany, gdy uruchamiasz skrypt Nextflow'a. Domyślnie Nextflow użyje nienazwanego workflow'a jako workflow'a wejściowego, gdy jest obecny, i to właśnie robiłeś do tej pory, z blokami `workflow` zaczynającymi się tak:

```groovy title="hello.nf" linenums="1"
workflow {
```

Ale nasz workflow powitań nie ma nienazwanego workflow'a, zamiast tego mamy nazwany workflow:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Dlatego Nextflow zgłosił błąd i nie zrobił tego, czego chcieliśmy.

Nie dodaliśmy składni `take:`/`emit:`, abyśmy mogli wywołać workflow bezpośrednio - zrobiliśmy to, abyśmy mogli go skomponować z innymi workflow'ami. Rozwiązaniem jest utworzenie głównego skryptu z nienazwanym workflow'em wejściowym, który importuje i wywołuje nasz nazwany workflow.

### 1.4. Utwórz i przetestuj główny workflow

Teraz utworzymy główny workflow, który importuje i używa workflow'a `greeting`.

Utwórz `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Zauważ, że nasz wpis workflow'a w tym pliku jest nienazwany, a to dlatego, że będziemy go używać jako workflow'a wejściowego.

Uruchom to i zobacz wyjście:

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
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Działa! Opakaliśmy nazwany workflow powitań w główny workflow z nienazwanym blokiem wejściowym `workflow`. Główny workflow używa workflow'a `GREETING_WORKFLOW` prawie (nie do końca) jak procesu i przekazuje kanał `names` jako argument.

### Podsumowanie

W tej sekcji nauczyłeś się kilku ważnych koncepcji:

- **Nazwane workflow'y**: Tworzenie nazwanego workflow'a (`GREETING_WORKFLOW`), który można importować i używać ponownie
- **Interfejsy workflow'ów**: Definiowanie jasnych wejść za pomocą `take:` i wyjść za pomocą `emit:`, aby stworzyć komponowalny workflow
- **Punkty wejścia**: Zrozumienie, że Nextflow potrzebuje nienazwanego workflow'a wejściowego, aby uruchomić skrypt
- **Komponowanie workflow'ów**: Importowanie i używanie nazwanego workflow'a wewnątrz innego workflow'a
- **Przestrzenie nazw workflow'ów**: Dostęp do wyjść workflow'a przy użyciu przestrzeni nazw `.out` (`GREETING_WORKFLOW.out.greetings`)

Masz teraz działający workflow powitań, który:

- Przyjmuje kanał nazw jako wejście
- Waliduje każdą nazwę
- Tworzy powitanie dla każdej prawidłowej nazwy
- Dodaje znaczniki czasu do powitań
- Udostępnia zarówno oryginalne, jak i powitania ze znacznikami czasu jako wyjścia

To modularne podejście pozwala testować workflow powitań niezależnie lub używać go jako komponentu w większych pipeline'ach.

---

## 2. Dodaj workflow transformacji

Teraz utwórzmy workflow, który stosuje transformacje tekstowe do powitań.

### 2.1. Utwórz plik workflow'a

```bash
touch workflows/transform.nf
```

### 2.2. Dodaj kod workflow'a

Dodaj ten kod do `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Input channel with messages

    main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase greetings
        reversed = reversed_ch  // Reversed uppercase greetings
}
```

Nie będziemy powtarzać wyjaśnienia składni komponowalnej tutaj, ale zauważ, że nazwany workflow jest ponownie zadeklarowany z blokami `take:` i `emit:`, a zawartość workflow'a jest umieszczona wewnątrz bloku `main:`.

### 2.3. Zaktualizuj główny workflow

Zaktualizuj `main.nf`, aby używał obu workflow'ów:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Uruchom kompletny pipeline:

```bash
nextflow run main.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Jeśli spojrzysz na jeden z tych odwróconych plików, zobaczysz, że jest to odwrócona wersja powitania wielkimi literami:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Podsumowanie

Powinieneś teraz mieć kompletny pipeline, który:

- Przetwarza nazwy przez workflow powitań
- Przekazuje powitania ze znacznikami czasu do workflow'a transformacji
- Produkuje zarówno wersje wielkimi literami, jak i odwrócone powitań

---

## Podsumowanie

W tej misji pobocznej zbadaliśmy potężną koncepcję komponowania workflow'ów w Nextflow'ie, która pozwala nam budować złożone pipeline'y z mniejszych, wielokrotnego użytku komponentów.

To modularne podejście oferuje kilka zalet w porównaniu z monolitycznymi pipeline'ami:

- Każdy workflow może być rozwijany, testowany i debugowany niezależnie
- Workflow'y mogą być używane ponownie w różnych pipeline'ach
- Ogólna struktura pipeline'a staje się bardziej czytelna i łatwiejsza w utrzymaniu
- Zmiany w jednym workflow'ie niekoniecznie wpływają na inne, jeśli interfejsy pozostają spójne
- Punkty wejścia mogą być konfigurowane, aby uruchamiać różne części Twojego pipeline'a według potrzeb

_Ważne jest jednak zauważenie, że chociaż wywoływanie workflow'ów jest trochę podobne do wywoływania procesów, to nie jest to dokładnie to samo. Nie możesz na przykład uruchomić workflow'a N razy, wywołując go z kanałem o rozmiarze N - musiałbyś przekazać kanał o rozmiarze N do workflow'a i iterować wewnętrznie._

Zastosowanie tych technik w Twojej własnej pracy umożliwi Ci budowanie bardziej zaawansowanych pipeline'ów Nextflow'a, które mogą obsługiwać złożone zadania bioinformatyczne, pozostając jednocześnie łatwymi w utrzymaniu i skalowalnymi.

### Kluczowe wzorce

1.  **Struktura workflow'a**: Zdefiniowaliśmy jasne wejścia i wyjścia dla każdego workflow'a przy użyciu składni `take:` i `emit:`, tworząc dobrze zdefiniowane interfejsy między komponentami, i opakaliśmy logikę workflow'a w bloku `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **Importy workflow'ów:** Zbudowaliśmy dwa niezależne moduły workflow'ów i zaimportowaliśmy je do głównego pipeline'u za pomocą instrukcji `include`.

    - Zaimportuj pojedynczy workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Zaimportuj wiele workflow'ów

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Zaimportuj z aliasem, aby uniknąć konfliktów nazw

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Punkty wejścia**: Nextflow wymaga nienazwanego workflow'a wejściowego, aby wiedzieć, gdzie rozpocząć wykonywanie. Ten workflow wejściowy wywołuje Twoje nazwane workflow'y.

    - Nienazwany workflow (punkt wejścia)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Nazwany workflow (wywoływany z workflow'a wejściowego)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **Zarządzanie przepływem danych:** Nauczyliśmy się, jak uzyskiwać dostęp do wyjść workflow'a przy użyciu notacji przestrzeni nazw (`WORKFLOW_NAME.out.channel_name`) i przekazywać je do innych workflow'ów.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Dodatkowe zasoby

- [Dokumentacja workflow'ów Nextflow'a](https://www.nextflow.io/docs/latest/workflow.html)
- [Dokumentacja operatorów kanałów](https://www.nextflow.io/docs/latest/operator.html)
- [Dokumentacja strategii błędów](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Co dalej?

Wróć do [menu misji pobocznych](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
