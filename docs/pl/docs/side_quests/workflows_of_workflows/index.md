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

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow'a (procesy, kanały, operatory, moduły).

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś/zrobiłaś, otwórz środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../envsetup/index.md).

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

#### Przejrzyj materiały

Znajdziesz katalog `modules` zawierający kilka definicji procesów, które rozwijają to, czego nauczyłeś/nauczyłaś się w 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

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

## 1. Tworzenie Greeting Workflow

Zacznijmy od stworzenia workflow'u, który waliduje nazwy i generuje powitania ze znacznikami czasu.

### 1.1. Utwórz strukturę workflow'u

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Dodaj kod pierwszego (pod)workflow'u

Dodaj ten kod do `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Łańcuch procesów: walidacja -> tworzenie powitania -> dodanie znacznika czasu
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Jest to kompletny workflow o strukturze podobnej do tych, które widziałeś/widziałaś w samouczku 'Hello Nextflow', który możemy testować niezależnie. Sprawdźmy to teraz:

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

Działa zgodnie z oczekiwaniami, ale żeby uczynić go kompozytowalnym, musimy wprowadzić kilka zmian.

### 1.3. Uczyń workflow kompozytowalnym

Kompozytowalne workflow'y różnią się nieco od tych, które widziałeś/widziałaś w samouczku 'Hello Nextflow':

- Blok workflow musi mieć nazwę
- Wejścia są deklarowane przy użyciu słowa kluczowego `take:`
- Zawartość workflow'u jest umieszczona wewnątrz bloku `main:`
- Wyjścia są deklarowane przy użyciu słowa kluczowego `emit:`

Zaktualizujmy greeting workflow, aby pasował do tej struktury. Zmień kod na następujący:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Kanał wejściowy z nazwami

    main:
        // Łańcuch procesów: walidacja -> tworzenie powitania -> dodanie znacznika czasu
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Oryginalne powitania
        timestamped = timestamped_ch  // Powitania ze znacznikami czasu
}
```

Widać, że workflow ma teraz nazwę oraz bloki `take:` i `emit:` — to właśnie przez nie będziemy łączyć go z workflow'em wyższego poziomu.
Zawartość workflow'u jest również umieszczona wewnątrz bloku `main:`. Zwróć uwagę, że usunęliśmy deklarację kanału wejściowego `names_ch`, ponieważ jest on teraz przekazywany jako argument do workflow'u.

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

Ten komunikat informuje o nowym koncepcie — „entry workflow". Entry workflow to workflow, który jest wywoływany podczas uruchamiania skryptu Nextflow. Domyślnie Nextflow używa nienazwanego workflow'u jako entry workflow, gdy taki istnieje — i właśnie to robiłeś/robiłaś do tej pory, z blokami workflow zaczynającymi się tak:

```groovy title="hello.nf" linenums="1"
workflow {
```

Nasz greeting workflow nie ma jednak nienazwanego workflow'u — mamy zamiast tego workflow z nazwą:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Dlatego Nextflow zgłosił błąd i nie wykonał tego, czego oczekiwaliśmy.

Nie dodaliśmy składni `take:`/`emit:` po to, żeby wywoływać workflow bezpośrednio — zrobiliśmy to, żeby móc go komponować z innymi workflow'ami. Rozwiązaniem jest stworzenie głównego skryptu z nienazwanym entry workflow, który importuje i wywołuje nasz nazwany workflow.

### 1.4. Utwórz i przetestuj główny workflow

Teraz stworzymy główny workflow, który importuje i używa workflow'u `greeting`.

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

Zwróć uwagę, że entry workflow w tym pliku jest nienazwany — właśnie dlatego, że będziemy go używać jako entry workflow.

Uruchom go i sprawdź wynik:

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

Działa! Opakowano nazwany greeting workflow w główny workflow z nienazwanym blokiem `workflow` jako entry workflow. Główny workflow używa `GREETING_WORKFLOW` niemal (choć nie do końca) jak procesu i przekazuje kanał `names` jako argument.

### Podsumowanie

W tej sekcji poznałeś/poznałaś kilka ważnych konceptów:

- **Nazwane workflow'y**: Tworzenie nazwanego workflow'u (`GREETING_WORKFLOW`), który można importować i ponownie używać
- **Interfejsy workflow'ów**: Definiowanie wyraźnych wejść za pomocą `take:` i wyjść za pomocą `emit:` w celu stworzenia kompozytowalnego workflow'u
- **Entry points**: Zrozumienie, że Nextflow potrzebuje nienazwanego entry workflow, aby wiedzieć, od czego zacząć wykonanie skryptu
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

## 2. Dodanie Transform Workflow

Teraz stwórzmy workflow, który stosuje transformacje tekstu do powitań.

### 2.1. Utwórz plik workflow'u

```bash
touch workflows/transform.nf
```

### 2.2. Dodaj kod workflow'u

Dodaj ten kod do `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Kanał wejściowy z wiadomościami

    main:
        // Zastosuj transformacje sekwencyjnie
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Powitania wielkimi literami
        reversed = reversed_ch  // Odwrócone powitania wielkimi literami
}
```

Nie będziemy tu powtarzać wyjaśnienia składni kompozytowalnej, ale zwróć uwagę, że nazwany workflow jest ponownie zadeklarowany z blokami `take:` i `emit:`, a zawartość workflow'u jest umieszczona wewnątrz bloku `main:`.

### 2.3. Zaktualizuj główny workflow

Zaktualizuj `main.nf`, aby używał obu workflow'ów:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Uruchom greeting workflow
    GREETING_WORKFLOW(names)

    // Uruchom transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Wyświetl wyniki
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

Jeśli zajrzysz do jednego z tych odwróconych plików, zobaczysz, że jest to wersja powitania wielkimi literami, zapisana od tyłu:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

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

Stosowanie tych technik we własnej pracy pozwoli Ci budować bardziej zaawansowane pipeline'y Nextflow, które mogą obsługiwać złożone zadania bioinformatyczne, pozostając przy tym łatwymi w utrzymaniu i skalowalnymi.

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

Wróć do [menu Side Quests](../) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
