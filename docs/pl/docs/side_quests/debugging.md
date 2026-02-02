# Debugowanie Workflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Debugowanie to kluczowa umiejętność, która może zaoszczędzić Ci godzin frustracji i pomóc Ci stać się bardziej efektywnym programistą Nextflow. Przez całą Swoją karierę, szczególnie na początku, będziesz napotykać błędy podczas budowania i utrzymywania Swoich workflow. Nauka systematycznych podejść do debugowania pomoże Ci szybko identyfikować i rozwiązywać problemy.

### Cele nauczania

W tym side queście zbadamy **systematyczne techniki debugowania** dla workflow Nextflow:

- **Debugowanie błędów składni**: Efektywne wykorzystywanie funkcji IDE i komunikatów błędów Nextflow
- **Debugowanie kanałów**: Diagnozowanie problemów z przepływem danych i strukturą kanałów
- **Debugowanie procesów**: Badanie błędów wykonania i problemów z zasobami
- **Wbudowane narzędzia debugowania**: Wykorzystanie trybu podglądu, uruchamiania zaślepek i katalogów roboczych Nextflow
- **Systematyczne podejścia**: Czterofazowa metodologia efektywnego debugowania

Pod koniec będziesz mieć solidną metodologię debugowania, która przekształca frustrujące komunikaty błędów w jasne mapy drogowe do rozwiązań.

### Wymagania wstępne

Przed podjęciem tego side questa powinieneś:

- Ukończyć tutorial [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory)

**Opcjonalnie:** Zalecamy najpierw ukończenie side questa [IDE Features for Nextflow Development](./ide_features.md).
Obejmuje on kompleksowe omówienie funkcji IDE wspierających debugowanie (podświetlanie składni, wykrywanie błędów itp.), których będziemy tutaj intensywnie używać.

---

## 0. Rozpoczęcie pracy

#### Otwórz codespace szkoleniowy

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otwierasz środowisko szkoleniowe zgodnie z opisem w [Konfiguracja Środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki dla tego tutoriala.

```bash
cd side-quests/debugging
```

Możesz ustawić VSCode, aby skupić się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz zestaw przykładowych workflow z różnymi typami błędów, których użyjemy do ćwiczeń:

??? abstract "Zawartość katalogu"

    ```console
    .
    ├── bad_bash_var.nf
    ├── bad_channel_shape.nf
    ├── bad_channel_shape_viewed_debug.nf
    ├── bad_channel_shape_viewed.nf
    ├── bad_number_inputs.nf
    ├── badpractice_syntax.nf
    ├── bad_resources.nf
    ├── bad_syntax.nf
    ├── buggy_workflow.nf
    ├── data
    │   ├── sample_001.fastq.gz
    │   ├── sample_002.fastq.gz
    │   ├── sample_003.fastq.gz
    │   ├── sample_004.fastq.gz
    │   ├── sample_005.fastq.gz
    │   └── sample_data.csv
    ├── exhausted.nf
    ├── invalid_process.nf
    ├── missing_output.nf
    ├── missing_software.nf
    ├── missing_software_with_stub.nf
    ├── nextflow.config
    └── no_such_var.nf
    ```

Te pliki reprezentują typowe scenariusze debugowania, z którymi spotkasz się w rzeczywistym rozwoju.

#### Przejrzyj zadanie

Twoim wyzwaniem jest uruchomienie każdego workflow, zidentyfikowanie błędu(ów) i naprawienie ich.

Dla każdego błędnego workflow:

1. **Uruchom workflow** i obserwuj błąd
2. **Przeanalizuj komunikat błędu**: co Nextflow Ci mówi?
3. **Zlokalizuj problem** w kodzie używając dostarczonych wskazówek
4. **Napraw błąd** i zweryfikuj, że Twoje rozwiązanie działa
5. **Zresetuj plik** przed przejściem do następnej sekcji (użyj `git checkout <filename>`)

Ćwiczenia postępują od prostych błędów składni do bardziej subtelnych problemów w czasie wykonania.
Rozwiązania są omawiane bezpośrednio, ale spróbuj rozwiązać każdy samodzielnie przed czytaniem dalej.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy, aby się zanurzyć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace jest uruchomiony
- [ ] Ustawiłem odpowiednio mój katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Błędy składni

Błędy składni to najczęstszy typ błędów, z którymi się spotkasz podczas pisania kodu Nextflow. Występują, gdy kod nie jest zgodny z oczekiwanymi regułami składni DSL Nextflow. Te błędy uniemożliwiają uruchomienie Twojego workflow w ogóle, więc ważne jest, aby nauczyć się je szybko identyfikować i naprawiać.

### 1.1. Brakujące nawiasy klamrowe

Jednym z najczęstszych błędów składni, a czasem jednym z bardziej złożonych do debugowania, są **brakujące lub niedopasowane nawiasy**.

Zacznijmy od praktycznego przykładu.

#### Uruchom pipeline

```bash
nextflow run bad_syntax.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [stupefied_bhabha] DSL2 - revision: ca6327fad2

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

**Kluczowe elementy komunikatów błędów składni:**

- **Plik i lokalizacja**: Pokazuje, który plik i linia/kolumna zawierają błąd (`bad_syntax.nf:24:1`)
- **Opis błędu**: Wyjaśnia, co parser znalazł, czego się nie spodziewał (`Unexpected input: '<EOF>'`)
- **Wskaźnik EOF**: Komunikat `<EOF>` (End Of File) wskazuje, że parser dotarł do końca pliku, wciąż oczekując więcej treści - klasyczny znak niedomkniętych nawiasów klamrowych

#### Sprawdź kod

Teraz przeanalizujmy `bad_syntax.nf`, aby zrozumieć, co powoduje błąd:

```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
// Brakujący nawias zamykający dla procesu

workflow {

    // Utwórz kanał wejściowy
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Wywołaj proces z kanałem wejściowym
    PROCESS_FILES(input_ch)
}
```

Na potrzeby tego przykładu zostawiliśmy komentarz pokazujący, gdzie jest błąd. Rozszerzenie Nextflow VSCode powinno również dawać wskazówki o tym, co może być nie tak, podświetlając niedopasowany nawias na czerwono i zaznaczając przedwczesny koniec pliku:

![Bad syntax](img/bad_syntax.png)

**Strategia debugowania błędów nawiasów:**

1. Użyj dopasowywania nawiasów VS Code (umieść kursor obok nawiasu)
2. Sprawdź panel Problemów dla komunikatów związanych z nawiasami
3. Upewnij się, że każdy otwierający `{` ma odpowiadający zamykający `}`

#### Napraw kod

Zastąp komentarz brakującym nawiasem zamykającym:

=== "Po"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }  // Dodaj brakujący nawias zamykający

    workflow {

        // Utwórz kanał wejściowy
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Wywołaj proces z kanałem wejściowym
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="bad_syntax.nf" hl_lines="14" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    // Brakujący nawias zamykający dla procesu

    workflow {

        // Utwórz kanał wejściowy
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Wywołaj proces z kanałem wejściowym
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline

Teraz uruchom workflow ponownie, aby potwierdzić, że działa:

```bash
nextflow run bad_syntax.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [insane_faggin] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [48/cd7f54] PROCESS_FILES (1) | 3 of 3 ✔
    ```

### 1.2. Używanie nieprawidłowych słów kluczowych lub dyrektyw procesów

Kolejnym częstym błędem składni jest **nieprawidłowa definicja procesu**. Może to się zdarzyć, jeśli zapomnisz zdefiniować wymagane bloki lub użyjesz nieprawidłowych dyrektyw w definicji procesu.

#### Uruchom pipeline

```bash
nextflow run invalid_process.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [nasty_jepsen] DSL2 - revision: da9758d614

    Error invalid_process.nf:3:1: Invalid process definition -- check for missing or out-of-order section labels
    │   3 | process PROCESS_FILES {
    │     | ^^^^^^^^^^^^^^^^^^^^^^^
    │   4 |     inputs:
    │   5 |     val sample_name
    │   6 |
    ╰   7 |     output:

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Sprawdź kod

Błąd wskazuje "Invalid process definition" i pokazuje kontekst wokół problemu. Patrząc na linie 3-7, możemy zobaczyć `inputs:` w linii 4, co jest problemem. Przeanalizujmy `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // BŁĄD: Powinno być 'input' nie 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Utwórz kanał wejściowy
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Wywołaj proces z kanałem wejściowym
    PROCESS_FILES(input_ch)
}
```

Patrząc na linię 4 w kontekście błędu, możemy dostrzec problem: używamy `inputs` zamiast poprawnego `input`. Rozszerzenie Nextflow VSCode również oznaczy to:

![Invalid process message](img/invalid_process_message.png)

#### Napraw kod

Zastąp nieprawidłowe słowo kluczowe poprawnym, odwołując się do [dokumentacji](https://www.nextflow.io/docs/latest/process.html#):

=== "Po"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Naprawiono: Zmieniono 'inputs' na 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Utwórz kanał wejściowy
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Wywołaj proces z kanałem wejściowym
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // BŁĄD: Powinno być 'input' nie 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Utwórz kanał wejściowy
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Wywołaj proces z kanałem wejściowym
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline

Teraz uruchom workflow ponownie, aby potwierdzić, że działa:

```bash
nextflow run invalid_process.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `invalid_process.nf` [silly_fermi] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [b7/76cd9d] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.3. Używanie nieprawidłowych nazw zmiennych

Nazwy zmiennych używane w blokach script muszą być prawidłowe, pochodzące albo z wejść, albo z kodu groovy wstawionego przed skryptem. Ale kiedy opanujesz złożoność na początku rozwoju pipeline, łatwo jest popełnić błędy w nazewnictwie zmiennych, a Nextflow szybko Ci o tym powie.

#### Uruchom pipeline

```bash
nextflow run no_such_var.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [gloomy_meninsky] DSL2 - revision: 0c4d3bc28c

    Error no_such_var.nf:17:39: `undefined_var` is not defined
    │  17 |     echo "Using undefined variable: ${undefined_var}" >> ${output_pref
    ╰     |                                       ^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Błąd jest wykrywany podczas kompilacji i wskazuje bezpośrednio na niezdefiniowaną zmienną w linii 17, ze znakiem karetki wskazującym dokładnie, gdzie jest problem.

#### Sprawdź kod

Przeanalizujmy `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Zdefiniuj zmienne w kodzie Groovy przed skryptem
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // BŁĄD: undefined_var nie zdefiniowana
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Komunikat błędu wskazuje, że zmienna nie jest rozpoznawana w szablonie skryptu, i tak - powinieneś móc zobaczyć `${undefined_var}` użyte w bloku script, ale nie zdefiniowane nigdzie indziej.

#### Napraw kod

Jeśli otrzymasz błąd 'No such variable', możesz go naprawić, definiując zmienną (poprawiając nazwy zmiennych wejściowych lub edytując kod groovy przed skryptem), lub usuwając ją z bloku script, jeśli nie jest potrzebna:

=== "Po"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Zdefiniuj zmienne w kodzie Groovy przed skryptem
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Usunięto linię z undefined_var
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Zdefiniuj zmienne w kodzie Groovy przed skryptem
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // BŁĄD: undefined_var nie zdefiniowana
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline

Teraz uruchom workflow ponownie, aby potwierdzić, że działa:

```bash
nextflow run no_such_var.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `no_such_var.nf` [suspicious_venter] DSL2 - revision: 6ba490f7c5

    executor >  local (3)
    [21/237300] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 1.4. Nieprawidłowe użycie zmiennych Bash

Na początku w Nextflow może być trudno zrozumieć różnicę między zmiennymi Nextflow (Groovy) a Bash. Może to wygenerować inną formę błędu nieprawidłowej zmiennej, który pojawia się podczas próby użycia zmiennych w treści Bash bloku script.

#### Uruchom pipeline

```bash
nextflow run bad_bash_var.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [infallible_mandelbrot] DSL2 - revision: 0853c11080

    Error bad_bash_var.nf:13:42: `prefix` is not defined
    │  13 |     echo "Processing ${sample_name}" > ${prefix}.txt
    ╰     |                                          ^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Sprawdź kod

Błąd wskazuje na linię 13, gdzie używane jest `${prefix}`. Przeanalizujmy `bad_bash_var.nf`, aby zobaczyć, co powoduje problem:

```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    prefix="${sample_name}_output"
    echo "Processing ${sample_name}" > ${prefix}.txt  # BŁĄD: ${prefix} to składnia Groovy, nie Bash
    """
}
```

W tym przykładzie definiujemy zmienną `prefix` w Bash, ale w procesie Nextflow składnia `$` użyta do odwołania się do niej (`${prefix}`) jest interpretowana jako zmienna Groovy, a nie Bash. Zmienna nie istnieje w kontekście Groovy, więc otrzymujemy błąd 'no such variable'.

#### Napraw kod

Jeśli chcesz użyć zmiennej Bash, musisz zmienić znaczenie znaku dolara w ten sposób:

=== "Po"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > \${prefix}.txt  # Naprawiono: Zmieniono znaczenie znaku dolara
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="bad_bash_var.nf" hl_lines="13" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        prefix="${sample_name}_output"
        echo "Processing ${sample_name}" > ${prefix}.txt  # BŁĄD: ${prefix} to składnia Groovy, nie Bash
        """
    }
    ```

To mówi Nextflow, aby interpretować to jako zmienną Bash.

#### Uruchom pipeline

Teraz uruchom workflow ponownie, aby potwierdzić, że działa:

```bash
nextflow run bad_bash_var.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_bash_var.nf` [naughty_franklin] DSL2 - revision: 58c1c83709

    executor >  local (3)
    [4e/560285] PROCESS_FILES (2) | 3 of 3 ✔
    ```

!!! tip "Zmienne Groovy vs Bash"

    Dla prostych manipulacji zmiennymi, takich jak konkatenacja ciągów lub operacje na prefiksach/sufiksach, zwykle bardziej czytelne jest użycie zmiennych Groovy w sekcji script niż zmiennych Bash w bloku script:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    To podejście unika konieczności zmiany znaczenia znaków dolara i sprawia, że kod jest łatwiejszy do czytania i utrzymania.

### 1.5. Instrukcje poza blokiem workflow

Rozszerzenie Nextflow VSCode podświetla problemy ze strukturą kodu, które spowodują błędy. Częstym przykładem jest definiowanie kanałów poza blokiem `workflow {}` - jest to teraz wymuszane jako błąd składni.

#### Uruchom pipeline

```bash
nextflow run badpractice_syntax.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [intergalactic_colden] DSL2 - revision: 5e4b291bde

    Error badpractice_syntax.nf:3:1: Statements cannot be mixed with script declarations -- move statements into a process or workflow
    │   3 | input_ch = channel.of('sample1', 'sample2', 'sample3')
    ╰     | ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Komunikat błędu jasno wskazuje problem: instrukcje (takie jak definicje kanałów) nie mogą być mieszane z deklaracjami skryptu poza blokiem workflow lub process.

#### Sprawdź kod

Przeanalizujmy `badpractice_syntax.nf`, aby zobaczyć, co powoduje błąd:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // BŁĄD: Kanał zdefiniowany poza workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Zdefiniuj zmienne w kodzie Groovy przed skryptem
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {
    PROCESS_FILES(input_ch)
}
```

Rozszerzenie VSCode również podświetli zmienną `input_ch` jako zdefiniowaną poza blokiem workflow:

![Non-lethal syntax error](img/nonlethal.png)

#### Napraw kod

Przenieś definicję kanału do wnętrza bloku workflow:

=== "Po"

    ```groovy title="badpractice_syntax.nf" hl_lines="21" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Zdefiniuj zmienne w kodzie Groovy przed skryptem
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Przeniesiono do wnętrza bloku workflow
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // BŁĄD: Kanał zdefiniowany poza workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Zdefiniuj zmienne w kodzie Groovy przed skryptem
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline

Uruchom workflow ponownie, aby potwierdzić, że poprawka działa:

```bash
nextflow run badpractice_syntax.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `badpractice_syntax.nf` [naughty_ochoa] DSL2 - revision: 5e4b291bde

    executor >  local (3)
    [6a/84a608] PROCESS_FILES (2) | 3 of 3 ✔
    ```

Trzymaj Swoje kanały wejściowe zdefiniowane w bloku workflow i ogólnie postępuj zgodnie z wszelkimi innymi zaleceniami, które daje rozszerzenie.

### Wnioski

Możesz systematycznie identyfikować i naprawiać błędy składni używając komunikatów błędów Nextflow i wizualnych wskaźników IDE. Typowe błędy składni obejmują brakujące nawiasy klamrowe, nieprawidłowe słowa kluczowe procesów, niezdefiniowane zmienne i niewłaściwe użycie zmiennych Bash vs. Nextflow. Rozszerzenie VSCode pomaga wychwycić wiele z nich przed uruchomieniem. Mając te umiejętności debugowania składni w Swoim zestawie narzędzi, będziesz w stanie szybko rozwiązywać najczęstsze błędy składni Nextflow i przejść do rozwiązywania bardziej złożonych problemów w czasie wykonania.

### Co dalej?

Naucz się debugować bardziej złożone błędy struktury kanałów, które występują nawet gdy składnia jest poprawna.

---

## 2. Błędy struktury kanałów

Błędy struktury kanałów są bardziej subtelne niż błędy składni, ponieważ kod jest syntaktycznie poprawny, ale kształty danych nie pasują do tego, czego oczekują procesy. Nextflow spróbuje uruchomić pipeline, ale może stwierdzić, że liczba wejść nie pasuje do tego, czego oczekuje i zawieść. Te błędy zwykle pojawiają się tylko w czasie wykonania i wymagają zrozumienia danych przepływających przez Twój workflow.

!!! tip "Debugowanie kanałów za pomocą `.view()`"

    W tej sekcji pamiętaj, że możesz używać operatora `.view()` do inspekcji zawartości kanału w dowolnym momencie w Swoim workflow. Jest to jedno z najpotężniejszych narzędzi debugowania do zrozumienia problemów ze strukturą kanałów. Zbadamy tę technikę szczegółowo w sekcji 2.4, ale śmiało używaj jej podczas pracy nad przykładami.

    ```groovy
    my_channel.view()  // Pokazuje, co przepływa przez kanał
    ```

### 2.1. Zła liczba kanałów wejściowych

Ten błąd występuje, gdy przekazujesz inną liczbę kanałów niż oczekuje proces.

#### Uruchom pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [happy_swartz] DSL2 - revision: d83e58dcd3

    Error bad_number_inputs.nf:23:5: Incorrect number of call arguments, expected 1 but received 2
    │  23 |     PROCESS_FILES(samples_ch, files_ch)
    ╰     |     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

#### Sprawdź kod

Komunikat błędu jasno stwierdza, że wywołanie oczekiwało 1 argumentu, ale otrzymało 2, i wskazuje na linię 23. Przeanalizujmy `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Proces oczekuje tylko 1 wejścia

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Utwórz dwa oddzielne kanały
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // BŁĄD: Przekazywanie 2 kanałów, ale proces oczekuje tylko 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Powinieneś zobaczyć niedopasowane wywołanie `PROCESS_FILES`, dostarczające wiele kanałów wejściowych, gdy proces definiuje tylko jeden. Rozszerzenie VSCode również podkreśli wywołanie procesu na czerwono i dostarczy komunikat diagnostyczny przy najechaniu myszką:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Napraw kod

Dla tego konkretnego przykładu proces oczekuje pojedynczego kanału i nie wymaga drugiego kanału, więc możemy to naprawić przekazując tylko kanał `samples_ch`:

=== "Po"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Proces oczekuje tylko 1 wejścia

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Utwórz dwa oddzielne kanały
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Naprawiono: Przekaż tylko kanał, którego oczekuje proces
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Przed"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Proces oczekuje tylko 1 wejścia

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Utwórz dwa oddzielne kanały
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // BŁĄD: Przekazywanie 2 kanałów, ale proces oczekuje tylko 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Uruchom pipeline

```bash
nextflow run bad_number_inputs.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_number_inputs.nf` [big_euler] DSL2 - revision: e302bd87be

    executor >  local (3)
    [48/497f7b] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Częściej niż w tym przykładzie możesz dodać dodatkowe wejścia do procesu i zapomnieć odpowiednio zaktualizować wywołanie workflow, co może prowadzić do tego typu błędu. Na szczęście jest to jeden z łatwiejszych do zrozumienia i naprawienia błędów, ponieważ komunikat błędu jest dość jasny co do niedopasowania.

### 2.2. Wyczerpanie kanału (Proces uruchamia się rzadziej niż oczekiwano)

Niektóre błędy struktury kanałów są znacznie bardziej subtelne i nie powodują żadnych błędów. Prawdopodobnie najczęstszym z nich jest wyzwanie, przed którym stają nowi użytkownicy Nextflow w zrozumieniu, że kanały kolejkowe mogą zostać wyczerpane i zabraknie elementów, co oznacza, że workflow kończy się przedwcześnie.

#### Uruchom pipeline

```bash
nextflow run exhausted.nf
```

??? success "Wyjście polecenia"

```console title="Exhausted channel output"
 N E X T F L O W   ~  version 25.10.2

Launching `exhausted.nf` [extravagant_gauss] DSL2 - revision: 08cff7ba2a

executor >  local (1)
[bd/f61fff] PROCESS_FILES (1) [100%] 1 of 1 ✔
```

Ten workflow kończy się bez błędu, ale przetwarza tylko jedną próbkę!

#### Sprawdź kod

Przeanalizujmy `exhausted.nf`, aby zobaczyć, czy to jest w porządku:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Zdefiniuj zmienne w kodzie Groovy przed skryptem
    output_prefix = "${reference}_${sample_name}"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    """
}

workflow {

    reference_ch = channel.of('baseline_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

Proces uruchamia się tylko raz zamiast trzy razy, ponieważ kanał `reference_ch` jest kanałem kolejkowym, który zostaje wyczerpany po pierwszym wykonaniu procesu. Gdy jeden kanał jest wyczerpany, cały proces zatrzymuje się, nawet jeśli inne kanały nadal mają elementy.

To jest powszechny wzorzec, w którym masz pojedynczy plik referencyjny, który musi być ponownie używany dla wielu próbek. Rozwiązaniem jest przekształcenie kanału referencyjnego w kanał wartości, który może być używany w nieskończoność.

#### Napraw kod

Istnieje kilka sposobów rozwiązania tego problemu w zależności od liczby dotkniętych plików.

**Opcja 1**: Masz pojedynczy plik referencyjny, którego wielokrotnie używasz. Możesz po prostu utworzyć typ kanału wartości, który może być używany wielokrotnie. Istnieją trzy sposoby na to:

**1a** Użyj `channel.value()`:

```groovy title="exhausted.nf (naprawiony - Opcja 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Kanał wartości może być ponownie używany
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Użyj operatora `first()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (naprawiony - Opcja 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Konwertuj na kanał wartości
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Użyj operatora `collect()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (naprawiony - Opcja 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Konwertuj na kanał wartości
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opcja 2**: W bardziej złożonych scenariuszach, być może gdzie masz wiele plików referencyjnych dla wszystkich próbek w kanale próbek, możesz użyć operatora `combine` do utworzenia nowego kanału, który łączy dwa kanały w krotki:

```groovy title="exhausted.nf (naprawiony - Opcja 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Tworzy iloczyn kartezjański

    PROCESS_FILES(combined_ch)
}
```

Operator `.combine()` generuje iloczyn kartezjański dwóch kanałów, więc każdy element w `reference_ch` zostanie sparowany z każdym elementem w `input_ch`. To pozwala procesowi uruchamiać się dla każdej próbki, wciąż używając referencji.

To wymaga dostosowania wejścia procesu. W naszym przykładzie początek definicji procesu musiałby zostać dostosowany w następujący sposób:

```groovy title="exhausted.nf (naprawiony - Opcja 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

To podejście może nie być odpowiednie we wszystkich sytuacjach.

#### Uruchom pipeline

Wypróbuj jedno z powyższych poprawek i uruchom workflow ponownie:

```bash
nextflow run exhausted.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `exhausted.nf` [maniac_leavitt] DSL2 - revision: f372a56a7d

    executor >  local (3)
    [80/0779e9] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Powinieneś teraz zobaczyć, że wszystkie trzy próbki są przetwarzane zamiast tylko jednej.

### 2.3. Zła struktura zawartości kanału

Gdy workflow osiągną pewien poziom złożoności, może być trudno śledzić wewnętrzne struktury każdego kanału, a ludzie często generują niedopasowania między tym, czego oczekuje proces, a tym, co faktycznie zawiera kanał. Jest to bardziej subtelne niż problem, który omówiliśmy wcześniej, gdzie liczba kanałów była nieprawidłowa. W tym przypadku możesz mieć poprawną liczbę kanałów wejściowych, ale wewnętrzna struktura jednego lub więcej z tych kanałów nie pasuje do tego, czego oczekuje proces.

#### Uruchom pipeline

```bash
nextflow run bad_channel_shape.nf
```

??? failure "Wyjście polecenia"

    ```console
    Launching `bad_channel_shape.nf` [hopeful_pare] DSL2 - revision: ffd66071a1

    executor >  local (3)
    executor >  local (3)
    [3f/c2dcb3] PROCESS_FILES (3) [  0%] 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `[sample1, file1.txt]_output.txt` expected by process `PROCESS_FILES (1)`


    Command executed:

      echo "Processing [sample1, file1.txt]" > [sample1, file1.txt]_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/d6/1fb69d1d93300bbc9d42f1875b981e

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Sprawdź kod

Nawiasy kwadratowe w komunikacie błędu dostarczają wskazówki - proces traktuje krotkę jako pojedynczą wartość, co nie jest tym, czego chcemy. Przeanalizujmy `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Oczekuje pojedynczej wartości, dostaje krotkę

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Możesz zobaczyć, że generujemy kanał złożony z krotek: `['sample1', 'file1.txt']`, ale proces oczekuje pojedynczej wartości, `val sample_name`. Wykonane polecenie pokazuje, że proces próbuje utworzyć plik o nazwie `[sample3, file3.txt]_output.txt`, co nie jest zamierzonym wyjściem.

#### Napraw kod

Aby to naprawić, jeśli proces wymaga obu wejść, możemy dostosować proces, aby akceptował krotkę:

=== "Opcja 1: Akceptuj krotkę w procesie"

    === "Po"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Naprawiono: Akceptuj krotkę

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

    === "Przed"

        ```groovy title="bad_channel_shape.nf" hl_lines="5" linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                val sample_name  // Oczekuje pojedynczej wartości, dostaje krotkę

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Opcja 2: Wydobądź pierwszy element"

    === "Po"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Naprawiono: Wydobądź pierwszy element
        }
        ```

    === "Przed"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Uruchom pipeline

Wybierz jedno z rozwiązań i uruchom ponownie workflow:

```bash
nextflow run bad_channel_shape.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape.nf` [clever_thompson] DSL2 - revision: 8cbcae3746

    executor >  local (3)
    [bb/80a958] PROCESS_FILES (2) | 3 of 3 ✔
    ```

### 2.4. Techniki debugowania kanałów

#### Używanie `.view()` do inspekcji kanałów

Najpotężniejszym narzędziem debugowania dla kanałów jest operator `.view()`. Dzięki `.view()` możesz zrozumieć kształt Swoich kanałów na wszystkich etapach, aby pomóc w debugowaniu.

#### Uruchom pipeline

Uruchom `bad_channel_shape_viewed.nf`, aby to zobaczyć w akcji:

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [maniac_poisson] DSL2 - revision: b4f24dc9da

    executor >  local (3)
    [c0/db76b3] PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

#### Sprawdź kod

Przeanalizujmy `bad_channel_shape_viewed.nf`, aby zobaczyć, jak używane jest `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debugowanie: Pokaż oryginalną zawartość kanału
    .map { tuple -> tuple[0] }        // Transformacja: Wydobądź pierwszy element
    .view { "After mapping: $it" }    // Debugowanie: Pokaż przekształconą zawartość kanału

    PROCESS_FILES(input_ch)
}
```

#### Napraw kod

Aby oszczędzić sobie nadmiernego używania operacji `.view()` w przyszłości do zrozumienia zawartości kanału, wskazane jest dodanie kilku komentarzy pomocnych:

```groovy title="bad_channel_shape_viewed.nf (z komentarzami)" linenums="16" hl_lines="8 9"
workflow {

    // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Będzie to ważniejsze w miarę jak Twoje workflow będą rosły w złożoności i struktura kanału stanie się bardziej nieprzejrzysta.

#### Uruchom pipeline

```bash
nextflow run bad_channel_shape_viewed.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed.nf` [marvelous_koch] DSL2 - revision: 03e79cdbad

    executor >  local (3)
    [ff/d67cec] PROCESS_FILES (2) | 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    ```

### Wnioski

Wiele błędów struktury kanałów może być utworzonych za pomocą prawidłowej składni Nextflow. Możesz debugować błędy struktury kanałów poprzez zrozumienie przepływu danych, używanie operatorów `.view()` do inspekcji i rozpoznawanie wzorców komunikatów błędów, takich jak nawiasy kwadratowe wskazujące nieoczekiwane struktury krotek.

### Co dalej?

Dowiedz się o błędach tworzonych przez definicje procesów.

---

## 3. Błędy struktury procesów

Większość błędów, które napotykasz związanych z procesami, będzie związana z błędami, które popełniłeś w formowaniu polecenia, lub z problemami związanymi z podstawowym oprogramowaniem. Mimo to, podobnie jak w przypadku problemów z kanałami powyżej, możesz popełnić błędy w definicji procesu, które nie kwalifikują się jako błędy składni, ale które spowodują błędy w czasie wykonania.

### 3.1. Brakujące pliki wyjściowe

Jednym z częstych błędów podczas pisania procesów jest zrobienie czegoś, co generuje niedopasowanie między tym, czego oczekuje proces, a tym, co jest generowane.

#### Uruchom pipeline

```bash
nextflow run missing_output.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [zen_stone] DSL2 - revision: 37ff61f926

    executor >  local (3)
    executor >  local (3)
    [fd/2642e9] process > PROCESS_FILES (2) [ 66%] 2 of 3, failed: 2
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Missing output file(s) `sample3.txt` expected by process `PROCESS_FILES (3)`


    Command executed:

      echo "Processing sample3" > sample3_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/02/9604d49fb8200a74d737c72a6c98ed

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

#### Sprawdź kod

Komunikat błędu wskazuje, że proces oczekiwał utworzenia pliku wyjściowego o nazwie `sample3.txt`, ale skrypt faktycznie tworzy `sample3_output.txt`. Przeanalizujmy definicję procesu w `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Oczekuje: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Tworzy: sample3_output.txt
    """
}
```

Powinieneś zobaczyć, że istnieje niedopasowanie między nazwą pliku wyjściowego w bloku `output:`, a tą użytą w skrypcie. To niedopasowanie powoduje, że proces zawodzi. Jeśli napotkasz tego rodzaju błąd, wróć i sprawdź, czy wyjścia pasują między Twoją definicją procesu a Twoim blokiem wyjściowym.

Jeśli problem nadal nie jest jasny, sprawdź sam katalog roboczy, aby zidentyfikować faktycznie utworzone pliki wyjściowe:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

Dla tego przykładu to podkreśliłoby dla nas, że sufiks `_output` jest włączany do nazwy pliku wyjściowego, wbrew naszej definicji `output:`.

#### Napraw kod

Napraw niedopasowanie, czyniąc nazwę pliku wyjściowego spójną:

=== "Po"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Naprawiono: Dopasuj wyjście skryptu

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }
    ```

=== "Przed"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_
