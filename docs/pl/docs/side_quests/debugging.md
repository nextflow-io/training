# Debugowanie Workflow'ów

Debugowanie to kluczowa umiejętność, która może zaoszczędzić Ci godzin frustracji i pomóc stać się bardziej efektywnym programistą Nextflow. W trakcie swojej kariery, szczególnie na początku, będziesz napotykać błędy podczas tworzenia i utrzymywania workflow'ów. Nauka systematycznych podejść do debugowania pomoże Ci szybko identyfikować i rozwiązywać problemy.

### Cele edukacyjne

W tym zadaniu pobocznym zbadamy **systematyczne techniki debugowania** dla workflow'ów Nextflow:

- **Debugowanie błędów składni**: Efektywne wykorzystanie funkcji IDE i komunikatów błędów Nextflow
- **Debugowanie kanałów**: Diagnozowanie problemów z przepływem danych i strukturą kanałów
- **Debugowanie procesów**: Badanie błędów wykonania i problemów z zasobami
- **Wbudowane narzędzia debugowania**: Wykorzystanie trybu podglądu, uruchamiania zastępczego i katalogów roboczych Nextflow
- **Systematyczne podejścia**: Czterofazowa metodologia efektywnego debugowania

Pod koniec będziesz dysponować solidną metodologią debugowania, która przekształca frustrujące komunikaty błędów w przejrzyste mapy drogowe do rozwiązań.

### Wymagania wstępne

Przed podjęciem tego zadania pobocznego powinieneś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących
- Swobodnie posługiwać się podstawowymi koncepcjami i mechanizmami Nextflow (procesy, kanały, operatory)

**Opcjonalnie:** Zalecamy najpierw ukończenie zadania pobocznego [IDE Features for Nextflow Development](./ide_features.md).
Obejmuje ono kompleksowe omówienie funkcji IDE wspierających debugowanie (podświetlanie składni, wykrywanie błędów itp.), które będziemy tu intensywnie wykorzystywać.

---

## 0. Rozpocznij pracę

#### Otwórz środowisko szkoleniowe w codespace

Jeśli jeszcze tego nie zrobiłeś, upewnij się, że otworzyłeś środowisko szkoleniowe zgodnie z opisem w [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki dla tego samouczka.

```bash
cd side-quests/debugging
```

Możesz ustawić VSCode, aby skupił się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz zestaw przykładowych workflow'ów z różnymi typami błędów, których użyjemy do ćwiczeń:

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

Te pliki reprezentują typowe scenariusze debugowania, które napotkasz w rzeczywistym rozwoju.

#### Przejrzyj zadanie

Twoim wyzwaniem jest uruchomienie każdego workflow'a, zidentyfikowanie błędu (błędów) i ich naprawienie.

Dla każdego błędnego workflow'a:

1. **Uruchom workflow'a** i obserwuj błąd
2. **Przeanalizuj komunikat błędu**: co mówi Ci Nextflow?
3. **Zlokalizuj problem** w kodzie, korzystając z dostarczonych wskazówek
4. **Napraw błąd** i zweryfikuj, że Twoje rozwiązanie działa
5. **Zresetuj plik** przed przejściem do następnej sekcji (użyj `git checkout <nazwa_pliku>`)

Ćwiczenia przechodzą od prostych błędów składni do bardziej subtelnych problemów w czasie wykonania.
Rozwiązania są omawiane w tekście, ale spróbuj rozwiązać każde samodzielnie przed czytaniem dalej.

#### Lista gotowości

Myślisz, że jesteś gotowy, aby zacząć?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Mój codespace działa
- [ ] Ustawiłem odpowiednio swój katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Błędy składni

Błędy składni to najczęstszy typ błędów, które napotkasz podczas pisania kodu Nextflow. Występują, gdy kod nie jest zgodny z oczekiwanymi regułami składni DSL Nextflow. Te błędy uniemożliwiają uruchomienie Twojego workflow'a, więc ważne jest, aby nauczyć się je szybko identyfikować i naprawiać.

### 1.1. Brakujące nawiasy klamrowe

Jednym z najczęstszych błędów składni, a czasem jednym z bardziej złożonych do debugowania, są **brakujące lub niedopasowane nawiasy**.

Zacznijmy od praktycznego przykładu.

#### Uruchom pipeline'a

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

- **Plik i lokalizacja**: Pokazuje, który plik i która linia/kolumna zawiera błąd (`bad_syntax.nf:24:1`)
- **Opis błędu**: Wyjaśnia, co parser znalazł, czego się nie spodziewał (`Unexpected input: '<EOF>'`)
- **Wskaźnik EOF**: Komunikat `<EOF>` (End Of File) wskazuje, że parser osiągnął koniec pliku, wciąż oczekując więcej treści - klasyczny znak niezamkniętych nawiasów klamrowych

#### Sprawdź kod

Teraz zbadajmy `bad_syntax.nf`, aby zrozumieć, co powoduje błąd:

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
// Missing closing brace for the process

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Na potrzeby tego przykładu zostawiliśmy komentarz pokazujący, gdzie jest błąd. Rozszerzenie Nextflow VSCode powinno również dawać Ci wskazówki, co może być nie tak, zaznaczając niedopasowany nawias na czerwono i podświetlając przedwczesny koniec pliku:

![Bad syntax](img/bad_syntax.png)

**Strategia debugowania błędów nawiasów:**

1. Użyj dopasowywania nawiasów VS Code (umieść kursor obok nawiasu)
2. Sprawdź panel Problemów pod kątem komunikatów związanych z nawiasami
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
    }  // Add the missing closing brace

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
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
    // Missing closing brace for the process

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline'a

Teraz uruchom workflow'a ponownie, aby potwierdzić, że działa:

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

### 1.2. Używanie nieprawidłowych słów kluczowych lub dyrektyw procesu

Innym częstym błędem składni jest **nieprawidłowa definicja procesu**. Może się to zdarzyć, jeśli zapomnisz zdefiniować wymagane bloki lub użyjesz nieprawidłowych dyrektyw w definicji procesu.

#### Uruchom pipeline'a

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

Błąd wskazuje "Invalid process definition" i pokazuje kontekst wokół problemu. Patrząc na linie 3-7, widzimy `inputs:` w linii 4, co jest problemem. Zbadajmy `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // ERROR: Should be 'input' not 'inputs'
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
```

Patrząc na linię 4 w kontekście błędu, możemy dostrzec problem: używamy `inputs` zamiast poprawnego `input`. Rozszerzenie Nextflow VSCode również to oznaczy:

![Invalid process message](img/invalid_process_message.png)

#### Napraw kod

Zastąp nieprawidłowe słowo kluczowe poprawnym, odwołując się do [dokumentacji](https://www.nextflow.io/docs/latest/process.html#):

=== "Po"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Fixed: Changed 'inputs' to 'input'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        inputs:  // ERROR: Should be 'input' not 'inputs'
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create input channel
        input_ch = channel.of('sample1', 'sample2', 'sample3')

        // Call the process with the input channel
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline'a

Teraz uruchom workflow'a ponownie, aby potwierdzić, że działa:

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

### 1.3. Używanie błędnych nazw zmiennych

Nazwy zmiennych używane w blokach skryptów muszą być prawidłowe, pochodzące albo z wejść, albo z kodu groovy wstawionego przed skryptem. Ale gdy zmagasz się ze złożonością na początku rozwoju pipeline'a, łatwo popełnić błędy w nazewnictwie zmiennych, a Nextflow szybko Cię o tym poinformuje.

#### Uruchom pipeline'a

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

Błąd jest wykrywany w czasie kompilacji i wskazuje bezpośrednio na niezdefiniowaną zmienną w linii 17, z daszkiem wskazującym dokładnie, gdzie jest problem.

#### Sprawdź kod

Zbadajmy `no_such_var.nf`:

```groovy title="no_such_var.nf" hl_lines="17" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
    def output_prefix = "${sample_name}_processed"
    def timestamp = new Date().format("yyyy-MM-dd")

    """
    echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Komunikat błędu wskazuje, że zmienna nie jest rozpoznawana w szablonie skryptu, i proszę bardzo - powinieneś widzieć `${undefined_var}` użyte w bloku skryptu, ale nie zdefiniowane gdzie indziej.

#### Napraw kod

Jeśli otrzymasz błąd 'No such variable', możesz go naprawić, definiując zmienną (poprawiając nazwy zmiennych wejściowych lub edytując kod groovy przed skryptem) lub usuwając ją z bloku skryptu, jeśli nie jest potrzebna:

=== "Po"

    ```groovy title="no_such_var.nf" hl_lines="15-17" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """  // Removed the line with undefined_var
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
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // ERROR: undefined_var not defined
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline'a

Teraz uruchom workflow'a ponownie, aby potwierdzić, że działa:

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

### 1.4. Złe użycie zmiennych Bash

Na początku pracy z Nextflow może być trudno zrozumieć różnicę między zmiennymi Nextflow (Groovy) a Bash. Może to generować inną formę błędu złej zmiennej, która pojawia się przy próbie użycia zmiennych w treści Bash bloku skryptu.

#### Uruchom pipeline'a

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

Błąd wskazuje na linię 13, gdzie używane jest `${prefix}`. Zbadajmy `bad_bash_var.nf`, aby zobaczyć, co powoduje problem:

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

W tym przykładzie definiujemy zmienną `prefix` w Bash, ale w procesie Nextflow składnia `$`, której użyliśmy do odwołania się do niej (`${prefix}`), jest interpretowana jako zmienna Groovy, a nie Bash. Zmienna nie istnieje w kontekście Groovy, więc otrzymujemy błąd 'no such variable'.

#### Napraw kod

Jeśli chcesz użyć zmiennej Bash, musisz poprzedzić znak dolara odwrotnym ukośnikiem w ten sposób:

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Naprawiono: Poprzedzono znak dolara odwrotnym ukośnikiem
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

To mówi Nextflow, aby interpretował to jako zmienną Bash.

#### Uruchom pipeline'a

Teraz uruchom workflow'a ponownie, aby potwierdzić, że działa:

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

    Dla prostych manipulacji zmiennymi, takich jak konkatenacja łańcuchów lub operacje prefiksu/sufiksu, zazwyczaj bardziej czytelne jest użycie zmiennych Groovy w sekcji skryptu niż zmiennych Bash w bloku skryptu:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    To podejście unika konieczności poprzedzania znaków dolara odwrotnym ukośnikiem i sprawia, że kod jest łatwiejszy do odczytania i utrzymania.

### 1.5. Instrukcje poza blokiem workflow

Rozszerzenie Nextflow VSCode podświetla problemy ze strukturą kodu, które spowodują błędy. Częstym przykładem jest definiowanie kanałów poza blokiem `workflow {}` - jest to teraz egzekwowane jako błąd składni.

#### Uruchom pipeline'a

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

Komunikat błędu jasno wskazuje problem: instrukcje (takie jak definicje kanałów) nie mogą być mieszane z deklaracjami skryptów poza blokiem workflow lub procesu.

#### Sprawdź kod

Zbadajmy `badpractice_syntax.nf`, aby zobaczyć, co powoduje błąd:

```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
#!/usr/bin/env nextflow

input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}_processed.txt"

    script:
    // Define variables in Groovy code before the script
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
        // Define variables in Groovy code before the script
        def output_prefix = "${sample_name}_processed"
        def timestamp = new Date().format("yyyy-MM-dd")

        """
        echo "Processing ${sample_name} on ${timestamp}" > ${output_prefix}.txt
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Moved inside workflow block
        PROCESS_FILES(input_ch)
    }
    ```

=== "Przed"

    ```groovy title="badpractice_syntax.nf" hl_lines="3" linenums="1"
    #!/usr/bin/env nextflow

    input_ch = channel.of('sample1', 'sample2', 'sample3')  // ERROR: Channel defined outside workflow

    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_processed.txt"

        script:
        // Define variables in Groovy code before the script
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

#### Uruchom pipeline'a

Uruchom workflow'a ponownie, aby potwierdzić, że poprawka działa:

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

Trzymaj swoje kanały wejściowe zdefiniowane w bloku workflow i ogólnie postępuj zgodnie z innymi zaleceniami rozszerzenia.

### Podsumowanie

Możesz systematycznie identyfikować i naprawiać błędy składni, używając komunikatów błędów Nextflow i wizualnych wskaźników IDE. Typowe błędy składni obejmują brakujące nawiasy klamrowe, nieprawidłowe słowa kluczowe procesu, niezdefiniowane zmienne i niewłaściwe użycie zmiennych Bash vs Nextflow. Rozszerzenie VSCode pomaga wychwycić wiele z nich przed uruchomieniem. Dzięki tym umiejętnościom debugowania składni będziesz w stanie szybko rozwiązywać najczęstsze błędy składni Nextflow i przejść do rozwiązywania bardziej złożonych problemów w czasie wykonania.

### Co dalej?

Naucz się debugować bardziej złożone błędy struktury kanałów, które występują nawet gdy składnia jest poprawna.

---

## 2. Błędy struktury kanałów

Błędy struktury kanałów są bardziej subtelne niż błędy składni, ponieważ kod jest składniowo poprawny, ale kształty danych nie pasują do tego, czego oczekują procesy. Nextflow spróbuje uruchomić pipeline'a, ale może stwierdzić, że liczba wejść nie pasuje do oczekiwań i zakończyć się niepowodzeniem. Te błędy zazwyczaj pojawiają się tylko w czasie wykonania i wymagają zrozumienia danych przepływających przez Twój workflow.

!!! tip "Debugowanie kanałów za pomocą `.view()`"

    W tej sekcji pamiętaj, że możesz użyć operatora `.view()` do inspekcji zawartości kanału w dowolnym punkcie Twojego workflow'a. To jedno z najpotężniejszych narzędzi debugowania do zrozumienia problemów ze strukturą kanałów. Zbadamy tę technikę szczegółowo w sekcji 2.4, ale możesz jej używać podczas pracy nad przykładami.

    ```groovy
    my_channel.view()  // Shows what's flowing through the channel
    ```

### 2.1. Zła liczba kanałów wejściowych

Ten błąd występuje, gdy przekazujesz inną liczbę kanałów niż oczekuje proces.

#### Uruchom pipeline'a

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

Komunikat błędu jasno stwierdza, że wywołanie oczekiwało 1 argumentu, ale otrzymało 2, i wskazuje na linię 23. Zbadajmy `bad_number_inputs.nf`:

```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
```

Powinieneś zobaczyć niedopasowane wywołanie `PROCESS_FILES`, dostarczające wiele kanałów wejściowych, gdy proces definiuje tylko jeden. Rozszerzenie VSCode również podkreśli wywołanie procesu na czerwono i dostarczy komunikat diagnostyczny po najechaniu myszą:

![Incorrect number of args message](img/incorrect_num_args.png)

#### Napraw kod

W tym konkretnym przykładzie proces oczekuje pojedynczego kanału i nie wymaga drugiego kanału, więc możemy to naprawić, przekazując tylko kanał `samples_ch`:

=== "Po"

    ```groovy title="bad_number_inputs.nf" hl_lines="23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // Fixed: Pass only the channel the process expects
        PROCESS_FILES(samples_ch)
    }
    ```

=== "Przed"

    ```groovy title="bad_number_inputs.nf" hl_lines="5 23" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:
            val sample_name  // Process expects only 1 input

        output:
            path "${sample_name}_output.txt"

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt
        """
    }

    workflow {

        // Create two separate channels
        samples_ch = channel.of('sample1', 'sample2', 'sample3')
        files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

        // ERROR: Passing 2 channels but process expects only 1
        PROCESS_FILES(samples_ch, files_ch)
    }
    ```

#### Uruchom pipeline'a

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

Częściej niż w tym przykładzie możesz dodać dodatkowe wejścia do procesu i zapomnieć odpowiednio zaktualizować wywołanie workflow'a, co może prowadzić do tego typu błędu. Na szczęście jest to jeden z łatwiejszych do zrozumienia i naprawienia błędów, ponieważ komunikat błędu jest dość jasny co do niedopasowania.

### 2.2. Wyczerpanie kanału (proces działa mniej razy niż oczekiwano)

Niektóre błędy struktury kanałów są znacznie bardziej subtelne i nie generują żadnych błędów. Prawdopodobnie najczęstszy z nich odzwierciedla wyzwanie, przed którym stają nowi użytkownicy Nextflow w zrozumieniu, że kanały kolejki mogą się wyczerpać i wyczerpać elementy, co oznacza, że workflow kończy się przedwcześnie.

#### Uruchom pipeline'a

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

Zbadajmy `exhausted.nf`, aby zobaczyć, czy to jest poprawne:

```groovy title="exhausted.nf" hl_lines="23 24" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
    val reference
    val sample_name

    output:
    path "${output_prefix}.txt"

    script:
    // Define variables in Groovy code before the script
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

Proces działa tylko raz zamiast trzy razy, ponieważ kanał `reference_ch` jest kanałem kolejki, który wyczerpuje się po pierwszym wykonaniu procesu. Gdy jeden kanał się wyczerpie, cały proces zatrzymuje się, nawet jeśli inne kanały wciąż mają elementy.

To jest typowy wzorzec, w którym masz pojedynczy plik referencyjny, który musi być ponownie używany dla wielu próbek. Rozwiązaniem jest przekształcenie kanału referencyjnego w kanał wartości, który może być używany w nieskończoność.

#### Napraw kod

Istnieje kilka sposobów rozwiązania tego problemu w zależności od tego, ile plików jest dotkniętych.

**Opcja 1**: Masz pojedynczy plik referencyjny, który wielokrotnie używasz. Możesz po prostu utworzyć typ kanału wartości, który może być używany wielokrotnie. Istnieją trzy sposoby, aby to zrobić:

**1a** Użyj `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Value channel can be reused
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Użyj operatora `first()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#first):

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Użyj operatora `collect()` [operator](https://www.nextflow.io/docs/latest/reference/operator.html#collect):

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Convert to value channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opcja 2**: W bardziej złożonych scenariuszach, być może gdy masz wiele plików referencyjnych dla wszystkich próbek w kanale próbek, możesz użyć operatora `combine`, aby utworzyć nowy kanał łączący dwa kanały w krotki:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Creates cartesian product

    PROCESS_FILES(combined_ch)
}
```

Operator `.combine()` generuje iloczyn kartezjański dwóch kanałów, więc każdy element w `reference_ch` zostanie sparowany z każdym elementem w `input_ch`. To pozwala procesowi działać dla każdej próbki, wciąż używając referencji.

Wymaga to dostosowania wejścia procesu. W naszym przykładzie początek definicji procesu musiałby zostać dostosowany w następujący sposób:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

To podejście może nie być odpowiednie we wszystkich sytuacjach.

#### Uruchom pipeline'a

Wypróbuj jedno z powyższych rozwiązań i uruchom workflow'a ponownie:

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

Powinieneś teraz zobaczyć wszystkie trzy próbki przetwarzane zamiast tylko jednej.

### 2.3. Zła struktura zawartości kanału

Gdy workflow'y osiągają pewien poziom złożoności, może być trochę trudno śledzić wewnętrzne struktury każdego kanału, a ludzie często generują niedopasowania między tym, czego oczekuje proces, a tym, co kanał faktycznie zawiera. Jest to bardziej subtelne niż problem, który omówiliśmy wcześniej, gdzie liczba kanałów była nieprawidłowa. W tym przypadku możesz mieć prawidłową liczbę kanałów wejściowych, ale wewnętrzna struktura jednego lub więcej z tych kanałów nie pasuje do tego, czego oczekuje proces.

#### Uruchom pipeline'a

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

Nawiasy kwadratowe w komunikacie błędu dostarczają wskazówki - proces traktuje krotkę jako pojedynczą wartość, co nie jest tym, czego chcemy. Zbadajmy `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Expects single value, gets tuple

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    PROCESS_FILES(input_ch)
}
```

Widzisz, że generujemy kanał złożony z krotek: `['sample1', 'file1.txt']`, ale proces oczekuje pojedynczej wartości, `val sample_name`. Wykonane polecenie pokazuje, że proces próbuje utworzyć plik o nazwie `[sample3, file3.txt]_output.txt`, co nie jest zamierzonym wyjściem.

#### Napraw kod

Aby to naprawić, jeśli proces wymaga obu wejść, moglibyśmy dostosować proces, aby akceptował krotkę:

=== "Opcja 1: Akceptuj krotkę w procesie"

    === "Po"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Fixed: Accept tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
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
                val sample_name  // Expects single value, gets tuple

            output:
                path "${sample_name}_output.txt"

            script:
            """
            echo "Processing ${sample_name}" > ${sample_name}_output.txt
            """
        }

        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

=== "Opcja 2: Wyodrębnij pierwszy element"

    === "Po"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Fixed: Extract first element
        }
        ```

    === "Przed"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Channel emits tuples, but process expects single values
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch)
        }
        ```

#### Uruchom pipeline'a

Wybierz jedno z rozwiązań i uruchom workflow'a ponownie:

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

Najpotężniejszym narzędziem debugowania dla kanałów jest operator `.view()`. Dzięki `.view()` możesz zrozumieć kształt swoich kanałów na wszystkich etapach, aby pomóc w debugowaniu.

#### Uruchom pipeline'a

Uruchom `bad_channel_shape_viewed.nf`, aby zobaczyć to w akcji:

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

Zbadajmy `bad_channel_shape_viewed.nf`, aby zobaczyć, jak używane jest `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Show original channel content
    .map { tuple -> tuple[0] }        // Transform: Extract first element
    .view { "After mapping: $it" }    // Debug: Show transformed channel content

    PROCESS_FILES(input_ch)
}
```

#### Napraw kod

Aby uchronić się przed nadmiernym używaniem operacji `.view()` w przyszłości do zrozumienia zawartości kanału, wskazane jest dodanie komentarzy pomocniczych:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
workflow {

    // Channel emits tuples, but process expects single values
    input_ch = channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        ) // [sample_name, file_name]
        .map { tuple -> tuple[0] } // sample_name

    PROCESS_FILES(input_ch)
}
```

Stanie się to ważniejsze, gdy Twoje workflow'y będą rosły w złożoności, a struktura kanałów stanie się bardziej nieprzejrzysta.

#### Uruchom pipeline'a

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

### Podsumowanie

Wiele błędów struktury kanałów można utworzyć przy użyciu prawidłowej składni Nextflow. Możesz debugować błędy struktury kanałów, rozumiejąc przepływ danych, używając operatorów `.view()` do inspekcji i rozpoznając wzorce błędów, takie jak nawiasy kwadratowe wskazujące nieoczekiwane struktury krotek.

### Co dalej?

Dowiedz się o błędach tworzonych przez definicje procesów.

---

## 3. Błędy struktury procesów

Większość błędów, które napotkasz związanych z procesami, będzie związana z błędami popełnionymi podczas formowania polecenia lub z problemami związanymi z podstawowym oprogramowaniem. To powiedziawszy, podobnie jak w przypadku problemów z kanałami powyżej, możesz popełnić błędy w definicji procesu, które nie kwalifikują się jako błędy składni, ale które spowodują błędy w czasie wykonania.

### 3.1. Brakujące pliki wyjściowe

Jednym z częstych błędów podczas pisania procesów jest zrobienie czegoś, co generuje niedopasowanie między tym, czego oczekuje proces, a tym, co jest generowane.

#### Uruchom pipeline'a

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

Komunikat błędu wskazuje, że proces oczekiwał utworzenia pliku wyjściowego o nazwie `sample3.txt`, ale skrypt faktycznie tworzy `sample3_output.txt`. Zbadajmy definicję procesu w `missing_output.nf`:

```groovy title="missing_output.nf" linenums="3" hl_lines="6 10"
process PROCESS_FILES {
    input:
    val sample_name

    output:
    path "${sample_name}.txt"  // Expects: sample3.txt

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
    """
}
```

Powinieneś zobaczyć, że istnieje niedopasowanie między nazwą pliku wyjściowego w bloku `output:` a tą użytą w skrypcie. To niedopasowanie powoduje niepowodzenie procesu. Jeśli napotkasz tego rodzaju błąd, wróć i sprawdź, czy wyjścia pasują między definicją procesu a blokiem wyjściowym.

Jeśli problem wciąż nie jest jasny, sprawdź sam katalog roboczy, aby zidentyfikować faktycznie utworzone pliki wyjściowe:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

W tym przykładzie podkreśliłoby nam to, że sufiks `_output` jest włączany do nazwy pliku wyjściowego, wbrew naszej definicji `output:`.

#### Napraw kod

Napraw niedopasowanie, czyniąc nazwę pliku wyjściowego spójną:

=== "Po"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Fixed: Match the script output

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
        path "${sample_name}.txt"  // Expects: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Creates: sample3_output.txt
        """
    }
    ```

#### Uruchom pipeline'a

```bash
nextflow run missing_output.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [elated_hamilton] DSL2 - revision: 961938ee2b

    executor >  local (3)
    [16/1c437c] PROCESS_FILES (3) | 3 of 3 ✔
    ```

### 3.2. Brakujące oprogramowanie

Inna klasa błędów występuje z powodu błędów w dostarczaniu oprogramowania. `missing_software.nf` to składniowo prawidłowy workflow, ale zależy od zewnętrznego oprogramowania, aby dostarczyć polecenie `cowpy`, którego używa.

#### Uruchom pipeline'a

```bash
nextflow run missing_software.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Proces nie ma dostępu do polecenia, które określamy. Czasami dzieje się tak, ponieważ skrypt jest obecny w katalogu `bin` workflow'a, ale nie został oznaczony jako wykonywalny. Innym razem dzieje się tak, ponieważ oprogramowanie nie jest zainstalowane w kontenerze lub środowisku, w którym działa workflow.

#### Sprawdź kod

Zwróć uwagę na kod wyjścia `127` - mówi Ci dokładnie, jaki jest problem. Zbadajmy `missing_software.nf`:

```groovy title="missing_software.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

#### Napraw kod

Byliśmy tu trochę nieuczciwi i tak naprawdę nie ma nic złego z kodem. Musimy tylko określić niezbędną konfigurację, aby uruchomić proces w taki sposób, aby miał dostęp do danego polecenia. W tym przypadku proces ma definicję kontenera, więc wszystko, co musimy zrobić, to uruchomić workflow'a z włączonym Dockerem.

#### Uruchom pipeline'a

Skonfigurowaliśmy dla Ciebie profil Docker w `nextflow.config`, więc możesz uruchomić workflow'a za pomocą:

```bash
nextflow run missing_software.nf -profile docker
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software.nf` [awesome_stonebraker] DSL2 - revision: 0296d12839

    executor >  local (3)
    [38/ab20d1] PROCESS_FILES (1) | 3 of 3 ✔
    ```

!!! note

    Aby dowiedzieć się więcej o tym, jak Nextflow używa kontenerów, zobacz [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Zła konfiguracja zasobów

W użyciu produkcyjnym będziesz konfigurować zasoby dla swoich procesów. Na przykład `memory` definiuje maksymalną ilość pamięci dostępnej dla Twojego procesu, a jeśli proces ją przekroczy, Twój scheduler zazwyczaj zabije proces i zwróci kod wyjścia `137`. Nie możemy tego tutaj zademonstrować, ponieważ używamy executora `local`, ale możemy pokazać coś podobnego z `time`.

#### Uruchom pipeline'a

`bad_resources.nf` ma konfigurację procesu z nierealistycznym ograniczeniem czasu wynoszącym 1 milisekundę:

```bash
nextflow run bad_resources.nf -profile docker
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [disturbed_elion] DSL2 - revision: 27d2066e86

    executor >  local (3)
    [c0/ded8e1] PROCESS_FILES (3) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (2)'

    Caused by:
      Process exceeded running time limit (1ms)

    Command executed:

      cowpy sample2 > sample2_output.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/53/f0a4cc56d6b3dc2a6754ff326f1349

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Sprawdź kod

Zbadajmy `bad_resources.nf`:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // ERROR: Unrealistic time limit

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Wiemy, że proces zajmie dłużej niż sekundę (dodaliśmy tam sleep, aby się upewnić), ale proces jest ustawiony na przekroczenie czasu po 1 milisekundzie. Ktoś był trochę nierealistyczny ze swoją konfiguracją!

#### Napraw kod

Zwiększ limit czasu do realistycznej wartości:

=== "Po"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Fixed: Realistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

=== "Przed"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '1 ms'  // ERROR: Unrealistic time limit

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Takes 1 second, but time limit is 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Uruchom pipeline'a

```bash
nextflow run bad_resources.nf -profile docker
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_resources.nf` [friendly_mcclintock] DSL2 - revision: 381567d2c1

    executor >  local (3)
    [c2/9b4c41] PROCESS_FILES (3) | 3 of 3 ✔
    ```

Jeśli upewnisz się, że czytasz swoje komunikaty błędów, niepowodzenia takie jak to nie powinny Cię zbyt długo zastanawiać. Ale upewnij się, że rozumiesz wymagania zasobowe poleceń, które uruchamiasz, abyś mógł odpowiednio skonfigurować swoje dyrektywy zasobów.

### 3.4. Techniki debugowania procesów

Gdy procesy zawodzą lub zachowują się nieoczekiwanie, potrzebujesz systematycznych technik do zbadania, co poszło nie tak. Katalog roboczy zawiera wszystkie informacje potrzebne do debugowania wykonania procesu.

#### Używanie inspekcji katalogu roboczego

Najpotężniejszym narzędziem debugowania dla procesów jest badanie katalogu roboczego. Gdy proces zawodzi, Nextflow tworzy katalog roboczy dla tego konkretnego wykonania procesu zawierający wszystkie pliki potrzebne do zrozumienia, co się stało.

#### Uruchom pipeline'a

Użyjmy przykładu `missing_output.nf` z wcześniejszych rozważań, aby zademonstrować inspekcję katalogu roboczego (wygeneruj ponownie niedopasowanie nazewnictwa wyjścia, jeśli musisz):

```bash
nextflow run missing_output.nf
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_output.nf` [irreverent_payne] DSL2 - revision: 3d5117f7e2

    executor >  local (3)
    [5d/d544a4] PROCESS_FILES (2) | 0 of 3 ✘
    ERROR ~ Error executing process > 'PROCESS_FILES (1)'

    Caused by:
      Missing output file(s) `sample1.txt` expected by process `PROCESS_FILES (1)`

    Command executed:

      echo "Processing sample1" > sample1_output.txt

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/side-quests/debugging/work/1e/2011154d0b0f001cd383d7364b5244

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

#### Sprawdź katalog roboczy

Gdy otrzymasz ten błąd, katalog roboczy zawiera wszystkie informacje debugowania. Znajdź ścieżkę katalogu roboczego z komunikatu błędu i zbadaj jego zawartość:

```bash
# Znajdź katalog roboczy z komunikatu błędu
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Następnie możesz zbadać kluczowe pliki:

##### Sprawdź skrypt polecenia

Plik `.command.sh` pokazuje dokładnie, jakie polecenie zostało wykonane:

```bash
# Wyświetl wykonane polecenie
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

To ujawnia:

- **Podstawienie zmiennych**: Czy zmienne Nextflow zostały prawidłowo rozwinięte
- **Ścieżki plików**: Czy pliki wejściowe zostały prawidłowo zlokalizowane
- **Struktura polecenia**: Czy składnia skryptu jest poprawna

Typowe problemy do wyszukania:

- **Brakujące cudzysłowy**: Zmienne zawierające spacje wymagają odpowiedniego cytowania
- **Złe ścieżki plików**: Pliki wejściowe, które nie istnieją lub są w złych lokalizacjach
- **Nieprawidłowe nazwy zmiennych**: Literówki w odwołaniach do zmiennych
- **Brakująca konfiguracja środowiska**: Polecenia zależne od konkretnych środowisk

##### Sprawdź wyjście błędu

Plik `.command.err` zawiera faktyczne komunikaty błędów:

```bash
# Wyświetl wyjście błędu
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Ten plik pokaże:

- **Kody wyjścia**: 127 (polecenie nie znalezione), 137 (zabity) itp.
- **Błędy uprawnień**: Problemy z dostępem do plików
- **Błędy oprogramowania**: Komunikaty błędów specyficzne dla aplikacji
- **Błędy zasobów**: Przekroczono limit pamięci/czasu

##### Sprawdź standardowe wyjście

Plik `.command.out` pokazuje, co wygenerowało Twoje polecenie:

```bash
# Wyświetl standardowe wyjście
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

To pomaga zweryfikować:

- **Oczekiwane wyjście**: Czy polecenie wygenerowało prawidłowe wyniki
- **Częściowe wykonanie**: Czy polecenie rozpoczęło się, ale zawiodło w połowie
- **Informacje debugowania**: Wszelkie wyjście diagnostyczne z Twojego skryptu

##### Sprawdź kod wyjścia

Plik `.exitcode` zawiera kod wyjścia dla procesu:

```bash
# Wyświetl kod wyjścia
cat work/*/*/.exitcode
```

Typowe kody wyjścia i ich znaczenia:

- **Kod wyjścia 127**: Polecenie nie znalezione - sprawdź instalację oprogramowania
- **Kod wyjścia 137**: Proces zabity - sprawdź limity pamięci/czasu

##### Sprawdź istnienie pliku

Gdy procesy zawodzą z powodu brakujących plików wyjściowych, sprawdź, jakie pliki faktycznie zostały utworzone:

```bash
# Wyświetl wszystkie pliki w katalogu roboczym
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

To pomaga zidentyfikować:

- **Niedopasowania nazw plików**: Pliki wyjściowe o innych nazwach niż oczekiwano
- **Problemy z uprawnieniami**: Pliki, których nie można było utworzyć
- **Problemy ze ścieżkami**: Pliki utworzone w złych katalogach

W naszym wcześniejszym przykładzie potwierdziło nam to, że podczas gdy nasz oczekiwany `sample3.txt` nie był obecny, `sample3_output.txt` był:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Podsumowanie

Debugowanie procesów wymaga badania katalogów roboczych, aby zrozumieć, co poszło nie tak. Kluczowe pliki obejmują `.command.sh` (wykonany skrypt), `.command.err` (komunikaty błędów) i `.command.out` (standardowe wyjście). Kody wyjścia takie jak 127 (polecenie nie znalezione) i 137 (proces zabity) dostarczają natychmiastowych wskazówek diagnostycznych o typie niepowodzenia.

### Co dalej?

Dowiedz się o wbudowanych narzędziach debugowania Nextflow i systematycznych podejściach do rozwiązywania problemów.

---

## 4. Wbudowane narzędzia debugowania i zaawansowane techniki

Nextflow dostarcza kilka potężnych wbudowanych narzędzi do debugowania i analizowania wykonania workflow'a. Te narzędzia pomagają zrozumieć, co poszło nie tak, gdzie poszło nie tak i jak to efektywnie naprawić.

### 4.1. Wyjście procesu w czasie rzeczywistym

Czasami musisz zobaczyć, co dzieje się wewnątrz działających procesów. Możesz włączyć wyjście procesu w czasie rzeczywistym, które pokazuje dokładnie, co robi każde zadanie podczas wykonywania.

#### Uruchom pipeline'a

`bad_channel_shape_viewed.nf` z naszych wcześniejszych przykładów wydrukował zawartość kanału używając `.view()`, ale możemy również użyć dyrektywy `debug`, aby wyświetlić zmienne z wnętrza samego procesu, co demonstrujemy w `bad_channel_shape_viewed_debug.nf`. Uruchom workflow'a:

```bash
nextflow run bad_channel_shape_viewed_debug.nf
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_channel_shape_viewed_debug.nf` [agitated_crick] DSL2 - revision: ea3676d9ec

    executor >  local (3)
    [c6/2dac51] process > PROCESS_FILES (3) [100%] 3 of 3 ✔
    Channel content: [sample1, file1.txt]
    Channel content: [sample2, file2.txt]
    Channel content: [sample3, file3.txt]
    After mapping: sample1
    After mapping: sample2
    After mapping: sample3
    Sample name inside process is sample2

    Sample name inside process is sample1

    Sample name inside process is sample3
    ```

#### Sprawdź kod

Zbadajmy `bad_channel_shape_viewed_debug.nf`, aby zobaczyć, jak działa dyrektywa `debug`:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Enable real-time output

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}
```

Dyrektywa `debug` może być szybkim i wygodnym sposobem na zrozumienie środowiska procesu.

### 4.2. Tryb podglądu

Czasami chcesz wychwycić problemy zanim jakiekolwiek procesy się uruchomią. Nextflow dostarcza flagę dla tego rodzaju proaktywnego debugowania: `-preview`.

#### Uruchom pipeline'a

Tryb podglądu pozwala testować logikę workflow'a bez wykonywania poleceń. Może to być bardzo przydatne do szybkiego sprawdzania struktury Twojego workflow'a i upewnienia się, że procesy są prawidłowo połączone bez uruchamiania jakichkolwiek faktycznych poleceń.

!!! note

    Jeśli naprawiłeś `bad_syntax.nf` wcześniej, przywróć błąd składni, usuwając nawias zamykający po bloku skryptu przed uruchomieniem tego polecenia.

Uruchom to polecenie:

```bash
nextflow run bad_syntax.nf -preview
```

??? failure "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `bad_syntax.nf` [magical_mercator] DSL2 - revision: 550b9a8873

    Error bad_syntax.nf:24:1: Unexpected input: '<EOF>'

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

Tryb podglądu jest szczególnie przydatny do wczesnego wychwytywania błędów składni bez uruchamiania jakichkolwiek procesów. Waliduje strukturę workflow'a i połączenia procesów przed wykonaniem.

### 4.3. Uruchamianie zastępcze do testowania logiki

Czasami błędy są trudne do debugowania, ponieważ polecenia trwają zbyt długo, wymagają specjalnego oprogramowania lub zawodzą z złożonych powodów. Uruchamianie zastępcze pozwala testować logikę workflow'a bez wykonywania faktycznych poleceń.

#### Uruchom pipeline'a

Podczas tworzenia procesu Nextflow możesz użyć dyrektywy `stub`, aby zdefiniować 'zastępcze' polecenia, które generują wyjścia o prawidłowej formie bez uruchamiania prawdziwego polecenia. To podejście jest szczególnie wartościowe, gdy chcesz zweryfikować, że logika Twojego workflow'a jest poprawna przed zajmowaniem się złożonościami faktycznego oprogramowania.

Na przykład, pamiętasz nasz `missing_software.nf` z wcześniejszych rozważań? Ten, gdzie mieliśmy brakujące oprogramowanie, które uniemożliwiało uruchomienie workflow'a, dopóki nie dodaliśmy `-profile docker`? `missing_software_with_stub.nf` to bardzo podobny workflow. Jeśli uruchomimy go w ten sam sposób, wygenerujemy ten sam błąd:

```bash
nextflow run missing_software_with_stub.nf
```

??? failure "Wyjście polecenia"

    ```console hl_lines="12 18"
    ERROR ~ Error executing process > 'PROCESS_FILES (3)'

    Caused by:
      Process `PROCESS_FILES (3)` terminated with an error exit status (127)


    Command executed:

      cowpy sample3 > sample3_output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/side-quests/debugging/work/82/42a5bfb60c9c6ee63ebdbc2d51aa6e

    Tip: you can try to figure out what's wrong by changing to the process work directory and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Jednak ten workflow nie wygeneruje błędów, jeśli uruchomimy go z `-stub-run`, nawet bez profilu `docker`:

```bash
nextflow run missing_software_with_stub.nf -stub-run
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `missing_software_with_stub.nf` [astonishing_shockley] DSL2 - revision: f1f4f05d7d

    executor >  local (3)
    [b5/2517a3] PROCESS_FILES (3) | 3 of 3 ✔
    ```

#### Sprawdź kod

Zbadajmy `missing_software_with_stub.nf`:

```groovy title="missing_software.nf (with stub)" hl_lines="16-19" linenums="3"
process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    cowpy ${sample_name} > ${sample_name}_output.txt
    """

    stub:
    """
    touch ${sample_name}_output.txt
    """
}
```

W stosunku do `missing_software.nf`, ten proces ma dyrektywę `stub:` określającą polecenie do użycia zamiast tego określonego w `script:`, w przypadku gdy Nextflow jest uruchamiany w trybie zastępczym.

Polecenie `touch`, którego tu używamy, nie zależy od żadnego oprogramowania ani odpowiednich wejść i będzie działać we wszystkich sytuacjach, pozwalając nam debugować logikę workflow'a bez martwienia się o wewnętrzne elementy procesu.

**Uruchamianie zastępcze pomaga debugować:**

- Strukturę kanałów i przepływ danych
- Połączenia i zależności procesów
- Propagację parametrów
- Logikę workflow'a bez zależności od oprogramowania

### 4.4. Systematyczne podejście do debugowania

Teraz, gdy nauczyłeś się indywidualnych technik debugowania - od plików śledzenia i katalogów roboczych po tryb podglądu, uruchamianie zastępcze i monitorowanie zasobów - połączmy je w systematyczną metodologię. Posiadanie ustrukturyzowanego podejścia zapobiega przytłoczeniu przez złożone błędy i zapewnia, że nie przegapisz ważnych wskazówek.

Ta metodologia łączy wszystkie narzędzia, które omówiliśmy, w efektywny workflow:

**Czterofazowa metoda debugowania:**

**Faza 1: Rozwiązywanie błędów składni (5 minut)**

1. Sprawdź czerwone podkreślenia w VSCode lub swoim IDE
2. Uruchom `nextflow run workflow.nf -preview`, aby zidentyfikować problemy składniowe
3. Napraw wszystkie błędy składni (brakujące nawiasy, końcowe przecinki itp.)
4. Upewnij się, że workflow parsuje się pomyślnie przed kontynuowaniem

**Faza 2: Szybka ocena (5 minut)**

1. Uważnie przeczytaj komunikaty błędów w czasie wykonania
2. Sprawdź, czy to błąd wykonania, logiki czy zasobów
3. Użyj trybu podglądu do testowania podstawowej logiki workflow'a

**Faza 3: Szczegółowe badanie (15-30 minut)**

1. Znajdź katalog roboczy nieudanego zadania
2. Zbadaj pliki dziennika
3. Dodaj operatory `.view()` do inspekcji kanałów
4. Użyj `-stub-run` do testowania logiki workflow'a bez wykonania

**Faza 4: Napraw i zwaliduj (15 minut)**

1. Wprowadź minimalne ukierunkowane poprawki
2. Testuj z wznowieniem: `nextflow run workflow.nf -resume`
3. Zweryfikuj kompletne wykonanie workflow'a

!!! tip "Używanie wznowienia do efektywnego debugowania"

    Gdy już zidentyfikujesz problem, potrzebujesz efektywnego sposobu testowania swoich poprawek bez marnowania czasu na ponowne uruchamianie udanych części Twojego workflow'a. Funkcjonalność `-resume` Nextflow jest nieoceniona podczas debugowania.

    Napotkałeś `-resume`, jeśli pracowałeś przez [Hello Nextflow](../hello_nextflow/), i ważne jest, abyś dobrze z niego korzystał podczas debugowania, aby zaoszczędzić sobie czekania, podczas gdy procesy przed Twoim problemowym procesem działają.

    **Strategia debugowania z wznowieniem:**

    1. Uruchom workflow do niepowodzenia
    2. Zbadaj katalog roboczy dla nieudanego zadania
    3. Napraw konkretny problem
    4. Wznów, aby przetestować tylko poprawkę
    5. Powtarzaj, aż workflow się zakończy

#### Profil konfiguracji debugowania

Aby uczynić to systematyczne podejście jeszcze bardziej efektywnym, możesz utworzyć dedykowaną konfigurację debugowania, która automatycznie włącza wszystkie potrzebne narzędzia:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Conservative resources for debugging
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Następnie możesz uruchomić pipeline'a z włączonym tym profilem:

```bash
nextflow run workflow.nf -profile debug
```

Ten profil włącza wyjście w czasie rzeczywistym, zachowuje katalogi robocze i ogranicza paralelizację dla łatwiejszego debugowania.

### 4.5. Praktyczne ćwiczenie debugowania

Teraz czas zastosować systematyczne podejście do debugowania w praktyce. Workflow `buggy_workflow.nf` zawiera kilka typowych błędów reprezentujących typy problemów, które napotkasz w rzeczywistym rozwoju.

!!! exercise

    Użyj systematycznego podejścia do debugowania, aby zidentyfikować i naprawić wszystkie błędy w `buggy_workflow.nf`. Ten workflow próbuje przetworzyć dane próbek z pliku CSV, ale zawiera wiele celowych błędów reprezentujących typowe scenariusze debugowania.

    Zacznij od uruchomienia workflow'a, aby zobaczyć pierwszy błąd:

    ```bash
    nextflow run buggy_workflow.nf
    ```

    ??? failure "Wyjście polecenia"

        ```console
        N E X T F L O W   ~  version 25.10.2

        Launching `buggy_workflow.nf` [wise_ramanujan] DSL2 - revision: d51a8e83fd

        ERROR ~ Range [11, 12) out of bounds for length 11

         -- Check '.nextflow.log' file for details
        ```

        Ten zagadkowy błąd wskazuje problem z parsowaniem wokół linii 11-12 w bloku `params{}`. Parser v2 wychwytuje problemy strukturalne wcześnie.

    Zastosuj czterofazową metodę debugowania, której się nauczyłeś:

    **Faza 1: Rozwiązywanie błędów składni**
    - Sprawdź czerwone podkreślenia w VSCode lub swoim IDE
    - Uruchom `nextflow run workflow.nf -preview`, aby zidentyfikować problemy składniowe
    - Napraw wszystkie błędy składni (brakujące nawiasy, końcowe przecinki itp.)
    - Upewnij się, że workflow parsuje się pomyślnie przed kontynuowaniem

    **Faza 2: Szybka ocena**
    - Uważnie przeczytaj komunikaty błędów w czasie wykonania
    - Zidentyfikuj, czy błędy są związane z wykonaniem, logiką czy zasobami
    - Użyj trybu `-preview` do testowania podstawowej logiki workflow'a

    **Faza 3: Szczegółowe badanie**
    - Zbadaj katalogi robocze dla nieudanych zadań
    - Dodaj operatory `.view()` do inspekcji kanałów
    - Sprawdź pliki dziennika w katalogach roboczych
    - Użyj `-stub-run` do testowania logiki workflow'a bez wykonania

    **Faza 4: Napraw i zwaliduj**
    - Wprowadź ukierunkowane poprawki
    - Użyj `-resume` do efektywnego testowania poprawek
    - Zweryfikuj kompletne wykonanie workflow'a

    **Narzędzia debugowania do Twojej dyspozycji:**
    ```bash
    # Tryb podglądu do sprawdzania składni
    nextflow run buggy_workflow.nf -preview

    # Profil debugowania dla szczegółowego wyjścia
    nextflow run buggy_workflow.nf -profile debug

    # Uruchamianie zastępcze do testowania logiki
    nextflow run buggy_workflow.nf -stub-run

    # Wznów po poprawkach
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution
        `buggy_workflow.nf` zawiera 9 lub 10 odrębnych błędów (w zależności od sposobu liczenia) obejmujących wszystkie główne kategorie debugowania. Oto systematyczne zestawienie każdego błędu i sposobu jego naprawienia

        Zacznijmy od tych błędów składni:

        **Błąd 1: Błąd składni - końcowy przecinek**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // ERROR: Trailing comma
        ```
        **Poprawka:** Usuń końcowy przecinek
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Błąd 2: Błąd składni - brakujący nawias zamykający**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // ERROR: Missing closing brace for processFiles process
        ```
        **Poprawka:** Dodaj brakujący nawias zamykający
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Add missing closing brace
        ```

        **Błąd 3: Błąd nazwy zmiennej**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // ERROR: should be sample_id
        cat ${input_file} > ${sample}_result.txt  // ERROR: should be sample_id
        ```
        **Poprawka:** Użyj prawidłowej nazwy zmiennej wejściowej
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Błąd 4: Błąd niezdefiniowanej zmiennej**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // ERROR: sample_ids undefined
        ```
        **Poprawka:** Użyj prawidłowego kanału i wyodrębnij identyfikatory próbek
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        W tym momencie workflow będzie działać, ale wciąż będziemy otrzymywać błędy (np. `Path value cannot be null` w `processFiles`), spowodowane złą strukturą kanału.

        **Błąd 5: Błąd struktury kanału - złe wyjście mapy**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // ERROR: processFiles expects tuple
        ```
        **Poprawka:** Zwróć strukturę krotki, której oczekuje processFiles
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Ale to zepsuje nasz sposób uruchamiania `heavyProcess()` powyżej, więc będziemy musieli użyć mapy, aby przekazać tylko identyfikatory próbek do tego procesu:

        **Błąd 6: Zła struktura kanału dla heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // ERROR: input_ch now has 2 elements per emission- heavyProcess only needs 1 (the first)
        ```
        **Poprawka:** Użyj prawidłowego kanału i wyodrębnij identyfikatory próbek
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Teraz dostajemy się trochę dalej, ale otrzymujemy błąd o `No such variable: i`, ponieważ nie poprzedziliśmy zmiennej Bash odwrotnym ukośnikiem.

        **Błąd 7: Błąd poprzedzania zmiennej Bash odwrotnym ukośnikiem**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // ERROR: $i not escaped
        ```
        **Poprawka:** Poprzedź zmienną bash odwrotnym ukośnikiem
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Teraz otrzymujemy `Process exceeded running time limit (1ms)`, więc naprawiamy limit czasu wykonania dla odpowiedniego procesu:

        **Błąd 8: Błąd konfiguracji zasobów**
        ```groovy linenums="36"
        time '1 ms'  // ERROR: Unrealistic time limit
        ```
        **Poprawka:** Zwiększ do realistycznego limitu czasu
        ```groovy linenums="36"
        time '100 s'
        ```

        Następnie mamy błąd `Missing output file(s)` do rozwiązania:

        **Błąd 9: Niedopasowanie nazwy pliku wyjściowego**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // ERROR: Wrong filename, should match output declaration
        ```
        **Poprawka:** Dopasuj deklarację wyjścia
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Pierwsze dwa procesy działały, ale nie trzeci.

        **Błąd 10: Niedopasowanie nazwy pliku wyjściowego**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Error: attempting to take input from the pwd rather than a process
        handleFiles(file_ch)
        ```
        **Poprawka:** Weź wyjście z poprzedniego procesu
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Dzięki temu cały workflow powinien działać.

        **Kompletny poprawiony workflow:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Błędny workflow do ćwiczeń debugowania
        * Ten workflow zawiera kilka celowych błędów do celów edukacyjnych
        */

        params{
            // Parameters with missing validation
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Proces z niedopasowaniem wejścia/wyjścia
        */
        process processFiles {
            publishDir "${params.output}/processed", mode: 'copy'

            input:
                tuple val(sample_id), path(input_file)

            output:
                path "${sample_id}_result.txt"

            script:
            """
            echo "Processing: ${sample_id}"
            cat ${input_file} > ${sample_id}_result.txt
            """
        }

        /*
        * Proces z problemami zasobów
        */
        process heavyProcess {
            publishDir "${params.output}/heavy", mode: 'copy'

            time '100 s'

            input:
                val sample_id

            output:
                path "${sample_id}_heavy.txt"

            script:
            """
            # Symuluj ciężkie obliczenia
            for i in {1..1000000}; do
                echo "Heavy computation \$i for ${sample_id}"
            done > ${sample_id}_heavy.txt
            """
        }

        /*
        * Proces z problemami obsługi plików
        */
        process handleFiles {
            publishDir "${params.output}/files", mode: 'copy'

            input:
                path input_file

            output:
                path "processed_${input_file}"

            script:
            """
            if [ -f "${input_file}" ]; then
                cp ${input_file} processed_${input_file}
            fi
            """
        }

        /*
        * Główny workflow z problemami kanałów
        */
        workflow {

            // Channel with incorrect usage
            input_ch = channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row -> [row.sample_id, file(row.fastq_path)] }

            processed_ch = processFiles(input_ch)

            heavy_ch = heavyProcess(input_ch.map{it[0]})

            handleFiles(heavyProcess.out)
        }
        ```

**Omówione kategorie błędów:**

- **Błędy składni**: Brakujące nawiasy, końcowe przecinki, niezdefiniowane zmienne
- **Błędy struktury kanałów**: Złe kształty danych, niezdefiniowane kanały
- **Błędy procesów**: Niedopasowania plików wyjściowych, poprzedzanie zmiennych odwrotnym ukośnikiem
- **Błędy zasobów**: Nierealistyczne limity czasu

**Kluczowe lekcje debugowania:**

1. **Uważnie czytaj komunikaty błędów** - często wskazują bezpośrednio na problem
2. **Używaj systematycznych podejść** - naprawiaj jeden błąd na raz i testuj z `-resume`
3. **Rozumiej przepływ danych** - błędy struktury kanałów są często najbardziej subtelne
4. **Sprawdzaj katalogi robocze** - gdy procesy zawodzą, dzienniki mówią dokładnie, co poszło nie tak

---

## Podsumowanie

W tym zadaniu pobocznym nauczyłeś się zestawu systematycznych technik debugowania workflow'ów Nextflow.
Zastosowanie tych technik w swojej pracy pozwoli Ci spędzać mniej czasu na walce z komputerem, szybciej rozwiązywać problemy i chronić się przed przyszłymi problemami.

### Kluczowe wzorce

**1. Jak identyfikować i naprawiać błędy składni**:

- Interpretowanie komunikatów błędów Nextflow i lokalizowanie problemów
- Typowe błędy składni: brakujące nawiasy, nieprawidłowe słowa kluczowe, niezdefiniowane zmienne
- Rozróżnianie między zmiennymi Nextflow (Groovy) a Bash
- Używanie funkcji rozszerzenia VS Code do wczesnego wykrywania błędów

```groovy
// Missing brace - look for red underlines in IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- missing!

// Wrong keyword
inputs:  // Should be 'input:'

// Undefined variable - escape with backslash for Bash variables
echo "${undefined_var}"      // Nextflow variable (error if not defined)
echo "\${bash_var}"          // Bash variable (escaped)
```

**2. Jak debugować problemy ze strukturą kanałów**:

- Rozumienie kardynalności kanałów i problemów z wyczerpaniem
- Debugowanie niedopasowań struktury zawartości kanału
- Używanie operatorów `.view()` do inspekcji kanałów
- Rozpoznawanie wzorców błędów, takich jak nawiasy kwadratowe w wyjściu

```groovy
// Inspect channel content
my_channel.view { "Content: $it" }

// Convert queue to value channel (prevents exhaustion)
reference_ch = channel.value('ref.fa')
// or
reference_ch = channel.of('ref.fa').first()
```

**3. Jak rozwiązywać problemy z wykonaniem procesów**:

- Diagnozowanie błędów brakujących plików wyjściowych
- Rozumienie kodów wyjścia (127 dla brakującego oprogramowania, 137 dla problemów z pamięcią)
- Badanie katalogów roboczych i plików poleceń
- Odpowiednie konfigurowanie zasobów

```bash
# Sprawdź, co faktycznie zostało wykonane
cat work/ab/cdef12/.command.sh

# Sprawdź wyjście błędu
cat work/ab/cdef12/.command.err

# Kod wyjścia 127 = polecenie nie znalezione
# Kod wyjścia 137 = zabity (limit pamięci/czasu)
```

**4. Jak używać wbudowanych narzędzi debugowania Nextflow**:

- Wykorzystywanie trybu podglądu i debugowania w czasie rzeczywistym
- Implementowanie uruchamiania zastępczego do testowania logiki
- Stosowanie wznowienia dla efektywnych cykli debugowania
- Przestrzeganie czterofazowej systematycznej metodologii debugowania

!!! tip "Szybka referencja debugowania"

    **Błędy składni?** → Sprawdź ostrzeżenia VSCode, uruchom `nextflow run workflow.nf -preview`

    **Problemy z kanałami?** → Użyj `.view()` do inspekcji zawartości: `my_channel.view()`

    **Niepowodzenia procesów?** → Sprawdź pliki katalogu roboczego:

    - `.command.sh` - wykonany skrypt
    - `.command.err` - komunikaty błędów
    - `.exitcode` - status wyjścia (127 = polecenie nie znalezione, 137 = zabity)

    **Tajemnicze zachowanie?** → Uruchom z `-stub-run`, aby przetestować logikę workflow'a

    **Wprowadzono poprawki?** → Użyj `-resume`, aby zaoszczędzić czas testowania: `nextflow run workflow.nf -resume`

---

### Dodatkowe zasoby

- [Przewodnik rozwiązywania problemów Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html): Oficjalna dokumentacja rozwiązywania problemów
- [Zrozumienie kanałów Nextflow](https://www.nextflow.io/docs/latest/channel.html): Głębokie zanurzenie w typy kanałów i zachowanie
- [Referencja dyrektyw procesów](https://www.nextflow.io/docs/latest/process.html#directives): Wszystkie dostępne opcje konfiguracji procesów
- [nf-test](https://www.nf-test.com/): Framework testowy dla pipeline'ów Nextflow
- [Społeczność Nextflow Slack](https://www.nextflow.io/slack-invite.html): Uzyskaj pomoc od społeczności

Dla workflow'ów produkcyjnych rozważ:

- Skonfigurowanie [Seqera Platform](https://seqera.io/platform/) do monitorowania i debugowania na dużą skalę
- Używanie [kontenerów Wave](https://seqera.io/wave/) dla odtwarzalnych środowisk oprogramowania

**Pamiętaj:** Efektywne debugowanie to umiejętność, która poprawia się z praktyką. Systematyczna metodologia i kompleksowy zestaw narzędzi, które tutaj nabyłeś, będą Ci dobrze służyć przez całą Twoją drogę rozwoju Nextflow.

---

## Co dalej?

Wróć do [menu zadań pobocznych](./index.md) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
