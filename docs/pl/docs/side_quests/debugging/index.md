# Debugowanie Workflow'ów

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Debugowanie to kluczowa umiejętność, która może zaoszczędzić Ci wiele godzin frustracji i pomóc stać się bardziej efektywnym programistą Nextflow. W trakcie swojej kariery, szczególnie na początku, będziesz napotykać błędy podczas budowania i utrzymywania workflow'ów. Nauka systematycznych podejść do debugowania pomoże Ci szybko identyfikować i rozwiązywać problemy.

### Cele szkolenia

W tym side queście poznamy **systematyczne techniki debugowania** workflow'ów Nextflow:

- **Debugowanie błędów składni**: Efektywne korzystanie z funkcji IDE i komunikatów o błędach Nextflow
- **Debugowanie kanałów**: Diagnozowanie problemów z przepływem danych i strukturą kanałów
- **Debugowanie procesów**: Badanie błędów wykonania i problemów z zasobami
- **Wbudowane narzędzia debugowania**: Korzystanie z trybu podglądu Nextflow, stub running i katalogów roboczych
- **Systematyczne podejścia**: Czterofazowa metodologia efektywnego debugowania

Po zakończeniu będziesz dysponować solidną metodologią debugowania, która przekształca frustrujące komunikaty o błędach w czytelne wskazówki prowadzące do rozwiązań.

### Wymagania wstępne

Przed podjęciem tego side questu powinieneś/powinnaś:

- Ukończyć samouczek [Hello Nextflow](../hello_nextflow/README.md) lub równoważny kurs dla początkujących.
- Swobodnie posługiwać się podstawowymi konceptami i mechanizmami Nextflow (procesy, kanały, operatory)

**Opcjonalnie:** Zalecamy wcześniejsze ukończenie side questu [IDE Features for Nextflow Development](../dev_environment/).
Obejmuje on kompleksowe omówienie funkcji IDE wspierających debugowanie (podświetlanie składni, wykrywanie błędów itp.), z których będziemy tu intensywnie korzystać.

---

## 0. Pierwsze kroki

#### Otwórz środowisko szkoleniowe

Jeśli jeszcze tego nie zrobiłeś/zrobiłaś, otwórz środowisko szkoleniowe zgodnie z opisem w sekcji [Konfiguracja środowiska](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Przejdź do katalogu projektu

Przejdźmy do katalogu, w którym znajdują się pliki tego samouczka.

```bash
cd side-quests/debugging
```

Możesz ustawić VSCode tak, aby skupiał się na tym katalogu:

```bash
code .
```

#### Przejrzyj materiały

Znajdziesz tu zestaw przykładowych workflow'ów z różnymi rodzajami błędów, których użyjemy do ćwiczeń:

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

Pliki te reprezentują typowe scenariusze debugowania, z którymi spotkasz się w rzeczywistym środowisku programistycznym.

#### Zapoznaj się z zadaniem

Twoim wyzwaniem jest uruchomienie każdego workflow'u, zidentyfikowanie błędów i ich naprawienie.

Dla każdego błędnego workflow'u:

1. **Uruchom workflow** i zaobserwuj błąd
2. **Przeanalizuj komunikat o błędzie**: co mówi Ci Nextflow?
3. **Zlokalizuj problem** w kodzie, korzystając z dostarczonych wskazówek
4. **Napraw błąd** i zweryfikuj, że rozwiązanie działa
5. **Zresetuj plik** przed przejściem do następnej sekcji (użyj `git checkout <filename>`)

Ćwiczenia przechodzą od prostych błędów składni do bardziej subtelnych problemów w czasie wykonania.
Rozwiązania są omawiane na bieżąco, ale spróbuj rozwiązać każde z nich samodzielnie przed przeczytaniem dalej.

#### Lista kontrolna gotowości

Myślisz, że jesteś gotowy/gotowa?

- [ ] Rozumiem cel tego kursu i jego wymagania wstępne
- [ ] Moje środowisko jest uruchomione i działa
- [ ] Ustawiłem/ustawiłam odpowiedni katalog roboczy
- [ ] Rozumiem zadanie

Jeśli możesz zaznaczyć wszystkie pola, możesz zaczynać.

---

## 1. Błędy składni

Błędy składni to najczęstszy rodzaj błędów, z jakimi spotkasz się podczas pisania kodu Nextflow. Występują, gdy kod nie jest zgodny z oczekiwanymi regułami składni DSL Nextflow. Uniemożliwiają one uruchomienie workflow'u, dlatego ważne jest, aby nauczyć się je szybko identyfikować i naprawiać.

### 1.1. Brakujące nawiasy klamrowe

Jednym z najczęstszych błędów składni, a zarazem jednym z trudniejszych do debugowania, są **brakujące lub niedopasowane nawiasy klamrowe**.

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

**Kluczowe elementy komunikatów o błędach składni:**

- **Plik i lokalizacja**: Wskazuje, który plik i wiersz/kolumna zawierają błąd (`bad_syntax.nf:24:1`)
- **Opis błędu**: Wyjaśnia, co parser napotkał, czego się nie spodziewał (`Unexpected input: '<EOF>'`)
- **Wskaźnik EOF**: Komunikat `<EOF>` (End Of File) oznacza, że parser dotarł do końca pliku, wciąż oczekując dalszej zawartości — klasyczny znak niezamkniętych nawiasów klamrowych

#### Sprawdź kod

Przyjrzyjmy się teraz `bad_syntax.nf`, aby zrozumieć przyczynę błędu:

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

Na potrzeby tego przykładu zostawiliśmy komentarz wskazujący miejsce błędu. Rozszerzenie Nextflow dla VSCode powinno również dawać Ci wskazówki dotyczące problemu — zaznacza niedopasowany nawias na czerwono i podświetla przedwczesny koniec pliku:

![Błędna składnia](img/bad_syntax.png)

**Strategia debugowania błędów nawiasów:**

1. Użyj dopasowywania nawiasów w VS Code (umieść kursor obok nawiasu)
2. Sprawdź panel Problems pod kątem komunikatów związanych z nawiasami
3. Upewnij się, że każdy otwierający `{` ma odpowiadający mu zamykający `}`

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

Uruchom workflow ponownie, aby potwierdzić, że działa:

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

Innym częstym błędem składni jest **nieprawidłowa definicja procesu**. Może się to zdarzyć, gdy zapomnisz zdefiniować wymagane bloki lub użyjesz nieprawidłowych dyrektyw w definicji procesu.

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

Komunikat o błędzie wskazuje na „Invalid process definition" i pokazuje kontekst wokół problemu. Patrząc na wiersze 3–7, widzimy `inputs:` w wierszu 4 — to właśnie jest problem. Przyjrzyjmy się `invalid_process.nf`:

```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    inputs:  // BŁĄD: Powinno być 'input', nie 'inputs'
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

Patrząc na wiersz 4 w kontekście błędu, możemy dostrzec problem: używamy `inputs` zamiast poprawnej dyrektywy `input`. Rozszerzenie Nextflow dla VSCode również to oznaczy:

![Komunikat o nieprawidłowym procesie](img/invalid_process_message.png)

#### Napraw kod

Zastąp nieprawidłowe słowo kluczowe poprawnym, korzystając z [dokumentacji](https://www.nextflow.io/docs/latest/process.html#):

=== "Po"

    ```groovy title="invalid_process.nf" hl_lines="4" linenums="1"
    #!/usr/bin/env nextflow

    process PROCESS_FILES {
        input:  // Naprawiono: zmieniono 'inputs' na 'input'
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
        inputs:  // BŁĄD: Powinno być 'input', nie 'inputs'
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

Uruchom workflow ponownie, aby potwierdzić, że działa:

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

Nazwy zmiennych używane w blokach skryptu muszą być prawidłowe — muszą pochodzić z wejść lub z kodu Groovy wstawionego przed skryptem. Jednak gdy na początku tworzenia pipeline'u zmagasz się ze złożonością, łatwo popełnić błędy w nazewnictwie zmiennych, a Nextflow szybko Cię o tym poinformuje.

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

Błąd jest wykrywany w czasie kompilacji i wskazuje bezpośrednio na niezdefiniowaną zmienną w wierszu 17, a daszek precyzyjnie pokazuje miejsce problemu.

#### Sprawdź kod

Przyjrzyjmy się `no_such_var.nf`:

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
    echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // BŁĄD: undefined_var nie jest zdefiniowana
    """
}

workflow {
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    PROCESS_FILES(input_ch)
}
```

Komunikat o błędzie wskazuje, że zmienna nie jest rozpoznawana w szablonie skryptu — i rzeczywiście, możesz zobaczyć `${undefined_var}` użyte w bloku skryptu, ale niezdefiniowane nigdzie indziej.

#### Napraw kod

Jeśli otrzymasz błąd „No such variable", możesz go naprawić, definiując zmienną (poprawiając nazwy zmiennych wejściowych lub edytując kod Groovy przed skryptem) albo usuwając ją z bloku skryptu, jeśli nie jest potrzebna:

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
        """  // Usunięto wiersz z undefined_var
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
        echo "Using undefined variable: ${undefined_var}" >> ${output_prefix}.txt  // BŁĄD: undefined_var nie jest zdefiniowana
        """
    }

    workflow {
        input_ch = channel.of('sample1', 'sample2', 'sample3')
        PROCESS_FILES(input_ch)
    }
    ```

#### Uruchom pipeline

Uruchom workflow ponownie, aby potwierdzić, że działa:

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

Na początku pracy z Nextflow trudno jest zrozumieć różnicę między zmiennymi Nextflow (Groovy) a zmiennymi Bash. Może to generować inną postać błędu złej zmiennej, który pojawia się przy próbie użycia zmiennych w zawartości Bash bloku skryptu.

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

Błąd wskazuje na wiersz 13, gdzie używane jest `${prefix}`. Przyjrzyjmy się `bad_bash_var.nf`, aby zobaczyć przyczynę problemu:

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

W tym przykładzie definiujemy zmienną `prefix` w Bash, ale w procesie Nextflow składnia `$`, której użyliśmy do jej odwołania (`${prefix}`), jest interpretowana jako zmienna Groovy, a nie Bash. Zmienna nie istnieje w kontekście Groovy, więc otrzymujemy błąd „no such variable".

#### Napraw kod

Jeśli chcesz użyć zmiennej Bash, musisz poprzedzić znak dolara ukośnikiem odwrotnym:

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
        echo "Processing ${sample_name}" > \${prefix}.txt  # Naprawiono: poprzedzono znak dolara ukośnikiem
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

To mówi Nextflow, aby interpretował tę zmienną jako zmienną Bash.

#### Uruchom pipeline

Uruchom workflow ponownie, aby potwierdzić, że działa:

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

!!! tip "Zmienne Groovy a zmienne Bash"

    W przypadku prostych operacji na zmiennych, takich jak konkatenacja ciągów znaków czy operacje na prefiksach/sufiksach, zazwyczaj bardziej czytelne jest używanie zmiennych Groovy w sekcji skryptu zamiast zmiennych Bash w bloku skryptu:

    ```groovy linenums="1"
    script:
    def output_prefix = "${sample_name}_processed"
    def output_file = "${output_prefix}.txt"
    """
    echo "Processing ${sample_name}" > ${output_file}
    """
    ```

    Takie podejście eliminuje potrzebę poprzedzania znaków dolara ukośnikiem i sprawia, że kod jest łatwiejszy do czytania i utrzymania.

### 1.5. Instrukcje poza blokiem workflow

Rozszerzenie Nextflow dla VSCode podświetla problemy ze strukturą kodu, które spowodują błędy. Typowym przykładem jest definiowanie kanałów poza blokiem `workflow {}` — jest to teraz egzekwowane jako błąd składni.

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

Komunikat o błędzie jasno wskazuje problem: instrukcje (takie jak definicje kanałów) nie mogą być mieszane z deklaracjami skryptu poza blokiem workflow lub process.

#### Sprawdź kod

Przyjrzyjmy się `badpractice_syntax.nf`, aby zobaczyć przyczynę błędu:

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

![Nieśmiertelny błąd składni](img/nonlethal.png)

#### Napraw kod

Przenieś definicję kanału do bloku workflow:

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
        input_ch = channel.of('sample1', 'sample2', 'sample3')  // Przeniesiono do bloku workflow
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

Definiuj kanały wejściowe wewnątrz bloku workflow i ogólnie stosuj się do innych zaleceń rozszerzenia.

### Podsumowanie

Możesz systematycznie identyfikować i naprawiać błędy składni, korzystając z komunikatów o błędach Nextflow i wizualnych wskaźników IDE. Typowe błędy składni obejmują brakujące nawiasy klamrowe, nieprawidłowe słowa kluczowe procesów, niezdefiniowane zmienne oraz nieprawidłowe użycie zmiennych Bash kontra Nextflow. Rozszerzenie VSCode pomaga wykryć wiele z nich przed uruchomieniem. Mając te umiejętności debugowania składni w swoim arsenale, będziesz w stanie szybko rozwiązywać najczęstsze błędy składni Nextflow i przejść do bardziej złożonych problemów w czasie wykonania.

### Co dalej?

Naucz się debugować bardziej złożone błędy struktury kanałów, które występują nawet gdy składnia jest poprawna.

---

## 2. Błędy struktury kanałów

Błędy struktury kanałów są bardziej subtelne niż błędy składni, ponieważ kod jest składniowo poprawny, ale kształty danych nie pasują do tego, czego oczekują procesy. Nextflow spróbuje uruchomić pipeline, ale może stwierdzić, że liczba wejść nie zgadza się z oczekiwaną i zakończyć działanie błędem. Tego rodzaju błędy zazwyczaj pojawiają się dopiero w czasie wykonania i wymagają zrozumienia danych przepływających przez workflow.

!!! tip "Debugowanie kanałów za pomocą `.view()`"

    W tej sekcji pamiętaj, że możesz używać operatora `.view()` do inspekcji zawartości kanału w dowolnym miejscu workflow'u. To jedno z najpotężniejszych narzędzi debugowania do rozumienia problemów ze strukturą kanałów. Szczegółowo omówimy tę technikę w sekcji 2.4, ale możesz jej używać już podczas pracy z przykładami.

    ```groovy
    my_channel.view()  // Pokazuje, co przepływa przez kanał
    ```

### 2.1. Nieprawidłowa liczba kanałów wejściowych

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

Komunikat o błędzie jasno stwierdza, że wywołanie oczekiwało 1 argumentu, ale otrzymało 2, i wskazuje na wiersz 23. Przyjrzyjmy się `bad_number_inputs.nf`:

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

Powinieneś/powinnaś zobaczyć niedopasowane wywołanie `PROCESS_FILES`, które dostarcza wiele kanałów wejściowych, podczas gdy proces definiuje tylko jeden. Rozszerzenie VSCode również podkreśli wywołanie procesu na czerwono i wyświetli komunikat diagnostyczny po najechaniu myszą:

![Komunikat o nieprawidłowej liczbie argumentów](img/incorrect_num_args.png)

#### Napraw kod

W tym konkretnym przykładzie proces oczekuje jednego kanału i nie potrzebuje drugiego, więc możemy to naprawić, przekazując tylko kanał `samples_ch`:

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

Częściej niż w tym przykładzie możesz dodać dodatkowe wejścia do procesu i zapomnieć odpowiednio zaktualizować wywołanie workflow'u, co może prowadzić do tego rodzaju błędu. Na szczęście jest to jeden z łatwiejszych do zrozumienia i naprawienia błędów, ponieważ komunikat o błędzie jest dość jasny co do niezgodności.

### 2.2. Wyczerpanie kanału (proces uruchamia się rzadziej niż oczekiwano)

Niektóre błędy struktury kanałów są znacznie bardziej subtelne i nie generują żadnych błędów. Prawdopodobnie najczęstszy z nich odzwierciedla wyzwanie, z którym borykają się nowi użytkownicy Nextflow: zrozumienie, że kanały kolejki mogą się wyczerpać i skończyć im się elementy, co powoduje przedwczesne zakończenie workflow'u.

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

Workflow kończy się bez błędu, ale przetwarza tylko jedną próbkę!

#### Sprawdź kod

Przyjrzyjmy się `exhausted.nf`, aby sprawdzić, czy to prawidłowe zachowanie:

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

Proces uruchamia się tylko raz zamiast trzy razy, ponieważ `reference_ch` jest kanałem kolejki, który wyczerpuje się po pierwszym wykonaniu procesu. Gdy jeden kanał się wyczerpie, cały proces zatrzymuje się, nawet jeśli inne kanały wciąż mają elementy.

Jest to typowy wzorzec, gdy masz jeden plik referencyjny, który musi być wielokrotnie używany dla wielu próbek. Rozwiązaniem jest przekształcenie kanału referencyjnego w kanał wartości, który może być używany wielokrotnie.

#### Napraw kod

Istnieje kilka sposobów rozwiązania tego problemu, w zależności od liczby plików, których dotyczy.

**Opcja 1**: Masz jeden plik referencyjny, który wielokrotnie używasz. Możesz po prostu utworzyć kanał wartości, który można używać wielokrotnie. Są trzy sposoby, aby to zrobić:

**1a** Użyj `channel.value()`:

```groovy title="exhausted.nf (fixed - Option 1a)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.value('baseline_reference')  // Kanał wartości można wielokrotnie używać
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1b** Użyj operatora `first()` ([dokumentacja](https://www.nextflow.io/docs/latest/reference/operator.html#first)):

```groovy title="exhausted.nf (fixed - Option 1b)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').first()  // Konwertuj na kanał wartości
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**1c.** Użyj operatora `collect()` ([dokumentacja](https://www.nextflow.io/docs/latest/reference/operator.html#collect)):

```groovy title="exhausted.nf (fixed - Option 1c)" hl_lines="2" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference').collect()  // Konwertuj na kanał wartości
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    PROCESS_FILES(reference_ch, input_ch)
}
```

**Opcja 2**: W bardziej złożonych scenariuszach, na przykład gdy masz wiele plików referencyjnych dla wszystkich próbek w kanale próbek, możesz użyć operatora `combine`, aby utworzyć nowy kanał łączący oba kanały w krotki:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="4" linenums="21"
workflow {
    reference_ch = channel.of('baseline_reference','other_reference')
    input_ch = channel.of('sample1', 'sample2', 'sample3')
    combined_ch = reference_ch.combine(input_ch)  // Tworzy iloczyn kartezjański

    PROCESS_FILES(combined_ch)
}
```

Operator `.combine()` generuje iloczyn kartezjański dwóch kanałów, więc każdy element `reference_ch` zostanie sparowany z każdym elementem `input_ch`. Pozwala to procesowi działać dla każdej próbki przy jednoczesnym korzystaniu z referencji.

Wymaga to dostosowania wejścia procesu. W naszym przykładzie początek definicji procesu należałoby zmienić następująco:

```groovy title="exhausted.nf (fixed - Option 2)" hl_lines="5" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        tuple val(reference), val(sample_name)
```

Takie podejście może nie być odpowiednie we wszystkich sytuacjach.

#### Uruchom pipeline

Wypróbuj jedną z powyższych poprawek i uruchom workflow ponownie:

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

Teraz powinieneś/powinnaś zobaczyć, że wszystkie trzy próbki są przetwarzane, a nie tylko jedna.

### 2.3. Nieprawidłowa struktura zawartości kanału

Gdy workflow'y osiągają pewien poziom złożoności, trudno jest śledzić wewnętrzne struktury każdego kanału, a ludzie często generują niezgodności między tym, czego oczekuje proces, a tym, co faktycznie zawiera kanał. Jest to bardziej subtelne niż problem omówiony wcześniej, gdzie liczba kanałów była nieprawidłowa. W tym przypadku możesz mieć prawidłową liczbę kanałów wejściowych, ale wewnętrzna struktura jednego lub więcej z nich nie odpowiada temu, czego oczekuje proces.

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

Nawiasy kwadratowe w komunikacie o błędzie dają wskazówkę — proces traktuje krotkę jako pojedynczą wartość, co nie jest zamierzonym zachowaniem. Przyjrzyjmy się `bad_channel_shape.nf`:

```groovy title="bad_channel_shape.nf" hl_lines="5 20-22" linenums="1"
#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Oczekuje pojedynczej wartości, otrzymuje krotkę

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

Widać, że generujemy kanał złożony z krotek: `['sample1', 'file1.txt']`, ale proces oczekuje pojedynczej wartości `val sample_name`. Wykonane polecenie pokazuje, że proces próbuje utworzyć plik o nazwie `[sample3, file3.txt]_output.txt`, co nie jest zamierzonym wyjściem.

#### Napraw kod

Aby to naprawić, jeśli proces wymaga obu wejść, możemy dostosować go do przyjmowania krotki:

=== "Opcja 1: Przyjmij krotkę w procesie"

    === "Po"

        ```groovy title="bad_channel_shape.nf" hl_lines="5"  linenums="1"
        #!/usr/bin/env nextflow

        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)  // Naprawiono: Przyjmij krotkę

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
                val sample_name  // Oczekuje pojedynczej wartości, otrzymuje krotkę

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

=== "Opcja 2: Wyodrębnij pierwszy element"

    === "Po"

        ```groovy title="bad_channel_shape.nf" hl_lines="9" linenums="16"
        workflow {

            // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
            input_ch = channel.of(
              ['sample1', 'file1.txt'],
              ['sample2', 'file2.txt'],
              ['sample3', 'file3.txt']
            )
            PROCESS_FILES(input_ch.map { it[0] })  // Naprawiono: Wyodrębnij pierwszy element
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

Wybierz jedno z rozwiązań i uruchom workflow ponownie:

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

Najpotężniejszym narzędziem debugowania kanałów jest operator `.view()`. Pozwala on zrozumieć kształt kanałów na wszystkich etapach, co pomaga w debugowaniu.

#### Uruchom pipeline

Uruchom `bad_channel_shape_viewed.nf`, aby zobaczyć to w działaniu:

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

Przyjrzyjmy się `bad_channel_shape_viewed.nf`, aby zobaczyć, jak używane jest `.view()`:

```groovy title="bad_channel_shape_viewed.nf" linenums="16" hl_lines="9 11"
workflow {

    // Kanał emituje krotki, ale proces oczekuje pojedynczych wartości
    input_ch = channel.of(
      ['sample1', 'file1.txt'],
      ['sample2', 'file2.txt'],
      ['sample3', 'file3.txt']
    )
    .view { "Channel content: $it" }  // Debug: Pokaż oryginalną zawartość kanału
    .map { tuple -> tuple[0] }        // Transformacja: Wyodrębnij pierwszy element
    .view { "After mapping: $it" }    // Debug: Pokaż przetransformowaną zawartość kanału

    PROCESS_FILES(input_ch)
}
```

#### Napraw kod

Aby uniknąć nadmiernego używania operacji `.view()` w przyszłości do rozumienia zawartości kanałów, warto dodać komentarze:

```groovy title="bad_channel_shape_viewed.nf (with comments)" linenums="16" hl_lines="8 9"
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

Stanie się to coraz ważniejsze w miarę jak Twoje workflow'y będą rosnąć w złożoności, a struktura kanałów stanie się mniej przejrzysta.

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

### Podsumowanie

Wiele błędów struktury kanałów można stworzyć przy użyciu prawidłowej składni Nextflow. Możesz je debugować, rozumiejąc przepływ danych, używając operatorów `.view()` do inspekcji i rozpoznając wzorce komunikatów o błędach, takie jak nawiasy kwadratowe wskazujące na nieoczekiwane struktury krotek.

### Co dalej?

Naucz się o błędach tworzonych przez definicje procesów.

---

## 3. Błędy struktury procesów

Większość błędów związanych z procesami dotyczy pomyłek w formułowaniu polecenia lub problemów z oprogramowaniem. Podobnie jak w przypadku problemów z kanałami, możesz jednak popełniać błędy w definicji procesu, które nie kwalifikują się jako błędy składni, ale spowodują błędy w czasie wykonania.

### 3.1. Brakujące pliki wyjściowe

Częstym błędem przy pisaniu procesów jest stworzenie niezgodności między tym, czego oczekuje proces, a tym, co jest generowane.

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

Komunikat o błędzie wskazuje, że proces oczekiwał pliku wyjściowego o nazwie `sample3.txt`, ale skrypt faktycznie tworzy `sample3_output.txt`. Przyjrzyjmy się definicji procesu w `missing_output.nf`:

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

Widać niezgodność między nazwą pliku wyjściowego w bloku `output:` a tą użytą w skrypcie. Powoduje to błąd procesu. Jeśli napotkasz taki błąd, wróć i sprawdź, czy wyjścia są zgodne między definicją procesu a blokiem output.

Jeśli problem nadal nie jest jasny, sprawdź sam katalog roboczy, aby zidentyfikować faktycznie utworzone pliki wyjściowe:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

W tym przykładzie pokazałoby nam to, że sufiks `_output` jest włączany do nazwy pliku wyjściowego, wbrew naszej definicji `output:`.

#### Napraw kod

Napraw niezgodność, ujednolicając nazwę pliku wyjściowego:

=== "Po"

    ```groovy title="missing_output.nf" hl_lines="6 10" linenums="3"
    process PROCESS_FILES {
        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"  // Naprawiono: Dopasuj do wyjścia skryptu

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
        path "${sample_name}.txt"  // Oczekuje: sample3.txt

        script:
        """
        echo "Processing ${sample_name}" > ${sample_name}_output.txt  // Tworzy: sample3_output.txt
        """
    }
    ```

#### Uruchom pipeline

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

Inną klasą błędów są pomyłki w dostarczaniu oprogramowania. `missing_software.nf` to składniowo poprawny workflow, ale zależy od zewnętrznego oprogramowania dostarczającego polecenie `cowpy`, którego używa.

#### Uruchom pipeline

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

Proces nie ma dostępu do polecenia, które określamy. Czasami dzieje się tak, ponieważ skrypt jest obecny w katalogu `bin` workflow'u, ale nie został uczyniony wykonywalnym. Innym razem oprogramowanie nie jest zainstalowane w kontenerze lub środowisku, w którym działa workflow.

#### Sprawdź kod

Zwróć uwagę na kod wyjścia `127` — mówi Ci dokładnie o problemie. Przyjrzyjmy się `missing_software.nf`:

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

Byliśmy tu trochę nieuczciwi — w kodzie tak naprawdę nie ma nic złego. Wystarczy określić niezbędną konfigurację, aby uruchomić proces w sposób zapewniający dostęp do danego polecenia. W tym przypadku proces ma definicję kontenera, więc wystarczy uruchomić workflow z włączonym Docker.

#### Uruchom pipeline

Przygotowaliśmy dla Ciebie profil Docker w `nextflow.config`, więc możesz uruchomić workflow za pomocą:

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

!!! note "Uwaga"

    Aby dowiedzieć się więcej o tym, jak Nextflow używa kontenerów, zobacz [Hello Nextflow](../hello_nextflow/05_hello_containers.md)

### 3.3. Nieprawidłowa konfiguracja zasobów

W środowisku produkcyjnym będziesz konfigurować zasoby dla swoich procesów. Na przykład `memory` definiuje maksymalną ilość pamięci dostępną dla procesu — jeśli ją przekroczy, harmonogram zazwyczaj zabija proces i zwraca kod wyjścia `137`. Nie możemy tego zademonstrować tutaj, ponieważ używamy executora `local`, ale możemy pokazać coś podobnego z `time`.

#### Uruchom pipeline

`bad_resources.nf` ma konfigurację procesu z nierealistycznym limitem czasu wynoszącym 1 milisekundę:

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

Przyjrzyjmy się `bad_resources.nf`:

```groovy title="bad_resources.nf" linenums="3" hl_lines="3"
process PROCESS_FILES {

    time '1 ms'  // BŁĄD: Nierealistyczny limit czasu

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  // Trwa 1 sekundę, ale limit czasu to 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}
```

Wiemy, że proces zajmie więcej niż sekundę (dodaliśmy sleep, aby mieć pewność), ale jest ustawiony na przekroczenie limitu czasu po 1 milisekundzie. Ktoś był trochę nierealistyczny w swojej konfiguracji!

#### Napraw kod

Zwiększ limit czasu do realistycznej wartości:

=== "Po"

    ```groovy title="bad_resources.nf" hl_lines="3" linenums="3"
    process PROCESS_FILES {

        time '100 s'  // Naprawiono: Realistyczny limit czasu

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

        time '1 ms'  // BŁĄD: Nierealistyczny limit czasu

        input:
        val sample_name

        output:
        path "${sample_name}_output.txt"

        script:
        """
        sleep 1  // Trwa 1 sekundę, ale limit czasu to 1ms
        cowpy ${sample_name} > ${sample_name}_output.txt
        """
    }
    ```

#### Uruchom pipeline

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

Jeśli będziesz uważnie czytać komunikaty o błędach, takie awarie nie powinny Cię długo zastanawiać. Upewnij się jednak, że rozumiesz wymagania zasobowe uruchamianych poleceń, aby móc odpowiednio skonfigurować dyrektywy zasobów.

### 3.4. Techniki debugowania procesów

Gdy procesy zawodzą lub zachowują się nieoczekiwanie, potrzebujesz systematycznych technik do zbadania, co poszło nie tak. Katalog roboczy zawiera wszystkie informacje potrzebne do debugowania wykonania procesu.

#### Używanie inspekcji katalogu roboczego

Najpotężniejszym narzędziem debugowania procesów jest badanie katalogu roboczego. Gdy proces zawodzi, Nextflow tworzy katalog roboczy dla tego konkretnego wykonania procesu, zawierający wszystkie pliki potrzebne do zrozumienia, co się stało.

#### Uruchom pipeline

Użyjmy przykładu `missing_output.nf` z wcześniejszego, aby zademonstrować inspekcję katalogu roboczego (wygeneruj ponownie niezgodność nazw plików wyjściowych, jeśli potrzebujesz):

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

Gdy otrzymasz ten błąd, katalog roboczy zawiera wszystkie informacje do debugowania. Znajdź ścieżkę katalogu roboczego z komunikatu o błędzie i sprawdź jego zawartość:

```bash
# Znajdź katalog roboczy z komunikatu o błędzie
ls work/02/9604d49fb8200a74d737c72a6c98ed/
```

Następnie możesz zbadać kluczowe pliki:

##### Sprawdź skrypt polecenia

Plik `.command.sh` pokazuje dokładnie, jakie polecenie zostało wykonane:

```bash
# Wyświetl wykonane polecenie
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.sh
```

Ujawnia to:

- **Podstawianie zmiennych**: Czy zmienne Nextflow zostały prawidłowo rozwinięte
- **Ścieżki plików**: Czy pliki wejściowe zostały prawidłowo zlokalizowane
- **Struktura polecenia**: Czy składnia skryptu jest poprawna

Typowe problemy do sprawdzenia:

- **Brakujące cudzysłowy**: Zmienne zawierające spacje wymagają odpowiedniego cytowania
- **Nieprawidłowe ścieżki plików**: Pliki wejściowe, które nie istnieją lub są w złych lokalizacjach
- **Nieprawidłowe nazwy zmiennych**: Literówki w odwołaniach do zmiennych
- **Brakująca konfiguracja środowiska**: Polecenia zależne od określonych środowisk

##### Sprawdź wyjście błędów

Plik `.command.err` zawiera faktyczne komunikaty o błędach:

```bash
# Wyświetl wyjście błędów
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.err
```

Ten plik pokaże:

- **Kody wyjścia**: 127 (polecenie nie znalezione), 137 (zabity), itp.
- **Błędy uprawnień**: Problemy z dostępem do plików
- **Błędy oprogramowania**: Komunikaty o błędach specyficzne dla aplikacji
- **Błędy zasobów**: Przekroczenie limitu pamięci/czasu

##### Sprawdź standardowe wyjście

Plik `.command.out` pokazuje, co wyprodukowało Twoje polecenie:

```bash
# Wyświetl standardowe wyjście
cat work/02/9604d49fb8200a74d737c72a6c98ed/.command.out
```

Pomaga to zweryfikować:

- **Oczekiwane wyjście**: Czy polecenie wyprodukowało właściwe wyniki
- **Częściowe wykonanie**: Czy polecenie zaczęło się, ale nie powiodło się w połowie
- **Informacje diagnostyczne**: Wszelkie dane diagnostyczne z Twojego skryptu

##### Sprawdź kod wyjścia

Plik `.exitcode` zawiera kod wyjścia procesu:

```bash
# Wyświetl kod wyjścia
cat work/*/*/.exitcode
```

Typowe kody wyjścia i ich znaczenie:

- **Kod wyjścia 127**: Polecenie nie znalezione — sprawdź instalację oprogramowania
- **Kod wyjścia 137**: Proces zabity — sprawdź limity pamięci/czasu

##### Sprawdź istnienie plików

Gdy procesy zawodzą z powodu brakujących plików wyjściowych, sprawdź, jakie pliki zostały faktycznie utworzone:

```bash
# Wylistuj wszystkie pliki w katalogu roboczym
ls -la work/02/9604d49fb8200a74d737c72a6c98ed/
```

Pomaga to zidentyfikować:

- **Niezgodności nazw plików**: Pliki wyjściowe o innych nazwach niż oczekiwano
- **Problemy z uprawnieniami**: Pliki, których nie można było utworzyć
- **Problemy ze ścieżkami**: Pliki utworzone w złych katalogach

We wcześniejszym przykładzie potwierdziło nam to, że oczekiwany `sample3.txt` nie był obecny, ale `sample3_output.txt` był:

```bash
❯ ls -h work/02/9604d49fb8200a74d737c72a6c98ed
sample3_output.txt
```

### Podsumowanie

Debugowanie procesów wymaga badania katalogów roboczych w celu zrozumienia, co poszło nie tak. Kluczowe pliki to `.command.sh` (wykonany skrypt), `.command.err` (komunikaty o błędach) i `.command.out` (standardowe wyjście). Kody wyjścia takie jak 127 (polecenie nie znalezione) i 137 (proces zabity) dostarczają natychmiastowych wskazówek diagnostycznych dotyczących rodzaju awarii.

### Co dalej?

Poznaj wbudowane narzędzia debugowania Nextflow i systematyczne podejścia do rozwiązywania problemów.

---

## 4. Wbudowane narzędzia debugowania i zaawansowane techniki

Nextflow udostępnia kilka potężnych wbudowanych narzędzi do debugowania i analizowania wykonania workflow'u. Pomagają one zrozumieć, co poszło nie tak, gdzie to się stało i jak to efektywnie naprawić.

### 4.1. Wyjście procesu w czasie rzeczywistym

Czasami musisz zobaczyć, co dzieje się wewnątrz uruchomionych procesów. Możesz włączyć wyjście procesu w czasie rzeczywistym, które pokazuje dokładnie, co robi każde zadanie podczas wykonywania.

#### Uruchom pipeline

`bad_channel_shape_viewed.nf` z naszych wcześniejszych przykładów drukował zawartość kanału za pomocą `.view()`, ale możemy również użyć dyrektywy `debug`, aby wyświetlać zmienne z wnętrza samego procesu — demonstrujemy to w `bad_channel_shape_viewed_debug.nf`. Uruchom workflow:

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

Przyjrzyjmy się `bad_channel_shape_viewed_debug.nf`, aby zobaczyć, jak działa dyrektywa `debug`:

```groovy title="bad_channel_shape_viewed_debug.nf" linenums="3" hl_lines="2"
process PROCESS_FILES {
    debug true  // Włącz wyjście w czasie rzeczywistym

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

Czasami chcesz wykryć problemy zanim uruchomią się jakiekolwiek procesy. Nextflow udostępnia flagę do tego rodzaju proaktywnego debugowania: `-preview`.

#### Uruchom pipeline

Tryb podglądu pozwala testować logikę workflow'u bez wykonywania poleceń. Może być bardzo przydatny do szybkiego sprawdzania struktury workflow'u i upewniania się, że procesy są prawidłowo połączone, bez uruchamiania żadnych rzeczywistych poleceń.

!!! note "Uwaga"

    Jeśli wcześniej naprawiłeś/naprawiłaś `bad_syntax.nf`, przywróć błąd składni, usuwając nawias zamykający po bloku skryptu przed uruchomieniem tego polecenia.

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

Tryb podglądu jest szczególnie przydatny do wczesnego wykrywania błędów składni bez uruchamiania żadnych procesów. Weryfikuje strukturę workflow'u i połączenia procesów przed wykonaniem.

### 4.3. Stub running do testowania logiki

Czasami błędy są trudne do debugowania, ponieważ polecenia trwają zbyt długo, wymagają specjalnego oprogramowania lub zawodzą z złożonych powodów. Stub running pozwala testować logikę workflow'u bez wykonywania rzeczywistych poleceń.

#### Uruchom pipeline

Podczas tworzenia procesu Nextflow możesz użyć dyrektywy `stub`, aby zdefiniować „fikcyjne" polecenia generujące wyjścia o właściwej formie bez uruchamiania rzeczywistego polecenia. Takie podejście jest szczególnie wartościowe, gdy chcesz zweryfikować poprawność logiki workflow'u przed zmierzeniem się ze złożonością rzeczywistego oprogramowania.

Pamiętasz nasz `missing_software.nf` z wcześniejszego? Ten, w którym brakowało oprogramowania uniemożliwiającego uruchomienie workflow'u, dopóki nie dodaliśmy `-profile docker`? `missing_software_with_stub.nf` to bardzo podobny workflow. Jeśli uruchomimy go w ten sam sposób, wygenerujemy ten sam błąd:

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

Przyjrzyjmy się `missing_software_with_stub.nf`:

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

W porównaniu z `missing_software.nf`, ten proces ma dyrektywę `stub:` określającą polecenie, które ma być użyte zamiast tego z `script:`, gdy Nextflow jest uruchamiany w trybie stub.

Polecenie `touch`, którego tu używamy, nie zależy od żadnego oprogramowania ani odpowiednich wejść i będzie działać w każdej sytuacji, pozwalając nam debugować logikę workflow'u bez martwienia się o wewnętrzne działanie procesu.

**Stub running pomaga debugować:**

- Strukturę kanałów i przepływ danych
- Połączenia procesów i zależności
- Propagację parametrów
- Logikę workflow'u bez zależności od oprogramowania

### 4.4. Systematyczne podejście do debugowania

Teraz, gdy poznałeś/poznałaś poszczególne techniki debugowania — od katalogów roboczych po tryb podglądu, stub running i monitorowanie zasobów — połączmy je w systematyczną metodologię. Ustrukturyzowane podejście zapobiega przytłoczeniu złożonymi błędami i zapewnia, że nie przeoczysz ważnych wskazówek.

Ta metodologia łączy wszystkie omówione narzędzia w efektywny workflow:

**Czterofazowa metoda debugowania:**

**Faza 1: Rozwiązywanie błędów składni (5 minut)**

1. Sprawdź czerwone podkreślenia w VSCode lub Twoim IDE
2. Uruchom `nextflow run workflow.nf -preview`, aby zidentyfikować problemy ze składnią
3. Napraw wszystkie błędy składni (brakujące nawiasy klamrowe, końcowe przecinki itp.)
4. Upewnij się, że workflow parsuje się pomyślnie przed kontynuowaniem

**Faza 2: Szybka ocena (5 minut)**

1. Uważnie przeczytaj komunikaty o błędach w czasie wykonania
2. Sprawdź, czy to błąd wykonania, logiki czy zasobów
3. Użyj trybu podglądu do testowania podstawowej logiki workflow'u

**Faza 3: Szczegółowe badanie (15–30 minut)**

1. Znajdź katalog roboczy nieudanego zadania
2. Zbadaj pliki dziennika
3. Dodaj operatory `.view()` do inspekcji kanałów
4. Użyj `-stub-run` do testowania logiki workflow'u bez wykonania

**Faza 4: Naprawa i walidacja (15 minut)**

1. Wprowadź minimalne, ukierunkowane poprawki
2. Przetestuj z resume: `nextflow run workflow.nf -resume`
3. Zweryfikuj kompletne wykonanie workflow'u

!!! tip "Używanie resume do efektywnego debugowania"

    Po zidentyfikowaniu problemu potrzebujesz efektywnego sposobu testowania poprawek bez marnowania czasu na ponowne uruchamianie pomyślnych części workflow'u. Funkcjonalność `-resume` Nextflow jest nieoceniona przy debugowaniu.

    Zetknąłeś/zetknęłaś się z `-resume` podczas pracy z [Hello Nextflow](../hello_nextflow/), i ważne jest, abyś z niego korzystał/korzystała podczas debugowania, aby zaoszczędzić czas oczekiwania na uruchomienie procesów poprzedzających problematyczny.

    **Strategia debugowania z resume:**

    1. Uruchom workflow do momentu awarii
    2. Zbadaj katalog roboczy nieudanego zadania
    3. Napraw konkretny problem
    4. Wznów, aby przetestować tylko poprawkę
    5. Powtarzaj, aż workflow zakończy się pomyślnie

#### Profil konfiguracji debugowania

Aby to systematyczne podejście było jeszcze bardziej efektywne, możesz utworzyć dedykowaną konfigurację debugowania, która automatycznie włącza wszystkie potrzebne narzędzia:

```groovy title="nextflow.config (debug profile)" linenums="1"
profiles {
    debug {
        process {
            debug = true
            cleanup = false

            // Konserwatywne zasoby do debugowania
            maxForks = 1
            memory = '2.GB'
            cpus = 1
        }
    }
}
```

Następnie możesz uruchomić pipeline z tym profilem:

```bash
nextflow run workflow.nf -profile debug
```

Ten profil włącza wyjście w czasie rzeczywistym, zachowuje katalogi robocze i ogranicza paralelizację dla łatwiejszego debugowania.

### 4.5. Praktyczne ćwiczenie debugowania

Czas zastosować systematyczne podejście do debugowania w praktyce. Workflow `buggy_workflow.nf` zawiera kilka typowych błędów reprezentujących rodzaje problemów, z którymi spotkasz się w rzeczywistym środowisku programistycznym.

!!! exercise "Ćwiczenie"

    Użyj systematycznego podejścia do debugowania, aby zidentyfikować i naprawić wszystkie błędy w `buggy_workflow.nf`. Ten workflow próbuje przetwarzać dane próbek z pliku CSV, ale zawiera wiele celowych błędów reprezentujących typowe scenariusze debugowania.

    Zacznij od uruchomienia workflow'u, aby zobaczyć pierwszy błąd:

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

        Ten tajemniczy błąd wskazuje na problem parsowania wokół wierszy 11–12 w bloku `params{}`. Parser v2 wykrywa problemy strukturalne wcześnie.

    Zastosuj czterofazową metodę debugowania, której się nauczyłeś/nauczyłaś:

    **Faza 1: Rozwiązywanie błędów składni**
    - Sprawdź czerwone podkreślenia w VSCode lub Twoim IDE
    - Uruchom `nextflow run workflow.nf -preview`, aby zidentyfikować problemy ze składnią
    - Napraw wszystkie błędy składni (brakujące nawiasy klamrowe, końcowe przecinki itp.)
    - Upewnij się, że workflow parsuje się pomyślnie przed kontynuowaniem

    **Faza 2: Szybka ocena**
    - Uważnie przeczytaj komunikaty o błędach w czasie wykonania
    - Zidentyfikuj, czy błędy są związane z wykonaniem, logiką czy zasobami
    - Użyj trybu `-preview` do testowania podstawowej logiki workflow'u

    **Faza 3: Szczegółowe badanie**
    - Zbadaj katalogi robocze nieudanych zadań
    - Dodaj operatory `.view()` do inspekcji kanałów
    - Sprawdź pliki dziennika w katalogach roboczych
    - Użyj `-stub-run` do testowania logiki workflow'u bez wykonania

    **Faza 4: Naprawa i walidacja**
    - Wprowadź ukierunkowane poprawki
    - Użyj `-resume` do efektywnego testowania poprawek
    - Zweryfikuj kompletne wykonanie workflow'u

    **Narzędzia debugowania do Twojej dyspozycji:**
    ```bash
    # Tryb podglądu do sprawdzania składni
    nextflow run buggy_workflow.nf -preview

    # Profil debug do szczegółowego wyjścia
    nextflow run buggy_workflow.nf -profile debug

    # Stub running do testowania logiki
    nextflow run buggy_workflow.nf -stub-run

    # Resume po poprawkach
    nextflow run buggy_workflow.nf -resume
    ```

    ??? solution "Rozwiązanie"
        `buggy_workflow.nf` zawiera 9 lub 10 odrębnych błędów (w zależności od sposobu liczenia) obejmujących wszystkie główne kategorie debugowania. Oto systematyczne omówienie każdego błędu i sposobu jego naprawienia.

        Zacznijmy od błędów składni:

        **Błąd 1: Błąd składni — końcowy przecinek**
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt",  // BŁĄD: Końcowy przecinek
        ```
        **Poprawka:** Usuń końcowy przecinek
        ```groovy linenums="21"
        output:
            path "${sample_id}_result.txt"
        ```

        **Błąd 2: Błąd składni — brakujący nawias zamykający**
        ```groovy linenums="24"
        script:
        """
        echo "Processing: ${sample}"
        cat ${input_file} > ${sample}_result.txt
        """
        // BŁĄD: Brakujący nawias zamykający dla procesu processFiles
        ```
        **Poprawka:** Dodaj brakujący nawias zamykający
        ```groovy linenums="29"
        """
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        """
        }  // Dodaj brakujący nawias zamykający
        ```

        **Błąd 3: Błąd nazwy zmiennej**
        ```groovy linenums="26"
        echo "Processing: ${sample}"     // BŁĄD: powinno być sample_id
        cat ${input_file} > ${sample}_result.txt  // BŁĄD: powinno być sample_id
        ```
        **Poprawka:** Użyj poprawnej nazwy zmiennej wejściowej
        ```groovy linenums="26"
        echo "Processing: ${sample_id}"
        cat ${input_file} > ${sample_id}_result.txt
        ```

        **Błąd 4: Błąd niezdefiniowanej zmiennej**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(sample_ids)  // BŁĄD: sample_ids niezdefiniowane
        ```
        **Poprawka:** Użyj poprawnego kanału i wyodrębnij identyfikatory próbek
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)
        ```

        W tym momencie workflow uruchomi się, ale nadal będziemy otrzymywać błędy (np. `Path value cannot be null` w `processFiles`), spowodowane złą strukturą kanału.

        **Błąd 5: Błąd struktury kanału — nieprawidłowe wyjście map**
        ```groovy linenums="83"
        .map { row -> row.sample_id }  // BŁĄD: processFiles oczekuje krotki
        ```
        **Poprawka:** Zwróć strukturę krotki, której oczekuje processFiles
        ```groovy linenums="83"
        .map { row -> [row.sample_id, file(row.fastq_path)] }
        ```

        Ale to zepsuje nasze wywołanie `heavyProcess()` powyżej, więc będziemy musieli użyć map, aby przekazać tylko identyfikatory próbek do tego procesu:

        **Błąd 6: Zła struktura kanału dla heavyProcess**
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch)  // BŁĄD: input_ch ma teraz 2 elementy na emisję — heavyProcess potrzebuje tylko 1 (pierwszego)
        ```
        **Poprawka:** Użyj poprawnego kanału i wyodrębnij identyfikatory próbek
        ```groovy linenums="87"
        heavy_ch = heavyProcess(input_ch.map{it[0]})
        ```

        Teraz docieramy dalej, ale otrzymujemy błąd `No such variable: i`, ponieważ nie poprzedziliśmy zmiennej Bash ukośnikiem.

        **Błąd 7: Błąd poprzedzania zmiennej Bash**
        ```groovy linenums="48"
        echo "Heavy computation $i for ${sample_id}"  // BŁĄD: $i nie jest poprzedzone ukośnikiem
        ```
        **Poprawka:** Poprzedź zmienną Bash ukośnikiem
        ```groovy linenums="48"
        echo "Heavy computation \${i} for ${sample_id}"
        ```

        Teraz otrzymujemy `Process exceeded running time limit (1ms)`, więc naprawiamy limit czasu dla odpowiedniego procesu:

        **Błąd 8: Błąd konfiguracji zasobów**
        ```groovy linenums="36"
        time '1 ms'  // BŁĄD: Nierealistyczny limit czasu
        ```
        **Poprawka:** Zwiększ do realistycznego limitu czasu
        ```groovy linenums="36"
        time '100 s'
        ```

        Następnie mamy błąd `Missing output file(s)` do rozwiązania:

        **Błąd 9: Niezgodność nazwy pliku wyjściowego**
        ```groovy linenums="49"
        done > ${sample_id}.txt  // BŁĄD: Nieprawidłowa nazwa pliku, powinna pasować do deklaracji output
        ```
        **Poprawka:** Dopasuj do deklaracji output
        ```groovy linenums="49"
        done > ${sample_id}_heavy.txt
        ```

        Pierwsze dwa procesy uruchomiły się, ale nie trzeci.

        **Błąd 10: Niezgodność nazwy pliku wyjściowego**
        ```groovy linenums="88"
        file_ch = channel.fromPath("*.txt") // Błąd: próba pobrania wejścia z bieżącego katalogu zamiast z procesu
        handleFiles(file_ch)
        ```
        **Poprawka:** Pobierz wyjście z poprzedniego procesu
        ```groovy linenums="88"
        handleFiles(heavyProcess.out)
        ```

        Po tym cały workflow powinien działać.

        **Kompletny poprawiony workflow:**
        ```groovy linenums="1"
        #!/usr/bin/env nextflow

        /*
        * Błędny workflow do ćwiczeń debugowania
        * Ten workflow zawiera kilka celowych błędów do celów edukacyjnych
        */

        params{
            // Parametry bez walidacji
            input: Path = 'data/sample_data.csv'
            output: String = 'results'
        }

        /*
        * Proces z niezgodnością wejścia/wyjścia
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

            // Kanał z nieprawidłowym użyciem
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

- **Błędy składni**: Brakujące nawiasy klamrowe, końcowe przecinki, niezdefiniowane zmienne
- **Błędy struktury kanałów**: Nieprawidłowe kształty danych, niezdefiniowane kanały
- **Błędy procesów**: Niezgodności nazw plików wyjściowych, poprzedzanie zmiennych
- **Błędy zasobów**: Nierealistyczne limity czasu

**Kluczowe lekcje debugowania:**

1. **Uważnie czytaj komunikaty o błędach** — często wskazują bezpośrednio na problem
2. **Używaj systematycznych podejść** — naprawiaj jeden błąd na raz i testuj z `-resume`
3. **Rozumiej przepływ danych** — błędy struktury kanałów są często najbardziej subtelne
4. **Sprawdzaj katalogi robocze** — gdy procesy zawodzą, dzienniki mówią Ci dokładnie, co poszło nie tak

---

## Podsumowanie

W tym side queście poznałeś/poznałaś zestaw systematycznych technik debugowania workflow'ów Nextflow.
Stosowanie tych technik w swojej pracy pozwoli Ci spędzać mniej czasu na walce z komputerem, szybciej rozwiązywać problemy i chronić się przed przyszłymi błędami.

### Kluczowe wzorce

**1. Jak identyfikować i naprawiać błędy składni**:

- Interpretowanie komunikatów o błędach Nextflow i lokalizowanie problemów
- Typowe błędy składni: brakujące nawiasy klamrowe, nieprawidłowe słowa kluczowe, niezdefiniowane zmienne
- Rozróżnianie między zmiennymi Nextflow (Groovy) a Bash
- Używanie funkcji rozszerzenia VS Code do wczesnego wykrywania błędów

```groovy
// Brakujący nawias — szukaj czerwonych podkreśleń w IDE
process FOO {
    script:
    """
    echo "hello"
    """
// } <-- brakuje!

// Nieprawidłowe słowo kluczowe
inputs:  // Powinno być 'input:'

// Niezdefiniowana zmienna — poprzedź ukośnikiem dla zmiennych Bash
echo "${undefined_var}"      // Zmienna Nextflow (błąd jeśli niezdefiniowana)
echo "\${bash_var}"          // Zmienna Bash (poprzedzona ukośnikiem)
```

**2. Jak debugować problemy ze strukturą kanałów**:

- Rozumienie kardynalności kanałów i problemów z wyczerpaniem
- Debugowanie niezgodności struktury zawartości kanałów
- Używanie operatorów `.view()` do inspekcji kanałów
- Rozpoznawanie wzorców błędów, takich jak nawiasy kwadratowe w wyjściu

```groovy
// Inspekcja zawartości kanału
my_channel.view { "Content: $it" }

// Konwersja kanału kolejki na kanał wartości (zapobiega wyczerpaniu)
reference_ch = channel.value('ref.fa')
// lub
reference_ch = channel.of('ref.fa').first()
```

**3. Jak rozwiązywać problemy z wykonaniem procesów**:

- Diagnozowanie błędów brakujących plików wyjściowych
- Rozumienie kodów wyjścia (127 dla brakującego oprogramowania, 137 dla problemów z pamięcią)
- Badanie katalogów roboczych i plików poleceń
- Odpowiednia konfiguracja zasobów

```bash
# Sprawdź, co faktycznie zostało wykonane
cat work/ab/cdef12/.command.sh

# Sprawdź wyjście błędów
cat work/ab/cdef12/.command.err

# Kod wyjścia 127 = polecenie nie znalezione
# Kod wyjścia 137 = zabity (limit pamięci/czasu)
```

**4. Jak używać wbudowanych narzędzi debugowania Nextflow**:

- Korzystanie z trybu podglądu i debugowania w czasie rzeczywistym
- Implementacja stub running do testowania logiki
- Stosowanie resume do efektywnych cykli debugowania
- Stosowanie czterofazowej systematycznej metodologii debugowania

!!! tip "Szybki przewodnik po debugowaniu"

    **Błędy składni?** → Sprawdź ostrzeżenia VSCode, uruchom `nextflow run workflow.nf -preview`

    **Problemy z kanałami?** → Użyj `.view()` do inspekcji zawartości: `my_channel.view()`

    **Awarie procesów?** → Sprawdź pliki w katalogu roboczym:

    - `.command.sh` — wykonany skrypt
    - `.command.err` — komunikaty o błędach
    - `.exitcode` — status wyjścia (127 = polecenie nie znalezione, 137 = zabity)

    **Tajemnicze zachowanie?** → Uruchom z `-stub-run`, aby przetestować logikę workflow'u

    **Wprowadziłeś/wprowadziłaś poprawki?** → Użyj `-resume`, aby zaoszczędzić czas: `nextflow run workflow.nf -resume`

---

### Dodatkowe zasoby

- [Przewodnik rozwiązywania problemów Nextflow](https://www.nextflow.io/docs/latest/troubleshooting.html): Oficjalna dokumentacja rozwiązywania problemów
- [Rozumienie kanałów Nextflow](https://www.nextflow.io/docs/latest/channel.html): Szczegółowe omówienie typów kanałów i ich zachowania
- [Dokumentacja dyrektyw procesów](https://www.nextflow.io/docs/latest/process.html#directives): Wszystkie dostępne opcje konfiguracji procesów
- [nf-test](https://www.nf-test.com/): Framework testowy dla pipeline'ów Nextflow
- [Społeczność Nextflow na Slack](https://www.nextflow.io/slack-invite.html): Uzyskaj pomoc od społeczności

W przypadku workflow'ów produkcyjnych rozważ:

- Skonfigurowanie [Seqera Platform](https://seqera.io/platform/) do monitorowania i debugowania na dużą skalę
- Używanie [kontenerów Wave](https://seqera.io/wave/) dla reprodukowalnych środowisk oprogramowania

**Pamiętaj:** Efektywne debugowanie to umiejętność, która doskonali się z praktyką. Systematyczna metodologia i kompleksowy zestaw narzędzi, które tu zdobyłeś/zdobyłaś, będą Ci dobrze służyć przez całą Twoją drogę z Nextflow.

---

## Co dalej?

Wróć do [menu Side Quests](../) lub kliknij przycisk w prawym dolnym rogu strony, aby przejść do następnego tematu na liście.
