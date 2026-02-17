# Część 2: Przetwarzanie pojedynczej próbki

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W Części 1 przetestowałeś polecenia {TOOL_A} i {TOOL_B} ręcznie w ich odpowiednich kontenerach.
Teraz opakujemy te same polecenia w workflow Nextflow'a.

## Zadanie

W tej części kursu stworzymy workflow, który wykonuje następujące czynności:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

To odtwarza kroki z Części 1, gdzie uruchamiałeś te polecenia ręcznie w ich kontenerach.

Jako punkt wyjścia dostarczamy Ci plik workflow'a, `{DOMAIN_DIR}.nf`, który przedstawia główne części workflow'a, a także dwa pliki modułów, {TOOL_A_MODULE}.nf i {TOOL_B_MODULE}.nf, które przedstawiają strukturę modułów.
Te pliki nie są funkcjonalne; ich celem jest jedynie służenie jako szkielety, które wypełnisz interesującymi częściami kodu.

## Plan lekcji

Aby proces tworzenia był bardziej edukacyjny, podzieliliśmy go na {N} kroków:

1. **Napisz jednostopniowy workflow, który uruchamia {TOOL_A_ACTION}.**
   To obejmuje tworzenie modułu, importowanie go i wywoływanie w workflow.
2. **Dodaj drugi proces, aby uruchomić {TOOL_B_ACTION}.**
   To wprowadza łączenie wyjść procesów z wejściami i obsługę plików pomocniczych.
3. **Dostosuj workflow do uruchamiania na zestawie próbek.**
   To obejmuje równoległe wykonywanie i wprowadza krotki, aby utrzymać powiązane pliki razem.
4. **Spraw, aby workflow akceptował arkusz próbek zawierający zestaw plików wejściowych.**
   To demonstruje powszechny wzorzec dostarczania wejść zbiorczo.

Każdy krok koncentruje się na konkretnym aspekcie tworzenia workflow'a.

---

## 1. Napisz jednostopniowy workflow, który uruchamia {TOOL_A_ACTION}

Ten pierwszy krok koncentruje się na podstawach: wczytywaniu {PRIMARY_INPUT_TYPE} i {TOOL_A_OUTPUT_DESCRIPTION}.

Przypomnij sobie polecenie `{TOOL_A_COMMAND_NAME}` z [Części 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

Polecenie przyjmuje {INPUT_DESCRIPTION} i produkuje {OUTPUT_DESCRIPTION}.
URI kontenera to `{TOOL_A_CONTAINER_URI}`.

Opakujemy te informacje w Nextflow'a w trzech etapach:

1. Skonfiguruj wejście
2. Napisz proces i wywołaj go w workflow
3. Skonfiguruj obsługę wyjścia

### 1.1. Skonfiguruj wejście

Musimy zadeklarować parametr wejściowy, utworzyć profil testowy, aby zapewnić wygodną wartość domyślną, i utworzyć kanał wejściowy.

#### 1.1.1. Dodaj deklarację parametru wejściowego

W głównym pliku workflow'a `{DOMAIN_DIR}.nf`, w sekcji `Pipeline parameters`, zadeklaruj parametr CLI o nazwie `{PRIMARY_PARAM_NAME}`.

=== "Po"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Przed"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

To konfiguruje parametr CLI, ale nie chcemy wpisywać ścieżki pliku za każdym razem, gdy uruchamiamy workflow podczas tworzenia.
Istnieje wiele opcji dostarczania wartości domyślnej; tutaj używamy profilu testowego.

#### 1.1.2. Utwórz profil testowy z wartością domyślną w `nextflow.config`

Profil testowy zapewnia wygodne wartości domyślne do wypróbowania workflow'a bez określania wejść w wierszu poleceń.
Jest to powszechna konwencja w ekosystemie Nextflow (zobacz [Hello Config](../../hello_nextflow/06_hello_config.md) po więcej szczegółów).

Dodaj blok `profiles` do `nextflow.config` z profilem `test`, który ustawia parametr `{PRIMARY_PARAM_NAME}` na jeden z testowych {PRIMARY_INPUT_TYPE}.

=== "Po"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Przed"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Tutaj używamy `${projectDir}`, wbudowanej zmiennej Nextflow'a, która wskazuje na katalog, w którym znajduje się skrypt workflow'a.
To ułatwia odwoływanie się do plików danych i innych zasobów bez kodowania bezwzględnych ścieżek na stałe.

#### 1.1.3. Skonfiguruj kanał wejściowy

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Napisz moduł {TOOL_A_NAME}

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Zaimportuj i wywołaj moduł w workflow

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Uruchom workflow'a

W tym momencie mamy jednostopniowy workflow, który powinien być w pełni funkcjonalny.

Możemy go uruchomić z `-profile test`, aby użyć wartości domyślnej skonfigurowanej w profilu testowym i uniknąć konieczności wpisywania ścieżki w wierszu poleceń.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Podsumowanie

Wiesz, jak utworzyć moduł zawierający proces, zaimportować go do workflow'a, wywołać go z kanałem wejściowym i opublikować wyniki.

### Co dalej?

Dodaj drugi proces, aby połączyć w łańcuch dodatkowe kroki analizy.

---

## 2. Dodaj drugi proces, aby uruchomić {TOOL_B_ACTION}

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Przypomnij sobie polecenie `{TOOL_B_COMMAND_NAME}` z [Części 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Napisz moduł {TOOL_B_NAME}

{MODULE_INSTRUCTIONS}

### 2.2. Zaimportuj i wywołaj moduł w workflow

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Uruchom workflow'a

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Podsumowanie

Wiesz, jak łączyć wyjścia procesów z wejściami i obsługiwać pliki pomocnicze w workflow.

### Co dalej?

Przeskaluj workflow, aby przetwarzać wiele próbek równolegle.

---

## 3. Dostosuj workflow do uruchamiania na zestawie próbek

Do tej pory workflow przetwarza pojedynczą próbkę.
Aby obsłużyć wiele próbek, musimy zmodyfikować sposób dostarczania wejść i wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia wykonywania.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Uruchom workflow na wielu próbkach

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Podsumowanie

Wiesz, jak wykorzystać paradygmat przepływu danych Nextflow'a do zrównoleglenia przetwarzania poszczególnych próbek na wielu próbkach wejściowych.

### Co dalej?

W [Części 3](03_multi_sample.md) dodasz agregację wielopróbkową, aby połączyć wyniki ze wszystkich próbek.
