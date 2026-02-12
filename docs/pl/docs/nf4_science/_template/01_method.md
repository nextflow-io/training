# Część 1: Przegląd metody i testy manualne

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![Przegląd pipeline'u](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Metody

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Zanim zaczniemy pisać jakikolwiek kod workflow'a, przetestujemy polecenia ręcznie na danych testowych.

### Zbiór danych

Udostępniamy następujące dane i powiązane zasoby:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Oprogramowanie

Głównymi narzędziami są [{TOOL_A}]({TOOL_A_URL}) oraz [{TOOL_B}]({TOOL_B_URL}).

Te narzędzia nie są zainstalowane w środowisku GitHub Codespaces, więc użyjemy ich przez kontenery (zobacz [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Uwaga"

    Upewnij się, że jesteś w katalogu `nf4-science/{DOMAIN_DIR}`, tak aby ostatnia część ścieżki wyświetlana po wpisaniu `pwd` to `{DOMAIN_DIR}`.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

W tej sekcji testujemy polecenia składające się na podejście przetwarzania pojedynczej próbki.
To właśnie te polecenia opakujemy w workflow Nextflow'a w części 2 tego szkolenia.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Zaczynamy od przetestowania poleceń na tylko jednej próbce.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Pobierz kontener

Uruchom polecenie `docker pull`, aby pobrać obraz kontenera:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Uruchom kontener interaktywnie

Uruchom kontener i zamontuj katalog `data`, aby narzędzia mogły uzyskać dostęp do plików wejściowych:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

Twój prompt zmieni się, wskazując, że jesteś wewnątrz kontenera.

#### 1.1.3. Uruchom polecenie

```bash
{TOOL_A_COMMAND}
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Wyjdź z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normalnego.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Pobierz kontener

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Uruchom kontener interaktywnie

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Uruchom polecenie

```bash
{TOOL_B_COMMAND}
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Wyjdź z kontenera

Aby wyjść z kontenera, wpisz `exit`.

```bash
exit
```

Twój prompt powinien wrócić do normalnego.
To kończy test przetwarzania pojedynczej próbki.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

W tej sekcji testujemy dodatkowe polecenia potrzebne do przetwarzania wielu próbek.
To właśnie te polecenia opakujemy w workflow Nextflow'a w części 3 tego szkolenia.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Podsumowanie

Wiesz, jak testować polecenia {TOOL_A} i {TOOL_B} w ich odpowiednich kontenerach, w tym jak {MULTI_SAMPLE_SUMMARY}.

### Co dalej?

Naucz się, jak opakować te same polecenia w workflow'y wykorzystujące kontenery do wykonywania pracy.
