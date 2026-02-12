# Część 3: Agregacja wielu próbek

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tłumaczenie wspomagane przez AI - [dowiedz się więcej i zasugeruj ulepszenia](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

W Części 2 zbudowałeś pipeline przetwarzający każdą próbkę niezależnie.
Teraz rozszerzymy go, aby zaimplementować {AGGREGATION_METHOD} dla wielu próbek, zgodnie z opisem w [Części 1](01_method.md).

## Zadanie

W tej części kursu rozszerzymy workflow'a, aby wykonywał następujące operacje:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Ta część bazuje bezpośrednio na workflow'ie stworzonym w Części 2.

??? info "Jak zacząć od tej sekcji"

    Ta sekcja kursu zakłada, że ukończyłeś [Część 2: Przetwarzanie pojedynczej próbki](./02_single_sample.md) i masz działający pipeline `{DOMAIN_DIR}.nf`.

    Jeśli nie ukończyłeś Części 2 lub chcesz zacząć od nowa w tej części, możesz użyć rozwiązania z Części 2 jako punktu wyjścia.
    Uruchom poniższe polecenia z katalogu `nf4-science/{DOMAIN_DIR}/`:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Otrzymasz w ten sposób kompletny workflow'a przetwarzający pojedyncze próbki.
    Możesz sprawdzić, czy działa poprawnie, uruchamiając następujące polecenie:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Plan lekcji

Podzieliliśmy to na dwa kroki:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Obejmuje to aktualizację poleceń i wyjść procesów.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Wprowadza operator `collect()` {AND_OTHER_CONCEPTS}.

!!! note "Uwaga"

     Upewnij się, że jesteś we właściwym katalogu roboczym:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Przypomnij sobie zmodyfikowane polecenie z [Części 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Po"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Przed"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Uruchom workflow'a, aby zweryfikować modyfikację

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Podsumowanie

Wiesz, jak modyfikować polecenia i wyjścia procesów, aby dostosować zachowanie workflow'a.

### Co dalej?

Dodaj krok agregacji wielu próbek.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Napisz moduł agregacji

{MODULE_INSTRUCTIONS}

### 2.2. Zbierz wyjścia z poszczególnych próbek i przekaż je do procesu agregacji

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Uruchom ukończony workflow'a

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Wyjście polecenia"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Podsumowanie

Masz kompletny pipeline, który przetwarza próbki indywidualnie i agreguje wyniki ze wszystkich próbek.
Wiesz, jak używać operatorów kanałów takich jak `collect()` do agregowania wyjść z poszczególnych próbek w celu przeprowadzenia analizy wielu próbek.

### Co dalej?

Gratulacje ukończenia tego kursu! Przejdź do [podsumowania kursu](next_steps.md), aby przejrzeć to, czego się nauczyłeś, i poznać kolejne kroki.
