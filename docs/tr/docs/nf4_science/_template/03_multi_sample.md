# Bölüm 3: Çoklu örnek toplama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bölüm 2'de, her örneği bağımsız olarak işleyen örnek başına bir işleme boru hattı oluşturdunuz.
Şimdi bunu genişleterek [Bölüm 1](01_method.md)'de ele alınan çoklu örnek {AGGREGATION_METHOD} işlemini uygulayacağız.

## Görev

Kursun bu bölümünde, iş akışını aşağıdakileri yapacak şekilde genişleteceğiz:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Bu bölüm doğrudan Bölüm 2'de üretilen iş akışı üzerine inşa edilmektedir.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Bölüm 2: Tek örnek işleme](./02_single_sample.md) bölümünü tamamladığınızı ve çalışan bir `{DOMAIN_DIR}.nf` boru hattınız olduğunu varsayar.

    Eğer Bölüm 2'yi tamamlamadıysanız veya bu bölüm için yeni başlamak istiyorsanız, Bölüm 2 çözümünü başlangıç noktanız olarak kullanabilirsiniz.
    `nf4-science/{DOMAIN_DIR}/` dizininin içinden şu komutları çalıştırın:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Bu size eksiksiz bir tek örnek işleme iş akışı sağlar.
    Aşağıdaki komutu çalıştırarak başarıyla çalıştığını test edebilirsiniz:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Ders planı

Bunu iki adıma ayırdık:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Bu, süreç komutlarını ve çıktılarını güncellemeyi kapsar.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Bu, `collect()` operatörünü {AND_OTHER_CONCEPTS} tanıtır.

!!! note "Not"

     Doğru çalışma dizininde olduğunuzdan emin olun:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

[Bölüm 1](01_method.md)'deki değiştirilmiş komutu hatırlayın:

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Sonra"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Önce"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Değişikliği doğrulamak için iş akışını çalıştırın

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Özet

Süreç komutlarını ve çıktılarını değiştirerek iş akışı davranışını nasıl uyarlayacağınızı biliyorsunuz.

### Sırada ne var?

Çoklu örnek toplama adımını ekleyin.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Toplama modülünü yazın

{MODULE_INSTRUCTIONS}

### 2.2. Örnek başına çıktıları toplayın ve toplama sürecine besleyin

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Tamamlanmış iş akışını çalıştırın

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Özet

Örnekleri ayrı ayrı işleyen ve tüm örnekler genelinde sonuçları toplayan eksiksiz bir boru hattınız var.
Örnek başına çıktıları çoklu örnek analizi için toplamak üzere `collect()` gibi kanal operatörlerini nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Bu kursu tamamladığınız için tebrikler! Öğrendiklerinizi gözden geçirmek ve sonraki adımları keşfetmek için [kurs özetine](next_steps.md) gidin.
