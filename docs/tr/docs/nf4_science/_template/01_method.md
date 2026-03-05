# Bölüm 1: Yönteme genel bakış ve manuel test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![Boru hattına genel bakış](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Yöntemler

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Herhangi bir iş akışı kodu yazmaya başlamadan önce, komutları bazı test verileri üzerinde manuel olarak deneyeceğiz.

### Veri seti

Aşağıdaki veri ve ilgili kaynakları sağlıyoruz:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Yazılım

İlgili ana araçlar [{TOOL_A}]({TOOL_A_URL}) ve [{TOOL_B}]({TOOL_B_URL})'dir.

Bu araçlar GitHub Codespaces ortamında yüklü değildir, bu nedenle onları konteynerlar aracılığıyla kullanacağız (bkz. [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Not"

    `nf4-science/{DOMAIN_DIR}` dizininde olduğunuzdan emin olun, böylece `pwd` yazdığınızda gösterilen yolun son kısmı `{DOMAIN_DIR}` olsun.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

Bu bölümde tek örnekli işleme yaklaşımını oluşturan komutları test ediyoruz.
Bunlar, bu kursun 2. Bölümünde bir Nextflow iş akışına sarmalayacağımız komutlardır.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Komutları sadece bir örnek üzerinde test ederek başlıyoruz.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Konteyneri çekin

Konteyner imajını indirmek için `docker pull` komutunu çalıştırın:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Konteyneri etkileşimli olarak başlatın

Konteyneri başlatın ve araçların girdi dosyalarına erişebilmesi için `data` dizinini bağlayın:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

İsteminiz, konteynerin içinde olduğunuzu belirtecek şekilde değişir.

#### 1.1.3. Komutu çalıştırın

```bash
{TOOL_A_COMMAND}
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmüş olmalıdır.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Konteyneri çekin

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Konteyneri etkileşimli olarak başlatın

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Komutu çalıştırın

```bash
{TOOL_B_COMMAND}
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Konteynerden çıkın

Konteynerden çıkmak için `exit` yazın.

```bash
exit
```

İsteminiz normale dönmüş olmalıdır.
Bu, tek örnekli işleme testini tamamlar.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

Bu bölümde çok örnekli işleme için gereken ek komutları test ediyoruz.
Bunlar, bu kursun 3. Bölümünde bir Nextflow iş akışına sarmalayacağımız komutlardır.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Özet

{TOOL_A} ve {TOOL_B} komutlarını kendi konteynerlarında nasıl test edeceğinizi, {MULTI_SAMPLE_SUMMARY} dahil olmak üzere biliyorsunuz.

### Sırada ne var?

Aynı komutları, işi yürütmek için konteynerlar kullanan iş akışlarına nasıl sarmalayacağınızı öğrenin.
