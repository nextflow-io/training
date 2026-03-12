# Bölüm 2: Tek örnek işleme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bölüm 1'de, {TOOL_A} ve {TOOL_B} komutlarını ilgili konteynerlerinde manuel olarak test ettiniz.
Şimdi aynı komutları bir Nextflow iş akışına dönüştüreceğiz.

## Görev

Kursun bu bölümünde, aşağıdaki işlemleri gerçekleştiren bir iş akışı geliştireceğiz:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Bu, Bölüm 1'deki adımları tekrarlar; orada bu komutları konteynerlerinde manuel olarak çalıştırmıştınız.

Başlangıç noktası olarak, size iş akışının ana bölümlerini özetleyen bir iş akışı dosyası (`{DOMAIN_DIR}.nf`) ve modüllerin yapısını özetleyen iki modül dosyası ({TOOL_A_MODULE}.nf ve {TOOL_B_MODULE}.nf) sağlıyoruz.
Bu dosyalar işlevsel değildir; amaçları sadece kodun ilginç kısımlarını doldurmanız için iskelet görevi görmektir.

## Ders planı

Geliştirme sürecini daha eğitici hale getirmek için bunu {N} adıma ayırdık:

1. **{TOOL_A_ACTION} çalıştıran tek aşamalı bir iş akışı yazın.**
   Bu, bir modül oluşturmayı, içe aktarmayı ve bir iş akışında çağırmayı kapsar.
2. **{TOOL_B_ACTION} çalıştırmak için ikinci bir süreç ekleyin.**
   Bu, süreç çıktılarını girdilere zincirlemeyi ve yardımcı dosyaları işlemeyi tanıtır.
3. **İş akışını bir örnek grubu üzerinde çalışacak şekilde uyarlayın.**
   Bu, paralel yürütmeyi kapsar ve ilişkili dosyaları bir arada tutmak için demetleri tanıtır.
4. **İş akışını toplu girdi dosyaları içeren bir örnek tablosu kabul edecek şekilde düzenleyin.**
   Bu, toplu olarak girdi sağlamak için yaygın bir deseni gösterir.

Her adım, iş akışı geliştirmenin belirli bir yönüne odaklanır.

---

## 1. {TOOL_A_ACTION} çalıştıran tek aşamalı bir iş akışı yazın

Bu ilk adım temellere odaklanır: {PRIMARY_INPUT_TYPE} yükleme ve {TOOL_A_OUTPUT_DESCRIPTION}.

[Bölüm 1](01_method.md)'deki `{TOOL_A_COMMAND_NAME}` komutunu hatırlayın:

```bash
{TOOL_A_COMMAND_SUMMARY}
```

Komut {INPUT_DESCRIPTION} alır ve {OUTPUT_DESCRIPTION} üretir.
Konteyner URI'si `{TOOL_A_CONTAINER_URI}` idi.

Bu bilgiyi alıp üç aşamada Nextflow'a dönüştüreceğiz:

1. Girdiyi ayarlayın
2. Süreci yazın ve iş akışında çağırın
3. Çıktı işlemeyi yapılandırın

### 1.1. Girdiyi ayarlayın

Bir girdi parametresi bildirmemiz, uygun bir varsayılan değer sağlamak için bir test profili oluşturmamız ve bir girdi kanalı oluşturmamız gerekiyor.

#### 1.1.1. Bir girdi parametresi bildirimi ekleyin

Ana iş akışı dosyası `{DOMAIN_DIR}.nf`'de, `Pipeline parameters` bölümü altında, `{PRIMARY_PARAM_NAME}` adlı bir CLI parametresi bildirin.

=== "Sonra"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Birincil girdi
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Önce"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Birincil girdi
    ```

Bu, CLI parametresini ayarlar; ancak geliştirme sırasında iş akışını her çalıştırdığımızda dosya yolunu yazmak istemiyoruz.
Varsayılan bir değer sağlamak için birden fazla seçenek vardır; burada bir test profili kullanıyoruz.

#### 1.1.2. `nextflow.config` dosyasında varsayılan değere sahip bir test profili oluşturun

Test profili, komut satırında girdi belirtmeden bir iş akışını denemek için uygun varsayılan değerler sağlar.
Bu, Nextflow ekosisteminde yaygın bir kuraldır (daha fazla ayrıntı için [Hello Config](../../hello_nextflow/06_hello_config.md)'e bakın).

`nextflow.config` dosyasına, `{PRIMARY_PARAM_NAME}` parametresini test {PRIMARY_INPUT_TYPE}lerinden birine ayarlayan bir `test` profili içeren bir `profiles` bloğu ekleyin.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Burada, iş akışı betiğinin bulunduğu dizini işaret eden yerleşik bir Nextflow değişkeni olan `${projectDir}` kullanıyoruz.
Bu, mutlak yolları sabit kodlamadan veri dosyalarına ve diğer kaynaklara başvurmayı kolaylaştırır.

#### 1.1.3. Girdi kanalını ayarlayın

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. {TOOL_A_NAME} modülünü yazın

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Modülü içe aktarın ve iş akışında çağırın

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. İş akışını çalıştırın

Bu noktada, tamamen işlevsel olması gereken tek adımlı bir iş akışımız var.

Test profilinde ayarlanan varsayılan değeri kullanmak ve komut satırına yol yazmak zorunda kalmamak için `-profile test` ile çalıştırabiliriz.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Özet

Bir süreç içeren bir modül oluşturmayı, bunu bir iş akışına içe aktarmayı, bir girdi kanalıyla çağırmayı ve sonuçları yayınlamayı biliyorsunuz.

### Sırada ne var?

Ek analiz adımlarını zincirlemek için ikinci bir süreç ekleyin.

---

## 2. {TOOL_B_ACTION} çalıştırmak için ikinci bir süreç ekleyin

{DESCRIBE_WHAT_THIS_STEP_ADDS}

[Bölüm 1](01_method.md)'deki `{TOOL_B_COMMAND_NAME}` komutunu hatırlayın:

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. {TOOL_B_NAME} modülünü yazın

{MODULE_INSTRUCTIONS}

### 2.2. Modülü içe aktarın ve iş akışında çağırın

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. İş akışını çalıştırın

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Özet

Süreç çıktılarını girdilere zincirlemeyi ve iş akışında yardımcı dosyaları işlemeyi biliyorsunuz.

### Sırada ne var?

İş akışını birden fazla örneği paralel olarak işleyecek şekilde ölçeklendirin.

---

## 3. İş akışını bir örnek grubu üzerinde çalışacak şekilde uyarlayın

Şu ana kadar, iş akışı tek bir örneği işliyor.
Birden fazla örneği işlemek için, girdilerin nasıl sağlandığını değiştirmemiz ve yürütmeyi paralelleştirmek için Nextflow'un veri akışı paradigmasından yararlanmamız gerekiyor.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. İş akışını birden fazla örnek üzerinde çalıştırın

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Komut çıktısı"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Özet

Birden fazla girdi örneği üzerinde örnek başına işlemeyi paralelleştirmek için Nextflow'un veri akışı paradigmasından nasıl yararlanacağınızı biliyorsunuz.

### Sırada ne var?

[Bölüm 3](03_multi_sample.md)'te, tüm örneklerdeki sonuçları birleştirmek için çok örnekli toplama ekleyeceksiniz.
