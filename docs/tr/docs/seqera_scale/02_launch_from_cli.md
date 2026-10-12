# Bölüm 2: Pipeline'ları komut satırından başlatma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[Bölüm 1](./01_run_with_seqera.md)'de nf-core/rnaseq'i Seqera web arayüzünden başlattınız.
Şimdi aynı işlemi `tw` CLI kullanarak komut satırından yapacak ve çalışma alanınıza yeni bir pipeline ekleyeceksiniz.

---

## 1. Pipeline'ları komut satırından başlatma

Çalıştırma görünümünde **Command line** sekmesine tıklayın.
Platform'un sizin adınıza oluşturup gönderdiği `nextflow run` komutunun tam halini göreceksiniz — Use nf-core kursunda manuel olarak çalıştırdığınız komutların aynısı.

Platform, Nextflow'un yerini almaz; onu düzenler.
Web arayüzü üzerinden yapabildiğiniz her şeyi, Platform API'siyle etkileşim kurmak için kullanılan komut satırı aracı olan `tw` CLI aracılığıyla bir terminalden de yapabilirsiniz.
Bu, betiklerden veya CI/CD pipeline'larından başlatmaları otomatikleştirmek için oldukça kullanışlıdır.

Bunu şimdi, önceki kurslarda kullandığınız codespace'ten yapacağız.

### 1.1. tw CLI'ı yükleme

`tw` ikili dosyasını indirip yüklemek için Codespace terminalinizde aşağıdaki komutları çalıştırın:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Kurulumu doğrulayın:

```bash
tw --version
```

??? success "Komut çıktısı"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

`tw` CLI yüklendi ve yapılandırılmaya hazır.

### 1.2. Erişim jetonu alma

`tw` CLI, Seqera ile kimlik doğrulamasını kişisel erişim jetonu kullanarak gerçekleştirir.

1. Seqera web arayüzünde sağ üst köşedeki avatarınıza tıklayın ve **Your tokens** seçeneğini seçin.
2. **Add token**'a tıklayın, bir isim verin (örneğin `training`) ve **Add**'e tıklayın.
3. Jeton değerini kopyalayın — yalnızca bir kez gösterilecektir.
   Hemen bir yere kaydetmezseniz yeni bir tane oluşturmanız gerekecektir.

### 1.3. CLI'ı yapılandırma

Kolaylık sağlamak amacıyla, az önce oluşturduğunuz erişim jetonunu ve çalışma alanı tanımlayıcısını içeren bir yapılandırma dosyası oluşturacağız.

Bu dizindeki `.seqera_config` dosyasını düzenleyicide açın ve iki değişkeni ayarlayın:

- **`TOWER_ACCESS_TOKEN`**: 1.2. bölümünde oluşturduğunuz jeton
- **`TOWER_WORKSPACE_ID`**: Çalışma alanınızın sayısal kimliği (1.4. bölümünde çalıştıracağınız `tw workspaces list` komutundaki `ID` sütunu)

Değerleri girdikten sonra yapılandırmayı yükleyin:

```bash
source .seqera_config
```

Bağlantıyı doğrulayın:

```bash
tw info
```

??? success "Komut çıktısı"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

`tw` CLI artık kimlik doğrulaması yapılmış ve Seqera hesabınıza bağlı durumda.
Yapılandırmayı yeniden yüklemek için her Codespace oturumunun başında `source .seqera_config` komutunu çalıştırın.

!!! tip "İpucu"

    Çalışma alanınızda birincil bir hesaplama ortamı tanımlı değilse, yapılandırma dosyanıza `export TOWER_COMPUTE_ENV=<compute-env-name>` satırını ekleyerek bir varsayılan belirleyebilirsiniz.
    Herhangi bir yapılandırma değeri, bayrağı açıkça geçirerek komut satırından geçersiz kılınabilir (örneğin `--compute-env other-env`).
    Seçeneklerin ve ortam değişkenlerinin tam listesi için [tw CLI referansına](https://docs.seqera.io/platform/latest/cli/reference) bakın.

### 1.4. Çalışma alanınızı CLI'dan keşfetme

Erişiminiz olan çalışma alanlarını listeleyin:

```bash
tw workspaces list
```

??? success "Komut çıktısı"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Az önce başlattığınız nf-core/rnaseq çalıştırması dahil, çalışma alanınızdaki çalıştırmaları görüntüleyin:

```bash
tw runs list
```

??? success "Komut çıktısı"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

Web arayüzünde izlediğiniz çalıştırmanın aynısı burada da görünmektedir.

!!! note "Not"

    `TOWER_WORKSPACE_ID`, `.seqera_config` dosyasında ayarlandığından tüm `tw` komutlarında `--workspace` parametresini atlayabilirsiniz.
    Yapılandırma olmadan bunu açıkça geçirmeniz gerekir:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Web arayüzünde görünen her şeye CLI'dan da erişilebilir.

### 1.5. nf-core/rnaseq'i CLI'dan başlatma

[Bölüm 1](./01_run_with_seqera.md)'de çalışma alanınıza eklediğiniz pipeline, CLI'da adıyla kullanılabilir.
`test` profiliyle başlatın:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Komut çıktısı"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Bağlantıyı tarayıcınızda açın ve çalıştırmanın **Runs** panelinde göründüğünü doğrulayın.

Çalıştığını gördükten sonra, CLI ile web arayüzünün aynı çalışma alanına ait iki farklı görünüm olduğunu doğrulamış olursunuz.

!!! note "Not"

    Pipeline'ı bir çalışma alanına eklemeden önce `tw launch` komutuna doğrudan bir GitHub URL'si de geçirebilirsiniz.
    Ancak pipeline'ı başlatmadan önce açıkça eklemek genellikle daha iyi bir yaklaşımdır: pipeline yapılandırmasını gelecekteki çalıştırmalar için kaydeder, adıyla erişilebilir kılar ve Launchpad'de tüm çalışma alanı üyelerine görünür hale getirir.

    `tw` kullanarak bir pipeline'ı doğrudan komut satırından çalışma alanına eklemek de mümkündür.
    Bir sonraki bölümde bunu nf-core/demo pipeline'ı ile nasıl yapacağınız gösterilmektedir.

### Özetle

`tw` CLI'ın kimlik doğrulamasını nasıl yapacağınızı, çalışma alanınızı nasıl inceleyeceğinizi ve terminalden kayıtlı bir pipeline'ı nasıl başlatacağınızı öğrendiniz.

### Sırada ne var?

Komut satırından çalışma alanınıza yeni bir pipeline ekleyin ve başlatın.

---

## 2. Yeni bir pipeline ekleme ve çalıştırma

GitHub'daki herhangi bir Nextflow pipeline'ı, kök dizininde `main.nf` giriş noktası ve `nextflow.config` dosyası bulunduğu sürece `tw pipelines add` komutuyla çalışma alanınıza eklenebilir.
nf-core/demo bu konuda pratik yapmak için iyi bir örnektir: Use nf-core kursunda zaten çalıştırdınız, dolayısıyla ne yaptığını ve ne beklemeniz gerektiğini biliyorsunuz.

### 2.1. nf-core/demo'yu çalışma alanınıza ekleme

Pipeline'ı çalışma alanınıza kaydetmek için aşağıdaki komutu çalıştırın:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Komut çıktısı"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

Pipeline artık kayıtlı ve Launchpad'de görünecektir.

### 2.2. Launchpad'de göründüğünü doğrulama

Eklendiğini onaylamak için çalışma alanınızdaki pipeline'ları listeleyin:

```bash
tw pipelines list
```

??? success "Komut çıktısı"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Tarayıcınızda çalışma alanınızı açın ve nf-core/demo'nun nf-core/rnaseq ile birlikte göründüğünü doğrulamak için **Launchpad**'e tıklayın.

!!! tip "İpucu"

    Pipeline'ları web arayüzü üzerinden de ekleyebilirsiniz: sol kenar çubuğunda **Launchpad**'e, ardından **Add pipeline**'a tıklayın ve formu uygun şekilde doldurun.

nf-core/demo girişindeki **Launch** düğmesine tıklayarak başlatma formunu açın.
`input` ve `outdir` parametrelerinin kırmızıyla vurgulandığını göreceksiniz — bunlar varsayılan değeri olmayan zorunlu alanlardır; çünkü `tw pipelines add` yalnızca pipeline kaynağını kaydeder, herhangi bir parametre ön yapılandırması yapmaz.
Sonraki iki bölümde bu değerlerin nasıl sağlanacağı anlatılmaktadır: önce web formu aracılığıyla, ardından komut satırından.

### 2.3. nf-core/demo'yu web arayüzünden başlatma

Başlatma formu açıkken iki zorunlu parametreyi doldurun.

`input` için nf-core/demo test profilindeki test örnek sayfası URL'sini girin.
Bunu, Use nf-core kursunda incelediğiniz pipeline deposundaki `conf/test.config` dosyasında bulabilirsiniz:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

`outdir` için pipeline'ın sonuçlarını yazabileceği bir bulut depolama yolu girin.
Çalışma alanınız için yapılandırılmış bucket'ı, çalıştırmaları düzenli tutmak amacıyla bir alt dizinle birlikte kullanın:

```
s3://my-bucket/demo-results
```

Her iki alan da doldurulduktan sonra mavi **Launch** düğmesine tıklayın.

Çalıştırma **Runs** panelinde görünür ve test veri kümesinde birkaç dakika içinde tamamlanır.
Görev tablosunu ve varsa yürütme raporlarını incelemek için çalıştırmaya tıklayın.

### 2.4. nf-core/demo'yu CLI'dan başlatma

`nextflow run` komutunun aksine, `tw launch` komutu `--input` veya `--outdir` gibi bireysel parametre bayraklarını kabul etmez.
Parametreler, `--params-file` ile geçirilen YAML veya JSON biçimindeki bir dosya aracılığıyla sağlanmalıdır.
Bu yaklaşım tekrarlanabilirliği teşvik eder: kaydedilmiş bir parametre dosyası, bir çalıştırma için kullanılan değerleri tam olarak belgeler ve bir çalıştırma yapılandırmasını tekrarlamayı ya da paylaşmayı kolaylaştırır.

Çalışma dizininizde bir parametre dosyası oluşturun:

```bash
touch params.yaml
```

Dosyayı düzenleyicide açın ve çıktı yolunu ekleyin:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Artık `input` örnek sayfasını sağlayan `test` profilini ve `outdir`'i sağlayan parametre dosyasını kullanarak pipeline'ı başlatabilirsiniz:

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Komut çıktısı"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Çalıştırmanın **Runs** panelinde göründüğünü doğrulamak için bağlantıyı açın.

!!! tip "İpucu"

    Bazı varsayılanlar belirlemek ve web formu aracılığıyla daha önce yaptıklarımızla eşleşmesi için bazı ek özellikler eklemek isterseniz, parametre dosyasını ilk kurulum adımına dahil edebilirsiniz:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Özetle

GitHub'da barındırılan herhangi bir Nextflow pipeline'ını çalışma alanınıza nasıl ekleyeceğinizi ve başlatacağınızı öğrendiniz; hem parametreleri manuel olarak doldurarak web arayüzünden, hem de bir profili parametre dosyasıyla birleştirerek `tw` CLI'dan.

---

## Özet

Bu bölümde şunları öğrendiniz:

- `tw` CLI'ın kimlik doğrulamasını yapma ve terminalden kayıtlı bir pipeline'ı başlatma
- CLI kullanarak GitHub'dan yeni bir pipeline ekleme ve Launchpad'de göründüğünü doğrulama
- Zorunlu parametreleri manuel olarak doldurarak Seqera web arayüzünden pipeline başlatma
- Bir Nextflow profili ve parametre dosyası kullanarak CLI'dan pipeline başlatma
