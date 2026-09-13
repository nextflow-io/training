# Bölüm 2: İşlem Kaynaklarını ve Hataları Yönetme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[Bölüm 1](./01_packaging_and_execution.md)'de, bir pipeline'ın görevlerinin nerede ve nasıl çalışacağını uyarladınız.
Bu bölümde, her görevin ne kadar işlem kaynağı alacağını ve bir görev başarısız olduğunda ne olacağını uyarlayacaksınız.

---

## 1. İşlem Kaynağı Tahsislerini Kontrol Etme

Nextflow, varsayılan olarak `cpus` yönergesi aracılığıyla her sürece tek bir CPU tahsis eder ve siz bir sınır belirlemedikçe bellek sınırı uygulamaz:

```groovy title="Built-in configuration"
process {
    cpus = 1
}
```

[Nextflow Run](../nextflow_run/index.md) bölümünden bu pipeline'ın yapılandırmasının tüm süreçler için `memory` değerini 1 GB olarak ayarladığını zaten biliyorsunuz.
Peki kendi pipeline'larınız için hangi değerleri kullanmanız gerektiğini nasıl anlarsınız?

### 1.1. Kaynak Kullanım Raporu Oluşturma

[Nextflow Run](../nextflow_run/02_configure_pipeline.md) bölümünde `-with-report` seçeneğiyle bir yürütme raporu oluşturdunuz.
Süreçlerinizin gerçekte ne kadar CPU ve belleğe ihtiyaç duyduğunu öğrenmek için de aynı raporu kullanırsınız: iş akışını bazı varsayılan tahsislerle çalıştırın, gerçek kullanımı kaydedin, ardından buna göre ayarlayın.

```bash
nextflow run main.nf -with-report report-config-1.html
```

Rapor, bir tarayıcıda açabileceğiniz bir HTML dosyasıdır.
Tahsis edilen kaynakların gerçekte ne kadarının kullanıldığı dahil olmak üzere, süreç başına çalışma süresi ve kaynak kullanımını ayrıntılı olarak gösterir.
Mevcut varsayılanlarla (1 CPU, 1 GB bellek) `cowpy` için gösterdikleri şunlardır:

| Metrik                | Değer  |
| --------------------- | ------ |
| CPU kullanımı         | 116%   |
| Kullanılan peak bellek | 6.4 MB |
| Tahsis edilen bellek  | 1 GB   |

`cowpy`, 1 GB tahsisinin %1'inden çok daha azını kullanır; `%cpu` değerinin 100'ün üzerinde olması, konteyner içinde kısa süreli artışlarla zaman zaman birden fazla CPU'nun işlem gücünü kullandığı anlamına gelir.

Mevcut özelliklerin tam listesi için [Reports](https://nextflow.io/docs/latest/reports.html) sayfasına bakın.

### 1.2. Belirli Bir Süreç İçin Kaynak Tahsisi Belirleme

Yukarıdaki rapor, `cowpy`'nin mevcut tahsisi dahilinde rahatça çalıştığını göstermektedir; ancak diyelim ki üretim ortamında daha büyük girdiler beklediğiniz için yine de daha fazla alan vermek istiyorsunuz.
`withName` kullanarak tek bir süreç için varsayılanları geçersiz kılabilirsiniz.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

Bu yapılandırmayla birlikte, [Bölüm 1](./01_packaging_and_execution.md)'deki `conda` ayarına ek olarak 2 GB ve 2 CPU talep eden `cowpy` dışındaki her süreç 1 GB bellek ve tek bir CPU talep eder.

!!! info "Bilgi"

    Makinenizde az sayıda CPU varsa ve süreç başına yüksek bir sayı tahsis ederseniz, Nextflow mevcut olandan daha fazla CPU talep etmeyeceğinden görev çağrıları birbirinin arkasında sıraya girebilir.

Önce ve sonrayı karşılaştırabilmek için farklı bir rapor dosya adıyla tekrar çalıştırın.

```bash
nextflow run main.nf -with-report report-config-2.html
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [voluminous_venter] revision: c3c85dec78

    executor >  local (8)
    [a1/0e96d4] sayHello (1)       | 3 of 3 ✔
    [a3/7173a3] convertToUpper (2) | 3 of 3 ✔
    [4f/a8ae3d] collectGreetings   | 1 of 1 ✔
    [91/3724f8] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

`cowpy` için iki raporu karşılaştırdığımızda:

| Metrik                | Önce (1 CPU, 1 GB) | Sonra (2 CPU, 2 GB) |
| --------------------- | ------------------- | -------------------- |
| Kullanılan peak bellek | 6.4 MB              | 6.4 MB               |
| CPU kullanımı         | 116%                | 118%                 |

Tahsisi iki katına çıkarmak gerçek kullanımı hiç değiştirmedi; bu da orijinal 1 GB / 1 CPU'nun bu basit iş yükü için zaten yeterli olduğunu gösteriyor.
Önemsiz olmayan veriler işleyen gerçek bir pipeline'da, sayıların süreçler arasında anlamlı biçimde farklılaşmasını beklersiniz; bu da tahsis etmek yerine tahmin yapmak yerine önce profil çıkarmanız gerektiğinin tam olarak nedenidir.

### 1.3. Kaynak Sınırları Ekleme

İşlem altyapınıza bağlı olarak, örneğin küme genelinde bir üst sınır gibi, talep edebilecekleriniz üzerinde katı kısıtlamalar olabilir.
`resourceLimits` yönergesi bu sınırları belirlemenizi sağlar:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow bunları hedef yürütücünün beklediği biçime çevirir.
Bir süreç sınırdan fazlasını talep ederse, talep reddedilmek yerine sınıra indirilir.

!!! warning "Uyarı"

    Bu, eğitim ortamında çalıştırabileceğiniz bir şey değildir; etkili olabilmesi için HPC altyapısı gerektirir.

??? info "Kurumsal referans yapılandırmaları"

    nf-core projesi, geniş bir HPC ve bulut yürütücüsü yelpazesini kapsayan, dünya genelindeki kurumlar tarafından paylaşılan bir [yapılandırma dosyaları koleksiyonu](https://nf-co.re/configs/) yönetmektedir.
    Kendi kurumunuz aralarında olsun ya da olmasın, bunlar faydalı bir başlangıç noktasıdır.

### Özetle

Kaynak kullanımını değerlendirmek için profil raporu oluşturmayı, belirli bir süreç için kaynak tahsislerini geçersiz kılmayı ve `resourceLimits` ile tahsisleri sınırlamayı öğrendiniz.

### Sırada ne var?

Kaynak tahsis tahmininiz doğru olsun ya da olmasın, bir görev başarısız olduğunda pipeline'ın otomatik olarak nasıl kurtarılacağını öğrenin.

---

## 2. Yeniden Denemelerle Görev Hatalarını Yönetme

Profil çıkarma, bir sürecin çoğu zaman neye ihtiyaç duyduğunu söyler; ancak gerçek iş yükleri değişkendir: çoğu girdi için yeterli olan bir tahsis, alışılmadık derecede büyük bir girdi için yetersiz kalabilir ve tahminler basitçe yanlış olabilir.
Tek bir başarısız görevin tüm çalışmayı durdurmasına izin vermek yerine, Nextflow başarısız bir görevi otomatik olarak yeniden deneyebilir ve isteğe bağlı olarak her denemede daha fazla kaynak sağlayabilir.

### 2.1. Başarısız Bir Görevi Otomatik Olarak Yeniden Deneme

Bunu uygulamada görmek için, `cowpy`'nin bellek tahsisini gerçekte ihtiyaç duyduğunun altına kasıtlı olarak ayarlayın: [1.1](#11-generate-a-resource-utilization-report) bölümünden peak değerinin yaklaşık 6.4 MB olduğunu hatırlayın, dolayısıyla 6 MB tam olarak yetmeyecektir.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5-7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 2.GB
            cpus = 2
        }
    }
    ```

`errorStrategy`, bir görev başarısız olduğunda Nextflow'un ne yapacağını belirtir: `'retry'`, tüm pipeline'ı durdurmak yerine görevi yeniden gönderir.
`maxRetries`, Nextflow vazgeçmeden önce kaç ek deneme yapılacağını sınırlar.

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı (kısaltılmış)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [desperate_brazil] revision: c3c85dec78

    executor >  local (10)
    [67/fe1f49] sayHello (1)       | 3 of 3 ✔
    [8a/f13335] convertToUpper (1) | 3 of 3 ✔
    [39/1b24ed] collectGreetings   | 1 of 1 ✔
    [7a/d5eb6f] cowpy              | 0 of 1, retries: 2 ✘
    [9d/b79eb3] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)
    [6d/1d9d84] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (2)
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (137)

    Command executed:
      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      /usr/local/bin/_activate_current_env.sh: line 35:    14 Killed                  micromamba activate "${ENV_NAME:-base}"

    Work dir:
      /workspaces/training/execution-config/work/7a/d5eb6feeac0eed18d95d3da7a7aeb4

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

137 çıkış kodu, yetersiz bellek nedeniyle sonlandırma için standart sinyaldir: konteynerin `cowpy`'yi çalıştırmaya yetecek kadar belleği yoktu.
Nextflow görevi iki kez yeniden denedi; `maxRetries = 2` ile eşleşen toplam üç deneme yapıldı.
Denemeler arasında bellek tahsisi hiç değişmediğinden, her deneme aynı sorunla karşılaştı; yeniden denemeler tükenince Nextflow hatayı tam olarak raporlar ve pipeline'ı sıfır olmayan bir durum koduyla durdurur.

Temel neden denemeler arasında değişmiyorsa, yeniden denemenin tek başına hiçbir şeyi düzeltmediğini unutmayın.

### 2.2. Her Yeniden Denemede Kaynakları Artırma

Bir süreç yönergesi içinde `task.attempt`, 1'den başlayan mevcut deneme numarasını tutar.
Her yeniden denemede bir kaynak tahsisini artırmak için bunu bir closure içinde kullanabilirsiniz.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="6" hl_lines="5 7"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = { 6.MB * task.attempt }
            errorStrategy = 'retry'
            maxRetries = 3
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
            memory = 6.MB
            errorStrategy = 'retry'
            maxRetries = 2
        }
    }
    ```

İş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı (kısaltılmış)"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [grave_joliot] revision: c3c85dec78

    executor >  local (9)
    [0f/211b8a] sayHello (2)       | 3 of 3 ✔
    [26/301ea2] convertToUpper (3) | 3 of 3 ✔
    [22/5a895b] collectGreetings   | 1 of 1 ✔
    [e1/beee86] cowpy              | 1 of 1, retries: 1 ✔
    [b6/7aed6a] NOTE: Process `cowpy` terminated with an error exit status (137) -- Execution is retried (1)

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hola-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

İlk deneme 6 MB'ta yine başarısız olur; ancak yeniden deneme 12 MB (`6.MB * 2`) ile çalışır ve başarılı olur; pipeline tüm çıktılar yayımlanarak tamamlanır.

!!! warning "Uyarı"

    Konsol çıktısı, pipeline bir bütün olarak başarılı olsa bile başarısız ilk denemeyi bildiren bir `NOTE:` satırı içermeye devam eder: Nextflow her yeniden denemeyi ayrı ayrı kaydeder; ancak yeniden denenen bir hata genel sonucu etkilemez.
    Çalışmanın gerçekten başarılı olup olmadığını doğrulamak için `Outputs:` özetini veya komutun çıkış durumunu kontrol edin.

Belirli bir hatanın hangisi olduğuna göre ölçeklendirme dahil daha gelişmiş yeniden deneme kalıpları için Nextflow belgelerindeki [Dinamik işlem kaynakları](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) sayfasına bakın.

### Özetle

Bir pipeline'ın başarısız görevleri otomatik olarak nasıl yeniden deneyeceğini ve `task.attempt` kullanarak her yeniden denemede kaynak tahsislerinin nasıl artırılacağını öğrendiniz.

### Sırada ne var?

Bu tür yapılandırmaları değiştirilebilir profillere nasıl paketleyeceğinizi öğreneceğiniz [Bölüm 3](./03_profiles.md)'e geçin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Kaynak kullanımını değerlendirmek için profil raporu oluşturma ve süreç başına kaynak tahsisi belirleme
- `resourceLimits` ile kaynak taleplerini sınırlama
- `errorStrategy` ve `maxRetries` ile başarısız bir görevi otomatik olarak yeniden deneme
- `task.attempt` kullanarak her yeniden denemede kaynak tahsisini artırma
