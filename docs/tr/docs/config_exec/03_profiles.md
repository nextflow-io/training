# Bölüm 3: Yapılandırmaları değiştirmek için profiller kullanın

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


[Bölüm 1](./01_packaging_and_execution.md) ve [Bölüm 2](./02_resources_and_retries.md) boyunca birkaç yapılandırma seçeneği biriktirdiniz: yazılım paketleme, yürütme platformu ve kaynak tahsisleri.
Pratikte, çalıştırdığınız ortama bağlı olarak bu seçeneklerin tüm kümelerini değiştirmek isteyeceksiniz; örneğin geliştirme için dizüstü bilgisayar, üretim için HPC kümesi.

Nextflow, farklı yapılandırmaları tanımlayan istediğiniz sayıda [profil](https://nextflow.io/docs/latest/config.html#profiles) oluşturmanıza ve çalışma zamanında tek bir bayrakla bir (veya birkaç) tanesini seçmenize olanak tanır.

Zaten birini kullandınız: [Nextflow Run](../nextflow_run/index.md) bölümündeki `test` profili, girdi parametrelerini küçük ve iyi tanımlanmış bir kümeyle geçersiz kılar.
Şimdi kendi altyapı profillerinizi oluşturacak ve bunları mevcut profille birleştireceksiniz.

---

## 1. Farklı ortamlar için profiller oluşturun

### 1.1. Profilleri ayarlayın

`nextflow.config` dosyasına iki profil ekleyin: biri Docker ile normal bir dizüstü bilgisayarda çalıştırmak için, diğeri Slurm zamanlayıcısı ve Conda kullanan bir üniversite HPC kümesi için.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="35"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

`univ_hpc` profili aynı zamanda kaynak sınırlarını da belirler; bu, paylaşımlı HPC altyapısında genellikle zorunludur.

### 1.2. İş akışını bir profille çalıştırın

Çalışma zamanında `-profile` ile bir profil seçin.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Uyarı"

    `univ_hpc` profili, eğitim ortamında çalışmayacaktır; zira burada Slurm zamanlayıcısı mevcut değildir.

Her zaman birlikte kullanılması gereken başka ayarlar bulursanız, bunları ilgili profile ekleyin.
İhtiyaç duyduğunuz diğer kombinasyonları gruplamak için ek profiller de oluşturabilirsiniz.

### 1.3. Birden fazla profille çalıştırın

Profiller birbirini dışlamaz.
`-profile <profil1>,<profil2>` ile birden fazlasını aynı anda etkinleştirebilirsiniz.
`my_laptop` profilini, Nextflow Run'dan zaten tanıdığınız `test` profiliyle birleştirin.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/config-exec/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

Tek tek dosya adları, `test` profilinden `batch = 'test'` değerini doğru şekilde alır (`COLLECTED-test-output.txt` ve benzerleri).

Aynı seçeneği belirleyen profilleri birleştirdiğinizde, Nextflow çakışmayı en son okuduğu değeri, yani dosyada daha sonra gelen değeri kullanarak çözer.
Çakışan ayarlar tamamen farklı yapılandırma kaynaklarından geliyorsa, standart [öncelik sırası](https://www.nextflow.io/docs/latest/config.html) geçerli olur.

### Özetle

Altyapıya özgü yapılandırmaları bir araya getiren profilleri nasıl tanımlayacağınızı, çalışma zamanında `-profile` ile nasıl seçeceğinizi, tek bir çalıştırmada birden fazla profili nasıl birleştireceğinizi ve birden fazla profil aynı seçeneği belirlediğinde Nextflow'un çakışmaları nasıl çözdüğünü öğrendiniz.

### Sırada ne var?

Her şeyi çalıştırmadan önce tam olarak çözümlenmiş yapılandırmayı nasıl inceleyeceğinizi öğrenin.

---

## 2. Çözümlenmiş yapılandırmayı inceleyin

[Nextflow Run](../nextflow_run/02_configure_pipeline.md) bölümünde tek bir profilin neye çözümlendiğini kontrol etmek için `nextflow config -profile test` komutunu zaten kullandınız.
Bu komut, birden fazla profili birleştirdiğinizde özellikle kullanışlı hale gelir: az önce gördüğünüz gibi, iki profil aynı seçeneği belirlediğinde hangi değerin gerçekten geçerli olduğunu elle hesaplamak güçleşebilir.
`nextflow config` komutu, pipeline'ı çalıştırmadan tüm bunları sizin için çözer.

### 2.1. Varsayılan yapılandırmayı çözümleyin

```bash
nextflow config
```

??? success "Komut çıktısı"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

Bu, pipeline'ı ekstra bayrak olmadan çalıştırsaydınız geçerli olacak yapılandırmanın tam olarak kendisidir.

### 2.2. Profiller etkinleştirilmiş halde yapılandırmayı çözümleyin

Gerçek bir çalıştırmada kullanacağınız profillerin aynısını ekleyin.

```bash
nextflow config -profile my_laptop,test
```

??? success "Komut çıktısı"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

İkisini karşılaştırmak neyin değiştiğini doğrular: `params.batch`, `params.character` ve `process.executor` değerlerinin tümü `my_laptop,test` profillerini yansıtır.
Bu, çözümlenmiş ayarları elle hesaplamanın sıkıcı ve hataya açık olacağı, çok katmanlı yapılandırmaya sahip pipeline'lar için özellikle değerli hale gelir.

### Özetle

Herhangi bir profil kombinasyonu için tam olarak çözümlenmiş yapılandırmayı, her şeyi çalıştırmadan önce incelemek amacıyla `nextflow config` komutunu nasıl kullanacağınızı öğrendiniz.

### Sırada ne var?

Nextflow pipeline'larını yapılandırmanın temellerini tamamladınız.
Buradan nereye gideceğinizi öğrenmek için [Kurs özeti](next_steps.md) sayfasına bakın.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Altyapıya özgü yapılandırmaları bir araya getiren profiller tanımlamak
- Tek bir çalıştırmada birden fazla profili birleştirmek ve aralarındaki çakışmaların nasıl çözüldüğünü anlamak
- Tam olarak çözümlenmiş yapılandırmayı incelemek için `nextflow config` komutunu kullanmak
