# Bölüm 2: Pipeline'ı Yapılandırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[Bölüm 1](./01_run_nextflow.md)'de, konteynerler kullanarak birden fazla girdiyi paralel olarak işleyen eksiksiz bir çok adımlı pipeline çalıştırdınız.
Şimdi `nextflow.config` kullanarak pipeline davranışının nasıl yapılandırılacağına bakacağız: önce size sağladığımız yapılandırma dosyasını inceleyerek, ardından yapılandırma sağlamanın birkaç farklı yolunu keşfederek ve son olarak çıktıların nasıl ve nereye yayımlanacağını kontrol ederek.

---

## 1. Ana Yapılandırma Dosyasını İnceleme

Nextflow, `nextflow.config` dosyasını çalışma dizininden otomatik olarak alır ve ayarlarını her çalıştırmaya uygular.

Dört alanı kapsayan bir yapılandırma dosyası sağlıyoruz: yazılım paketleme, süreç ayarları, pipeline parametreleri ve yürütme profilleri.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Yazılım paketleme
     */
    docker.enabled = true

    /*
     * Süreç ayarları
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Pipeline parametreleri
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Profiller
     */
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

Her birini tek tek inceleyelim, ardından pipeline'ı bir profille çalıştırarak profilleri pratikte kullanalım.

!!! note "Not"

    Bu yapılandırma, tek bir makinede yerel yürütmeyi kapsar.
    Nextflow aynı zamanda HPC zamanlayıcılarını (SLURM, PBS, LSF) ve bulut yürütücülerini (AWS Batch, Google Cloud Batch, Azure Batch) de destekler; bunların tümü aynı `nextflow.config` mekanizmasıyla yapılandırılır.
    Bu seçeneklerin tam bir açıklaması için [Configure Execution](../config_exec/01_packaging_and_execution.md) kursundaki [Bölüm 1: Hesaplama Ortamınıza Uyum Sağlama](../config_exec/index.md) sayfasına bakın.

### 1.1. Yazılım Paketleme

Yazılım paketleme, Nextflow'un süreçlerinizin ihtiyaç duyduğu araçları nasıl sağladığını belirler; bu bir konteyner imajı, bir Conda ortamı veya başka bir şey olabilir.

```groovy title="nextflow.config" linenums="1"
/*
 * Yazılım paketleme
 */
docker.enabled = true
```

Bu satır, her süreç için Docker'ı etkinleştirir.
`container` yönergesi tanımlayan her süreç, belirtilen imajın içinde çalışır.

### 1.2. Süreç Ayarları

Bir sürecin, pipeline'ınızdaki `sayHello` veya `cowpy` gibi tek bir adım olduğunu hatırlayın.
Nextflow, her birinin gerçekte nasıl çalışacağına dair pek çok şeyi yapılandırmanıza olanak tanır: ne kadar CPU ve bellek alacağı, hangi konteyneri veya Conda ortamını kullanacağı ve daha fazlası.

```groovy title="nextflow.config" linenums="6"
/*
 * Süreç ayarları
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Bu, her süreci tek bir CPU ve 1 GB bellekle sınırlar.

Nextflow ayrıca tek tek adlandırılmış süreçler veya süreç grupları için farklı değerler belirlemenize de olanak tanır; bunu nasıl yapacağınızı [Configure Execution](../config_exec/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) kursunun [Bölüm 2: Hesaplama Kaynaklarını ve Hataları Yönetme](../config_exec/index.md) bölümünde öğreneceksiniz.

### 1.3. Pipeline Parametreleri

Parametreler, pipeline'ın komut satırı girdileridir; daha önce doğrudan komut satırında belirlediğiniz `--input`, `--batch` ve `--character` bayraklarıdır.
Bunlar için burada varsayılan değerler belirlemek, her seferinde bunları yazmak zorunda kalmamanızı sağlar; ancak bu bölümün ilerleyen kısımlarında göreceğiniz gibi, bunları sağlamanın birkaç farklı yolu daha vardır.

```groovy title="nextflow.config" linenums="14"
/*
 * Pipeline parametreleri
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Bu varsayılanlar, komut satırında bir parametre sağlanmadığında devreye girer; dolayısıyla `nextflow run main.nf` komutunu hiçbir bayrak olmadan çalıştırmak yine de işe yarar.

### 1.4. Profiller

Profiller, bir dizi ayarı tek bir ad altında gruplandırmanıza olanak tanır; böylece her seferinde değerleri elle değiştirmek yerine tek bir bayrakla tüm yapılandırmalar arasında geçiş yapabilirsiniz.

```groovy title="nextflow.config" linenums="23"
/*
 * Profiller
 */
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

`test` profili, pipeline'ı küçük ve iyi tanımlanmış bir girdi setiyle çalıştırmak için üç parametreyi geçersiz kılar; her nf-core pipeline'ı hızlı doğrulama için böyle bir profille birlikte gelir ve bu, kendi pipeline'larınızda da uygulamaya değer bir kuraldır.

`conda` profili, yazılım paketlemeyi Docker'dan Conda'ya geçirir.

Bir profili etkinleştirmek için komut satırında `-profile <isim>` parametresini kullanırsınız.

`test` profilini pratikte kullanalım.

```bash
nextflow run main.nf -profile test
```

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

Pipeline, `batch = 'test'` ve `character = 'tux'` değerleriyle çalışır.
`results/test/` dizinine bakın: toplu iş adı artık dizin yolunun kendisinin bir parçasıdır ve ASCII sanatı hindi yerine tux penguenini göstermektedir.

!!! note "Not"

    Aynı anda birden fazla profili etkinleştirebilir ve herhangi bir şey çalıştırmadan önce tam olarak çözümlenmiş sonucu görmek için `nextflow config -profile <isim>,<isim>` komutunu kullanabilirsiniz.
    Profillerin birleştirilmesi ve Nextflow'un aralarındaki çakışmaları nasıl çözdüğü, [Configure Execution](../config_exec/03_profiles.md) kursunun [Bölüm 3: Yapılandırmaları Değiştirmek için Profilleri Kullanma](../config_exec/index.md) bölümünde ayrıntılı olarak ele alınmaktadır.

### Özetle

Bir `nextflow.config` dosyasının en yaygın öğelerinin ne işe yaradığını ve bir profilin nasıl etkinleştirileceğini öğrendiniz.

### Sırada ne var?

Ana `nextflow.config` dosyasını değiştirmeden yapılandırma değerleri sağlamanın birkaç farklı yolunu öğrenin; bu yöntemler bireysel çalıştırmaları yapılandırmak ve tam bir ayar setini başkasıyla paylaşmak için kullanışlıdır.

---

## 2. Ek Dosyalar Aracılığıyla Yapılandırma Sağlama

`nextflow.config` dosyasında varsayılan değerleri belirlemek, nadiren değişen değerler için iyi çalışır.
Nextflow ayrıca iki daha hedefli mekanizma sunar: yürütmeyi belirli bir ortama uyarlamak için çalıştırmaya özgü bir yapılandırma dosyası ve bir iş ortağıyla tam bir girdi değerleri setini paylaşmak için bir parametre dosyası.

### 2.1. Çalıştırmaya Özgü Bir Yapılandırma Dosyası Kullanma

Diyelim ki pipeline'ı Docker'ın bulunmadığı bir makineye taşıyorsunuz ve her sürece daha fazla çalışma alanı vermek istiyorsunuz.
Yalnızca ihtiyacınız olan geçersiz kılmaları içeren yeni bir yapılandırma dosyası oluşturun:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Ana pipeline'ınızla birlikte `-c` parametresiyle geçirin:

```bash
nextflow run main.nf -c custom.config
```

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow, `custom.config` dosyasını pipeline'ın kendi `nextflow.config` dosyasının üzerine birleştirir; böylece her süreç artık varsayılanlar yerine 2 CPU ve 2 GB bellek alır ve Docker yerine Conda üzerinden çalışır.
`cowpy`, konteynerinin yanı sıra bir Conda paketi de tanımlanmış tek süreçtir; dolayısıyla Nextflow'un gerçekten bir ortam oluşturduğunu göreceğiniz süreç bu olacaktır:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Yalnızca kaynak tahsisini ve paketlemeyi geçersiz kılan, pipeline parametrelerine dokunmayan küçük bir dosya, nf-core pipeline'larının kurumsal yapılandırmalardan beklediği tam da bu kalıptır.
Gerçek dünya örnekleri için [nf-core/configs](https://github.com/nf-core/configs) deposuna göz atın.

Bu, normal yapılandırmanıza dokunmadan bir pipeline'ı yeni bir ortama uyarlamanın tek kullanımlık bir yolunu sunar.

### 2.2. Parametre Dosyası Kullanma

Diyelim ki bir iş ortağıyla tam bir çalıştırma parametreleri setini paylaşmanız veya bunları bir yayın için kaydetmeniz gerekiyor.

Nextflow, YAML veya JSON formatında [parametre dosyaları](https://nextflow.io/docs/latest/config.html#parameter-file) sağlamanıza olanak tanır; bunlar tam ve tekrarlanabilir bir değer setini dağıtmanın daha basit bir yoludur.

`test-params.yaml` adlı bir parametre dosyası çalışma dizininizde zaten mevcuttur:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Bu dosya düz YAML olduğundan, `nextflow.config` dosyasında kullanılan eşittir işaretleri (`=`) yerine iki nokta üst üste (`:`) kullanır.

!!! info "Bilgi"

    JSON sürümü olan `test-params.json` da sağlanmıştır. Kendiniz deneyebilirsiniz; geçirme sözdizimi aynıdır.

Dosyayı `-params-file` ile geçirin:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Dosya içeriği"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Bir pipeline'ın birden fazla parametresi olduğunda parametre dosyası özellikle değerlidir: iş akışı betiğinde herhangi bir değişiklik yapmadan veya uzun bir komut satırı oluşturmadan tümünü tek seferde sağlamanıza olanak tanır ve sonuçlarınızla birlikte dağıtmak da kolaydır.

### Özetle

Yapılandırma sağlamanın iki farklı yolunu daha öğrendiniz: yürütmeyi yeni bir ortama uyarlamak için çalıştırmaya özgü bir yapılandırma dosyası ve tam, tekrarlanabilir girdi değerlerini paylaşmak için bir parametre dosyası.

### Sırada ne var?

Pipeline'ınızın çıktılarının nasıl ve nereye yayımlanacağını nasıl kontrol edeceğinizi öğrenin.

---

## 3. Pipeline Çıktılarını Yönetme

Bir pipeline yazarı, çıktıların kodda nasıl düzenleneceğine karar verir; ancak nereye gideceklerini veya nasıl aktarılacaklarını kontrol etmek için bu koda dokunmanız gerekmez.
Nextflow bunun yerine yapılandırma düzeyinde yollar sunar: temel bir çıktı dizini belirleyin ve dosyaların kopyalanıp kopyalanmayacağını ya da sembolik bağlantı olarak mı oluşturulacağını seçin.

### 3.1. Çıktı Dizinini Özelleştirme

Nextflow varsayılan olarak çıktıları `results/` altında yayımlar.
`-output-dir` (veya kısa formu `-o`) ile farklı bir konuma yönlendirin:

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Dizin içeriği"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

Çıktılar artık yerleşik `results/batch/` varsayılanı yerine `outputs/batch/` altına yerleşir.
Pipeline'ın kendi kodu, `batch/` ve `intermediates/` alt dizinleri gibi bu temel dizin içindeki yapıyı belirlemeye devam eder; `-output-dir` yalnızca bu yapının nerede başlayacağını kontrol eder.

`-output-dir`, aslında `outputDir` yapılandırma seçeneği için yalnızca bir komut satırı kısayoludur; dolayısıyla yapılandırmanın gidebileceği her yere yerleştirilebilir: doğrudan `nextflow.config` içine, bir profil içine veya bu bölümün başında kullandığınız gibi bir `-c` katman dosyasına.
Örneğin, aşağıdaki kod parçacığı aynı ayarın komut satırında geçirilmek yerine doğrudan `nextflow.config` içine yerleştirildiğini göstermektedir:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Bu tür bir yapılandırma seçeneğinin bulunabileceği tam liste için Nextflow referansındaki [Yapılandırma dosyası](https://nextflow.io/docs/latest/config.html) sayfasına bakın.

### 3.2. Çıktıların Nasıl Yayımlanacağını Seçme

Nextflow varsayılan olarak çıktıları, gerçek kopyalar olarak değil, `work/` altındaki çıktı konumlarına işaret eden sembolik bağlantılar olarak yayımlar:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Pipeline yazarları, iş akışı kodundaki her bir süreç için 'yayımlama modunu' `'copy'` veya `'move'` olarak ayarlayabilir.
Bunu genellikle pipeline'ın son çıktıları için yaparlar; tam pipeline çalıştırıldıktan sonra silinebilecek ara dosyalar için ise varsayılan `'symlink'` davranışını bırakırlar.

Bu yaklaşım diskte veri çoğaltmayı önler; ancak `-resume` kullanma yeteneğini kaybetmeden `work/` altındaki görev dizinlerini silemeyeceğiniz anlamına gelir.
Tüm çıktı dosyalarının düzgün şekilde kopyalanmasını istiyorsanız, pipeline yapılandırmanızda [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) değerini `'copy'` olarak ayarlayın. (`-output-dir`'den farklı olarak bunun için komut satırı bayrağı yoktur; yalnızca yapılandırma dosyasıyla ayarlanabilir.)

`nextflow.config` dosyasında ayarlamayı deneyin:

=== "Sonra"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Önce"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Ardından pipeline'ı çalıştırın; çıktılardaki farkı görebilmek için toplu iş adını değiştirin:

```bash
nextflow run main.nf --batch withmode
```

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Daha önce yaptığınız gibi çıktı dosyalarından birine bakın:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Artık `work/` dizini temizlense bile erişilebilir kalacak gerçek, bağımsız bir dosyadır.

!!! warning "Uyarı"

    `workflow.output.mode` ayarı yalnızca pipeline kodunda henüz bir mod belirlenmemiş çıktılar için varsayılan değeri doldurur.
    Ne ayarlarsanız ayarlayın, yazarın sabit olarak kodladığı bir modu geçersiz kılamaz.

### Özetle

Pipeline koduna dokunmadan temel çıktı dizinini nasıl özelleştireceğinizi ve kopyalanmış ile sembolik bağlantılı çıktılar arasında nasıl seçim yapacağınızı öğrendiniz.

### Sırada ne var?

Geçmiş çalıştırmaların geçmişini nasıl inceleyeceğinizi, yürütme raporlarını nasıl oluşturacağınızı ve eski çalışma dizinlerini nasıl temizleyeceğinizi öğreneceğiniz [Bölüm 3](./03_manage_executions.md)'e geçin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- `nextflow.config` ve profiller kullanarak pipeline davranışını yapılandırma
- Çalıştırmaya özgü bir yapılandırma dosyası veya parametre dosyası aracılığıyla yapılandırma sağlama
- Çıktı dizinini özelleştirme ve kopyalanmış ile sembolik bağlantılı çıktılar arasında seçim yapma
