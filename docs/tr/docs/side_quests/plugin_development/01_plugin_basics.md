# Bölüm 1: Eklenti Temelleri

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu bölümde, eklentilerin Nextflow'u nasıl genişlettiğini öğrenecek, ardından üç farklı eklentiyi çalışırken görmek için deneyeceksiniz.

---

## 1. Eklentiler nasıl çalışır?

Eklentiler, Nextflow'u çeşitli uzantı türleri aracılığıyla genişletir:

| Uzantı Türü          | Ne yapar                                             | Örnek                            |
| -------------------- | ---------------------------------------------------- | -------------------------------- |
| Fonksiyonlar         | İş akışlarından çağrılabilir özel fonksiyonlar ekler | `samplesheetToList()`            |
| İş akışı monitörleri | Görev tamamlama gibi olaylara yanıt verir            | Özel günlükleme, Slack uyarıları |
| Yürütücüler          | Görev yürütme arka uçları ekler                      | AWS Batch, Kubernetes            |
| Dosya sistemleri     | Depolama arka uçları ekler                           | S3, Azure Blob                   |

Fonksiyonlar ve iş akışı monitörleri (Nextflow API'sinde "trace observer" olarak adlandırılır), eklenti yazarları için en yaygın türlerdir.
Yürütücüler ve dosya sistemleri genellikle platform satıcıları tarafından oluşturulur.

Sonraki alıştırmalar, her iki türü de çalışırken görebilmeniz için fonksiyon eklentilerini ve bir observer eklentisini göstermektedir.

---

## 2. Fonksiyon eklentilerini kullanma

Fonksiyon eklentileri, iş akışlarınıza `include` ile içe aktardığınız çağrılabilir fonksiyonlar ekler.
İki tanesini deneyeceksiniz: nf-hello (basit bir örnek) ve nf-schema (yaygın olarak kullanılan gerçek dünya eklentisi).
Her iki alıştırma da aynı `hello.nf` pipeline'ını değiştirir; bu sayede eklentilerin mevcut bir iş akışını nasıl geliştirdiğini görebilirsiniz.

### 2.1. nf-hello: elle yazılmış kodu değiştirme

[nf-hello](https://github.com/nextflow-io/nf-hello) eklentisi, rastgele dizeler üreten bir `randomString` fonksiyonu sağlar.
Pipeline, bu fonksiyonun kendi satır içi sürümünü zaten tanımlamaktadır; bunu eklentideki sürümle değiştireceksiniz.

#### 2.1.1. Başlangıç noktasını inceleme

Pipeline'a bakın:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Rastgele bir alfanümerik dize üret
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

Pipeline, kendi `randomString` fonksiyonunu satır içinde tanımlar ve ardından her selamlamaya rastgele bir kimlik eklemek için kullanır.

Çalıştırın:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

Çıktı sıranız ve rastgele dizeleriniz farklı olacaktır; betiği tekrar çalıştırırsanız farklı bir rastgele selamlama kümesi elde edersiniz.

#### 2.1.2. Eklentiyi yapılandırma

Satır içi fonksiyonu eklentideki bir fonksiyonla değiştirin. `nextflow.config` dosyanıza şunu ekleyin:

```groovy title="nextflow.config"
// Eklenti geliştirme alıştırmaları için yapılandırma
plugins {
    id 'nf-hello@0.5.0'
}
```

Eklentiler, `nextflow.config` dosyasında `plugins {}` bloğu kullanılarak tanımlanır.
Nextflow, bunları topluluk ve resmi eklentilerin merkezi deposu olan [Nextflow Plugin Registry](https://registry.nextflow.io/)'den otomatik olarak indirir.

#### 2.1.3. Eklenti fonksiyonunu kullanma

Satır içi `randomString` fonksiyonunu eklenti sürümüyle değiştirin:

=== "Sonra"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Önce"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Rastgele bir alfanümerik dize üret
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

`include` ifadesi, `randomString` fonksiyonunu daha geniş bir katkıda bulunan havuzu tarafından kanıtlanmış, test edilmiş ve bakımı yapılmış bir kütüphaneden içe aktarır; bu kişiler hataları tespit edip düzeltebilir.
Her pipeline kendi fonksiyon kopyasını korumak yerine, eklentiyi kullanan her pipeline aynı denetlenmiş uygulamayı alır.
Bu durum, yinelenen kodu ve beraberinde gelen bakım yükünü azaltır.
`#!groovy include { function } from 'plugin/plugin-id'` sözdizimi, Nextflow modülleri için kullanılan `include` ifadesinin aynısıdır; yalnızca `plugin/` öneki eklenir.
[`randomString` için kaynak kodunu](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) GitHub'daki nf-hello deposunda inceleyebilirsiniz.

#### 2.1.4. Çalıştırma

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Rastgele dizeleriniz farklı olacaktır.)

Çıktıda hâlâ rastgele sonekler bulunmaktadır; ancak artık `randomString` satır içi koddan değil, nf-hello eklentisinden gelmektedir.
"Pipeline is starting!" ve "Pipeline complete!" mesajları yenidir.
Bunlar, Bölüm 5'te inceleyeceğiniz eklentinin observer bileşeninden gelmektedir.

Nextflow, eklentileri ilk kullanıldıklarında otomatik olarak indirir; bu nedenle `nf-hello@0.5.0` tanımlayan her pipeline, projeler arasında kod kopyalamadan aynı test edilmiş `randomString` fonksiyonunu alır.

Artık bir fonksiyon eklentisi kullanmanın üç adımını gördünüz: `nextflow.config` dosyasında tanımlayın, fonksiyonu `include` ile içe aktarın ve iş akışınızda çağırın.
Bir sonraki alıştırma, bu adımların aynısını gerçek dünya eklentisine uygular.

### 2.2. nf-schema: doğrulamalı CSV ayrıştırma

[nf-schema](https://github.com/nextflow-io/nf-schema) eklentisi, en yaygın kullanılan Nextflow eklentilerinden biridir.
Beklenen sütunları ve türleri tanımlayan bir JSON şeması kullanarak CSV/TSV dosyalarını ayrıştıran `samplesheetToList` fonksiyonunu sağlar.

Pipeline şu anda `greetings.csv` dosyasını `splitCsv` ve manuel bir `map` kullanarak okumaktadır; ancak nf-schema bunu doğrulanmış, şema tabanlı ayrıştırmayla değiştirebilir.
Alıştırma dizininde bir JSON şema dosyası (`greetings_schema.json`) zaten sağlanmıştır.

??? info "Şema nedir?"

    Şema, geçerli verinin nasıl göründüğünün biçimsel bir tanımıdır.
    Hangi sütunların beklendiğini, her değerin hangi türde olması gerektiğini (string, number vb.) ve hangi alanların zorunlu olduğunu tanımlar.

    Bunu bir sözleşme olarak düşünebilirsiniz: girdi verisi şemayla eşleşmezse, araç sorunu pipeline'da ilerleyen aşamalarda kafa karıştırıcı hatalara yol açmak yerine erken aşamada yakalayabilir.

#### 2.2.1. Şemayı inceleme

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

Şema iki sütun tanımlar (`greeting` ve `language`) ve `greeting` sütununu zorunlu olarak işaretler.
Birisi `greeting` sütunu eksik bir CSV geçirirse, nf-schema pipeline çalışmadan önce hatayı yakalar.

#### 2.2.2. nf-schema'yı yapılandırmaya ekleme

`nextflow.config` dosyasını her iki eklentiyi de içerecek şekilde güncelleyin:

=== "Sonra"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. hello.nf dosyasını samplesheetToList kullanacak şekilde güncelleme

`splitCsv` girdisini `samplesheetToList` ile değiştirin:

=== "Sonra"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Önce"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Özel `splitCsv` ve `map` ayrıştırma kodu, pipeline çalışmadan önce örnek sayfayı şemaya göre doğrulayan kanıtlanmış ve test edilmiş bir fonksiyon olan `samplesheetToList` ile değiştirilir.
Bu durum, elle yazılmış ayrıştırma mantığının bakım yükünü azaltırken pipeline kullanıcılarının deneyimini de iyileştirir; kullanıcılar girdileri beklenen biçimle eşleşmediğinde açık hata mesajları alır.
Her satır, sütun sırasına göre bir değerler listesi hâline gelir; dolayısıyla `row[0]` selamlama, `row[1]` ise dildir.

#### 2.2.4. Çalıştırma

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Rastgele dizeleriniz farklı olacaktır.)

Çıktı aynıdır; ancak artık şema, pipeline çalışmadan önce CSV yapısını doğrulamaktadır.
Karmaşık örnek sayfaları ve çok sayıda sütun içeren gerçek pipeline'larda bu tür doğrulama, manuel `splitCsv` + `map` kullanımının gözden kaçıracağı hataları önler.

#### 2.2.5. Doğrulamayı çalışırken görme

Şema doğrulamasının neler yakaladığını görmek için `greetings.csv` dosyasına hatalar ekleyin.

Zorunlu `greeting` sütununu `message` olarak yeniden adlandırın:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Pipeline'ı çalıştırın:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

Şema bir `greeting` sütunu gerektirdiğinden ve bunu bulamadığından pipeline çalışmayı reddeder.

Şimdi zorunlu sütunu geri yükleyin; ancak isteğe bağlı `language` sütununu `lang` olarak yeniden adlandırın:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Bu sefer pipeline çalışır; ancak bir uyarı yazdırır:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Zorunlu sütunlar kesin hatalara yol açar; isteğe bağlı sütunlar ise uyarılara neden olur.
Bu tür erken geri bildirim, onlarca sütun içeren gerçek pipeline'larda hata ayıklama süresini kısaltır.

#### 2.2.6. Doğrulama davranışını yapılandırma

`lang` hakkındaki uyarı yararlıdır; ancak önem derecesini yapılandırma aracılığıyla kontrol edebilirsiniz.
Eklentiler, davranışlarını kontrol eden kendi yapılandırma kapsamlarını içerebilir.
nf-schema eklentisi `validation` yapılandırma kapsamını içerir; buradaki ayarları değiştirerek nf-schema'nın nasıl davranacağını değiştirebilirsiniz.

Tanınmayan başlıkların uyarı yerine hata oluşturmasını sağlamak için `nextflow.config` dosyasına bir `validation` bloğu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Aynı `lang` sütunu yerindeyken pipeline'ı tekrar çalıştırın:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

Pipeline artık uyarı vermek yerine başarısız olur.
Pipeline kodu değişmedi; yalnızca yapılandırma değişti.

Devam etmeden önce `greetings.csv` dosyasını özgün hâline geri yükleyin ve `validation` bloğunu kaldırın:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

Hem nf-hello hem de nf-schema, `include` ile içe aktardığınız ve iş akışı kodunuzda çağırdığınız fonksiyonlar sağlayan fonksiyon eklentileridir.
Bir sonraki alıştırma, hiçbir `include` ifadesi gerektirmeden çalışan farklı bir eklenti türünü göstermektedir.

---

## 3. Bir observer eklentisi kullanma: nf-co2footprint

Tüm eklentiler içe aktarılacak fonksiyonlar sağlamaz.
[nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) eklentisi, pipeline'ınızın kaynak kullanımını izlemek ve karbon ayak izini tahmin etmek için bir **trace observer** kullanır.
Herhangi bir pipeline kodunu değiştirmenize gerek yoktur; yalnızca yapılandırmaya ekleyin.

### 3.1. nf-co2footprint'i yapılandırmaya ekleme

`nextflow.config` dosyasını güncelleyin:

=== "Sonra"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Pipeline'ı çalıştırma

```bash
nextflow run hello.nf
```

Eklenti, yürütme sırasında çeşitli INFO ve WARN mesajları üretir.
Bunlar, yerel bir makinede çalışan küçük bir örnek için normaldir:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Bölge, yürütücü, CPU modeli ve bellek hakkındaki uyarılar, eklentinin yerel bir eğitim ortamının tam donanım ayrıntılarını algılayamamasından kaynaklanmaktadır.
Üretim ortamında (örneğin bir HPC kümesi veya bulut), bu değerler mevcut olacak ve tahminler daha doğru olacaktır.

Sonunda şuna benzer bir satır arayın:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Sayılarınız farklı olacaktır.)

### 3.3. Raporu görüntüleme

Eklenti, çalışma dizininizde çıktı dosyaları oluşturur:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Özete bakın:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(Sayılarınız farklı olacaktır.)

İlk bölüm, ham enerji ve emisyon rakamlarını göstermektedir.
"Which equals" bölümü, bu sayıları tanıdık eşdeğerlere dönüştürerek bağlama oturtmaktadır.
Özet ayrıca eklentinin yapılandırma seçeneklerini listeleyen bir bölüm ve hesaplama yönteminin dayandığı [Green Algorithms](https://doi.org/10.1002/advs.202100707) araştırma makalesine bir atıf içermektedir.

### 3.4. Eklentiyi yapılandırma

Bölüm 3.2'deki "Target zone null" uyarısı, eklentide yapılandırılmış bir konum olmadığı için görüntülendi.
nf-co2footprint eklentisi, coğrafi konumunuzu ayarlayabileceğiniz bir `co2footprint` yapılandırma kapsamı tanımlar.

`nextflow.config` dosyasına bir `co2footprint` bloğu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "İpucu"

    Tercih ederseniz kendi ülke kodunuzu kullanabilirsiniz (örneğin `'US'`, `'DE'`, `'FR'`).

Pipeline'ı çalıştırın:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

Bölge uyarısı artık görünmüyor.
Eklenti artık küresel yedek değer (480.0 gCO₂eq/kWh) yerine GB'ye özgü karbon yoğunluğunu (163.92 gCO₂eq/kWh) kullanmaktadır.

!!! note "Not"

    `WARN: Unrecognized config option 'co2footprint.location'` mesajını da görebilirsiniz.
    Bu görsel bir uyarıdır ve güvenle yoksayılabilir; eklenti değeri yine de doğru şekilde okur.

Bölüm 6'da kendi eklentiniz için bir yapılandırma kapsamı oluşturacaksınız.

Bu eklenti tamamen observer mekanizması aracılığıyla çalışır; iş akışı yaşam döngüsü olaylarına bağlanarak kaynak ölçümlerini toplar ve pipeline tamamlandığında raporunu oluşturur.

Artık fonksiyon eklentilerini (`include` ile içe aktarılan) ve bir observer eklentisini (yalnızca yapılandırma aracılığıyla etkinleştirilen) denediniz.
Bunlar en yaygın iki uzantı türüdür; ancak bölüm 1'deki tablonun gösterdiği gibi eklentiler yürütücüler ve dosya sistemleri de ekleyebilir.

---

## 4. Eklentileri keşfetme

[Nextflow Plugin Registry](https://registry.nextflow.io/), mevcut eklentileri bulmak için merkezi bir merkezdir.

![registry.nextflow.io'daki nf-hello eklenti sayfası](img/plugin-registry-nf-hello.png)

Her eklenti sayfası, açıklamasını, mevcut sürümlerini, kurulum talimatlarını ve belgelere bağlantıları gösterir.

---

## 5. Eklenti geliştirmeye hazırlanma

Aşağıdaki bölümler (Bölüm 2-6), nf-schema'ya dayanan ancak nf-hello veya nf-co2footprint'e dayanmayan ayrı bir pipeline dosyası olan `greet.nf` kullanır.

`nextflow.config` dosyasını yalnızca nf-schema'yı tutacak şekilde güncelleyin:

```groovy title="nextflow.config"
// Eklenti geliştirme alıştırmaları için yapılandırma
plugins {
    id 'nf-schema@2.6.1'
}
```

co2footprint çıktı dosyalarını kaldırın:

```bash
rm -f co2footprint_*
```

`hello.nf` dosyası, Bölüm 1 çalışmanızı referans olarak saklar; bundan sonra `greet.nf` ile çalışacaksınız.

---

## Özetle

Üç farklı eklenti kullandınız:

- **nf-hello**: `include` ile içe aktarılan `randomString` fonksiyonunu sağlayan bir fonksiyon eklentisi
- **nf-schema**: Şema doğrulamalı CSV ayrıştırması için `samplesheetToList` fonksiyonunu sağlayan bir fonksiyon eklentisi
- **nf-co2footprint**: Kaynak kullanımını otomatik olarak izleyen, `include` gerektirmeyen bir observer eklentisi

Temel kalıplar:

- Eklentiler, `nextflow.config` dosyasında `#!groovy plugins { id 'plugin-name@version' }` ile tanımlanır
- Fonksiyon eklentileri `#!groovy include { function } from 'plugin/plugin-id'` gerektirir
- Observer eklentileri, yapılandırmada tanımlandıktan sonra otomatik olarak çalışır
- Eklentiler, davranışı özelleştirmek için yapılandırma kapsamları tanımlayabilir (örneğin `#!groovy validation {}`, `#!groovy co2footprint {}`)
- [Nextflow Plugin Registry](https://registry.nextflow.io/), mevcut eklentileri listeler

---

## Sırada ne var?

Aşağıdaki bölümler, kendi eklentinizi nasıl oluşturacağınızı göstermektedir.
Eklenti geliştirmeyle ilgilenmiyorsanız burada durabilir veya [Özet](summary.md) bölümüne atlayabilirsiniz.

[Bölüm 2'ye devam edin :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
