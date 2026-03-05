# Temel Nextflow Kodlama Kalıpları

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, Java Sanal Makinesi üzerinde çalışan bir programlama dilidir. Nextflow, [Groovy](http://groovy-lang.org/) üzerine kurulu olsa ve sözdiziminin çoğunu onunla paylaşsa da, Nextflow "uzantılı Groovy"den fazlasıdır -- tam olarak belirtilmiş bir [sözdizimi](https://nextflow.io/docs/latest/reference/syntax.html) ve [standart kütüphane](https://nextflow.io/docs/latest/reference/stdlib.html) ile bağımsız bir dildir.

Değişkenler, map'ler ve listeler için temel sözdiziminin ötesine geçmeden çok sayıda Nextflow kodu yazabilirsiniz. Çoğu Nextflow eğitimi iş akışı düzenlemesine (channel'lar, process'ler ve veri akışı) odaklanır ve sadece bununla bile şaşırtıcı derecede ileri gidebilirsiniz.

Ancak, veri manipüle etmeniz, karmaşık dosya adlarını ayrıştırmanız, koşullu mantık uygulamanız veya sağlam üretim iş akışları oluşturmanız gerektiğinde, kodunuzun iki ayrı yönünü düşünmek yardımcı olur: **dataflow** (channel'lar, operatörler, process'ler ve workflow'lar) ve **scripting** (closure'lar, fonksiyonlar ve process script'leri içindeki kod). Bu ayrım bir nebze keyfi olsa da (hepsi Nextflow kodu), ne zaman pipeline'ınızı düzenleyip ne zaman veriyi manipüle ettiğinizi anlamak için yararlı bir zihinsel model sağlar. Her ikisinde de ustalaşmak, açık ve bakımı yapılabilir iş akışları yazma yeteneğinizi önemli ölçüde artırır.

### Öğrenme hedefleri

Bu yan görev, sizi temel kavramlardan üretime hazır kalıplara kadar uygulamalı bir yolculuğa çıkarır.
Basit bir CSV okuyan iş akışını karmaşık bir biyoinformatik pipeline'a dönüştüreceğiz; adım adım gerçekçi zorluklarla geliştireceğiz:

- **Sınırları anlamak:** Dataflow işlemleri ile scripting arasında ayrım yapmak ve birlikte nasıl çalıştıklarını anlamak
- **Veri manipülasyonu:** Güçlü operatörler kullanarak map'leri ve koleksiyonları çıkartmak, dönüştürmek ve alt kümelemek
- **String işleme:** Regex kalıpları ile karmaşık dosya adlandırma şemalarını ayrıştırmak ve değişken enterpolasyonunda ustalaşmak
- **Yeniden kullanılabilir fonksiyonlar:** Daha temiz, daha bakımı kolay iş akışları için karmaşık mantığı adlandırılmış fonksiyonlara çıkartmak
- **Dinamik mantık:** Farklı girdi türlerine uyum sağlayan process'ler oluşturmak ve dinamik kaynak tahsisi için closure'lar kullanmak
- **Koşullu yönlendirme:** Metadata özelliklerine göre örnekleri farklı process'lere akıllıca yönlendirmek
- **Güvenli işlemler:** Eksik verileri null-safe operatörlerle zarif bir şekilde ele almak ve girdileri net hata mesajlarıyla doğrulamak
- **Yapılandırma tabanlı işleyiciler:** Günlükleme, bildirimler ve yaşam döngüsü yönetimi için workflow olay işleyicilerini kullanmak

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmanız gerekir.
- Temel Nextflow kavramlarını ve mekanizmalarını (process'ler, channel'lar, operatörler, dosyalarla çalışma, metadata) rahatça kullanıyor olmanız gerekir.
- Yaygın programlama yapılarına (değişkenler, map'ler, listeler) temel düzeyde aşina olmanız gerekir.

Bu eğitim, programlama kavramlarını karşılaştıkça açıklayacaktır; bu nedenle kapsamlı programlama deneyimine ihtiyacınız yoktur.
Temel kavramlardan başlayıp ileri düzey kalıplara doğru ilerleyeceğiz.

---

## 0. Başlarken

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitimin dosyalarının bulunduğu dizine geçelim.

```bash
cd side-quests/essential_scripting_patterns
```

#### Materyalleri inceleyin

Bir ana iş akışı dosyası ve örnek veri dosyalarını içeren bir `data` dizini bulacaksınız.

```console title="Directory contents"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Örnek CSV dosyamız, özelliklerine göre farklı işleme gerektiren biyolojik örnekler hakkında bilgi içerir:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Bu gerçekçi veri setini kullanarak, gerçek biyoinformatik iş akışlarında karşılaşacağınız pratik programlama tekniklerini keşfedeceğiz.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım
<!-- - [ ] I understand the assignment -->

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Dataflow ve Scripting: Sınırları Anlamak

### 1.1. Neyin Ne Olduğunu Belirlemek

Nextflow iş akışları yazarken, **dataflow** (verinin channel'lar ve process'ler arasında nasıl hareket ettiği) ile **scripting** (veriyi manipüle eden ve kararlar veren kod) arasında ayrım yapmak önemlidir. Birlikte nasıl çalıştıklarını gösteren bir iş akışı oluşturalım.

#### 1.1.1. Temel Nextflow Workflow'u

CSV dosyasını okuyan basit bir iş akışı ile başlayın (bunu sizin için `main.nf` dosyasında zaten hazırladık):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` bloğu pipeline yapımızı tanımlarken, `channel.fromPath()` bir dosya yolundan channel oluşturur. `.splitCsv()` operatörü CSV dosyasını işler ve her satırı bir map veri yapısına dönüştürür.

Ham CSV verisini görmek için bu iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Map Operatörünü Ekleme

Şimdi veriyi dönüştürmek için scripting ekleyeceğiz; muhtemelen zaten aşina olduğunuz `.map()` operatörünü kullanarak. Bu operatör, her öğeyi dönüştürmek için kod yazabileceğiniz bir 'closure' alır.

!!! note

    **Closure**, etrafta taşınabilen ve daha sonra çalıştırılabilen bir kod bloğudur. Satır içi tanımladığınız bir fonksiyon olarak düşünün. Closure'lar süslü parantezler `{ }` ile yazılır ve parametre alabilir. Nextflow operatörlerinin çalışma şeklinin temelini oluştururlar ve bir süredir Nextflow yazıyorsanız, farkında olmadan zaten kullanıyor olabilirsiniz!

Map işlemi şöyle görünür:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Bu ilk **closure**'ımız - argüman olarak geçirebileceğiniz anonim bir fonksiyon (Python'daki lambda'lara veya JavaScript'teki arrow fonksiyonlarına benzer). Closure'lar, Nextflow operatörleriyle çalışmak için gereklidir.

`{ row -> return row }` closure'u `row` adında bir parametre alır (herhangi bir isim olabilir: `item`, `sample`, vb.).

`.map()` operatörü her channel öğesini işlediğinde, o öğeyi closure'unuza geçirir. Burada `row`, bir seferde bir CSV satırını tutar.

Bu değişikliği uygulayın ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

Girdiyi değiştirmeden döndürdüğümüz için öncekiyle aynı çıktıyı göreceksiniz. Bu, map operatörünün doğru çalıştığını onaylar. Şimdi veriyi dönüştürmeye başlayalım.

#### 1.1.3. Map Veri Yapısı Oluşturma

Şimdi her veri satırını dönüştürmek için closure'umuzun içine **scripting** mantığı yazacağız. Bu, veri akışını düzenlemek yerine bireysel veri öğelerini işlediğimiz yerdir.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Veri dönüşümü için scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

`sample_meta` map'i, ilgili bilgileri depolayan bir anahtar-değer veri yapısıdır (Python'daki sözlükler, JavaScript'teki nesneler veya Ruby'deki hash'ler gibi): örnek ID'si, organizma, doku tipi, dizileme derinliği ve kalite puanı.

Verimizi temizlemek için `.toLowerCase()` ve `.replaceAll()` gibi string manipülasyon metotlarını, CSV'den gelen string verileri uygun sayısal tiplere dönüştürmek için ise `.toInteger()` ve `.toDouble()` gibi tip dönüşüm metotlarını kullanırız.

Bu değişikliği uygulayın ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Koşullu Mantık Ekleme

Şimdi daha fazla scripting ekleyelim - bu sefer veri değerlerine göre kararlar vermek için ternary operatörü kullanarak.

Aşağıdaki değişikliği yapın:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

Ternary operatörü, `koşul ? doğruysa_değer : yanlışsa_değer` kalıbını izleyen bir if/else ifadesinin kısayoludur. Bu satırın anlamı: "Kalite 40'tan büyükse 'high' kullan, değilse 'normal' kullan". Kuzeni olan **Elvis operatörü** (`?:`), bir şey null veya boş olduğunda varsayılan değerler sağlar - bu kalıbı eğitimin ilerleyen bölümlerinde inceleyeceğiz.

Map toplama operatörü `+`, mevcut olanı değiştirmek yerine **yeni bir map** oluşturur. Bu satır, `sample_meta`'daki tüm anahtar-değer çiftlerini ve yeni `priority` anahtarını içeren yeni bir map oluşturur.

!!! Note

    Closure'lara geçirilen map'leri asla değiştirmeyin - her zaman `+` kullanarak (örneğin) yenilerini oluşturun. Nextflow'da aynı veri genellikle aynı anda birden fazla işlemden geçer. Bir map'i yerinde değiştirmek, aynı nesneye başvuran diğer işlemlerde öngörülemeyen yan etkilere neden olabilir. Yeni map'ler oluşturmak, her işlemin kendi temiz kopyasına sahip olmasını sağlar.

Değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Kalite puanlarına dayalı bir öncelik seviyesiyle metadata'mızı zenginleştirmek için koşullu mantığı başarıyla ekledik.

#### 1.1.5. `.subMap()` ile Map'leri Alt Kümeleme

`+` operatörü bir map'e anahtar eklerken, bazen tam tersini yapmanız gerekir - yalnızca belirli anahtarları çıkarmak. `.subMap()` metodu bunun için mükemmeldir.

Metadata'mızın yalnızca tanımlama alanlarını içeren basitleştirilmiş bir versiyonunu oluşturmak için bir satır ekleyelim:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Bu, `view()` işlemi tarafından görüntülenen tam metadata'yı ve `println` ile yazdırdığımız çıkarılmış alt kümeyi gösterir.

`.subMap()` metodu bir anahtar listesi alır ve yalnızca bu anahtarları içeren yeni bir map döndürür. Orijinal map'te bir anahtar yoksa, sonuca dahil edilmez.

Bu, farklı process'ler için farklı metadata versiyonları oluşturmanız gerektiğinde özellikle yararlıdır - bazıları tam metadata'ya ihtiyaç duyabilirken, diğerleri yalnızca minimum tanımlama alanlarına ihtiyaç duyabilir.

Artık bunlara ihtiyacımız olmadığı için, iş akışınızı önceki durumuna geri döndürmek için println ifadelerini kaldırın.

!!! tip "Map İşlemleri Özeti"

    - **Anahtar ekleme**: `map1 + [new_key: value]` - Ek anahtarlarla yeni map oluşturur
    - **Anahtar çıkarma**: `map1.subMap(['key1', 'key2'])` - Yalnızca belirtilen anahtarlarla yeni map oluşturur
    - **Her iki işlem de yeni map'ler oluşturur** - Orijinal map'ler değişmez

#### 1.1.6. Map'leri Birleştirme ve Sonuçları Döndürme

Şimdiye kadar yalnızca Nextflow topluluğunun 'meta map' dediği şeyi döndürüyorduk ve bu metadata'nın ilişkili olduğu dosyaları görmezden geliyorduk. Ancak Nextflow iş akışları yazıyorsanız, muhtemelen bu dosyalarla bir şeyler yapmak istiyorsunuzdur.

2 elemanlı bir tuple'dan oluşan bir channel yapısı çıktılayalım: zenginleştirilmiş metadata map'i ve karşılık gelen dosya yolu. Bu, Nextflow'da process'lere veri geçirmek için yaygın bir kalıptır.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Bu değişikliği uygulayın ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Bu `[meta, file]` tuple yapısı, Nextflow'da process'lere hem metadata hem de ilişkili dosyaları geçirmek için yaygın bir kalıptır.

!!! note

    **Map'ler ve Metadata**: Map'ler, Nextflow'da metadata ile çalışmanın temelidir. Metadata map'leriyle çalışma hakkında daha ayrıntılı bir açıklama için [Metadata ile Çalışma](./metadata.md) yan görevine bakın.

İş akışımız temel kalıbı göstermektedir: **dataflow işlemleri** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) verinin pipeline boyunca nasıl hareket ettiğini düzenlerken, `.map()` closure'u içindeki **scripting** (map'ler `[key: value]`, string metotları, tip dönüşümleri, ternary operatörleri) bireysel veri öğelerinin dönüşümünü gerçekleştirir.

### 1.2. Farklı Türleri Anlamak: Channel ve List

Şimdiye kadar iyi gidiyoruz; dataflow işlemleri ile scripting arasında ayrım yapabiliyoruz. Peki aynı metot adı her iki bağlamda da mevcut olduğunda ne olur?

Mükemmel bir örnek, Nextflow standart kütüphanesinde hem channel türleri hem de List türleri için mevcut olan `collect` metodudur. Bir List üzerindeki `collect()` metodu her elemanı dönüştürürken, bir channel üzerindeki `collect()` operatörü tüm channel emisyonlarını tek öğeli bir channel'a toplar.

Bunu bazı örnek verilerle gösterelim; önce channel `collect()` operatörünün ne yaptığını hatırlayarak başlayalım. `collect.nf` dosyasına bakın:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - birden fazla channel emisyonunu tek bir grupta toplar
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Adımlar:

- Örnek ID'lerinin bir List'ini tanımlayın
- Her örnek ID'sini ayrı ayrı yayan `fromList()` ile bir channel oluşturun
- Her öğeyi `view()` ile akış sırasında yazdırın
- Channel'ın `collect()` operatörü ile tüm öğeleri tek bir listede toplayın
- Toplanan sonucu (tüm örnek ID'lerini içeren tek öğe) ikinci bir `view()` ile yazdırın

Channel'ın yapısını değiştirdik, ancak verinin kendisini değiştirmedik.

Bunu doğrulamak için iş akışını çalıştırın:

```bash
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` her channel emisyonu için bir çıktı döndürür, bu nedenle bu tek çıktının orijinal 3 öğeyi tek bir listede grupladığını biliyoruz.

Şimdi bir List üzerinde `collect` metodunun çalışmasını görelim. Orijinal örnek ID'leri listesine List'in `collect` metodunu uygulamak için `collect.nf` dosyasını değiştirin:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla channel emisyonunu tek bir grupta toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - her elemanı dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla channel emisyonunu tek bir grupta toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

Bu yeni kod parçasında:

- Orijinal listedeki her örnek ID'sini dönüştürmek için List'in `collect` metodunu kullanan yeni bir `formatted_ids` değişkeni tanımlıyoruz
- Sonucu `println` ile yazdırıyoruz

Değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Bu sefer verinin yapısını DEĞİŞTİRMEDİK; listede hala 3 öğemiz var, ancak her öğeyi List'in `collect` metodu kullanarak dönüştürerek değiştirilmiş değerlerle yeni bir liste ürettik. Bu, bir channel üzerinde `map` operatörü kullanmaya benzer, ancak channel yerine bir List veri yapısı üzerinde çalışır.

`collect`, bir noktayı vurgulamak için kullandığımız uç bir örnektir. Temel ders şudur: iş akışları yazarken, her zaman **veri yapıları** (List, Map, vb.) ile **channel'lar** (dataflow yapıları) arasında ayrım yapın. İşlemler aynı adı paylaşabilir ancak çağrıldıkları türe bağlı olarak tamamen farklı davranabilir.

### 1.3. Spread Operatörü (`*.`) - Özellik Çıkarma Kısayolu

List'in `collect` metodu ile ilişkili olan spread operatörü (`*.`), koleksiyonlardan özellikleri çıkarmak için kısa bir yol sağlar. Temelde yaygın bir `collect` kalıbı için sözdizimsel bir kolaylıktır.

`collect.nf` dosyamıza bir gösteri ekleyelim:

=== "After"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla channel emisyonunu tek bir grupta toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - her elemanı dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Spread operatörü - kısa özellik erişimi
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Before"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla channel emisyonunu tek bir grupta toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - her elemanı dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Güncellenmiş iş akışını çalıştırın:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Spread operatörü `*.`, yaygın bir collect kalıbı için kısayoldur:

```groovy
// Bunlar eşdeğerdir:
def ids = samples*.id
def ids = samples.collect { it.id }

// Metot çağrılarıyla da çalışır:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Spread operatörü, bir nesne listesinden tek bir özellik çıkarmanız gerektiğinde özellikle yararlıdır - tam `collect` closure'unu yazmaktan daha okunabilirdir.

!!! tip "Spread ve Collect Ne Zaman Kullanılmalı"

    - **Spread (`*.`)** basit özellik erişimi için: `samples*.id`, `files*.name`
    - **Collect** dönüşümler veya karmaşık mantık için: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Çıkarımlar

Bu bölümde öğrendikleriniz:

- **Dataflow ve scripting**: Channel operatörleri verinin pipeline'ınız boyunca nasıl aktığını düzenlerken, scripting bireysel veri öğelerini dönüştürür
- **Türleri anlamak**: Aynı metot adı (örneğin `collect`), çağrıldığı türe bağlı olarak farklı davranabilir (Channel ve List)
- **Bağlam önemlidir**: Channel'larla (dataflow) mı yoksa veri yapılarıyla (scripting) mı çalıştığınızın her zaman farkında olun

Bu sınırları anlamak, hata ayıklama, dokümantasyon ve bakımı kolay iş akışları yazmak için gereklidir.

Sıradaki bölümde, gerçek dünya verilerini işlemek için gerekli olan string işleme yeteneklerini daha derinlemesine inceleyeceğiz.

---

## 2. String İşleme ve Dinamik Script Üretimi

String işlemede ustalaşmak, kırılgan iş akışlarını sağlam pipeline'lardan ayırır. Bu bölüm karmaşık dosya adlarını ayrıştırma, dinamik script üretimi ve değişken enterpolasyonunu kapsar.

### 2.1. Desen Eşleştirme ve Düzenli İfadeler

Biyoinformatik dosyaları genellikle metadata kodlayan karmaşık adlandırma kurallarına sahiptir. Bunu düzenli ifadelerle desen eşleştirme kullanarak otomatik olarak çıkaralım.

`main.nf` iş akışımıza geri dönecek ve dosya adlarından ek örnek bilgisi çıkarmak için desen eşleştirme mantığı ekleyeceğiz. Veri setimizdeki FASTQ dosyaları, `SAMPLE_001_S1_L001_R1_001.fastq.gz` gibi adlarla Illumina tarzı adlandırma kurallarını takip eder. Bunlar kriptik görünebilir, ancak aslında örnek ID'si, şerit numarası ve okuma yönü gibi yararlı metadata'yı kodlar. Bu adları ayrıştırmak için regex yeteneklerini kullanacağız.

Mevcut `main.nf` iş akışınızda aşağıdaki değişikliği yapın:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Veri dönüşümü için scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Veri dönüşümü için scripting
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Bu, temel **string işleme kavramlarını** gösterir:

1. `~/pattern/` sözdizimi ile **düzenli ifade sabitleri** - bu, ters eğik çizgileri escape etmeye gerek kalmadan bir regex deseni oluşturur
2. `=~` operatörü ile **desen eşleştirme** - bu, bir string'i bir regex desenine karşı eşleştirmeye çalışır
3. `[0][1]`, `[0][2]`, vb. ile grupları yakalayan **Matcher nesneleri** - `[0]` tüm eşleşmeyi, `[1]`, `[2]`, vb. parantez içindeki yakalanan grupları ifade eder

`^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` regex desenini inceleyelim:

| Desen               | Eşleşir                                  | Yakalar                          |
| ------------------- | ---------------------------------------- | -------------------------------- |
| `^(.+)`             | Baştan örnek adı                         | Grup 1: örnek adı                |
| `_S(\d+)`           | Örnek numarası `_S1`, `_S2`, vb.         | Grup 2: örnek numarası           |
| `_L(\d{3})`         | Şerit numarası `_L001`                   | Grup 3: şerit (3 rakam)          |
| `_(R[12])`          | Okuma yönü `_R1` veya `_R2`              | Grup 4: okuma yönü               |
| `_(\d{3})`          | Parça numarası `_001`                    | Grup 5: parça (3 rakam)          |
| `\.fastq(?:\.gz)?$` | Dosya uzantısı `.fastq` veya `.fastq.gz` | Yakalanmaz (?: yakalamayan grup) |

Bu, metadata'yı otomatik olarak çıkarmak için Illumina tarzı adlandırma kurallarını ayrıştırır.

Değiştirilmiş iş akışını çalıştırın:

```bash title="Test pattern matching"
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Bu, dosya adlarından zenginleştirilmiş metadata'yı gösterir.

### 2.2. Process'lerde Dinamik Script Üretimi

Process script blokları temelde shell'e geçirilen çok satırlı string'lerdir. Girdi özelliklerine göre farklı script string'leri dinamik olarak üretmek için **koşullu mantık** (if/else, ternary operatörleri) kullanabilirsiniz. Bu, process tanımlarını çoğaltmadan farklı girdi türlerini (tek uçlu ve çift uçlu dizileme okumaları gibi) ele almak için gereklidir.

Bu kalıbı gösteren bir process ekleyelim. `modules/fastp.nf` dosyasını açın ve inceleyin:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

Process, FASTQ dosyalarını girdi olarak alır ve adaptörleri kesmek ve düşük kaliteli okumaları filtrelemek için `fastp` aracını çalıştırır. Ne yazık ki, bu process'i yazan kişi örnek veri setimizdeki tek uçlu okumaları hesaba katmamış. İş akışımıza ekleyip ne olacağını görelim:

Önce, `main.nf` iş akışınızın en başına modülü dahil edin:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Ardından `ch_samples` channel'ını `FASTP` process'ine bağlamak için `workflow` bloğunu değiştirin:

=== "After"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Bu değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Process'in ikinci girdi dosyası için `null` değeriyle `fastp` çalıştırmaya çalıştığını görebilirsiniz, bu da başarısız olmasına neden oluyor. Bunun nedeni, veri setimizin tek uçlu okumalar içermesi, ancak process'in çift uçlu okumaları (aynı anda iki girdi dosyası) bekleyecek şekilde sabit kodlanmış olmasıdır.

`FASTP` process'inin `script:` bloğuna koşullu mantık ekleyerek bunu düzeltin. Bir if/else ifadesi okuma dosya sayısını kontrol eder ve komutu buna göre ayarlar.

=== "After"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Basit tek uçlu ve çift uçlu tespit
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Artık iş akışı hem tek uçlu hem de çift uçlu okumaları zarif bir şekilde ele alabilir. Koşullu mantık girdi dosyalarının sayısını kontrol eder ve `fastp` için uygun komutu oluşturur. Çalışıp çalışmadığını görelim:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

İyi görünüyor! Çalıştırılan gerçek komutları kontrol edersek (görev hash'inize göre özelleştirin):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Nextflow'un tek uçlu okumalar için doğru komutu seçtiğini görebiliriz:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Dinamik script mantığının bir diğer yaygın kullanımı [Nextflow for Science Genomics modülünde](../../nf4science/genomics/02_joint_calling) görülebilir. O modülde, çağrılan GATK process'i birden fazla girdi dosyası alabilir, ancak her birinin doğru bir komut satırı oluşturmak için `-V` ile ön ekli olması gerekir. Process, girdi dosyaları koleksiyonunu (`all_gvcfs`) doğru komut argümanlarına dönüştürmek için scripting kullanır:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Process script bloklarında scripting kullanmanın bu kalıpları son derece güçlüdür ve birçok senaryoda uygulanabilir - değişken girdi türlerini ele almaktan dosya koleksiyonlarından karmaşık komut satırı argümanları oluşturmaya kadar, process'lerinizi gerçek dünya verilerinin çeşitli gereksinimlerine gerçekten uyarlanabilir hale getirir.

### 2.3. Değişken Enterpolasyonu: Nextflow ve Shell Değişkenleri

Process script'leri Nextflow değişkenlerini, shell değişkenlerini ve komut yerine koymaları karıştırır; her birinin farklı enterpolasyon sözdizimi vardır. Yanlış sözdizimi kullanmak hatalara yol açar. Bunları, işleme raporu oluşturan bir process ile inceleyelim.

`modules/generate_report.nf` modül dosyasına bakın:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Bu process, örnek ID'si ve dosya adıyla basit bir rapor yazar. Şimdi farklı türde değişkenleri karıştırmamız gerektiğinde ne olduğunu görmek için çalıştıralım.

Process'i `main.nf` dosyanıza dahil edin ve iş akışına ekleyin:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

Şimdi iş akışını çalıştırın ve `results/reports/` dizinindeki oluşturulan raporları kontrol edin. Her örnek hakkında temel bilgiler içermelidirler.

<!-- TODO: add the run command -->

??? success "Komut çıktısı"

    ```console
    <!-- TODO: output -->
    ```

Peki işlemenin ne zaman ve nerede yapıldığı hakkında bilgi eklemek istersek? Process'i, rapora geçerli kullanıcı, hostname ve tarihi dahil etmek için **shell** değişkenleri ve biraz komut yerine koyma kullanacak şekilde değiştirelim:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Bunu çalıştırırsanız bir hata fark edeceksiniz - Nextflow `${USER}` ifadesini mevcut olmayan bir Nextflow değişkeni olarak yorumlamaya çalışır.

??? failure "Komut çıktısı"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Bash'in bunları ele alabilmesi için escape etmemiz gerekir.

Shell değişkenlerini ve komut yerine koymaları ters eğik çizgi (`\`) ile escape ederek düzeltin:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

Artık çalışıyor! Ters eğik çizgi (`\`), Nextflow'a "bunu yorumlama, Bash'e geçir" der.

### Çıkarımlar

Bu bölümde **string işleme** tekniklerini öğrendiniz:

- **Dosya ayrıştırma için düzenli ifadeler**: Karmaşık dosya adlandırma kurallarından metadata çıkarmak için `=~` operatörü ve regex desenleri (`~/pattern/`) kullanma
- **Dinamik script üretimi**: Girdi özelliklerine göre farklı script string'leri üretmek için koşullu mantık (if/else, ternary operatörleri) kullanma
- **Değişken enterpolasyonu**: Nextflow'un string'leri ne zaman yorumladığını ve shell'in ne zaman yorumladığını anlama
  - `${var}` - Nextflow değişkenleri (iş akışı derleme zamanında Nextflow tarafından yorumlanır)
  - `\${var}` - Shell ortam değişkenleri (escape edilmiş, çalışma zamanında bash'e geçirilir)
  - `\$(cmd)` - Shell komut yerine koyma (escape edilmiş, çalışma zamanında bash tarafından yürütülür)

Bu string işleme ve üretim kalıpları, gerçek dünya biyoinformatik iş akışlarında karşılaşacağınız çeşitli dosya formatlarını ve adlandırma kurallarını ele almak için gereklidir.

---

## 3. Yeniden Kullanılabilir Fonksiyonlar Oluşturma

Channel operatörlerinde veya process tanımlarında satır içi karmaşık iş akışı mantığı, okunabilirliği ve bakımı azaltır. **Fonksiyonlar**, bu mantığı adlandırılmış, yeniden kullanılabilir bileşenlere çıkarmanızı sağlar.

Map işlemimiz uzun ve karmaşık hale geldi. `def` anahtar sözcüğünü kullanarak yeniden kullanılabilir bir fonksiyona çıkaralım.

Mevcut iş akışımızla bunun nasıl göründüğünü göstermek için, `separateMetadata` adlı yeniden kullanılabilir bir fonksiyon tanımlamak üzere `def` kullanarak aşağıdaki değişikliği yapın:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

Bu mantığı bir fonksiyona çıkararak, gerçek iş akışı mantığını çok daha temiz bir hale getirdik:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Bu, iş akışı mantığını bir bakışta okumayı ve anlamayı çok daha kolay hale getirir. `separateMetadata` fonksiyonu, metadata'yı ayrıştırma ve zenginleştirme için tüm karmaşık mantığı kapsüller, bu da onu yeniden kullanılabilir ve test edilebilir kılar.

Hala çalıştığından emin olmak için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

Çıktı, her iki process'in de başarıyla tamamlandığını göstermelidir. İş akışı artık çok daha temiz ve bakımı daha kolaydır; tüm karmaşık metadata işleme mantığı `separateMetadata` fonksiyonunda kapsüllenmiştir.

### Çıkarımlar

Bu bölümde **fonksiyon oluşturmayı** öğrendiniz:

- **`def` ile fonksiyon tanımlama**: Adlandırılmış fonksiyonlar oluşturmak için anahtar sözcük (Python'daki `def` veya JavaScript'teki `function` gibi)
- **Fonksiyon kapsamı**: Script düzeyinde tanımlanan fonksiyonlara Nextflow iş akışınız boyunca erişilebilir
- **Dönüş değerleri**: Fonksiyonlar otomatik olarak son ifadeyi döndürür veya açık `return` kullanır
- **Daha temiz kod**: Karmaşık mantığı fonksiyonlara çıkarmak, herhangi bir dilde temel bir yazılım mühendisliği uygulamasıdır

Sıradaki bölümde, dinamik kaynak tahsisi için process directive'lerinde closure'ları nasıl kullanacağımızı inceleyeceğiz.

---

## 4. Closure'lar ile Dinamik Kaynak Directive'leri

Şimdiye kadar scripting'i process'lerin `script` bloğunda kullandık. Ancak **closure'lar** (Bölüm 1.1'de tanıtılan), özellikle dinamik kaynak tahsisi için process directive'lerinde de son derece yararlıdır. FASTP process'imize örnek özelliklerine göre uyum sağlayan kaynak directive'leri ekleyelim.

### 4.1. Örneğe özgü kaynak tahsisi

Şu anda FASTP process'imiz varsayılan kaynakları kullanıyor. Yüksek derinlikli örnekler için daha fazla CPU tahsis ederek daha akıllı hale getirelim. Dinamik bir `cpus` directive'i ve statik bir `memory` directive'i eklemek için `modules/fastp.nf` dosyasını düzenleyin:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

`{ meta.depth > 40000000 ? 2 : 1 }` closure'u **ternary operatörünü** (Bölüm 1.1'de ele alınan) kullanır ve her görev için değerlendirilir, böylece örnek başına kaynak tahsisi yapılabilir. Yüksek derinlikli örnekler (>40M okuma) 2 CPU alırken, diğerleri 1 CPU alır.

!!! note "Directive'lerde Girdi Değişkenlerine Erişim"

    Closure, herhangi bir girdi değişkenine (burada `meta` gibi) erişebilir çünkü Nextflow bu closure'ları her görev yürütme bağlamında değerlendirir.

İş akışını tekrar `-ansi-log false` seçeneğiyle çalıştırarak görev hash'lerini daha kolay görebilirsiniz.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Herhangi bir görev için CPU tahsisini görmek üzere çalıştırılan `docker` komutunu kontrol edebilirsiniz:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Şöyle bir şey görmelisiniz:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Bu örnekte 2 CPU isteyen bir örnek seçtik (`--cpu-shares 2048`), çünkü yüksek derinlikli bir örnekti, ancak örnek derinliğine bağlı olarak farklı CPU tahsisleri görmelisiniz. Diğer görevler için de bunu deneyin.

### 4.2. Yeniden deneme stratejileri

Bir diğer güçlü kalıp, yeniden deneme stratejileri için `task.attempt` kullanmaktır. Bunun neden yararlı olduğunu göstermek için, FASTP'ye ayrılan bellek tahsisini ihtiyaç duyduğundan daha azına düşüreceğiz. `modules/fastp.nf` dosyasındaki `memory` directive'ini `1.GB` olarak değiştirin:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... ve iş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Bu, process'in bellek sınırlarını aştığı için sonlandırıldığını gösterir.

Bu, gerçek dünya iş akışlarında çok yaygın bir senaryodur - bazen bir görevin ne kadar belleğe ihtiyaç duyacağını çalıştırana kadar bilemezsiniz.

İş akışımızı daha dayanıklı hale getirmek için, yine bir Groovy closure'u kullanarak her denemede bellek tahsisini artıran bir yeniden deneme stratejisi uygulayabiliriz. `memory` directive'ini temel belleği `task.attempt` ile çarpacak şekilde değiştirin ve `errorStrategy 'retry'` ile `maxRetries 2` directive'lerini ekleyin:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Artık process yetersiz bellek nedeniyle başarısız olursa, Nextflow daha fazla bellekle yeniden deneyecektir:

- İlk deneme: 1 GB (task.attempt = 1)
- İkinci deneme: 2.GB (task.attempt = 2)

... ve `maxRetries` sınırına kadar devam eder.

### Çıkarımlar

Closure'larla dinamik directive'ler şunları yapmanıza olanak tanır:

- Girdi özelliklerine göre kaynak tahsis etme
- Artan kaynaklarla otomatik yeniden deneme stratejileri uygulama
- Birden fazla faktörü birleştirme (metadata, deneme sayısı, öncelikler)
- Karmaşık kaynak hesaplamaları için koşullu mantık kullanma

Bu, iş akışlarınızı hem daha verimli (fazla tahsis yapmama) hem de daha dayanıklı (daha fazla kaynakla otomatik yeniden deneme) hale getirir.

---

## 5. Koşullu Mantık ve Process Kontrolü

Daha önce channel verilerini dönüştürmek için scripting ile `.map()` kullandık. Şimdi hangi process'lerin verilere göre çalıştığını kontrol etmek için koşullu mantık kullanacağız - farklı örnek türlerine uyum sağlayan esnek iş akışları için gerekli.

Nextflow'un [dataflow operatörleri](https://www.nextflow.io/docs/latest/reference/operator.html), çalışma zamanında değerlendirilen closure'lar alır ve koşullu mantığın channel içeriğine göre iş akışı kararlarını yönlendirmesini sağlar.

### 5.1. `.branch()` ile Yönlendirme

Örneğin, dizileme örneklerimizin yalnızca belirli bir kapsam eşiğinin üzerindeki insan örnekleri ise FASTP ile kırpılması gerektiğini varsayalım. Fare örnekleri veya düşük kapsamlı örnekler bunun yerine Trimgalore ile çalıştırılmalıdır (bu yapay bir örnek, ancak konuyu açıklıyor).

`modules/trimgalore.nf` dosyasında basit bir Trimgalore process'i sağladık; isterseniz bakabilirsiniz, ancak bu alıştırma için detaylar önemli değil. Asıl nokta, örnekleri metadata'larına göre yönlendirmek istememizdir.

`modules/trimgalore.nf` dosyasından yeni modülü dahil edin:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... ve ardından örnekleri metadata'larına göre dallandırmak ve uygun kırpma process'ine yönlendirmek için `main.nf` iş akışınızı şu şekilde değiştirin:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Bu değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Burada, örnekleri metadata'larına göre yönlendirmek için `.branch{}` operatörünün içinde küçük ama güçlü koşullu ifadeler kullandık. Yüksek kapsamlı insan örnekleri `FASTP`'den geçerken, diğer tüm örnekler `TRIMGALORE`'dan geçer.

### 5.2. `.filter()` ile Doğruluk Değeri Kullanımı

İş akışı yürütmesini kontrol etmek için bir diğer güçlü kalıp, pipeline'dan hangi öğelerin devam edeceğini belirlemek için bir closure kullanan `.filter()` operatörüdür. Filter closure'unun içinde, hangi öğelerin geçeceğine karar veren **boolean ifadeleri** yazarsınız.

Nextflow (birçok dinamik dil gibi), boolean bağlamlarında hangi değerlerin `true` veya `false` olarak değerlendirildiğini belirleyen bir **"doğruluk değeri" (truthiness)** kavramına sahiptir:

- **Doğru (Truthy)**: Null olmayan değerler, boş olmayan string'ler, sıfır olmayan sayılar, boş olmayan koleksiyonlar
- **Yanlış (Falsy)**: `null`, boş string'ler `""`, sıfır `0`, boş koleksiyonlar `[]` veya `[:]`, `false`

Bu, `meta.id`'nin tek başına (açık `!= null` olmadan) ID'nin var olup olmadığını ve boş olmadığını kontrol ettiği anlamına gelir. Kalite gereksinimlerimizi karşılamayan örnekleri filtrelemek için bunu kullanalım.

Branch işleminden önce aşağıdakileri ekleyin:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Geçersiz veya düşük kaliteli örnekleri filtrele
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

İş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Bazı örnekleri hariç tutan bir filtre seçtiğimiz için daha az görev yürütüldü.

`meta.id && meta.organism && meta.depth >= 25000000` filtre ifadesi, doğruluk değerini açık karşılaştırmalarla birleştirir:

- `meta.id && meta.organism` her iki alanın da var olduğunu ve boş olmadığını kontrol eder (doğruluk değeri kullanarak)
- `meta.depth >= 25000000` açık bir karşılaştırma ile yeterli dizileme derinliğini sağlar

!!! note "Pratikte Doğruluk Değeri"

    `meta.id && meta.organism` ifadesi şunu yazmaktan daha kısadır:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Bu, filtreleme mantığını çok daha temiz ve okunması kolay hale getirir.

### Çıkarımlar

Bu bölümde, `.branch{}` ve `.filter{}` gibi Nextflow operatörlerinin closure arayüzlerini kullanarak iş akışı yürütmesini kontrol etmek için koşullu mantık kullanmayı ve kısa koşullu ifadeler yazmak için doğruluk değerinden yararlanmayı öğrendiniz.

Pipeline'ımız artık örnekleri uygun process'lere akıllıca yönlendiriyor, ancak üretim iş akışlarının geçersiz verileri zarif bir şekilde ele alması gerekir. İş akışımızı eksik veya null değerlere karşı dayanıklı hale getirelim.

---

## 6. Güvenli Gezinme ve Elvis Operatörleri

`separateMetadata` fonksiyonumuz şu anda tüm CSV alanlarının mevcut ve geçerli olduğunu varsayıyor. Peki eksik verilerle ne olur? Hadi öğrenelim.

### 6.1. Problem: Var Olmayan Özelliklere Erişim

Diyelim ki isteğe bağlı dizileme çalıştırma bilgisi için destek eklemek istiyoruz. Bazı laboratuvarlarda, örneklerin dizileme çalıştırma ID'si veya parti numarası için ek bir alanı olabilir, ancak mevcut CSV'mizde bu sütun yok. Yine de erişmeye çalışalım.

`separateMetadata` fonksiyonunu bir run_id alanı içerecek şekilde değiştirin:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Şimdi iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Bu, bir NullPointerException ile çöker.

Sorun şu ki `row.run_id`, CSV'mizde `run_id` sütunu olmadığı için `null` döndürür. `null` üzerinde `.toUpperCase()` çağırmaya çalıştığımızda çöker. İşte güvenli gezinme operatörü burada devreye girer.

### 6.2. Güvenli Gezinme Operatörü (`?.`)

Güvenli gezinme operatörü (`?.`), bir `null` değer üzerinde çağrıldığında istisna fırlatmak yerine `null` döndürür. `?.`'dan önceki nesne `null` ise, tüm ifade metodu çalıştırmadan `null` olarak değerlendirilir.

Fonksiyonu güvenli gezinme kullanacak şekilde güncelleyin:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    <!-- TODO: output -->
    ```

Çökme yok! İş akışı artık eksik alanı zarif bir şekilde ele alıyor. `row.run_id` `null` olduğunda, `?.` operatörü `.toUpperCase()` çağrısını engeller ve `run_id`, istisna fırlatmak yerine `null` olur.

### 6.3. Varsayılanlar için Elvis Operatörü (`?:`)

Elvis operatörü (`?:`), sol taraf "yanlış" (daha önce açıklandığı gibi) olduğunda varsayılan değerler sağlar. Adını Elvis Presley'den alır çünkü `?:` yandan bakıldığında onun ünlü saçına ve gözlerine benzer!

Güvenli gezinme kullandığımıza göre, `run_id` bu alana sahip olmayan örnekler için `null` olacaktır. Bir varsayılan değer sağlamak ve `sample_meta` map'imize eklemek için Elvis operatörünü kullanalım:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Sonuçları görmek için iş akışına bir `view()` operatörü de ekleyin:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Mükemmel! Artık tüm örneklerin gerçek çalıştırma ID'leri (büyük harfle) veya varsayılan 'UNSPECIFIED' değeriyle bir `run` alanı var. `?.` ve `?:` kombinasyonu hem güvenlik (çökme yok) hem de mantıklı varsayılanlar sağlar.

Çalıştığını onayladığımıza göre `.view()` operatörünü kaldırın.

!!! tip "Güvenli Gezinme ve Elvis'i Birleştirme"

    `value?.method() ?: 'default'` kalıbı, üretim iş akışlarında yaygındır:

    - `value?.method()` - Metodu güvenle çağırır, `value` `null` ise `null` döndürür
    - `?: 'default'` - Sonuç `null` ise yedek değer sağlar

    Bu kalıp, eksik/tamamlanmamış verileri zarif bir şekilde ele alır.

Bu operatörleri fonksiyonlarda, operatör closure'larında (`.map{}`, `.filter{}`), process script'lerinde ve yapılandırma dosyalarında tutarlı bir şekilde kullanın. Gerçek dünya verileriyle çalışırken çökmeleri önlerler.

### Çıkarımlar

- **Güvenli gezinme (`?.`)**: Null değerlerde çökmeleri önler - istisna fırlatmak yerine null döndürür
- **Elvis operatörü (`?:`)**: Varsayılanlar sağlar - `value ?: 'default'`
- **Birleştirme**: `value?.method() ?: 'default'` yaygın kalıptır

Bu operatörler, iş akışlarını eksik verilere karşı dayanıklı hale getirir - gerçek dünya çalışmaları için gereklidir.

---

## 7. `error()` ve `log.warn` ile Doğrulama

Bazen girdi parametreleri geçersiz olduğunda iş akışını hemen durdurmanız gerekir. Nextflow'da doğrulama mantığı uygulamak için `error()` ve `log.warn` gibi yerleşik fonksiyonların yanı sıra `if` ifadeleri ve boolean mantığı gibi standart programlama yapılarını kullanabilirsiniz. İş akışımıza doğrulama ekleyelim.

Workflow bloğunuzdan önce bir doğrulama fonksiyonu oluşturun, bunu workflow'dan çağırın ve channel oluşturmayı CSV dosya yolu için bir parametre kullanacak şekilde değiştirin. Parametre eksikse veya dosya mevcut değilse, net bir mesajla yürütmeyi durdurmak için `error()` çağırın.

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Girdi parametresinin sağlandığını kontrol et
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // CSV dosyasının var olduğunu kontrol et
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Şimdi CSV dosyası olmadan çalıştırmayı deneyin:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

İş akışı, daha sonra gizemli bir şekilde başarısız olmak yerine hemen net bir hata mesajıyla durur.

Şimdi var olmayan bir dosyayla çalıştırın:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Son olarak, doğru dosyayla çalıştırın:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Komut çıktısı"

    ```console
    <!-- TODO: output -->
    ```

Bu sefer başarıyla çalışır.

`separateMetadata` fonksiyonunun içine de doğrulama ekleyebilirsiniz. Düşük dizileme derinliğine sahip örnekler için uyarı vermek, ancak yine de iş akışının devam etmesine izin vermek için ölümcül olmayan `log.warn` kullanalım:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Verinin anlamlı olduğunu doğrula
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

İş akışını orijinal CSV ile tekrar çalıştırın:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Örneklerden biri için düşük dizileme derinliği hakkında bir uyarı görüyoruz.

### Çıkarımlar

- **`error()`**: İş akışını net bir mesajla hemen durdurur
- **`log.warn`**: İş akışını durdurmadan uyarı verir
- **Erken doğrulama**: Yardımcı hatalarla hızlı başarısızlık için işlemeden önce girdileri kontrol edin
- **Doğrulama fonksiyonları**: İş akışı başlangıcında çağrılabilecek yeniden kullanılabilir doğrulama mantığı oluşturun

Doğru doğrulama, sorunları erken yakalayarak ve net hata mesajları vererek iş akışlarını daha dayanıklı ve kullanıcı dostu hale getirir.

---

## 8. Workflow Olay İşleyicileri

Şimdiye kadar iş akışı script'lerimizde ve process tanımlarımızda kod yazıyorduk. Ancak bilmeniz gereken önemli bir özellik daha var: workflow olay işleyicileri.

Olay işleyicileri, iş akışınızın yaşam döngüsünde belirli noktalarda çalışan closure'lardır. Günlükleme, bildirimler veya temizlik işlemleri eklemek için mükemmeldirler. Bu işleyiciler, workflow tanımınızın yanında iş akışı script'inizde tanımlanmalıdır.

### 8.1. `onComplete` İşleyicisi

En yaygın kullanılan olay işleyicisi, iş akışınız tamamlandığında (başarılı veya başarısız olsun) çalışan `onComplete`'dir. Pipeline sonuçlarımızı özetlemek için bir tane ekleyelim.

Olay işleyicisini `main.nf` dosyanıza, workflow tanımınızın içine ekleyin:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Bu closure, iş akışı tamamlandığında çalışır. İçinde, yürütme hakkında yararlı özellikler sağlayan `workflow` nesnesine erişebilirsiniz.

İş akışınızı çalıştırın, sonunda bu özeti göreceksiniz!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

Koşullu mantık ekleyerek daha kullanışlı hale getirelim:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Artık başarı/başarısızlık mesajı ve belirtildiyse çıktı dizini dahil daha bilgilendirici bir özet alıyoruz:

<!-- TODO: add run command -->

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

Dosya işlemleri kullanarak özeti bir dosyaya da yazabilirsiniz:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // Bir log dosyasına yaz
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. `onError` İşleyicisi

`onComplete`'in yanı sıra, kullanabileceğiniz bir olay işleyicisi daha var: yalnızca iş akışı başarısız olursa çalışan `onError`:

```groovy title="main.nf - onError handler"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Detaylı hata günlüğü yaz
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

İş akışı script'inizde birden fazla işleyiciyi birlikte kullanabilirsiniz:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### Çıkarımlar

Bu bölümde öğrendikleriniz:

- **Olay işleyici closure'ları**: İş akışı script'inizde farklı yaşam döngüsü noktalarında çalışan closure'lar
- **`onComplete` işleyicisi**: Yürütme özetleri ve sonuç raporlama için
- **`onError` işleyicisi**: Hata işleme ve başarısızlık günlükleme için
- **Workflow nesne özellikleri**: `workflow.success`, `workflow.duration`, `workflow.errorMessage`, vb. erişim

Olay işleyicileri, iş akışı script'lerinize sofistike günlükleme ve bildirim yetenekleri eklemek için Nextflow dilinin tüm gücünü nasıl kullanabileceğinizi gösterir.

---

## Özet

Tebrikler, başardınız!

Bu yan görev boyunca, temel metadata işlemeden sofistike, üretime hazır bir iş akışına evrilen kapsamlı bir örnek işleme pipeline'ı oluşturdunuz.
Her bölüm bir öncekinin üzerine inşa edilerek, programlama yapılarının basit iş akışlarını güçlü veri işleme sistemlerine nasıl dönüştürdüğünü aşağıdaki faydalarla gösterdi:

- **Daha temiz kod**: Dataflow ve scripting arasındaki farkı anlamak, daha düzenli iş akışları yazmanıza yardımcı olur
- **Sağlam işleme**: Güvenli gezinme ve Elvis operatörleri, iş akışlarını eksik verilere karşı dayanıklı hale getirir
- **Esnek işleme**: Koşullu mantık, iş akışlarınızın farklı örnek türlerini uygun şekilde işlemesini sağlar
- **Uyarlanabilir kaynaklar**: Dinamik directive'ler, girdi özelliklerine göre kaynak kullanımını optimize eder

Bu ilerleme, biyoinformatik pipeline'ların gerçek dünya evrimini yansıtır: birkaç örneği işleyen araştırma prototiplerinden, laboratuvarlar ve kurumlar genelinde binlerce örneği işleyen üretim sistemlerine.
Çözdüğünüz her zorluk ve öğrendiğiniz her kalıp, geliştiricilerin Nextflow iş akışlarını ölçeklendirirken karşılaştıkları gerçek sorunları yansıtır.

Bu kalıpları kendi çalışmanızda uygulamak, sağlam, üretime hazır iş akışları oluşturmanıza olanak tanıyacaktır.

### Temel kalıplar

1.  **Dataflow ve Scripting:** Dataflow işlemleri (channel düzenleme) ile scripting (veriyi manipüle eden kod) arasında ayrım yapmayı ve Channel ve List üzerindeki `collect` gibi farklı türlerdeki işlemler arasındaki kritik farkları öğrendiniz.

    - Dataflow: channel düzenleme

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: koleksiyonlar üzerinde veri işleme

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **İleri Düzey String İşleme**: Dosya adlarını ayrıştırmak için düzenli ifadelerde, process'lerde dinamik script üretiminde ve değişken enterpolasyonunda (Nextflow, Bash ve Shell) ustalaştınız.

    - Desen eşleştirme

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Koşullu dönüşlü fonksiyon

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Dosya koleksiyonundan komut argümanlarına (process script bloğunda)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Yeniden Kullanılabilir Fonksiyonlar Oluşturma**: Karmaşık mantığı channel operatörlerinden çağrılabilecek adlandırılmış fonksiyonlara çıkarmayı öğrenerek iş akışlarını daha okunabilir ve bakımı kolay hale getirdiniz.

    - Adlandırılmış fonksiyon tanımlama

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* kısalık için kod gizlendi */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* kısalık için kod gizlendi */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Adlandırılmış fonksiyonu workflow'da çağırma

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closure'lar ile Dinamik Kaynak Directive'leri**: Girdi özelliklerine dayalı uyarlanabilir kaynak tahsisi için process directive'lerinde closure'ları kullanmayı keşfettiniz.

    - Adlandırılmış closure'lar ve kompozisyon

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Kapsam erişimli closure'lar

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Koşullu Mantık ve Process Kontrolü**: `.branch()` ve `.filter()` operatörlerini kullanarak akıllı yönlendirme eklediniz, kısa koşullu ifadeler için doğruluk değerinden yararlandınız.

    - Verileri farklı iş akışı dallarına yönlendirmek için `.branch()` kullanma

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Groovy Truth ile boolean değerlendirme

    ```groovy
    if (sample.files) println "Has files"
    ```

    - Doğruluk değeri ile verileri alt kümelemek için `filter()` kullanma

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Güvenli Gezinme ve Elvis Operatörleri**: Null-safe özellik erişimi için `?.` ve varsayılan değerler sağlamak için `?:` kullanarak pipeline'ı eksik verilere karşı dayanıklı hale getirdiniz.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **error() ve log.warn ile Doğrulama**: Girdileri erken doğrulamayı ve net hata mesajlarıyla hızlı başarısızlık sağlamayı öğrendiniz.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Yapılandırma Olay İşleyicileri**: Günlükleme, bildirimler ve yaşam döngüsü yönetimi için workflow olay işleyicilerini (`onComplete` ve `onError`) kullanmayı öğrendiniz.

    - Günlükleme ve bildirim için `onComplete` kullanma

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Özellikle başarısızlık durumunda işlem yapmak için `onError` kullanma

    ```groovy
    workflow.onError = {
        // Detaylı hata günlüğü yaz
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Ek kaynaklar

- [Nextflow Dil Referansı](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operatörleri](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Sözdizimi](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standart Kütüphane](https://nextflow.io/docs/latest/reference/stdlib.html)

Daha ileri düzey özellikleri keşfetmeniz gerektiğinde bu kaynakları incelediğinizden emin olun.

Becerilerinizi geliştirmek ve genişletmek için pratik yapmaktan fayda göreceksiniz:

- Dataflow ve scripting arasında uygun ayrımla daha temiz iş akışları yazma
- Nextflow, Bash ve shell değişkenleriyle yaygın tuzaklardan kaçınmak için değişken enterpolasyonunda ustalaşma
- Verimli, uyarlanabilir iş akışları için dinamik kaynak directive'leri kullanma
- Dosya koleksiyonlarını uygun şekilde biçimlendirilmiş komut satırı argümanlarına dönüştürme
- Regex ve string işleme kullanarak farklı dosya adlandırma kurallarını ve girdi formatlarını zarif bir şekilde ele alma
- İleri düzey closure kalıpları ve fonksiyonel programlama kullanarak yeniden kullanılabilir, bakımı kolay kod oluşturma
- Koleksiyon işlemlerini kullanarak karmaşık veri setlerini işleme ve düzenleme
- İş akışlarınızı üretime hazır hale getirmek için doğrulama, hata işleme ve günlükleme ekleme
- Olay işleyicileri ile iş akışı yaşam döngüsü yönetimi uygulama

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) geri dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ altındaki düğmeye tıklayın.
