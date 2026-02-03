# Temel Nextflow Kodlama Kalıpları

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, Java Sanal Makinesi üzerinde çalışan bir programlama dilidir. Nextflow, [Groovy](http://groovy-lang.org/) üzerine kurulu olsa ve sözdiziminin çoğunu onunla paylaşsa da, Nextflow "uzantılı Groovy"den fazlasıdır -- tam olarak belirtilmiş bir [sözdizimi](https://nextflow.io/docs/latest/reference/syntax.html) ve [standart kütüphane](https://nextflow.io/docs/latest/reference/stdlib.html) ile bağımsız bir dildir.

Değişkenler, map'ler ve listeler için temel sözdiziminin ötesine geçmeden çok sayıda Nextflow kodu yazabilirsiniz. Çoğu Nextflow eğitimi iş akışı düzenlemesine (kanallar, process'ler ve veri akışı) odaklanır ve sadece bununla bile şaşırtıcı derecede ileri gidebilirsiniz.

Ancak, veri manipüle etmeniz, karmaşık dosya adlarını ayrıştırmanız, koşullu mantık uygulamanız veya sağlam üretim iş akışları oluşturmanız gerektiğinde, kodunuzun iki ayrı yönünü düşünmek yardımcı olur: **dataflow** (kanallar, operatörler, process'ler ve workflow'lar) ve **scripting** (closure'lar, fonksiyonlar ve process script'leri içindeki kod). Bu ayrım bir nebze keyfi olsa da -hepsi Nextflow kodu- ne zaman pipeline'ınızı düzenleyip ne zaman veriyi manipüle ettiğinizi anlamak için yararlı bir zihinsel model sağlar. Her ikisinde de ustalaşmak, açık ve bakımı yapılabilir iş akışları yazma yeteneğinizi önemli ölçüde artırır.

### Öğrenme hedefleri

Bu yan görev, sizi temel kavramlardan üretime hazır kalıplara kadar uygulamalı bir yolculuğa çıkarır.
Basit bir CSV okuyan iş akışını gerçekçi zorluklarla adım adım geliştirerek karmaşık bir biyoinformatik pipeline'a dönüştüreceğiz:

- **Sınırları anlamak:** Dataflow işlemleri ile scripting arasında ayrım yapmak ve bunların birlikte nasıl çalıştığını anlamak
- **Veri manipülasyonu:** Güçlü operatörler kullanarak map'leri ve koleksiyonları çıkartmak, dönüştürmek ve alt kümelemek
- **String işleme:** Regex kalıpları ile karmaşık dosya adlandırma şemalarını ayrıştırmak ve değişken enterpolasyonunda ustalaşmak
- **Yeniden kullanılabilir fonksiyonlar:** Daha temiz, daha bakımı kolay iş akışları için karmaşık mantığı adlandırılmış fonksiyonlara çıkartmak
- **Dinamik mantık:** Farklı girdi türlerine uyum sağlayan process'ler oluşturmak ve dinamik kaynak tahsisi için closure'lar kullanmak
- **Koşullu yönlendirme:** Örnekleri meta veri özelliklerine göre farklı process'lerden akıllıca yönlendirmek
- **Güvenli işlemler:** Eksik verileri null-safe operatörlerle zarif bir şekilde işlemek ve girdileri açık hata mesajlarıyla doğrulamak
- **Yapılandırma tabanlı işleyiciler:** Günlükleme, bildirimler ve yaşam döngüsü yönetimi için workflow olay işleyicileri kullanmak

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (process'ler, kanallar, operatörler, dosyalarla çalışma, meta veri) kullanmakta rahat olmalısınız
- Yaygın programlama yapılarına (değişkenler, map'ler, listeler) temel aşinalığınız olmalı

Bu eğitim, karşılaştığımız programlama kavramlarını açıklayacaktır, bu nedenle kapsamlı programlama deneyimine ihtiyacınız yoktur.
Temel kavramlarla başlayıp gelişmiş kalıplara doğru ilerleyeceğiz.

---

## 0. Başlayın

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/essential_scripting_patterns
```

#### Materyalleri gözden geçirin

Ana iş akışı dosyasını ve örnek veri dosyalarını içeren bir `data` dizini bulacaksınız.

```console title="Dizin içeriği"
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

Örnek CSV dosyamız, özelliklerine göre farklı işleme ihtiyacı olan biyolojik örnekler hakkında bilgi içerir:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Gerçek biyoinformatik iş akışlarında karşılaşacağınız pratik programlama tekniklerini keşfetmek için bu gerçekçi veri setini kullanacağız.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım
<!-- - [ ] I understand the assignment -->

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Dataflow vs Scripting: Sınırları Anlamak

### 1.1. Neyin Ne Olduğunu Belirleme

Nextflow iş akışları yazarken, **dataflow** (verinin kanallar ve process'ler arasında nasıl hareket ettiği) ile **scripting** (veriyi manipüle eden ve kararlar veren kod) arasında ayrım yapmak önemlidir. Bunların birlikte nasıl çalıştığını gösteren bir iş akışı oluşturalım.

#### 1.1.1. Temel Nextflow İş Akışı

Sadece CSV dosyasını okuyan basit bir iş akışıyla başlayın (bunu `main.nf` dosyasında sizin için zaten yaptık):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` bloğu pipeline yapımızı tanımlarken, `channel.fromPath()` bir dosya yolundan bir kanal oluşturur. `.splitCsv()` operatörü CSV dosyasını işler ve her satırı bir map veri yapısına dönüştürür.

Ham CSV verilerini görmek için bu iş akışını çalıştırın:

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

Şimdi veriyi dönüştürmek için scripting ekleyeceğiz, muhtemelen zaten aşina olduğunuz `.map()` operatörünü kullanarak. Bu operatör, her öğeyi dönüştürmek için kod yazabileceğimiz bir 'closure' alır.

!!! note

    **Closure**, etrafta taşınabilen ve daha sonra çalıştırılabilen bir kod bloğudur. Bunu satır içinde tanımladığınız bir fonksiyon gibi düşünün. Closure'lar süslü parantezlerle `{ }` yazılır ve parametre alabilir. Nextflow operatörlerinin nasıl çalıştığının temelindedirler ve bir süredir Nextflow yazıyorsanız, farkında bile olmadan onları kullanıyor olabilirsiniz!

İşte o map işlemi nasıl görünüyor:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Bu bizim ilk **closure**'ımız - argüman olarak geçebileceğiniz anonim bir fonksiyon (Python'daki lambda'lara veya JavaScript'teki ok fonksiyonlarına benzer). Closure'lar Nextflow operatörleriyle çalışmak için esastır.

`{ row -> return row }` closure'ı bir parametre alır: `row` (herhangi bir isim olabilir: `item`, `sample`, vb.).

`.map()` operatörü her kanal öğesini işlediğinde, o öğeyi closure'ınıza geçirir. Burada `row`, her seferinde bir CSV satırını tutar.

Bu değişikliği uygulayın ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

Girdinin değişmeden geri döndüğü için öncekiyle aynı çıktıyı göreceksiniz. Bu, map operatörünün doğru çalıştığını doğrular. Şimdi veriyi dönüştürmeye başlayalım.

#### 1.1.3. Bir Map Veri Yapısı Oluşturma

Şimdi her veri satırını dönüştürmek için closure'ımızın içinde **scripting** mantığı yazacağız. Veri akışını düzenlemek yerine bireysel veri öğelerini işlediğimiz yer burasıdır.

=== "Sonra"

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

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

`sample_meta` map'i, ilgili bilgileri depolayan bir anahtar-değer veri yapısıdır (Python'daki sözlükler, JavaScript'teki nesneler veya Ruby'deki hash'ler gibi): örnek kimliği, organizma, doku tipi, dizileme derinliği ve kalite skoru.

Verimizi temizlemek için `.toLowerCase()` ve `.replaceAll()` gibi string manipülasyon yöntemlerini, CSV'den gelen string verileri uygun sayısal türlere dönüştürmek için de `.toInteger()` ve `.toDouble()` gibi tür dönüştürme yöntemlerini kullanıyoruz.

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

Şimdi daha fazla scripting ekleyelim - bu sefer veri değerlerine göre karar vermek için üçlü operatör kullanarak.

Aşağıdaki değişikliği yapın:

=== "Sonra"

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

=== "Önce"

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

Üçlü operatör, `koşul ? doğruysa_değer : yanlışsa_değer` kalıbını izleyen bir if/else ifadesinin kısaltmasıdır. Bu satır şu anlama gelir: "Kalite 40'tan büyükse 'high' kullan, aksi takdirde 'normal' kullan". Kuzeni **Elvis operatörü** (`?:`), bir şey null veya boş olduğunda varsayılan değerler sağlar - bu kalıbı bu eğitimde daha sonra keşfedeceğiz.

Map toplama operatörü `+`, mevcut olanı değiştirmek yerine **yeni bir map** oluşturur. Bu satır, `sample_meta`daki tüm anahtar-değer çiftlerini artı yeni `priority` anahtarını içeren yeni bir map oluşturur.

!!! Note

    Closure'lara geçirilen map'leri asla değiştirmeyin - her zaman `+` kullanarak yenilerini oluşturun (örneğin). Nextflow'da, aynı veri genellikle aynı anda birden fazla işlemden geçer. Bir map'i yerinde değiştirmek, diğer işlemler aynı nesneye başvurduğunda öngörülemeyen yan etkilere neden olabilir. Yeni map'ler oluşturmak, her işlemin kendi temiz kopyasına sahip olmasını sağlar.

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

Meta verimizi kalite skorlarına göre bir öncelik seviyesiyle zenginleştirmek için başarıyla koşullu mantık ekledik.

#### 1.1.5. `.subMap()` ile Map'leri Alt Kümeleme

`+` operatörü bir map'e anahtar eklerken, bazen tersini yapmanız gerekir - sadece belirli anahtarları çıkartmak. `.subMap()` yöntemi bunun için mükemmeldir.

Meta verimizin sadece tanımlama alanlarını içeren basitleştirilmiş bir versiyonunu oluşturmak için bir satır ekleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
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
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "Sadece ID alanları: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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

    Sadece ID alanları: [id:sample_001, organism:human, tissue:liver]
    Sadece ID alanları: [id:sample_002, organism:mouse, tissue:brain]
    Sadece ID alanları: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Bu, hem `view()` işlemi tarafından görüntülenen tam meta veriyi hem de `println` ile yazdırdığımız çıkarılmış alt kümeyi gösterir.

`.subMap()` yöntemi bir anahtar listesi alır ve sadece o anahtarları içeren yeni bir map döndürür. Eğer bir anahtar orijinal map'te yoksa, sonuçta yer almaz.

Bu özellikle farklı process'ler için farklı meta veri versiyonları oluşturmanız gerektiğinde kullanışlıdır - bazıları tam meta veriye ihtiyaç duyabilirken diğerleri sadece minimal tanımlama alanlarına ihtiyaç duyar.

Şimdi ileriye dönük ihtiyacımız olmadığı için iş akışınızı önceki durumuna döndürmek üzere o println ifadelerini kaldırın.

!!! tip "Map İşlemleri Özeti"

    - **Anahtar ekle**: `map1 + [new_key: value]` - Ek anahtarlarla yeni map oluşturur
    - **Anahtar çıkart**: `map1.subMap(['key1', 'key2'])` - Sadece belirtilen anahtarlarla yeni map oluşturur
    - **Her iki işlem de yeni map'ler oluşturur** - Orijinal map'ler değişmeden kalır

#### 1.1.6. Map'leri Birleştirme ve Sonuçları Döndürme

Şimdiye kadar sadece Nextflow topluluğunun 'meta map' dediği şeyi döndürüyorduk ve o meta verilerin ilişkili olduğu dosyaları göz ardı ediyorduk. Ancak Nextflow iş akışları yazıyorsanız, muhtemelen o dosyalarla bir şeyler yapmak istersiniz.

2 elementten oluşan bir tuple içeren bir kanal yapısı çıkaralım: zenginleştirilmiş meta veri map'i ve ilgili dosya yolu. Bu, Nextflow'da process'lere veri geçirmek için yaygın bir kalıptır.

=== "Sonra"

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

=== "Önce"

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

Bu `[meta, file]` tuple yapısı, Nextflow'da hem meta veriyi hem de ilişkili dosyaları process'lere geçirmek için yaygın bir kalıptır.

!!! note

    **Map'ler ve Meta Veri**: Map'ler Nextflow'da meta veriyle çalışmanın temelidir. Meta veri map'leriyle çalışmanın daha ayrıntılı açıklaması için [Meta veriyle çalışma](./metadata.md) yan görevine bakın.

İş akışımız temel kalıbı göstermektedir: **dataflow işlemleri** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) verinin pipeline boyunca nasıl hareket ettiğini düzenlerken, `.map()` closure'ı içindeki **scripting** (map'ler `[key: value]`, string yöntemleri, tür dönüşümleri, üçlü operatörler) bireysel veri öğelerinin dönüşümünü ele alır.

### 1.2. Farklı Türleri Anlamak: Channel vs List

Şimdiye kadar çok iyi, dataflow işlemleri ile scripting arasında ayrım yapabiliyoruz. Peki ya aynı yöntem adı hem dataflow hem de scripting bağlamlarında bulunduğunda?

Mükemmel bir örnek, Nextflow standart kütüphanesinde hem kanal türleri hem de List türleri için var olan `collect` yöntemidir. List üzerindeki `collect()` yöntemi her elementi dönüştürürken, kanal üzerindeki `collect()` operatörü tüm kanal emisyonlarını tek öğeli bir kanalda toplar.

Bunu bazı örnek verilerle gösterelim, kanal `collect()` operatörünün ne yaptığını hatırlayarak başlayalım. `collect.nf` dosyasına bakın:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - birden fazla kanal emisyonunu bir araya toplar
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de toplandı)" }
```

Adımlar:

- Örnek kimliklerin bir listesini tanımla
- Her örnek kimliğini ayrı ayrı yayan `fromList()` ile bir kanal oluştur
- Akıştan geçerken her öğeyi `view()` ile yazdır
- Tüm öğeleri kanalın `collect()` operatörü ile tek bir listede topla
- Toplanan sonucu (tüm örnek kimliklerini içeren tek öğe) ikinci bir `view()` ile yazdır

Kanalın yapısını değiştirdik ama verinin kendisini değiştirmedik.

Bunu onaylamak için iş akışını çalıştırın:

```bash
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Bireysel kanal öğesi: sample_001
    Bireysel kanal öğesi: sample_002
    Bireysel kanal öğesi: sample_003
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe 1'de toplandı)
    ```

`view()` her kanal emisyonu için bir çıktı döndürür, bu yüzden bu tek çıktının 3 orijinal öğenin hepsini bir listede topladığını biliyoruz.

Şimdi List üzerindeki `collect` yöntemini iş başında görelim. Orijinal örnek kimlik listesine List'in `collect` yöntemini uygulamak için `collect.nf` dosyasını değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de toplandı)" }

    // List.collect() - her elementi dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()}'e dönüştürüldü)"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de toplandı)" }
    ```

Bu yeni kod parçasında:

- Orijinal listedeki her örnek kimliğini dönüştürmek için List'in `collect` yöntemini kullanan yeni bir `formatted_ids` değişkeni tanımla
- `println` kullanarak sonucu yazdır

Değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() sonucu: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 öğe 3'e dönüştürüldü)
    Bireysel kanal öğesi: sample_001
    Bireysel kanal öğesi: sample_002
    Bireysel kanal öğesi: sample_003
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe 1'de toplandı)
    ```

Bu sefer, bazı örnekleri dışlayan bir filtre seçtiğimiz için verinin yapısını değiştirmedik, hala listede 3 öğemiz var, ancak değiştirilmiş değerlere sahip yeni bir liste üretmek için List'in `collect` yöntemini kullanarak her öğeyi dönüştürdük. Bu, bir kanal üzerinde `map` operatörü kullanmaya benzer, ancak kanal yerine bir List veri yapısı üzerinde çalışıyor.

`collect` burada bir nokta belirtmek için kullandığımız uç bir durumdur. Önemli ders şudur: iş akışları yazarken her zaman **veri yapıları** (Listeler, Map'ler, vb.) ile **kanallar** (dataflow yapıları) arasında ayrım yapın. İşlemler aynı adları paylaşabilir ancak üzerinde çağrıldıkları türe bağlı olarak tamamen farklı davranabilir.

### 1.3. Yayılma Operatörü (`*.`) - Özellik Çıkarımı için Kısayol

List'in `collect` yöntemiyle ilgili olan yayılma operatörü (`*.`), koleksiyonlardan özellikleri çıkartmanın özlü bir yolunu sağlar. Esasen yaygın bir `collect` kalıbı için sözdizimsel şekerdir.

Bir gösterimi `collect.nf` dosyamıza ekleyelim:

=== "Sonra"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de toplandı)" }

    // List.collect() - her elementi dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()}'e dönüştürüldü)"

    // Yayılma operatörü - özlü özellik erişimi
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Yayılma operatörü sonucu: ${all_ids}"
    ```

=== "Önce"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de toplandı)" }

    // List.collect() - her elementi dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()}'e dönüştürüldü)"
    ```

Güncellenmiş iş akışını çalıştırın:

```bash title="Yayılma operatörünü test et"
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() sonucu: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 öğe 3'e dönüştürüldü)
    Yayılma operatörü sonucu: [s1, s2, s3]
    Bireysel kanal öğesi: sample_001
    Bireysel kanal öğesi: sample_002
    Bireysel kanal öğesi: sample_003
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe 1'de toplandı)
    ```

Yayılma operatörü `*.`, yaygın bir collect kalıbı için kısayoldur:

```groovy
// Bunlar eşdeğerdir:
def ids = samples*.id
def ids = samples.collect { it.id }

// Yöntem çağrılarıyla da çalışır:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Yayılma operatörü, özellikle bir nesne listesinden tek bir özelliği çıkartmanız gerektiğinde kullanışlıdır - tam `collect` closure'ını yazmaktan daha okunabilirdir.

!!! tip "Yayılma vs Collect Ne Zaman Kullanılır"

    - **Yayılma (`*.`) kullan** basit özellik erişimi için: `samples*.id`, `files*.name`
    - **Collect kullan** dönüşümler veya karmaşık mantık için: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Çıkarımlar

Bu bölümde şunları öğrendiniz:

- **Dataflow vs scripting**: Kanal operatörleri verinin pipeline'ınız boyunca nasıl aktığını düzenlerken, scripting bireysel veri öğelerini dönüştürür
- **Türleri anlamak**: Aynı yöntem adı (örneğin `collect`) üzerinde çağrıldığı türe bağlı olarak farklı davranabilir (Channel vs List)
- **Bağlam önemlidir**: Kanallarla mı (dataflow) yoksa veri yapılarıyla mı (scripting) çalıştığınızın her zaman farkında olun

Bu sınırları anlamak hata ayıklama, dokümantasyon ve bakımı yapılabilir iş akışları yazmak için esastır.

Sonrasında, gerçek dünya verileriyle başa çıkmak için gerekli olan string işleme yeteneklerine daha derinlemesine dalacağız.

---

## 2. String İşleme ve Dinamik Script Oluşturma

String işlemede ustalaşmak, kırılgan iş akışlarını sağlam pipeline'lardan ayırır. Bu bölüm karmaşık dosya adlarını ayrıştırmayı, dinamik script oluşturmayı ve değişken enterpolasyonunu kapsar.

### 2.1. Kalıp Eşleme ve Düzenli İfadeler

Biyoinformatik dosyaları genellikle meta veri kodlayan karmaşık adlandırma kurallarına sahiptir. Bunları düzenli ifadelerle kalıp eşleme kullanarak otomatik olarak çıkaralım.

`main.nf` iş akışımıza geri döneceğiz ve dosya adlarından ek örnek bilgilerini çıkartmak için bazı kalıp eşleme mantığı ekleyeceğiz. Veri setimizdeki FASTQ dosyaları `SAMPLE_001_S1_L001_R1_001.fastq.gz` gibi isimlerle Illumina tarzı adlandırma kurallarını takip eder. Bunlar gizemli görünebilir, ancak aslında örnek kimliği, şerit numarası ve okuma yönü gibi yararlı meta verileri kodlarlar. Bu isimleri ayrıştırmak için regex yeteneklerini kullanacağız.

Mevcut `main.nf` iş akışınıza aşağıdaki değişikliği yapın:

=== "Sonra"

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

=== "Önce"

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

Bu, anahtar **string işleme kavramlarını** gösterir:

1. **Düzenli ifade değişmezleri** `~/pattern/` sözdizimi kullanarak - bu ters eğik çizgileri kaçırmaya gerek kalmadan bir regex kalıbı oluşturur
2. **Kalıp eşleme** `=~` operatörü ile - bu bir string'i bir regex kalıbıyla eşleştirmeye çalışır
3. **Eşleştirici nesneler** `[0][1]`, `[0][2]`, vb. ile grupları yakalar - `[0]` tüm eşleşmeyi ifade eder, `[1]`, `[2]`, vb. parantez içindeki yakalanan gruplara işaret eder

Regex kalıbını `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` parçalayalım:

| Kalıp               | Eşleşir                                  | Yakalar                           |
| ------------------- | ---------------------------------------- | --------------------------------- |
| `^(.+)`             | Başlangıçtan örnek adı                   | Grup 1: örnek adı                 |
| `_S(\d+)`           | Örnek numarası `_S1`, `_S2`, vb.         | Grup 2: örnek numarası            |
| `_L(\d{3})`         | Şerit numarası `_L001`                   | Grup 3: şerit (3 basamak)         |
| `_(R[12])`          | Okuma yönü `_R1` veya `_R2`              | Grup 4: okuma yönü                |
| `_(\d{3})`          | Parça numarası `_001`                    | Grup 5: parça (3 basamak)         |
| `\.fastq(?:\.gz)?$` | Dosya uzantısı `.fastq` veya `.fastq.gz` | Yakalanmadı (?: yakalamama grubu) |

Bu, meta verileri otomatik olarak çıkartmak için Illumina tarzı adlandırma kurallarını ayrıştırır.

Değiştirilmiş iş akışını çalıştırın:

```bash title="Kalıp eşlemeyi test et"
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

Bu, dosya adlarından zenginleştirilmiş meta verileri gösterir.

### 2.2. Process'lerde Dinamik Script Oluşturma

Process script blokları esasen shell'e geçirilen çok satırlı string'lerdir. Girdi özelliklerine göre dinamik olarak farklı script string'leri oluşturmak için **koşullu mantık** (if/else, üçlü operatörler) kullanabilirsiniz. Bu, process tanımlarını çoğaltmadan tek uçlu ve çift uçlu dizileme okumaları gibi çeşitli girdi türlerini işlemek için esastır.

Bu kalıbı gösteren iş akışımıza bir process ekleyelim. `modules/fastp.nf` dosyasını açın ve bakın:

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

Process girdi olarak FASTQ dosyaları alır ve adaptörleri kesmek ve düşük kaliteli okumaları filtrelemek için `fastp` aracını çalıştırır. Ne yazık ki, bu process'i yazan kişi örnek veri setimizde sahip olduğumuz tek uçlu okumalara izin vermemiş. Bunu iş akışımıza ekleyelim ve ne olduğunu görelim:

İlk olarak, `main.nf` iş akışınızın en ilk satırına modülü dahil edin:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Sonra `ch_samples` kanalını `FASTP` process'ine bağlamak için `workflow` bloğunu değiştirin:

=== "Sonra"

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

=== "Önce"

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

Process'in ikinci girdi dosyası için `null` değeriyle `fastp` çalıştırmaya çalıştığını görebilirsiniz, bu da başarısız olmasına neden oluyor. Bunun nedeni, veri setimizin tek uçlu okumalar içermesi, ancak process'in çift uçlu okumaları (bir seferde iki girdi dosyası) bekleyecek şekilde sabit kodlanmış olmasıdır.

Bunu, `FASTP` process'inin `script:` bloğuna koşullu mantık ekleyerek düzeltin. Bir if/else ifadesi okuma dosya sayısını kontrol eder ve komutu buna göre ayarlar.

=== "Sonra"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Basit tek uçlu vs çift uçlu tespit
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

=== "Önce"

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

Şimdi iş akışı hem tek uçlu hem de çift uçlu okumaları zarif bir şekilde işleyebilir. Koşullu mantık girdi dosya sayısını kontrol eder ve `fastp` için uygun komutu oluşturur. Çalışıp çalışmadığını görelim:

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

İyi görünüyor! Gerçekten çalıştırılan komutları kontrol edersek (görev hash'inizi özelleştirin):

```console title="Çalıştırılan komutları kontrol et"
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

Dinamik script mantığının bir başka yaygın kullanımı [Bilim için Nextflow Genomik modülünde](../../nf4science/genomics/02_joint_calling) görülebilir. O modülde, çağrılan GATK process'i birden fazla girdi dosyası alabilir, ancak her birinin doğru bir komut satırı oluşturmak için `-V` ile öneklenmesi gerekir. Process, bir girdi dosyaları koleksiyonunu (`all_gvcfs`) doğru komut argümanlarına dönüştürmek için scripting kullanır:

```groovy title="GATK için komut satırı manipülasyonu" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Process script bloklarında scripting kullanmanın bu kalıpları son derece güçlüdür ve birçok senaryoda uygulanabilir - değişken girdi türlerini işlemekten dosya koleksiyonlarından karmaşık komut satırı argü
