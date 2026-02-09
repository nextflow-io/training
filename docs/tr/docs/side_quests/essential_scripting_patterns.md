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

Bu sefer, verinin yapısını değiştirmedik, hala listede 3 öğemiz var, ancak değiştirilmiş değerlere sahip yeni bir liste üretmek için List'in `collect` yöntemini kullanarak her öğeyi dönüştürdük. Bu, bir kanal üzerinde `map` operatörü kullanmaya benzer, ancak kanal yerine bir List veri yapısı üzerinde çalışıyor.

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

Process script bloklarında scripting kullanmanın bu kalıpları son derece güçlüdür ve birçok senaryoda uygulanabilir - değişken girdi türlerini işlemekten dosya koleksiyonlarından karmaşık komut satırı argümanları oluşturmaya kadar, process'lerinizi gerçek dünya verilerinin çeşitli gereksinimlerine gerçekten uyarlanabilir hale getirir.

### 2.3. Değişken Enterpolasyonu: Nextflow ve Shell Değişkenleri

Process script'leri Nextflow değişkenlerini, shell değişkenlerini ve komut ikamelerini karıştırır, her biri farklı enterpolasyon sözdizimi ile. Yanlış sözdizimi kullanmak hatalara neden olur. Bunları bir rapor oluşturan bir process ile keşfedelim.

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

Bu process örnek kimliği ve dosya adı ile basit bir rapor yazar. Şimdi farklı değişken türlerini karıştırmamız gerektiğinde ne olduğunu görmek için çalıştıralım.

Process'i `main.nf`'nize dahil edin ve iş akışına ekleyin:

=== "Sonra"

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

=== "Önce"

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

Şimdi iş akışını çalıştırın ve `results/reports/` içinde oluşturulan raporları kontrol edin. Her örnek hakkında temel bilgiler içermelidirler.

<!-- TODO: add the run command -->

??? success "Komut çıktısı"

    ```console
    <!-- TODO: output -->
    ```

Peki ya işlemenin ne zaman ve nerede gerçekleştiğine dair bilgi eklemek istersek? Process'i **shell** değişkenlerini ve biraz komut ikamesini kullanarak rapora mevcut kullanıcı, hostname ve tarihi dahil edecek şekilde değiştirelim:

=== "Sonra"

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

=== "Önce"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Bunu çalıştırırsanız, bir hata göreceksiniz - Nextflow `${USER}`'ı var olmayan bir Nextflow değişkeni olarak yorumlamaya çalışır.

??? failure "Komut çıktısı"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Bash'in bunu işleyebilmesi için kaçırmamız gerekiyor.

Shell değişkenlerini ve komut ikamelerini ters eğik çizgi (`\`) ile kaçırarak düzeltin:

=== "Sonra"

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

=== "Önce"

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

Şimdi çalışıyor! Ters eğik çizgi (`\`) Nextflow'a "bunu yorumlama, Bash'e geçir" der.

### Çıkarımlar

Bu bölümde **string işleme** tekniklerini öğrendiniz:

- **Dosya ayrıştırma için düzenli ifadeler**: Karmaşık dosya adlandırma kurallarından meta veri çıkartmak için `=~` operatörünü ve regex kalıplarını (`~/pattern/`) kullanma
- **Dinamik script oluşturma**: Girdi özelliklerine göre farklı script string'leri oluşturmak için koşullu mantık (if/else, üçlü operatörler) kullanma
- **Değişken enterpolasyonu**: Nextflow'un string'leri ne zaman yorumladığını ve shell'in ne zaman yorumladığını anlama
  - `${var}` - Nextflow değişkenleri (iş akışı derleme zamanında Nextflow tarafından enterpolasyon yapılır)
  - `\${var}` - Shell ortam değişkenleri (kaçırılmış, çalışma zamanında bash'e geçirilir)
  - `\$(cmd)` - Shell komut ikamesi (kaçırılmış, çalışma zamanında bash tarafından çalıştırılır)

Bu string işleme ve oluşturma kalıpları, gerçek dünya biyoinformatik iş akışlarında karşılaşacağınız çeşitli dosya formatları ve adlandırma kurallarını işlemek için esastır.

---

## 3. Yeniden Kullanılabilir Fonksiyonlar Oluşturma

Kanal operatörlerinde veya process tanımlarında satır içi karmaşık iş akışı mantığı okunabilirliği ve bakımı azaltır. **Fonksiyonlar**, bu mantığı adlandırılmış, yeniden kullanılabilir bileşenlere çıkartmanıza olanak tanır.

Map işlemimiz uzun ve karmaşık hale geldi. `def` anahtar kelimesini kullanarak bunu yeniden kullanılabilir bir fonksiyona çıkaralım.

Mevcut iş akışımızla bunun nasıl göründüğünü göstermek için, `separateMetadata` adlı yeniden kullanılabilir bir fonksiyon tanımlamak üzere `def` kullanarak aşağıdaki değişikliği yapın:

=== "Sonra"

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

=== "Önce"

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

Bu mantığı bir fonksiyona çıkartarak, gerçek iş akışı mantığını çok daha temiz bir şeye indirgedik:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Bu, iş akışı mantığını bir bakışta okumayı ve anlamayı çok daha kolay hale getirir. `separateMetadata` fonksiyonu, meta verileri ayrıştırma ve zenginleştirme için tüm karmaşık mantığı kapsüller, onu yeniden kullanılabilir ve test edilebilir hale getirir.

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

Çıktı her iki process'in de başarıyla tamamlandığını göstermelidir. İş akışı artık çok daha temiz ve bakımı kolay, tüm karmaşık meta veri işleme mantığı `separateMetadata` fonksiyonunda kapsüllenmiş durumda.

### Çıkarımlar

Bu bölümde **fonksiyon oluşturmayı** öğrendiniz:

- **`def` ile fonksiyon tanımlama**: Adlandırılmış fonksiyonlar oluşturmak için anahtar kelime (Python'daki `def` veya JavaScript'teki `function` gibi)
- **Fonksiyon kapsamı**: Script seviyesinde tanımlanan fonksiyonlar Nextflow iş akışınız boyunca erişilebilir
- **Dönüş değerleri**: Fonksiyonlar otomatik olarak son ifadeyi döndürür veya açık `return` kullanır
- **Daha temiz kod**: Karmaşık mantığı fonksiyonlara çıkartmak herhangi bir dilde temel bir yazılım mühendisliği pratiğidir

Sonrasında, dinamik kaynak tahsisi için process yönergelerinde closure'ları nasıl kullanacağımızı keşfedeceğiz.

---

## 4. Closure'larla Dinamik Kaynak Yönergeleri

Şimdiye kadar process'lerin `script` bloğunda scripting kullandık. Ancak **closure'lar** (Bölüm 1.1'de tanıtıldı) özellikle dinamik kaynak tahsisi için process yönergelerinde de inanılmaz derecede kullanışlıdır. FASTP process'imize örnek özelliklerine göre uyum sağlayan kaynak yönergeleri ekleyelim.

### 4.1. Örneğe özgü kaynak tahsisi

Şu anda FASTP process'imiz varsayılan kaynakları kullanıyor. Yüksek derinlikli örnekler için daha fazla CPU tahsis ederek daha akıllı hale getirelim. Dinamik bir `cpus` yönergesi ve statik bir `memory` yönergesi içerecek şekilde `modules/fastp.nf`'yi düzenleyin:

=== "Sonra"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Önce"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

`{ meta.depth > 40000000 ? 2 : 1 }` closure'ı **üçlü operatörü** (Bölüm 1.1'de ele alındı) kullanır ve her görev için değerlendirilir, örnek başına kaynak tahsisine olanak tanır. Yüksek derinlikli örnekler (>40M okuma) 2 CPU alırken, diğerleri 1 CPU alır.

!!! note "Yönergelerde Girdi Değişkenlerine Erişim"

    Closure, herhangi bir girdi değişkenine (burada `meta` gibi) erişebilir çünkü Nextflow bu closure'ları her görev yürütmesinin bağlamında değerlendirir.

Görev hash'lerini görmeyi kolaylaştırmak için `-ansi-log false` seçeneğiyle iş akışını tekrar çalıştırın.

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

Herhangi bir görev için CPU tahsisini görmek üzere çalıştırılan tam `docker` komutunu kontrol edebilirsiniz:

```console title="Docker komutunu kontrol et"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Şuna benzer bir şey görmelisiniz:

```bash title="docker komutu"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Bu örnekte yüksek derinlikli bir örnek olduğu için 2 CPU (`--cpu-shares 2048`) isteyen bir örnek seçtik, ancak örnek derinliğine bağlı olarak farklı CPU tahsisleri görmelisiniz. Bunu diğer görevler için de deneyin.

### 4.2. Yeniden deneme stratejileri

Bir başka güçlü kalıp, yeniden deneme stratejileri için `task.attempt` kullanmaktır. Bunun neden yararlı olduğunu göstermek için, FASTP'ye bellek tahsisini ihtiyaç duyduğundan daha aza indirerek başlayacağız. `modules/fastp.nf`'deki `memory` yönergesini `1.GB` olarak değiştirin:

=== "Sonra"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Önce"

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

Bu, process'in bellek limitlerini aştığı için öldürüldüğünü gösterir.

Bu, gerçek dünya iş akışlarında çok yaygın bir senaryodur - bazen bir görevin ne kadar belleğe ihtiyaç duyacağını çalıştırana kadar bilemezsiniz.

İş akışımızı daha sağlam hale getirmek için, her denemede bellek tahsisini artıran bir yeniden deneme stratejisi uygulayabiliriz, bir kez daha Groovy closure kullanarak. `memory` yönergesini temel belleği `task.attempt` ile çarpmak üzere değiştirin ve `errorStrategy 'retry'` ve `maxRetries 2` yönergelerini ekleyin:

=== "Sonra"

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

=== "Önce"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Şimdi process yetersiz bellek nedeniyle başarısız olursa, Nextflow daha fazla bellek ile yeniden deneyecektir:

- İlk deneme: 1 GB (task.attempt = 1)
- İkinci deneme: 2.GB (task.attempt = 2)

... ve böyle devam eder, `maxRetries` limitine kadar.

### Çıkarımlar

Closure'larla dinamik yönergeler şunları yapmanıza olanak tanır:

- Girdi özelliklerine göre kaynak tahsis etme
- Artan kaynaklarla otomatik yeniden deneme stratejileri uygulama
- Birden fazla faktörü birleştirme (meta veri, deneme sayısı, öncelikler)
- Karmaşık kaynak hesaplamaları için koşullu mantık kullanma

Bu, iş akışlarınızı hem daha verimli (aşırı tahsis etmeme) hem de daha sağlam (daha fazla kaynakla otomatik yeniden deneme) hale getirir.

---

## 5. Koşullu Mantık ve Process Kontrolü

Daha önce, kanal verilerini dönüştürmek için scripting ile `.map()` kullandık. Şimdi hangi process'lerin veriye göre çalıştırılacağını kontrol etmek için koşullu mantık kullanacağız - farklı örnek türlerine uyum sağlayan esnek iş akışları için esastır.

Nextflow'un [dataflow operatörleri](https://www.nextflow.io/docs/latest/reference/operator.html) çalışma zamanında değerlendirilen closure'lar alır, koşullu mantığın kanal içeriğine göre iş akışı kararlarını yönlendirmesini sağlar.

### 5.1. `.branch()` ile Yönlendirme

Örneğin, dizileme örneklerimizin sadece belirli bir eşiğin üzerinde kapsama sahip insan örnekleri olmaları durumunda FASTP ile kırpılması gerektiğini varsayalım. Fare örnekleri veya düşük kapsama sahip örnekler bunun yerine Trimgalore ile çalıştırılmalıdır (bu uydurulmuş bir örnek, ancak noktayı gösteriyor).

`modules/trimgalore.nf`'de basit bir Trimgalore process'i sağladık, isterseniz bakabilirsiniz, ancak detaylar bu alıştırma için önemli değil. Önemli nokta, örnekleri meta verilerine göre yönlendirmek istememizdir.

`modules/trimgalore.nf`'den yeni formu dahil edin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... ve sonra `main.nf` iş akışınızı örnekleri meta verilerine göre dallandırmak ve uygun kırpma process'inden yönlendirmek için şu şekilde değiştirin:

=== "Sonra"

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

=== "Önce"

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

Burada, örnekleri meta verilerine göre yönlendirmek için `.branch{}` operatörü içinde küçük ama güçlü koşullu ifadeler kullandık. Yüksek kapsama sahip insan örnekleri `FASTP`'den geçerken, diğer tüm örnekler `TRIMGALORE`'den geçer.

### 5.2. Doğruluk ile `.filter()` Kullanma

İş akışı yürütmesini kontrol etmek için bir başka güçlü kalıp, hangi öğelerin pipeline'da devam etmesi gerektiğini belirlemek için bir closure kullanan `.filter()` operatörüdür. Filtre closure'ı içinde, hangi öğelerin geçeceğine karar veren **boolean ifadeleri** yazacaksınız.

Nextflow'un (birçok dinamik dil gibi) boolean bağlamlarda hangi değerlerin `true` veya `false` olarak değerlendirileceğini belirleyen bir **"doğruluk"** kavramı vardır:

- **Doğru**: Null olmayan değerler, boş olmayan string'ler, sıfır olmayan sayılar, boş olmayan koleksiyonlar
- **Yanlış**: `null`, boş string'ler `""`, sıfır `0`, boş koleksiyonlar `[]` veya `[:]`, `false`

Bu, `meta.id`'nin tek başına (açık `!= null` olmadan) kimliğin var olup olmadığını ve boş olmadığını kontrol ettiği anlamına gelir. Kalite gereksinimlerimizi karşılamayan örnekleri filtrelemek için bunu kullanalım.

Dal işleminden önce aşağıdakini ekleyin:

=== "Sonra"

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

=== "Önce"

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

Bazı örnekleri dışlayan bir filtre seçtiğimiz için, daha az görev yürütüldü.

Filtre ifadesi `meta.id && meta.organism && meta.depth >= 25000000` doğruluğu açık karşılaştırmalarla birleştirir:

- `meta.id && meta.organism` her iki alanın da var olduğunu ve boş olmadığını kontrol eder (doğruluk kullanarak)
- `meta.depth >= 25000000` açık bir karşılaştırma ile yeterli dizileme derinliğini sağlar

!!! note "Pratikte Doğruluk"

    `meta.id && meta.organism` ifadesi şunu yazmaktan daha özlüdür:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Bu, filtreleme mantığını çok daha temiz ve okunması kolay hale getirir.

### Çıkarımlar

Bu bölümde, `.branch{}` ve `.filter{}` gibi Nextflow operatörlerinin closure arayüzlerini kullanarak iş akışı yürütmesini kontrol etmek için koşullu mantık kullanmayı, özlü koşullu ifadeler yazmak için doğruluğu kullanmayı öğrendiniz.

Pipeline'ımız artık örnekleri uygun process'lerden akıllıca yönlendiriyor, ancak üretim iş akışlarının geçersiz verileri zarif bir şekilde işlemesi gerekir. İş akışımızı eksik veya null değerlere karşı sağlam hale getirelim.

---

## 6. Güvenli Navigasyon ve Elvis Operatörleri

`separateMetadata` fonksiyonumuz şu anda tüm CSV alanlarının mevcut ve geçerli olduğunu varsayıyor. Peki eksik verilerle ne olur? Öğrenelim.

### 6.1. Problem: Var Olmayan Özelliklere Erişim

Diyelim ki isteğe bağlı dizileme çalıştırma bilgisi için destek eklemek istiyoruz. Bazı laboratuvarlarda, örneklerin dizileme çalıştırma kimliği veya parti numarası için ek bir alanı olabilir, ancak mevcut CSV'mizde bu sütun yok. Yine de ona erişmeyi deneyelim.

`separateMetadata` fonksiyonunu bir run_id alanı içerecek şekilde değiştirin:

=== "Sonra"

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

=== "Önce"

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

Bu bir NullPointerException ile çöküyor.

Problem, `row.run_id`'nin `null` döndürmesidir çünkü `run_id` sütunu CSV'mizde mevcut değildir. `null` üzerinde `.toUpperCase()` çağırmaya çalıştığımızda çöküyor. Güvenli navigasyon operatörünün günü kurtardığı yer burasıdır.

### 6.2. Güvenli Navigasyon Operatörü (`?.`)

Güvenli navigasyon operatörü (`?.`), `null` bir değer üzerinde çağrıldığında istisna atmak yerine `null` döndürür. `?.` öncesindeki nesne `null` ise, tüm ifade yöntemi çalıştırmadan `null` olarak değerlendirilir.

Güvenli navigasyon kullanmak için fonksiyonu güncelleyin:

=== "Sonra"

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

=== "Önce"

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

Çökme yok! İş akışı artık eksik alanı zarif bir şekilde işliyor. `row.run_id` `null` olduğunda, `?.` operatörü `.toUpperCase()` çağrısını önler ve `run_id` istisna oluşturmak yerine `null` olur.

### 6.3. Varsayılanlar için Elvis Operatörü (`?:`)

Elvis operatörü (`?:`), sol taraf "yanlış" olduğunda (daha önce açıklandığı gibi) varsayılan değerler sağlar. Yan tarafından bakıldığında `?:` Elvis Presley'nin ünlü saçı ve gözleri gibi göründüğü için Elvis operatörü olarak adlandırılır!

Artık güvenli navigasyon kullandığımıza göre, `run_id` o alan olmayan örnekler için `null` olacaktır. Varsayılan bir değer sağlamak ve `sample_meta` map'imize eklemek için Elvis operatörünü kullanalım:

=== "Sonra"

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

=== "Önce"

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

Ayrıca sonuçları görmek için iş akışına bir `view()` operatörü ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Önce"

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

Mükemmel! Şimdi tüm örneklerin gerçek çalıştırma kimliği (büyük harfle) veya varsayılan değer 'UNSPECIFIED' ile bir `run` alanı var. `?.` ve `?:` kombinasyonu hem güvenlik (çökme yok) hem de mantıklı varsayılanlar sağlar.

Çalıştığını onayladığımıza göre `.view()` operatörünü çıkarın.

!!! tip "Güvenli Navigasyon ve Elvis'i Birleştirme"

    `value?.method() ?: 'default'` kalıbı üretim iş akışlarında yaygındır:

    - `value?.method()` - Yöntemi güvenli bir şekilde çağırır, `value` `null` ise `null` döndürür
    - `?: 'default'` - Sonuç `null` ise yedek sağlar

    Bu kalıp eksik/eksik verileri zarif bir şekilde işler.

Bu operatörleri fonksiyonlarda, operatör closure'larında (`.map{}`, `.filter{}`), process script'lerinde ve yapılandırma dosyalarında tutarlı bir şekilde kullanın. Gerçek dünya verileriyle çalışırken çökmeleri önlerler.

### Çıkarımlar

- **Güvenli navigasyon (`?.`)**: Null değerlerde çökmeleri önler - istisna atmak yerine null döndürür
- **Elvis operatörü (`?:`)**: Varsayılanlar sağlar - `value ?: 'default'`
- **Birleştirme**: `value?.method() ?: 'default'` yaygın kalıptır

Bu operatörler iş akışlarını eksik verilere karşı dirençli hale getirir - gerçek dünya çalışması için esastır.

---

## 7. `error()` ve `log.warn` ile Doğrulama

Bazen girdi parametreleri geçersizse iş akışını hemen durdurmanız gerekir. Nextflow'da, doğrulama mantığı uygulamak için `error()` ve `log.warn` gibi yerleşik fonksiyonların yanı sıra `if` ifadeleri ve boolean mantığı gibi standart programlama yapılarını kullanabilirsiniz. İş akışımıza doğrulama ekleyelim.

İş akışı bloğunuzdan önce bir doğrulama fonksiyonu oluşturun, iş akışından çağırın ve kanal oluşturmayı CSV dosya yolu için bir parametre kullanacak şekilde değiştirin. Parametre eksikse veya dosya mevcut değilse, açık bir mesajla yürütmeyi durdurmak için `error()` çağırın.

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Girdi parametresinin sağlandığını kontrol et
        if (!params.input) {
            error("Girdi CSV dosya yolu sağlanmadı. Lütfen --input <file.csv> belirtin")
        }

        // CSV dosyasının var olduğunu kontrol et
        if (!file(params.input).exists()) {
            error("Girdi CSV dosyası bulunamadı: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Önce"

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
    Girdi CSV dosya yolu sağlanmadı. Lütfen --input <file.csv> belirtin
    ```

İş akışı daha sonra gizemli bir şekilde başarısız olmak yerine açık bir hata mesajıyla hemen durur

Şimdi var olmayan bir dosyayla çalıştırın:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Girdi CSV dosyası bulunamadı: ./data/nonexistent.csv
    ```

Son olarak, doğru dosyayla çalıştırın:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Komut çıktısı"

    ```console
    <!-- TODO: output -->
    ```

Bu sefer başarıyla çalışıyor.

`separateMetadata` fonksiyonu içinde de doğrulama ekleyebilirsiniz. Düşük dizileme derinliğine sahip örnekler için uyarılar vermek üzere ölümcül olmayan `log.warn` kullanalım, ancak iş akışının devam etmesine izin verelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Verinin mantıklı olduğunu doğrula
        if (sample_meta.depth < 30000000) {
            log.warn "${sample_meta.id} için düşük dizileme derinliği: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Önce"

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
    WARN: sample_002 için düşük dizileme derinliği: 25000000
    ```

Örneklerden biri için düşük dizileme derinliği hakkında bir uyarı görüyoruz.

### Çıkarımlar

- **`error()`**: Açık mesajla iş akışını hemen durdurur
- **`log.warn`**: İş akışını durdurmadan uyarılar verir
- **Erken doğrulama**: Yararlı hatalarla hızlı başarısızlık için işlemeden önce girdileri kontrol edin
- **Doğrulama fonksiyonları**: İş akışı başlangıcında çağrılabilecek yeniden kullanılabilir doğrulama mantığı oluşturun

Uygun doğrulama, sorunları erken yakalayarak ve açık hata mesajlarıyla iş akışlarını daha sağlam ve kullanıcı dostu hale getirir.

---

## 8. Workflow Olay İşleyicileri

Şimdiye kadar, iş akışı script'lerimizde ve process tanımlarımızda kod yazıyorduk. Ancak bilmeniz gereken bir önemli özellik daha var: workflow olay işleyicileri.

Olay işleyicileri, iş akışınızın yaşam döngüsünde belirli noktalarda çalışan closure'lardır. Günlükleme, bildirimler veya temizleme işlemleri eklemek için mükemmeldir. Bu işleyiciler, iş akışı tanımınızın yanında iş akışı script'inizde tanımlanmalıdır.

### 8.1. `onComplete` İşleyicisi

En yaygın kullanılan olay işleyicisi, iş akışınız bittiğinde (başarılı veya başarısız olsun) çalışan `onComplete`'dir. Pipeline sonuçlarımızı özetlemek için bir tane ekleyelim.

Olay işleyicisini `main.nf` dosyanıza, iş akışı tanımınızın içine ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline yürütme özeti:"
            println "=========================="
            println "Tamamlandı: ${workflow.complete}"
            println "Süre      : ${workflow.duration}"
            println "Başarı    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "çıkış durumu : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Bu closure iş akışı tamamlandığında çalışır. İçinde, yürütme hakkında yararlı özellikler sağlayan `workflow` nesnesine erişiminiz vardır.

İş akışınızı çalıştırın ve bu özetin sonunda göründüğünü göreceksiniz!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: sample_002 için düşük dizileme derinliği: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline yürütme özeti:
    ==========================
    Tamamlandı: 2025-10-10T12:14:24.885384+01:00
    Süre      : 2.9s
    Başarı    : true
    workDir   : /workspaces/training/side-quests/essential_scripting_patterns/work
    çıkış durumu : 0
    ```

Koşullu mantık ekleyerek daha kullanışlı hale getirelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline yürütme özeti:"
            println "=========================="
            println "Tamamlandı: ${workflow.complete}"
            println "Süre      : ${workflow.duration}"
            println "Başarı    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "çıkış durumu : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline başarıyla tamamlandı!"
            } else {
                println "❌ Pipeline başarısız oldu!"
                println "Hata: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline yürütme özeti:"
            println "=========================="
            println "Tamamlandı: ${workflow.complete}"
            println "Süre      : ${workflow.duration}"
            println "Başarı    : ${workflow.success}"
            println "workDir   : ${workflow.workDir}"
            println "çıkış durumu : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Şimdi belirtilmişse başarı/başarısızlık mesajı ve çıktı dizini dahil olmak üzere daha da bilgilendirici bir özet alıyoruz:

<!-- TODO: add run command -->

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: sample_002 için düşük dizileme derinliği:25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline yürütme özeti:
    ==========================
    Tamamlandı: 2025-10-10T12:16:00.522569+01:00
    Süre      : 3.6s
    Başarı    : true
    workDir   : /workspaces/training/side-quests/essential_scripting_patterns/work
    çıkış durumu : 0

    ✅ Pipeline başarıyla tamamlandı!
    ```

Dosya işlemlerini kullanarak özeti bir dosyaya da yazabilirsiniz:

```groovy title="main.nf - Özeti dosyaya yazma"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onComplete = {
        def summary = """
        Pipeline Yürütme Özeti
        ===========================
        Tamamlandı: ${workflow.complete}
        Süre      : ${workflow.duration}
        Başarı    : ${workflow.success}
        Komut     : ${workflow.commandLine}
        """

        println summary

        // Bir log dosyasına yaz
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. `onError` İşleyicisi

`onComplete`'in yanı sıra, kullanabileceğiniz bir başka olay işleyicisi daha var: `onError`, sadece iş akışı başarısız olursa çalışır:

```groovy title="main.nf - onError işleyicisi"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "="* 50
        println "Pipeline yürütmesi başarısız oldu!"
        println "Hata mesajı: ${workflow.errorMessage}"
        println "="* 50

        // Detaylı hata günlüğü yaz
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Hata Raporu
        =====================
        Zaman: ${new Date()}
        Hata: ${workflow.errorMessage}
        Hata raporu: ${workflow.errorReport ?: 'Detaylı rapor mevcut değil'}
        """

        println "Hata detayları şuraya yazıldı: ${error_file}"
    }
}
```

İş akışı script'inizde birden fazla işleyiciyi birlikte kullanabilirsiniz:

```groovy title="main.nf - Birleştirilmiş işleyiciler"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "Workflow başarısız oldu: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "BAŞARI ✅" : "BAŞARISIZ ❌"

        println """
        Pipeline tamamlandı: ${status}
        Süre: ${duration_mins} dakika
        """
    }
}
```

### Çıkarımlar

Bu bölümde şunları öğrendiniz:

- **Olay işleyici closure'ları**: Farklı yaşam döngüsü noktalarında çalışan iş akışı script'inizdeki closure'lar
- **`onComplete` işleyicisi**: Yürütme özetleri ve sonuç raporlaması için
- **`onError` işleyicisi**: Hata işleme ve başarısızlıkları günlükleme için
- **Workflow nesne özellikleri**: `workflow.success`, `workflow.duration`, `workflow.errorMessage`, vb.'ye erişim

Olay işleyicileri, gelişmiş günlükleme ve bildirim yetenekleri eklemek için iş akışı script'lerinizde Nextflow dilinin tam gücünü nasıl kullanabileceğinizi gösterir.

---

## Özet

Tebrikler, başardınız!

Bu yan görev boyunca, temel meta veri işlemeden gelişmiş, üretime hazır bir iş akışına evrilen kapsamlı bir örnek işleme pipeline'ı oluşturdunuz.
Her bölüm bir öncekinin üzerine inşa edildi, programlama yapılarının basit iş akışlarını güçlü veri işleme sistemlerine nasıl dönüştürdüğünü gösterdi, aşağıdaki faydalarla:

- **Daha açık kod**: Dataflow vs scripting'i anlamak daha organize iş akışları yazmanıza yardımcı olur
- **Sağlam işleme**: Güvenli navigasyon ve Elvis operatörleri iş akışlarını eksik verilere karşı dirençli hale getirir
- **Esnek işleme**: Koşullu mantık iş akışlarınızın farklı örnek türlerini uygun şekilde işlemesini sağlar
- **Uyarlanabilir kaynaklar**: Dinamik yönergeler girdi özelliklerine göre kaynak kullanımını optimize eder

Bu ilerleme, birkaç örneği işleyen araştırma prototiplerinden laboratuvarlar ve kurumlar arasında binlerce örneği işleyen üretim sistemlerine kadar gerçek dünya biyoinformatik pipeline'larının evrimini yansıtır.
Çözdüğünüz her zorluk ve öğrendiğiniz her kalıp, geliştiricilerin Nextflow iş akışlarını ölçeklendirirken karşılaştıkları gerçek sorunları yansıtır.

Bu kalıpları kendi çalışmanızda uygulamak, sağlam, üretime hazır iş akışları oluşturmanızı sağlayacaktır.

### Anahtar kalıplar

1.  **Dataflow vs Scripting:** Dataflow işlemleri (kanal düzenlemesi) ile scripting (veriyi manipüle eden kod) arasında ayrım yapmayı öğrendiniz, Channel vs List üzerindeki `collect` gibi farklı türlerdeki işlemler arasındaki önemli farkları da içerecek şekilde.

    - Dataflow: kanal düzenlemesi

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: koleksiyonlar üzerinde veri işleme

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Gelişmiş String İşleme**: Dosya adlarını ayrıştırmak için düzenli ifadelerde, process'lerde dinamik script oluşturmada ve değişken enterpolasyonunda (Nextflow vs Bash vs Shell) ustalaştınız.

    - Kalıp eşleme

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Koşullu dönüşle fonksiyon

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Dosya koleksiyonunu komut argümanlarına (process script bloğunda)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Yeniden Kullanılabilir Fonksiyonlar Oluşturma**: Karmaşık mantığı kanal operatörlerinden çağrılabilen adlandırılmış fonksiyonlara çıkartmayı öğrendiniz, iş akışlarını daha okunabilir ve bakımı kolay hale getirdiniz.

    - Adlandırılmış bir fonksiyon tanımla

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

    - Adlandırılmış fonksiyonu bir iş akışında çağır

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closure'larla Dinamik Kaynak Yönergeleri**: Girdi özelliklerine göre uyarlanabilir kaynak tahsisi için process yönergelerinde closure'ları kullanmayı keşfettiniz.

    - Adlandırılmış closure'lar ve kompozisyon

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Kapsam erişimi olan closure'lar

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Koşullu Mantık ve Process Kontrolü**: `.branch()` ve `.filter()` operatörlerini kullanarak akıllı yönlendirme eklediniz, özlü koşullu ifadeler için doğruluğu kullandınız.

    - Veriyi farklı iş akışı dallarından yönlendirmek için `.branch()` kullan

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
    if (sample.files) println "Dosyalar var"
    ```

    - Veriyi 'doğruluk' ile alt kümelemek için `filter()` kullan

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Güvenli Navigasyon ve Elvis Operatörleri**: Null-safe özellik erişimi için `?.` ve varsayılan değerler sağlamak için `?:` kullanarak pipeline'ı eksik verilere karşı sağlam hale getirdiniz.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **error() ve log.warn ile Doğrulama**: Girdileri erken doğrulamayı ve açık hata mesajlarıyla hızlı başarısızlık sağlamayı öğrendiniz.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Geçersiz: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Hata: ${e.message}"
    }
    ```

8.  **Yapılandırma Olay İşleyicileri**: Günlükleme, bildirimler ve yaşam döngüsü yönetimi için workflow olay işleyicilerini (`onComplete` ve `onError`) kullanmayı öğrendiniz.

    - Günlükleme ve bildirim için `onComplete` kullanma

    ```groovy
    workflow.onComplete = {
        println "Başarı     : ${workflow.success}"
        println "çıkış durumu : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline başarıyla tamamlandı!"
        } else {
            println "❌ Pipeline başarısız oldu!"
            println "Hata: ${workflow.errorMessage}"
        }
    }
    ```

    - Özellikle başarısızlık durumunda işlem yapmak için `onError` kullanma

    ```groovy
    workflow.onError = {
        // Detaylı hata günlüğü yaz
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Zaman: ${new Date()}
        Hata: ${workflow.errorMessage}
        Hata raporu: ${workflow.errorReport ?: 'Detaylı rapor mevcut değil'}
        """

        println "Hata detayları şuraya yazıldı: ${error_file}"
    }
    ```

### Ek kaynaklar

- [Nextflow Dil Referansı](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operatörleri](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Sözdizimi](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standart Kütüphanesi](https://nextflow.io/docs/latest/reference/stdlib.html)

Daha gelişmiş özellikleri keşfetmeniz gerektiğinde bu kaynaklara mutlaka göz atın.

Becerilerinizi pratik yaparak ve genişleterek geliştirmekten faydalanacaksınız:

- Dataflow ve scripting arasında uygun ayrım ile daha temiz iş akışları yazma
- Nextflow, Bash ve shell değişkenleriyle yaygın tuzaklardan kaçınmak için değişken enterpolasyonunda ustalaşma
- Verimli, uyarlanabilir iş akışları için dinamik kaynak yönergeleri kullanma
- Dosya koleksiyonlarını düzgün biçimlendirilmiş komut satırı argümanlarına dönüştürme
- Regex ve string işleme kullanarak farklı dosya adlandırma kurallarını ve girdi formatlarını zarif bir şekilde işleme
- Gelişmiş closure kalıpları ve fonksiyonel programlama kullanarak yeniden kullanılabilir, bakımı yapılabilir kod oluşturma
- Koleksiyon işlemlerini kullanarak karmaşık veri setlerini işleme ve organize etme
- İş akışlarınızı üretime hazır hale getirmek için doğrulama, hata işleme ve günlükleme ekleme
- Olay işleyicileriyle workflow yaşam döngüsü yönetimi uygulama

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ altındaki düğmeye tıklayın.
