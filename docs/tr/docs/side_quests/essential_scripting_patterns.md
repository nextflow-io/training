# Temel Nextflow Betik Kalıpları

Nextflow, Java Sanal Makinesi üzerinde çalışan bir programlama dilidir. Nextflow [Groovy](http://groovy-lang.org/) üzerine inşa edilmiş ve sözdiziminin çoğunu paylaşsa da, Nextflow sadece "uzantılı Groovy"den daha fazlasıdır -- tam olarak belirtilmiş bir [sözdizimi](https://nextflow.io/docs/latest/reference/syntax.html) ve [standart kütüphane](https://nextflow.io/docs/latest/reference/stdlib.html) ile bağımsız bir dildir.

Değişkenler, map'ler ve listeler için temel sözdiziminin ötesine geçmeden çok sayıda Nextflow yazabilirsiniz. Çoğu Nextflow eğitimi iş akışı düzenlemesine (kanallar, süreçler ve veri akışı) odaklanır ve sadece bununla şaşırtıcı derecede ileri gidebilirsiniz.

Ancak, veri manipüle etmeniz, karmaşık dosya adlarını ayrıştırmanız, koşullu mantık uygulamanız veya sağlam üretim iş akışları oluşturmanız gerektiğinde, kodunuzun iki farklı yönü hakkında düşünmek yardımcı olur: **veri akışı** (kanallar, operatörler, süreçler ve iş akışları) ve **betik yazma** (closure'lar, fonksiyonlar ve süreç betikleri içindeki kod). Bu ayrım biraz keyfi olsa da --hepsi Nextflow kodu-- pipeline'ınızı ne zaman düzenlediğiniz ile ne zaman veri manipüle ettiğinizi anlamak için yararlı bir zihinsel model sağlar. Her ikisinde de ustalaşmak, açık ve sürdürülebilir iş akışları yazma yeteneğinizi önemli ölçüde geliştirir.

### Öğrenme hedefleri

Bu yan görev, sizi temel kavramlardan üretime hazır kalıplara kadar uygulamalı bir yolculuğa çıkarır.
Basit bir CSV okuma iş akışını, gerçekçi zorluklar aracılığıyla adım adım geliştirerek sofistike bir biyoinformatik pipeline'ına dönüştüreceğiz:

- **Sınırları anlama:** Veri akışı işlemleri ile betik yazma arasında ayrım yapın ve bunların birlikte nasıl çalıştığını anlayın
- **Veri manipülasyonu:** Güçlü operatörler kullanarak map'leri ve koleksiyonları çıkarın, dönüştürün ve alt kümelere ayırın
- **String işleme:** Regex kalıpları ile karmaşık dosya adlandırma şemalarını ayrıştırın ve değişken interpolasyonunda ustalaşın
- **Yeniden kullanılabilir fonksiyonlar:** Daha temiz, daha sürdürülebilir iş akışları için karmaşık mantığı adlandırılmış fonksiyonlara çıkarın
- **Dinamik mantık:** Farklı girdi türlerine uyum sağlayan süreçler oluşturun ve dinamik kaynak tahsisi için closure'lar kullanın
- **Koşullu yönlendirme:** Örnekleri meta veri özelliklerine göre farklı süreçlerden akıllıca geçirin
- **Güvenli işlemler:** Null-safe operatörlerle eksik verileri zarif bir şekilde ele alın ve açık hata mesajlarıyla girdileri doğrulayın
- **Yapılandırma tabanlı işleyiciler:** Loglama, bildirimler ve yaşam döngüsü yönetimi için iş akışı olay işleyicilerini kullanın

### Ön koşullar

Bu yan görevi üstlenmeden önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, dosyalarla çalışma, meta veriler) kullanmakta rahat olmalısınız
- Yaygın programlama yapılarına (değişkenler, map'ler, listeler) temel düzeyde aşina olmalısınız

Bu eğitim, programlama kavramlarını karşılaştıkça açıklayacaktır, bu nedenle kapsamlı programlama deneyimine ihtiyacınız yoktur.
Temel kavramlarla başlayıp gelişmiş kalıplara doğru ilerleyeceğiz.

---

## 0. Başlayın

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

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

Örnek CSV'miz, özelliklerine göre farklı işleme ihtiyaç duyan biyolojik örnekler hakkında bilgi içerir:

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

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
<!-- - [ ] I understand the assignment -->

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Veri Akışı vs Betik Yazma: Sınırları Anlamak

### 1.1. Neyin Ne Olduğunu Belirleme

Nextflow iş akışları yazarken, **veri akışı** (verilerin kanallar ve süreçler arasında nasıl hareket ettiği) ile **betik yazma** (verileri manipüle eden ve kararlar veren kod) arasında ayrım yapmak önemlidir. Bunların birlikte nasıl çalıştığını gösteren bir iş akışı oluşturalım.

#### 1.1.1. Temel Nextflow İş Akışı

Sadece CSV dosyasını okuyan basit bir iş akışıyla başlayın (bunu sizin için `main.nf` dosyasında zaten yaptık):

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

Şimdi verileri dönüştürmek için betik yazma ekleyeceğiz, muhtemelen zaten aşina olduğunuz `.map()` operatörünü kullanarak. Bu operatör, her öğeyi dönüştürmek için kod yazabileceğimiz bir 'closure' alır.

!!! note

    **Closure**, etrafta taşınabilen ve daha sonra çalıştırılabilen bir kod bloğudur. Bunu satır içi tanımladığınız bir fonksiyon olarak düşünün. Closure'lar süslü parantezlerle `{ }` yazılır ve parametre alabilir. Nextflow operatörlerinin nasıl çalıştığının temelini oluştururlar ve bir süredir Nextflow yazıyorsanız, farkında olmadan zaten kullanıyor olabilirsiniz!

İşte o map işlemi şöyle görünür:

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

Bu bizim ilk **closure**'ımız - argüman olarak geçirebileceğiniz anonim bir fonksiyon (Python'daki lambda'lara veya JavaScript'teki ok fonksiyonlarına benzer). Closure'lar Nextflow operatörleriyle çalışmak için esastır.

`{ row -> return row }` closure'ı bir `row` parametresi alır (herhangi bir isim olabilir: `item`, `sample`, vb.).

`.map()` operatörü her kanal öğesini işlediğinde, o öğeyi closure'ınıza geçirir. Burada, `row` her seferinde bir CSV satırını tutar.

Bu değişikliği uygulayın ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

Girdinin değişmeden döndürüldüğü için öncekiyle aynı çıktıyı göreceksiniz. Bu, map operatörünün doğru çalıştığını doğrular. Şimdi verileri dönüştürmeye başlayalım.

#### 1.1.3. Bir Map Veri Yapısı Oluşturma

Şimdi her veri satırını dönüştürmek için closure'ımızın içine **betik yazma** mantığı yazacağız. Burada veri akışını düzenlemek yerine bireysel veri öğelerini işliyoruz.

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Veri dönüşümü için betik yazma
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

`sample_meta` map'i, ilgili bilgileri saklayan bir anahtar-değer veri yapısıdır (Python'daki sözlükler, JavaScript'teki nesneler veya Ruby'deki hash'ler gibi): örnek ID, organizma, doku tipi, dizileme derinliği ve kalite skoru.

Verilerimizi temizlemek için `.toLowerCase()` ve `.replaceAll()` gibi string manipülasyon metodlarını ve CSV'den gelen string verileri uygun sayısal türlere dönüştürmek için `.toInteger()` ve `.toDouble()` gibi tür dönüştürme metodlarını kullanıyoruz.

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

Şimdi daha fazla betik yazma ekleyelim - bu sefer veri değerlerine göre kararlar almak için üçlü operatör kullanarak.

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

Üçlü operatör, `koşul ? doğruysa_değer : yanlışsa_değer` kalıbını izleyen bir if/else ifadesinin kısaltmasıdır. Bu satır şu anlama gelir: "Kalite 40'tan büyükse, 'high' kullan, aksi takdirde 'normal' kullan". Kuzeni olan **Elvis operatörü** (`?:`), bir şey null veya boş olduğunda varsayılan değerler sağlar - bu kalıbı bu eğitimde daha sonra keşfedeceğiz.

Map toplama operatörü `+`, mevcut olanı değiştirmek yerine **yeni bir map** oluşturur. Bu satır, `sample_meta`'dan tüm anahtar-değer çiftlerini artı yeni `priority` anahtarını içeren yeni bir map oluşturur.

!!! Note

    Closure'lara geçirilen map'leri asla değiştirmeyin - her zaman `+` kullanarak (örneğin) yenilerini oluşturun. Nextflow'da, aynı veriler genellikle aynı anda birden fazla işlemden geçer. Bir map'i yerinde değiştirmek, diğer işlemler aynı nesneye referans verdiğinde öngörülemeyen yan etkilere neden olabilir. Yeni map'ler oluşturmak, her işlemin kendi temiz kopyasına sahip olmasını sağlar.

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

Kalite skorlarına göre bir öncelik seviyesi ile meta verilerimizi zenginleştirmek için başarıyla koşullu mantık ekledik.

#### 1.1.5. `.subMap()` ile Map'leri Alt Kümelere Ayırma

`+` operatörü bir map'e anahtarlar eklerken, bazen tersini yapmanız gerekir - yalnızca belirli anahtarları çıkarmak. `.subMap()` metodu bunun için mükemmeldir.

Meta verilerimizin yalnızca tanımlama alanlarını içeren basitleştirilmiş bir versiyonunu oluşturmak için bir satır ekleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Veri dönüşümü için betik yazma
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
                // Veri dönüşümü için betik yazma
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

Bu, hem `view()` işlemi tarafından görüntülenen tam meta verileri hem de `println` ile yazdırdığımız çıkarılmış alt kümeyi gösterir.

`.subMap()` metodu bir anahtar listesi alır ve yalnızca o anahtarları içeren yeni bir map döndürür. Bir anahtar orijinal map'te yoksa, sonuca dahil edilmez.

Bu, farklı süreçler için farklı meta veri versiyonları oluşturmanız gerektiğinde özellikle yararlıdır - bazıları tam meta verilere ihtiyaç duyarken diğerleri yalnızca minimal tanımlama alanlarına ihtiyaç duyabilir.

Şimdi ileriye dönük olarak bunlara ihtiyacımız olmadığı için iş akışınızı önceki durumuna geri yüklemek için o println ifadelerini kaldırın.

!!! tip "Map İşlemleri Özeti"

    - **Anahtar ekleme**: `map1 + [new_key: value]` - Ek anahtarlarla yeni map oluşturur
    - **Anahtar çıkarma**: `map1.subMap(['key1', 'key2'])` - Yalnızca belirtilen anahtarlarla yeni map oluşturur
    - **Her iki işlem de yeni map'ler oluşturur** - Orijinal map'ler değişmeden kalır

#### 1.1.6. Map'leri Birleştirme ve Sonuçları Döndürme

Şimdiye kadar, yalnızca Nextflow topluluğunun 'meta map' dediği şeyi döndürüyorduk ve bu meta verilerin ilişkili olduğu dosyaları görmezden geliyorduk. Ancak Nextflow iş akışları yazıyorsanız, muhtemelen bu dosyalarla bir şeyler yapmak istersiniz.

2 elemanlı bir tuple'dan oluşan bir kanal yapısı çıkaralım: zenginleştirilmiş meta veri map'i ve karşılık gelen dosya yolu. Bu, Nextflow'da süreçlere veri geçirmek için yaygın bir kalıptır.

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

Bu `[meta, file]` tuple yapısı, Nextflow'da hem meta verileri hem de ilişkili dosyaları süreçlere geçirmek için yaygın bir kalıptır.

!!! note

    **Map'ler ve Meta Veriler**: Map'ler Nextflow'da meta verilerle çalışmanın temelidir. Meta veri map'leriyle çalışma hakkında daha ayrıntılı bir açıklama için [Meta verilerle çalışma](./metadata.md) yan görevine bakın.

İş akışımız temel kalıbı gösterir: **veri akışı işlemleri** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) verilerin pipeline boyunca nasıl hareket ettiğini düzenlerken, `.map()` closure'ı içindeki **betik yazma** (map'ler `[key: value]`, string metodları, tür dönüşümleri, üçlü operatörler) bireysel veri öğelerinin dönüşümünü ele alır.

### 1.2. Farklı Türleri Anlama: Kanal vs Liste

Şimdiye kadar çok iyi, veri akışı işlemleri ile betik yazma arasında ayrım yapabiliriz. Peki aynı metod adı hem kanal türleri hem de Liste türleri için Nextflow standart kütüphanesinde mevcut olduğunda ne olur?

Mükemmel bir örnek, hem kanal türleri hem de Liste türleri için var olan `collect` metodudur. Bir Liste üzerindeki `collect()` metodu her elemanı dönüştürürken, bir kanal üzerindeki `collect()` operatörü tüm kanal emisyonlarını tek öğeli bir kanalda toplar.

Bunu bazı örnek verilerle gösterelim, kanal `collect()` operatörünün ne yaptığını hatırlayarak başlayalım. `collect.nf` dosyasına bakın:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - birden fazla kanal emisyonunu bir araya toplar
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de gruplandı)" }
```

Adımlar:

- Örnek ID'lerden oluşan bir Liste tanımlayın
- Her örnek ID'yi ayrı ayrı yayan `fromList()` ile bir kanal oluşturun
- Akıştan geçerken her öğeyi `view()` ile yazdırın
- Tüm öğeleri kanalın `collect()` operatörü ile tek bir listede toplayın
- Toplanan sonucu (tüm örnek ID'leri içeren tek öğe) ikinci bir `view()` ile yazdırın

Kanalın yapısını değiştirdik, ancak verinin kendisini değiştirmedik.

Bunu doğrulamak için iş akışını çalıştırın:

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
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe 1'de gruplandı)
    ```

`view()` her kanal emisyonu için bir çıktı döndürür, bu nedenle bu tek çıktının 3 orijinal öğenin tümünü bir listede gruplandırdığını biliyoruz.

Şimdi bir Liste üzerindeki `collect` metodunu iş başında görelim. Orijinal örnek ID listesine Listenin `collect` metodunu uygulamak için `collect.nf` dosyasını değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de gruplandı)" }

    // List.collect() - her elemanı dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()} öğeye dönüştürüldü)"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de gruplandı)" }
    ```

Bu yeni snippet'te:

- Orijinal listedeki her örnek ID'yi dönüştürmek için Listenin `collect` metodunu kullanan yeni bir `formatted_ids` değişkeni tanımlıyoruz
- Sonucu `println` kullanarak yazdırıyoruz

Değiştirilmiş iş akışını çalıştırın:

```bash
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() sonucu: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 öğe 3 öğeye dönüştürüldü)
    Bireysel kanal öğesi: sample_001
    Bireysel kanal öğesi: sample_002
    Bireysel kanal öğesi: sample_003
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe 1'de gruplandı)
    ```

Bu sefer, bazı örnekleri hariç tutan bir filtre seçtiğimiz için verinin yapısını DEĞİŞTİRMEDİK, listede hala 3 öğemiz var, ancak değiştirilmiş değerlerle yeni bir liste üretmek için Listenin `collect` metodunu kullanarak her öğeyi DÖNÜŞTÜRDÜK. Bu, bir kanal üzerinde `map` operatörü kullanmaya benzer, ancak bir kanal yerine bir Liste veri yapısı üzerinde çalışıyor.

`collect` burada bir noktayı vurgulamak için kullandığımız aşırı bir durumdur. Temel ders şudur: iş akışları yazarken her zaman **veri yapıları** (Listeler, Map'ler, vb.) ile **kanallar** (veri akışı yapıları) arasında ayrım yapın. İşlemler aynı isimleri paylaşabilir ancak çağrıldıkları türe bağlı olarak tamamen farklı davranabilir.

### 1.3. Yayılma Operatörü (`*.`) - Özellik Çıkarma için Kısayol

Listenin `collect` metoduyla ilgili olan yayılma operatörüdür (`*.`), koleksiyonlardan özellikleri çıkarmak için kısa bir yol sağlar. Esasen yaygın bir `collect` kalıbı için sözdizimsel şekerdir.

`collect.nf` dosyamıza bir gösteri ekleyelim:

=== "Sonra"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal emisyonunu bir araya toplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de gruplandı)" }

    // List.collect() - her elemanı dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()} öğeye dönüştürüldü)"

    // Yayılma operatörü - kısa özellik erişimi
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
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe 1'de gruplandı)" }

    // List.collect() - her elemanı dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()} öğeye dönüştürüldü)"
    ```

Güncellenmiş iş akışını çalıştırın:

```bash title="Yayılma operatörünü test et"
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() sonucu: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 öğe 3 öğeye dönüştürüldü)
    Yayılma operatörü sonucu: [s1, s2, s3]
    Bireysel kanal öğesi: sample_001
    Bireysel kanal öğesi: sample_002
    Bireysel kanal öğesi: sample_003
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe 1'de gruplandı)
    ```

Yayılma operatörü `*.`, yaygın bir collect kalıbı için kısayoldur:

```groovy
// Bunlar eşdeğerdir:
def ids = samples*.id
def ids = samples.collect { it.id }

// Metod çağrılarıyla da çalışır:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Yayılma operatörü, bir nesne listesinden tek bir özelliği çıkarmanız gerektiğinde özellikle yararlıdır - tam `collect` closure'ını yazmaktan daha okunabilirdir.

!!! tip "Yayılma vs Collect Ne Zaman Kullanılır"

    - **Yayılma (`*.`) kullanın** basit özellik erişimi için: `samples*.id`, `files*.name`
    - **Collect kullanın** dönüşümler veya karmaşık mantık için: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Özet

Bu bölümde şunları öğrendiniz:

- **Veri akışı vs betik yazma**: Kanal operatörleri verilerin pipeline'ınız boyunca nasıl aktığını düzenlerken, betik yazma bireysel veri öğelerini dönüştürür
- **Türleri anlama**: Aynı metod adı (örneğin `collect`) çağrıldığı türe (Kanal vs Liste) bağlı olarak farklı davranabilir
- **Bağlam önemlidir**: Kanallarla (veri akışı) mı yoksa veri yapılarıyla (betik yazma) mı çalıştığınızın her zaman farkında olun

Bu sınırları anlamak, hata ayıklama, dokümantasyon ve sürdürülebilir iş akışları yazmak için esastır.

Sırada, gerçek dünya verileriyle başa çıkmak için gerekli olan string işleme yeteneklerine daha derinlemesine dalacağız.

---

## 2. String İşleme ve Dinamik Betik Oluşturma

String işlemede ustalaşmak, kırılgan iş akışlarını sağlam pipeline'lardan ayırır. Bu bölüm, karmaşık dosya adlarını ayrıştırmayı, dinamik betik oluşturmayı ve değişken interpolasyonunu kapsar.

### 2.1. Kalıp Eşleştirme ve Düzenli İfadeler

Biyoinformatik dosyaları genellikle meta verileri kodlayan karmaşık adlandırma kurallarına sahiptir. Bunu düzenli ifadelerle kalıp eşleştirme kullanarak otomatik olarak çıkaralım.

`main.nf` iş akışımıza geri döneceğiz ve dosya adlarından ek örnek bilgilerini çıkarmak için bazı kalıp eşleştirme mantığı ekleyeceğiz. Veri setimizdeki FASTQ dosyaları, `SAMPLE_001_S1_L001_R1_001.fastq.gz` gibi adlara sahip Illumina tarzı adlandırma kurallarını takip eder. Bunlar şifreli görünebilir, ancak aslında örnek ID, şerit numarası ve okuma yönü gibi yararlı meta verileri kodlarlar. Bu adları ayrıştırmak için regex yeteneklerini kullanacağız.

Mevcut `main.nf` iş akışınızda aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Veri dönüşümü için betik yazma
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
                // Veri dönüşümü için betik yazma
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

1. `~/kalıp/` sözdizimini kullanan **düzenli ifade literalleri** - bu, ters eğik çizgileri kaçırmaya gerek kalmadan bir regex kalıbı oluşturur
2. `=~` operatörü ile **kalıp eşleştirme** - bu, bir string'i bir regex kalıbıyla eşleştirmeye çalışır
3. `[0][1]`, `[0][2]`, vb. ile grupları yakalayan **eşleştirici nesneler** - `[0]` tüm eşleşmeyi ifade eder, `[1]`, `[2]`, vb. parantez içindeki yakalanan grupları ifade eder

Regex kalıbını `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` parçalayalım:

| Kalıp               | Eşleşir                                    | Yakalar                                |
| ------------------- | ------------------------------------------ | -------------------------------------- |
| `^(.+)`             | Baştan örnek adı                           | Grup 1: örnek adı                      |
| `_S(\d+)`           | Örnek numarası `_S1`, `_S2`, vb.           | Grup 2: örnek numarası                 |
| `_L(\d{3})`         | Şerit numarası `_L001`                     | Grup 3: şerit (3 basamak)              |
| `_(R[12])`          | Okuma yönü `_R1` veya `_R2`                | Grup 4: okuma yönü                     |
| `_(\d{3})`          | Parça numarası `_001`                      | Grup 5: parça (3 basamak)              |
| `\.fastq(?:\.gz)?$` | Dosya uzantısı `.fastq` veya `.fastq.gz`   | Yakalanmadı (?: yakalamayan)           |

Bu, meta verileri otomatik olarak çıkarmak için Illumina tarzı adlandırma kurallarını ayrıştırır.

Değiştirilmiş iş akışını çalıştırın:

```bash title="Kalıp eşleştirmeyi test et"
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

### 2.2. Süreçlerde Dinamik Betik Oluşturma

Süreç betik blokları esasen shell'e geçirilen çok satırlı string'lerdir. Girdi özelliklerine göre farklı betik string'leri dinamik olarak oluşturmak için **koşullu mantık** (if/else, üçlü operatörler) kullanabilirsiniz. Bu, tek uçlu ve çift uçlu dizileme okumaları gibi çeşitli girdi türlerini ele almak için süreç tanımlarını çoğaltmadan gereklidir.

Bu kalıbı gösteren iş akışımıza bir süreç ekleyelim. `modules/fastp.nf` dosyasını açın ve bir göz atın:

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

Süreç, girdi olarak FASTQ dosyalarını alır ve adaptörleri kesmek ve düşük kaliteli okumaları filtrelemek için `fastp` aracını çalıştırır. Ne yazık ki, bu süreci yazan kişi örnek veri setimizdeki tek uçlu okumalara izin vermedi. Bunu iş akışımıza ekleyelim ve ne olduğunu görelim:

İlk olarak, `main.nf` iş akışınızın ilk satırına modülü dahil edin:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Ardından `ch_samples` kanalını `FASTP` sürecine bağlamak için `workflow` bloğunu değiştirin:

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

Sürecin ikinci girdi dosyası için `null` değeriyle `fastp` çalıştırmaya çalıştığını ve bunun başarısız olmasına neden olduğunu görebilirsiniz. Bunun nedeni, veri setimizin tek uçlu okumalar içermesi, ancak sürecin çift uçlu okumaları (aynı anda iki girdi dosyası) bekleyecek şekilde sabit kodlanmış olmasıdır.

Bunu, `FASTP` süreci `script:` bloğuna koşullu mantık ekleyerek düzeltin. Bir if/else ifadesi okuma dosyası sayısını kontrol eder ve komutu buna göre ayarlar.

=== "Sonra"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Basit tek uçlu vs çift uçlu algılama
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

Artık iş akışı hem tek uçlu hem de çift uçlu okumaları zarif bir şekilde ele alabilir. Koşullu mantık, girdi dosyalarının sayısını kontrol eder ve `fastp` için uygun komutu oluşturur. Çalışıp çalışmadığını görelim:

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

İyi görünüyor! Gerçekten çalıştırılan komutları kontrol edersek (görev hash'iniz için özelleştirin):

```console title="Çalıştırılan komutları kontrol et"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Nextflow'un tek uçlu okumalar için doğru komutu doğru bir şekilde seçtiğini görebiliriz:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Dinamik betik mantığının başka bir yaygın kullanımı [Nextflow for Science Genomics modülünde](../../nf4science/genomics/02_joint_calling) görülebilir. Bu modülde, çağrılan GATK süreci birden fazla girdi dosyası alabilir, ancak her birinin doğru bir komut satırı oluşturmak için `-V` ile öneklenmesi gerekir. Süreç, bir girdi dosyaları koleksiyonunu (`all_gvcfs`) doğru komut argümanlarına dönüştürmek için betik yazma kullanır:

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

Süreç betik bloklarında betik yazma kullanmanın bu kalıpları son derece güçlüdür ve birçok senaryoda uygulanabilir - değişken girdi türlerini ele almaktan dosya koleksiyonlarından karmaşık komut satırı argümanları oluşturmaya kadar, süreçlerinizi gerçek dünya verilerinin çeşitli gereksinimlerine gerçekten uyarlanabilir hale getirir.

### 2.3. Değişken İnterpolasyonu: Nextflow ve Shell Değişkenleri

Süreç betikleri Nextflow değişkenlerini, shell değişkenlerini ve komut ikamelerini karıştırır, her biri farklı interpolasyon sözdizimi ile. Yanlış sözdizimi kullanmak hatalara neden olur. Bunları bir işleme raporu oluşturan bir süreçle keşfedelim.

Modül dosyası `modules/generate_report.nf` dosyasına bir göz atın:

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

Bu süreç, örnek ID ve dosya adı ile basit bir rapor yazar. Şimdi farklı türde değişkenleri karıştırmamız gerektiğinde ne olduğunu görmek için çalıştıralım.

Süreci `main.nf` dosyanıza dahil edin ve iş akışına ekleyin:

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

Şimdi iş akışını çalıştırın ve `results/reports/` içinde oluşturulan raporları kontrol edin. Her örnek hakkında temel bilgiler içermelidir.

<!-- TODO: add the run command -->

??? success "Komut çıktısı"

    ```console
    <!-- TODO: output -->
    ```

Peki işlemenin ne zaman ve nerede gerçekleştiğine dair bilgi eklemek istersek? Rapora mevcut kullanıcıyı, ana bilgisayar adını ve tarihi dahil etmek için **shell** değişkenlerini ve biraz komut ikamesini kullanmak üzere süreci değiştirelim:

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

Bunu çalıştırırsanız, bir hata göreceksiniz - Nextflow `${USER}` değişkenini var olmayan bir Nextflow değişkeni olarak yorumlamaya çalışır.

??? failure "Komut çıktısı"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Bash'in bunu ele alabilmesi için kaçırmamız gerekiyor.

Shell değişkenlerini ve komut ikamelerini ters eğik çizgi (`\`) ile kaçırarak bunu düzeltin:

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

### Özet

Bu bölümde **string işleme** tekniklerini öğrendiniz:

- **Dosya ayrıştırma için düzenli ifadeler**: Karmaşık dosya adlandırma kurallarından meta verileri çıkarmak için `=~` operatörünü ve regex kalıplarını (`~/kalıp/`) kullanma
- **Dinamik betik oluşturma**: Girdi özelliklerine göre farklı betik string'leri oluşturmak için koşullu mantık (if/else, üçlü operatörler) kullanma
- **Değişken interpolasyonu**: Nextflow'un string'leri ne zaman yorumladığını ve shell'in ne zaman yorumladığını anlama
  - `${var}` - Nextflow değişkenleri (iş akışı derleme zamanında Nextflow tarafından interpolasyona uğrar)
  - `\${var}` - Shell ortam değişkenleri (kaçırılmış, çalışma zamanında bash'e geçirilir)
  - `\$(cmd)` - Shell komut ikamesi (kaçırılmış, çalışma zamanında bash tarafından çalıştırılır)

Bu string işleme ve oluşturma kalıpları, gerçek dünya biyoinformatik iş akışlarında karşılaşacağınız çeşitli dosya formatlarını ve adlandırma kurallarını ele almak için esastır.

---

## 3. Yeniden Kullanılabilir Fonksiyonlar Oluşturma

Kanal operatörlerinde veya süreç tanımlarında satır içi karmaşık iş akışı mantığı, okunabilirliği ve sürdürülebilirliği azaltır. **Fonksiyonlar**, bu mantığı adlandırılmış, yeniden kullanılabilir bileşenlere çıkarmanıza olanak tanır.

Map işlemimiz uzun ve karmaşık hale geldi. Bunu `def` anahtar kelimesini kullanarak yeniden kullanılabilir bir fonksiyona çıkaralım.

Mevcut iş akışımızla bunun nasıl göründüğünü göstermek için, `separateMetadata` adlı yeniden kullanılabilir bir fonksiyon tanımlamak için `def` kullanarak aşağıdaki değişikliği yapın:

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

Bu mantığı bir fonksiyona çıkararak, gerçek iş akışı mantığını çok daha temiz bir şeye indirgedik:

```groovy title="minimal iş akışı"
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

Çıktı, her iki sürecin de başarıyla tamamlandığını göstermelidir. İş akışı artık çok daha temiz ve sürdürülebilir, tüm karmaşık meta veri işleme mantığı `separateMetadata` fonksiyonunda kapsüllenmiş.

### Özet

Bu bölümde **fonksiyon oluşturmayı** öğrendiniz:

- **`def` ile fonksiyon tanımlama**: Adlandırılmış fonksiyonlar oluşturmak için anahtar kelime (Python'daki `def` veya JavaScript'teki `function` gibi)
- **Fonksiyon kapsamı**: Betik düzeyinde tanımlanan fonksiyonlar Nextflow iş akışınız boyunca erişilebilir
- **Dönüş değerleri**: Fonksiyonlar otomatik olarak son ifadeyi döndürür veya açık `return` kullanır
- **Daha temiz kod**: Karmaşık mantığı fonksiyonlara çıkarmak, herhangi bir dilde temel bir yazılım mühendisliği uygulamasıdır

Sırada, dinamik kaynak tahsisi için süreç yönergelerinde closure'ların nasıl kullanılacağını keşfedeceğiz.

---

## 4. Closure'larla Dinamik Kaynak Yönergeleri

Şimdiye kadar süreçlerin `script` bloğunda betik yazma kullandık. Ancak **closure'lar** (Bölüm 1.1'de tanıtıldı) süreç yönergelerinde de, özellikle dinamik kaynak tahsisi için inanılmaz derecede yararlıdır. FASTP sürecimize örnek özelliklerine göre uyum sağlayan kaynak yönergeleri ekleyelim.

### 4.1. Örneğe özgü kaynak tahsisi

Şu anda, FASTP sürecimiz varsayılan kaynakları kullanıyor. Yüksek derinlikli örnekler için daha fazla CPU tahsis ederek daha akıllı hale getirelim. Dinamik bir `cpus` yönergesi ve statik bir `memory` yönergesi içerecek şekilde `modules/fastp.nf` dosyasını düzenleyin:

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

`{ meta.depth > 40000000 ? 2 : 1 }` closure'ı **üçlü operatörü** (Bölüm 1.1'de ele alındı) kullanır ve her görev için değerlendirilir, örnek başına kaynak tahsisine izin verir. Yüksek derinlikli örnekler (>40M okuma) 2 CPU alırken, diğerleri 1 CPU alır.

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

Herhangi bir görev için CPU tahsisini görmek için çalıştırılan tam `docker` komutunu kontrol edebilirsiniz:

```console title="Docker komutunu kontrol et"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Şöyle bir şey görmelisiniz:

```bash title="docker komutu"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Bu örnekte, yüksek derinlikli bir örnek olduğu için 2 CPU (`--cpu-shares 2048`) isteyen bir örnek seçtik, ancak örnek derinliğine bağlı olarak farklı CPU tahsisleri görmelisiniz. Bunu diğer görevler için de deneyin.

### 4.2. Yeniden deneme stratejileri

Başka bir güçlü kalıp, yeniden deneme stratejileri için `task.attempt` kullanmaktır. Bunun neden yararlı olduğunu göstermek için, FASTP'ye bellek tahsisini ihtiyaç duyduğundan daha aza indirerek başlayacağız. `modules/fastp.nf` dosyasındaki `memory` yönergesini `1.GB` olarak değiştirin:

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

Bu, sürecin bellek sınırlarını aştığı için öldürüldüğünü gösterir.

Bu, gerçek dünya iş akışlarında çok yaygın bir senaryodur - bazen bir görevin ne kadar belleğe ihtiyaç duyacağını çalıştırana kadar bilemezsiniz.

İş akışımızı daha sağlam hale getirmek için, her denemede bellek tahsisini artıran bir yeniden deneme stratejisi uygulayabiliriz, bir kez daha bir Groovy closure kullanarak. `memory` yönergesini temel belleği `task.attempt` ile çarpmak için değiştirin ve `errorStrategy 'retry'` ve `maxRetries 2` yönergelerini ekleyin:

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

Şimdi süreç yetersiz bellek nedeniyle başarısız olursa, Nextflow daha fazla bellekle yeniden deneyecektir:

- İlk deneme: 1 GB (task.attempt = 1)
- İkinci deneme: 2.GB (task.attempt = 2)

... ve böyle devam eder, `maxRetries` sınırına kadar.

### Özet

Closure'larla dinamik yönergeler şunları yapmanıza olanak tanır:

- Girdi özelliklerine göre kaynakları tahsis etme
- Artan kaynaklarla otomatik yeniden deneme stratejileri uygulama
- Birden fazla faktörü (meta veriler, deneme numarası, öncelikler) birleştirme
- Karmaşık kaynak hesaplamaları için koşullu mantık kullanma

Bu, iş akışlarınızı hem daha verimli (aşırı tahsis etmeme) hem de daha sağlam (daha fazla kaynakla otomatik yeniden deneme) hale getirir.

---

## 5. Koşullu Mantık ve Süreç Kontrolü

Daha önce, kanal verilerini dönüştürmek için betik yazma ile `.map()` kullandık. Şimdi hangi süreçlerin verilere göre çalıştırılacağını kontrol etmek için koşullu mantık kullanacağız - farklı örnek türlerine uyum sağlayan esnek iş akışları için gerekli.

Nextflow'un [veri akışı operatörleri](https://www.nextflow.io/docs/latest/reference/operator.html), çalışma zamanında değerlendirilen closure'lar alır, koşullu mantığın kanal içeriğine göre iş akışı kararlarını yönlendirmesini sağlar.

### 5.1. `.branch()` ile Yönlendirme

Örneğin, dizileme örneklerimizin yalnızca belirli bir eşiğin üzerinde kapsama sahip insan örnekleri olmaları durumunda FASTP ile kesilmesi gerektiğini varsayalım. Fare örnekleri veya düşük kapsama örnekleri bunun yerine Trimgalore ile çalıştırılmalıdır (bu uydurulmuş bir örnektir, ancak noktayı gösterir).

`modules/trimgalore.nf` dosyasında basit bir Trimgalore süreci sağladık, isterseniz bir göz atın, ancak ayrıntılar bu alıştırma için önemli değil. Temel nokta, örnekleri meta verilerine göre yönlendirmek istememizdir.

`modules/trimgalore.nf` dosyasından yeni formu dahil edin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... ve ardından örnekleri meta verilerine göre dallandırmak ve uygun kesme sürecinden geçirmek için `main.nf` iş akışınızı şu şekilde değiştirin:

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

İş akışı yürütmesini kontrol etmek için başka bir güçlü kalıp, hangi öğelerin pipeline'dan devam etmesi gerektiğini belirlemek için bir closure kullanan `.filter()` operatörüdür. Filtre closure'ı içinde, hangi öğelerin geçeceğine karar veren **boolean ifadeleri** yazacaksınız.

Nextflow'un (birçok dinamik dil gibi) boolean bağlamlarda hangi değerlerin `true` veya `false` olarak değerlendirileceğini belirleyen bir **"doğruluk"** kavramı vardır:

- **Doğru**: Null olmayan değerler, boş olmayan string'ler, sıfır olmayan sayılar, boş olmayan koleksiyonlar
- **Yanlış**: `null`, boş string'ler `""`, sıfır `0`, boş koleksiyonlar `[]` veya `[:]`, `false`

Bu, `meta.id` tek başına (açık `!= null` olmadan) ID'nin var olup olmadığını ve boş olmadığını kontrol eder anlamına gelir. Kalite gereksinimlerimizi karşılamayan örnekleri filtrelemek için bunu kullanalım.

Dal işleminden önce aşağıdakileri ekleyin:

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

Bazı örnekleri hariç tutan bir filtre seçtiğimiz için, daha az görev yürütüldü.

Filtre ifadesi `meta.id && meta.organism && meta.depth >= 25000000`, doğruluğu açık karşılaştırmalarla birleştirir:

- `meta.id && meta.organism` her iki alanın da var olduğunu ve boş olmadığını kontrol eder (doğruluk kullanarak)
- `meta.depth >= 25000000` açık bir karşılaştırma ile yeterli dizileme derinliğini sağlar

!!! note "Pratikte Doğruluk"

    `meta.id && meta.organism` ifadesi şunu yazmaktan daha kısadır:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Bu, filtreleme mantığını çok daha temiz ve okunması daha kolay hale getirir.

### Özet

Bu bölümde, Nextflow operatörlerinin closure arayüzlerini (`.branch{}` ve `.filter{}` gibi) kullanarak iş akışı yürütmesini kontrol etmek için koşullu mantık kullanmayı, kısa koşullu ifadeler yazmak için doğruluğu kullanmayı öğrendiniz.

Pipeline'ımız artık örnekleri uygun süreçlerden akıllıca geçiriyor, ancak üretim iş akışlarının geçersiz verileri zarif bir şekilde ele alması gerekir. İş akışımızı eksik veya null değerlere karşı sağlam hale getirelim.

---

## 6. Güvenli Navigasyon ve Elvis Operatörleri

`separateMetadata` fonksiyonumuz şu anda tüm CSV alanlarının mevcut ve geçerli olduğunu varsayıyor. Peki eksik verilerle ne olur? Öğrenelim.

### 6.1. Sorun: Var Olmayan Özelliklere Erişim

Diyelim ki isteğe bağlı dizileme çalıştırma bilgisi için destek eklemek istiyoruz. Bazı laboratuvarlarda, örneklerin dizileme çalıştırma ID'si veya parti numarası için ek bir alanı olabilir, ancak mevcut CSV'mizde bu sütun yok. Yine de ona erişmeyi deneyelim.

Bir run_id alanı içerecek şekilde `separateMetadata` fonksiyonunu değiştirin:

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

Bu, NullPointerException ile çöküyor.

Sorun, `row.run_id`'nin `null` döndürmesidir çünkü `run_id` sütunu CSV'mizde mevcut değildir. `null` üzerinde `.toUpperCase()` çağırmaya çalıştığımızda çöküyor. Güvenli navigasyon operatörünün günü kurtardığı yer burasıdır.

### 6.2. Güvenli Navigasyon Operatörü (`?.`)

Güvenli navigasyon operatörü (`?.`), `null` bir değer üzerinde çağrıldığında istisna atmak yerine `null` döndürür. `?.` öncesindeki nesne `null` ise, tüm ifade metodu çalıştırmadan `null` olarak değerlendirilir.

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

Çökme yok! İş akışı artık eksik alanı zarif bir şekilde ele alıyor. `row.run_id` `null` olduğunda, `?.` operatörü `.toUpperCase()` çağrısını önler ve `run_id`, istisna oluşturmak yerine `null` olur.

### 6.3. Varsayılanlar için Elvis Operatörü (`?:`)

Elvis operatörü (`?:`), sol taraf "yanlış" olduğunda (daha önce açıklandığı gibi) varsayılan değerler sağlar. Yan taraftan bakıldığında `?:` Elvis Presley'in ünlü saçı ve gözleri gibi göründüğü için Elvis operatörü olarak adlandırılır!

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

Mükemmel! Artık tüm örneklerin gerçek çalıştırma ID'leri (büyük harfle) veya varsayılan 'UNSPECIFIED' değeri ile bir `run` alanı var. `?.` ve `?:` kombinasyonu hem güvenlik (çökme yok) hem de mantıklı varsayılanlar sağlar.

Çalıştığını doğruladığımıza göre `.view()` operatörünü çıkarın.

!!! tip "Güvenli Navigasyon ve Elvis'i Birleştirme"

    `value?.method() ?: 'default'` kalıbı üretim iş akışlarında yaygındır:

    - `value?.method()` - Metodu güvenli bir şekilde çağırır, `value` `null` ise `null` döndürür
    - `?: 'default'` - Sonuç `null` ise yedek sağlar

    Bu kalıp eksik/eksik verileri zarif bir şekilde ele alır.

Bu operatörleri fonksiyonlarda, operatör closure'larında (`.map{}`, `.filter{}`), süreç betiklerinde ve yapılandırma dosyalarında tutarlı bir şekilde kullanın. Gerçek dünya verileriyle başa çıkarken çökmeleri önlerler.

### Özet

- **Güvenli navigasyon (`?.`)**: Null değerlerde çökmeleri önler - istisna atmak yerine null döndürür
- **Elvis operatörü (`?:`)**: Varsayılanlar sağlar - `value ?: 'default'`
- **Birleştirme**: `value?.method() ?: 'default'` yaygın kalıptır

Bu operatörler iş akışlarını eksik verilere karşı dirençli hale getirir - gerçek dünya çalışması için gereklidir.

---

## 7. `error()` ve `log.warn` ile Doğrulama

Bazen girdi parametreleri geçersizse iş akışını hemen durdurmanız gerekir. Nextflow'da, doğrulama mantığı uygulamak için `error()` ve `log.warn` gibi yerleşik fonksiyonların yanı sıra `if` ifadeleri ve boolean mantığı gibi standart programlama yapılarını kullanabilirsiniz. İş akışımıza doğrulama ekleyelim.

İş akışı bloğunuzdan önce bir doğrulama fonksiyonu oluşturun, iş akışından çağırın ve kanal oluşturmayı CSV dosya yolu için bir parametre kullanacak şekilde değiştirin. Parametre eksikse veya dosya yoksa, açık bir mesajla yürütmeyi durdurmak için `error()` çağırın.

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

İş akışı, daha sonra gizemli bir şekilde başarısız olmak yerine açık bir hata mesajıyla hemen durur

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

`separateMetadata` fonksiyonu içine de doğrulama ekleyebilirsiniz. Düşük dizileme derinliğine sahip örnekler için uyarılar vermek, ancak iş akışının devam etmesine izin vermek için ölümcül olmayan `log.warn` kullanalım:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Verilerin mantıklı olduğunu doğrula
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

### Özet

- **`error()`**: Açık mesajla iş akışını hemen durdurur
- **`log.warn`**: İş akışını durdurmadan uyarılar verir
- **Erken doğrulama**: İşlemeden önce girdileri kontrol ederek yararlı hatalarla hızlı başarısızlık sağlar
- **Doğrulama fonksiyonları**: İş akışı başlangıcında çağrılabilen yeniden kullanılabilir doğrulama mantığı oluşturur

Uygun doğrulama, sorunları erken yakalayarak ve açık hata mesajlarıyla iş akışlarını daha sağlam ve kullanıcı dostu hale getirir.

---

## 8. İş Akışı Olay İşleyicileri

Şimdiye kadar, iş akışı betiklerimizde ve süreç tanımlarımızda kod yazıyorduk. Ancak bilmeniz gereken bir önemli özellik daha var: iş akışı olay işleyicileri.

Olay işleyicileri, iş akışınızın yaşam döngüsünde belirli noktalarda çalışan closure'lardır. Loglama, bildirimler veya temizleme işlemleri eklemek için mükemmeldirler. Bu işleyiciler, iş akışı tanımınızın yanında iş akışı betiğinizde tanımlanmalıdır.

### 8.1. `onComplete` İşleyicisi

En yaygın kullanılan olay işleyicisi, iş akışınız bittiğinde (başarılı veya başarısız olsun) çalışan `onComplete`'tir. Pipeline sonuçlarımızı özetlemek için bir tane ekleyelim.

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

Bu closure, iş akışı tamamlandığında çalışır. İçinde, yürütme hakkında yararlı özellikler sağlayan `workflow` nesnesine erişiminiz vardır.

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
            println "workDir   : ${workflow.workDir}"```groovy title="main.nf" linenums="66" hl_lines="5-16"
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
    WARN: sample_002 için düşük dizileme derinliği: 25000000
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

`onComplete`'in yanı sıra, kullanabileceğiniz bir başka olay işleyicisi daha var: `onError`, yalnızca iş akışı başarısız olursa çalışır:

```groovy title="main.nf - onError işleyicisi"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "="* 50
        println "Pipeline yürütmesi başarısız oldu!"
        println "Hata mesajı: ${workflow.errorMessage}"
        println "="* 50

        // Ayrıntılı hata logu yaz
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        İş Akışı Hata Raporu
        =====================
        Zaman: ${new Date()}
        Hata: ${workflow.errorMessage}
        Hata raporu: ${workflow.errorReport ?: 'Ayrıntılı rapor mevcut değil'}
        """

        println "Hata ayrıntıları şuraya yazıldı: ${error_file}"
    }
}
```

İş akışı betiğinizde birden fazla işleyiciyi birlikte kullanabilirsiniz:

```groovy title="main.nf - Birleştirilmiş işleyiciler"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "İş akışı başarısız oldu: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "BAŞARILI ✅" : "BAŞARISIZ ❌"

        println """
        Pipeline tamamlandı: ${status}
        Süre: ${duration_mins} dakika
        """
    }
}
```

### Özet

Bu bölümde şunları öğrendiniz:

- **Olay işleyici closure'ları**: Farklı yaşam döngüsü noktalarında çalışan iş akışı betiğinizdeki closure'lar
- **`onComplete` işleyicisi**: Yürütme özetleri ve sonuç raporlaması için
- **`onError` işleyicisi**: Hata işleme ve başarısızlıkları loglama için
- **Workflow nesne özellikleri**: `workflow.success`, `workflow.duration`, `workflow.errorMessage`, vb.'ye erişim

Olay işleyicileri, sofistike loglama ve bildirim yetenekleri eklemek için iş akışı betiklerinizde Nextflow dilinin tam gücünü nasıl kullanabileceğinizi gösterir.

---

## Özet

Tebrikler, başardınız!

Bu yan görev boyunca, temel meta veri işlemeden sofistike, üretime hazır bir iş akışına evrilen kapsamlı bir örnek işleme pipeline'ı oluşturdunuz.
Her bölüm bir öncekinin üzerine inşa edildi ve programlama yapılarının basit iş akışlarını güçlü veri işleme sistemlerine nasıl dönüştürdüğünü gösterdi, aşağıdaki faydalarla:

- **Daha net kod**: Veri akışı vs betik yazmayı anlamak daha organize iş akışları yazmanıza yardımcı olur
- **Sağlam işleme**: Güvenli navigasyon ve Elvis operatörleri iş akışlarını eksik verilere karşı dirençli hale getirir
- **Esnek işleme**: Koşullu mantık, iş akışlarınızın farklı örnek türlerini uygun şekilde işlemesine olanak tanır
- **Uyarlanabilir kaynaklar**: Dinamik yönergeler, girdi özelliklerine göre kaynak kullanımını optimize eder

Bu ilerleme, birkaç örneği işleyen araştırma prototiplerinden laboratuvarlar ve kurumlar arasında binlerce örneği işleyen üretim sistemlerine kadar gerçek dünya biyoinformatik pipeline'larının evrimini yansıtır.
Çözdüğünüz her zorluk ve öğrendiğiniz her kalıp, geliştiricilerin Nextflow iş akışlarını ölçeklendirirken karşılaştıkları gerçek sorunları yansıtır.

Bu kalıpları kendi çalışmanızda uygulamak, sağlam, üretime hazır iş akışları oluşturmanızı sağlayacaktır.

### Temel kalıplar

1.  **Veri Akışı vs Betik Yazma:** Veri akışı işlemleri (kanal düzenlemesi) ile betik yazma (verileri manipüle eden kod) arasında ayrım yapmayı öğrendiniz, Kanal vs Liste üzerindeki `collect` gibi farklı türlerdeki işlemler arasındaki önemli farklar dahil.

    - Veri akışı: kanal düzenlemesi

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Betik yazma: koleksiyonlarda veri işleme

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Gelişmiş String İşleme**: Dosya adlarını ayrıştırmak için düzenli ifadelerde, süreçlerde dinamik betik oluşturmada ve değişken interpolasyonunda (Nextflow vs Bash vs Shell) ustalaştınız.

    - Kalıp eşleştirme

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

    - Dosya koleksiyonunu komut argümanlarına (süreç betik bloğunda)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Yeniden Kullanılabilir Fonksiyonlar Oluşturma**: Karmaşık mantığı kanal operatörlerinden çağrılabilen adlandırılmış fonksiyonlara çıkarmayı öğrendiniz, iş akışlarını daha okunabilir ve sürdürülebilir hale getirdiniz.

    - Adlandırılmış bir fonksiyon tanımlayın

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

    - Adlandırılmış fonksiyonu bir iş akışında çağırın

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closure'larla Dinamik Kaynak Yönergeleri**: Girdi özelliklerine göre uyarlanabilir kaynak tahsisi için süreç yönergelerinde closure'ları kullanmayı keşfettiniz.

    - Adlandırılmış closure'lar ve kompozisyon

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Kapsam erişimli closure'lar

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Koşullu Mantık ve Süreç Kontrolü**: `.branch()` ve `.filter()` operatörlerini kullanarak akıllı yönlendirme eklediniz, kısa koşullu ifadeler için doğruluğu kullandınız.

    - Verileri farklı iş akışı dallarından yönlendirmek için `.branch()` kullanın

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Groovy Truth ile Boolean değerlendirme

    ```groovy
    if (sample.files) println "Dosyalar var"
    ```

    - 'Doğruluk' ile verileri alt kümelere ayırmak için `filter()` kullanın

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

8.  **Yapılandırma Olay İşleyicileri**: Loglama, bildirimler ve yaşam döngüsü yönetimi için iş akışı olay işleyicilerini (`onComplete` ve `onError`) kullanmayı öğrendiniz.

    - Loglamak ve bildirimde bulunmak için `onComplete` kullanma

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
        // Ayrıntılı hata logu yaz
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Zaman: ${new Date()}
        Hata: ${workflow.errorMessage}
        Hata raporu: ${workflow.errorReport ?: 'Ayrıntılı rapor mevcut değil'}
        """

        println "Hata ayrıntıları şuraya yazıldı: ${error_file}"
    }
    ```

### Ek kaynaklar

- [Nextflow Dil Referansı](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operatörleri](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Betik Sözdizimi](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standart Kütüphanesi](https://nextflow.io/docs/latest/reference/stdlib.html)

Daha gelişmiş özellikleri keşfetmeniz gerektiğinde bu kaynaklara mutlaka göz atın.

Becerilerinizi geliştirmek için pratik yapmaktan ve genişletmekten faydalanacaksınız:

- Veri akışı ile betik yazma arasında uygun ayrım yaparak daha temiz iş akışları yazın
- Nextflow, Bash ve shell değişkenleriyle yaygın tuzaklardan kaçınmak için değişken interpolasyonunda ustalaşın
- Verimli, uyarlanabilir iş akışları için dinamik kaynak yönergeleri kullanın
- Dosya koleksiyonlarını düzgün biçimlendirilmiş komut satırı argümanlarına dönüştürün
- Regex ve string işleme kullanarak farklı dosya adlandırma kurallarını ve girdi formatlarını zarif bir şekilde ele alın
- Gelişmiş closure kalıpları ve fonksiyonel programlama kullanarak yeniden kullanılabilir, sürdürülebilir kod oluşturun
- Koleksiyon işlemlerini kullanarak karmaşık veri setlerini işleyin ve düzenleyin
- İş akışlarınızı üretime hazır hale getirmek için doğrulama, hata işleme ve loglama ekleyin
- Olay işleyicileriyle iş akışı yaşam döngüsü yönetimini uygulayın

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ altındaki düğmeye tıklayın.
