# Temel Nextflow Betik Kalıpları

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, Java Sanal Makinesi üzerinde çalışan bir programlama dilidir. Nextflow, [Groovy](http://groovy-lang.org/) üzerine inşa edilmiş olup sözdiziminin büyük bölümünü Groovy ile paylaşsa da yalnızca "uzantılı Groovy" değildir; tam olarak belirlenmiş bir [sözdizimine](https://nextflow.io/docs/latest/reference/syntax.html) ve [standart kütüphaneye](https://nextflow.io/docs/latest/reference/stdlib.html) sahip bağımsız bir dildir.

Değişkenler, map'ler ve listeler için temel sözdiziminin ötesine geçmeden Nextflow ile çok şey yazabilirsiniz. Nextflow eğitimlerinin büyük çoğunluğu iş akışı orkestrasyonuna (kanallar, süreçler ve veri akışı) odaklanır; yalnızca bunlarla da şaşırtıcı ölçüde ilerleyebilirsiniz.

Ancak veri işlemeniz, karmaşık dosya adlarını ayrıştırmanız, koşullu mantık uygulamanız veya dayanıklı üretim iş akışları oluşturmanız gerektiğinde, kodunuzun iki farklı boyutunu düşünmek işe yarar: **veri akışı** (kanallar, operatörler, süreçler ve iş akışları) ve **betik yazımı** (closure'lar, fonksiyonlar ve süreç betikleri içindeki kod). Bu ayrım bir ölçüde keyfi olsa da —sonuçta hepsi Nextflow kodudur— pipeline'ınızı ne zaman orkestre ettiğinizi, ne zaman veri işlediğinizi anlamak için kullanışlı bir zihinsel model sunar. Her ikisinde de ustalaşmak, açık ve sürdürülebilir iş akışları yazma yetkinliğinizi önemli ölçüde artırır.

### Öğrenme hedefleri

Bu yan görev, temel kavramlardan üretime hazır kalıplara uzanan uygulamalı bir yolculuğa çıkarır.
Basit bir CSV okuma iş akışını, gerçekçi zorluklarla adım adım geliştirerek sofistike bir biyoinformatik pipeline'a dönüştüreceğiz:

- **Sınırları anlamak:** Veri akışı işlemleri ile betik yazımını birbirinden ayırt etmek ve birlikte nasıl çalıştıklarını kavramak
- **Veri işleme:** Güçlü operatörler kullanarak map'leri ve koleksiyonları çıkarmak, dönüştürmek ve alt kümelere ayırmak
- **Dize işleme:** Regex kalıplarıyla karmaşık dosya adlandırma şemalarını ayrıştırmak ve değişken enterpolasyonunda ustalaşmak
- **Yeniden kullanılabilir fonksiyonlar:** Karmaşık mantığı, daha temiz ve sürdürülebilir iş akışları için adlandırılmış fonksiyonlara taşımak
- **Dinamik mantık:** Farklı girdi türlerine uyum sağlayan süreçler oluşturmak ve dinamik kaynak tahsisi için closure'lar kullanmak
- **Koşullu yönlendirme:** Örnekleri, meta veri özelliklerine göre farklı süreçlere akıllıca yönlendirmek
- **Güvenli işlemler:** Null-safe operatörlerle eksik verileri zarif biçimde ele almak ve açık hata mesajlarıyla girdileri doğrulamak
- **Yapılandırma tabanlı işleyiciler:** Günlükleme, bildirimler ve yaşam döngüsü yönetimi için iş akışı olay işleyicileri kullanmak

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmanız gerekir:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, dosyalarla çalışma, meta veri) rahatça kullanabilmek
- Yaygın programlama yapılarına (değişkenler, map'ler, listeler) temel düzeyde aşinalık

Bu eğitim, programlama kavramlarını karşılaştıkça açıklayacaktır; dolayısıyla kapsamlı programlama deneyimine ihtiyacınız yoktur.
Temel kavramlardan başlayarak ileri düzey kalıplara doğru ilerleyeceğiz.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitime ait dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/essential_scripting_patterns
```

#### Materyalleri inceleyin

Ana iş akışı dosyasını ve örnek veri dosyalarını içeren bir `data` dizinini bulacaksınız.

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

Örnek CSV dosyamız, özelliklerine göre farklı işlem gerektiren biyolojik örneklere ilişkin bilgiler içermektedir:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Gerçek biyoinformatik iş akışlarında karşılaşacağınız pratik programlama tekniklerini keşfetmek için bu gerçekçi veri kümesini kullanacağız.

<!-- TODO: Bunu daha alana özgüsüz hale getirebilir miyiz? -->

<!-- TODO: bir atama ifadesi ekle? #### Görevi inceleyin -->

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
<!-- - [ ] Görevi anlıyorum -->

Tüm kutuları işaretleyebildiyseniz başlayabilirsiniz.

---

## 1. Veri Akışı ve Betik Yazımı: Sınırları Anlamak

### 1.1. Neyin Ne Olduğunu Belirlemek

Nextflow iş akışları yazarken **veri akışı** (verinin kanallar ve süreçler arasında nasıl hareket ettiği) ile **betik yazımı** (veriyi işleyen ve kararlar alan kod) arasındaki farkı gözetmek önemlidir. Birlikte nasıl çalıştıklarını gösteren bir iş akışı oluşturalım.

#### 1.1.1. Temel Nextflow İş Akışı

CSV dosyasını yalnızca okuyan basit bir iş akışıyla başlayın (`main.nf` dosyasında bunu sizin için hazırladık):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

`workflow` bloğu pipeline yapımızı tanımlarken `channel.fromPath()` bir dosya yolundan kanal oluşturur. `.splitCsv()` operatörü CSV dosyasını işler ve her satırı bir map veri yapısına dönüştürür.

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

#### 1.1.2. Map Operatörünü Eklemek

Şimdi veriyi dönüştürmek için betik yazımı ekleyeceğiz; muhtemelen zaten aşina olduğunuz `.map()` operatörünü kullanacağız. Bu operatör, her öğeyi dönüştürmek için kod yazabileceğimiz bir 'closure' alır.

!!! note "Not"

    **Closure**, etrafta dolaştırılıp daha sonra çalıştırılabilen bir kod bloğudur. Satır içinde tanımladığınız bir fonksiyon olarak düşünebilirsiniz. Closure'lar küme parantezleri `{ }` ile yazılır ve parametre alabilir. Nextflow operatörlerinin çalışma biçiminin temelini oluştururlar; bir süredir Nextflow yazıyorsanız, farkında olmadan zaten kullanıyor olabilirsiniz!

Map işleminin nasıl göründüğü aşağıdadır:

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

Bu, ilk **closure**'ımızdır; argüman olarak geçirebileceğiniz anonim bir fonksiyondur (Python'daki lambda'lara veya JavaScript'teki ok fonksiyonlarına benzer). Closure'lar, Nextflow operatörleriyle çalışmak için vazgeçilmezdir.

`{ row -> return row }` closure'ı bir `row` parametresi alır (herhangi bir isim olabilir: `item`, `sample` vb.).

`.map()` operatörü her kanal öğesini işlerken onu closure'ınıza iletir. Burada `row`, her seferinde bir CSV satırını tutar.

Bu değişikliği uygulayın ve iş akışını çalıştırın:

```bash
nextflow run main.nf
```

Girdiyi değiştirmeden döndürdüğümüz için öncekiyle aynı çıktıyı göreceksiniz. Bu, map operatörünün doğru çalıştığını doğrular. Şimdi veriyi dönüştürmeye başlayalım.

#### 1.1.3. Map Veri Yapısı Oluşturmak

Şimdi her veri satırını dönüştürmek için closure'ımızın içine **betik yazımı** mantığı ekleyeceğiz. Veri akışını orkestre etmek yerine bireysel veri öğelerini işlediğimiz yer burasıdır.

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Veri dönüşümü için betik yazımı
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

`sample_meta` map'i, ilgili bilgileri depolayan bir anahtar-değer veri yapısıdır (Python'daki sözlüklere, JavaScript'teki nesnelere veya Ruby'deki hash'lere benzer): örnek kimliği, organizma, doku türü, dizileme derinliği ve kalite skoru.

Verilerimizi temizlemek için `.toLowerCase()` ve `.replaceAll()` gibi dize işleme yöntemlerini, CSV'den gelen dize verilerini uygun sayısal türlere dönüştürmek için ise `.toInteger()` ve `.toDouble()` gibi tür dönüştürme yöntemlerini kullanıyoruz.

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

#### 1.1.4. Koşullu Mantık Eklemek

Şimdi daha fazla betik yazımı ekleyelim; bu sefer veri değerlerine göre kararlar almak için üçlü operatörü kullanacağız.

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

Üçlü operatör, `koşul ? doğruysa_değer : yanlışsa_değer` kalıbını izleyen bir if/else ifadesinin kısaltmasıdır. Bu satır şu anlama gelir: "Kalite 40'tan büyükse 'high', değilse 'normal' kullan." Yakın akrabası olan **Elvis operatörü** (`?:`), bir şey null veya boş olduğunda varsayılan değerler sağlar; bu kalıbı eğitimin ilerleyen bölümlerinde inceleyeceğiz.

Map toplama operatörü `+`, mevcut map'i değiştirmek yerine **yeni bir map** oluşturur. Bu satır, `sample_meta`'daki tüm anahtar-değer çiftlerini ve yeni `priority` anahtarını içeren yeni bir map oluşturur.

!!! Note "Not"

    Closure'lara iletilen map'leri hiçbir zaman değiştirmeyin; her zaman `+` gibi operatörler kullanarak yenilerini oluşturun. Nextflow'da aynı veri çoğunlukla birden fazla işlemden eş zamanlı olarak geçer. Bir map'i yerinde değiştirmek, aynı nesneye başvuran diğer işlemlerde öngörülemeyen yan etkilere yol açabilir. Yeni map'ler oluşturmak, her işlemin kendi temiz kopyasına sahip olmasını sağlar.

Değiştirilen iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Meta verilerimizi kalite skorlarına göre öncelik düzeyi ile zenginleştirmek için koşullu mantığı başarıyla ekledik.

#### 1.1.5. `.subMap()` ile Map'leri Alt Kümelere Ayırmak

`+` operatörü bir map'e anahtar eklerken, bazen tam tersini yapmanız gerekir: yalnızca belirli anahtarları çıkarmak. `.subMap()` yöntemi bunun için mükemmeldir.

Yalnızca kimlik alanlarını içeren, meta verilerimizin basitleştirilmiş bir sürümünü oluşturmak için bir satır ekleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Veri dönüşümü için betik yazımı
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "Yalnızca kimlik alanları: ${id_only}"

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
                // Veri dönüşümü için betik yazımı
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

Değiştirilen iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    Yalnızca kimlik alanları: [id:sample_001, organism:human, tissue:liver]
    Yalnızca kimlik alanları: [id:sample_002, organism:mouse, tissue:brain]
    Yalnızca kimlik alanları: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Bu çıktı, `view()` işlemiyle görüntülenen tam meta veriyi ve `println` ile yazdırdığımız çıkarılmış alt kümeyi birlikte göstermektedir.

`.subMap()` yöntemi bir anahtar listesi alır ve yalnızca bu anahtarları içeren yeni bir map döndürür. Orijinal map'te bir anahtar yoksa sonuca dahil edilmez.

Bu yöntem, farklı süreçler için farklı meta veri sürümleri oluşturmanız gerektiğinde özellikle kullanışlıdır; bazı süreçler tam meta veriye ihtiyaç duyarken diğerleri yalnızca minimal kimlik alanlarına ihtiyaç duyabilir.

Şimdi bu println ifadelerini kaldırarak iş akışınızı önceki durumuna döndürün; ilerleyen bölümlerde bunlara ihtiyacımız olmayacak.

!!! tip "İpucu: Map İşlemleri Özeti"

    - **Anahtar eklemek**: `map1 + [new_key: value]` - Ek anahtarlarla yeni map oluşturur
    - **Anahtar çıkarmak**: `map1.subMap(['key1', 'key2'])` - Yalnızca belirtilen anahtarlarla yeni map oluşturur
    - **Her iki işlem de yeni map'ler oluşturur** - Orijinal map'ler değişmeden kalır

#### 1.1.6. Map'leri Birleştirmek ve Sonuçları Döndürmek

Şimdiye kadar yalnızca Nextflow topluluğunun 'meta map' olarak adlandırdığı yapıyı döndürdük ve bu meta verilerin ilişkili olduğu dosyaları görmezden geldik. Ancak Nextflow iş akışları yazıyorsanız, büyük olasılıkla bu dosyalarla bir şeyler yapmak isteyeceksiniz.

İki öğeden oluşan bir demet içeren bir kanal yapısı çıktısı oluşturalım: zenginleştirilmiş meta veri map'i ve karşılık gelen dosya yolu. Bu, süreçlere veri iletmek için Nextflow'da yaygın bir kalıptır.

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

Bu `[meta, file]` demet yapısı, hem meta veriyi hem de ilişkili dosyaları süreçlere iletmek için Nextflow'da yaygın bir kalıptır.

!!! note "Not"

    **Map'ler ve Meta Veri**: Map'ler, Nextflow'da meta veriyle çalışmanın temelini oluşturur. Meta veri map'leriyle çalışma hakkında daha ayrıntılı bir açıklama için [Meta Veriyle Çalışma](../metadata/) yan görevine bakın.

İş akışımız temel kalıbı göstermektedir: **veri akışı işlemleri** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) verinin pipeline boyunca nasıl hareket ettiğini orkestre ederken, `.map()` closure'ı içindeki **betik yazımı** (map'ler `[key: value]`, dize yöntemleri, tür dönüşümleri, üçlü operatörler) bireysel veri öğelerinin dönüşümünü gerçekleştirir.

### 1.2. Farklı Türleri Anlamak: Kanal ve Liste

Şimdiye kadar iyi gidiyoruz; veri akışı işlemleri ile betik yazımını birbirinden ayırt edebiliyoruz. Peki ya aynı yöntem adı her iki bağlamda da mevcutsa?

Bunun mükemmel bir örneği `collect` yöntemidir; Nextflow standart kütüphanesinde hem kanal türleri hem de List türleri için mevcuttur. Bir List üzerindeki `collect()` yöntemi her öğeyi dönüştürürken, bir kanal üzerindeki `collect()` operatörü tüm kanal yayınlarını tek öğeli bir kanalda toplar.

Bunu örnek verilerle gösterelim; önce kanal `collect()` operatörünün ne yaptığını hatırlayalım. `collect.nf` dosyasına bakın:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - birden fazla kanal yayınını tek bir yayında gruplar
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe tek bir öğede gruplandı)" }
```

Adımlar:

- Örnek kimliklerinden oluşan bir List tanımlayın
- Her örnek kimliğini ayrı ayrı yayınlayan `fromList()` ile bir kanal oluşturun
- Her öğeyi geçerken `view()` ile yazdırın
- Tüm öğeleri kanalın `collect()` operatörüyle tek bir listede toplayın
- Toplanan sonucu (tüm örnek kimliklerini içeren tek öğe) ikinci bir `view()` ile yazdırın

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
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe tek bir öğede gruplandı)
    ```

`view()` her kanal yayını için bir çıktı döndürür; dolayısıyla bu tek çıktının orijinal 3 öğeyi tek bir listede grupladığını biliyoruz.

Şimdi bir List üzerindeki `collect` yöntemini görelim. `collect.nf` dosyasını, orijinal örnek kimliği listesine List'in `collect` yöntemini uygulayacak şekilde değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal yayınını tek bir yayında gruplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe tek bir öğede gruplandı)" }

    // List.collect() - her öğeyi dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()} öğeye dönüştürüldü)"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal yayınını tek bir yayında gruplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe tek bir öğede gruplandı)" }
    ```

Bu yeni kod parçasında:

- Orijinal listedeki her örnek kimliğini dönüştürmek için List'in `collect` yöntemini kullanan yeni bir `formatted_ids` değişkeni tanımlıyoruz
- Sonucu `println` ile yazdırıyoruz

Değiştirilen iş akışını çalıştırın:

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
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe tek bir öğede gruplandı)
    ```

Bu sefer verinin yapısını DEĞİŞTİRMEDİK; listede hâlâ 3 öğe var. Ancak List'in `collect` yöntemiyle her öğeyi DÖNÜŞTÜRDÜKve değiştirilmiş değerlere sahip yeni bir liste elde ettik. Bu, bir kanal üzerinde `map` operatörü kullanmaya benzer; ancak kanal yerine bir List veri yapısı üzerinde çalışmaktadır.

Burada bir noktayı vurgulamak için `collect`'in aşırı bir örneğini kullandık. Temel ders şudur: İş akışları yazarken her zaman **veri yapıları** (List'ler, Map'ler vb.) ile **kanallar** (veri akışı yapıları) arasındaki farkı gözetin. İşlemler aynı adı paylaşabilir, ancak çağrıldıkları türe bağlı olarak tamamen farklı davranabilir.

### 1.3. Spread Operatörü (`*.`) - Özellik Çıkarma için Kısayol

List'in `collect` yöntemiyle ilişkili olan spread operatörü (`*.`), koleksiyonlardan özellik çıkarmak için kısa ve öz bir yol sunar. Özünde yaygın bir `collect` kalıbı için sözdizimsel bir şekerdir.

`collect.nf` dosyamıza bir gösterim ekleyelim:

=== "Sonra"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal yayınını tek bir yayında gruplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe tek bir öğede gruplandı)" }

    // List.collect() - her öğeyi dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()} öğeye dönüştürüldü)"

    // Spread operatörü - kısa ve öz özellik erişimi
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operatörü sonucu: ${all_ids}"
    ```

=== "Önce"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - birden fazla kanal yayınını tek bir yayında gruplar
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Bireysel kanal öğesi: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() sonucu: ${list} (${list.size()} öğe tek bir öğede gruplandı)" }

    // List.collect() - her öğeyi dönüştürür, yapıyı korur
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() sonucu: ${formatted_ids} (${sample_ids.size()} öğe ${formatted_ids.size()} öğeye dönüştürüldü)"
    ```

Güncellenmiş iş akışını çalıştırın:

```bash title="Spread operatörünü test edin"
nextflow run collect.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() sonucu: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 öğe 3 öğeye dönüştürüldü)
    Spread operatörü sonucu: [s1, s2, s3]
    Bireysel kanal öğesi: sample_001
    Bireysel kanal öğesi: sample_002
    Bireysel kanal öğesi: sample_003
    channel.collect() sonucu: [sample_001, sample_002, sample_003] (3 öğe tek bir öğede gruplandı)
    ```

Spread operatörü `*.`, yaygın bir collect kalıbı için kısayoldur:

```groovy
// Bunlar eşdeğerdir:
def ids = samples*.id
def ids = samples.collect { it.id }

// Yöntem çağrılarıyla da çalışır:
def names = files*.getName()
def names = files.collect { it.getName() }
```

Spread operatörü, bir nesne listesinden tek bir özellik çıkarmanız gerektiğinde özellikle kullanışlıdır; tam `collect` closure'ını yazmaktan daha okunabilirdir.

!!! tip "İpucu: Spread ve Collect'i Ne Zaman Kullanmalı"

    - **Spread (`*.`) kullanın** basit özellik erişimi için: `samples*.id`, `files*.name`
    - **Collect kullanın** dönüşümler veya karmaşık mantık için: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Özetle

Bu bölümde şunları öğrendiniz:

- **Veri akışı ve betik yazımı**: Kanal operatörleri verinin pipeline'ınızda nasıl aktığını orkestre ederken, betik yazımı bireysel veri öğelerini dönüştürür
- **Türleri anlamak**: Aynı yöntem adı (örneğin `collect`), çağrıldığı türe bağlı olarak (Kanal ve List) farklı davranabilir
- **Bağlam önemlidir**: Kanallarla (veri akışı) mı yoksa veri yapılarıyla (betik yazımı) mı çalıştığınızın her zaman farkında olun

Bu sınırları anlamak, hata ayıklama, belgeleme ve sürdürülebilir iş akışları yazma açısından temel öneme sahiptir.

Sırada, gerçek dünya verilerini işlemek için vazgeçilmez olan dize işleme yeteneklerini daha ayrıntılı inceleyeceğiz.

---

## 2. Dize İşleme ve Dinamik Betik Oluşturma

Dize işlemede ustalaşmak, kırılgan iş akışlarını dayanıklı pipeline'lardan ayırır. Bu bölüm, karmaşık dosya adlarını ayrıştırmayı, dinamik betik oluşturmayı ve değişken enterpolasyonunu kapsamaktadır.

### 2.1. Kalıp Eşleştirme ve Düzenli İfadeler

Biyoinformatik dosyaları çoğunlukla meta veriyi kodlayan karmaşık adlandırma kurallarına sahiptir. Düzenli ifadelerle kalıp eşleştirme kullanarak bu bilgileri otomatik olarak çıkaralım.

`main.nf` iş akışımıza dönerek dosya adlarından ek örnek bilgisi çıkarmak için kalıp eşleştirme mantığı ekleyeceğiz. Veri kümemizdeki FASTQ dosyaları, `SAMPLE_001_S1_L001_R1_001.fastq.gz` gibi adlara sahip Illumina tarzı adlandırma kurallarını izlemektedir. Bu adlar şifreli görünebilir, ancak örnek kimliği, şerit numarası ve okuma yönü gibi kullanışlı meta verileri kodlamaktadır. Bu adları ayrıştırmak için regex yeteneklerini kullanacağız.

Mevcut `main.nf` iş akışınızda aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Veri dönüşümü için betik yazımı
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
                // Veri dönüşümü için betik yazımı
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

Bu, temel **dize işleme kavramlarını** göstermektedir:

1. `~/pattern/` sözdizimini kullanan **düzenli ifade değişmezleri** - ters eğik çizgilerden kaçmaya gerek kalmadan regex kalıbı oluşturur
2. `=~` operatörüyle **kalıp eşleştirme** - bir dizeyi regex kalıbıyla eşleştirmeye çalışır
3. `[0][1]`, `[0][2]` vb. ile grupları yakalayan **Matcher nesneleri** - `[0]` tüm eşleşmeyi, `[1]`, `[2]` vb. parantez içindeki yakalanan grupları ifade eder

`^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` regex kalıbını inceleyelim:

| Kalıp               | Eşleşir                                    | Yakalar                                    |
| ------------------- | ------------------------------------------ | ------------------------------------------ |
| `^(.+)`             | Baştan örnek adı                           | Grup 1: örnek adı                          |
| `_S(\d+)`           | Örnek numarası `_S1`, `_S2` vb.            | Grup 2: örnek numarası                     |
| `_L(\d{3})`         | Şerit numarası `_L001`                     | Grup 3: şerit (3 basamak)                  |
| `_(R[12])`          | Okuma yönü `_R1` veya `_R2`               | Grup 4: okuma yönü                         |
| `_(\d{3})`          | Parça numarası `_001`                      | Grup 5: parça (3 basamak)                  |
| `\.fastq(?:\.gz)?$` | Dosya uzantısı `.fastq` veya `.fastq.gz`   | Yakalanmaz (?:, yakalamayan gruptur)       |

Bu, meta veriyi otomatik olarak çıkarmak için Illumina tarzı adlandırma kurallarını ayrıştırır.

Değiştirilen iş akışını çalıştırın:

```bash title="Kalıp eşleştirmeyi test edin"
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

Bu çıktı, dosya adlarından zenginleştirilmiş meta veriyi göstermektedir.

### 2.2. Süreçlerde Dinamik Betik Oluşturma

Süreç betik blokları, özünde kabuğa iletilen çok satırlı dizelerdir. Girdi özelliklerine göre farklı betik dizeleri dinamik olarak oluşturmak için **koşullu mantık** (if/else, üçlü operatörler) kullanabilirsiniz. Bu, süreç tanımlarını çoğaltmadan tek uçlu ve çift uçlu dizileme okumaları gibi farklı girdi türlerini işlemek için vazgeçilmezdir.

İş akışımıza bu kalıbı gösteren bir süreç ekleyelim. `modules/fastp.nf` dosyasını açın ve inceleyin:

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

Süreç, FASTQ dosyalarını girdi olarak alır ve adaptörleri kırpmak ve düşük kaliteli okumaları filtrelemek için `fastp` aracını çalıştırır. Ne yazık ki bu süreci yazan kişi, örnek veri kümemizdeki tek uçlu okumaları dikkate almamış. Süreci iş akışımıza ekleyelim ve ne olduğunu görelim:

Önce modülü `main.nf` iş akışınızın en başına dahil edin:

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

Bu değiştirilen iş akışını çalıştırın:

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

Sürecin ikinci girdi dosyası için `null` değeriyle `fastp` çalıştırmaya çalıştığını ve bu nedenle başarısız olduğunu görebilirsiniz. Bunun nedeni, veri kümemizin tek uçlu okumalar içermesi, ancak sürecin çift uçlu okumalar (aynı anda iki girdi dosyası) beklediği şekilde sabit kodlanmış olmasıdır.

Bunu düzeltmek için `FASTP` sürecinin `script:` bloğuna koşullu mantık ekleyin. Bir if/else ifadesi, okuma dosyası sayısını kontrol eder ve komutu buna göre ayarlar.

=== "Sonra"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Basit tek uçlu ve çift uçlu algılama
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

Artık iş akışı hem tek uçlu hem de çift uçlu okumaları sorunsuz biçimde işleyebilir. Koşullu mantık, girdi dosyalarının sayısını kontrol eder ve `fastp` için uygun komutu oluşturur. Çalışıp çalışmadığını görelim:

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

Harika görünüyor! Çalıştırılan gerçek komutları kontrol edersek (görev hash'inize göre özelleştirin):

```console title="Çalıştırılan komutları kontrol edin"
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

Dinamik betik mantığının başka bir yaygın kullanımı [Nextflow for Science Genomics modülünde](../../nf4science/genomics/02_joint_calling) görülebilir. O modülde çağrılan GATK süreci birden fazla girdi dosyası alabilir, ancak doğru bir komut satırı oluşturmak için her birinin başına `-V` eklenmelidir. Süreç, bir girdi dosyaları koleksiyonunu (`all_gvcfs`) doğru komut argümanlarına dönüştürmek için betik yazımını kullanır:

```groovy title="GATK için komut satırı işleme" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Süreç betik bloklarında betik yazımı kullanmanın bu kalıpları son derece güçlüdür ve pek çok senaryoda uygulanabilir; değişken girdi türlerini işlemekten dosya koleksiyonlarından karmaşık komut satırı argümanları oluşturmaya kadar, süreçlerinizi gerçek dünyanın çeşitli gereksinimlerine gerçek anlamda uyarlanabilir kılar.

### 2.3. Değişken Enterpolasyonu: Nextflow ve Kabuk Değişkenleri

Süreç betikleri, Nextflow değişkenlerini, kabuk değişkenlerini ve komut ikamelerini bir arada kullanır; her birinin farklı enterpolasyon sözdizimi vardır. Yanlış sözdizimi kullanmak hatalara yol açar. Bunu, bir işleme raporu oluşturan bir süreçle inceleyelim.

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
    echo "İşleniyor: ${reads}" > ${meta.id}_report.txt
    echo "Örnek: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Bu süreç, örnek kimliği ve dosya adıyla basit bir rapor yazar. Şimdi farklı türde değişkenleri bir arada kullanmamız gerektiğinde ne olduğunu görmek için çalıştıralım.

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

Şimdi iş akışını çalıştırın ve `results/reports/` dizinindeki oluşturulan raporları kontrol edin. Her örnek hakkında temel bilgiler içermeleri gerekir.

<!-- TODO: çalıştırma komutunu ekle -->

??? success "Komut çıktısı"

    ```console
    <!-- TODO: çıktı -->
    ```

Peki işlemenin ne zaman ve nerede gerçekleştiğine dair bilgi eklemek istersek? Rapora geçerli kullanıcıyı, ana bilgisayar adını ve tarihi dahil etmek için **kabuk** değişkenlerini ve komut ikamesini kullanalım:

=== "Sonra"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "İşleniyor: ${reads}" > ${meta.id}_report.txt
        echo "Örnek: ${meta.id}" >> ${meta.id}_report.txt
        echo "İşleyen: ${USER}" >> ${meta.id}_report.txt
        echo "Ana bilgisayar: $(hostname)" >> ${meta.id}_report.txt
        echo "Tarih: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Önce"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "İşleniyor: ${reads}" > ${meta.id}_report.txt
        echo "Örnek: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Bunu çalıştırırsanız bir hata fark edeceksiniz; Nextflow, `${USER}`'ı mevcut olmayan bir Nextflow değişkeni olarak yorumlamaya çalışır.

??? failure "Komut çıktısı"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Bash'in bunu işleyebilmesi için kaçış karakteri eklememiz gerekiyor.

Kabuk değişkenlerini ve komut ikamelerini ters eğik çizgi (`\`) ile kaçırarak düzeltin:

=== "Sonra"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "İşleniyor: ${reads}" > ${meta.id}_report.txt
        echo "Örnek: ${meta.id}" >> ${meta.id}_report.txt
        echo "İşleyen: \${USER}" >> ${meta.id}_report.txt
        echo "Ana bilgisayar: \$(hostname)" >> ${meta.id}_report.txt
        echo "Tarih: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Önce"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "İşleniyor: ${reads}" > ${meta.id}_report.txt
        echo "Örnek: ${meta.id}" >> ${meta.id}_report.txt
        echo "İşleyen: ${USER}" >> ${meta.id}_report.txt
        echo "Ana bilgisayar: $(hostname)" >> ${meta.id}_report.txt
        echo "Tarih: $(date)" >> ${meta.id}_report.txt
        """
    ```

Artık çalışıyor! Ters eğik çizgi (`\`), Nextflow'a "bunu yorumlama, Bash'e ilet" der.

### Özetle

Bu bölümde **dize işleme** tekniklerini öğrendiniz:

- **Dosya ayrıştırma için düzenli ifadeler**: Karmaşık dosya adlandırma kurallarından meta veri çıkarmak için `=~` operatörünü ve regex kalıplarını (`~/pattern/`) kullanmak
- **Dinamik betik oluşturma**: Girdi özelliklerine göre farklı betik dizeleri oluşturmak için koşullu mantık (if/else, üçlü operatörler) kullanmak
- **Değişken enterpolasyonu**: Nextflow'un dizeleri ne zaman yorumladığını, kabuğun ne zaman yorumladığını anlamak
  - `${var}` - Nextflow değişkenleri (iş akışı derleme zamanında Nextflow tarafından enterpolasyon yapılır)
  - `\${var}` - Kabuk ortam değişkenleri (kaçırılmış, çalışma zamanında bash'e iletilir)
  - `\$(cmd)` - Kabuk komut ikamesi (kaçırılmış, çalışma zamanında bash tarafından çalıştırılır)

Bu dize işleme ve oluşturma kalıpları, gerçek dünya biyoinformatik iş akışlarında karşılaşacağınız çeşitli dosya formatlarını ve adlandırma kurallarını işlemek için vazgeçilmezdir.

---

## 3. Yeniden Kullanılabilir Fonksiyonlar Oluşturmak

Kanal operatörlerinde veya süreç tanımlarında satır içi karmaşık iş akışı mantığı, okunabilirliği ve sürdürülebilirliği azaltır. **Fonksiyonlar**, bu mantığı adlandırılmış, yeniden kullanılabilir bileşenlere taşımanızı sağlar.

Map işlemimiz uzun ve karmaşık bir hale geldi. `def` anahtar sözcüğünü kullanarak bunu yeniden kullanılabilir bir fonksiyona çıkaralım.

Mevcut iş akışımızla bunun nasıl göründüğünü göstermek için aşağıdaki değişikliği yapın; `separateMetadata` adlı yeniden kullanılabilir bir fonksiyon tanımlamak için `def` kullanın:

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

Bu mantığı bir fonksiyona çıkararak gerçek iş akışı mantığını çok daha temiz bir hale getirdik:

```groovy title="minimal iş akışı"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Bu, iş akışı mantığını bir bakışta okumayı ve anlamayı çok daha kolay kılar. `separateMetadata` fonksiyonu, meta veriyi ayrıştırma ve zenginleştirmeye yönelik tüm karmaşık mantığı kapsüller; böylece yeniden kullanılabilir ve test edilebilir hale gelir.

İş akışının hâlâ çalıştığından emin olmak için çalıştırın:

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

Çıktı, her iki sürecin de başarıyla tamamlandığını göstermelidir. İş akışı artık çok daha temiz ve sürdürülebilir; tüm karmaşık meta veri işleme mantığı `separateMetadata` fonksiyonunda kapsüllenmiştir.

### Özetle

Bu bölümde **fonksiyon oluşturmayı** öğrendiniz:

- **`def` ile fonksiyon tanımlamak**: Adlandırılmış fonksiyonlar oluşturmak için kullanılan anahtar sözcük (Python'daki `def` veya JavaScript'teki `function` gibi)
- **Fonksiyon kapsamı**: Betik düzeyinde tanımlanan fonksiyonlar, Nextflow iş akışınız boyunca erişilebilirdir
- **Dönüş değerleri**: Fonksiyonlar son ifadeyi otomatik olarak döndürür ya da açık `return` kullanılabilir
- **Daha temiz kod**: Karmaşık mantığı fonksiyonlara çıkarmak, herhangi bir dilde temel bir yazılım mühendisliği uygulamasıdır

Sırada, dinamik kaynak tahsisi için süreç yönergelerinde closure'ların nasıl kullanılacağını inceleyeceğiz.

---

## 4. Closure'larla Dinamik Kaynak Yönergeleri

Şimdiye kadar süreçlerin `script` bloğunda betik yazımı kullandık. Ancak **closure'lar** (Bölüm 1.1'de tanıtıldı), özellikle dinamik kaynak tahsisi için süreç yönergelerinde de son derece kullanışlıdır. FASTP sürecimize, örnek özelliklerine göre uyum sağlayan kaynak yönergeleri ekleyelim.

### 4.1. Örneğe Özgü Kaynak Tahsisi

Şu anda FASTP sürecimiz varsayılan kaynakları kullanmaktadır. Yüksek derinlikli örnekler için daha fazla CPU tahsis ederek daha akıllı hale getirelim. `modules/fastp.nf` dosyasını düzenleyerek dinamik bir `cpus` yönergesi ve statik bir `memory` yönergesi ekleyin:

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

`{ meta.depth > 40000000 ? 2 : 1 }` closure'ı **üçlü operatörü** (Bölüm 1.1'de ele alındı) kullanır ve her görev için değerlendirilerek örnek başına kaynak tahsisine olanak tanır. Yüksek derinlikli örnekler (>40M okuma) 2 CPU alırken diğerleri 1 CPU alır.

!!! note "Not: Yönergelerde Girdi Değişkenlerine Erişim"

    Closure, Nextflow bu closure'ları her görev yürütme bağlamında değerlendirdiğinden, herhangi bir girdi değişkenine (burada `meta` gibi) erişebilir.

Görev hash'lerini daha kolay görmek için `-ansi-log false` seçeneğiyle iş akışını yeniden çalıştırın.

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

Herhangi bir görev için CPU tahsisini görmek amacıyla çalıştırılan `docker` komutunu kontrol edebilirsiniz:

```console title="Docker komutunu kontrol edin"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Şuna benzer bir şey görmelisiniz:

```bash title="docker komutu"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Bu örnekte yüksek derinlikli bir örnek olduğu için 2 CPU talep eden bir örnek seçtik (`--cpu-shares 2048`); ancak örnek derinliğine bağlı olarak farklı CPU tahsisleri görmelisiniz. Diğer görevler için de bunu deneyin.

### 4.2. Yeniden Deneme Stratejileri

Bir diğer güçlü kalıp, yeniden deneme stratejileri için `task.attempt` kullanmaktır. Bunun neden yararlı olduğunu göstermek için önce FASTP'a ihtiyaç duyduğundan daha az bellek tahsis edeceğiz. `modules/fastp.nf` dosyasındaki `memory` yönergesini `1.GB` olarak değiştirin:

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

... ve iş akışını yeniden çalıştırın:

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

Bu, sürecin bellek sınırını aştığı için sonlandırıldığını göstermektedir.

Bu, gerçek dünya iş akışlarında çok yaygın bir senaryodur; bazen bir görevin ne kadar belleğe ihtiyaç duyacağını çalıştırana kadar bilemezsiniz.

İş akışımızı daha dayanıklı hale getirmek için, yine bir Groovy closure kullanarak her denemede bellek tahsisini artıran bir yeniden deneme stratejisi uygulayabiliriz. `memory` yönergesini temel belleği `task.attempt` ile çarpacak şekilde değiştirin ve `errorStrategy 'retry'` ile `maxRetries 2` yönergelerini ekleyin:

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

Artık süreç yetersiz bellek nedeniyle başarısız olursa Nextflow daha fazla bellekle yeniden dener:

- Birinci deneme: 1 GB (task.attempt = 1)
- İkinci deneme: 2 GB (task.attempt = 2)

... ve `maxRetries` sınırına kadar bu şekilde devam eder.

### Özetle

Closure'larla dinamik yönergeler şunları yapmanızı sağlar:

- Girdi özelliklerine göre kaynak tahsis etmek
- Artan kaynaklarla otomatik yeniden deneme stratejileri uygulamak
- Birden fazla faktörü (meta veri, deneme sayısı, öncelikler) bir arada kullanmak
- Karmaşık kaynak hesaplamaları için koşullu mantık kullanmak

Bu, iş akışlarınızı hem daha verimli (aşırı tahsis yapmadan) hem de daha dayanıklı (daha fazla kaynakla otomatik yeniden deneme) kılar.

---

## 5. Koşullu Mantık ve Süreç Kontrolü

Daha önce kanal verilerini dönüştürmek için betik yazımıyla `.map()` kullandık. Şimdi verilere göre hangi süreçlerin çalışacağını kontrol etmek için koşullu mantık kullanacağız; bu, farklı örnek türlerine uyum sağlayan esnek iş akışları için vazgeçilmezdir.

Nextflow'un [veri akışı operatörleri](https://www.nextflow.io/docs/latest/reference/operator.html), çalışma zamanında değerlendirilen closure'lar alır; bu sayede koşullu mantık, kanal içeriğine göre iş akışı kararlarını yönlendirebilir.

### 5.1. `.branch()` ile Yönlendirme

Örneğin, dizileme örneklerimizin yalnızca belirli bir eşiğin üzerinde kapsama sahip insan örnekleri olması durumunda FASTP ile kırpılması gerektiğini varsayalım. Fare örnekleri veya düşük kapsama sahip örnekler bunun yerine Trimgalore ile çalıştırılmalıdır (bu yapay bir örnektir, ancak konuyu iyi açıklamaktadır).

`modules/trimgalore.nf` dosyasında basit bir Trimgalore süreci sağladık; isterseniz inceleyebilirsiniz, ancak ayrıntılar bu alıştırma için önemli değildir. Önemli olan, örnekleri meta verilerine göre yönlendirmek istediğimizdir.

`modules/trimgalore.nf` dosyasından yeni modülü dahil edin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... ardından `main.nf` iş akışınızı, örnekleri meta verilerine göre dallandıracak ve uygun kırpma sürecine yönlendirecek şekilde değiştirin:

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

Bu değiştirilen iş akışını çalıştırın:

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

Burada, örnekleri meta verilerine göre yönlendirmek için `.branch{}` operatörü içinde küçük ama güçlü koşullu ifadeler kullandık. Yüksek kapsama sahip insan örnekleri `FASTP`'tan geçerken diğer tüm örnekler `TRIMGALORE`'dan geçer.

### 5.2. Doğruluk Değeriyle `.filter()` Kullanmak

İş akışı yürütmesini kontrol etmek için bir diğer güçlü kalıp, pipeline boyunca hangi öğelerin devam etmesi gerektiğini belirlemek için closure kullanan `.filter()` operatörüdür. Filter closure'ının içine, hangi öğelerin geçeceğine karar veren **boolean ifadeler** yazarsınız.

Nextflow (pek çok dinamik dil gibi), boolean bağlamlarda hangi değerlerin `true` veya `false` olarak değerlendirileceğini belirleyen bir **"doğruluk değeri"** kavramına sahiptir:

- **Doğru (truthy)**: Null olmayan değerler, boş olmayan dizeler, sıfır olmayan sayılar, boş olmayan koleksiyonlar
- **Yanlış (falsy)**: `null`, boş dizeler `""`, sıfır `0`, boş koleksiyonlar `[]` veya `[:]`, `false`

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

İş akışını yeniden çalıştırın:

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

Bazı örnekleri dışlayan bir filtre seçtiğimiz için daha az görev çalıştırıldı.

`meta.id && meta.organism && meta.depth >= 25000000` filtre ifadesi, doğruluk değerini açık karşılaştırmalarla birleştirir:

- `meta.id && meta.organism`, her iki alanın da var olduğunu ve boş olmadığını kontrol eder (doğruluk değeri kullanarak)
- `meta.depth >= 25000000`, açık bir karşılaştırmayla yeterli dizileme derinliğini sağlar

!!! note "Not: Pratikte Doğruluk Değeri"

    `meta.id && meta.organism` ifadesi şunu yazmaktan çok daha kısadır:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Bu, filtreleme mantığını çok daha temiz ve okunması kolay hale getirir.

### Özetle

Bu bölümde, `.branch{}` ve `.filter{}` gibi Nextflow operatörlerinin closure arayüzlerini kullanarak iş akışı yürütmesini kontrol etmek için koşullu mantığı nasıl kullanacağınızı ve kısa koşullu ifadeler yazmak için doğruluk değerinden nasıl yararlanacağınızı öğrendiniz.

Pipeline'ımız artık örnekleri uygun süreçlere akıllıca yönlendiriyor; ancak üretim iş akışlarının geçersiz verileri zarif biçimde işlemesi gerekir. İş akışımızı eksik veya null değerlere karşı dayanıklı hale getirelim.

---

## 6. Güvenli Gezinme ve Elvis Operatörleri

`separateMetadata` fonksiyonumuz şu anda tüm CSV alanlarının mevcut ve geçerli olduğunu varsaymaktadır. Peki eksik verilerle ne olur? Bunu öğrenelim.

### 6.1. Sorun: Var Olmayan Özelliklere Erişmek

Diyelim ki isteğe bağlı dizileme çalışması bilgisi için destek eklemek istiyoruz. Bazı laboratuvarlarda örneklerin dizileme çalışması kimliği veya toplu numarası için ek bir alanı olabilir, ancak mevcut CSV'mizde bu sütun yok. Yine de erişmeye çalışalım.

`separateMetadata` fonksiyonunu bir `run_id` alanı içerecek şekilde değiştirin:

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

Sorun şu: `row.run_id`, CSV'mizde `run_id` sütunu olmadığı için `null` döndürüyor. `null` üzerinde `.toUpperCase()` çağırmaya çalıştığımızda çöküyor. İşte burada güvenli gezinme operatörü devreye giriyor.

### 6.2. Güvenli Gezinme Operatörü (`?.`)

Güvenli gezinme operatörü (`?.`), `null` bir değer üzerinde çağrıldığında istisna fırlatmak yerine `null` döndürür. `?.` öncesindeki nesne `null` ise, yöntem çalıştırılmadan tüm ifade `null` olarak değerlendirilir.

Güvenli gezinme kullanmak için fonksiyonu güncelleyin:

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

Yeniden çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    <!-- TODO: çıktı -->
    ```

Çökme yok! İş akışı artık eksik alanı zarif biçimde işliyor. `row.run_id` `null` olduğunda, `?.` operatörü `.toUpperCase()` çağrısını engeller ve `run_id`, istisna fırlatmak yerine `null` olur.

### 6.3. Varsayılanlar için Elvis Operatörü (`?:`)

Elvis operatörü (`?:`), sol taraf "falsy" (daha önce açıklandığı gibi) olduğunda varsayılan değerler sağlar. Adını Elvis Presley'den alır; çünkü `?:` yana çevrildiğinde onun ünlü saçına ve gözlerine benzer!

Artık güvenli gezinme kullandığımıza göre, bu alana sahip olmayan örnekler için `run_id` `null` olacak. Varsayılan bir değer sağlamak ve bunu `sample_meta` map'imize eklemek için Elvis operatörünü kullanalım:

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

Sonuçları görmek için iş akışına bir `view()` operatörü de ekleyin:

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

Mükemmel! Artık tüm örnekler, gerçek çalışma kimliklerini (büyük harfle) veya varsayılan 'UNSPECIFIED' değerini içeren bir `run` alanına sahip. `?.` ve `?:` kombinasyonu hem güvenlik (çökme yok) hem de mantıklı varsayılanlar sağlar.

Çalıştığını doğruladıktan sonra `.view()` operatörünü kaldırın.

!!! tip "İpucu: Güvenli Gezinme ve Elvis'i Birleştirmek"

    `value?.method() ?: 'default'` kalıbı üretim iş akışlarında yaygındır:

    - `value?.method()` - Yöntemi güvenle çağırır, `value` `null` ise `null` döndürür
    - `?: 'default'` - Sonuç `null` ise yedek değer sağlar

    Bu kalıp, eksik/tamamlanmamış verileri zarif biçimde işler.

Bu operatörleri fonksiyonlarda, operatör closure'larında (`.map{}`, `.filter{}`), süreç betiklerinde ve yapılandırma dosyalarında tutarlı biçimde kullanın. Gerçek dünya verilerini işlerken çökmeleri önlerler.

### Özetle

- **Güvenli gezinme (`?.`)**: Null değerlerde çökmeleri önler; istisna fırlatmak yerine null döndürür
- **Elvis operatörü (`?:`)**: Varsayılanlar sağlar - `value ?: 'default'`
- **Birleştirme**: `value?.method() ?: 'default'` yaygın kalıptır

Bu operatörler, iş akışlarını eksik verilere karşı dayanıklı kılar; gerçek dünya çalışmaları için vazgeçilmezdir.

---

## 7. `error()` ve `log.warn` ile Doğrulama

Bazen girdi parametreleri geçersizse iş akışını hemen durdurmak gerekir. Nextflow'da doğrulama mantığı uygulamak için `error()` ve `log.warn` gibi yerleşik fonksiyonların yanı sıra `if` ifadeleri ve boolean mantığı gibi standart programlama yapılarını kullanabilirsiniz. İş akışımıza doğrulama ekleyelim.

İş akışı bloğunuzdan önce bir doğrulama fonksiyonu oluşturun, iş akışından çağırın ve CSV dosya yolu için bir parametre kullanmak üzere kanal oluşturmayı değiştirin. Parametre eksikse veya dosya yoksa, açık bir mesajla yürütmeyi durdurmak için `error()` çağırın.

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Girdi parametresinin sağlandığını kontrol et
        if (!params.input) {
            error("Girdi CSV dosya yolu belirtilmedi. Lütfen --input <file.csv> ile belirtin")
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
    Girdi CSV dosya yolu belirtilmedi. Lütfen --input <file.csv> ile belirtin
    ```

İş akışı, daha sonra gizemli bir şekilde başarısız olmak yerine açık bir hata mesajıyla hemen durur.

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
    <!-- TODO: çıktı -->
    ```

Bu sefer başarıyla çalışır.

`separateMetadata` fonksiyonu içine de doğrulama ekleyebilirsiniz. Düşük dizileme derinliğine sahip örnekler için uyarı vermek, ancak iş akışının devam etmesine izin vermek için önemli olmayan `log.warn` kullanın:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Verinin mantıklı olduğunu doğrula
        if (sample_meta.depth < 30000000) {
            log.warn "Düşük dizileme derinliği ${sample_meta.id} için: ${sample_meta.depth}"
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

İş akışını orijinal CSV ile yeniden çalıştırın:

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
    WARN: Düşük dizileme derinliği sample_002 için: 25000000
    ```

Örneklerden biri için düşük dizileme derinliği uyarısı görüyoruz.

### Özetle

- **`error()`**: İş akışını açık mesajla hemen durdurur
- **`log.warn`**: İş akışını durdurmadan uyarı verir
- **Erken doğrulama**: Hızlı başarısız olmak ve yardımcı hatalar vermek için işlemeden önce girdileri kontrol edin
- **Doğrulama fonksiyonları**: İş akışı başlangıcında çağrılabilen yeniden kullanılabilir doğrulama mantığı oluşturun

Doğru doğrulama, sorunları erken ve açık hata mesajlarıyla yakalayarak iş akışlarını daha dayanıklı ve kullanıcı dostu hale getirir.

---

## 8. İş Akışı Olay İşleyicileri

Şimdiye kadar iş akışı betiklerimizde ve süreç tanımlarımızda kod yazdık. Ancak bilmeniz gereken bir özellik daha var: iş akışı olay işleyicileri.

Olay işleyicileri, iş akışınızın yaşam döngüsündeki belirli noktalarda çalışan closure'lardır. Günlükleme, bildirimler veya temizleme işlemleri eklemek için mükemmeldir. Bu işleyiciler, iş akışı tanımınızın yanında iş akışı betiğinizde tanımlanmalıdır.

### 8.1. `onComplete` İşleyicisi

En yaygın kullanılan olay işleyicisi, iş akışınız tamamlandığında (başarılı veya başarısız olsun) çalışan `onComplete`'dir. Pipeline sonuçlarımızı özetlemek için bir tane ekleyelim.

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
            println "Süre       : ${workflow.duration}"
            println "Başarı     : ${workflow.success}"
            println "workDir    : ${workflow.workDir}"
            println "Çıkış kodu : ${workflow.exitStatus}"
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

Bu closure, iş akışı tamamlandığında çalışır. İçinde, yürütme hakkında yararlı özellikler sağlayan `workflow` nesnesine erişebilirsiniz.

İş akışınızı çalıştırın; sonunda bu özeti göreceksiniz!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Düşük dizileme derinliği sample_002 için: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline yürütme özeti:
    ==========================
    Tamamlandı: 2025-10-10T12:14:24.885384+01:00
    Süre       : 2.9s
    Başarı     : true
    workDir    : /workspaces/training/side-quests/essential_scripting_patterns/work
    Çıkış kodu : 0
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
            println "Süre       : ${workflow.duration}"
            println "Başarı     :```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline yürütme özeti:"
            println "=========================="
            println "Tamamlandı: ${workflow.complete}"
            println "Süre       : ${workflow.duration}"
            println "Başarı     : ${workflow.success}"
            println "workDir    : ${workflow.workDir}"
            println "Çıkış kodu : ${workflow.exitStatus}"
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
            println "Süre       : ${workflow.duration}"
            println "Başarı     : ${workflow.success}"
            println "workDir    : ${workflow.workDir}"
            println "Çıkış kodu : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Artık belirtilmişse çıktı dizinini ve başarı/başarısızlık mesajını da içeren çok daha bilgilendirici bir özet elde ediyoruz:

<!-- TODO: çalıştırma komutunu ekle -->

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Düşük dizileme derinliği sample_002 için: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline yürütme özeti:
    ==========================
    Tamamlandı: 2025-10-10T12:16:00.522569+01:00
    Süre       : 3.6s
    Başarı     : true
    workDir    : /workspaces/training/side-quests/essential_scripting_patterns/work
    Çıkış kodu : 0

    ✅ Pipeline başarıyla tamamlandı!
    ```

Dosya işlemlerini kullanarak özeti bir dosyaya da yazabilirsiniz:

```groovy title="main.nf - Özeti dosyaya yazmak"
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

        // Günlük dosyasına yaz
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. `onError` İşleyicisi

`onComplete`'in yanı sıra kullanabileceğiniz bir olay işleyicisi daha vardır: yalnızca iş akışı başarısız olduğunda çalışan `onError`:

```groovy title="main.nf - onError işleyicisi"
workflow {
    // ... iş akışı kodunuz ...

    workflow.onError = {
        println "="* 50
        println "Pipeline yürütmesi başarısız oldu!"
        println "Hata mesajı: ${workflow.errorMessage}"
        println "="* 50

        // Ayrıntılı hata günlüğü yaz
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

```groovy title="main.nf - Birleşik işleyiciler"
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

### Özetle

Bu bölümde şunları öğrendiniz:

- **Olay işleyici closure'ları**: İş akışı betiğinizdeki, farklı yaşam döngüsü noktalarında çalışan closure'lar
- **`onComplete` işleyicisi**: Yürütme özetleri ve sonuç raporlaması için
- **`onError` işleyicisi**: Hata işleme ve başarısızlıkları günlükleme için
- **Workflow nesnesi özellikleri**: `workflow.success`, `workflow.duration`, `workflow.errorMessage` vb. erişimi

Olay işleyicileri, iş akışı betiklerinize sofistike günlükleme ve bildirim yetenekleri eklemek için Nextflow dilinin tüm gücünü nasıl kullanabileceğinizi göstermektedir.

---

## Özet

Tebrikler, başardınız!

Bu yan görev boyunca, temel meta veri işlemeden başlayarak sofistike, üretime hazır bir iş akışına evrilen kapsamlı bir örnek işleme pipeline'ı oluşturdunuz.
Her bölüm bir öncekinin üzerine inşa edildi; programlama yapılarının basit iş akışlarını güçlü veri işleme sistemlerine nasıl dönüştürdüğünü şu faydalarla gösterdi:

- **Daha net kod**: Veri akışı ve betik yazımını anlamak, daha düzenli iş akışları yazmanıza yardımcı olur
- **Dayanıklı işleme**: Güvenli gezinme ve Elvis operatörleri, iş akışlarını eksik verilere karşı dayanıklı kılar
- **Esnek işleme**: Koşullu mantık, iş akışlarınızın farklı örnek türlerini uygun şekilde işlemesini sağlar
- **Uyarlanabilir kaynaklar**: Dinamik yönergeler, girdi özelliklerine göre kaynak kullanımını optimize eder

Bu ilerleme, biyoinformatik pipeline'larının gerçek dünya evrimini yansıtmaktadır; birkaç örneği işleyen araştırma prototiplerinden laboratuvarlar ve kurumlar genelinde binlerce örneği işleyen üretim sistemlerine. Çözdüğünüz her zorluk ve öğrendiğiniz her kalıp, geliştiricilerin Nextflow iş akışlarını ölçeklendirirken karşılaştığı gerçek sorunları yansıtmaktadır.

Bu kalıpları kendi çalışmalarınızda uygulamak, dayanıklı, üretime hazır iş akışları oluşturmanızı sağlayacaktır.

### Temel kalıplar

1.  **Veri Akışı ve Betik Yazımı:** Veri akışı işlemleri (kanal orkestrasyonu) ile betik yazımı (veriyi işleyen kod) arasındaki farkı öğrendiniz; Channel ve List üzerindeki `collect` gibi farklı türlerdeki işlemler arasındaki kritik farklar dahil.

    - Veri akışı: kanal orkestrasyonu

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Betik yazımı: koleksiyonlar üzerinde veri işleme

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Gelişmiş Dize İşleme**: Dosya adlarını ayrıştırmak için düzenli ifadelerde, süreçlerde dinamik betik oluşturmada ve değişken enterpolasyonunda (Nextflow, Bash ve Kabuk) ustalaştınız.

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

    - Dosya koleksiyonundan komut argümanlarına (süreç betik bloğunda)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Yeniden Kullanılabilir Fonksiyonlar Oluşturmak**: Karmaşık mantığı, kanal operatörlerinden çağrılabilen adlandırılmış fonksiyonlara çıkarmayı öğrendiniz; iş akışlarını daha okunabilir ve sürdürülebilir hale getirdiniz.

    - Adlandırılmış fonksiyon tanımlamak

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

    - Adlandırılmış fonksiyonu iş akışında çağırmak

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Closure'larla Dinamik Kaynak Yönergeleri**: Girdi özelliklerine göre uyarlanabilir kaynak tahsisi için süreç yönergelerinde closure kullanmayı keşfettiniz.

    - Adlandırılmış closure'lar ve bileşim

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Kapsam erişimli closure'lar

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Koşullu Mantık ve Süreç Kontrolü**: `.branch()` ve `.filter()` operatörlerini kullanarak akıllı yönlendirme eklediniz; kısa koşullu ifadeler için doğruluk değerinden yararlandınız.

    - Veriyi farklı iş akışı dallarına yönlendirmek için `.branch()` kullanmak

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Groovy doğruluk değeriyle boolean değerlendirme

    ```groovy
    if (sample.files) println "Dosyalar mevcut"
    ```

    - 'Doğruluk değeri' ile veriyi alt kümelere ayırmak için `filter()` kullanmak

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

7.  **`error()` ve `log.warn` ile Doğrulama**: Girdileri erken doğrulamayı ve açık hata mesajlarıyla hızlı başarısız olmayı öğrendiniz.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Geçersiz: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Hata: ${e.message}"
    }
    ```

8.  **Yapılandırma Olay İşleyicileri**: Günlükleme, bildirimler ve yaşam döngüsü yönetimi için iş akışı olay işleyicilerini (`onComplete` ve `onError`) kullanmayı öğrendiniz.

    - Günlükleme ve bildirim için `onComplete` kullanmak

    ```groovy
    workflow.onComplete = {
        println "Başarı     : ${workflow.success}"
        println "Çıkış kodu : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline başarıyla tamamlandı!"
        } else {
            println "❌ Pipeline başarısız oldu!"
            println "Hata: ${workflow.errorMessage}"
        }
    }
    ```

    - Özellikle başarısızlık durumunda harekete geçmek için `onError` kullanmak

    ```groovy
    workflow.onError = {
        // Ayrıntılı hata günlüğü yaz
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

Daha gelişmiş özellikleri keşfetmeniz gerektiğinde bu kaynaklara başvurmayı unutmayın.

Becerilerinizi geliştirmek ve genişletmek için pratik yapmanız şu konularda fayda sağlayacaktır:

- Veri akışı ve betik yazımı arasında uygun ayrımla daha temiz iş akışları yazmak
- Nextflow, Bash ve kabuk değişkenlerindeki yaygın tuzaklardan kaçınmak için değişken enterpolasyonunda ustalaşmak
- Verimli, uyarlanabilir iş akışları için dinamik kaynak yönergeleri kullanmak
- Dosya koleksiyonlarını doğru biçimlendirilmiş komut satırı argümanlarına dönüştürmek
- Regex ve dize işleme kullanarak farklı dosya adlandırma kurallarını ve girdi formatlarını zarif biçimde işlemek
- Gelişmiş closure kalıpları ve fonksiyonel programlama kullanarak yeniden kullanılabilir, sürdürülebilir kod oluşturmak
- Koleksiyon işlemlerini kullanarak karmaşık veri kümelerini işlemek ve düzenlemek
- İş akışlarınızı üretime hazır hale getirmek için doğrulama, hata işleme ve günlükleme eklemek
- Olay işleyicileriyle iş akışı yaşam döngüsü yönetimi uygulamak

---

## Sırada ne var?

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
