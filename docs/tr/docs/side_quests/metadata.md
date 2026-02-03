# Metadata ve meta map'ler

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herhangi bir bilimsel analizde, sadece ham veri dosyalarıyla çalışmayız.
Her dosya kendi ek bilgileriyle birlikte gelir: ne olduğu, nereden geldiği ve onu özel kılan nedir.
Bu ek bilgilere metadata (üst veri) diyoruz.

Metadata, diğer verileri tanımlayan veridir.
Metadata, dosyalar ve deneysel koşullar hakkındaki önemli ayrıntıları takip eder ve analizlerin her veri setinin benzersiz özelliklerine göre uyarlanmasına yardımcı olur.

Bunu bir kütüphane kataloğu gibi düşünün: kitaplar gerçek içeriği (ham veri) içerirken, katalog kartları her kitap hakkında temel bilgiler sağlar—ne zaman yayınlandığı, kimin yazdığı, nerede bulunacağı (metadata).
Nextflow pipeline'larında metadata şunlar için kullanılabilir:

- Dosyaya özgü bilgileri workflow boyunca takip etmek
- Dosya özelliklerine göre process'leri yapılandırmak
- İlgili dosyaları ortak analiz için gruplamak

### Öğrenme hedefleri

Bu yan görevde, workflow'larda metadata'nın nasıl ele alınacağını keşfedeceğiz.
Temel dosya bilgilerini içeren basit bir veri tablosundan (biyoinformatikte genellikle samplesheet olarak adlandırılır) başlayarak, şunları öğreneceksiniz:

- CSV dosyalarından dosya metadata'sını okuma ve ayrıştırma
- Metadata map'lerini oluşturma ve manipüle etme
- Workflow yürütme sırasında yeni metadata alanları ekleme
- Process davranışını özelleştirmek için metadata kullanma

Bu beceriler, karmaşık dosya ilişkilerini ve işleme gereksinimlerini işleyebilen daha sağlam ve esnek pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramlarını ve mekanizmalarını (process'ler, channel'lar, operatörler) rahatça kullanabiliyor olmalısınız

---

## 0. Başlangıç

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

Bu eğitimin dosyalarının bulunduğu dizine geçelim.

```bash
cd side-quests/metadata
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

Bir ana workflow dosyası ve bir veri tablosu ile birkaç veri dosyası içeren bir `data` dizini bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

`main.nf` dosyasındaki workflow, kademeli olarak tam işlevli bir workflow'a genişleteceğiniz bir taslaktır.

Veri tablosu, veri dosyalarına giden yolları ve ilişkili bazı metadata'yı 3 sütunda düzenlenmiş şekilde listeler:

- `id`: açıklayıcı, dosyaya verilen bir ID
- `character`: bir karakter adı, daha sonra farklı yaratıklar çizmek için kullanacağız
- `data`: farklı dillerde selamlaşmalar içeren `.txt` dosyalarının yolları

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Her veri dosyası, beş dilden birinde (fr: Fransızca, de: Almanca, es: İspanyolca, it: İtalyanca, en: İngilizce) selamlaşma metni içerir.

Ayrıca size `langid` adında konteynerleştirilmiş bir dil analiz aracı sağlayacağız.

#### Görevi inceleyin

Meydan okumanız, şunları yapacak bir Nextflow workflow'u yazmaktır:

1. Her dosyadaki dili otomatik olarak **Tanımlama**
2. Dosyaları dil ailesine göre **Gruplama** (Cermen dilleri vs Roman dilleri)
3. Her dosya için dilini ve metadata'sına göre işlemeyi **Özelleştirme**
4. Çıktıları dil grubuna göre **Organize etme**

Bu, dosyaya özgü metadata'nın işleme kararlarını yönlendirdiği tipik bir workflow modelini temsil eder; tam olarak metadata map'lerinin zarif bir şekilde çözdüğü sorun türü.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Bir veri tablosundan metadata yükleme

Başlangıç noktası olarak size verdiğimiz workflow taslağını incelemek için `main.nf` workflow dosyasını açın.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Örnek veri tablosunu bir dosya olarak yüklemek için temel bir channel factory kurduğumuzu görebilirsiniz, ancak bu henüz dosyanın içeriğini okumayacaktır.
Bunu ekleyerek başlayalım.

### 1.1. `splitCsv` ile içeriği okuma

Dosya içeriğini minimal çabayla uygun şekilde ayrıştıracak bir operatör seçmemiz gerekiyor.
Veri tablomuz CSV formatında olduğundan, bu iş [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörü içindir, bu operatör dosyadaki her satırı channel'da bir öğe olarak yükler.

Channel oluşturma koduna bir `splitCsv()` işlemi eklemek için aşağıdaki değişiklikleri yapın, artı dosyanın içeriğinin channel'a doğru şekilde yüklendiğini kontrol etmek için bir `view()` işlemi.

=== "Sonra"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Nextflow'a CSV dosyasının ilk satırını başlık satırı olarak okumasını söylemek için `header: true` seçeneğini kullandığımıza dikkat edin.

Bundan ne çıktığını görelim mi?
Workflow'u çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Operatörün CSV dosyasındaki her satır için, sütun başlıklarını karşılık gelen değerlerin anahtarları olarak kullanarak bir anahtar-değer çifti map'i oluşturduğunu görebiliriz.

Her map girişi veri tablomuzun bir sütununa karşılık gelir:

- `id`
- `character`
- `recording`

Bu harika! Belirli alanlara her dosyadan erişmeyi kolaylaştırıyor.
Örneğin, dosya ID'sine `id` ile veya txt dosya yoluna `recording` ile erişebiliriz.

??? info "(İsteğe bağlı) Map'ler hakkında daha fazla bilgi"

    Nextflow'un üzerine kurulu olduğu programlama dili olan Groovy'de, map, Python'daki dictionary'lere, JavaScript'teki object'lere veya Ruby'deki hash'lere benzer bir anahtar-değer veri yapısıdır.

    İşte pratikte bir map'i nasıl tanımlayabileceğinizi ve içeriğine nasıl erişebileceğinizi gösteren çalıştırılabilir bir script:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Basit bir map oluştur
    def my_map = [id:'sampleA', character:'squirrel']

    // Tüm map'i yazdır
    println "map: ${my_map}"

    // Nokta notasyonu kullanarak bireysel değerlere eriş
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Düzgün bir `workflow` bloğu olmasa bile, Nextflow bunu bir workflow gibi çalıştırabilir:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Ve çıktıda görmeyi bekleyebileceğiniz şey:

    ```console title="Çıktı"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` ile belirli alanları seçme

Diyelim ki veri tablosundan `character` sütununa erişmek ve yazdırmak istiyoruz.
Channel'ımızdaki her öğe üzerinde yineleme yapmak ve map nesnesinden özellikle `character` girişini seçmek için Nextflow `map` operatörünü kullanabiliriz.

Workflow'da aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Şimdi workflow'u tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Başarılı! Veri tablomuzdan türetilen map yapısından yararlanarak her satır için bireysel sütunlardan değerlere eriştik.

Şimdi veri tablosunu başarıyla okuduk ve her satırdaki verilere erişimimiz var, pipeline mantığımızı uygulamaya başlayabiliriz.

### 1.3. Metadata'yı bir 'meta map' olarak organize etme

Workflow'un mevcut durumunda, girdi dosyaları (`recording` anahtarı altında) ve ilişkili metadata (`id`, `character`) hepsi aynı düzeydedir, sanki hepsi aynı büyük çantadalar.
Pratik sonuç, bu channel'ı tüketen her process'in bu yapıyı göz önünde bulundurarak yapılandırılması gerektiğidir:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Veri tablosundaki sütun sayısı değişmediği sürece bu iyidir.
Ancak, veri tablosuna sadece bir sütun daha eklerseniz, channel'ın şekli artık process'in beklediğiyle eşleşmeyecek ve workflow hata üretecektir.
Ayrıca process'i, script bloğu tarafından gerekli olmayan değişkenleri sabit kodlamak zorunda kalabileceğiniz için, biraz farklı girdi verilerine sahip olabilecek başkalarıyla paylaşmayı zorlaştırır.

Bu sorunu önlemek için, veri tablosunun kaç sütun içerdiğinden bağımsız olarak channel yapısını tutarlı tutmanın bir yolunu bulmamız gerekiyor.

Bunu, tüm metadata'yı tuple içindeki bir öğede toplayarak yapabiliriz, buna metadata map veya daha basit ifadeyle 'meta map' diyeceğiz.

`map` işleminde aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Channel öğelerimizi, meta map ve karşılık gelen dosya nesnesi olmak üzere iki öğeden oluşan bir tuple olarak yeniden yapılandırdık.

Workflow'u çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console title="Meta map'i görüntüle"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Şimdi, channel'daki her öğe önce metadata map'i, sonra da karşılık gelen dosya nesnesini içerir:

```console title="Örnek çıktı yapısı"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Sonuç olarak, veri tablosuna daha fazla sütun eklemek `meta` map'inde daha fazla metadata kullanılabilir hale getirecek, ancak channel şeklini değiştirmeyecektir.
Bu, metadata öğelerini girdi belirtimine sabit kodlamak zorunda kalmadan channel'ı tüketen process'ler yazmamızı sağlar:

```groovy title="Sözdizimi örneği"
    input:
    tuple val(meta), file(recording)
```

Bu, Nextflow workflow'larında metadata'yı organize etmek için yaygın olarak kullanılan bir kuraldır.

### Çıkarımlar

Bu bölümde, şunları öğrendiniz:

- **Metadata neden önemlidir:** Metadata'yı verinizle birlikte tutmak, workflow boyunca önemli dosya bilgilerini korur.
- **Veri tablolarını nasıl okuyacağınız:** Başlık bilgisine sahip CSV dosyalarını okumak ve satırları yapılandırılmış verilere dönüştürmek için `splitCsv` kullanımı
- **Meta map nasıl oluşturulur:** Tuple yapısı `[ [id:value, ...], file ]` kullanarak metadata'yı dosya verisinden ayırma

---

## 2. Metadata'yı manipüle etme

Artık metadata'mız yüklendiğine göre, onunla bir şeyler yapalım!

Her yaratığın kayıt dosyasında bulunan dili tanımlamak için [`langid`](https://github.com/saffsd/langid.py) adlı bir araç kullanacağız.
Araç bir dizi dil üzerinde önceden eğitilmiş olarak gelir ve bir metin parçası verildiğinde, `stdout`'a bir dil tahmini ve ilişkili bir olasılık puanı çıktısı verir.

### 2.1. Process'i içe aktarın ve kodu inceleyin

Size `langid` aracını saran `IDENTIFY_LANGUAGE` adlı önceden yazılmış bir process modülü sağlıyoruz, bu nedenle sadece workflow bloğundan önce bir include ifadesi eklemeniz gerekiyor.

Workflow'da aşağıdaki düzenlemeyi yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Kodunu incelemek için modül dosyasını açabilirsiniz:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Her girdi dosyasının dilini tahmin etmek için langid kullan
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Gördüğünüz gibi, girdi tanımı girdi channel'ımıza uyguladığımız aynı `tuple val(meta), path(file)` yapısını kullanıyor.

Çıktı tanımı, girdi yapısına benzer bir yapıya sahip bir tuple olarak yapılandırılmıştır, ancak üçüncü öğe olarak `stdout` da içerir.
Bu `tuple val(meta), path(file), <output>` modeli, metadata'yı hem girdi verisiyle hem de çıktılarla ilişkili tutar ve pipeline boyunca akar.

Burada Nextflow'un [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) çıktı niteleyicisini kullandığımızı unutmayın çünkü araç çıktısını dosya yazmak yerine doğrudan konsola yazdırır; ve komut satırında olasılık puanını kaldırmak, yeni satır karakterlerini kaldırarak dizeyi temizlemek ve yalnızca dil tahminini döndürmek için `sed` kullanırız.

### 2.2. `IDENTIFY_LANGUAGE`'e bir çağrı ekleyin

Artık process workflow için kullanılabilir olduğuna göre, veri channel'ında çalıştırmak için `IDENTIFY_LANGUAGE` process'ine bir çağrı ekleyebiliriz.

Workflow'da aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Orijinal `.view()` işlemini channel oluşturmada kaldırdığımızı unutmayın.

Artık workflow'u çalıştırabiliriz:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Mükemmel! Şimdi her karakterin hangi dilde konuştuğuna dair bir tahminimiz var.

Ve daha önce belirtildiği gibi, çıktıya girdi dosyasını ve meta map'i de dahil ettik, bu da her ikisinin de az önce ürettiğimiz yeni bilgiyle ilişkili kalması anlamına gelir.
Bu, bir sonraki adımda faydalı olacaktır.

!!! note

    Daha genel olarak, meta map'i sonuçlarla ilişkili tutmanın bu modeli, aynı tanımlayıcıları paylaşan ilgili sonuçları ilişkilendirmeyi kolaylaştırır.

    Zaten öğrenmiş olduğunuz gibi, sonuçları bunlar arasında eşleştirmek için channel'lardaki öğelerin sırasına güvenemezsiniz.
    Bunun yerine, verileri doğru şekilde ilişkilendirmek için anahtarlar kullanmalısınız ve meta map'ler bu amaç için ideal bir yapı sağlar.

    Bu kullanım senaryosunu [Splitting & Grouping](./splitting_and_grouping.md) yan görevinde ayrıntılı olarak inceliyoruz.

### 2.3. Process çıktılarıyla metadata'yı genişletme

Az önce ürettiğimiz sonuçların kendilerinin dosyaların içeriği hakkında bir metadata türü olduğu göz önüne alındığında, bunları meta map'imize eklemek faydalı olacaktır.

Ancak, mevcut meta map'i yerinde değiştirmek istemiyoruz.
Teknik bir bakış açısından, bunu yapmak _mümkündür_, ancak güvenli değildir.

Bu nedenle, bunun yerine `+` operatörünü (bir Groovy özelliği) kullanarak mevcut meta map'in içeriğini artı yeni bilgileri tutan yeni bir `lang: lang_id` anahtar-değer çiftini içeren yeni bir meta map oluşturacağız.
Ve bunu eski map'i yenisiyle değiştirmek için bir [`map`](https://www.nextflow.io/docs/latest/operator.html#map) işlemiyle birleştireceğiz.

Workflow'da yapmanız gereken düzenlemeler şunlardır:

=== "Sonra"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

`+` operatörüne henüz aşina değilseniz veya bu kafa karıştırıcı geliyorsa, aşağıdaki ayrıntılı açıklamayı gözden geçirmek için birkaç dakika ayırın.

??? info "`+` operatörü kullanarak yeni meta map'in oluşturulması"

    **İlk olarak, Groovy `+` operatörünü kullanarak iki map'in içeriğini birleştirebileceğimizi bilmeniz gerekir.**

    Diyelim ki aşağıdaki map'lere sahibiz:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Bunları şu şekilde birleştirebiliriz:

    ```groovy
    new_map = map1 + map2
    ```

    `new_map`'in içeriği şu şekilde olacaktır:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Harika!

    **Peki zaten bir map'in parçası olmayan bir alan eklemeniz gerekirse ne olur?**

    Diyelim ki `map1`'den tekrar başlıyorsunuz, ancak dil tahmini kendi map'inde değil (bir `map2` yok).
    Bunun yerine, `lang_id` adlı bir değişkende tutuluyor ve değerini (`'fr'`) `lang` anahtarıyla saklamak istediğinizi biliyorsunuz.

    Aslında şunu yapabilirsiniz:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Burada, `[lang: new_info]` anında yeni bir adsız map oluşturur ve `map1 + ` `map1`'i yeni adsız map ile birleştirerek daha önce olduğu gibi aynı `new_map` içeriğini üretir.

    Düzgün, değil mi?

    **Şimdi bunu bir Nextflow `channel.map()` işlemi bağlamına aktaralım.**

    Kod şöyle olur:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Bu şunları yapar:

    - `map1, lang_id ->` tuple'daki iki öğeyi alır
    - `[map1 + [lang: lang_id]]` yukarıda detaylandırıldığı gibi yeni map'i oluşturur

    Çıktı, yukarıdaki örneğimizdeki `new_map` ile aynı içeriğe sahip tek bir adsız map'tir.
    Yani etkin bir şekilde şunu dönüştürdük:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    şuna:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Umarım `map1`'i `meta` olarak değiştirirsek, bunun temelde workflow'umuzda meta map'imize dil tahminini eklemek için ihtiyacımız olan her şey olduğunu görebilirsiniz.

    Bir şey dışında!

    Workflow'umuz söz konusu olduğunda, **`meta, file, lang_id`'den oluşan tuple'da `file` nesnesinin varlığını da hesaba katmamız gerekir**.

    Yani buradaki kod şöyle olur:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    `file`'ın `map` işleminde neden dolaştığını anlamakta zorlanıyorsanız, `[meta + [lang: lang_id], file]` yerine o satırın `[new_map, file]` olduğunu hayal edin.
    Bu, basitçe `file`'ı tuple'daki ikinci pozisyondaki orijinal yerinde bıraktığımızı daha açık hale getirmelidir. Sadece `new_info` değerini aldık ve onu birinci pozisyondaki map'e katladık.

    **Ve bu bizi `tuple val(meta), path(file)` channel yapısına geri getiriyor!**

Bu kodun ne yaptığını anladığınızdan emin olduğunuzda, çalışıp çalışmadığını görmek için workflow'u çalıştırın:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Evet, bu kontrol edildi!
Process'in çıktısını `meta, file, lang_id`'den düzgün bir şekilde yeniden düzenledik, böylece `lang_id` artık meta map'teki anahtarlardan biri ve channel'ın tuple'ları tekrar `meta, file` modeline uyuyor.

<!-- TODO (future) subMap kullanarak bir anahtarın nasıl kaldırılacağını da göstermeli miyiz?! Veya bunu nerede bulacağımızı not edelim. -->

### 2.4. Koşullu ifadeler kullanarak bir dil grubu atama

Artık dil tahminlerimiz olduğuna göre, yeni gruplamalar atamak için bilgileri kullanalım.

Örnek verilerimizde, karakterlerimiz tarafından kullanılan diller cermen dilleri (İngilizce, Almanca) ve roman dilleri (Fransızca, İspanyolca, İtalyanca) olarak gruplandırılabilir.
Pipeline'da daha sonra bu sınıflandırmanın hazır bir şekilde mevcut olması faydalı olabilir, bu yüzden o bilgiyi meta map'e ekleyelim.

Ve, iyi haber, bu map operatörünü kullanmaya mükemmel bir şekilde uyan bir başka durumdur!

Özellikle, `lang_group` adında bir değişken tanımlayacağız, her veri parçası için `lang_group`'a hangi değerin atanacağını belirlemek için bazı basit koşullu mantık kullanacağız.

Genel sözdizimi şöyle görünecek:

```groovy
.map { meta, file ->

    // lang_group'u tanımlayan koşullu mantık buraya gelir

    [meta + [lang_group: lang_group], file]
}
```

Gördüğünüz gibi bu, önceki adımda kullandığımız anında map birleştirme işlemine çok benziyor.
Sadece koşullu ifadeleri yazmamız gerekiyor.

İşte uygulamak istediğimiz koşullu mantık:

- Varsayılan değeri `'unknown'` olan `lang_group` adında bir değişken tanımlayın.
- Eğer `lang` Almanca (`'de'`) veya İngilizce (`'en'`) ise, `lang_group`'u `germanic` olarak değiştirin.
- Değilse eğer `lang` Fransızca (`'fr'`), İspanyolca (`'es'`) ve İtalyanca (`'it'`) içeren bir listede yer alıyorsa, `lang_group`'u `romance` olarak değiştirin.

Nextflow'da koşullu ifadeleri nasıl yazacağınızı zaten biliyorsanız kendiniz yazmayı deneyin.

!!! tip

    Map işlemi içinde `lang` değerine `meta.lang` ile erişebilirsiniz.

Workflow'da aşağıdaki değişiklikleri yapmayı tamamlamalısınız:

=== "Sonra"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

İşte ana noktalar:

- Varsayılan değeri `unknown` olarak ayarlanmış `lang_group` değişkenini oluşturmak için `def lang_group = "unknown"` kullanırız.
- İki cermen dili için alternatif `.equals()` testleri ve üç roman dili için bir listede varlık testi ile koşullu mantık için bir `if {} else if {}` yapısı kullanırız.
- Güncellenmiş meta map'i oluşturmak için daha önce olduğu gibi `meta + [lang_group:lang_group]` birleştirme işlemini kullanırız.

<!-- TODO (future) Ek kaynaklar bölümünde ilgili belgelere not/bağlantılar ekle -->

Hepsi mantıklı olduğunda, sonucu görmek için workflow'u tekrar çalıştırın:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Gördüğünüz gibi, channel öğeleri `[meta, file]` yapısını korur, ancak meta map artık bu yeni sınıflandırmayı içerir.

### Çıkarımlar

Bu bölümde, nasıl yapılacağını öğrendiniz:

- **Girdi metadata'sını çıktı channel'larına uygulama**: Metadata'yı bu şekilde kopyalamak, sonuçları daha sonra metadata içeriğine göre ilişkilendirmemize olanak tanır.
- **Özel anahtarlar oluşturma**: Meta map'inizde iki yeni anahtar oluşturdunuz, bunları `meta + [new_key:value]` ile mevcut meta map'e birleştirerek. Biri bir process'ten hesaplanan bir değere, diğeri `map` operatöründe belirlediğiniz bir koşula dayalı olarak.

Bunlar, pipeline boyunca ilerlerken yeni ve mevcut metadata'yı dosyalarla ilişkilendirmenize olanak tanır.
Bir process'in parçası olarak metadata kullanmasanız bile, meta map'i bu şekilde veriyle ilişkili tutmak tüm ilgili bilgileri bir arada tutmayı kolaylaştırır.

---

## 3. Bir process'te meta map bilgisini kullanma

Artık meta map'i nasıl oluşturup güncelleyeceğinizi bildiğinize göre, gerçekten eğlenceli kısma geçebiliriz: metadata'yı bir process'te gerçekten kullanma.

Daha spesifik olarak, workflow'umuza her hayvanı ASCII art olarak çizmek ve kaydedilen metni bir konuşma balonunda söyletmek için ikinci bir adım ekleyeceğiz.
Bunu [`cowpy`](https://github.com/jeffbuttars/cowpy) adlı bir araç kullanarak yapacağız.

??? info "`cowpy` ne yapar?"

    `cowpy`, keyfi metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII art üreten bir komut satırı aracıdır.
    Tony Monroe'nun klasik [cowsay](https://en.wikipedia.org/wiki/Cowsay) aracının python implementasyonudur.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    İsteğe bağlı olarak, varsayılan inek yerine kullanılacak bir karakter (veya 'cowacter') seçebilirsiniz.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Hello Nextflow kursunu tamamladıysanız, bu aracı zaten çalışırken gördünüz.
Eğer görmediyseniz, endişelenmeyin; ilerlerken bilmeniz gereken her şeyi ele alacağız.

### 3.1. Process'i içe aktarın ve kodu inceleyin

Size `cowpy` aracını saran `COWPY` adlı önceden yazılmış bir process modülü sağlıyoruz, bu nedenle sadece workflow bloğundan önce bir include ifadesi eklemeniz gerekiyor.

Workflow'da aşağıdaki düzenlemeyi yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Kodunu incelemek için modül dosyasını açabilirsiniz:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy ile ASCII art oluştur
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Gördüğünüz gibi, bu process şu anda bir girdi dosyası (görüntülenecek metni içeren) ve ASCII art'ta çizilmesi gereken karakteri belirten bir değer almak üzere tasarlanmıştır, genellikle workflow seviyesinde bir komut satırı parametresi ile sağlanır.

### 3.2. Bir meta map alanını girdi olarak geçirme

Hello Nextflow kursunda `cowpy` aracını kullandığımızda, nihai görüntüyü çizmek için hangi karakterin kullanılacağını belirlemek için bir komut satırı parametresi kullandık.
Bu mantıklıydı, çünkü pipeline'ın her çalıştırması başına sadece bir görüntü oluşturuyorduk.

Ancak, bu eğitimde, işlediğimiz her konu için uygun bir görüntü oluşturmak istiyoruz, bu nedenle bir komut satırı parametresi kullanmak çok sınırlayıcı olurdu.

İyi haber: veri tablomuzda ve dolayısıyla meta map'imizde bir `character` sütunumuz var.
Process'in her girdi için kullanması gereken karakteri ayarlamak için bunu kullanalım.

Bu amaçla, üç şey yapmamız gerekecek:

1. Üzerinde daha rahat çalışabilmemiz için önceki process'ten çıkan çıktı channel'ına bir isim verin
2. İlgilendiğimiz bilgilere nasıl erişeceğimizi belirleyin
3. İkinci process'e bir çağrı ekleyin ve bilgiyi uygun şekilde besleyin

Başlayalım.

#### 3.2.1. Önceki çıktı channel'ını adlandırma

Önceki manipülasyonları doğrudan ilk process'in çıktı channel'ı üzerinde uyguladık, `IDENTIFY_LANGUAGE.out`.
Bir sonraki process'e channel'ın içeriğini beslemek için (ve bunu açık ve kolay okunur bir şekilde yapmak için) ona kendi adını vermek istiyoruz, `ch_languages`.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana workflow'da, `.view()` operatörünü `.set { ch_languages }` ile değiştirin ve channel'a ismiyle başvurabileceğimizi test eden bir satır ekleyin.

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // Geçici: ch_languages'a göz at
        ch_languages.view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Her selamlaşmanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

Bunu çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Bu, artık channel'a ismiyle başvurabileceğimizi doğruluyor.

#### 3.2.2. Dosya ve karakter metadata'sına erişim

Modül koduna bakarak `COWPY` process'inin bir metin dosyası ve bir `character` değeri verilmesini beklediğini biliyoruz.
`COWPY` process çağrısını yazmak için, channel'daki her öğeden karşılık gelen dosya nesnesini ve metadata'yı nasıl çıkaracağımızı bilmemiz gerekiyor.

Genellikle olduğu gibi, bunu yapmanın en basit yolu bir `map` işlemi kullanmaktır.

Channel'ımız `[meta, file]` olarak yapılandırılmış tuple'lar içerir, bu nedenle `file` nesnesine doğrudan erişebiliriz ve meta map içinde saklanan `character` değerine `meta.character` olarak başvurarak erişebiliriz.

Ana workflow'da, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34"
        // Geçici: dosya ve karaktere eriş
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="34"
        // Geçici: ch_languages'a göz at
        ch_languages.view()
    ```

`.view` işlemlerinin çıktısını daha okunabilir hale getirmek için closure'ları (örneğin `{ file -> "File: " + file }`) kullandığımızı unutmayın.

Bunu çalıştıralım:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Dosya yolları ve karakter değerleri çıktınızda farklı bir sırayla gelebilir._

Bu, channel'daki her öğe için dosyaya ve karaktere erişebildiğimizi doğruluyor.

#### 3.2.3. `COWPY` process'ini çağırma

Şimdi hepsini bir araya getirelim ve gerçekten `ch_languages` channel'ında `COWPY` process'ini çağıralım.

Ana workflow'da, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34"
        // ASCII art oluşturmak için cowpy çalıştır
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="34"
        // Geçici: dosya ve karaktere eriş
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Sadece iki
