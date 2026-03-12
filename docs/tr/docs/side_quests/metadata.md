# Metadata ve meta map'ler

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herhangi bir bilimsel analizde, sadece ham veri dosyalarıyla çalışmayız.
Her dosya kendi ek bilgileriyle birlikte gelir: ne olduğu, nereden geldiği ve onu özel kılan nedir.
Bu ek bilgilere metadata (üst veri) diyoruz.

Metadata, diğer verileri tanımlayan veridir.
Metadata, dosyalar ve deneysel koşullar hakkındaki önemli ayrıntıları takip eder ve analizlerin her veri setinin benzersiz özelliklerine göre uyarlanmasına yardımcı olur.

Bunu bir kütüphane kataloğu gibi düşünün: kitaplar gerçek içeriği (ham veri) içerirken, katalog kartları her kitap hakkında temel bilgiler sağlar—ne zaman yayınlandığı, kimin yazdığı, nerede bulunacağı (metadata).
Nextflow pipeline'larında metadata şunlar için kullanılabilir:

- Dosyaya özgü bilgileri iş akışı boyunca takip etmek
- Dosya özelliklerine göre süreçleri yapılandırmak
- İlgili dosyaları ortak analiz için gruplamak

### Öğrenme hedefleri

Bu yan görevde, iş akışlarında metadata'nın nasıl ele alınacağını keşfedeceğiz.
Temel dosya bilgilerini içeren basit bir veri tablosundan (biyoinformatikte genellikle samplesheet olarak adlandırılır) başlayarak şunları öğreneceksiniz:

- CSV dosyalarından dosya metadata'sını okuma ve ayrıştırma
- Metadata map'lerini oluşturma ve manipüle etme
- İş akışı yürütme sırasında yeni metadata alanları ekleme
- Süreç davranışını özelleştirmek için metadata kullanma

Bu beceriler, karmaşık dosya ilişkilerini ve işleme gereksinimlerini işleyebilen daha sağlam ve esnek pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramlarını ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabiliyor olmalısınız.

---

## 0. Başlarken

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

Bir ana iş akışı dosyası ve bir veri tablosu ile birkaç veri dosyası içeren bir `data` dizini bulacaksınız.

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

`main.nf` dosyasındaki iş akışı, kademeli olarak tam işlevli bir iş akışına genişleteceğiniz bir taslaktır.

Veri tablosu, veri dosyalarına giden yolları ve ilişkili bazı metadata'yı 3 sütunda düzenlenmiş şekilde listeler:

- `id`: açıklayıcı, dosyaya verilen bir kimlik
- `character`: bir karakter adı; daha sonra farklı yaratıklar çizmek için kullanacağız
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

Meydan okumanız, şunları yapacak bir Nextflow iş akışı yazmaktır:

1. Her dosyadaki dili otomatik olarak **Tanımlama**
2. Dosyaları dil ailesine göre **Gruplama** (Cermen dilleri vs Roman dilleri)
3. Her dosya için dilini ve metadata'sına göre işlemeyi **Özelleştirme**
4. Çıktıları dil grubuna göre **Organize etme**

Bu, dosyaya özgü metadata'nın işleme kararlarını yönlendirdiği tipik bir iş akışı modelini temsil eder; tam olarak metadata map'lerinin zarif bir şekilde çözdüğü sorun türü.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışır durumda
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Bir veri tablosundan metadata yükleme

Başlangıç noktası olarak size verdiğimiz iş akışı taslağını incelemek için `main.nf` dosyasını açın.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Örnek veri tablosunu bir dosya olarak yüklemek için temel bir channel factory kurduğumuzu görebilirsiniz; ancak bu henüz dosyanın içeriğini okumayacaktır.
Bunu ekleyerek başlayalım.

### 1.1. `splitCsv` ile içeriği okuma

Dosya içeriğini minimal çabayla uygun şekilde ayrıştıracak bir operatör seçmemiz gerekiyor.
Veri tablomuz CSV formatında olduğundan, bu iş [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörü içindir; bu operatör dosyadaki her satırı kanalda bir öğe olarak yükler.

Kanal oluşturma koduna bir `splitCsv()` işlemi eklemek için aşağıdaki değişiklikleri yapın. Ayrıca dosyanın içeriğinin kanala doğru şekilde yüklendiğini kontrol etmek için bir `view()` işlemi de ekleyin.

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
İş akışını çalıştırın:

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
Örneğin, dosya kimliğine `id` ile veya txt dosya yoluna `recording` ile erişebiliriz.

??? info "(İsteğe bağlı) Map'ler hakkında daha fazla bilgi"

    Nextflow'un üzerine kurulu olduğu programlama dili olan Groovy'de, map; Python'daki dictionary'lere, JavaScript'teki object'lere veya Ruby'deki hash'lere benzer bir anahtar-değer veri yapısıdır.

    İşte pratikte bir map'i nasıl tanımlayabileceğinizi ve içeriğine nasıl erişebileceğinizi gösteren çalıştırılabilir bir betik:

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

    Düzgün bir `workflow` bloğu olmasa bile, Nextflow bunu bir iş akışı gibi çalıştırabilir:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Çıktıda görmeyi bekleyebileceğiniz şey:

    ```console title="Çıktı"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` ile belirli alanları seçme

Diyelim ki veri tablosundan `character` sütununa erişmek ve yazdırmak istiyoruz.
Kanalımızdaki her öğe üzerinde yineleme yapmak ve map nesnesinden özellikle `character` girişini seçmek için Nextflow `map` operatörünü kullanabiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

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

Şimdi iş akışını tekrar çalıştırın:

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

Veri tablosunu başarıyla okuduk ve her satırdaki verilere erişimimiz var; artık pipeline mantığımızı uygulamaya başlayabiliriz.

### 1.3. Metadata'yı bir 'meta map' olarak organize etme

İş akışının mevcut durumunda, girdi dosyaları (`recording` anahtarı altında) ve ilişkili metadata (`id`, `character`) hepsi aynı düzeydedir; sanki hepsi aynı büyük çantadadır.
Bunun pratik sonucu, bu kanalı tüketen her sürecin bu yapıyı göz önünde bulundurarak yapılandırılması gerektiğidir:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Veri tablosundaki sütun sayısı değişmediği sürece bu iyidir.
Ancak veri tablosuna sadece bir sütun daha eklerseniz, kanalın şekli artık sürecin beklediğiyle eşleşmeyecek ve iş akışı hata üretecektir.
Ayrıca süreç, betik bloğu tarafından gerekli olmayan değişkenleri sabit kodlamak zorunda kalabileceğiniz için biraz farklı girdi verilerine sahip olabilecek başkalarıyla paylaşmayı zorlaştırır.

Bu sorunu önlemek için, veri tablosunun kaç sütun içerdiğinden bağımsız olarak kanal yapısını tutarlı tutmanın bir yolunu bulmamız gerekiyor.

Bunu, tüm metadata'yı tuple içindeki bir öğede toplayarak yapabiliriz; buna metadata map veya daha basit ifadeyle 'meta map' diyeceğiz.

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

Kanal öğelerimizi, meta map ve karşılık gelen dosya nesnesi olmak üzere iki öğeden oluşan bir tuple olarak yeniden yapılandırdık.

İş akışını çalıştıralım:

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

Artık kanaldaki her öğe önce metadata map'i, sonra da karşılık gelen dosya nesnesini içerir:

```console title="Örnek çıktı yapısı"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Sonuç olarak, veri tablosuna daha fazla sütun eklemek `meta` map'inde daha fazla metadata kullanılabilir hale getirecek; ancak kanal şeklini değiştirmeyecektir.
Bu, metadata öğelerini girdi belirtimine sabit kodlamak zorunda kalmadan kanalı tüketen süreçler yazmamızı sağlar:

```groovy title="Sözdizimi örneği"
    input:
    tuple val(meta), file(recording)
```

Bu, Nextflow iş akışlarında metadata'yı organize etmek için yaygın olarak kullanılan bir kuraldır.

### Özet

Bu bölümde şunları öğrendiniz:

- **Metadata neden önemlidir:** Metadata'yı verinizle birlikte tutmak, iş akışı boyunca önemli dosya bilgilerini korur.
- **Veri tablolarını nasıl okuyacağınız:** Başlık bilgisine sahip CSV dosyalarını okumak ve satırları yapılandırılmış verilere dönüştürmek için `splitCsv` kullanımı
- **Meta map nasıl oluşturulur:** Tuple yapısı `[ [id:value, ...], file ]` kullanarak metadata'yı dosya verisinden ayırma

---

## 2. Metadata'yı manipüle etme

Artık metadata'mız yüklendiğine göre, onunla bir şeyler yapalım!

Her yaratığın kayıt dosyasında bulunan dili tanımlamak için [`langid`](https://github.com/saffsd/langid.py) adlı bir araç kullanacağız.
Araç bir dizi dil üzerinde önceden eğitilmiş olarak gelir; bir metin parçası verildiğinde `stdout`'a bir dil tahmini ve ilişkili bir olasılık puanı çıktısı verir.

### 2.1. Süreci içe aktarın ve kodu inceleyin

Size `langid` aracını saran `IDENTIFY_LANGUAGE` adlı önceden yazılmış bir süreç modülü sağlıyoruz; bu nedenle sadece iş akışı bloğundan önce bir include ifadesi eklemeniz gerekiyor.

İş akışında aşağıdaki düzenlemeyi yapın:

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

Gördüğünüz gibi, girdi tanımı girdi kanalımıza uyguladığımız aynı `tuple val(meta), path(file)` yapısını kullanıyor.

Çıktı tanımı, girdi yapısına benzer bir yapıya sahip bir tuple olarak yapılandırılmıştır; ancak üçüncü öğe olarak `stdout` da içerir.
Bu `tuple val(meta), path(file), <output>` modeli, metadata'yı hem girdi verisiyle hem de çıktılarla ilişkili tutar ve pipeline boyunca akar.

Burada Nextflow'un [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) çıktı niteleyicisini kullandığımızı unutmayın; çünkü araç çıktısını dosya yazmak yerine doğrudan konsola yazdırır. Komut satırında olasılık puanını kaldırmak, yeni satır karakterlerini kaldırarak dizeyi temizlemek ve yalnızca dil tahminini döndürmek için `sed` kullanırız.

### 2.2. `IDENTIFY_LANGUAGE`'e bir çağrı ekleyin

Artık süreç iş akışı için kullanılabilir olduğuna göre, veri kanalında çalıştırmak için `IDENTIFY_LANGUAGE` sürecine bir çağrı ekleyebiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

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

Orijinal `.view()` işlemini kanal oluşturmadan kaldırdığımıza dikkat edin.

Artık iş akışını çalıştırabiliriz:

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

Daha önce belirtildiği gibi, çıktıya girdi dosyasını ve meta map'i de dahil ettik; bu da her ikisinin de az önce ürettiğimiz yeni bilgiyle ilişkili kalması anlamına gelir.
Bu, bir sonraki adımda faydalı olacaktır.

!!! note "Not"

    Daha genel olarak, meta map'i sonuçlarla ilişkili tutmanın bu modeli, aynı tanımlayıcıları paylaşan ilgili sonuçları ilişkilendirmeyi kolaylaştırır.

    Zaten öğrenmiş olduğunuz gibi, sonuçları eşleştirmek için kanallardaki öğelerin sırasına güvenemezsiniz.
    Bunun yerine, verileri doğru şekilde ilişkilendirmek için anahtarlar kullanmalısınız; meta map'ler bu amaç için ideal bir yapı sağlar.

    Bu kullanım senaryosunu [Splitting & Grouping](./splitting_and_grouping.md) yan görevinde ayrıntılı olarak inceliyoruz.

### 2.3. Süreç çıktılarıyla metadata'yı genişletme

Az önce ürettiğimiz sonuçların kendilerinin dosyaların içeriği hakkında bir metadata türü olduğu göz önüne alındığında, bunları meta map'imize eklemek faydalı olacaktır.

Ancak mevcut meta map'i yerinde değiştirmek istemiyoruz.
Teknik bir bakış açısından bunu yapmak _mümkündür_, ancak güvenli değildir.

Bu nedenle, bunun yerine `+` operatörünü (bir Groovy özelliği) kullanarak mevcut meta map'in içeriğini artı yeni bilgileri tutan yeni bir `lang: lang_id` anahtar-değer çiftini içeren yeni bir meta map oluşturacağız.
Bunu eski map'i yenisiyle değiştirmek için bir [`map`](https://www.nextflow.io/docs/latest/operator.html#map) işlemiyle birleştireceğiz.

İş akışında yapmanız gereken düzenlemeler şunlardır:

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

    Diyelim ki `map1`'den tekrar başlıyorsunuz, ancak dil tahmini kendi map'inde değil (`map2` yok).
    Bunun yerine, `lang_id` adlı bir değişkende tutuluyor ve değerini (`'fr'`) `lang` anahtarıyla saklamak istediğinizi biliyorsunuz.

    Aslında şunu yapabilirsiniz:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Burada `[lang: new_info]` anında yeni bir adsız map oluşturur; `map1 + ` ise `map1`'i yeni adsız map ile birleştirerek daha önce olduğu gibi aynı `new_map` içeriğini üretir.

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
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    şuna:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    `map1`'i `meta` olarak değiştirirsek, bunun temelde iş akışımızda meta map'imize dil tahminini eklemek için ihtiyacımız olan her şey olduğunu görebilirsiniz.

    Bir şey dışında!

    İş akışımız söz konusu olduğunda, **`meta, file, lang_id`'den oluşan tuple'da `file` nesnesinin varlığını da hesaba katmamız gerekir**.

    Yani buradaki kod şöyle olur:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    `file`'ın `map` işleminde neden dolaştığını anlamakta zorlanıyorsanız, `[meta + [lang: lang_id], file]` yerine o satırın `[new_map, file]` olduğunu hayal edin.
    Bu, `file`'ı tuple'daki ikinci pozisyondaki orijinal yerinde bıraktığımızı daha açık hale getirir. Sadece `lang_id` değerini aldık ve onu birinci pozisyondaki map'e kattık.

    **Ve bu bizi `tuple val(meta), path(file)` kanal yapısına geri getiriyor!**

Bu kodun ne yaptığını anladığınızdan emin olduğunuzda, çalışıp çalışmadığını görmek için iş akışını çalıştırın:

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
Sürecin çıktısını `meta, file, lang_id`'den düzgün bir şekilde yeniden düzenledik; böylece `lang_id` artık meta map'teki anahtarlardan biri ve kanalın tuple'ları tekrar `meta, file` modeline uyuyor.

<!-- TODO (future) subMap kullanarak bir anahtarın nasıl kaldırılacağını da göstermeli miyiz?! Veya bunu nerede bulacağımızı not edelim. -->

### 2.4. Koşullu ifadeler kullanarak bir dil grubu atama

Artık dil tahminlerimiz olduğuna göre, yeni gruplamalar atamak için bu bilgileri kullanalım.

Örnek verilerimizde, karakterlerimiz tarafından kullanılan diller cermen dilleri (İngilizce, Almanca) ve roman dilleri (Fransızca, İspanyolca, İtalyanca) olarak gruplandırılabilir.
Pipeline'da daha sonra bu sınıflandırmanın hazır bir şekilde mevcut olması faydalı olabilir; bu yüzden o bilgiyi meta map'e ekleyelim.

İyi haber: bu, `map` operatörünü kullanmaya mükemmel şekilde uyan bir başka durumdur!

Özellikle, `lang_group` adında bir değişken tanımlayacağız ve her veri parçası için `lang_group`'a hangi değerin atanacağını belirlemek için bazı basit koşullu mantık kullanacağız.

Genel sözdizimi şöyle görünecek:

```groovy
.map { meta, file ->

    // lang_group'u tanımlayan koşullu mantık buraya gelir

    [meta + [lang_group: lang_group], file]
}
```

Gördüğünüz gibi bu, önceki adımda kullandığımız anında map birleştirme işlemine çok benziyor.
Sadece koşullu ifadeleri yazmamız gerekiyor.

Uygulamak istediğimiz koşullu mantık şudur:

- Varsayılan değeri `'unknown'` olan `lang_group` adında bir değişken tanımlayın.
- Eğer `lang` Almanca (`'de'`) veya İngilizce (`'en'`) ise, `lang_group`'u `germanic` olarak değiştirin.
- Değilse eğer `lang` Fransızca (`'fr'`), İspanyolca (`'es'`) ve İtalyanca (`'it'`) içeren bir listede yer alıyorsa, `lang_group`'u `romance` olarak değiştirin.

Nextflow'da koşullu ifadeleri nasıl yazacağınızı zaten biliyorsanız kendiniz yazmayı deneyin.

!!! tip "İpucu"

    Map işlemi içinde `lang` değerine `meta.lang` ile erişebilirsiniz.

İş akışında aşağıdaki değişiklikleri yapmayı tamamlamalısınız:

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

Hepsi mantıklı olduğunda, sonucu görmek için iş akışını tekrar çalıştırın:

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

Gördüğünüz gibi, kanal öğeleri `[meta, file]` yapısını korur; ancak meta map artık bu yeni sınıflandırmayı içerir.

### Özet

Bu bölümde şunları öğrendiniz:

- **Girdi metadata'sını çıktı kanallarına uygulama**: Metadata'yı bu şekilde kopyalamak, sonuçları daha sonra metadata içeriğine göre ilişkilendirmemize olanak tanır.
- **Özel anahtarlar oluşturma**: Meta map'inizde iki yeni anahtar oluşturdunuz; bunları `meta + [new_key:value]` ile mevcut meta map'e birleştirerek. Biri bir süreçten hesaplanan bir değere, diğeri `map` operatöründe belirlediğiniz bir koşula dayalı olarak.

Bunlar, pipeline boyunca ilerlerken yeni ve mevcut metadata'yı dosyalarla ilişkilendirmenize olanak tanır.
Bir sürecin parçası olarak metadata kullanmasanız bile, meta map'i bu şekilde veriyle ilişkili tutmak tüm ilgili bilgileri bir arada tutmayı kolaylaştırır.

---

## 3. Bir süreçte meta map bilgisini kullanma

Artık meta map'i nasıl oluşturup güncelleyeceğinizi bildiğinize göre, gerçekten eğlenceli kısma geçebiliriz: metadata'yı bir süreçte gerçekten kullanma.

Daha spesifik olarak, iş akışımıza her hayvanı ASCII art olarak çizmek ve kaydedilen metni bir konuşma balonunda söyletmek için ikinci bir adım ekleyeceğiz.
Bunu [`cowpy`](https://github.com/jeffbuttars/cowpy) adlı bir araç kullanarak yapacağız.

??? info "`cowpy` ne yapar?"

    `cowpy`, keyfi metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII art üreten bir komut satırı aracıdır.
    Tony Monroe'nun klasik [cowsay](https://en.wikipedia.org/wiki/Cowsay) aracının Python implementasyonudur.

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
Görmediyseniz endişelenmeyin; ilerlerken bilmeniz gereken her şeyi ele alacağız.

### 3.1. Süreci içe aktarın ve kodu inceleyin

Size `cowpy` aracını saran `COWPY` adlı önceden yazılmış bir süreç modülü sağlıyoruz; bu nedenle sadece iş akışı bloğundan önce bir include ifadesi eklemeniz gerekiyor.

İş akışında aşağıdaki düzenlemeyi yapın:

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

Gördüğünüz gibi, bu süreç şu anda bir girdi dosyası (görüntülenecek metni içeren) ve ASCII art'ta çizilmesi gereken karakteri belirten bir değer almak üzere tasarlanmıştır; genellikle iş akışı seviyesinde bir komut satırı parametresi ile sağlanır.

### 3.2. Bir meta map alanını girdi olarak geçirme

Hello Nextflow kursunda `cowpy` aracını kullandığımızda, nihai görüntüyü çizmek için hangi karakterin kullanılacağını belirlemek için bir komut satırı parametresi kullandık.
Bu mantıklıydı; çünkü pipeline'ın her çalıştırması başına sadece bir görüntü oluşturuyorduk.

Ancak bu eğitimde, işlediğimiz her konu için uygun bir görüntü oluşturmak istiyoruz; bu nedenle bir komut satırı parametresi kullanmak çok sınırlayıcı olurdu.

İyi haber: veri tablomuzda ve dolayısıyla meta map'imizde bir `character` sütunumuz var.
Sürecin her girdi için kullanması gereken karakteri ayarlamak için bunu kullanalım.

Bu amaçla üç şey yapmamız gerekecek:

1. Üzerinde daha rahat çalışabilmemiz için önceki süreçten çıkan çıktı kanalına bir isim verin
2. İlgilendiğimiz bilgilere nasıl erişeceğimizi belirleyin
3. İkinci sürece bir çağrı ekleyin ve bilgiyi uygun şekilde besleyin

Başlayalım.

#### 3.2.1. Önceki çıktı kanalını adlandırma

Önceki manipülasyonları doğrudan ilk sürecin çıktı kanalı üzerinde uyguladık: `IDENTIFY_LANGUAGE.out`.
Bir sonraki sürece kanalın içeriğini beslemek için (ve bunu açık ve kolay okunur bir şekilde yapmak için) ona kendi adını vermek istiyoruz: `ch_languages`.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana iş akışında, `.view()` operatörünü `.set { ch_languages }` ile değiştirin ve kanala ismiyle başvurabildiğimizi test eden bir satır ekleyin.

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

Bu, artık kanala ismiyle başvurabildiğimizi doğruluyor.

#### 3.2.2. Dosya ve karakter metadata'sına erişim

Modül koduna bakarak `COWPY` sürecinin bir metin dosyası ve bir `character` değeri verilmesini beklediğini biliyoruz.
`COWPY` süreç çağrısını yazmak için, kanaldaki her öğeden karşılık gelen dosya nesnesini ve metadata'yı nasıl çıkaracağımızı bilmemiz gerekiyor.

Genellikle olduğu gibi, bunu yapmanın en basit yolu bir `map` işlemi kullanmaktır.

Kanalımız `[meta, file]` olarak yapılandırılmış tuple'lar içerir; bu nedenle `file` nesnesine doğrudan erişebiliriz ve meta map içinde saklanan `character` değerine `meta.character` olarak başvurarak erişebiliriz.

Ana iş akışında aşağıdaki kod değişikliklerini yapın:

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

`.view` işlemlerinin çıktısını daha okunabilir hale getirmek için closure'ları (örneğin `{ file -> "File: " + file }`) kullandığımıza dikkat edin.

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

Bu, kanaldaki her öğe için dosyaya ve karaktere erişebildiğimizi doğruluyor.

#### 3.2.3. `COWPY` sürecini çağırma

Şimdi hepsini bir araya getirelim ve gerçekten `ch_languages` kanalında `COWPY` sürecini çağıralım.

Ana iş akışında aşağıdaki kod değişikliklerini yapın:

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

Sadece iki map işlemini (`.view()` ifadeleri hariç) süreç çağrısına girdi olarak kopyaladığımızı görüyorsunuz.
Aralarındaki virgülü unutmadığınızdan emin olun!

Biraz hantal; ancak bir sonraki bölümde bunu nasıl iyileştireceğimizi göreceğiz.

Bunu çalıştıralım:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Results dizinine bakarsanız, her selamlaşmanın karşılık gelen karakter tarafından söylendiği ASCII art'ı içeren bireysel dosyaları görmelisiniz.

??? abstract "Dizin ve örnek dosya içeriği"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Bu, meta map'teki bilgileri pipeline'ın ikinci adımında komutu parametrelemek için kullanabildığimizi gösteriyor.

Ancak yukarıda belirtildiği gibi, ilgili kodun bir kısmı biraz hantaldı; çünkü metadata'yı hâlâ iş akışı gövdesi bağlamındayken açmamız gerekti.
Bu yaklaşım meta map'ten az sayıda alan kullanmak için iyi çalışır; ancak çok daha fazlasını kullanmak istesek kötü ölçeklenirdi.

`multiMap()` adlı başka bir operatör var ki bu biraz basitleştirmemize olanak tanır; ancak o zaman bile ideal değildir.

??? info "(İsteğe bağlı) `multiMap()` ile alternatif versiyon"

    Merak ediyorsanız, hem `file` hem de `character`'ı çıkaran tek bir `map()` işlemi yazamadık; çünkü bu onları bir tuple olarak döndürürdü.
    `file` ve `character` öğelerini sürece ayrı ayrı beslemek için iki ayrı `map()` işlemi yazmak zorunda kaldık.

    Teknik olarak bunu tek bir eşleme işlemiyle yapmanın başka bir yolu var: birden fazla kanal yayınlayabilen `multiMap()` operatörünü kullanarak.
    Örneğin, yukarıdaki `COWPY` çağrısını aşağıdaki kodla değiştirebilirsiniz:

    === "Sonra"

        ```groovy title="main.nf" linenums="34"
            // ASCII art oluşturmak için cowpy çalıştır
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Önce"

        ```groovy title="main.nf" linenums="34"
            // ASCII art oluşturmak için cowpy çalıştır
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Bu tam olarak aynı sonucu üretir.

Her iki durumda da, iş akışı seviyesinde açma işlemi yapmak zorunda olmak garip.

Tüm meta map'i sürece besleyebilsek ve ihtiyacımız olanı orada seçebilsek daha iyi olurdu.

### 3.3. Tüm meta map'i geçirme ve kullanma

Meta map'in amacı sonuçta tüm metadata'yı bir paket olarak birlikte geçirmektir.
Yukarıda bunu yapamamanın tek nedeni, sürecin bir meta map kabul edecek şekilde ayarlanmamış olmasıydı.
Ancak süreç kodunu kontrol ettiğimiz için bunu değiştirebiliriz.

İş akışını basitleştirebilmemiz için `COWPY` sürecini, ilk süreçte kullandığımız `[meta, file]` tuple yapısını kabul edecek şekilde değiştirelim.

Bu amaçla üç şey yapmamız gerekecek:

1. `COWPY` süreç modülünün girdi tanımlarını değiştirmek
2. Meta map'i kullanmak için süreç komutunu güncellemek
3. İş akışı gövdesindeki süreç çağrısını güncellemek

Hazır mısınız? Başlayalım!

#### 3.3.1. `COWPY` modül girdisini değiştirme

`cowpy.nf` modül dosyasında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Önce"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Bu, eğitimde daha önce ele aldığımız `[meta, file]` tuple yapısını kullanmamızı sağlar.

Eğitimi kısa tutmak için süreç çıktı tanımını meta map'i çıktılayacak şekilde güncellemediğimizi; ancak `IDENTIFY_LANGUAGE` sürecinin modelini izleyerek bunu kendiniz bir alıştırma olarak yapmakta özgür olduğunuzu unutmayın.

#### 3.3.2. Meta map alanını kullanmak için komutu güncelleme

Tüm meta map artık süreç içinde kullanılabilir; bu nedenle içerdiği bilgilere doğrudan komut bloğunun içinden başvurabiliriz.

`cowpy.nf` modül dosyasında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Önce"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Daha önce bağımsız bir girdi olarak geçirilen `character` değerine yapılan referansı, `meta.character` kullanarak başvurduğumuz meta map'teki değerle değiştirdik.

Şimdi süreç çağrısını buna göre güncelleyelim.

#### 3.3.3. Süreç çağrısını güncelleme ve çalıştırma

Süreç artık girdisinin `[meta, file]` tuple yapısını kullanmasını bekliyor; bu da önceki sürecin çıktıladığı şeydir. Bu nedenle `ch_languages` kanalının tamamını `COWPY` sürecine basitçe geçirebiliriz.

Ana iş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // ASCII art oluşturmak için cowpy çalıştır
    COWPY(ch_languages)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // ASCII art oluşturmak için cowpy çalıştır
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Bu, çağrıyı önemli ölçüde basitleştirir!

Önceki yürütmenin sonuçlarını silelim ve çalıştıralım:

```bash
rm -r results
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Results dizinine bakarsanız, öncekiyle aynı çıktıları görmelisiniz; _yani_ her selamlaşmanın karşılık gelen karakter tarafından söylendiği ASCII art'ı içeren bireysel dosyalar.

??? abstract "Dizin içeriği"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Bu, daha basit kodla öncekiyle aynı sonuçları üretir.

Elbette bu, süreç kodunu değiştirebildiğinizi varsayar.
Bazı durumlarda, değiştirme özgürlüğünüz olmayan mevcut süreçlere güvenmek zorunda kalabilirsiniz; bu da seçeneklerinizi sınırlar.
İyi haber: [nf-core](https://nf-co.re/) projesinden modüller kullanmayı planlıyorsanız, nf-core modüllerinin hepsinin standart olarak `[meta, file]` tuple yapısını kullanacak şekilde ayarlanmış olmasıdır.

### 3.4. Eksik gerekli girdilerde sorun giderme

`character` değeri, `COWPY` sürecinin başarılı bir şekilde çalışması için gereklidir.
Bir yapılandırma dosyasında bunun için varsayılan bir değer belirlemiyorsak, veri tablosunda bunun için bir değer SAĞLAMALIYIZ.

**Sağlamazsak ne olur?**
Bu, girdi veri tablosunun ne içerdiğine ve iş akışının hangi versiyonunu çalıştırdığımıza bağlıdır.

#### 3.4.1. Character sütunu var ama boş

Bir veri toplama hatasını simüle etmek için veri tablomuzda girdilerden birinin character değerini sildiğimizi varsayalım:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Yukarıda kullandığımız iş akışının her iki versiyonu için de, veri tablosu okunduğunda tüm girdiler için `character` anahtarı oluşturulacak; ancak `sampleA` için değer boş bir dize olacaktır.

Bu bir hataya neden olacaktır.

??? failure "Komut çıktısı"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Nextflow o örnek için `cowpy` komut satırını çalıştırdığında, `${meta.character}` boş bir dizeyle doldurulur; bu nedenle `cowpy` aracı `-c` argümanı için değer sağlanmadığını söyleyen bir hata fırlatır.

#### 3.4.2. Character sütunu veri tablosunda mevcut değil

Şimdi `character` sütununu veri tablomuzdan tamamen sildiğimizi varsayalım:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Bu durumda veri tablosu okunduğunda `character` anahtarı hiç oluşturulmayacaktır.

##### 3.4.2.1. İş akışı seviyesinde erişilen değer

Bölüm 3.2'de yazdığımız kodun versiyonunu kullanıyorsak, Nextflow `COWPY` sürecini çağırmadan ÖNCE meta map'teki `character` anahtarına erişmeye çalışacaktır.

Talimatla eşleşen öğe bulamayacak; bu nedenle `COWPY`'yi hiç çalıştırmayacaktır.

??? success "Komut çıktısı"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Nextflow açısından bakıldığında, bu iş akışı başarıyla çalıştırıldı!
Ancak istediğimiz çıktıların hiçbiri üretilmeyecektir.

##### 3.4.2.2. Süreç seviyesinde erişilen değer

Bölüm 3.3'teki versiyonu kullanıyorsak, Nextflow tüm meta map'i `COWPY` sürecine geçirecek ve komutu çalıştırmaya çalışacaktır.

Bu bir hataya neden olacaktır; ancak ilk durumla karşılaştırıldığında farklı bir hata.

??? failure "Komut çıktısı"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Bu, `meta.character` mevcut olmadığı için olur; bu nedenle ona erişme girişimimiz `null` döndürür. Sonuç olarak, Nextflow komut satırına tam anlamıyla `null` yerleştirir; bu da tabii ki `cowpy` aracı tarafından tanınmaz.

#### 3.4.3. Çözümler

İş akışı yapılandırmasının bir parçası olarak varsayılan bir değer sağlamanın yanı sıra, bunu daha sağlam bir şekilde ele almak için yapabileceğimiz iki şey vardır:

1. Veri tablosunun gerekli tüm bilgileri içerdiğinden emin olmak için iş akışınıza girdi doğrulaması uygulayın. Hello nf-core eğitim kursunda [girdi doğrulamasına giriş](../hello_nf-core/05_input_validation.md) bulabilirsiniz. <!-- TODO (future) pending a proper Validation side quest -->

2. Süreç modülünüzü kullanan herkesin gerekli girdileri hemen tanımlayabilmesini sağlamak istiyorsanız, gerekli metadata özelliğini açık bir girdi haline de getirebilirsiniz.

İşte bunun nasıl çalışacağının bir örneği.

İlk olarak, süreç seviyesinde girdi tanımını şu şekilde güncelleyin:

=== "Sonra"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Önce"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Ardından, iş akışı seviyesinde `character` özelliğini metadata'dan çıkarmak ve onu girdi tuple'ının açık bir bileşeni yapmak için bir eşleme işlemi kullanın:

=== "Sonra"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Bu yaklaşımın `character`'ın gerekli olduğunu açıkça gösterme avantajı vardır ve sürecin diğer bağlamlarda yeniden kullanılmasını kolaylaştırır.

Bu önemli bir tasarım ilkesini vurgular:

**Meta map'i isteğe bağlı, açıklayıcı bilgiler için kullanın; ancak gerekli değerleri açık girdiler olarak çıkarın.**

Meta map, kanal yapılarını temiz tutmak ve keyfi kanal yapılarını önlemek için mükemmeldir; ancak bir süreçte doğrudan referans verilen zorunlu öğeler için bunları açık girdiler olarak çıkarmak daha sağlam ve bakımı kolay kod oluşturur.

### Özet

Bu bölümde, metadata'yı bir sürecin yürütülmesini özelleştirmek için nasıl kullanacağınızı; hem iş akışı seviyesinde hem de süreç seviyesinde erişerek öğrendiniz.

---

## Ek alıştırma

Bir sürecin içinden meta map bilgilerini kullanma pratiği yapmak istiyorsanız, çıktıların nasıl adlandırıldığını ve/veya organize edildiğini özelleştirmek için meta map'ten `lang` ve `lang_group` gibi diğer bilgileri kullanmayı deneyin.

Örneğin, bu sonucu üretmek için kodu değiştirmeyi deneyin:

```console title="Results dizini içeriği"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Özet

Bu yan görevde, Nextflow iş akışlarında metadata ile etkili bir şekilde nasıl çalışılacağını keşfettiniz.

Metadata'yı açık ve veriye bağlı tutmanın bu modeli, Nextflow'da temel bir en iyi uygulamadır ve dosya bilgilerini sabit kodlamaya göre çeşitli avantajlar sunar:

- Dosya metadata'sı iş akışı boyunca dosyalarla ilişkili kalır
- Süreç davranışı dosya başına özelleştirilebilir
- Çıktı organizasyonu dosya metadata'sını yansıtabilir
- Dosya bilgileri pipeline yürütme sırasında genişletilebilir

Bu modeli kendi çalışmanızda uygulamak, sağlam ve bakımı kolay biyoinformatik iş akışları oluşturmanızı sağlayacaktır.

### Temel modeller

1.  **Metadata'yı Okuma ve Yapılandırma:** CSV dosyalarını okumak ve veri dosyalarınızla ilişkili kalan organize metadata map'leri oluşturmak.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **İş Akışı Sırasında Metadata'yı Genişletme:** Pipeline'ınız ilerledikçe süreç çıktıları ekleyerek ve koşullu mantıkla değerler türeterek metadata'nıza yeni bilgiler ekleme.

    - Süreç çıktısına dayalı yeni anahtarlar ekleme

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Koşullu bir ifade kullanarak yeni anahtarlar ekleme

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Süreç Davranışını Özelleştirme:** Metadata'yı süreç içinde kullanma.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Ek kaynaklar

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
