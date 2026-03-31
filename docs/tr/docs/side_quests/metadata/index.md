# Metadata ve Meta Map'ler

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herhangi bir bilimsel analizde, yalnızca ham veri dosyalarıyla çalışmak nadiren yeterlidir.
Her dosya, kendine özgü ek bilgilerle birlikte gelir: ne olduğu, nereden geldiği ve onu özel kılan nedir.
Bu ekstra bilgiye metadata diyoruz.

Metadata, diğer verileri tanımlayan veridir.
Metadata, dosyalar ve deneysel koşullar hakkındaki önemli ayrıntıları takip eder ve analizlerin her veri setinin kendine özgü özelliklerine göre uyarlanmasına yardımcı olur.

Bunu bir kütüphane kataloğu gibi düşünebilirsiniz: kitaplar gerçek içeriği (ham veri) barındırırken, katalog kartları her kitap hakkında temel bilgileri sağlar; ne zaman yayımlandığı, kimin yazdığı, nerede bulunacağı (metadata).
Nextflow pipeline'larında metadata şu amaçlarla kullanılabilir:

- İş akışı boyunca dosyaya özgü bilgileri takip etmek
- Dosya özelliklerine göre süreçleri yapılandırmak
- İlgili dosyaları ortak analiz için gruplamak

### Öğrenme hedefleri

Bu yan görevde, iş akışlarında metadata'nın nasıl ele alınacağını keşfedeceğiz.
Temel dosya bilgilerini içeren basit bir veri sayfasından (biyoinformatikte genellikle samplesheet olarak adlandırılır) başlayarak şunları öğreneceksiniz:

- CSV dosyalarından dosya metadata'sını okumak ve ayrıştırmak
- Metadata map'leri oluşturmak ve düzenlemek
- İş akışı yürütme sırasında yeni metadata alanları eklemek
- Süreç davranışını özelleştirmek için metadata kullanmak

Bu beceriler, karmaşık dosya ilişkilerini ve işleme gereksinimlerini yönetebilen daha dayanıklı ve esnek pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmanız gerekir:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabilmek.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için gerekli dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/metadata
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

Bir ana iş akışı dosyası ve bir veri sayfası ile birkaç veri dosyası içeren bir `data` dizini bulacaksınız.

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

`main.nf` dosyasındaki iş akışı, kademeli olarak tam işlevli bir iş akışına dönüştüreceğiniz bir taslaktır.

Veri sayfası, veri dosyalarının yollarını ve ilgili metadata'yı 3 sütun halinde listeler:

- `id`: açıklaması kendinden belli, dosyaya verilen bir kimlik
- `character`: daha sonra farklı yaratıklar çizmek için kullanacağımız bir karakter adı
- `data`: farklı dillerde selamlaşma ifadeleri içeren `.txt` dosyalarının yolları

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

Her veri dosyası, beş dilden birinde selamlama metni içerir (fr: Fransızca, de: Almanca, es: İspanyolca, it: İtalyanca, en: İngilizce).

Ayrıca size `langid` adlı konteynerize bir dil analizi aracı da sunacağız.

#### Görevi inceleyin

Göreviniz, aşağıdakileri yapacak bir Nextflow iş akışı yazmaktır:

1. Her dosyadaki dili otomatik olarak **tanımlamak**
2. Dosyaları dil ailesine göre **gruplamak** (Cermen dilleri ve Roman dilleri)
3. Her dosyanın diline ve metadata'sına göre işlemeyi **özelleştirmek**
4. Çıktıları dil grubuna göre **düzenlemek**

Bu, dosyaya özgü metadata'nın işleme kararlarını yönlendirdiği tipik bir iş akışı örüntüsünü temsil eder; tam da metadata map'lerin zarif bir şekilde çözdüğü türden bir problem.

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebildiyseniz, başlayabilirsiniz.

---

## 1. Veri sayfasından metadata yükleme

Başlangıç noktası olarak size verdiğimiz iş akışı taslağını incelemek için `main.nf` iş akışı dosyasını açın.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Örnek veri sayfasını bir dosya olarak yüklemek için temel bir kanal fabrikası kurduğumuzu görebilirsiniz; ancak bu henüz dosyanın içeriğini okumayacaktır.
Bunu ekleyerek başlayalım.

### 1.1. `splitCsv` ile içerikleri okuma

Dosya içeriğini minimum çabayla uygun şekilde ayrıştıracak bir operatör seçmemiz gerekiyor.
Veri sayfamız CSV formatında olduğundan, bu iş için [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörü biçilmiş kaftandır; dosyadaki her satırı kanalda bir eleman olarak yükler.

Kanal oluşturma koduna bir `splitCsv()` işlemi ve dosyanın içeriğinin kanala doğru şekilde yüklenip yüklenmediğini kontrol etmek için bir `view()` işlemi eklemek üzere aşağıdaki değişiklikleri yapın.

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

CSV dosyasının ilk satırını başlık satırı olarak okumak için `header: true` seçeneğini kullandığımıza dikkat edin.

Bunun ne ürettiğine bakalım, değil mi?
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

Operatörün, CSV dosyasındaki her satır için sütun başlıklarını anahtar olarak kullanan anahtar-değer çiftlerinden oluşan bir map oluşturduğunu görebiliriz.

Her map girişi, veri sayfamızdaki bir sütuna karşılık gelir:

- `id`
- `character`
- `recording`

Harika! Bu, her dosyadan belirli alanlara erişmeyi kolaylaştırır.
Örneğin, dosya kimliğine `id` ile veya txt dosya yoluna `recording` ile erişebiliriz.

??? info "(İsteğe bağlı) Map'ler hakkında daha fazla bilgi"

    Nextflow'un üzerine inşa edildiği programlama dili olan Groovy'de, map; Python'daki sözlükler, JavaScript'teki nesneler veya Ruby'deki hash'ler gibi bir anahtar-değer veri yapısıdır.

    Bir map'i nasıl tanımlayabileceğinizi ve içeriğine pratikte nasıl erişebileceğinizi gösteren çalıştırılabilir bir betik:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Basit bir map oluştur
    def my_map = [id:'sampleA', character:'squirrel']

    // Map'in tamamını yazdır
    println "map: ${my_map}"

    // Nokta gösterimi kullanarak tek tek değerlere eriş
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Düzgün bir `workflow` bloğu olmasa da Nextflow bunu bir iş akışıymış gibi çalıştırabilir:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Çıktıda şunları görmeyi bekleyebilirsiniz:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` ile belirli alanları seçme

Diyelim ki veri sayfasındaki `character` sütununa erişmek ve onu yazdırmak istiyoruz.
Kanalımızdaki her öğe üzerinde yinelemek ve map nesnesinden özellikle `character` girişini seçmek için Nextflow `map` operatörünü kullanabiliriz.

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

Başarılı! Her satır için veri sayfamızdan türetilen map yapısından yararlanarak tek tek sütunlardaki değerlere eriştik.

Veri sayfasını başarıyla okuyup her satırdaki verilere erişebildiğimize göre, pipeline mantığımızı uygulamaya başlayabiliriz.

### 1.3. Metadata'yı bir 'meta map' içinde düzenleme

İş akışının mevcut durumunda, girdi dosyaları (`recording` anahtarı altında) ve ilgili metadata (`id`, `character`) hepsi aynı düzeyde, sanki tek büyük bir torbada gibi.
Bunun pratik sonucu, bu kanalı kullanan her sürecin bu yapıyı göz önünde bulundurarak yapılandırılması gerektiğidir:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Veri sayfasındaki sütun sayısı değişmediği sürece bu sorun değil.
Ancak veri sayfasına tek bir sütun bile eklerseniz, kanalın şekli artık sürecin beklediğiyle eşleşmeyecek ve iş akışı hata üretecektir.
Ayrıca süreci, biraz farklı girdi verilerine sahip olabilecek başkalarıyla paylaşmayı zorlaştırır ve betik bloğu tarafından gereksinim duyulmayan değişkenleri sürece sabit kodlamak zorunda kalabilirsiniz.

Bu sorunu önlemek için, veri sayfasının kaç sütun içerdiğinden bağımsız olarak kanal yapısını tutarlı tutmanın bir yolunu bulmamız gerekiyor.

Bunu, tüm metadata'yı demet içindeki bir öğede toplayarak yapabiliriz; buna metadata map veya daha kısaca 'meta map' diyeceğiz.

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

Kanal elemanlarını, meta map ve karşılık gelen dosya nesnesinden oluşan iki elemanlı bir demete yeniden yapılandırdık.

İş akışını çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console title="View meta map"
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

Artık kanaldaki her eleman, önce metadata map'i, ardından karşılık gelen dosya nesnesini içeriyor:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Sonuç olarak, veri sayfasına daha fazla sütun eklemek `meta` map'inde daha fazla metadata kullanılabilir hale getirecek, ancak kanal şeklini değiştirmeyecektir.
Bu, metadata öğelerini girdi tanımına sabit kodlamak zorunda kalmadan kanalı kullanan süreçler yazmamızı sağlar:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Bu, Nextflow iş akışlarında metadata'yı düzenlemek için yaygın olarak kullanılan bir kuraldır.

### Özetle

Bu bölümde şunları öğrendiniz:

- **Metadata neden önemlidir:** Metadata'yı verilerinizle birlikte tutmak, iş akışı boyunca önemli dosya bilgilerini korur.
- **Veri sayfaları nasıl okunur:** Başlık bilgisiyle CSV dosyalarını okumak ve satırları yapılandırılmış veriye dönüştürmek için `splitCsv` kullanmak.
- **Meta map nasıl oluşturulur:** `[ [id:değer, ...], dosya ]` demet yapısını kullanarak metadata'yı dosya verisinden ayırmak.

---

## 2. Metadata'yı düzenleme

Metadata'mızı yüklediğimize göre, şimdi onunla bir şeyler yapalım!

Her yaratığın kayıt dosyasında hangi dilin bulunduğunu tanımlamak için [`langid`](https://github.com/saffsd/langid.py) adlı bir araç kullanacağız.
Araç, bir dil kümesi üzerinde önceden eğitilmiş olarak gelir ve bir metin parçası verildiğinde, her ikisini de `stdout`'a yazan bir dil tahmini ve ilgili olasılık skoru üretir.

### 2.1. Süreci içe aktarma ve kodu inceleme

`langid` aracını saran `IDENTIFY_LANGUAGE` adlı önceden yazılmış bir süreç modülü sağlıyoruz; bu nedenle yalnızca iş akışı bloğundan önce bir include ifadesi eklemeniz gerekiyor.

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

Gördüğünüz gibi, girdi tanımı girdi kanalımıza az önce uyguladığımız `tuple val(meta), path(file)` yapısını kullanıyor.

Çıktı tanımı, girdinkine benzer bir yapıya sahip bir demet olarak yapılandırılmıştır; ancak üçüncü eleman olarak `stdout` da içermektedir.
Bu `tuple val(meta), path(file), <çıktı>` örüntüsü, metadata'nın pipeline boyunca hem girdi verisiyle hem de çıktılarla ilişkili kalmasını sağlar.

Araç çıktısını bir dosyaya yazmak yerine doğrudan konsola yazdırdığından Nextflow'un [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) çıktı niteleyicisini kullandığımıza dikkat edin; olasılık skorunu kaldırmak, yeni satır karakterlerini temizleyerek dizeyi düzenlemek ve yalnızca dil tahminini döndürmek için komut satırında `sed` kullanıyoruz.

### 2.2. `IDENTIFY_LANGUAGE` çağrısı ekleme

Süreç artık iş akışı için kullanılabilir olduğuna göre, veri kanalı üzerinde çalıştırmak için `IDENTIFY_LANGUAGE` sürecine bir çağrı ekleyebiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
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

Kanal oluşturmadaki orijinal `.view()` işlemini kaldırdığımıza dikkat edin.

Şimdi iş akışını çalıştırabiliriz:

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

Mükemmel! Artık her karakterin hangi dili konuştuğuna dair bir tahminimiz var.

Daha önce belirtildiği gibi, çıktıya girdi dosyasını ve meta map'i de dahil ettik; bu da her ikisinin az önce ürettiğimiz yeni bilgiyle ilişkili kalmaya devam ettiği anlamına gelir.
Bu, bir sonraki adımda işe yarayacaktır.

!!! note "Not"

    Daha genel olarak, meta map'i sonuçlarla ilişkili tutma örüntüsü, aynı tanımlayıcıları paylaşan ilgili sonuçları birbirleriyle ilişkilendirmeyi kolaylaştırır.

    Daha önce öğrendiğiniz gibi, sonuçları kanallar arasında eşleştirmek için öğelerin sırasına güvenemezsiniz.
    Bunun yerine, verileri doğru şekilde ilişkilendirmek için anahtarlar kullanmanız gerekir; meta map'ler bu amaç için ideal bir yapı sağlar.

    Bu kullanım senaryosunu [Splitting & Grouping](../splitting_and_grouping/) yan görevinde ayrıntılı olarak inceliyoruz.

### 2.3. Süreç çıktılarıyla metadata'yı zenginleştirme

Az önce ürettiğimiz sonuçlar, dosyaların içeriği hakkında bir tür metadata olduğundan, bunları meta map'imize eklemek yararlı olacaktır.

Ancak mevcut meta map'i yerinde değiştirmek istemiyoruz.
Teknik açıdan bakıldığında bunu yapmak _mümkün_, ancak güvenli değil.

Bu nedenle, `+` operatörünü (bir Groovy özelliği) kullanarak mevcut meta map'in içeriğini ve yeni bilgileri tutan yeni bir `lang: lang_id` anahtar-değer çiftini içeren yeni bir meta map oluşturacağız.
Bunu, eski map'i yenisiyle değiştirmek için bir [`map`](https://www.nextflow.io/docs/latest/operator.html#map) işlemiyle birleştireceğiz.

İş akışında yapmanız gereken düzenlemeler şunlardır:

=== "Sonra"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

`+` operatörüne henüz aşina değilseniz veya bu kafa karıştırıcı geliyorsa, aşağıdaki ayrıntılı açıklamayı incelemek için birkaç dakika ayırın.

??? info "`+` operatörü kullanılarak yeni meta map'in oluşturulması"

    **Öncelikle, iki map'in içeriğini Groovy `+` operatörünü kullanarak birleştirebileceğimizi bilmeniz gerekiyor.**

    Diyelim ki şu map'lere sahibiz:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Bunları şu şekilde birleştirebiliriz:

    ```groovy
    new_map = map1 + map2
    ```

    `new_map`'in içeriği şöyle olacaktır:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Harika!

    **Peki ya bir map'in parçası olmayan bir alan eklemeniz gerekiyorsa?**

    Diyelim ki `map1`'den yeniden başlıyorsunuz, ancak dil tahmini kendi map'inde değil (yani `map2` yok).
    Bunun yerine, `lang_id` adlı bir değişkende tutuluyor ve değerini (`'fr'`) `lang` anahtarıyla saklamak istiyorsunuz.

    Aslında şunu yapabilirsiniz:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Burada `[lang: new_info]` anında yeni isimsiz bir map oluşturur ve `map1 + ` ifadesi `map1`'i yeni isimsiz map ile birleştirerek daha önce olduğu gibi aynı `new_map` içeriğini üretir.

    Güzel, değil mi?

    **Şimdi bunu bir Nextflow `channel.map()` işlemi bağlamına taşıyalım.**

    Kod şu hale gelir:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Bu şunları yapar:

    - `map1, lang_id ->` demetteki iki öğeyi alır
    - `[map1 + [lang: lang_id]]` yukarıda ayrıntılandırıldığı gibi yeni map'i oluşturur

    Çıktı, yukarıdaki örneğimizdeki `new_map` ile aynı içeriğe sahip tek bir isimsiz map'tir.
    Yani şunu dönüştürdük:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    buna:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    `map1`'i `meta` olarak değiştirirsek, iş akışımızdaki meta map'e dil tahminini eklemek için temelde ihtiyacımız olan tek şeyin bu olduğunu görebilirsiniz.

    Tek bir şey dışında!

    İş akışımız söz konusu olduğunda, **`meta, file, lang_id`'den oluşan demetteki `file` nesnesinin varlığını da hesaba katmamız gerekiyor**.

    Bu nedenle buradaki kod şu hale gelir:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    `map` işleminde `file`'ın neden oradan oraya taşınıyor gibi göründüğünü anlamakta güçlük çekiyorsanız, `[meta + [lang: lang_id], file]` satırının `[new_map, file]` olarak okunduğunu hayal edin.
    Bu, `file`'ı demetteki ikinci konumunda bıraktığımızı daha net göstermelidir. Sadece `new_info` değerini alıp birinci konumdaki map'e katlıyoruz.

    **Ve bu bizi `tuple val(meta), path(file)` kanal yapısına geri getiriyor!**

Bu kodun ne yaptığını anladığınızdan emin olduktan sonra, çalışıp çalışmadığını görmek için iş akışını çalıştırın:

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

Evet, bu doğru!
Sürecin çıktısını `meta, file, lang_id`'den düzgün bir şekilde yeniden düzenledik; böylece `lang_id` artık meta map'teki anahtarlardan biri oldu ve kanalın demetleri yeniden `meta, file` modeline uyuyor.

<!-- TODO (gelecek) subMap kullanarak bir anahtarın nasıl kaldırılacağını da göstermeli miyiz?! Veya nerede bulunacağını belirtmeli miyiz. -->

### 2.4. Koşullu ifadeler kullanarak dil grubu atama

Artık dil tahminlerimize sahip olduğumuza göre, bu bilgiyi yeni gruplamalar atamak için kullanalım.

Örnek verilerimizde, karakterlerimizin kullandığı diller Cermen dilleri (İngilizce, Almanca) ve Roman dilleri (Fransızca, İspanyolca, İtalyanca) olarak gruplandırılabilir.
Bu sınıflandırmanın pipeline'ın ilerleyen aşamalarında hazır bulunması yararlı olabilir; bu nedenle bu bilgiyi meta map'e ekleyelim.

Ve iyi haber şu ki, bu da `map` operatörünü kullanmaya mükemmel şekilde uygun bir durum!

Özellikle, `lang_group` adlı bir değişken tanımlayacak ve her veri parçası için `lang_group`'a hangi değerin atanacağını belirlemek için basit koşullu mantık kullanacağız.

Genel sözdizimi şöyle görünecek:

```groovy
.map { meta, file ->

    // lang_group'u tanımlayan koşullu mantık buraya gelir

    [meta + [lang_group: lang_group], file]
}
```

Bunun önceki adımda kullandığımız anında map birleştirme işlemine çok benzediğini görebilirsiniz.
Sadece koşullu ifadeleri yazmamız gerekiyor.

Uygulamak istediğimiz koşullu mantık şu:

- Varsayılan değeri `'unknown'` olan `lang_group` adlı bir değişken tanımlayın.
- `lang` Almanca (`'de'`) veya İngilizce (`'en'`) ise `lang_group`'u `germanic` olarak değiştirin.
- Aksi takdirde `lang`, Fransızca (`'fr'`), İspanyolca (`'es'`) ve İtalyanca (`'it'`) içeren bir listede yer alıyorsa `lang_group`'u `romance` olarak değiştirin.

Nextflow'da koşullu ifadeleri nasıl yazacağınızı zaten biliyorsanız, kendiniz yazmayı deneyin.

!!! tip "İpucu"

    `lang` değerine map işlemi içinde `meta.lang` ile erişebilirsiniz.

İş akışında aşağıdaki değişiklikleri yapmanız gerekecek:

=== "Sonra"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
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
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Temel noktalar şunlardır:

- `lang_group` değişkenini `unknown` varsayılan değeriyle oluşturmak için `def lang_group = "unknown"` kullanıyoruz.
- Koşullu mantık için `if {} else if {}` yapısını kullanıyoruz; iki Cermen dili için alternatif `.equals()` testleri ve üç Roman dili için bir listede varlık testi.
- Güncellenmiş meta map'i oluşturmak için daha önce olduğu gibi `meta + [lang_group:lang_group]` birleştirme işlemini kullanıyoruz.

<!-- TODO (gelecek) İlgili belgelere ek kaynaklar bölümünde not/bağlantı ekle -->

Her şey anlamlı geldiğinde, sonucu görmek için iş akışını tekrar çalıştırın:

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

Gördüğünüz gibi, kanal elemanları `[meta, file]` yapısını korurken meta map artık bu yeni sınıflandırmayı da içeriyor.

### Özetle

Bu bölümde şunları öğrendiniz:

- **Girdi metadata'sını çıktı kanallarına uygulamak**: Metadata'yı bu şekilde kopyalamak, sonuçları daha sonra metadata içeriğine göre ilişkilendirmemizi sağlar.
- **Özel anahtarlar oluşturmak**: Meta map'inizde iki yeni anahtar oluşturdunuz; bunları `meta + [yeni_anahtar:değer]` ile mevcut meta map'e birleştirdiniz. Biri bir süreçten hesaplanan değere, diğeri `map` operatöründe belirlediğiniz bir koşula dayalı.

Bunlar, pipeline'ınızda ilerledikçe yeni ve mevcut metadata'yı dosyalarla ilişkilendirmenizi sağlar.
Metadata'yı bir sürecin parçası olarak kullanmasanız bile, meta map'i bu şekilde verilerle ilişkili tutmak tüm ilgili bilgileri bir arada tutmayı kolaylaştırır.

---

## 3. Bir süreçte meta map bilgisini kullanma

Meta map'i nasıl oluşturacağınızı ve güncelleyeceğinizi öğrendiğinize göre, gerçekten eğlenceli kısma geçebiliriz: metadata'yı bir süreçte gerçekten kullanmak.

Daha özellikle, iş akışımıza her hayvanı ASCII sanatı olarak çizen ve konuşma balonunda kaydedilen metni söyleten ikinci bir adım ekleyeceğiz.
Bunu [`cowpy`](https://github.com/jeffbuttars/cowpy) adlı bir araç kullanarak yapacağız.

??? info "`cowpy` ne yapar?"

    `cowpy`, rastgele metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII sanatı üreten bir komut satırı aracıdır.
    Tony Monroe'nun klasik [cowsay](https://en.wikipedia.org/wiki/Cowsay) aracının Python uygulamasıdır.

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

Hello Nextflow kursunu tamamladıysanız, bu aracı daha önce görmüşsünüzdür.
Görmediyseniz endişelenmeyin; ilerledikçe bilmeniz gereken her şeyi ele alacağız.

### 3.1. Süreci içe aktarma ve kodu inceleme

`cowpy` aracını saran `COWPY` adlı önceden yazılmış bir süreç modülü sağlıyoruz; bu nedenle yalnızca iş akışı bloğundan önce bir include ifadesi eklemeniz gerekiyor.

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

// cowpy ile ASCII sanatı oluştur
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

Gördüğünüz gibi, bu süreç şu anda bir girdi dosyası (görüntülenecek metni içeren) ve ASCII sanatında çizilmesi gereken karakteri belirten bir değer almak üzere tasarlanmıştır; bu değer genellikle iş akışı düzeyinde bir komut satırı parametresiyle sağlanır.

### 3.2. Bir meta map alanını girdi olarak geçirme

Hello Nextflow kursunda `cowpy` aracını kullandığımızda, son görüntüyü çizmek için hangi karakterin kullanılacağını belirlemek amacıyla bir komut satırı parametresi kullandık.
Bu mantıklıydı, çünkü pipeline'ın her çalıştırmasında yalnızca bir görüntü üretiyorduk.

Ancak bu eğitimde, işlediğimiz her özne için uygun bir görüntü oluşturmak istiyoruz; bu nedenle bir komut satırı parametresi kullanmak çok kısıtlayıcı olurdu.

İyi haber: veri sayfamızda ve dolayısıyla meta map'imizde bir `character` sütunu var.
Bunu, sürecin her giriş için kullanması gereken karakteri ayarlamak için kullanalım.

Bu amaçla üç şey yapmamız gerekecek:

1. Önceki süreçten çıkan çıktı kanalına bir isim vermek, böylece üzerinde daha rahat çalışabiliriz.
2. İlgilendiğimiz bilgiye nasıl erişeceğimizi belirlemek.
3. İkinci sürece bir çağrı eklemek ve bilgiyi uygun şekilde beslemek.

Başlayalım.

#### 3.2.1. Önceki çıktı kanalını adlandırma

Önceki düzenlemeleri doğrudan ilk sürecin çıktı kanalı olan `IDENTIFY_LANGUAGE.out` üzerinde uyguladık.
Bu kanalın içeriğini bir sonraki sürece beslemek (ve bunu açık ve okunması kolay bir şekilde yapmak) için ona kendi adını vermek istiyoruz: `ch_languages`.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana iş akışında `.view()` operatörünü `.set { ch_languages }` ile değiştirin ve kanala adıyla başvurabildiğimizi test eden bir satır ekleyin.

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
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
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
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

Bu, kanala artık adıyla başvurabildiğimizi doğrular.

#### 3.2.2. Dosya ve karakter metadata'sına erişme

Modül koduna bakarak `COWPY` sürecinin bir metin dosyası ve bir `character` değeri beklendiğini biliyoruz.
`COWPY` sürecine yapılacak çağrıyı yazmak için, kanaldaki her elemandan karşılık gelen dosya nesnesini ve metadata'yı nasıl çıkaracağımızı bilmemiz yeterli.

Her zaman olduğu gibi, bunu yapmanın en basit yolu bir `map` işlemi kullanmaktır.

Kanalımız `[meta, file]` yapısında demetler içeriyor; bu nedenle `file` nesnesine doğrudan erişebilir ve meta map içinde saklanan `character` değerine `meta.character` olarak başvurabiliriz.

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

`.view` işlemlerinin çıktısını daha okunabilir hale getirmek için closure'lar (örneğin `{ file -> "File: " + file }`) kullandığımıza dikkat edin.

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

_Dosya yolları ve karakter değerleri çıktınızda farklı bir sırada gelebilir._

Bu, kanaldaki her eleman için dosyaya ve karaktere erişebildiğimizi doğrular.

#### 3.2.3. `COWPY` sürecini çağırma

Şimdi her şeyi bir araya getirelim ve `ch_languages` kanalı üzerinde `COWPY` sürecini gerçekten çağıralım.

Ana iş akışında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34"
        // ASCII sanatı oluşturmak için cowpy'yi çalıştır
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

İki map işlemini (`.view()` ifadeleri olmadan) süreç çağrısına girdi olarak kopyaladığımızı görüyorsunuz.
Aralarındaki virgülü unutmamaya dikkat edin!

Biraz hantal, ama bunu bir sonraki bölümde nasıl daha iyi hale getireceğimizi göreceğiz.

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

Sonuçlar dizinine bakarsanız, her selamlamanın karşılık gelen karakter tarafından söylendiği ASCII sanatını içeren tek tek dosyaları görmelisiniz.

??? abstract "Dizin ve örnek dosya içerikleri"

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

Bu, pipeline'ın ikinci adımındaki komutu parametreleştirmek için meta map'teki bilgileri kullanabildğimizi gösteriyor.

Ancak yukarıda belirtildiği gibi, meta verileri iş akışı gövdesi bağlamında açmak zorunda kaldığımızdan, kullanılan kodun bir kısmı biraz hantaldı.
Bu yaklaşım, meta map'ten az sayıda alan kullanmak için gayet iyi çalışır; ancak çok daha fazlasını kullanmak istesek ölçeklenmez.

`multiMap()` adlı başka bir operatör var ve bu işlemi biraz kolaylaştırmamıza izin veriyor; ancak o bile ideal değil.

??? info "(İsteğe bağlı) `multiMap()` ile alternatif sürüm"

    Merak ediyorsanız, hem `file` hem de `character`'ı çıkaran tek bir `map()` işlemi yazamadık; çünkü bu onları bir demet olarak döndürürdü.
    `file` ve `character` elemanlarını sürece ayrı ayrı beslemek için iki ayrı `map()` işlemi yazmak zorunda kaldık.

    Teknik olarak bunu tek bir eşleme işlemiyle yapmanın başka bir yolu var: birden fazla kanal yayınlayabilen `multiMap()` operatörünü kullanmak.
    Örneğin, yukarıdaki `COWPY` çağrısını aşağıdaki kodla değiştirebilirsiniz:

    === "Sonra"

        ```groovy title="main.nf" linenums="34"
            // ASCII sanatı oluşturmak için cowpy'yi çalıştır
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Önce"

        ```groovy title="main.nf" linenums="34"
            // ASCII sanatı oluşturmak için cowpy'yi çalıştır
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Bu tam olarak aynı sonucu üretir.

Her iki durumda da, iş akışı düzeyinde bir miktar açma işlemi yapmak zorunda kalmak hantaldır.

Tüm meta map'i sürece besleyip ihtiyacımız olanı orada seçebilsek daha iyi olurdu.

### 3.3. Tüm meta map'i geçirme ve kullanma

Meta map'in amacı, sonuçta tüm metadata'yı bir paket olarak birlikte geçirmektir.
Yukarıda bunu yapamamanın tek nedeni, sürecin bir meta map kabul edecek şekilde ayarlanmamış olmasıydı.
Ancak süreç kodunu kontrol ettiğimizden, bunu değiştirebiliriz.

İş akışını kolaylaştırmak için `COWPY` sürecini ilk süreçte kullandığımız `[meta, file]` demet yapısını kabul edecek şekilde değiştirelim.

Bu amaçla üç şey yapmamız gerekecek:

1. `COWPY` süreç modülünün girdi tanımlarını değiştirmek
2. Süreç komutunu meta map'i kullanacak şekilde güncellemek
3. İş akışı gövdesindeki süreç çağrısını güncellemek

Hazır mısınız? Başlayalım!

#### 3.3.1. `COWPY` modülü girdisini değiştirme

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

Bu, eğitimin önceki bölümlerinde ele aldığımız `[meta, file]` demet yapısını kullanmamızı sağlar.

Eğitimi sade tutmak için süreç çıktı tanımını meta map'i çıkaracak şekilde güncellemedik; ancak `IDENTIFY_LANGUAGE` sürecini model alarak bunu kendiniz bir alıştırma olarak yapabilirsiniz.

#### 3.3.2. Komutu meta map alanını kullanacak şekilde güncelleme

Tüm meta map artık sürecin içinde kullanılabilir; bu nedenle içerdiği bilgilere doğrudan komut bloğunun içinden başvurabiliriz.

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

Daha önce bağımsız bir girdi olarak geçirilen `character` değerine yapılan başvuruyu, `meta.character` kullanarak meta map'te tutulan değerle değiştirdik.

Şimdi süreç çağrısını buna göre güncelleyelim.

#### 3.3.3. Süreç çağrısını güncelleme ve çalıştırma

Süreç artık girdisinin `[meta, file]` demet yapısını kullanmasını bekliyor; bu da önceki sürecin çıktısıyla aynı yapı. Bu nedenle `ch_languages` kanalının tamamını `COWPY` sürecine geçirebiliriz.

Ana iş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // ASCII sanatı oluşturmak için cowpy'yi çalıştır
    COWPY(ch_languages)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // ASCII sanatı oluşturmak için cowpy'yi çalıştır
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Bu, çağrıyı önemli ölçüde basitleştiriyor!

Önceki çalıştırmanın sonuçlarını silelim ve çalıştıralım:

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

Sonuçlar dizinine bakarsanız, öncekiyle aynı çıktıları görmelisiniz; yani her selamlamanın karşılık gelen karakter tarafından söylendiği ASCII sanatını içeren tek tek dosyalar.

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

Bu, daha basit kodla öncekiyle aynı sonuçları üretiyor.

Elbette bu, süreç kodunu değiştirebildiğinizi varsayar.
Bazı durumlarda, değiştirme iznine sahip olmadığınız mevcut süreçlere güvenmek zorunda kalabilirsiniz; bu da seçeneklerinizi kısıtlar.
[nf-core](https://nf-co.re/) projesinden modüller kullanmayı planlıyorsanız iyi haber şu: nf-core modüllerinin tamamı standart olarak `[meta, file]` demet yapısını kullanacak şekilde ayarlanmıştır.

### 3.4. Eksik zorunlu girdilerin sorun giderimi

`character` değeri, `COWPY` sürecinin başarıyla çalışması için zorunludur.
Bir yapılandırma dosyasında varsayılan bir değer belirlemediğimiz sürece, veri sayfasında bu değeri MUTLAKA sağlamamız gerekir.

**Sağlamazsak ne olur?**
Bu, girdi veri sayfasının ne içerdiğine ve iş akışının hangi sürümünü çalıştırdığımıza bağlıdır.

#### 3.4.1. Karakter sütunu var ama boş

Diyelim ki bir veri toplama hatasını simüle etmek için veri sayfasındaki girişlerden birinin karakter değerini siliyoruz:

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

Yukarıda kullandığımız iş akışının her iki sürümünde de, veri sayfası okunduğunda tüm girişler için `character` anahtarı oluşturulacak; ancak `sampleA` için değer boş bir dize olacaktır.

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

Nextflow bu örnek için `cowpy` komut satırını çalıştırdığında, `${meta.character}` `cowpy` komut satırında boş bir dizeyle doldurulur; bu nedenle `cowpy` aracı `-c` argümanı için değer sağlanmadığını belirten bir hata fırlatır.

#### 3.4.2. Karakter sütunu veri sayfasında hiç yok

Şimdi diyelim ki `character` sütununu veri sayfasından tamamen siliyoruz:

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

Bu durumda, veri sayfası okunduğunda `character` anahtarı hiç oluşturulmayacaktır.

##### 3.4.2.1. İş akışı düzeyinde erişilen değer

3.2. bölümünde yazdığımız kod sürümünü kullanıyorsak, Nextflow `COWPY` sürecini çağırmadan ÖNCE meta map'teki `character` anahtarına erişmeye çalışacaktır.

Talimatla eşleşen herhangi bir eleman bulamayacak; bu nedenle `COWPY`'yi hiç çalıştırmayacaktır.

??? success "Komut çıktısı"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Nextflow açısından bakıldığında, bu iş akışı başarıyla çalıştı!
Ancak istediğimiz çıktıların hiçbiri üretilmeyecektir.

##### 3.4.2.2. Süreç düzeyinde erişilen değer

3.3. bölümündeki sürümü kullanıyorsak, Nextflow tüm meta map'i `COWPY` sürecine geçirecek ve komutu çalıştırmaya çalışacaktır.

Bu bir hataya neden olacak; ancak ilk durumdan farklı bir hata.

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

Bu, `meta.character` mevcut olmadığından erişmeye çalışmamızın `null` döndürmesi nedeniyle oluşur. Sonuç olarak Nextflow, komut satırına tam anlamıyla `null` ekler; bu da `cowpy` aracı tarafından tanınmaz.

#### 3.4.3. Çözümler

Yapılandırma dosyasının bir parçası olarak varsayılan bir değer sağlamanın yanı sıra, bunu daha sağlam bir şekilde ele almak için yapabileceğimiz iki şey var:

1. Veri sayfasının gerekli tüm bilgileri içerdiğinden emin olmak için iş akışınıza girdi doğrulaması uygulayın. Hello nf-core eğitim kursunda [girdi doğrulamaya giriş](../hello_nf-core/05_input_validation.md) bölümünü bulabilirsiniz. <!-- TODO (gelecek) uygun bir Doğrulama yan görevi bekliyor -->

2. Süreç modülünüzü kullanan herkesin gerekli girdileri hemen tanımlayabilmesini sağlamak istiyorsanız, gerekli metadata özelliğini açık bir girdi haline de getirebilirsiniz.

Bunun nasıl çalışacağına dair bir örnek:

Önce süreç düzeyinde, girdi tanımını şu şekilde güncelleyin:

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

Ardından, iş akışı düzeyinde, `character` özelliğini metadata'dan çıkarmak ve girdi demetinin açık bir bileşeni haline getirmek için bir eşleme işlemi kullanın:

=== "Sonra"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Bu yaklaşımın avantajı, `character`'ın zorunlu olduğunu açıkça göstermesi ve süreci diğer bağlamlarda yeniden kullanmayı kolaylaştırmasıdır.

Bu, önemli bir tasarım ilkesini vurgular:

**Meta map'i isteğe bağlı, açıklayıcı bilgiler için kullanın; ancak zorunlu değerleri açık girdiler olarak çıkarın.**

Meta map, kanal yapılarını temiz tutmak ve rastgele kanal yapılarını önlemek için mükemmeldir; ancak bir süreçte doğrudan başvurulan zorunlu elemanlar için bunları açık girdiler olarak çıkarmak daha sağlam ve sürdürülebilir kod oluşturur.

### Özetle

Bu bölümde, metadata'yı bir sürecin yürütülmesini özelleştirmek için nasıl kullanacağınızı öğrendiniz; buna iş akışı düzeyinde veya süreç düzeyinde erişim de dahildir.

---

## Ek alıştırma

Meta map bilgisini bir sürecin içinden kullanmayı pratik yapmak istiyorsanız, çıktıların nasıl adlandırıldığını ve/veya düzenlendiğini özelleştirmek için meta map'teki `lang` ve `lang_group` gibi diğer bilgileri kullanmayı deneyin.

Örneğin, şu sonucu üretmek için kodu değiştirmeyi deneyin:

```console title="Results directory contents"
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

<!-- TODO (gelecek) Çalışılmış çözüm sağla -->
<!-- yeniden adlandırma, süreç içindeki meta'yı kullanmalı -->
<!-- çıktı düzenlemesi, iş akışı çıktılarındaki meta'yı kullanmalı -->

---

## Özet

Bu yan görevde, Nextflow iş akışlarında metadata ile etkili bir şekilde çalışmayı keşfettiniz.

Metadata'yı açık tutma ve verilerle ilişkilendirme örüntüsü, Nextflow'da temel bir en iyi uygulama olup dosya bilgilerini sabit kodlamaya kıyasla çeşitli avantajlar sunar:

- Dosya metadata'sı, iş akışı boyunca dosyalarla ilişkili kalır
- Süreç davranışı dosya başına özelleştirilebilir
- Çıktı düzenlemesi dosya metadata'sını yansıtabilir
- Dosya bilgileri pipeline yürütme sırasında genişletilebilir

Bu örüntüyü kendi çalışmalarınızda uygulamak, sağlam ve sürdürülebilir biyoinformatik iş akışları oluşturmanızı sağlayacaktır.

### Temel örüntüler

1.  **Metadata'yı Okuma ve Yapılandırma:** CSV dosyalarını okumak ve veri dosyalarınızla ilişkili kalan düzenli metadata map'leri oluşturmak.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **İş Akışı Sırasında Metadata'yı Genişletme:** Süreç çıktıları ekleyerek ve koşullu mantık aracılığıyla değerler türeterek pipeline ilerledikçe metadata'ya yeni bilgiler eklemek.

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

3.  **Süreç Davranışını Özelleştirme:** Metadata'yı sürecin içinde kullanmak.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Ek kaynaklar

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Sırada ne var?

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
