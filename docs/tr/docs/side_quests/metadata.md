# Metadata ve meta map'ler

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herhangi bir bilimsel analizde, nadiren sadece ham veri dosyalarıyla çalışırız.
Her dosya kendi ek bilgileriyle birlikte gelir: ne olduğu, nereden geldiği ve onu özel kılan nedir.
Bu ekstra bilgiye metadata (üstveri) diyoruz.

Metadata, diğer verileri tanımlayan veridir.
Metadata, dosyalar ve deneysel koşullar hakkında önemli detayları takip eder ve analizleri her veri setinin benzersiz özelliklerine göre uyarlamaya yardımcı olur.

Bunu bir kütüphane kataloğu gibi düşünün: kitaplar gerçek içeriği (ham veri) içerirken, katalog kartları her kitap hakkında temel bilgiler sağlar—ne zaman yayınlandığı, kim yazdığı, nerede bulunacağı (metadata).
Nextflow pipeline'larında metadata şunlar için kullanılabilir:

- Dosyaya özgü bilgileri iş akışı boyunca takip etmek
- Süreçleri dosya özelliklerine göre yapılandırmak
- İlgili dosyaları ortak analiz için gruplamak

### Öğrenme hedefleri

Bu yan görevde, iş akışlarında metadata'yı nasıl ele alacağımızı keşfedeceğiz.
Temel dosya bilgilerini içeren basit bir veri tablosuyla (biyoinformatikte genellikle samplesheet olarak adlandırılır) başlayarak, şunları öğreneceksiniz:

- CSV dosyalarından dosya metadata'sını okuma ve ayrıştırma
- Metadata map'leri oluşturma ve manipüle etme
- İş akışı yürütme sırasında yeni metadata alanları ekleme
- Süreç davranışını özelleştirmek için metadata kullanma

Bu beceriler, karmaşık dosya ilişkilerini ve işleme gereksinimlerini ele alabilen daha sağlam ve esnek pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan görevi üstlenmeden önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramlarını ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabiliyor olmalısınız

---

## 0. Başlarken

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/metadata
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri gözden geçirin

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

Veri tablosu, veri dosyalarının yollarını ve ilişkili bazı metadata'ları 3 sütunda düzenlenmiş şekilde listeler:

- `id`: açıklayıcı, dosyaya verilen bir ID
- `character`: daha sonra farklı yaratıklar çizmek için kullanacağımız bir karakter adı
- `data`: farklı dillerde selamlamalar içeren `.txt` dosyalarının yolları

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

Her veri dosyası, beş dilden birinde (fr: Fransızca, de: Almanca, es: İspanyolca, it: İtalyanca, en: İngilizce) bir selamlama metni içerir.

Ayrıca size `langid` adlı konteynerleştirilmiş bir dil analiz aracı sağlayacağız.

#### Görevi gözden geçirin

Göreviniz şunları yapacak bir Nextflow iş akışı yazmaktır:

1. Her dosyadaki dili otomatik olarak **tanımlayın**
2. Dosyaları dil ailesine göre **gruplayın** (Cermen dilleri vs Roman dilleri)
3. Her dosya için işlemeyi diline ve metadata'sına göre **özelleştirin**
4. Çıktıları dil grubuna göre **düzenleyin**

Bu, dosyaya özgü metadata'nın işleme kararlarını yönlendirdiği tipik bir iş akışı modelini temsil eder; tam olarak metadata map'lerinin zarif bir şekilde çözdüğü türden bir problem.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Bir veri tablosundan metadata yükleme

Başlangıç noktası olarak size verdiğimiz iş akışı taslağını incelemek için `main.nf` iş akışı dosyasını açın.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Örnek veri tablosunu bir dosya olarak yüklemek için temel bir kanal fabrikası kurduğumuzu görebilirsiniz, ancak bu henüz dosyanın içeriğini okumayacaktır.
Bunu ekleyerek başlayalım.

### 1.1. `splitCsv` ile içeriği okuma

Dosya içeriğini bizim tarafımızdan minimum çabayla uygun şekilde ayrıştıracak bir operatör seçmemiz gerekiyor.
Veri tablomuz CSV formatında olduğundan, bu [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörü için bir iştir; bu operatör dosyadaki her satırı kanalda bir eleman olarak yükler.

Kanal oluşturma koduna bir `splitCsv()` işlemi eklemek için aşağıdaki değişiklikleri yapın, ayrıca dosyanın içeriğinin kanala doğru şekilde yüklendiğini kontrol etmek için bir `view()` işlemi ekleyin.

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

Bundan ne çıktığına bakalım, olur mu?
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

Operatörün CSV dosyasındaki her satır için, sütun başlıklarını karşılık gelen değerler için anahtar olarak kullanarak bir anahtar-değer çiftleri map'i oluşturduğunu görebiliriz.

Her map girişi veri tablomuzda bir sütuna karşılık gelir:

- `id`
- `character`
- `recording`

Bu harika! Her dosyadan belirli alanlara erişmeyi kolaylaştırır.
Örneğin, dosya ID'sine `id` ile veya txt dosya yoluna `recording` ile erişebiliriz.

??? info "(İsteğe bağlı) Map'ler hakkında daha fazla bilgi"

    Nextflow'un üzerine inşa edildiği programlama dili olan Groovy'de, map, Python'daki sözlüklere, JavaScript'teki nesnelere veya Ruby'deki hash'lere benzer bir anahtar-değer veri yapısıdır.

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

    Uygun bir `workflow` bloğu olmasa bile, Nextflow bunu bir iş akışıymış gibi çalıştırabilir:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    Ve işte çıktıda görmeyi bekleyebileceğiniz şey:

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

Artık veri tablosunu başarıyla okuduk ve her satırdaki verilere erişimimiz var, pipeline mantığımızı uygulamaya başlayabiliriz.

### 1.3. Metadata'yı bir 'meta map' içinde düzenleme

İş akışının mevcut durumunda, girdi dosyaları (`recording` anahtarı altında) ve ilişkili metadata (`id`, `character`) hepsi aynı düzeydedir, sanki hepsi büyük bir çantadaymış gibi.
Pratik sonuç, bu kanalı tüketen her sürecin bu yapı göz önünde bulundurularak yapılandırılması gerektiğidir:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Veri tablosundaki sütun sayısı değişmediği sürece bu iyidir.
Ancak, veri tablosuna sadece bir sütun bile eklerseniz, kanalın şekli artık sürecin beklediğiyle eşleşmeyecek ve iş akışı hatalar üretecektir.
Ayrıca süreci, biraz farklı girdi verilerine sahip olabilecek başkalarıyla paylaşmayı zorlaştırır ve betik bloğu tarafından gerekmeyen değişkenleri sürece sabit kodlamak zorunda kalabilirsiniz.

Bu sorunu önlemek için, veri tablosunun kaç sütun içerdiğinden bağımsız olarak kanal yapısını tutarlı tutmanın bir yolunu bulmamız gerekiyor.

Bunu, tüm metadata'yı tuple içinde bir öğede toplayarak yapabiliriz, buna metadata map veya daha basitçe 'meta map' diyeceğiz.

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

Kanal elemanlarımızı iki elemandan oluşan bir tuple olarak yeniden yapılandırdık: meta map ve karşılık gelen dosya nesnesi.

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

Şimdi, kanaldaki her eleman önce metadata map'i ve ikinci olarak karşılık gelen dosya nesnesini içerir:

```console title="Örnek çıktı yapısı"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Sonuç olarak, veri tablosuna daha fazla sütun eklemek `meta` map'inde daha fazla metadata kullanılabilir hale getirecek, ancak kanal şeklini değiştirmeyecektir.
Bu, metadata öğelerini girdi belirtimine sabit kodlamak zorunda kalmadan kanalı tüketen süreçler yazmamızı sağlar:

```groovy title="Sözdizimi örneği"
    input:
    tuple val(meta), file(recording)
```

Bu, Nextflow iş akışlarında metadata düzenlemek için yaygın olarak kullanılan bir kuraldır.

### Özet

Bu bölümde şunları öğrendiniz:

- **Metadata neden önemlidir:** Metadata'yı verilerinizle birlikte tutmak, iş akışı boyunca önemli dosya bilgilerini korur.
- **Veri tablolarını nasıl okuyacağınız:** Başlık bilgisi içeren CSV dosyalarını okumak ve satırları yapılandırılmış verilere dönüştürmek için `splitCsv` kullanımı
- **Meta map nasıl oluşturulur:** `[ [id:value, ...], file ]` tuple yapısını kullanarak metadata'yı dosya verisinden ayırma

---

## 2. Metadata'yı manipüle etme

Artık metadata'mızı yüklediğimize göre, onunla bir şeyler yapalım!

Her yaratığın kayıt dosyasında bulunan dili tanımlamak için [`langid`](https://github.com/saffsd/langid.py) adlı bir araç kullanacağız.
Araç bir dizi dil üzerinde önceden eğitilmiş olarak gelir ve bir metin parçası verildiğinde, hem bir dil tahmini hem de ilişkili bir olasılık skoru çıktısı verir, her ikisi de `stdout`'a.

### 2.1. Süreci içe aktarın ve kodu inceleyin

Size `langid` aracını sarmalayan `IDENTIFY_LANGUAGE` adlı önceden yazılmış bir süreç modülü sağlıyoruz, bu nedenle sadece workflow bloğundan önce bir include ifadesi eklemeniz gerekiyor.

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

Gördüğünüz gibi, girdi tanımı girdi kanalımıza uyguladığımız `tuple val(meta), path(file)` yapısını kullanıyor.

Çıktı tanımı, girdininkilere benzer bir yapıya sahip bir tuple olarak yapılandırılmıştır, ancak üçüncü eleman olarak `stdout`'u da içerir.
Bu `tuple val(meta), path(file), <output>` modeli, metadata'yı hem girdi verileri hem de çıktılarla ilişkili tutar ve pipeline boyunca akar.

Burada Nextflow'un [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) çıktı niteleyicisini kullandığımıza dikkat edin çünkü araç çıktısını bir dosya yazmak yerine doğrudan konsola yazdırır; ve olasılık skorunu kaldırmak, yeni satır karakterlerini kaldırarak dizeyi temizlemek ve yalnızca dil tahminini döndürmek için komut satırında `sed` kullanırız.

### 2.2. `IDENTIFY_LANGUAGE`'e bir çağrı ekleyin

Artık süreç iş akışı için kullanılabilir olduğuna göre, veri kanalı üzerinde çalıştırmak için `IDENTIFY_LANGUAGE` sürecine bir çağrı ekleyebiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Her selamlamanın dilini tanımlamak için langid çalıştır
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

Ve daha önce belirtildiği gibi, çıktıya girdi dosyasını ve meta map'i de dahil ettik, bu da her ikisinin de yeni ürettiğimiz bilgiyle ilişkili kalması anlamına gelir.
Bu, bir sonraki adımda faydalı olacaktır.

!!! note

    Daha genel olarak, meta map'i sonuçlarla ilişkili tutma modeli, aynı tanımlayıcıları paylaşan ilgili sonuçları ilişkilendirmeyi kolaylaştırır.

    Zaten öğrenmiş olduğunuz gibi, sonuçları bunlar arasında eşleştirmek için kanallardaki öğelerin sırasına güvenemezsiniz.
    Bunun yerine, verileri doğru şekilde ilişkilendirmek için anahtarlar kullanmalısınız ve meta map'ler bu amaç için ideal bir yapı sağlar.

    Bu kullanım durumunu [Splitting & Grouping](./splitting_and_grouping.md) yan görevinde ayrıntılı olarak inceliyoruz.

### 2.3. Süreç çıktılarıyla metadata'yı genişletme

Az önce ürettiğimiz sonuçlar dosyaların içeriği hakkında bir metadata biçimi olduğundan, bunları meta map'imize eklemek faydalı olacaktır.

Ancak, mevcut meta map'i yerinde değiştirmek istemiyoruz.
Teknik açıdan bakıldığında, bunu yapmak _mümkündür_, ancak güvenli değildir.

Bu nedenle, bunun yerine mevcut meta map'in içeriğini artı yeni bilgiyi tutan yeni bir `lang: lang_id` anahtar-değer çiftini içeren yeni bir meta map oluşturacağız, `+` operatörünü (bir Groovy özelliği) kullanarak.
Ve bunu eski map'i yenisiyle değiştirmek için bir [`map`](https://www.nextflow.io/docs/latest/operator.html#map) işlemiyle birleştireceğiz.

İş akışında yapmanız gereken düzenlemeler şunlardır:

=== "Sonra"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Her selamlamanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Her selamlamanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

`+` operatörüne henüz aşina değilseniz veya bu kafa karıştırıcı görünüyorsa, aşağıdaki ayrıntılı açıklamayı gözden geçirmek için birkaç dakika ayırın.

??? info "`+` operatörü kullanarak yeni meta map oluşturma"

    **İlk olarak, iki map'in içeriğini Groovy operatörü `+` kullanarak birleştirebileceğimizi bilmeniz gerekir.**

    Diyelim ki aşağıdaki map'lere sahibiz:

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

    **Peki zaten bir map'in parçası olmayan bir alan eklemeniz gerekirse ne olur?**

    Diyelim ki `map1`'den tekrar başlıyorsunuz, ancak dil tahmini kendi map'inde değil (`map2` yok).
    Bunun yerine, `lang_id` adlı bir değişkende tutuluyor ve değerini (`'fr'`) `lang` anahtarıyla saklamak istediğinizi biliyorsunuz.

    Aslında şunu yapabilirsiniz:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Burada, `[lang: new_info]` anında yeni bir isimsiz map oluşturur ve `map1 + ` `map1`'i yeni isimsiz map ile birleştirir, daha öncekiyle aynı `new_map` içeriğini üretir.

    Güzel, değil mi?

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

    Çıktı, yukarıdaki örneğimizdeki `new_map` ile aynı içeriğe sahip tek bir isimsiz map'tir.
    Yani etkili bir şekilde şunu dönüştürdük:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    şuna:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Umarım `map1`'i `meta` olarak değiştirirsek, bunun temelde iş akışımızdaki meta map'imize dil tahminini eklemek için ihtiyacımız olan her şey olduğunu görebilirsiniz.

    Bir şey hariç!

    İş akışımız durumunda, **tuple'da `file` nesnesinin varlığını da hesaba katmamız gerekir**, bu `meta, file, lang_id`'den oluşur.

    Yani buradaki kod şöyle olur:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    `file`'ın `map` işleminde neden hareket ediyor gibi göründüğünü takip etmekte zorlanıyorsanız, `[meta + [lang: lang_id], file]` yerine o satırın `[new_map, file]` okuduğunu hayal edin.
    Bu, `file`'ı tuple'daki ikinci konumdaki orijinal yerinde bıraktığımızı daha net hale getirmelidir. Sadece `new_info` değerini aldık ve tuple'daki ilk konumdaki map'e katladık.

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

Evet, bu doğru!
Sürecin çıktısını `meta, file, lang_id`'den düzgün bir şekilde yeniden düzenledik, böylece `lang_id` artık meta map'teki anahtarlardan biri ve kanalın tuple'ları tekrar `meta, file` modeline uyuyor.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Koşullu ifadeler kullanarak bir dil grubu atama

Artık dil tahminlerimize sahip olduğumuza göre, yeni gruplamalar atamak için bilgiyi kullanalım.

Örnek verilerimizde, karakterlerimiz tarafından kullanılan diller cermen dilleri (İngilizce, Almanca) ve roman dilleri (Fransızca, İspanyolca, İtalyanca) olarak gruplandırılabilir.
Bu sınıflandırmanın pipeline'ın ilerleyen bölümlerinde bir yerde hazır bulunması faydalı olabilir, bu yüzden bu bilgiyi meta map'e ekleyelim.

Ve iyi haber, bu yine `map` operatörünü kullanmaya mükemmel şekilde uygun başka bir durum!

Özellikle, `lang_group` adlı bir değişken tanımlayacağız, her veri parçası için `lang_group`'a hangi değerin atanacağını belirlemek için bazı basit koşullu mantık kullanacağız.

Genel sözdizimi şöyle görünecek:

```groovy
.map { meta, file ->

    // lang_group'u tanımlayan koşullu mantık buraya gelir

    [meta + [lang_group: lang_group], file]
}
```

Bunun önceki adımda kullandığımız anında map birleştirme işlemine çok benzediğini görebilirsiniz.
Sadece koşullu ifadeleri yazmamız gerekiyor.

İşte uygulamak istediğimiz koşullu mantık:

- Varsayılan değeri `'unknown'` olan `lang_group` adlı bir değişken tanımlayın.
- Eğer `lang` Almanca (`'de'`) veya İngilizce (`'en'`) ise, `lang_group`'u `germanic` olarak değiştirin.
- Aksi takdirde eğer `lang` Fransızca (`'fr'`), İspanyolca (`'es'`) ve İtalyanca (`'it'`) içeren bir listede yer alıyorsa, `lang_group`'u `romance` olarak değiştirin.

Nextflow'da koşullu ifadeleri nasıl yazacağınızı zaten biliyorsanız, kendiniz yazmayı deneyin.

!!! tip

    Map işlemi içinde `lang` değerine `meta.lang` ile erişebilirsiniz.

İş akışında aşağıdaki değişiklikleri yapmalısınız:

=== "Sonra"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Her selamlamanın dilini tanımlamak için langid çalıştır
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
        // Her selamlamanın dilini tanımlamak için langid çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

İşte ana noktalar:

- `lang_group` değişkenini varsayılan değeri `unknown` olarak ayarlanmış şekilde oluşturmak için `def lang_group = "unknown"` kullanırız.
- Koşullu mantık için bir `if {} else if {}` yapısı kullanırız, iki cermen dili için alternatif `.equals()` testleri ve üç roman dili için bir listede varlık testi ile.
- Güncellenmiş meta map'i oluşturmak için daha önce olduğu gibi `meta + [lang_group:lang_group]` birleştirme işlemini kullanırız.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Bunların hepsi anlamlı olduğunda, sonucu görmek için iş akışını tekrar çalıştırın:

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

Gördüğünüz gibi, kanal elemanları `[meta, file]` yapılarını korur, ancak meta map artık bu yeni sınıflandırmayı içerir.

### Özet

Bu bölümde şunları nasıl yapacağınızı öğrendiniz:

- **Girdi metadata'sını çıktı kanallarına uygulama**: Metadata'yı bu şekilde kopyalamak, sonuçları daha sonra metadata içeriğine göre ilişkilendirmemize olanak tanır.
- **Özel anahtarlar oluşturma**: Meta map'inizde iki yeni anahtar oluşturdunuz, bunları `meta + [new_key:value]` ile mevcut meta map'e birleştirerek. Biri bir süreçten hesaplanan bir değere dayalı, diğeri `map` operatöründe ayarladığınız bir koşula dayalı.

Bunlar, pipeline'ınızda ilerlerken dosyalarla yeni ve mevcut metadata'yı ilişkilendirmenize olanak tanır.
Bir sürecin parçası olarak metadata kullanmasanız bile, meta map'i bu şekilde verilerle ilişkili tutmak tüm ilgili bilgileri bir arada tutmayı kolaylaştırır.

---

## 3. Bir süreçte meta map bilgisini kullanma

Artık meta map'i nasıl oluşturacağınızı ve güncelleyeceğinizi bildiğinize göre, gerçekten eğlenceli kısma geçebiliriz: metadata'yı bir süreçte gerçekten kullanmak.

Daha spesifik olarak, her hayvanı ASCII sanatı olarak çizmek ve kaydedilen metni bir konuşma balonunda söyletmek için iş akışımıza ikinci bir adım ekleyeceğiz.
Bunu [`cowpy`](https://github.com/jeffbuttars/cowpy) adlı bir araç kullanarak yapacağız.

??? info "`cowpy` ne yapar?"

    `cowpy`, keyifli bir şekilde rastgele metin girdilerini görüntülemek için ASCII sanatı üreten bir komut satırı aracıdır.
    Tony Monroe'nun klasik [cowsay](https://en.wikipedia.org/wiki/Cowsay) aracının bir python uygulamasıdır.

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

Hello Nextflow kursunu tamamladıysanız, bu aracı zaten çalışırken görmüşsünüzdür.
Değilse, endişelenmeyin; ilerledikçe bilmeniz gereken her şeyi ele alacağız.

### 3.1. Süreci içe aktarın ve kodu inceleyin

Size `cowpy` aracını sarmalayan `COWPY` adlı önceden yazılmış bir süreç modülü sağlıyoruz, bu nedenle sadece workflow bloğundan önce bir include ifadesi eklemeniz gerekiyor.

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

Gördüğünüz gibi, bu süreç şu anda bir girdi dosyası (görüntülenecek metni içeren) ve ASCII sanatında çizilmesi gereken karakteri belirten bir değer almak üzere tasarlanmıştır, genellikle iş akışı düzeyinde bir komut satırı parametresi tarafından sağlanır.

### 3.2. Bir meta map alanını girdi olarak geçirme

Hello Nextflow kursunda `cowpy` aracını kullandığımızda, hangi karakterin kullanılacağını belirlemek için bir komut satırı parametresi kullandık.
Bu mantıklıydı, çünkü pipeline'ın her çalıştırmasında sadece bir görüntü oluşturuyorduk.

Ancak, bu eğitimde, işlediğimiz her konu için uygun bir görüntü oluşturmak istiyoruz, bu nedenle bir komut satırı parametresi kullanmak çok sınırlayıcı olurdu.

İyi haber: veri tablomuzda ve dolayısıyla meta map'imizde bir `character` sütunumuz var.
Sürecin her giriş için kullanması gereken karakteri ayarlamak için bunu kullanalım.

Bu amaçla, üç şey yapmamız gerekecek:

1. Önceki süreçten çıkan çıktı kanalına bir isim verin, böylece üzerinde daha rahat çalışabiliriz.
2. İlgilenilen bilgiye nasıl erişeceğimizi belirleyin
3. İkinci sürece bir çağrı ekleyin ve bilgiyi uygun şekilde besleyin.

Hadi başlayalım.

#### 3.2.1. Önceki çıktı kanalını adlandırın

Önceki manipülasyonları doğrudan ilk sürecin çıktı kanalı olan `IDENTIFY_LANGUAGE.out` üzerinde uyguladık.
Bu kanalın içeriğini bir sonraki sürece beslemek için (ve bunu açık ve okunması kolay bir şekilde yapmak için) ona kendi adını vermek istiyoruz, `ch_languages`.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana iş akışında, `.view()` operatörünü `.set { ch_languages }` ile değiştirin ve kanala ismiyle başvurabileceğimizi test eden bir satır ekleyin.

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Her selamlamanın dilini tanımlamak için langid çalıştır
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

        // Geçici: ch_languages'e göz atın
        ch_languages.view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Her selamlamanın dilini tanımlamak için langid çalıştır
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

Bu, artık kanala ismiyle başvurabileceğimizi doğrular.

#### 3.2.2. Dosya ve karakter metadata'sına erişin

Modül koduna baktığımızda `COWPY` sürecinin bir metin dosyası ve bir `character` değeri verilmesini beklediğini biliyoruz.
`COWPY` süreç çağrısını yazmak için, kanaldaki her elemandan karşılık gelen dosya nesnesini ve metadata'yı nasıl çıkaracağımızı bilmemiz gerekiyor.

Genellikle olduğu gibi, bunu yapmanın en basit yolu bir `map` işlemi kullanmaktır.

Kanalımız `[meta, file]` olarak yapılandırılmış tuple'lar içerir, bu nedenle `file` nesnesine doğrudan erişebiliriz ve meta map içinde saklanan `character` değerine `meta.character` olarak başvurarak erişebiliriz.

Ana iş akışında, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34"
        // Geçici: dosya ve karaktere eriş
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="34"
        // Geçici: ch_languages'e göz atın
        ch_languages.view()
    ```

`.view` işlemlerinin çıktısını daha okunabilir hale getirmek için closure'lar (`{ file -> "File: " + file }` gibi) kullandığımıza dikkat edin.

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

#### 3.2.3. `COWPY` sürecini çağırın

Şimdi hepsini bir araya getirelim ve `ch_languages` kanalında `COWPY` sürecini gerçekten çağıralım.

Ana iş akışında, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34"
        // ASCII sanatı oluşturmak için cowpy çalıştır
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
Sadece aralarındaki virgülü unutmadığınızdan emin olun!

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

Sonuçlar dizinine bakarsanız, karşılık gelen karakter tarafından konuşulan her selamlamanın ASCII sanatını içeren bireysel dosyaları görmelisiniz.

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

Bu, pipeline'ın ikinci adımındaki komutu parametrize etmek için meta map'teki bilgiyi kullanabildiğimizi gösterir.

Ancak, yukarıda belirtildiği gibi, dahil olan kodun bir kısmı biraz hantaldı, çünkü hala iş akışı gövdesi bağlamındayken meta verileri açmak zorunda kaldık.
Bu yaklaşım, meta map'ten az sayıda alan kullanmak için iyi çalışır, ancak çok daha fazlasını kullanmak istersek kötü ölçeklenirdi.

Bunu biraz düzene sokmamıza izin veren `multiMap()` adlı başka bir operatör var, ancak o zaman bile ideal değil.

??? info "(İsteğe bağlı) `multiMap()` ile alternatif versiyon"

    Merak ediyorsanız, hem `file` hem de `character`'ı çıktılayan tek bir `map()` işlemi yazamadık, çünkü bu onları bir tuple olarak döndürürdü.
    `file` ve `character` elemanlarını sürece ayrı ayrı beslemek için iki ayrı `map()` işlemi yazmak zorunda kaldık.

    Teknik olarak bunu tek bir eşleme işlemi aracılığıyla yapmanın başka bir yolu var, birden fazla kanal yayabilen `multiMap()` operatörünü kullanarak.
    Örneğin, yukarıdaki `COWPY` çağrısını aşağıdaki kodla değiştirebilirsiniz:

    === "Sonra"

        ```groovy title="main.nf" linenums="34"
            // ASCII sanatı oluşturmak için cowpy çalıştır
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Önce"

        ```groovy title="main.nf" linenums="34"
            // ASCII sanatı oluşturmak için cowpy çalıştır
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Bu tamamen aynı sonucu üretir.

Her iki durumda da, iş akışı düzeyinde bir miktar açma yapmak zorunda kalmamız gariptir.

Tüm meta map'i sürece besleyebilsek ve orada ihtiyacımız olanı seçebilsek daha iyi olurdu.

### 3.3. Tüm meta map'i geçirme ve kullanma

Meta map'in amacı sonuçta tüm metadata'yı bir paket olarak geçirmektir.
Yukarıda bunu yapamamızın tek nedeni, sürecin bir meta map kabul edecek şekilde ayarlanmamış olmasıdır.
Ancak süreç kodunu kontrol ettiğimiz için bunu değiştirebiliriz.

İş akışını düzene sokmak için `COWPY` sürecini ilk süreçte kullandığımız `[meta, file]` tuple yapısını kabul edecek şekilde değiştirelim.

Bu amaçla, üç şey yapmamız gerekecek:

1. `COWPY` süreç modülünün girdi tanımlarını değiştirin
2. Meta map'i kullanmak için süreç komutunu güncelleyin
3. İş akışı gövdesindeki süreç çağrısını güncelleyin

Hazır mısınız? Hadi gidelim!

#### 3.3.1. `COWPY` modül girdisini değiştirin

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

Eğitimi düzenli tutmak için süreç çıktı tanımını meta map'i çıktılamak üzere güncellemedik, ancak `IDENTIFY_LANGUAGE` sürecinin modelini takip ederek bunu kendiniz bir alıştırma olarak yapmaktan çekinmeyin.

#### 3.3.2. Meta map alanını kullanmak için komutu güncelleyin

Tüm meta map artık süreç içinde mevcut, bu nedenle içerdiği bilgilere doğrudan komut bloğunun içinden başvurabiliriz.

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

Daha önce bağımsız bir girdi olarak geçirilen `character` değerine yapılan referansı, `meta.character` kullanarak meta map'te tutulan değerle değiştirdik.

Şimdi süreç çağrısını buna göre güncelleyelim.

#### 3.3.3. Süreç çağrısını güncelleyin ve çalıştırın

Süreç artık girdisinin `[meta, file]` tuple yapısını kullanmasını bekliyor, bu da önceki sürecin çıktısıdır, bu nedenle tüm `ch_languages` kanalını `COWPY` sürecine basitçe geçirebiliriz.

Ana iş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // ASCII sanatı oluşturmak için cowpy çalıştır
    COWPY(ch_languages)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // ASCII sanatı oluşturmak için cowpy çalıştır
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

Sonuçlar dizinine bakarsanız, daha önce olduğu gibi aynı çıktıları görmelisiniz, _yani_ karşılık gelen karakter tarafından konuşulan her selamlamanın ASCII sanatını içeren bireysel dosyalar.

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

Yani bu, daha basit kodla daha önce olduğu gibi aynı sonuçları üretir.

Tabii ki, bu süreç kodunu değiştirebileceğinizi varsayar.
Bazı durumlarda, değiştirme özgürlüğünüz olmayan mevcut süreçlere güvenmek zorunda kalabilirsiniz, bu da seçeneklerinizi sınırlar.
İyi haber, [nf-core](https://nf-co.re/) projesinden modüller kullanmayı planlıyorsanız, nf-core modüllerinin hepsinin standart olarak `[meta, file]` tuple yapısını kullanacak şekilde ayarlanmış olmasıdır.

### 3.4. Eksik gerekli girdilerin sorun giderme

`character` değeri `COWPY` sürecinin başarıyla çalışması için gereklidir.
Bir yapılandırma dosyasında bunun için varsayılan bir değer ayarlamazsak, veri tablosunda bunun için bir değer sağlamalıyız.

**Sağlamazsak ne olur?**
Girdi veri tablosunun ne içerdiğine ve iş akışının hangi versiyonunu çalıştırdığımıza bağlıdır.

#### 3.4.1. Karakter sütunu var ancak boş

Diyelim ki bir veri toplama hatasını simüle etmek için veri tablomuzda girişlerden birinin karakter değerini siliyoruz:

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

Yukarıda kullandığımız iş akışının her iki versiyonu için de, veri tablosu okunduğunda tüm girişler için `character` anahtarı oluşturulacaktır, ancak `sampleA` için değer boş bir dize olacaktır.

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

Nextflow o örnek için `cowpy` komut satırını çalıştırdığında, `${meta.character}` `cowpy` komut satırında boş bir dizeyle doldurulur, bu nedenle `cowpy` aracı `-c` argümanı için değer sağlanmadığını söyleyen bir hata verir.

#### 3.4.2. Karakter sütunu veri tablosunda mevcut değil

Şimdi veri tablomuzdan `character` sütununu tamamen sildiğimizi varsayalım:

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

##### 3.4.2.1. Değere iş akışı düzeyinde erişildi

Bölüm 3.2'de yazdığımız kodun versiyonunu kullanıyorsak, Nextflow `COWPY` sürecini çağırmadan ÖNCE meta map'teki `character` anahtarına erişmeye çalışacaktır.

Talimatla eşleşen hiçbir eleman bulamayacak, bu nedenle `COWPY`'yi hiç çalıştırmayacaktır.

??? success "Komut çıktısı"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Nextflow açısından, bu iş akışı başarıyla çalıştı!
Ancak, istediğimiz çıktıların hiçbiri üretilmeyecektir.

##### 3.4.2.2. Değere süreç düzeyinde erişildi

Bölüm 3.3'teki versiyonu kullanıyorsak, Nextflow tüm meta map'i `COWPY` sürecine geçirecek ve komutu çalıştırmaya çalışacaktır.

Bu bir hataya neden olacak, ancak ilk duruma kıyasla farklı bir hata.

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

Bu olur çünkü `meta.character` mevcut değildir, bu nedenle ona erişme girişimimiz `null` döndürür. Sonuç olarak, Nextflow kelimenin tam anlamıyla komut satırına `null` ekler, bu da tabii ki `cowpy` aracı tarafından tanınmaz.

#### 3.4.3. Çözümler

İş akışı yapılandırmasının bir parçası olarak varsayılan bir değer sağlamanın yanı sıra, bunu daha sağlam bir şekilde ele almak için yapabileceğimiz iki şey var:

1. Veri tablosunun tüm gerekli bilgileri içerdiğinden emin olmak için iş akışınıza girdi doğrulaması uygulayın. Hello nf-core eğitim kursunda [girdi doğrulamasına giriş](../hello_nf-core/05_input_validation.md) bulabilirsiniz. <!-- TODO (future) pending a proper Validation side quest -->

2. Süreç modülünüzü kullanan herkesin gerekli girdileri hemen tanımlayabilmesini sağlamak istiyorsanız, gerekli metadata özelliğini açık bir girdi yapabilirsiniz.

İşte bunun nasıl çalışacağına dair bir örnek.

İlk olarak, süreç düzeyinde, girdi tanımını şu şekilde güncelleyin:

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

Ardından, iş akışı düzeyinde, `character` özelliğini metadata'dan çıkarmak ve girdi tuple'ının açık bir bileşeni yapmak için bir eşleme işlemi kullanın:

=== "Sonra"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Bu yaklaşım, `character`'ın gerekli olduğunu açıkça gösterme avantajına sahiptir ve süreci diğer bağlamlarda yeniden dağıtmayı kolaylaştırır.

Bu önemli bir tasarım ilkesini vurgular:

**İsteğe bağlı, açıklayıcı bilgiler için meta map'i kullanın, ancak gerekli değerleri açık girdiler olarak çıkarın.**

Meta map, kanal yapılarını temiz tutmak ve rastgele kanal yapılarını önlemek için mükemmeldir, ancak bir süreçte doğrudan başvurulan zorunlu elemanlar için, bunları açık girdiler olarak çıkarmak daha sağlam ve sürdürülebilir kod oluşturur.

### Özet

Bu bölümde, bir sürecin yürütülmesini özelleştirmek için metadata'yı nasıl kullanacağınızı, iş akışı düzeyinde veya süreç düzeyinde erişerek öğrendiniz.

---

## Ek alıştırma

Bir süreç içinden meta map bilgisini kullanma pratiği yapmak isterseniz, çıktıların nasıl adlandırıldığını ve/veya düzenlendiğini özelleştirmek için meta map'ten `lang` ve `lang_group` gibi diğer bilgi parçalarını kullanmayı deneyin.

Örneğin, bu sonucu üretmek için kodu değiştirmeyi deneyin:

```console title="Sonuçlar dizini içeriği"
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

Bu yan görevde, Nextflow iş akışlarında metadata ile etkili bir şekilde nasıl çalışacağınızı keşfettiniz.

Metadata'yı açık ve veriye ekli tutma modeli, Nextflow'da temel bir en iyi uygulamadır ve dosya bilgilerini sabit kodlamaya göre çeşitli avantajlar sunar:

- Dosya metadata'sı iş akışı boyunca dosyalarla ilişkili kalır
- Süreç davranışı dosya başına özelleştirilebilir
- Çıktı organizasyonu dosya metadata'sını yansıtabilir
- Dosya bilgileri pipeline yürütme sırasında genişletilebilir

Bu modeli kendi çalışmanızda uygulamak, sağlam, sürdürülebilir biyoinformatik iş akışları oluşturmanızı sağlayacaktır.

### Ana modeller

1.  **Metadata'yı Okuma ve Yapılandırma:** CSV dosyalarını okuma ve veri dosyalarınızla ilişkili kalan düzenli metadata map'leri oluşturma.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **İş Akışı Sırasında Metadata'yı Genişletme** Pipeline'ınız ilerledikçe metadata'nıza yeni bilgiler ekleme, süreç çıktıları ekleyerek ve koşullu mantık yoluyla değerler türeterek.

    - Süreç çıktısına dayalı yeni anahtarlar ekleme

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Koşullu bir cümle kullanarak yeni anahtarlar ekleme

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

3.  **Süreç Davranışını Özelleştirme:** Süreç içinde metadata kullanma.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Ek kaynaklar

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
