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
- "Meta map + veri dosyası" arayüzünün neden yaygın olarak kullanılan bir kural olduğunu anlamak
- İş akışı yürütme sırasında yeni metadata alanları eklemek
- Süreç davranışını özelleştirmek ve çıktıları düzenlemek için metadata kullanmak

Bu beceriler, karmaşık dosya ilişkilerini ve işleme gereksinimlerini yönetebilen daha dayanıklı ve esnek pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmanız gerekir:

- [Hello Nextflow](../../hello_nextflow/index.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabilmek.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

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

Editör, proje dizinine odaklanmış şekilde açılır.

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
- `data`: farklı dillerde selamlama ifadeleri içeren `.txt` dosyalarının yolları

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

Her karakterin kayıt dosyasındaki selamlamayı söylediği ASCII sanatını oluşturmak için [`COWPY`](https://github.com/jeffbuttars/cowpy) adlı bir araç kullanacağız.

??? info "`COWPY` ne yapar?"

    `COWPY`, rastgele metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII sanatı üreten bir komut satırı aracıdır.
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

Bunun yanı sıra, her karakterin hangi dili konuştuğunu belirlemek ve pipeline çıktılarını buna göre düzenlemek için `langid` adlı bir dil analizi aracı kullanacağız.

#### Görevi inceleyin

Göreviniz, aşağıdakileri yapacak bir Nextflow iş akışı yazmaktır:

1. Her karakterin **ASCII sanatını oluşturmak**
2. Çıktıları dil grubuna göre **düzenlemek** (Cermen dilleri ve Roman dilleri)

Bu, dosyaya özgü metadata'nın işleme kararlarını yönlendirdiği tipik bir iş akışı örüntüsünü temsil eder; tam da metadata map'lerin zarif bir şekilde çözdüğü türden bir problem.

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebildiyseniz, başlayabilirsiniz.

---

## 1. Metadata yükleme ve kullanma için temel seçenekler

Başlangıç noktası olarak size verdiğimiz iş akışı taslağını incelemek için `main.nf` iş akışı dosyasını açın.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

[`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörü, dosyadaki her satırı bir kanal elemanı olarak okur.
Bu, başlangıç kursumuz Hello Nextflow'da CSV verilerini yüklemek için kullandığımız yaklaşımın aynısıdır.
Bunun nasıl çalıştığını hatırlamak için [bu bölüme](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) bakabilirsiniz.

`header: true` seçeneğiyle ilk satır sütun başlıkları olarak işlenir; böylece her eleman, sütun adlarını anahtar olarak kullanan anahtar-değer çiftlerinden oluşan bir map'e dönüşür.

Henüz veriler üzerinde herhangi bir süreç çalıştırmadığımızdan, `publish` ve `output` bloklarının yalnızca birer taslak olduğuna dikkat edin.

### 1.1. İş akışını çalıştırın

Her şey yüklendikten sonra kanal içeriğinin nasıl yapılandırıldığını görmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Gördüğünüz gibi, operatör CSV dosyasındaki her satır için sütun başlıklarını anahtar olarak kullanan anahtar-değer çiftlerinden oluşan bir map oluşturmuştur.

Her map girişi, veri sayfamızdaki bir sütuna karşılık gelir:

- `id`
- `character`
- `recording`

Bu, her satırdan belirli alanlara erişmeyi kolaylaştırır.
Örneğin, dosya kimliğine `id` ile veya txt dosya yoluna `recording` ile erişebiliriz.

??? info "(İsteğe bağlı) Groovy map'leri hakkında daha fazla bilgi"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. `map` ile belirli bir alanı seçme

Kanaldaki her öğe üzerinde yinelemek ve nokta gösterimi kullanarak adıyla erişebileceğimiz `character` alanını seçmek için `map` operatörünü kullanacağız.

#### 1.2.1. Map işlemini ekleyin

`character` sütununa erişmek için `.view()` işleminden önce `map` işlemini şu şekilde ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Belirli bir alana bu şekilde erişmek, bunun nasıl çalıştığını hatırlamak istiyorsanız Hello Nextflow'un [bu bölümünde](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) daha ayrıntılı açıklanmaktadır.

#### 1.2.2. İş akışını çalıştırın

Çıkarılan karakter adlarını görüntüleyebildiğinizi doğrulamak için iş akışını çalıştırın.

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Bu, her satır için `character` sütunundaki değerlere erişebildiğimizi gösteriyor.

Şimdi bu verilerle bir şeyler yapalım: `COWPY` kullanarak ASCII sanatı oluşturmak için `character` ve `recording` alanlarını birlikte kullanalım.

### 1.3. `multiMap` ile alt kanallar yayınlama

Size önceden yazılmış bir `COWPY` süreç modülü sağlıyoruz; bu nedenle önce sürecin girdi gereksinimlerini incelemeniz gerekiyor.

Sürecin nasıl göründüğünü görmek için dosyayı açabilirsiniz:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// cowpy ile ASCII sanatı oluştur
process COWPY {

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

Gördüğünüz gibi, süreç iki ayrı girdi alıyor: bir kayıt dosyası ve bir karakter adı.
Her ikisinin de değerleri mevcut; ancak şu anda kanaldaki her elemanın içinde bir arada bulunuyorlar.

Birden fazla alanı ayrı kanallara çıkarmanın bir yolu, tek bir işlemde bir kanalı birden fazla adlandırılmış alt kanala bölen [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap) operatörüdür.

#### 1.3.1. multiMap işlemini ekleyin

`map` işlemini `multiMap` ile değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

`multiMap` bloğu, her satırdan iki adlandırılmış alt kanal (`file` ve `character`) tanımlar; bunlara `ch_datasheet.file` ve `ch_datasheet.character` olarak erişebiliriz.

#### 1.3.2. Alt kanallar üzerinde COWPY'yi çağırın

Şimdi `COWPY` sürecini dahil edin ve her alt kanalı ayrı bir argüman olarak verin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Bu, `COWPY`'nin gerektirdiği şekilde iki alanı ayrı ayrı geçirmemizi sağlar.

#### 1.3.3. Çıktı yayımlamayı ayarlayın

Son olarak, `COWPY`'nin çıktısını `publish:` bloğuna ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Bu, iş akışı tarafından üretilen çıktıları kolayca görüntülememizi sağlayacaktır.

#### 1.3.4. İş akışını çalıştırın

`COWPY`'nin sağladığımız girdiler üzerinde çalıştığını doğrulamak için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Gördüğünüz gibi, `COWPY` her dosya için doğru karakteri kullanarak çalıştı.

??? abstract "Sonuçlar dizini içeriği"

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

??? example "results/cowpy-guten_tag.txt dosyasının içeriği"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Bu yaklaşım işe yarıyor; ancak bir sınırlaması var: kanalı iki ayrı alt kanala bölmek zorunda kaldık.
Sürece daha fazla alan geçirmek isteseydik, bunları daha fazla alt kanala bölmemiz gerekirdi.
Bu durum giderek zahmetli ve karmaşık bir hal alabilir.

İyi haber: bunu yapmanın daha basit bir yolu var.

### 1.4. Her şeyi sürece tek bir girdi olarak gruplamak

Alanları ayrı kanallara bölmek yerine, süreci tüm girdileri tek bir demet olarak alacak şekilde güncelleyebiliriz; bu da sürece yapılan çağrıyı basitleştirir.

#### 1.4.1. COWPY sürecini güncelleyin

`COWPY`'yi her satırdaki üç elemana karşılık gelen bir demet kabul edecek şekilde güncelleyin:

=== "Sonra"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy ile ASCII sanatı oluştur
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Önce"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // cowpy ile ASCII sanatı oluştur
    process COWPY {

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

Artık süreç, ihtiyaç duyabileceğimiz her şeyi içeren tek bir girdi alıyor.

#### 1.4.2. Girdi demetini oluşturmak için `map()` kullanın

Sürece demette geçirmek istediğimiz elemanları belirtmek için yine bir eşleme işlemi kullanmamız gerekiyor:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

`splitCsv`'den gelen Groovy map'ini olduğu gibi neden geçiremediğimizi merak edebilirsiniz.
Bunun nedeni, Nextflow'a kayıt dosyasının bir path olarak işlenmesi gerektiğini (yani düzgün şekilde hazırlanması gerektiğini) açıkça belirtmemiz gerektiğidir.
Bu, `COWPY`'nin girdi arayüzü düzeyinde gerçekleşir; burada `recording` elemanı açıkça `path` olarak tanımlanır.

#### 1.4.3. Süreç çağrısını güncelleyin

Son olarak, süreç çağrısındaki iki ayrı girdiyi az önce oluşturduğumuz tek demetle değiştirelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Bu, sürece yapılan çağrıyı biraz basitleştiriyor.

#### 1.4.4. İş akışını çalıştırın

`COWPY`'nin verileri hâlâ doğru şekilde işleyebildiğini doğrulamak için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

Çıktı, öncekiyle aynı yedi `cowpy-*.txt` dosyasıdır; artık `COWPY`'ye daha basit bir çağrıyla üretilmektedir.

??? abstract "Sonuçlar dizini içeriği"

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

??? example "results/cowpy-guten_tag.txt dosyasının içeriği"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Bu, `multiMap` yaklaşımına göre küçük bir iyileştirmedir.
Ancak girdi demetini oluşturmak için orijinal Groovy map'ini hâlâ açmak zorunda kaldık; üstelik süreç ile veri sayfası arasında sıkı bir bağlantı var: `COWPY` girdi tanımı artık `id`, `character` ve `recording` sütun adlarına doğrudan başvuruyor.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Bir iş arkadaşı farklı yapıda bir veri sayfası kullanırsa (ek sütunlar veya farklı sırada sütunlar içeren), bu süreç değiştirilmeden çalışmayacaktır.
Bu durum, sürecin girdi yapısının veri sayfasının tam bileşimine bağlı olması nedeniyle onu kırılgan kılar.

Bunu çözmek için, tam yapısını süreç arayüzüne sabit kodlamadan tüm metadata'yı bir paket olarak geçirmenin bir yolunu bulmamız gerekiyor.

### 1.5. Meta map + dosya arayüzü kullanma

Çözüm, kanaldaki iki ayrı kaygıyı birbirinden ayırmaktır: **bir örneğe ait metadata** ve **veri dosyasının** kendisi.
Tüm metadata'yı tek bir map'te (yani "meta map"te) toplayarak, veri sayfasının kaç metadata sütunu içerdiğinden bağımsız olarak tutarlı iki elemanlı bir demet elde ederiz:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Veri sayfasına sütun eklemek veya çıkarmak `meta`'nın içeriğini değiştirir; ancak `[meta, file]` demet yapısı sabit kalır.
Bu yapıyı kabul eden süreçlerin kaç metadata alanı bulunduğunu bilmesine gerek yoktur.

#### 1.5.1. Demet içeriğini meta map olarak yeniden düzenleyin

`map` işlemini `[meta, file]` demeti üretecek şekilde yeniden yapılandıralım:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Bir sonraki adımda güncellenecek

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Bir `view()` ifadesi eklediğimize, `COWPY` çağrısını yorum satırına aldığımıza ve süreç girdi tanımı henüz yeni yapıyla eşleşmediğinden `COWPY.out`'u `channel.empty()` ile değiştirdiğimize dikkat edin.

#### 1.5.2. Yeniden düzenlenmiş içeriği incelemek için iş akışını çalıştırın

Yeni kanal yapısını görmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Kanaldaki her eleman artık iki elemanlı bir demettir: önce meta map, ardından dosya.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Veri sayfasına daha sonra bir `language` sütunu eklersek, süreç girdi tanımında herhangi bir değişiklik yapmadan `meta.language` olarak erişilebilir hale gelecektir.

#### 1.5.3. `COWPY` sürecini meta map kullanacak şekilde güncelleyin

`COWPY`'yi `[meta, file]` demet yapısını kabul edecek şekilde güncelleyin:

=== "Sonra"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy ile ASCII sanatı oluştur
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Önce"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // cowpy ile ASCII sanatı oluştur
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Betik bloğunun içinde `meta.character`, meta map'teki `character` alanına erişir.
Meta map'teki herhangi bir alana aynı şekilde erişilebilir.

#### 1.5.4. Süreç çağrısını güncelleyin

`COWPY` çağrısını geri yükleyin ve çıktısını yayımlama için bağlayın:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Bir sonraki adımda güncellenecek

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

Çıktı yayımlamayı da geri yükledik.

#### 1.5.5. İş akışını çalıştırın

Her şeyin çalıştığını doğrulamak için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

Sonuçlar dizini artık ASCII sanatı dosyalarını içeriyor.

??? abstract "Dizin içeriği"

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

??? example "results/cowpy-guten_tag.txt dosyasının içeriği"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Süreç artık tüm metadata'yı `meta` aracılığıyla bir paket olarak alıyor, ihtiyaç duyduğunu (`meta.character`) kullanıyor ve geri kalanını yok sayıyor.

Bu, tüm [nf-core](https://nf-co.re/) modülleri tarafından kullanılan standart arayüzdür.
`tuple val(meta), path(file)` örüntüsü nf-core modül kütüphanesinde tutarlı biçimde yer alır; bu nedenle bu kuralı benimseyen iş akışları, nf-core modüllerini minimum sürtüşmeyle değiştirebilir.

### Özetle

Bu bölümde şunları öğrendiniz:

- **Veri sayfaları nasıl okunur:** Başlık bilgisiyle CSV dosyalarını ayrıştırmak için `splitCsv` kullanmak.
- **Meta map kuralı neden var:** Metadata'yı veri dosyalarından ayırarak `[meta, file]` demetleri oluşturmak, veri sayfası geliştikçe kanal yapısını kararlı tutar.
- **Meta map alanları bir sürecin içinde nasıl kullanılır:** Meta map'teki herhangi bir alana betik bloğunda nokta gösterimiyle erişilebilir.

---

## 2. Ek metadata düzenlemeleri

Meta map arayüzü yerli yerine oturduğuna göre, veriler pipeline boyunca akarken onu zenginleştirebiliriz.

Her kayıt dosyasındaki dili belirlemek için [`langid`](https://github.com/saffsd/langid.py) adlı bir araç kullanacağız.
Bir metin parçası verildiğinde, `stdout`'a bir dil tahmini ve olasılık skoru yazdırır.

### 2.1. Dil tanımlama adımı ekleme

`langid` aracını saran `IDENTIFY_LANGUAGE` adlı önceden yazılmış bir süreç modülü sağlıyoruz.

Kodunu incelemek için modül dosyasını açın:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

Girdi tanımı, 1. Bölümde oluşturduğumuz `tuple val(meta), path(file)` yapısını kullanıyor; bu nedenle `ch_datasheet` herhangi bir uyarlama yapılmadan doğrudan bu sürece beslenebilir.

Çıktı, üçüncü eleman olarak `stdout` ekler: bu, `langid`'in konsola yazdırdığı dil tahminini yakalar.
`sed` komutu olasılık skorunu ve sondaki yeni satır karakterini kaldırarak yalnızca iki harfli dil kodunu bırakır.

#### 2.1.1. `IDENTIFY_LANGUAGE`'a bir çağrı ekleyin

`IDENTIFY_LANGUAGE` süreç modülünü dahil edin ve veri sayfası kanalı üzerinde çağırın:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

Bu sürecin ana çıktısı yalnızca bir dizedir; dolayısıyla yayımlanacak çıktı dosyası yoktur.
Bunun yerine, işlemin sonuçlarını görüntülemek için `IDENTIFY_LANGUAGE.out.view()` kullanıyoruz.

#### 2.1.2. İş akışını çalıştırın

`COWPY` görevlerini yeniden çalıştırmaktan kaçınmak için `-resume` kullanarak dil tanımlamasını üretmek amacıyla iş akışını çalıştırın:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Artık veri setindeki her dosya için bir dil tahminine sahibiz.

Çıktı demetinin `[meta, file, lang_id]`'den oluştuğuna dikkat edin; yani meta map ve dosya, yeni sonucun yanı sıra taşınmaya devam ediyor.

!!! note "Not"

    Meta map'i sonuçlarla ilişkili tutma örüntüsü, sonuçları daha sonra kanallar arasında birleştirmeyi kolaylaştırır.
    Verileri doğru şekilde ilişkilendirmek için kanallardaki öğelerin sırasına güvenemezsiniz.
    Bunun yerine anahtar kullanmanız gerekir.
    Meta map'ler bu amaç için ideal bir yapı sağlar.

    Bu kullanım senaryosu [Splitting & Grouping](../splitting_and_grouping/index.md) yan görevinde ayrıntılı olarak incelenmektedir.

### 2.2. Metadata'yı süreç çıktılarıyla zenginleştirme

Dil tahmini, dosyadaki veriler hakkında bir tür metadata'dır.
Bunu ayrı bir eleman olarak tutmak yerine meta map'e katlayalım.

#### 2.2.1. Yeni ve genişletilmiş bir meta map oluşturun

Groovy `+` operatörünü kullanarak orijinalinin yerini alacak yeni bir meta map oluşturabiliriz:

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Bu işlemin özü `#!groovy meta + [lang: lang_id]`'dir.

Bu kod, dil kodunu içeren tek anahtar-değer çiftinden oluşan geçici bir map oluşturur (`[lang: lang_id]`), ardından Groovy `+` operatörünü kullanarak bunu önceden var olan metadata'yı içeren orijinal `meta` map'iyle birleştirir ve yeni, genişletilmiş bir meta map üretir.

Daha ayrıntılı bir açıklama için aşağıdaki kutuya bakın.

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

    Diyelim ki `map1`'den yeniden başlıyorsunuz, ancak dil tahmini kendi map'inde değil (`map2` yok).
    Bunun yerine, `lang_id` adlı bir değişkende tutuluyor ve değerini (`'fr'`) `lang` anahtarıyla saklamak istiyorsunuz.

    Aslında şunu yapabilirsiniz:

    ```groovy
    new_map = map1 + [lang: lang_id]
    ```

    Burada `[lang: lang_id]` anında yeni isimsiz bir map oluşturur ve `map1 + ` ifadesi `map1`'i yeni isimsiz map ile birleştirerek daha önce olduğu gibi aynı `new_map` içeriğini üretir.

    Güzel, değil mi?

    **Şimdi bunu bir Nextflow `channel.map()` işlemi bağlamına taşıyalım.**

    Kod şu hale gelir:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Bu şunları yapar:

    - `#!groovy map1, lang_id ->` demetteki iki öğeyi alır
    - `#!groovy map1 + [lang: lang_id]` yukarıda ayrıntılandırıldığı gibi yeni map'i oluşturur

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

    `map` işleminde `file`'ın neden oradan oraya taşınıyor gibi göründüğünü anlamakta güçlük çekiyorsanız, `#!groovy [meta + [lang: lang_id], file]` satırının `[new_map, file]` olarak okunduğunu hayal edin.
    Bu, `file`'ı demetteki ikinci konumunda bıraktığımızı daha net göstermelidir. Sadece `new_info` değerini alıp birinci konumdaki map'e katlıyoruz.

    **Ve bu bizi `tuple val(meta), path(file)` kanal yapısına geri getiriyor!**

#### 2.2.2. İş akışını çalıştırın

Kodun ne yaptığını anladığınızdan emin olduktan sonra, çalışıp çalışmadığını görmek için iş akışını çalıştırın:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "Meta map'ten anahtar kaldırma"

    Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) metodunu kullanarak meta map'ten bir anahtarı kaldırabilirsiniz; bu metot yalnızca belirttiğiniz anahtarları içeren yeni bir map döndürür:

    ```groovy
    meta.subMap(['id', 'character'])  // yalnızca 'id' ve 'character' içeren bir map döndürür
    ```

    Bu, bir aşağı akış süreci veya modülünün meta map'te biriken tüm alanlara ihtiyaç duymadığı durumlarda kullanışlıdır.

### 2.3. Koşullu ifadeler kullanarak dil grubu atama

Dil tahmini meta map'te yer aldığına göre, bundan daha fazla metadata türetebiliriz.
Veri setimizdeki diller iki aileye ayrılır: Cermen dilleri (İngilizce, Almanca) ve Roman dilleri (Fransızca, İspanyolca, İtalyanca).
`lang_group` alanı eklemek, bu sınıflandırmayı aşağı akışta kullanılabilir kılacaktır.

#### 2.3.1. Koşullu mantıkla bir `map` işlemi ekleyin

Dil ailesini atamak için koşullu mantık içeren ikinci bir `map` işlemi kullanacağız:

```groovy
.map { meta, file ->

    // lang_group'u tanımlayan koşullu mantık buraya gelir

    [meta + [lang_group: lang_group], file]
}
```

Uygulanacak mantık şu:

- Varsayılan değer olarak `lang_group = 'unknown'` ile başlayın.
- `meta.lang` `'de'` veya `'en'` ise `lang_group`'u `'germanic'` olarak ayarlayın.
- Aksi takdirde `meta.lang` `['fr', 'es', 'it']` listesinde yer alıyorsa `lang_group`'u `'romance'` olarak ayarlayın.

!!! tip "İpucu"

    `lang` değerine map işlemi içinde `meta.lang` ile erişebilirsiniz.

İş akışında aşağıdaki değişiklikleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Her selamlamanın dilini tanımlamak için langid'i çalıştır
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Temel noktalar şunlardır:

- `def lang_group = "unknown"` değişkeni güvenli bir varsayılan değerle başlatır.
- `if / else if` yapısı iki dil ailesini ele alır; diğer her şey `'unknown'` olarak kalır.
- `#!groovy .set { ch_languages }` elde edilen kanala bir sonraki adımda kullanılmak üzere bir isim verir.

<!-- TODO (gelecek) İlgili belgelere ek kaynaklar bölümünde not/bağlantı ekle -->

#### 2.3.2. İş akışını çalıştırın:

Çalıştığını doğrulamak için iş akışını çalıştırın:

```bash
nextflow run main.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Meta map artık dört alan taşıyor: `id`, `character`, `lang` ve `lang_group`.
Kanal yapısı hâlâ `[meta, file]`.

### 2.4. Çıktıları adlandırmak ve düzenlemek için metadata kullanma

`lang` ve `lang_group` artık meta map'te mevcut olduğuna göre, bunları çıktı dosya adlarına dil kodu eklemek ve dosyaları dil ailesine göre alt dizinlere düzenlemek için kullanabiliriz.

Bu üç değişiklik gerektirir: `COWPY` sürecini çıktısını yeniden adlandıracak ve `meta`'yı yaydığı şeye dahil edecek şekilde güncellemek, `COWPY` çağrısını `ch_languages` üzerinde çalışacak şekilde güncellemek ve alt dizin yolunu belirtmek için çıktı bloğunu güncellemek.

#### 2.4.1. `COWPY` sürecini güncelleyin

Meta map'teki dil kodunu kullanarak çıktı dosyasını yeniden adlandırın ve alt dizin yönlendirmesi için `lang_group`'a erişebilmesi amacıyla çıktıya `meta`'yı ekleyin:

=== "Sonra"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Önce"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Bu, girdi tanımını hiç değiştirmeden bir sürecin davranışını özelleştirmek için diğer metadata alanlarından nasıl yararlanabileceğimizi gösteriyor.

#### 2.4.2. `COWPY` çağrısını `ch_languages` üzerinde çalışacak şekilde güncelleyin

`COWPY(ch_datasheet)` yerine `COWPY(ch_languages)` kullanın:

=== "Sonra"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

Kanal içeriğini artık incelememize gerek olmadığından `ch_languages.view()` satırını da kaldırıyoruz.

#### 2.4.3. Çıktı bloğunu güncelleyin

Her dosyayı dil grubu alt dizinine yönlendirmek için `output {}` bloğuna bir `path` closure'ı ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Bu, çıktıları büyük bir esneklikle düzenlemek için metadata'yı nasıl kullanabileceğimizi gösteriyor.

#### 2.4.4. Tam pipeline'ı çalıştırın

Önceki sonuçları silin ve tam pipeline'ı çalıştırın:

```bash
rm -r results
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

Sonuçlar dizini artık dil ailesine göre düzenlenmiş olup her dosya tespit edilen dile göre adlandırılmıştır:

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

`output {}` bloğundaki `path` closure'ı her `[meta, file]` demetini alır ve alt dizin adı olarak `meta.lang_group`'u döndürür.
Dosya adının kendisi sürecin çıktısından gelir (`#!groovy "${meta.lang}-${input_file}"`).
Her iki metadata parçası da (dil kodu ve dil grubu) bu bölümde oluşturulan zenginleştirilmiş meta map'ten gelmektedir.

### Özetle

Bu bölümde şunları öğrendiniz:

- **Meta map'i süreç çıktılarıyla nasıl zenginleştirirsiniz:** `#!groovy meta + [anahtar: değer]` ile yeni anahtarlar eklemek, metadata'yı zenginleştirirken `[meta, file]` kanal yapısını korur.
- **Metadata'dan metadata nasıl türetilir:** `map` işlemi içindeki koşullu mantık, mevcut alanlardan yeni alanlar hesaplayabilir.
- **Çıktıları düzenlemek için metadata nasıl kullanılır:** `output {}` bloğundaki `path` closure'ı, dosyaları alt dizinlere yönlendirmek için meta map'ten okuyabilir.

---

## 3. Dayanıklılık değerlendirmeleri

Metadata değerleri süreç davranışını yönlendirdiğinde, eksik veya tamamlanmamış veriler teşhis edilmesi güç sorunlara yol açabilir.
Ne bekleyeceğinizi ve nasıl ele alacağınızı aşağıda bulabilirsiniz.

### 3.1. Zorunlu bir metadata alanı eksik olduğunda ne olur

`character` değeri, `COWPY` sürecinin geçerli bir sonuç üretmesi için zorunludur.
Hata modu, sütunun veri sayfasında mevcut olup olmadığına ve değerin boş mu yoksa tamamen yoksa mu olduğuna bağlıdır.

#### 3.1.1. Sütun var ama değer boş

Diyelim ki veri sayfasındaki bir girişin `character` alanı boş:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

Veri sayfası ayrıştırıldığında tüm girişler için `character` anahtarı oluşturulur; ancak `sampleA` için `meta.character` boş bir dize olacaktır.
Nextflow `#!groovy ${meta.character}` ifadesini komuta yerleştirdiğinde, `COWPY` aracı `-c` için boş bir argüman alır ve hata verir:

??? failure "Komut çıktısı"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

Hata mesajı (`expected one argument`), boş `-c` bayrağına işaret ediyor.
Work dizinindeki `.command.sh` dosyasını incelemek, komutun boş bir değerle çalıştırıldığını doğrular.

#### 3.1.2. Sütun veri sayfasında hiç yok

`character` sütunu tamamen yoksa:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

`character` anahtarı meta map'te hiç oluşturulmaz.
Süreç betiği `#!groovy ${meta.character}` ifadesini değerlendirdiğinde, eksik anahtar `null` döndürür ve Nextflow komuta tam anlamıyla `null` dizesini yerleştirir:

??? failure "Komut çıktısı"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

Çalıştırılan komuttaki `cowpy -c null` ifadesi, tanısal ipucudur.

### 3.2. Eksik metadata'yı ele alma stratejileri

İş akışlarını eksik metadata'ya karşı daha dayanıklı hale getirmek için birbirini tamamlayan iki yaklaşım vardır.

**1. Girdi doğrulaması**

En güvenilir çözüm, herhangi bir işlem başlamadan önce veri sayfasını doğrulamaktır; böylece sorunlar, çalışmanın ortasında anlaşılması güç bir süreç hatası olarak ortaya çıkmak yerine erken ve net bir hata mesajıyla yakalanır.
[Hello nf-core](../../hello_nf-core/05_input_validation.md) eğitimi, nf-schema eklentisini kullanarak girdi doğrulamasının nasıl ekleneceğini ele almaktadır. <!-- TODO (gelecek) uygun bir Doğrulama yan görevi bekliyor -->

**2. Zorunlu değerler için açık süreç girdileri**

Süreç arayüzünün belirli bir değerin zorunlu olduğunu iletmesini istiyorsanız, bunu meta map'ten açık bir girdi olarak çıkarmayı düşünebilirsiniz:

=== "Süreç tanımı"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "İş akışı çağrısı"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Bu yaklaşım, `character`'ı süreç sözleşmesinin görünür ve zorunlu bir parçası haline getirir.
Modülü okuyan herkes, bir karakter değerinin sağlanması gerektiğini hemen görebilir.
Alan yoksa, iş akışı süreç çalışmadan önce kanal düzeyinde açıkça hata verir.

Bu, yararlı bir tasarım ilkesini vurgular:

**Meta map'i isteğe bağlı veya açıklayıcı bilgiler için kullanın; zorunlu değerleri açık girdiler olarak çıkarın.**

Meta map, kanal yapılarını temiz ve kararlı tutar; ancak bir süreç tarafından gerçekten zorunlu olan değerler için bunları adlandırılmış girdiler olarak ortaya çıkarmak, netliği artırır ve modülün diğer bağlamlarda doğru şekilde kullanılmasını kolaylaştırır.

### Özetle

Bu bölümde şunları gördünüz:

- **Eksik metadata nasıl kendini gösterir:** Boş bir alan boş bir argüman üretir; olmayan bir alan komuta tam anlamıyla `null` olarak yerleştirilir.
- **İki tamamlayıcı strateji:** Sorunları erken yakalamak için girdi doğrulaması ve gereksinimleri açıkça iletmek için açık süreç girdileri.

---

## Özet

Bu yan görevde, Nextflow iş akışlarında metadata ile etkili bir şekilde çalışmayı keşfettiniz.

"Meta map + veri dosyası" demet örüntüsü, Nextflow'da temel bir kuraldır ve metadata'yı tek tek değerler olarak geçirmeye kıyasla çeşitli avantajlar sunar:

- Veri sayfası geliştikçe kanal yapısı kararlı kalır
- Süreç davranışı, alan adlarını sabit kodlamadan örnek başına özelleştirilebilir
- Metadata, çıktıları adlandırmak, gruplamak ve düzenlemek için pipeline boyunca kullanılabilir
- Bu arayüze göre yazılmış modüller, nf-core modülleri dahil, birbirleriyle değiştirilebilir

### Temel örüntüler

1.  **Metadata'yı okuma ve yapılandırma:** Bir CSV veri sayfasını ayrıştırın ve meta map oluşturun.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **İş akışı sırasında metadata'yı genişletme:** Süreç çıktılarından veya türetilmiş mantıktan yeni anahtarlar ekleyin.

    ```groovy
    // Süreç çıktısından
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // Koşullu mantıktan
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Bir sürecin içinde metadata kullanma:** Betik bloğunda nokta gösterimiyle herhangi bir alana erişin.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Çıktıları metadata değerine göre düzenleme:** `output {}` bloğunda `path` closure'ı kullanın.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Ek kaynaklar

- [map operatörü](https://www.nextflow.io/docs/latest/operator.html#map)
- [multiMap operatörü](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [stdout çıktı niteleyicisi](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Sırada ne var?

[Yan Görevler menüsüne](../index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
