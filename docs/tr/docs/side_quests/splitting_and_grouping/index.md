# Bölme ve Gruplama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, verilerle esnek bir şekilde çalışmak için güçlü araçlar sunar. Temel yeteneklerden biri, verileri farklı akışlara bölmek ve ardından ilgili öğeleri yeniden bir araya getirmektir. Bu özellik, analiz için sonuçları birleştirmeden önce farklı örnek türlerini ayrı ayrı işlemeniz gereken biyoinformatik iş akışlarında özellikle değerlidir.

Bunu posta sıralama işlemine benzetebilirsiniz: mektupları hedefe göre ayırır, her yığını farklı şekilde işler, ardından aynı kişiye gidecek öğeleri yeniden birleştirirsiniz. Nextflow, bilimsel verilerle bunu gerçekleştirmek için özel operatörler kullanır. Bu yaklaşım, dağıtık hesaplama ve biyoinformatik iş akışlarında **scatter/gather** deseni olarak da yaygın biçimde bilinir.

Nextflow'un kanal sistemi bu esnekliğin merkezinde yer alır. Kanallar, iş akışınızın farklı bölümlerini birbirine bağlayarak verinin analiziniz boyunca akmasını sağlar. Tek bir veri kaynağından birden fazla kanal oluşturabilir, her kanalı farklı şekilde işleyebilir ve gerektiğinde kanalları yeniden birleştirebilirsiniz. Bu yaklaşım, karmaşık biyoinformatik analizlerin dallanma ve birleşme yollarını doğal olarak yansıtan iş akışları tasarlamanıza olanak tanır.

### Öğrenme hedefleri

Bu yan görevde, Nextflow'un kanal operatörlerini kullanarak verileri bölmeyi ve gruplamayı öğreneceksiniz.
Örnek bilgileri ve ilişkili veri dosyaları içeren bir CSV dosyasıyla başlayıp bu verileri işleyecek ve yeniden düzenleyeceksiniz.

Bu yan görevin sonunda, aşağıdaki teknikleri kullanarak veri akışlarını etkili biçimde ayırabilecek ve birleştirebileceksiniz:

- `splitCsv` kullanarak dosyalardan veri okuma
- `filter` ve `map` ile verileri filtreleme ve dönüştürme
- `join` ve `groupTuple` kullanarak ilgili verileri birleştirme
- Paralel işleme için `combine` ile veri kombinasyonları oluşturma
- `subMap` ve tekilleştirme stratejileri kullanarak veri yapısını optimize etme
- Kanal yapılarını işlemenize yardımcı olmak için adlandırılmış closure'larla yeniden kullanılabilir fonksiyonlar oluşturma

Bu beceriler, temiz ve sürdürülebilir bir kod yapısını korurken birden fazla girdi dosyasını ve farklı veri türlerini verimli şekilde işleyebilen iş akışları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, dosyalarla çalışma, meta veri) rahatça kullanabiliyor olmalısınız.

**İsteğe bağlı:** Önce [İş akışlarında meta veri](../metadata/) yan görevini tamamlamanızı öneririz.
Bu görev, `splitCsv` ile CSV dosyalarını okuma ve burada yoğun biçimde kullanacağımız meta map'leri oluşturma konularının temellerini kapsar.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı şekilde eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitime ait dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/splitting_and_grouping
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

Bir ana iş akışı dosyası ve `samplesheet.csv` adlı bir örnek sayfası içeren `data` dizini bulacaksınız.

```console title="Directory contents"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Örnek sayfası, farklı hastalardan alınan örneklere ilişkin bilgileri içermektedir; hasta kimliği, örnek tekrar numarası, tür (normal veya tümör) ve varsayımsal veri dosyalarının yolları (aslında mevcut değiller, ancak varmış gibi davranacağız) yer almaktadır.

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Bu örnek sayfası, üç hastadan (A, B, C) alınan sekiz örneği listelemektedir.

Her hasta için `tumor` türünde (genellikle tümör biyopsilerinden elde edilen) veya `normal` türünde (sağlıklı doku ya da kandan alınan) örnekler mevcuttur.
Kanser analiziyle tanışık değilseniz, bunun eşleştirilmiş tümör/normal örnekler kullanarak karşılaştırmalı analizler gerçekleştiren deneysel bir modele karşılık geldiğini bilmeniz yeterlidir.

Özellikle A hastası için iki teknik replikat (tekrar) setimiz bulunmaktadır.

!!! note "Not"

    Bu deneysel tasarıma aşina değilseniz endişelenmeyin; bu eğitimi anlamak için kritik değildir.

#### Görevi inceleyin

Göreviniz, aşağıdakileri gerçekleştirecek bir Nextflow iş akışı yazmaktır:

1. **Okuma**: CSV dosyasından örnek verilerini okuyun ve meta map'lerle yapılandırın
2. **Ayırma**: Örnekleri türe göre (normal ve tümör) farklı kanallara ayırın
3. **Birleştirme**: Eşleştirilmiş tümör/normal çiftlerini hasta kimliği ve replikat numarasına göre birleştirin
4. **Dağıtma**: Paralel işleme için örnekleri genomik aralıklara dağıtın
5. **Gruplama**: İlgili örnekleri aşağı akış analizi için yeniden bir araya getirin

Bu, bağımsız işleme için verileri bölmeniz ve ardından karşılaştırmalı analiz için ilgili öğeleri yeniden birleştirmeniz gereken yaygın bir biyoinformatik desenini temsil etmektedir.

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebildiyseniz başlamaya hazırsınız.

---

## 1. Örnek verileri okuma

### 1.1. `splitCsv` ile örnek verileri okuma ve meta map'ler oluşturma

`splitCsv` ile örnek verileri okuyarak meta map deseninde düzenleyerek başlayalım. `main.nf` dosyasında iş akışına zaten başlandığını göreceksiniz.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Not"

    Bu eğitim boyunca, tüm kanal değişkenleri için Nextflow kanalları olduklarını açıkça belirtmek amacıyla `ch_` önekini kullanacağız.

[İş akışlarında meta veri](../metadata/) yan görevini tamamladıysanız bu deseni tanıyacaksınız. CSV'yi okumak için `splitCsv` kullanacak ve meta veriyi dosya yollarından ayırmak amacıyla verileri hemen bir meta map ile yapılandıracağız.

!!! info "Bilgi"

    Bu eğitimde `map` adı verilen iki farklı kavramla karşılaşacağız:

    - **Veri yapısı**: Anahtar-değer çiftlerini depolayan Groovy map'i (diğer dillerdeki sözlük/hash yapılarına eşdeğer)
    - **Kanal operatörü**: Bir kanaldaki öğeleri dönüştüren `.map()` operatörü

    Hangisini kastettiğimizi bağlamda açıklayacağız; ancak Nextflow ile çalışırken bu ayrımı anlamak önemlidir.

`main.nf` dosyasına şu değişiklikleri uygulayın:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Bu, `splitCsv` işlemini (CSV'yi başlıklarla okuma) ve `map` işlemini (verileri `[meta, dosya]` demetleri olarak yapılandırma) tek adımda birleştirir. Bu değişikliği uygulayın ve pipeline'ı çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Artık her öğenin bir `[meta, dosya]` demeti olduğu bir kanalımız var; meta veri, dosya yollarından ayrılmış durumda. Bu yapı, iş yükümüzü meta veri alanlarına göre bölmemize ve gruplamamıza olanak tanır.

---

## 2. Verileri filtreleme ve dönüştürme

### 2.1. `filter` ile verileri filtreleme

Verileri bir koşula göre filtrelemek için [`filter` operatörünü](https://www.nextflow.io/docs/latest/operator.html#filter) kullanabiliriz. Yalnızca normal örnekleri işlemek istediğimizi varsayalım. Bunu, verileri `type` alanına göre filtreleyerek yapabiliriz. Bunu `view` operatöründen önce ekleyelim.

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Filtrelenmiş sonucu görmek için iş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Verileri yalnızca normal örnekleri içerecek şekilde başarıyla filtreledik. Bunun nasıl çalıştığını özetleyelim.

`filter` operatörü, kanaldaki her öğeye uygulanan bir closure alır. Closure `true` döndürürse öğe dahil edilir; `false` döndürürse öğe dışlanır.

Bizim durumumuzda yalnızca `meta.type == 'normal'` olan örnekleri tutmak istiyoruz. Closure, her örneğe atıfta bulunmak için `meta,file` demetini kullanır, `meta.type` ile örnek türüne erişir ve bunun `'normal'`'e eşit olup olmadığını kontrol eder.

Bu, yukarıda tanıttığımız tek closure ile gerçekleştirilir:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Ayrı filtrelenmiş kanallar oluşturma

Şu anda filtreyi doğrudan CSV'den oluşturulan kanala uyguluyoruz; ancak bunu birden fazla şekilde filtrelemek istiyoruz, bu nedenle normal örnekler için ayrı bir filtrelenmiş kanal oluşturmak üzere mantığı yeniden yazalım:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Sonuçları görmek için pipeline'ı çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Verileri başarıyla filtreledik ve normal örnekler için ayrı bir kanal oluşturduk.

Tümör örnekleri için de filtrelenmiş bir kanal oluşturalım:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Normal ve tümör örneklerini iki farklı kanala ayırdık ve çıktıda farklı şekilde etiketlemek için `view()` fonksiyonuna bir closure sağladık: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Özetle

Bu bölümde şunları öğrendiniz:

- **Verileri filtreleme**: `filter` ile verilerin nasıl filtreleneceği
- **Verileri bölme**: Verilerin bir koşula göre farklı kanallara nasıl bölüneceği
- **Verileri görüntüleme**: Verileri yazdırmak ve farklı kanallardan gelen çıktıları etiketlemek için `view`'ın nasıl kullanılacağı

Normal ve tümör örneklerini iki farklı kanala ayırdık. Sonraki adımda, normal ve tümör örneklerini `id` alanına göre birleştireceğiz.

---

## 3. Kanalları tanımlayıcılara göre birleştirme

Önceki bölümde normal ve tümör örneklerini iki farklı kanala ayırdık. Bunlar, türlerine göre belirli süreçler veya iş akışları kullanılarak bağımsız olarak işlenebilir. Peki aynı hastadan alınan normal ve tümör örneklerini karşılaştırmak istediğimizde ne olur? Bu noktada, örnekleri `id` alanına göre eşleştirerek yeniden birleştirmemiz gerekir.

Nextflow, kanalları birleştirmek için birçok yöntem içerir; ancak bu durumda en uygun operatör [`join`](https://www.nextflow.io/docs/latest/operator.html#join)'dir. SQL'e aşinaysanız, birleştirme anahtarını ve gerçekleştirilecek birleştirme türünü belirttiğimiz `JOIN` işlemi gibi davrandığını söyleyebiliriz.

### 3.1. Hasta kimliğine göre birleştirmek için `map` ve `join` kullanma

[`join`](https://www.nextflow.io/docs/latest/operator.html#join) belgelerine bakarsak, varsayılan olarak iki kanalı her demetin ilk öğesine göre birleştirdiğini görebiliriz.

#### 3.1.1. Veri yapısını kontrol etme

Konsol çıktısı hâlâ mevcut değilse, veri yapımızı kontrol etmek ve `id` alanına göre birleştirmek için nasıl değiştirmemiz gerektiğini görmek amacıyla pipeline'ı çalıştıralım.

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

`id` alanının her meta map'teki ilk öğe olduğunu görebiliriz. `join`'in çalışması için her demetteki `id` alanını yalıtmamız gerekir. Bundan sonra, iki kanalı birleştirmek için `join` operatörünü kullanabiliriz.

#### 3.1.2. `id` alanını yalıtma

`id` alanını yalıtmak için, `id` alanını ilk öğe olarak içeren yeni bir demet oluşturmak amacıyla [`map` operatörünü](https://www.nextflow.io/docs/latest/operator.html#map) kullanabiliriz.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Fark ince olabilir, ancak her demetin ilk öğesinin `id` alanı olduğunu görebilmelisiniz.

#### 3.1.3. İki kanalı birleştirme

Artık `id` alanına göre iki kanalı birleştirmek için `join` operatörünü kullanabiliriz.

Yine, birleştirilmiş çıktıları yazdırmak için `view` kullanacağız.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Çok geniş olduğu için okumak biraz zor olabilir; ancak örneklerin `id` alanına göre birleştirildiğini görebilmelisiniz. Her demet artık şu biçimdedir:

- `id`: Örnek kimliği
- `normal_meta_map`: Tür, replikat ve BAM dosyasının yolunu içeren normal örnek meta verisi
- `normal_sample_file`: Normal örnek dosyası
- `tumor_meta_map`: Tür, replikat ve BAM dosyasının yolunu içeren tümör örnek meta verisi
- `tumor_sample`: Tür, replikat ve BAM dosyasının yolunu içeren tümör örneği

!!! warning "Uyarı"

    `join` operatörü, eşleşmeyen demetleri atacaktır. Bu örnekte tüm örneklerin tümör ve normal için eşleştirildiğinden emin olduk; ancak bu durum geçerli değilse eşleşmeyen demetleri korumak için `remainder: true` parametresini kullanmanız gerekir. Daha fazla ayrıntı için [belgelere](https://www.nextflow.io/docs/latest/operator.html#join) bakın.

Artık bir demetteki bir alanı yalıtmak için `map`'i ve demetleri ilk alana göre birleştirmek için `join`'i nasıl kullanacağınızı biliyorsunuz.
Bu bilgiyle, kanalları paylaşılan bir alana göre başarıyla birleştirebilirsiniz.

Sonraki adımda, birden fazla alana göre birleştirmek istediğiniz durumu ele alacağız.

### 3.2. Birden fazla alana göre birleştirme

sampleA için 2 replikatımız var, ancak sampleB ve sampleC için yalnızca 1 replikat bulunuyor. Bu durumda `id` alanını kullanarak bunları etkili biçimde birleştirebildik; peki senkronizasyon bozulursa ne olur? Farklı replikatlardan normal ve tümör örneklerini karıştırabiliriz!

Bunu önlemek için birden fazla alana göre birleştirebiliriz. Bunu başarmanın aslında birden fazla yolu vardır; ancak hem örnek `id`'sini hem de `replicate` numarasını içeren yeni bir birleştirme anahtarı oluşturmaya odaklanacağız.

Yeni bir birleştirme anahtarı oluşturarak başlayalım. Bunu, `id` ve `repeat` alanlarını ilk öğe olarak içeren yeni bir demet oluşturmak için [`map` operatörünü](https://www.nextflow.io/docs/latest/operator.html#map) kullanarak daha önce yaptığımız gibi yapabiliriz.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Artık birleştirmenin hem `id` hem de `repeat` alanlarını kullanarak gerçekleştiğini görmeliyiz. İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Her birleştirilmiş sonucun ilk öğesi olarak iki elemanlı bir demetin (`id` ve `repeat` alanları) nasıl yer aldığına dikkat edin. Bu, karmaşık öğelerin birleştirme anahtarı olarak kullanılabileceğini ve aynı koşullardan gelen örnekler arasında oldukça karmaşık eşleştirmeler yapılabileceğini göstermektedir.

Farklı anahtarlarla birleştirmenin daha fazla yolunu keşfetmek istiyorsanız, ek seçenekler ve örnekler için [join operatörü belgelerine](https://www.nextflow.io/docs/latest/operator.html#join) bakın.

### 3.3. Yeni bir birleştirme anahtarı oluşturmak için `subMap` kullanma

Önceki yaklaşım, birleştirme anahtarımızdaki alan adlarını kaybettirir; `id` ve `repeat` alanları yalnızca bir değer listesine dönüşür. Alan adlarını daha sonra erişmek üzere korumak için [`subMap` metodunu](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) kullanabiliriz.

`subMap` metodu, bir map'ten yalnızca belirtilen anahtar-değer çiftlerini çıkarır. Burada birleştirme anahtarımızı oluşturmak için yalnızca `id` ve `repeat` alanlarını çıkaracağız.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Artık yalnızca `id` ve `repeat` alanlarını içermekle kalmayıp alan adlarını da koruyan yeni bir birleştirme anahtarımız var; böylece bunlara daha sonra `meta.id` ve `meta.repeat` gibi isimle erişebiliriz.

### 3.4. Map'te adlandırılmış closure kullanma

Tekrarı önlemek ve hataları azaltmak için adlandırılmış bir closure kullanabiliriz. Adlandırılmış closure, birden fazla yerde çağırabileceğimiz yeniden kullanılabilir bir fonksiyon oluşturmamıza olanak tanır.

Bunu yapmak için önce closure'ı yeni bir değişken olarak tanımlıyoruz:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Map dönüşümünü yeniden kullanabileceğimiz adlandırılmış bir değişken olarak tanımladık.

Ayrıca, bu kanalı alan herhangi bir sürecin dosyayı doğru şekilde işleyebilmesi için `file()` kullanarak dosya yolunu bir Path nesnesine dönüştürdüğümüze dikkat edin (daha fazla bilgi için [Dosyalarla çalışma](../working_with_files/) bölümüne bakın).

Closure'ı iş akışımızda uygulayalım:

=== "Sonra"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Önce"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "Not"

    `map` operatörü, closure'ı argüman olarak iletmek için `{ }` yerine `( )` kullanmaya geçti. Bunun nedeni, `map` operatörünün argüman olarak bir closure beklemesi ve `{ }` ifadesinin anonim bir closure tanımlamak için kullanılmasıdır. Adlandırılmış bir closure çağırırken `( )` sözdizimini kullanın.

Her şeyin hâlâ çalıştığını doğrulamak için iş akışını bir kez daha çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Adlandırılmış closure kullanmak, aynı dönüşümü birden fazla yerde yeniden kullanmamıza olanak tanır; bu da hata riskini azaltır ve kodu daha okunabilir ve sürdürülebilir hale getirir.

### 3.5. Veri tekrarını azaltma

İş akışımızda çok fazla tekrarlanan veri var. Birleştirilmiş örneklerdeki her öğe, `id` ve `repeat` alanlarını tekrarlıyor. Bu bilgi zaten gruplama anahtarında mevcut olduğundan bu fazlalıktan kaçınabiliriz. Hatırlatmak gerekirse, mevcut veri yapımız şöyle görünüyor:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

`id` ve `repeat` alanları gruplama anahtarında mevcut olduğundan, tekrarı önlemek için bunları her kanal öğesinin geri kalanından kaldıralım. Bunu, yalnızca `type` alanını içeren yeni bir map oluşturmak için `subMap` metodunu kullanarak yapabiliriz. Bu yaklaşım, veri yapımızdaki fazlalığı ortadan kaldırırken gerekli tüm bilgileri korumamıza olanak tanır.

=== "Sonra"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Artık closure, ilk öğenin `id` ve `repeat` alanlarını, ikinci öğenin ise yalnızca `type` alanını içerdiği bir demet döndürüyor. Bu, `id` ve `repeat` bilgilerini gruplama anahtarında bir kez depolayarak fazlalığı ortadan kaldırırken gerekli tüm bilgileri korur.

Bunun nasıl göründüğünü görmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

`id` ve `repeat` alanlarını gruplama anahtarında yalnızca bir kez belirttiğimizi ve örnek verisinde `type` alanına sahip olduğumuzu görebiliriz. Hiçbir bilgiyi kaybetmedik; ancak kanal içeriklerimizi daha özlü hale getirmeyi başardık.

### 3.6. Gereksiz bilgileri kaldırma

Yukarıda tekrarlanan bilgileri kaldırdık; ancak kanallarımızda hâlâ bazı gereksiz bilgiler var.

Başlangıçta normal ve tümör örneklerini `filter` kullanarak ayırdık, ardından bunları `id` ve `repeat` anahtarlarına göre birleştirdik. `join` operatörü, demetlerin birleştirilme sırasını korur; bu nedenle bizim durumumuzda, sol tarafta normal örnekler ve sağ tarafta tümör örnekleriyle, elde edilen kanal şu yapıyı korur: `id, <normal öğeler>, <tümör öğeleri>`.

Kanaldaki her öğenin konumunu bildiğimizden, `[type:normal]` ve `[type:tumor]` meta verilerini düşürerek yapıyı daha da basitleştirebiliriz.

=== "Sonra"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Sonucu görmek için tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Özetle

Bu bölümde şunları öğrendiniz:

- **Demetleri işleme**: Bir demetteki bir alanı yalıtmak için `map`'in nasıl kullanılacağı
- **Demetleri birleştirme**: Demetleri ilk alana göre birleştirmek için `join`'in nasıl kullanılacağı
- **Birleştirme anahtarları oluşturma**: Yeni bir birleştirme anahtarı oluşturmak için `subMap`'in nasıl kullanılacağı
- **Adlandırılmış closure'lar**: Map'te adlandırılmış closure'ın nasıl kullanılacağı
- **Birden fazla alana göre birleştirme**: Daha hassas eşleştirme için birden fazla alana göre nasıl birleştirileceği
- **Veri yapısı optimizasyonu**: Gereksiz bilgileri kaldırarak kanal yapısının nasıl düzenleneceği

Artık bir örnek sayfasını bölebilen, normal ve tümör örneklerini filtreleyebilen, bunları örnek kimliği ve replikat numarasına göre birleştirebilen ve sonuçları yazdırabilen bir iş akışınız var.

Bu, bağımsız olarak işledikten sonra örnekleri veya diğer veri türlerini eşleştirmeniz gereken biyoinformatik iş akışlarında yaygın bir desendir; dolayısıyla kullanışlı bir beceridir. Sonraki adımda, bir örneği birden fazla kez tekrarlamayı ele alacağız.

## 4. Örnekleri aralıklara dağıtma

Biyoinformatik iş akışlarındaki temel bir desen, analizleri genomik bölgelere dağıtmaktır. Örneğin, varyant çağırma işlemi genomu aralıklara (kromozomlar veya daha küçük bölgeler gibi) bölerek paralelleştirilebilir. Bu paralelleştirme stratejisi, hesaplama yükünü birden fazla çekirdek veya düğüme dağıtarak genel yürütme süresini önemli ölçüde azaltır ve pipeline verimliliğini artırır.

Aşağıdaki bölümde, örnek verilerimizi birden fazla genomik aralığa nasıl dağıtacağımızı göstereceğiz. Her örneği her aralıkla eşleştirerek farklı genomik bölgelerin paralel işlenmesine olanak tanıyacağız. Bu, veri kümemizin boyutunu aralık sayısıyla çarparak daha sonra bir araya getirilebilecek birden fazla bağımsız analiz birimi oluşturacaktır.

### 4.1. `combine` kullanarak örnekleri aralıklara yayma

Bir aralık kanalı oluşturarak başlayalım. İşleri basit tutmak için manuel olarak tanımlayacağımız 3 aralık kullanacağız. Gerçek bir iş akışında bunları bir dosya girdisinden okuyabilir veya hatta çok sayıda aralık dosyası içeren bir kanal oluşturabilirsiniz.

=== "Sonra"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Şimdi, her örneği her aralık için tekrarlamak istediğimizi hatırlayın. Bu, örnekler ve aralıkların Kartezyen çarpımı olarak da adlandırılır. Bunu [`combine` operatörünü](https://www.nextflow.io/docs/latest/operator.html#combine) kullanarak başarabiliriz. Bu operatör, kanal 1'deki her öğeyi alır ve kanal 2'deki her öğe için tekrarlar. İş akışımıza bir combine operatörü ekleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Şimdi çalıştıralım ve ne olduğunu görelim:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Başarılı! 3 aralıklı listemizde her örneği her aralık için tekrarladık. Kanaldaki öğe sayısını etkili biçimde üçe katladık.

Okumak biraz zor; bu nedenle bir sonraki bölümde düzenleyeceğiz.

### 4.2. Kanalı düzenleme

Örnek verilerimizi daha anlaşılır hale getirmek için `map` operatörünü kullanabiliriz. Aralık dizesini ilk öğedeki birleştirme map'ine taşıyalım.

=== "Sonra"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Bu map işleminin adım adım ne yaptığını açıklayalım.

Önce, kodu daha okunabilir kılmak için adlandırılmış parametreler kullanıyoruz. `grouping_key`, `normal`, `tumor` ve `interval` adlarını kullanarak demetteki öğelere indeks yerine isimle başvurabiliyoruz:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Ardından, `grouping_key`'i `interval` alanıyla birleştiriyoruz. `grouping_key`, `id` ve `repeat` alanlarını içeren bir map'tir. Groovy'nin map toplama işlemini (`+`) kullanarak `interval` ile yeni bir map oluşturuyoruz:

```groovy
                grouping_key + [interval: interval],
```

Son olarak, bunu üç öğeli bir demet olarak döndürüyoruz: birleştirilmiş meta veri map'i, normal örnek dosyası ve tümör örnek dosyası:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Tekrar çalıştıralım ve kanal içeriklerini kontrol edelim:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Verilerinizi doğru yapıya dönüştürmek için `map` kullanmak zor olabilir; ancak etkili veri işleme için kritik öneme sahiptir.

Artık her örneği tüm genomik aralıklara yayarak paralel olarak işlenebilecek birden fazla bağımsız analiz birimi oluşturduk. Peki ilgili örnekleri yeniden bir araya getirmek istersek? Bir sonraki bölümde, ortak öznitelikleri paylaşan örnekleri nasıl gruplayacağımızı öğreneceğiz.

### Özetle

Bu bölümde şunları öğrendiniz:

- **Örnekleri aralıklara yayma**: Örnekleri aralıklar üzerinde tekrarlamak için `combine`'ın nasıl kullanılacağı
- **Kartezyen çarpım oluşturma**: Örneklerin ve aralıkların tüm kombinasyonlarının nasıl üretileceği
- **Kanal yapısını düzenleme**: Daha iyi okunabilirlik için verileri yeniden yapılandırmak amacıyla `map`'in nasıl kullanılacağı
- **Paralel işleme hazırlığı**: Dağıtık analiz için verilerin nasıl hazırlanacağı

## 5. `groupTuple` kullanarak örnekleri toplama

Önceki bölümlerde, bir girdi dosyasından verileri nasıl böleceğimizi ve belirli alanlara göre (bizim durumumuzda normal ve tümör örnekleri) nasıl filtreleyeceğimizi öğrendik. Ancak bu yalnızca tek bir birleştirme türünü kapsar. Örnekleri belirli bir özniteliğe göre gruplamak istersek ne olur? Örneğin, eşleştirilmiş normal-tümör çiftlerini birleştirmek yerine, türlerinden bağımsız olarak "sampleA"dan gelen tüm örnekleri birlikte işlemek isteyebiliriz. Bu desen, sonunda sonuçları karşılaştırmadan veya birleştirmeden önce verimlilik nedeniyle ilgili örnekleri ayrı ayrı işlemek isteyebileceğiniz biyoinformatik iş akışlarında yaygındır.

Nextflow bunu yapmak için yerleşik yöntemler içerir; inceleyeceğimiz ana yöntem `groupTuple`'dır.

Aynı `id` ve `interval` alanlarına sahip tüm örneklerimizi gruplayarak başlayalım; bu, teknik replikatları gruplamak ancak anlamlı biçimde farklı örnekleri ayrı tutmak istediğimiz bir analize tipik olarak karşılık gelir.

Bunu yapmak için gruplama değişkenlerimizi yalıtmamız gerekir; böylece bunları tek başına kullanabiliriz.

İlk adım, önceki bölümde yaptığımıza benzer. Gruplama değişkenimizi demetin ilk öğesi olarak yalıtmamız gerekir. Hatırlayın, ilk öğemiz şu anda `id`, `repeat` ve `interval` alanlarından oluşan bir map'tir:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Map'ten `id` ve `interval` alanlarımızı yalıtmak için daha önce kullandığımız `subMap` metodunu yeniden kullanabiliriz. Daha önce olduğu gibi, her örnek için demetin ilk öğesine `subMap` metodunu uygulamak amacıyla `map` operatörünü kullanacağız.

=== "Sonra"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Tekrar çalıştıralım ve kanal içeriklerini kontrol edelim:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

`id` ve `interval` alanlarını başarıyla yalıttığımızı, ancak örnekleri henüz gruplamadığımızı görebiliriz.

!!! note "Not"

    Burada `replicate` alanını atıyoruz. Bunun nedeni, aşağı akış işleme için buna ihtiyaç duymamamızdır. Bu eğitimi tamamladıktan sonra, sonraki gruplamayı etkilemeden bunu dahil edip edemeyeceğinizi deneyin!

Şimdi [`groupTuple` operatörünü](https://www.nextflow.io/docs/latest/operator.html#grouptuple) kullanarak örnekleri bu yeni gruplama öğesine göre gruplayalım.

=== "Sonra"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

Hepsi bu kadar! Tek bir satır kod ekledik. Çalıştırdığımızda ne olduğunu görelim:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Verilerimizin yapısının değiştiğine ve her kanal öğesinde dosyaların artık `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` gibi demetler içinde yer aldığına dikkat edin. Bunun nedeni, `groupTuple` kullandığımızda Nextflow'un bir grubun her örneği için tek dosyaları birleştirmesidir. Verileri aşağı akışta işlemeye çalışırken bunu hatırlamak önemlidir.

!!! note "Not"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose), groupTuple'ın tersidir. Bir kanaldaki öğeleri açar ve düzleştirir. `transpose` ekleyerek yukarıda gerçekleştirdiğimiz gruplamayı geri almayı deneyin!

### Özetle

Bu bölümde şunları öğrendiniz:

- **İlgili örnekleri gruplama**: Örnekleri ortak özniteliklere göre toplamak için `groupTuple`'ın nasıl kullanılacağı
- **Gruplama anahtarlarını yalıtma**: Gruplama için belirli alanları çıkarmak amacıyla `subMap`'in nasıl kullanılacağı
- **Gruplandırılmış veri yapılarını işleme**: `groupTuple` tarafından oluşturulan iç içe yapıyla nasıl çalışılacağı
- **Teknik replikat işleme**: Aynı deneysel koşulları paylaşan örneklerin nasıl gruplandırılacağı

---

## Özet

Bu yan görevde, kanalları kullanarak verileri bölmeyi ve gruplamayı öğrendiniz.

Veriler pipeline boyunca akarken bunları değiştirerek, döngüler veya while ifadeleri kullanmadan ölçeklenebilir bir pipeline oluşturabilirsiniz; bu da daha geleneksel yaklaşımlara göre çeşitli avantajlar sunar:

- Ek kod yazmadan istediğimiz kadar az veya çok girdiyle ölçeklenebiliriz
- Yineleme yerine pipeline boyunca veri akışını yönetmeye odaklanırız
- Gerektiği kadar karmaşık veya basit olabiliriz
- Pipeline daha bildirimsel hale gelir; nasıl yapılması gerektiği yerine ne yapılması gerektiğine odaklanır
- Nextflow, bağımsız işlemleri paralel olarak çalıştırarak yürütmeyi bizim için optimize eder

Bu kanal işlemlerinde ustalaşmak, döngülere veya yinelemeli programlamaya başvurmadan karmaşık veri ilişkilerini işleyebilen esnek ve ölçeklenebilir pipeline'lar oluşturmanıza olanak tanır; Nextflow'un yürütmeyi optimize etmesine ve bağımsız işlemleri otomatik olarak paralelleştirmesine izin verir.

### Temel desenler

1.  **Yapılandırılmış girdi verisi oluşturma:** Meta map'lerle bir CSV dosyasından başlama ([İş akışlarında meta veri](../metadata/) bölümündeki desenleri temel alarak)

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Verileri ayrı kanallara bölme:** `type` alanına göre verileri bağımsız akışlara bölmek için `filter` kullandık

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Eşleştirilmiş örnekleri birleştirme:** `id` ve `repeat` alanlarına göre ilgili örnekleri yeniden birleştirmek için `join` kullandık

    - İki kanalı anahtara göre birleştirme (demetin ilk öğesi)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Birleştirme anahtarını çıkarma ve bu değere göre birleştirme

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - subMap kullanarak birden fazla alana göre birleştirme

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Aralıklara dağıtma:** Paralel işleme için örneklerin genomik aralıklarla Kartezyen çarpımlarını oluşturmak amacıyla `combine` kullandık.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Gruplama anahtarlarına göre toplama:** Her demetteki ilk öğeye göre gruplamak, böylece `id` ve `interval` alanlarını paylaşan örnekleri toplayarak teknik replikatları birleştirmek için `groupTuple` kullandık.

    ```groovy
    channel.groupTuple()
    ```

6.  **Veri yapısını optimize etme:** Belirli alanları çıkarmak için `subMap` kullandık ve dönüşümleri yeniden kullanılabilir kılmak için adlandırılmış closure oluşturduk.

    - Bir map'ten belirli alanları çıkarma

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Yeniden kullanılabilir dönüşümler için adlandırılmış closure kullanma

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Ek kaynaklar

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Sırada ne var?

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
