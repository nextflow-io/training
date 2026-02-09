# Bölme ve Gruplama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, verilerle esnek bir şekilde çalışmak için güçlü araçlar sağlar. Temel bir yetenek, verileri farklı akışlara bölmek ve ardından ilgili öğeleri tekrar bir araya getirmektir. Bu, özellikle farklı türdeki örnekleri ayrı ayrı işlemeden önce analiz için sonuçları birleştirmeniz gereken biyoinformatik iş akışlarında değerlidir.

Bunu posta sıralama gibi düşünün: mektupları varış yerine göre ayırırsınız, her yığını farklı şekilde işlersiniz, ardından aynı kişiye giden öğeleri yeniden birleştirirsiniz. Nextflow bunu bilimsel verilerle gerçekleştirmek için özel operatörler kullanır. Bu yaklaşım aynı zamanda dağıtık hesaplama ve biyoinformatik iş akışlarında yaygın olarak **scatter/gather** deseni olarak bilinir.

Nextflow'un kanal sistemi bu esnekliğin merkezindedir. Kanallar, iş akışınızın farklı bölümlerini birbirine bağlayarak verilerin analiziniz boyunca akmasını sağlar. Tek bir veri kaynağından birden fazla kanal oluşturabilir, her kanalı farklı şekilde işleyebilir ve ardından gerektiğinde kanalları tekrar birleştirebilirsiniz. Bu yaklaşım, karmaşık biyoinformatik analizlerinin dallanma ve birleşme yollarını doğal olarak yansıtan iş akışları tasarlamanıza olanak tanır.

### Öğrenme hedefleri

Bu yan görevde, Nextflow'un kanal operatörlerini kullanarak verileri bölmeyi ve gruplamayı öğreneceksiniz.
Örnek bilgilerini ve ilişkili veri dosyalarını içeren bir CSV dosyasıyla başlayacağız, ardından bu verileri manipüle edip yeniden düzenleyeceğiz.

Bu yan görevin sonunda, aşağıdaki teknikleri kullanarak veri akışlarını etkili bir şekilde ayırabilecek ve birleştirebileceksiniz:

- `splitCsv` kullanarak dosyalardan veri okuma
- `filter` ve `map` ile verileri filtreleme ve dönüştürme
- `join` ve `groupTuple` kullanarak ilgili verileri birleştirme
- Paralel işleme için `combine` ile veri kombinasyonları oluşturma
- `subMap` ve tekilleştirme stratejileri kullanarak veri yapısını optimize etme
- Kanal yapılarını manipüle etmenize yardımcı olacak adlandırılmış closure'lar ile yeniden kullanılabilir fonksiyonlar oluşturma

Bu beceriler, temiz ve sürdürülebilir kod yapısını korurken birden fazla girdi dosyasını ve farklı veri türlerini verimli bir şekilde işleyebilen iş akışları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmalısınız:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olun.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, dosyalarla çalışma, meta veriler) rahatça kullanabiliyor olun

**İsteğe bağlı:** Önce [İş akışlarında meta veriler](./metadata.md) yan görevini tamamlamanızı öneririz.
Bu, burada yoğun olarak kullanacağımız `splitCsv` ile CSV dosyalarını okuma ve meta haritaları oluşturma temellerini kapsar.

---

## 0. Başlarken

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/splitting_and_grouping
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri gözden geçirin

Bir ana iş akışı dosyası ve `samplesheet.csv` adlı bir örnek listesi içeren bir `data` dizini bulacaksınız.

```console title="Dizin içeriği"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Örnek listesi, hasta kimliği, örnek tekrar numarası, tür (normal veya tümör) ve varsayımsal veri dosyalarının yolları (gerçekte mevcut değiller, ancak varmış gibi davranacağız) dahil olmak üzere farklı hastalardan alınan örnekler hakkında bilgi içerir.

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

Bu örnek listesi, üç hastadan (A, B, C) sekiz örneği listeler.

Her hasta için, `tumor` (tipik olarak tümör biyopsilerinden kaynaklanan) veya `normal` (sağlıklı doku veya kandan alınan) türünde örneklerimiz var.
Kanser analizi konusunda bilginiz yoksa, bunun karşılaştırmalı analizler yapmak için eşleştirilmiş tümör/normal örnekleri kullanan bir deneysel modele karşılık geldiğini bilin.

Özellikle hasta A için iki set teknik tekrarımız (repeats) var.

!!! note

    Bu deneysel tasarıma aşina değilseniz endişelenmeyin, bu eğitimi anlamak için kritik değil.

#### Görevi gözden geçirin

Meydan okumanız şunları yapacak bir Nextflow iş akışı yazmaktır:

1. Bir CSV dosyasından örnek verilerini **okumak** ve meta haritalarla yapılandırmak
2. Örnekleri türe göre (normal vs tümör) farklı kanallara **ayırmak**
3. Eşleşen tümör/normal çiftlerini hasta kimliği ve tekrar numarasına göre **birleştirmek**
4. Paralel işleme için örnekleri genomik aralıklara **dağıtmak**
5. İlgili örnekleri alt akış analizi için tekrar **gruplamak**

Bu, bağımsız işleme için verileri bölmeniz, ardından karşılaştırmalı analiz için ilgili öğeleri yeniden birleştirmeniz gereken yaygın bir biyoinformatik desenini temsil eder.

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Örnek verileri okuma

### 1.1. `splitCsv` ile örnek verileri okuma ve meta haritaları oluşturma

Örnek verileri `splitCsv` ile okuyarak ve meta harita desenine göre düzenleyerek başlayalım. `main.nf` dosyasında, iş akışını zaten başlattığımızı göreceksiniz.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note

    Bu eğitim boyunca, Nextflow kanalları olduklarını açıkça belirtmek için tüm kanal değişkenleri için `ch_` önekini kullanacağız.

[İş akışlarında meta veriler](./metadata.md) yan görevini tamamladıysanız, bu deseni tanıyacaksınız. CSV'yi okumak için `splitCsv` kullanacağız ve meta verileri dosya yollarından ayırmak için verileri hemen bir meta haritasıyla yapılandıracağız.

!!! info

    Bu eğitimde `map` adında iki farklı kavramla karşılaşacağız:

    - **Veri yapısı**: Anahtar-değer çiftlerini saklayan Groovy map (diğer dillerdeki sözlükler/hash'lere eşdeğer)
    - **Kanal operatörü**: Bir kanaldaki öğeleri dönüştüren `.map()` operatörü

    Bağlama göre hangisini kastettiğimizi netleştireceğiz, ancak bu ayrım Nextflow ile çalışırken anlaşılması önemlidir.

Bu değişiklikleri `main.nf` dosyasına uygulayın:

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

Bu, `splitCsv` işlemini (CSV'yi başlıklarla okuma) ve `map` işlemini (verileri `[meta, file]` demetleri olarak yapılandırma) tek adımda birleştirir. Bu değişikliği uygulayın ve pipeline'ı çalıştırın:

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

Artık her öğenin bir `[meta, file]` demeti olduğu bir kanalımız var - meta veriler dosya yollarından ayrılmış. Bu yapı, meta veri alanlarına göre iş yükümüzü bölmemize ve gruplamamıza olanak tanır.

---

## 2. Verileri filtreleme ve dönüştürme

### 2.1. `filter` ile verileri filtreleme

Verileri bir koşula göre filtrelemek için [`filter` operatörünü](https://www.nextflow.io/docs/latest/operator.html#filter) kullanabiliriz. Diyelim ki sadece normal örnekleri işlemek istiyoruz. Bunu `type` alanına göre verileri filtreleyerek yapabiliriz. Bunu `view` operatöründen önce ekleyelim.

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

Verileri başarıyla filtreledik ve sadece normal örnekleri dahil ettik. Bunun nasıl çalıştığını özetleyelim.

`filter` operatörü, kanaldaki her öğeye uygulanan bir closure alır. Closure `true` döndürürse, öğe dahil edilir; `false` döndürürse, öğe hariç tutulur.

Bizim durumumuzda, sadece `meta.type == 'normal'` olan örnekleri tutmak istiyoruz. Closure, her örneğe atıfta bulunmak için `meta,file` demetini kullanır, `meta.type` ile örnek türüne erişir ve `'normal'` ile eşit olup olmadığını kontrol eder.

Bu, yukarıda tanıttığımız tek closure ile gerçekleştirilir:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Ayrı filtrelenmiş kanallar oluşturma

Şu anda filtreyi doğrudan CSV'den oluşturulan kanala uyguluyoruz, ancak bunu birden fazla şekilde filtrelemek istiyoruz, bu yüzden normal örnekler için ayrı bir filtrelenmiş kanal oluşturmak üzere mantığı yeniden yazalım:

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

Normal ve tümör örneklerini iki farklı kanala ayırdık ve çıktıda farklı etiketlemek için `view()` fonksiyonuna bir closure sağladık: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Özet

Bu bölümde şunları öğrendiniz:

- **Verileri filtreleme**: `filter` ile verileri nasıl filtreleyeceğinizi
- **Verileri bölme**: Bir koşula göre verileri farklı kanallara nasıl böleceğinizi
- **Verileri görüntüleme**: Verileri yazdırmak ve farklı kanallardan çıktıyı etiketlemek için `view` nasıl kullanılır

Artık normal ve tümör örneklerini iki farklı kanala ayırdık. Sırada, normal ve tümör örneklerini `id` alanına göre birleştireceğiz.

---

## 3. Tanımlayıcılara göre kanalları birleştirme

Önceki bölümde, normal ve tümör örneklerini iki farklı kanala ayırdık. Bunlar, türlerine göre belirli süreçler veya iş akışları kullanılarak bağımsız olarak işlenebilir. Peki aynı hastadan normal ve tümör örneklerini karşılaştırmak istediğimizde ne olur? Bu noktada, `id` alanlarına göre eşleştirdiğimizden emin olarak onları tekrar bir araya getirmemiz gerekir.

Nextflow, kanalları birleştirmek için birçok yöntem içerir, ancak bu durumda en uygun operatör [`join`](https://www.nextflow.io/docs/latest/operator.html#join)'dir. SQL'e aşinaysanız, `JOIN` işlemi gibi davranır; burada birleştirme anahtarını ve gerçekleştirilecek birleştirme türünü belirtiriz.

### 3.1. Hasta kimliğine göre birleştirmek için `map` ve `join` kullanma

[`join`](https://www.nextflow.io/docs/latest/operator.html#join) belgelerini kontrol edersek, varsayılan olarak iki kanalı her demetteki ilk öğeye göre birleştirdiğini görebiliriz.

#### 3.1.1. Veri yapısını kontrol etme

Konsol çıktısı hala mevcut değilse, veri yapımızı kontrol etmek ve `id` alanına göre birleştirmek için nasıl değiştirmemiz gerektiğini görmek için pipeline'ı çalıştıralım.

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

`id` alanının her meta haritasındaki ilk öğe olduğunu görebiliriz. `join`'in çalışması için, her demetteki `id` alanını izole etmeliyiz. Bundan sonra, iki kanalı birleştirmek için basitçe `join` operatörünü kullanabiliriz.

#### 3.1.2. `id` alanını izole etme

`id` alanını izole etmek için, [`map` operatörünü](https://www.nextflow.io/docs/latest/operator.html#map) kullanarak `id` alanını ilk öğe olarak içeren yeni bir demet oluşturabiliriz.

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

İnce olabilir, ancak her demetteki ilk öğenin `id` alanı olduğunu görebilmelisiniz.

#### 3.1.3. İki kanalı birleştirme

Şimdi `id` alanına göre iki kanalı birleştirmek için `join` operatörünü kullanabiliriz.

Bir kez daha, birleştirilmiş çıktıları yazdırmak için `view` kullanacağız.

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

Çok geniş olduğu için söylemek biraz zor, ancak örneklerin `id` alanına göre birleştirildiğini görebilmelisiniz. Her demet artık şu formata sahip:

- `id`: Örnek kimliği
- `normal_meta_map`: Tür, tekrar ve bam dosyasının yolu dahil normal örnek meta verileri
- `normal_sample_file`: Normal örnek dosyası
- `tumor_meta_map`: Tür, tekrar ve bam dosyasının yolu dahil tümör örnek meta verileri
- `tumor_sample`: Tür, tekrar ve bam dosyasının yolu dahil tümör örneği

!!! warning

    `join` operatörü eşleşmeyen demetleri atacaktır. Bu örnekte, tüm örneklerin tümör ve normal için eşleştiğinden emin olduk, ancak bu doğru değilse eşleşmeyen demetleri tutmak için `remainder: true` parametresini kullanmalısınız. Daha fazla ayrıntı için [belgelere](https://www.nextflow.io/docs/latest/operator.html#join) bakın.

Artık bir demetteki bir alanı izole etmek için `map` nasıl kullanılacağını ve ilk alana göre demetleri birleştirmek için `join` nasıl kullanılacağını biliyorsunuz.
Bu bilgiyle, paylaşılan bir alana göre kanalları başarıyla birleştirebiliriz.

Sırada, birden fazla alana göre birleştirmek istediğiniz durumu ele alacağız.

### 3.2. Birden fazla alana göre birleştirme

SampleA için 2 tekrarımız var, ancak sampleB ve sampleC için sadece 1 tane. Bu durumda `id` alanını kullanarak onları etkili bir şekilde birleştirebildik, ancak senkronize olmasalardı ne olurdu? Farklı tekrarlardan normal ve tümör örneklerini karıştırabilirdik!

Bunu önlemek için birden fazla alana göre birleştirebiliriz. Bunu başarmanın aslında birden fazla yolu var, ancak hem örnek `id` hem de `replicate` numarasını içeren yeni bir birleştirme anahtarı oluşturmaya odaklanacağız.

Yeni bir birleştirme anahtarı oluşturarak başlayalım. Bunu daha önce olduğu gibi [`map` operatörünü](https://www.nextflow.io/docs/latest/operator.html#map) kullanarak `id` ve `repeat` alanlarını ilk öğe olarak içeren yeni bir demet oluşturarak yapabiliriz.

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

Şimdi birleştirmenin hem `id` hem de `repeat` alanlarını kullanarak gerçekleştiğini görmeliyiz. İş akışını çalıştırın:

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

Her birleştirilmiş sonucun ilk öğesi olarak iki öğeli bir demetin (`id` ve `repeat` alanları) nasıl olduğuna dikkat edin. Bu, karmaşık öğelerin birleştirme anahtarı olarak nasıl kullanılabileceğini gösterir ve aynı koşullardan örnekler arasında oldukça karmaşık eşleştirmeyi mümkün kılar.

Farklı anahtarlara göre birleştirmenin daha fazla yolunu keşfetmek istiyorsanız, ek seçenekler ve örnekler için [join operatörü belgelerine](https://www.nextflow.io/docs/latest/operator.html#join) bakın.

### 3.3. Yeni bir birleştirme anahtarı oluşturmak için `subMap` kullanma

Önceki yaklaşım, birleştirme anahtarımızdan alan adlarını kaybeder - `id` ve `repeat` alanları sadece bir değer listesi haline gelir. Daha sonra erişim için alan adlarını korumak için [`subMap` metodunu](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) kullanabiliriz.

`subMap` metodu, bir haritadan yalnızca belirtilen anahtar-değer çiftlerini çıkarır. Burada birleştirme anahtarımızı oluşturmak için sadece `id` ve `repeat` alanlarını çıkaracağız.

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

Artık sadece `id` ve `repeat` alanlarını içeren değil, aynı zamanda alan adlarını da koruyan yeni bir birleştirme anahtarımız var, böylece daha sonra onlara isimle erişebiliriz, örneğin `meta.id` ve `meta.repeat`.

### 3.4. Map'te adlandırılmış closure kullanma

Tekrarı önlemek ve hataları azaltmak için adlandırılmış bir closure kullanabiliriz. Adlandırılmış bir closure, birden fazla yerde çağırabileceğimiz yeniden kullanılabilir bir fonksiyon oluşturmamıza olanak tanır.

Bunu yapmak için önce closure'ı yeni bir değişken olarak tanımlarız:

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

Yeniden kullanabileceğimiz adlandırılmış bir değişken olarak map dönüşümünü tanımladık.

Ayrıca dosya yolunu `file()` kullanarak bir Path nesnesine dönüştürdüğümüze dikkat edin, böylece bu kanalı alan herhangi bir süreç dosyayı doğru şekilde işleyebilir (daha fazla bilgi için [Dosyalarla çalışma](./working_with_files.md) bölümüne bakın).

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

!!! note

    `map` operatörü, closure'ı argüman olarak geçirmek için `{ }` yerine `( )` kullanmaya geçti. Bunun nedeni, `map` operatörünün argüman olarak bir closure beklediği ve `{ }` anonim bir closure tanımlamak için kullanıldığıdır. Adlandırılmış bir closure çağırırken `( )` sözdizimini kullanın.

Her şeyin hala çalıştığını kontrol etmek için iş akışını bir kez daha çalıştırın:

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

Adlandırılmış bir closure kullanmak, aynı dönüşümü birden fazla yerde yeniden kullanmamıza, hata riskini azaltmamıza ve kodu daha okunabilir ve sürdürülebilir hale getirmemize olanak tanır.

### 3.5. Veri tekrarını azaltma

İş akışımızda çok fazla tekrarlanan veri var. Birleştirilmiş örneklerdeki her öğe `id` ve `repeat` alanlarını tekrarlıyor. Bu bilgi zaten gruplama anahtarında mevcut olduğundan, bu fazlalıktan kaçınabiliriz. Hatırlatma olarak, mevcut veri yapımız şöyle görünüyor:

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

`id` ve `repeat` alanları gruplama anahtarında mevcut olduğundan, tekrardan kaçınmak için bunları her kanal öğesinin geri kalanından kaldıralım. Bunu sadece `type` alanını içeren yeni bir harita oluşturmak için `subMap` metodunu kullanarak yapabiliriz. Bu yaklaşım, veri yapımızdaki fazlalığı ortadan kaldırırken tüm gerekli bilgileri korumamıza olanak tanır.

=== "Sonra"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Şimdi closure, ilk öğenin `id` ve `repeat` alanlarını içerdiği ve ikinci öğenin sadece `type` alanını içerdiği bir demet döndürür. Bu, `id` ve `repeat` bilgilerini gruplama anahtarında bir kez saklayarak fazlalığı ortadan kaldırırken tüm gerekli bilgileri korur.

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

`id` ve `repeat` alanlarını gruplama anahtarında sadece bir kez belirttiğimizi ve örnek verilerinde `type` alanına sahip olduğumuzu görebiliriz. Herhangi bir bilgi kaybetmedik ama kanal içeriklerimizi daha özlü hale getirmeyi başardık.

### 3.6. Gereksiz bilgileri kaldırma

Yukarıda tekrarlanan bilgileri kaldırdık, ancak kanallarımızda hala başka gereksiz bilgiler var.

Başlangıçta, `filter` kullanarak normal ve tümör örneklerini ayırdık, ardından bunları `id` ve `repeat` anahtarlarına göre birleştirdik. `join` operatörü, demetlerin birleştirildiği sırayı korur, bu nedenle bizim durumumuzda, sol tarafta normal örnekler ve sağ tarafta tümör örnekleri ile, ortaya çıkan kanal bu yapıyı korur: `id, <normal öğeler>, <tümör öğeler>`.

Kanalımızdaki her öğenin konumunu bildiğimiz için, `[type:normal]` ve `[type:tumor]` meta verilerini bırakarak yapıyı daha da basitleştirebiliriz.

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

### Özet

Bu bölümde şunları öğrendiniz:

- **Demetleri Manipüle Etme**: Bir demetteki bir alanı izole etmek için `map` nasıl kullanılır
- **Demetleri Birleştirme**: İlk alana göre demetleri birleştirmek için `join` nasıl kullanılır
- **Birleştirme Anahtarları Oluşturma**: Yeni bir birleştirme anahtarı oluşturmak için `subMap` nasıl kullanılır
- **Adlandırılmış Closure'lar**: Map'te adlandırılmış bir closure nasıl kullanılır
- **Birden Fazla Alan Birleştirme**: Daha hassas eşleştirme için birden fazla alana göre nasıl birleştirilir
- **Veri Yapısı Optimizasyonu**: Gereksiz bilgileri kaldırarak kanal yapısı nasıl düzenlenir

Artık bir örnek listesini bölebilen, normal ve tümör örneklerini filtreleyebilen, bunları örnek kimliği ve tekrar numarasına göre birleştirebilen ve ardından sonuçları yazdırabilen bir iş akışınız var.

Bu, bağımsız işlemeden sonra örnekleri veya diğer veri türlerini eşleştirmeniz gereken biyoinformatik iş akışlarında yaygın bir desendir, bu nedenle yararlı bir beceridir. Sırada, bir örneği birden çok kez tekrarlamaya bakacağız.

## 4. Örnekleri aralıklara yayma

Biyoinformatik iş akışlarında temel bir desen, analizi genomik bölgelere dağıtmaktır. Örneğin, varyant çağırma, genomu aralıklara (kromozomlar veya daha küçük bölgeler gibi) bölerek paralelleştirilebilir. Bu paralelleştirme stratejisi, hesaplama yükünü birden fazla çekirdek veya düğüme dağıtarak pipeline verimliliğini önemli ölçüde artırır ve genel yürütme süresini azaltır.

Aşağıdaki bölümde, örnek verilerimizi birden fazla genomik aralığa nasıl dağıtacağımızı göstereceğiz. Her örneği her aralıkla eşleştireceğiz, böylece farklı genomik bölgelerin paralel işlenmesine izin vereceğiz. Bu, veri kümesi boyutumuzu aralık sayısıyla çarpacak ve daha sonra bir araya getirilebilecek birden fazla bağımsız analiz birimi oluşturacaktır.

### 4.1. `combine` kullanarak örnekleri aralıklara yayma

Bir aralık kanalı oluşturarak başlayalım. Hayatı basit tutmak için, manuel olarak tanımlayacağımız sadece 3 aralık kullanacağız. Gerçek bir iş akışında, bunları bir dosya girdisinden okuyabilir veya hatta birçok aralık dosyası içeren bir kanal oluşturabilirsiniz.

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

Şimdi unutmayın, her örneği her aralık için tekrarlamak istiyoruz. Bu bazen örneklerin ve aralıkların Kartezyen çarpımı olarak adlandırılır. Bunu [`combine` operatörünü](https://www.nextflow.io/docs/latest/operator.html#combine) kullanarak başarabiliriz. Bu, kanal 1'deki her öğeyi alacak ve kanal 2'deki her öğe için tekrarlayacaktır. İş akışımıza bir combine operatörü ekleyelim:

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

Başarılı! Her örneği 3 aralık listemizde her tek aralık için tekrarladık. Kanalımızdaki öğe sayısını etkili bir şekilde üçe katladık.

Yine de okumak biraz zor, bu yüzden bir sonraki bölümde düzenleyeceğiz.

### 4.2. Kanalı düzenleme

Örnek verilerimizi düzenlemek ve yeniden yapılandırmak için `map` operatörünü kullanabiliriz, böylece anlaşılması daha kolay olur. Aralık dizesini ilk öğedeki birleştirme haritasına taşıyalım.

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

Bu map işleminin adım adım ne yaptığını inceleyelim.

İlk olarak, kodu daha okunabilir hale getirmek için adlandırılmış parametreler kullanıyoruz. `grouping_key`, `normal`, `tumor` ve `interval` adlarını kullanarak, demetteki öğelere indeks yerine isimle başvurabiliriz:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Sonra, `grouping_key`'i `interval` alanıyla birleştiriyoruz. `grouping_key`, `id` ve `repeat` alanlarını içeren bir haritadır. `interval` ile yeni bir harita oluşturuyoruz ve bunları Groovy'nin harita toplama (`+`) kullanarak birleştiriyoruz:

```groovy
                grouping_key + [interval: interval],
```

Son olarak, bunu üç öğeli bir demet olarak döndürüyoruz: birleştirilmiş meta veri haritası, normal örnek dosyası ve tümör örnek dosyası:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Tekrar çalıştıralım ve kanal içeriğini kontrol edelim:

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

Verilerinizi doğru yapıya zorlamak için `map` kullanmak zor olabilir, ancak etkili veri manipülasyonu için çok önemlidir.

Artık her örneği tüm genomik aralıklara yayarak paralel olarak işlenebilecek birden fazla bağımsız analiz birimi oluşturduk. Peki ilgili örnekleri tekrar bir araya getirmek istersek ne olur? Bir sonraki bölümde, ortak özellikleri paylaşan örnekleri nasıl gruplayacağımızı öğreneceğiz.

### Özet

Bu bölümde şunları öğrendiniz:

- **Örnekleri aralıklara yayma**: Örnekleri aralıklara yaymak için `combine` nasıl kullanılır
- **Kartezyen çarpımlar oluşturma**: Örneklerin ve aralıkların tüm kombinasyonları nasıl oluşturulur
- **Kanal yapısını düzenleme**: Daha iyi okunabilirlik için verileri yeniden yapılandırmak için `map` nasıl kullanılır
- **Paralel işleme hazırlığı**: Dağıtık analiz için veriler nasıl ayarlanır

## 5. `groupTuple` kullanarak örnekleri toplama

Önceki bölümlerde, bir girdi dosyasından verileri nasıl böleceğimizi ve belirli alanlara göre (bizim durumumuzda normal ve tümör örnekleri) nasıl filtreleyeceğimizi öğrendik. Ancak bu sadece tek bir birleştirme türünü kapsar. Peki örnekleri belirli bir özelliğe göre gruplamak istersek ne olur? Örneğin, eşleşen normal-tümör çiftlerini birleştirmek yerine, türlerinden bağımsız olarak "sampleA"dan tüm örnekleri birlikte işlemek isteyebiliriz. Bu desen, sonunda sonuçları karşılaştırmadan veya birleştirmeden önce ilgili örnekleri verimlilik nedenleriyle ayrı ayrı işlemek isteyebileceğiniz biyoinformatik iş akışlarında yaygındır.

Nextflow bunu yapmak için yerleşik yöntemler içerir, bakacağımız ana yöntem `groupTuple`'dır.

Aynı `id` ve `interval` alanlarına sahip tüm örneklerimizi gruplayarak başlayalım, bu, teknik tekrarları gruplamak istediğimiz ancak anlamlı olarak farklı örnekleri ayrı tutmak istediğimiz bir analizin tipik bir örneği olacaktır.

Bunu yapmak için, gruplama değişkenimizi izole etmeliyiz, böylece onları tek başına kullanabiliriz.

İlk adım önceki bölümde yaptığımıza benzer. Gruplama değişkenimizi demetin ilk öğesi olarak izole etmeliyiz. Unutmayın, ilk öğemiz şu anda `id`, `repeat` ve `interval` alanlarının bir haritasıdır:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Haritadan `id` ve `interval` alanlarımızı izole etmek için daha önce kullandığımız `subMap` metodunu yeniden kullanabiliriz. Daha önce olduğu gibi, her örnek için demetin ilk öğesine `subMap` metodunu uygulamak için `map` operatörünü kullanacağız.

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

Tekrar çalıştıralım ve kanal içeriğini kontrol edelim:

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

`id` ve `interval` alanlarını başarıyla izole ettiğimizi, ancak örnekleri henüz gruplamadığımızı görebiliriz.

!!! note

    Burada `replicate` alanını atıyoruz. Bunun nedeni, daha fazla alt akış işleme için buna ihtiyacımız olmamasıdır. Bu eğitimi tamamladıktan sonra, daha sonraki gruplamayı etkilemeden onu dahil edip edemeyeceğinizi görün!

Şimdi örnekleri bu yeni gruplama öğesine göre [`groupTuple` operatörünü](https://www.nextflow.io/docs/latest/operator.html#grouptuple) kullanarak gruplayalım.

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

İşte bu kadar! Sadece tek bir satır kod ekledik. Çalıştırdığımızda ne olduğunu görelim:

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

Verilerimizin yapısının değiştiğine ve her kanal öğesi içinde dosyaların artık `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` gibi demetlerde bulunduğuna dikkat edin. Bunun nedeni, `groupTuple` kullandığımızda, Nextflow'un bir grubun her örneği için tek dosyaları birleştirmesidir. Verileri alt akışta işlemeye çalışırken bunu hatırlamak önemlidir.

!!! note

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose), groupTuple'ın tersidir. Bir kanaldaki öğeleri açar ve düzleştirir. `transpose` eklemeyi deneyin ve yukarıda gerçekleştirdiğimiz gruplamayı geri alın!

### Özet

Bu bölümde şunları öğrendiniz:

- **İlgili örnekleri gruplama**: Ortak özelliklere göre örnekleri toplamak için `groupTuple` nasıl kullanılır
- **Gruplama anahtarlarını izole etme**: Gruplama için belirli alanları çıkarmak için `subMap` nasıl kullanılır
- **Gruplanmış veri yapılarını işleme**: `groupTuple` tarafından oluşturulan iç içe yapıyla nasıl çalışılır
- **Teknik tekrar işleme**: Aynı deneysel koşulları paylaşan örnekler nasıl gruplandırılır

---

## Özet

Bu yan görevde, kanalları kullanarak verileri nasıl böleceğinizi ve gruplayacağınızı öğrendiniz.

Verileri pipeline boyunca akarken değiştirerek, döngüler veya while ifadeleri kullanmadan ölçeklenebilir bir pipeline oluşturabilirsiniz, bu da daha geleneksel yaklaşımlara göre çeşitli avantajlar sunar:

- Ek kod olmadan istediğimiz kadar çok veya az girdiye ölçeklenebiliriz
- İterasyon yerine verilerin pipeline boyunca akışını işlemeye odaklanırız
- Gerektiği kadar karmaşık veya basit olabiliriz
- Pipeline daha bildirimsel hale gelir, nasıl olması gerektiğinden ziyade ne olması gerektiğine odaklanır
- Nextflow, bağımsız işlemleri paralel olarak çalıştırarak yürütmeyi bizim için optimize edecektir

Bu kanal işlemlerinde ustalaşmak, döngülere veya yinelemeli programlamaya başvurmadan karmaşık veri ilişkilerini işleyen esnek, ölçeklenebilir pipeline'lar oluşturmanıza olanak tanır ve Nextflow'un yürütmeyi optimize etmesine ve bağımsız işlemleri otomatik olarak paralelleştirmesine izin verir.

### Temel desenler

1.  **Yapılandırılmış girdi verileri oluşturma:** Meta haritalarla bir CSV dosyasından başlama ([İş akışlarında meta veriler](./metadata.md) desenlerini temel alarak)

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

3.  **Eşleşen örnekleri birleştirme:** `id` ve `repeat` alanlarına göre ilgili örnekleri yeniden birleştirmek için `join` kullandık

    - Anahtara göre iki kanalı birleştirme (demetin ilk öğesi)

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

4.  **Aralıklara dağıtma:** Paralel işleme için örneklerin genomik aralıklarla Kartezyen çarpımlarını oluşturmak için `combine` kullandık.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Gruplama anahtarlarına göre toplama:** Her demetteki ilk öğeye göre gruplamak için `groupTuple` kullandık, böylece `id` ve `interval` alanlarını paylaşan örnekleri toplayarak teknik tekrarları birleştirdik.

    ```groovy
    channel.groupTuple()
    ```

6.  **Veri yapısını optimize etme:** Belirli alanları çıkarmak için `subMap` kullandık ve dönüşümleri yeniden kullanılabilir hale getirmek için adlandırılmış bir closure oluşturduk.

    - Bir haritadan belirli alanları çıkarma

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

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
