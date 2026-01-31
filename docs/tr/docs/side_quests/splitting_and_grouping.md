# Bölme ve Gruplama

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow, verilerle esnek bir şekilde çalışmak için güçlü araçlar sağlar. Temel bir özellik, verileri farklı akışlara bölmek ve ardından ilişkili öğeleri tekrar birleştirmektir. Bu, özellikle analiz için sonuçları birleştirmeden önce farklı örnek türlerini ayrı ayrı işlemeniz gereken biyoinformatik iş akışlarında değerlidir.

Bunu mektup sıralama gibi düşünün: mektupları varış yerine göre ayırırsınız, her yığını farklı şekilde işlersiniz, ardından aynı kişiye giden öğeleri yeniden birleştirirsiniz. Nextflow, bilimsel verilerle bunu başarmak için özel operatörler kullanır. Bu yaklaşım aynı zamanda dağıtık hesaplama ve biyoinformatik iş akışlarında yaygın olarak **scatter/gather** deseni olarak bilinir.

Nextflow'un channel sistemi bu esnekliğin merkezindedir. Channel'lar iş akışınızın farklı bölümlerini bağlayarak verilerin analiziniz boyunca akmasını sağlar. Tek bir veri kaynağından birden fazla channel oluşturabilir, her channel'ı farklı şekilde işleyebilir ve gerektiğinde channel'ları tekrar birleştirebilirsiniz. Bu yaklaşım, karmaşık biyoinformatik analizlerin dallanma ve birleşme yollarını doğal olarak yansıtan iş akışları tasarlamanıza olanak tanır.

### Öğrenme hedefleri

Bu yan görevde, Nextflow'un channel operatörlerini kullanarak verileri bölmeyi ve gruplamayı öğreneceksiniz.
Örnek bilgileri ve ilişkili veri dosyalarını içeren bir CSV dosyasıyla başlayacak, ardından bu verileri manipüle edip yeniden düzenleyeceğiz.

Bu yan görevin sonunda, aşağıdaki teknikleri kullanarak veri akışlarını etkili bir şekilde ayırıp birleştirebileceksiniz:

- `splitCsv` kullanarak dosyalardan veri okuma
- `filter` ve `map` ile verileri filtreleme ve dönüştürme
- `join` ve `groupTuple` kullanarak ilişkili verileri birleştirme
- Paralel işleme için `combine` ile veri kombinasyonları oluşturma
- `subMap` ve tekilleştirme stratejileri kullanarak veri yapısını optimize etme
- Channel yapılarını manipüle etmenize yardımcı olacak adlandırılmış closure'lar ile yeniden kullanılabilir fonksiyonlar oluşturma

Bu beceriler, temiz ve bakımı kolay kod yapısını korurken birden fazla girdi dosyası ve farklı veri türlerini verimli bir şekilde işleyebilen iş akışları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (process'ler, channel'lar, operatörler, dosyalarla çalışma, meta veri) rahatça kullanabiliyor olmalısınız

**İsteğe bağlı:** Önce [İş akışlarında meta veri](./metadata.md) yan görevini tamamlamanızı öneririz.
Bu, burada yoğun şekilde kullanacağımız `splitCsv` ile CSV dosyalarını okuma ve meta map oluşturma temellerini kapsar.

---

## 0. Başlangıç

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md)'nda açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/splitting_and_grouping
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

Bir ana workflow dosyası ve `samplesheet.csv` adlı bir örnek sayfası içeren bir `data` dizini bulacaksınız.

```console title="Dizin içeriği"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

Örnek sayfası, hasta kimliği, örnek tekrar numarası, tür (normal veya tümör) ve varsayımsal veri dosyalarına giden yollar (gerçekte mevcut değil, ancak varmış gibi davranacağız) dahil olmak üzere farklı hastalardan alınan örnekler hakkında bilgi içerir.

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

Bu örnek sayfası, üç hastadan (A, B, C) sekiz örnek listeler.

Her hasta için `tumor` (genellikle tümör biyopsilerinden kaynaklanan) veya `normal` (sağlıklı doku veya kandan alınan) türünde örneklerimiz var.
Kanser analiziyle tanışık değilseniz, bunun karşılaştırmalı analizler gerçekleştirmek için eşleştirilmiş tümör/normal örnekleri kullanan deneysel bir modele karşılık geldiğini bilin.

Özellikle A hastası için iki set teknik replik (tekrar) var.

!!! note "Not"

    Bu deneysel tasarıma aşina değilseniz endişelenmeyin, bu eğitimi anlamak için kritik değil.

#### Görevi inceleyin

Göreviniz şunları yapacak bir Nextflow iş akışı yazmaktır:

1. Bir CSV dosyasından örnek verilerini **okuma** ve meta map'lerle yapılandırma
2. Türe göre (normal vs tümör) örnekleri farklı channel'lara **ayırma**
3. Hasta kimliği ve replik numarasına göre eşleşen tümör/normal çiftlerini **birleştirme**
4. Paralel işleme için örnekleri genomik aralıklara **dağıtma**
5. Aşağı akış analizi için ilişkili örnekleri tekrar **gruplama**

Bu, verileri bağımsız işleme için bölmeniz ve ardından karşılaştırmalı analiz için ilişkili öğeleri yeniden birleştirmeniz gereken yaygın bir biyoinformatik desenidir.

#### Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Örnek verilerini okuma

### 1.1. `splitCsv` ile örnek verilerini okuma ve meta map'ler oluşturma

`splitCsv` ile örnek verilerini okuyarak ve meta map deseninde organize ederek başlayalım. `main.nf` dosyasında, iş akışını zaten başlattığımızı göreceksiniz.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Not"

    Bu eğitim boyunca, Nextflow channel'ları olduklarını açıkça belirtmek için tüm channel değişkenleri için `ch_` önekini kullanacağız.

[İş akışlarında meta veri](./metadata.md) yan görevini tamamladıysanız, bu deseni tanıyacaksınız. Meta veriyi dosya yollarından ayırmak için CSV'yi okumak ve verileri hemen bir meta map ile yapılandırmak için `splitCsv` kullanacağız.

!!! info "Bilgi"

    Bu eğitimde `map` olarak adlandırılan iki farklı kavramla karşılaşacağız:

    - **Veri yapısı**: Anahtar-değer çiftlerini saklayan Groovy map (diğer dillerdeki sözlükler/hash'ler ile eşdeğer)
    - **Channel operatörü**: Bir channel'daki öğeleri dönüştüren `.map()` operatörü

    Bağlamda hangisini kastettiğimizi açıklayacağız, ancak Nextflow ile çalışırken bu ayrımı anlamak önemlidir.

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

Bu, `splitCsv` işlemini (başlıklarla CSV'yi okuma) ve `map` işlemini (verileri `[meta, file]` tuple'ları olarak yapılandırma) tek adımda birleştirir. Bu değişikliği uygulayın ve pipeline'ı çalıştırın:

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

Artık her öğenin bir `[meta, file]` tuple'ı olduğu bir channel'ımız var - meta veri dosya yollarından ayrılmış durumda. Bu yapı, meta veri alanlarına göre iş yükümüzü bölmemize ve gruplamamıza olanak tanır.

---

## 2. Verileri filtreleme ve dönüştürme

### 2.1. `filter` ile verileri filtreleme

Bir koşula göre verileri filtrelemek için [`filter` operatörünü](https://www.nextflow.io/docs/latest/operator.html#filter) kullanabiliriz. Diyelim ki sadece normal örnekleri işlemek istiyoruz. Bunu `type` alanına göre verileri filtreleyerek yapabiliriz. Bunu `view` operatöründen önce ekleyelim.

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

Verileri başarıyla sadece normal örnekleri içerecek şekilde filtreledik. Bunun nasıl çalıştığını özetleyelim.

`filter` operatörü, channel'daki her öğeye uygulanan bir closure alır. Closure `true` döndürürse, öğe dahil edilir; `false` döndürürse, öğe hariç tutulur.

Bizim durumumuzda, sadece `meta.type == 'normal'` olan örnekleri tutmak istiyoruz. Closure, her örneğe başvurmak için `meta,file` tuple'ını kullanır, `meta.type` ile örnek türüne erişir ve `'normal'`'e eşit olup olmadığını kontrol eder.

Bu, yukarıda tanıttığımız tek closure ile gerçekleştirilir:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Ayrı filtrelenmiş channel'lar oluşturma

Şu anda filtreyi doğrudan CSV'den oluşturulan channel'a uyguluyoruz, ancak bunu birden fazla şekilde filtrelemek istiyoruz, bu yüzden mantığı normal örnekler için ayrı bir filtrelenmiş channel oluşturmak üzere yeniden yazalım:

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

Verileri başarıyla filtreledik ve normal örnekler için ayrı bir channel oluşturduk.

Tümör örnekleri için de filtrelenmiş bir channel oluşturalım:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal örnek: ' + it}
        ch_tumor_samples
            .view{'Tümör örnek: ' + it}
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

    Tümör örnek: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tümör örnek: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal örnek: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal örnek: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal örnek: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal örnek: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tümör örnek: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tümör örnek: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Normal ve tümör örneklerini iki farklı channel'a ayırdık ve çıktıda farklı şekilde etiketlemek için `view()`'a bir closure sağladık: `ch_tumor_samples.view{'Tümör örnek: ' + it}`.

### Özet

Bu bölümde öğrendikleriniz:

- **Verileri filtreleme**: `filter` ile verileri nasıl filtreleyeceğiniz
- **Verileri bölme**: Bir koşula göre verileri farklı channel'lara nasıl böleceğiniz
- **Verileri görüntüleme**: `view` kullanarak verileri yazdırma ve farklı channel'lardan çıktıları nasıl etiketleyeceğiniz

Artık normal ve tümör örneklerini iki farklı channel'a ayırdık. Şimdi normal ve tümör örneklerini `id` alanına göre birleştireceğiz.

---

## 3. Channel'ları tanımlayıcılara göre birleştirme

Önceki bölümde, normal ve tümör örneklerini iki farklı channel'a ayırdık. Bunlar türlerine göre belirli process'ler veya iş akışları kullanılarak bağımsız olarak işlenebilir. Ancak aynı hastadan normal ve tümör örneklerini karşılaştırmak istediğimizde ne olur? Bu noktada, örnekleri `id` alanına göre eşleştirerek tekrar birleştirmemiz gerekir.

Nextflow, channel'ları birleştirmek için birçok yöntem içerir, ancak bu durumda en uygun operatör [`join`](https://www.nextflow.io/docs/latest/operator.html#join)'dir. SQL'e aşinaysanız, birleştirme anahtarını ve gerçekleştirilecek birleştirme türünü belirttiğimiz `JOIN` işlemi gibi çalışır.

### 3.1. Hasta kimliğine göre birleştirmek için `map` ve `join` kullanma

[`join`](https://www.nextflow.io/docs/latest/operator.html#join) belgelerine bakarsak, varsayılan olarak her tuple'daki ilk öğeye göre iki channel'ı birleştirdiğini görebiliriz.

#### 3.1.1. Veri yapısını kontrol etme

Konsol çıktısı hala mevcut değilse, veri yapımızı kontrol etmek ve `id` alanına göre birleştirmek için nasıl değiştirmemiz gerektiğini görmek için pipeline'ı çalıştıralım.

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tümör örnek: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tümör örnek: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal örnek: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal örnek: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal örnek: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal örnek: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tümör örnek: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tümör örnek: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

`id` alanının her meta map'teki ilk öğe olduğunu görebiliriz. `join`'in çalışması için her tuple'daki `id` alanını izole etmeliyiz. Bundan sonra, iki channel'ı birleştirmek için `join` operatörünü kolayca kullanabiliriz.

#### 3.1.2. `id` alanını izole etme

`id` alanını izole etmek için, ilk öğe olarak `id` alanıyla yeni bir tuple oluşturmak üzere [`map` operatörünü](https://www.nextflow.io/docs/latest/operator.html#map) kullanabiliriz.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal örnek: ' + it}
        ch_tumor_samples
            .view{'Tümör örnek: ' + it}
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal örnek: ' + it}
        ch_tumor_samples
            .view{'Tümör örnek: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tümör örnek: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tümör örnek: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal örnek: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal örnek: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tümör örnek: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tümör örnek: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal örnek: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal örnek: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Belki de ince bir fark, ancak her tuple'daki ilk öğenin `id` alanı olduğunu görmelisiniz.

#### 3.1.3. İki channel'ı birleştirme

Artık `id` alanına göre iki channel'ı birleştirmek için `join` operatörünü kullanabiliriz.

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
            .view{'Normal örnek: ' + it}
        ch_tumor_samples
            .view{'Tümör örnek: ' + it}
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

Çok geniş olduğu için söylemek biraz zor, ancak örneklerin `id` alanına göre birleştirildiğini görmelisiniz. Her tuple artık şu formata sahip:

- `id`: Örnek kimliği
- `normal_meta_map`: Tür, replik ve bam dosyası yolu dahil normal örnek meta verisi
- `normal_sample_file`: Normal örnek dosyası
- `tumor_meta_map`: Tür, replik ve bam dosyası yolu dahil tümör örnek meta verisi
- `tumor_sample`: Tür, replik ve bam dosyası yolu dahil tümör örneği

!!! warning "Uyarı"

    `join` operatörü eşleşmeyen tuple'ları atacaktır. Bu örnekte, tüm örneklerin tümör ve normal için eşleştirildiğinden emin olduk, ancak bu doğru değilse eşleşmeyen tuple'ları tutmak için `remainder: true` parametresini kullanmanız gerekir. Daha fazla ayrıntı için [belgelere](https://www.nextflow.io/docs/latest/operator.html#join) bakın.

Artık bir tuple'daki bir alanı izole etmek için `map`'i ve tuple'ları ilk alana göre birleştirmek için `join`'i nasıl kullanacağınızı biliyorsunuz.
Bu bilgiyle, paylaşılan bir alana göre channel'ları başarıyla birleştirebilirsiniz.

Şimdi, birden fazla alanda birleştirmek istediğiniz durumu ele alacağız.

### 3.2. Birden fazla alanda birleştirme

ÖrnekA için 2 replikimiz var, ancak örnekB ve örnekC için sadece 1'er tane. Bu durumda `id` alanını kullanarak bunları etkili bir şekilde birleştirebildik, ancak eşzamanlı olmasalardı ne olurdu? Farklı repliklerden normal ve tümör örneklerini karıştırabilirdik!

Bunu önlemek için birden fazla alanda birleştirme yapabiliriz. Bunu başarmanın aslında birden fazla yolu var, ancak hem örnek `id`'sini hem de `replicate` numarasını içeren yeni bir birleştirme anahtarı oluşturmaya odaklanacağız.

Yeni bir birleştirme anahtarı oluşturarak başlayalım. Bunu, ilk öğe olarak `id` ve `repeat` alanlarıyla yeni bir tuple oluşturmak için [`map` operatörünü](https://www.nextflow.io/docs/latest/operator.html#map) kullanarak daha önce yaptığımız gibi yapabiliriz.

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

Her birleştirilmiş sonucun ilk öğesi olarak iki öğeli bir tuple (`id` ve `repeat` alanları) olduğunu not edin. Bu, karmaşık öğelerin birleştirme anahtarı olarak nasıl kullanılabileceğini gösterir ve aynı koşullardan gelen örnekler arasında oldukça karmaşık eşleştirme sağlar.

Farklı anahtarlarda birleştirmenin daha fazla yolunu keşfetmek istiyorsanız, ek seçenekler ve örnekler için [join operatörü belgelerine](https://www.nextflow.io/docs/latest/operator.html#join) bakın.

### 3.3. Yeni bir birleştirme anahtarı oluşturmak için `subMap` kullanma

Önceki yaklaşım, birleştirme anahtarımızdaki alan adlarını kaybeder - `id` ve `repeat` alanları sadece bir değerler listesi haline gelir. Daha sonra erişim için alan adlarını korumak için [`subMap` metodunu](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>) kullanabiliriz.

`subMap` metodu, bir map'ten yalnızca belirtilen anahtar-değer çiftlerini çıkarır. Burada birleştirme anahtarımızı oluşturmak için sadece `id` ve `repeat` alanlarını çıkaracağız.

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

Artık yalnızca `id` ve `repeat` alanlarını içeren değil, aynı zamanda alan adlarını koruyan yeni bir birleştirme anahtarımız var, böylece daha sonra bunlara `meta.id` ve `meta.repeat` gibi adla erişebiliriz.

### 3.4. map'te adlandırılmış closure kullanma

Çoğaltmayı önlemek ve hataları azaltmak için adlandırılmış bir closure kullanabiliriz. Adlandırılmış bir closure, birden fazla yerde çağırabileceğimiz yeniden kullanılabilir bir fonksiyon oluşturmamıza olanak tanır.

Bunu yapmak için, önce closure'ı yeni bir değişken olarak tanımlarız:

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

Ayrıca dosya yolunu `file()` kullanarak bir Path nesnesine dönüştürdüğümüzü not edin, böylece bu channel'ı alan herhangi bir process dosyayı doğru şekilde işleyebilir (daha fazla bilgi için [Dosyalarla çalışma](./working_with_files.md)'ya bakın).

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

    `map` operatörü closure'ı argüman olarak geçirmek için `{ }` yerine `( )` kullanmaya geçti. Bunun nedeni, `map` operatörünün argüman olarak bir closure beklemesi ve `{ }` işaretinin anonim bir closure tanımlamak için kullanılmasıdır. Adlandırılmış bir closure çağırırken `( )` sözdizimini kullanın.

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

Adlandırılmış bir closure kullanmak, aynı dönüşümü birden fazla yerde yeniden kullanmamıza olanak tanır, hata riskini azaltır ve kodu daha okunabilir ve bakımı kolay hale getirir.

### 3.5. Veri çoğaltmasını azaltma

İş akışımızda çok fazla çoğaltılmış veri var. Birleştirilmiş örneklerdeki her öğe `id` ve `repeat` alanlarını tekrarlar. Bu bilgi zaten gruplama anahtarında mevcut olduğundan, bu fazlalıktan kaçınabiliriz. Hatırlatma olarak, mevcut veri yapımız şöyle görünüyor:

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

`id` ve `repeat` alanları gruplama anahtarında mevcut olduğundan, çoğaltmayı önlemek için bunları her channel öğesinin geri kalanından kaldıralım. Bunu, yalnızca `type` alanıyla yeni bir map oluşturmak için `subMap` metodunu kullanarak yapabiliriz. Bu yaklaşım, veri yapımızdaki fazlalığı ortadan kaldırırken tüm gerekli bilgileri korumamıza olanak tanır.

=== "Sonra"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Şimdi closure, ilk öğenin `id` ve `repeat` alanlarını, ikinci öğenin sadece `type` alanını içerdiği bir tuple döndürüyor. Bu, `id` ve `repeat` bilgilerini gruplama anahtarında bir kez saklayarak fazlalığı ortadan kaldırır ve tüm gerekli bilgileri korur.

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

`id` ve `repeat` alanlarını gruplama anahtarında yalnızca bir kez belirttiğimizi ve örnek verisinde `type` alanına sahip olduğumuzu görebiliriz. Hiçbir bilgi kaybetmedik ama channel içeriklerimizi daha öz hale getirmeyi başardık.

### 3.6. Gereksiz bilgileri kaldırma

Yukarıda çoğaltılmış bilgileri kaldırdık, ancak channel'larımızda hala başka gereksiz bilgiler var.

Başlangıçta, `filter` kullanarak normal ve tümör örneklerini ayırdık, ardından `id` ve `repeat` anahtarlarına göre birleştirdik. `join` operatörü tuple'ların birleştirildiği sırayı korur, bu nedenle bizim durumumuzda, normal örnekler sol tarafta ve tümör örnekleri sağ tarafta, ortaya çıkan channel bu yapıyı korur: `id, <normal öğeleri>, <tümör öğeleri>`.

Channel'ımızdaki her öğenin konumunu bildiğimiz için, `[type:normal]` ve `[type:tumor]` meta verilerini bırakarak yapıyı daha da basitleştirebiliriz.

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

Bu bölümde öğrendikleriniz:

- **Tuple'ları manipüle etme**: Bir tuple'daki bir alanı izole etmek için `map` nasıl kullanılır
- **Tuple'ları birleştirme**: Tuple'ları ilk alana göre birleştirmek için `join` nasıl kullanılır
- **Birleştirme anahtarları oluşturma**: Yeni bir birleştirme anahtarı oluşturmak için `subMap` nasıl kullanılır
- **Adlandırılmış Closure'lar**: map'te adlandırılmış bir closure nasıl kullanılır
- **Çoklu alan birleştirme**: Daha hassas eşleştirme için birden fazla alanda nasıl birleştirme yapılır
- **Veri yapısı optimizasyonu**: Gereksiz bilgileri kaldırarak channel yapısı nasıl düzenlenir

Artık bir örnek sayfasını bölebilen, normal ve tümör örneklerini filtreleyebilen, bunları örnek kimliği ve replik numarasına göre birleştirebilen ve sonuçları yazdırabilen bir iş akışınız var.

Bu, bağımsız olarak işledikten sonra örnekleri veya diğer veri türlerini eşleştirmeniz gereken biyoinformatik iş akışlarında yaygın bir desendir, bu yüzden yararlı bir beceridir. Şimdi, bir örneği birden fazla kez tekrarlamaya bakacağız.

## 4. Örnekleri aralıklara yayma

Biyoinformatik iş akışlarında temel bir desen, analizi genomik bölgelere dağıtmaktır. Örneğin, varyant çağırma, genomu aralıklara (kromozomlar veya daha küçük bölgeler gibi) bölerek paralelleştirilebilir. Bu paralelleştirme stratejisi, hesaplama yükünü birden fazla çekirdeğe veya düğüme dağıtarak pipeline verimliliğini önemli ölçüde artırır ve toplam yürütme süresini azaltır.

Aşağıdaki bölümde, örnek verilerimizi birden fazla genomik aralığa nasıl dağıtacağımızı göstereceğiz. Her örneği her aralıkla eşleştireceğiz, bu da farklı genomik bölgelerin paralel işlenmesine olanak tanır. Bu, veri kümemizin boyutunu aralık sayısıyla çarpacak ve daha sonra bir araya getirilebilecek birden fazla bağımsız analiz birimi oluşturacaktır.

### 4.1. `combine` kullanarak örnekleri aralıklara yayma

Bir aralık channel'ı oluşturarak başlayalım. İşleri basit tutmak için, manuel olarak tanımlayacağımız sadece 3 aralık kullanacağız. Gerçek bir iş akışında, bunları bir dosya girdisinden okuyabilir veya hatta çok sayıda aralık dosyasıyla bir channel oluşturabilirsiniz.

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

Şimdi hatırlayın, her örneği her aralık için tekrarlamak istiyoruz. Bu bazen örneklerin ve aralıkların Kartezyen çarpımı olarak adlandırılır. Bunu [`combine` operatörünü](https://www.nextflow.io/docs/latest/operator.html#combine) kullanarak başarabiliriz. Bu, channel 1'den her öğeyi alacak ve channel 2'deki her öğe için tekrarlayacaktır. İş akışımıza bir combine operatörü ekleyelim:

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

Başarılı! Her örneği 3 aralıklık listemizdeki her bir aralık için tekrarladık. Channel'ımızdaki öğe sayısını etkili bir şekilde üçe katladık.

Okumak biraz zor olsa da, bir sonraki bölümde düzenleyeceğiz.

### 4.2. Channel'ı düzenleme

Örnek verilerimizi daha kolay anlaşılır hale getirmek için `map` operatörünü kullanarak düzenleyip yeniden yapılandırabiliriz. Aralık dizesini ilk öğedeki birleştirme map'ine taşıyalım.

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

İlk olarak, kodu daha okunabilir hale getirmek için adlandırılmış parametreler kullanıyoruz. `grouping_key`, `normal`, `tumor` ve `interval` adlarını kullanarak, tuple'daki öğelere dizin yerine adla başvurabiliriz:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

Ardından, `grouping_key`'i `interval` alanıyla birleştiriyoruz. `grouping_key`, `id` ve `repeat` alanlarını içeren bir map'tir. `interval` ile yeni bir map oluşturuyoruz ve Groovy'nin map toplama (`+`) işlemini kullanarak bunları birleştiriyoruz:

```groovy
                grouping_key + [interval: interval],
```

Son olarak, bunu üç öğeli bir tuple olarak döndürüyoruz: birleştirilmiş meta veri map'i, normal örnek dosyası ve tümör örnek dosyası:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Tekrar çalıştıralım ve channel içeriklerini kontrol edelim:

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

Artık her örneği tüm genomik aralıklarda tekrarladık, paralel olarak işlenebilecek birden fazla bağımsız analiz birimi oluşturduk. Peki ya ilişkili örnekleri tekrar bir araya getirmek istersek? Bir sonraki bölümde, ortak özellikleri paylaşan örnekleri nasıl gruplayacağımızı öğreneceğiz.

### Özet

Bu bölümde öğrendikleriniz:

- **Örnekleri aralıklara yayma**: Örnekleri aralıklara yaymak için `combine` nasıl kullanılır
- **Kartezyen çarpımlar oluşturma**: Örneklerin ve aralıkların tüm kombinasyonları nasıl oluşturulur
- **Channel yapısını düzenleme**: Daha iyi okunabilirlik için verileri yeniden yapılandırmak üzere `map` nasıl kullanılır
- **Paralel işleme hazırlığı**: Dağıtık analiz için veriler nasıl ayarlanır

## 5. `groupTuple` kullanarak örnekleri toplama

Önceki bölümlerde, bir girdi dosyasından verileri bölmeyi ve belirli alanlara göre filtrelemeyi (bizim durumumuzda normal ve tümör örnekleri) öğrendik. Ancak bu yalnızca tek bir birleştirme türünü kapsar. Ya örnekleri belirli bir özelliğe göre gruplamak istersek? Örneğin, eşleşen normal-tümör çiftlerini birleştirmek yerine, türlerine bakılmaksızın "sampleA"dan tüm örnekleri birlikte işlemek isteyebiliriz. Bu desen, verimlilik nedenleriyle ilişkili örnekleri ayrı ayrı işlemek ve sonunda sonuçları karşılaştırmak veya birleştirmek istediğiniz biyoinformatik iş akışlarında yaygındır.

Nextflow bunu yapmak için yerleşik yöntemler içerir, bakacağımız ana yöntem `groupTuple`'dır.

Aynı `id` ve `interval` alanlarına sahip tüm örneklerimizi gruplayarak başlayalım, bu teknik replikleri gruplamak istediğimiz ancak anlamlı olarak farklı örnekleri ayrı tutmak istediğimiz bir analizin tipik örneği olur.

Bunu yapmak için, gruplama değişkenlerimizi izole etmek için ayırmalıyız.

İlk adım, önceki bölümde yaptığımıza benzer. Gruplama değişkenimizi tuple'ın ilk öğesi olarak izole etmeliyiz. Hatırlayın, ilk öğemiz şu anda `id`, `repeat` ve `interval` alanlarının bir map'i:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

`id` ve `interval` alanlarımızı map'ten izole etmek için daha önceki `subMap` metodunu yeniden kullanabiliriz. Daha önce olduğu gibi, her örnek için tuple'ın ilk öğesine `subMap` metodunu uygulamak için `map` operatörünü kullanacağız.

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

Channel içeriklerini kontrol etmek için tekrar çalıştıralım:

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

`id` ve `interval` alanlarını başarıyla izole ettiğimizi ancak örnekleri henüz gruplamadığımızı görebiliriz.

!!! note "Not"

    Burada `replicate` alanını atıyoruz. Bunun nedeni, daha aşağı akış işleme için buna ihtiyacımız olmamasıdır. Bu eğitimi tamamladıktan sonra, sonraki gruplamayı etkilemeden bunu dahil edip edemeyeceğinize bakın!

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

Hepsi bu kadar! Sadece tek bir satır kod ekledik. Çalıştırdığımızda ne olacağını görelim:

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

Verilerimizin yapısının değiştiğini ve her channel öğesinde dosyaların artık `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]` gibi tuple'lar içinde yer aldığını not edin. Bunun nedeni, `groupTuple` kullandığımızda, Nextflow'un bir grubun her örneği için tekil dosyaları birleştirmesidir. Aşağı akışta verileri işlemeye çalışırken bunu hatırlamak önemlidir.

!!! note "Not"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) groupTuple'ın tersidir. Bir channel'daki öğeleri açar ve düzleştirir. Yukarıda gerçekleştirdiğimiz gruplamayı geri almak için `transpose` eklemeyi deneyin!

### Özet

Bu bölümde öğrendikleriniz:

- **Gruplama anahtarlarını izole etme**: Gruplama için belirli alanları çıkarmak üzere `subMap` nasıl kullanılır
- **Gruplanmış veri yapılarını işleme**: `groupTuple` tarafından oluşturulan iç içe yapıyla nasıl çalışılır
- **Teknik replik işleme**: Aynı deneysel koşulları paylaşan örnekler nasıl gruplandırılır

---

## Özet

Bu yan görevde, channel'ları kullanarak verileri bölmeyi ve gruplamayı öğrendiniz.

Verileri pipeline boyunca akarken değiştirerek, döngüler veya while ifadeleri kullanmadan ölçeklenebilir bir pipeline oluşturabilirsiniz; bu, daha geleneksel yaklaşımlara göre çeşitli avantajlar sunar:

- Ek kod olmadan istediğimiz kadar çok veya az girdiyle ölçeklenebiliriz
- Yineleme yerine verilerin pipeline boyunca akışını yönetmeye odaklanıyoruz
- Gerektiği kadar karmaşık veya basit olabiliriz
- Pipeline daha bildirimsel hale gelir, nasıl yapılacağı yerine ne yapılması gerektiğine odaklanır
- Nextflow, bağımsız işlemleri paralel olarak çalıştırarak yürütmeyi bizim için optimize eder

Bu channel işlemlerinde ustalaşmak, karmaşık veri ilişkilerini döngülere veya yinelemeli programlamaya başvurmadan ele alan esnek, ölçeklenebilir pipeline'lar oluşturmanıza olanak tanır ve Nextflow'un yürütmeyi optimize etmesine ve bağımsız işlemleri otomatik olarak paralelleştirmesine izin verir.

### Temel desenler

1.  **Yapılandırılmış girdi verileri oluşturma:** Meta map'lerle bir CSV dosyasından başlama ([İş akışlarında meta veri](./metadata.md)'den desenlere dayanarak)

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Verileri ayrı channel'lara bölme:** `type` alanına göre verileri bağımsız akışlara bölmek için `filter` kullandık

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Eşleşen örnekleri birleştirme:** `id` ve `repeat` alanlarına göre ilişkili örnekleri yeniden birleştirmek için `join` kullandık

    - İki channel'ı anahtara göre birleştirme (tuple'ın ilk öğesi)

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

    - subMap kullanarak birden fazla alanda birleştirme

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Aralıklara dağıtma:** Paralel işleme için örneklerin genomik aralıklarla Kartezyen çarpımlarını oluşturmak üzere `combine` kullandık.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Gruplama anahtarlarına göre toplama:** Her tuple'daki ilk öğeye göre gruplamak için `groupTuple` kullandık, böylece `id` ve `interval` alanlarını paylaşan örnekleri topladık ve teknik replikleri birleştirdik.

    ```groovy
    channel.groupTuple()
    ```

6.  **Veri yapısını optimize etme:** Belirli alanları çıkarmak için `subMap` kullandık ve dönüşümleri yeniden kullanılabilir hale getirmek için adlandırılmış bir closure oluşturduk.

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

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
