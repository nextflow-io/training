# Geliştirme Ortamı

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Modern Entegre Geliştirme Ortamları (IDE'ler) Nextflow geliştirme deneyiminizi önemli ölçüde dönüştürebilir. Bu yan görev, özellikle VS Code ve Nextflow eklentisini kullanarak daha hızlı kod yazma, hataları erken yakalama ve karmaşık iş akışlarında verimli gezinme konularına odaklanır.

!!! note "Bu geleneksel bir eğitim değildir"

    Diğer eğitim modüllerinin aksine, bu kılavuz adım adım bir eğitimden ziyade hızlı ipuçları, tavsiyeler ve pratik örnekler koleksiyonu olarak düzenlenmiştir. Her bölüm, mevcut geliştirme ihtiyaçlarınıza ve ilgi alanlarınıza göre bağımsız olarak keşfedilebilir. İstediğiniz gibi gezinin ve iş akışı geliştirmeniz için en faydalı olacak özelliklere odaklanın.

## Önce bilmeniz gerekenler

Bu kılavuz, [Hello Nextflow](../hello_nextflow/) eğitim kursunu tamamladığınızı ve aşağıdakiler dahil temel Nextflow kavramlarında rahat olduğunuzu varsayar:

- **Temel iş akışı yapısı**: Süreçleri, iş akışlarını ve bunların nasıl birbirine bağlandığını anlama
- **Kanal işlemleri**: Kanal oluşturma, süreçler arasında veri aktarma ve temel operatörleri kullanma
- **Modüller ve organizasyon**: Yeniden kullanılabilir modüller oluşturma ve include ifadelerini kullanma
- **Yapılandırma temelleri**: Parametreler, süreç yönergeleri ve profiller için `nextflow.config` kullanma

## Burada neler öğreneceksiniz

Bu kılavuz, sizi daha verimli bir Nextflow geliştiricisi yapacak **IDE üretkenlik özelliklerine** odaklanır:

- **Gelişmiş sözdizimi vurgulama**: VS Code'un kod yapınız hakkında size ne gösterdiğini anlama
- **Akıllı otomatik tamamlama**: Daha hızlı kod yazımı için içeriğe duyarlı önerilerden yararlanma
- **Hata algılama ve tanılama**: İş akışınızı çalıştırmadan önce sözdizimi hatalarını yakalama
- **Kod gezinmesi**: Süreçler, modüller ve tanımlar arasında hızlıca geçiş yapma
- **Biçimlendirme ve organizasyon**: Tutarlı, okunabilir kod stilini koruma
- **Yapay zeka destekli geliştirme** (opsiyonel): IDE'nizle entegre modern yapay zeka araçlarını kullanma

!!! info "Neden şimdi IDE özellikleri?"

    [Hello Nextflow](../hello_nextflow/) kursu sırasında muhtemelen VS Code'u zaten kullanıyordunuz, ancak odağı IDE özelliklerinden ziyade Nextflow temellerini öğrenmeye yönlendirdik. Artık süreçler, iş akışları, kanallar ve modüller gibi temel Nextflow kavramlarında rahatsınız, sizi daha verimli bir geliştirici yapacak gelişmiş IDE özelliklerinden yararlanmaya hazırsınız.

    Bunu geliştirme ortamınızı "seviye atlama" olarak düşünün - kullandığınız aynı editör, ne konusunda size yardımcı olduklarını anladıktan sonra gerçekten değerli hale gelen çok daha güçlü yeteneklere sahip.

---

## 0. Kurulum ve Isınma

IDE özelliklerini keşfetmek için özel olarak bir çalışma alanı oluşturalım:

```bash title="IDE özellikleri dizinine gidin"
cd side-quests/ide_features
```

Bu dizini VS Code'da açın:

```bash title="Geçerli dizinde VS Code'u açın"
code .
```

`ide_features` dizini çeşitli IDE özelliklerini gösteren örnek iş akışlarını içerir:

```bash title="Dizin yapısını göster"
tree .
```

```console title="Proje yapısı"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "Örnek Dosyalar Hakkında"

    - `basic_workflow.nf` çalıştırabileceğiniz ve değiştirebileceğiniz çalışan temel bir iş akışıdır
    - `complex_workflow.nf` sadece gösterim amaçlı tasarlanmıştır ve gezinme özelliklerini gösterir - başarıyla çalışmayabilir ancak gerçekçi çok dosyalı iş akışı yapısını gösterir

### Klavye Kısayolları

Bu kılavuzdaki bazı özellikler opsiyonel klavye kısayollarını kullanacaktır. Bu materyale tarayıcı üzerinden GitHub Codespaces aracılığıyla erişiyor olabilirsiniz ve bu durumda bazen kısayollar beklendiği gibi çalışmayabilir çünkü sisteminizde başka şeyler için kullanılıyorlar.

VS Code'u yerel olarak çalıştırıyorsanız, muhtemelen gerçekte iş akışları yazarken olacağınız gibi, kısayollar açıklandığı gibi çalışacaktır.

Mac kullanıyorsanız, bazı (hepsi değil) klavye kısayolları "ctrl" yerine "cmd" kullanacaktır ve bunu metinde `Ctrl/Cmd` şeklinde belirteceğiz.

### 0.1. Nextflow Eklentisini Kurma

!!! note "Zaten Devcontainer Kullanıyor musunuz?"

    **GitHub Codespaces**'te çalışıyorsanız veya **yerel bir devcontainer** kullanıyorsanız, Nextflow eklentisi muhtemelen sizin için zaten kurulu ve yapılandırılmış durumdadır. Aşağıdaki manuel kurulum adımlarını atlayabilir ve doğrudan eklenti özelliklerini keşfetmeye devam edebilirsiniz.

Eklentiyi manuel olarak kurmak için:

1. VS Code'u açın
2. Sol taraftaki eklentiler simgesine tıklayarak Extensions görünümüne gidin: ![eklentiler simgesi](img/extensions_icon.png) (VSCode'u yerel olarak çalıştırıyorsanız kısayol `Ctrl/Cmd+Shift+X`)
3. "Nextflow" arayın
4. Resmi Nextflow eklentisini kurun

![Nextflow Eklentisini Kur](img/install_extension.png)

### 0.2. Çalışma Alanı Düzeni

Hello Nextflow boyunca VS Code kullandığınız için temellere zaten aşinasınız. Bu oturum için çalışma alanınızı verimli bir şekilde nasıl düzenleyeceğiniz aşağıda açıklanmıştır:

- **Editör Alanı**: Dosyaları görüntülemek ve düzenlemek için. Dosyaları yan yana karşılaştırmak için bunu birden fazla bölmeye ayırabilirsiniz.
- **Dosya Gezgini** tıklayın (![dosya gezgini simgesi](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Sisteminizdeki yerel dosyalar ve klasörler. Dosyalar arasında gezinmek için bunu solda açık tutun
- **Entegre Terminal** (`Ctrl+Shift+` ters tırnak hem Windows hem de MacOS için): Alt kısımda bilgisayarla etkileşim için bir terminal. Nextflow veya diğer komutları çalıştırmak için bunu kullanın.
- **Sorunlar Paneli** (`Ctrl+Shift+M`): VS Code algıladığı hataları ve sorunları burada gösterecektir. Bu, sorunları bir bakışta vurgulamak için kullanışlıdır.

Örnekler üzerinde çalışırken düzeninizi özelleştirmek için panelleri sürükleyebilir veya gizleyebilirsiniz (kenar çubuğunu açıp kapatmak için `Ctrl/Cmd+B`).

### Çıkarım

VS Code'u Nextflow eklentisi ile kurdunuz ve verimli geliştirme için çalışma alanı düzenini anlıyorsunuz.

### Sırada ne var?

Sözdizimi vurgulamanın Nextflow kod yapısını bir bakışta anlamanıza nasıl yardımcı olduğunu öğrenin.

---

## 1. Sözdizimi Vurgulama ve Kod Yapısı

Artık çalışma alanınız hazır olduğuna göre, VS Code'un sözdizimi vurgulamasının Nextflow kodunu daha etkili bir şekilde okumanıza ve yazmanıza nasıl yardımcı olduğunu keşfedelim.

### 1.1. Nextflow Sözdizimi Öğeleri

Sözdizimi vurgulamayı çalışırken görmek için `basic_workflow.nf`'i açın:

![Sözdizimi Vitrini](img/syntax_showcase.png)

VS Code'un nasıl vurguladığına dikkat edin:

- **Anahtar kelimeler** (`process`, `workflow`, `input`, `output`, `script`) farklı renklerde
- **Dize değişmezleri** ve **parametreler** farklı stillerde
- **Yorumlar** sönük bir renkte
- **Değişkenler** ve **fonksiyon çağrıları** uygun vurguyla
- **Kod blokları** uygun girinti kılavuzlarıyla

!!! note "Temaya Bağlı Renkler"

    Gördüğünüz belirli renkler VS Code temanıza (karanlık/açık mod), renk ayarlarınıza ve yaptığınız özelleştirmelere bağlı olacaktır. Önemli olan, kod yapısının seçtiğiniz renk şemasından bağımsız olarak anlaşılmasını kolaylaştırmak için farklı sözdizimi öğelerinin görsel olarak birbirinden ayırt edilmesidir.

### 1.2. Kod Yapısını Anlama

Sözdizimi vurgulaması size şunları hızlıca belirlemenize yardımcı olur:

- **Süreç sınırları**: Farklı süreçler arasında net ayrım
- **Girdi/çıktı blokları**: Veri akışı tanımlarını bulmayı kolaylaştırır
- **Script blokları**: Yürütülen gerçek komutlar
- **Kanal işlemleri**: Veri dönüştürme adımları
- **Yapılandırma yönergeleri**: Sürece özgü ayarlar

Bu görsel organizasyon, birden fazla süreç ve karmaşık veri akışları içeren karmaşık iş akışlarıyla çalışırken paha biçilmez hale gelir.

### Çıkarım

VS Code'un sözdizimi vurgulamasının Nextflow kod yapısını okumanıza ve daha hızlı geliştirme için farklı dil öğelerini belirlemenize nasıl yardımcı olduğunu anlıyorsunuz.

### Sırada ne var?

Akıllı otomatik tamamlamanın içeriğe duyarlı önerilerle kod yazmayı nasıl hızlandırdığını öğrenin.

---

## 2. Akıllı Otomatik Tamamlama

VS Code'un otomatik tamamlama özellikleri, içeriğe göre uygun seçenekler önererek daha hızlı ve daha az hatayla kod yazmanıza yardımcı olur.

### 2.1. İçeriğe Duyarlı Öneriler

Otomatik tamamlama seçenekleri kodunuzda nerede olduğunuza bağlı olarak değişir:

#### Kanal İşlemleri

`basic_workflow.nf`'i tekrar açın ve workflow bloğunda `channel.` yazmayı deneyin:

![Kanal otomatik tamamlama](img/autocomplete_channel.png)

Şunlar için öneriler göreceksiniz:

- `fromPath()` - Dosya yollarından kanal oluştur
- `fromFilePairs()` - Eşleştirilmiş dosyalardan kanal oluştur
- `of()` - Değerlerden kanal oluştur
- `fromSRA()` - SRA erişim numaralarından kanal oluştur
- Ve daha fazlası...

Bu, tam metod adlarını hatırlamak zorunda kalmadan kullanılacak doğru kanal fabrikasını hızlıca bulmanıza yardımcı olur.

Ayrıca kanallara uygulanabilecek operatörleri de keşfedebilirsiniz. Örneğin, mevcut işlemleri görmek için `FASTQC.out.html.` yazın:

![Kanal işlemleri otomatik tamamlama](img/autocomplete_operators.png)

#### Süreç Yönergeleri

Bir süreç script bloğunun içinde, mevcut çalışma zamanı özelliklerini görmek için `task.` yazın:

![Görev özellikleri otomatik tamamlama](img/autocomplete_task.png)

#### Yapılandırma

nextflow.config dosyasını açın ve mevcut süreç yönergelerini görmek için herhangi bir yere `process.` yazın:

![Yapılandırma otomatik tamamlama](img/autocomplete_config.png)

Şunlar için öneriler göreceksiniz:

- `executor`
- `memory`
- `cpus`

Bu, süreçleri yapılandırırken zamandan tasarruf sağlar ve farklı yapılandırma kapsamlarında çalışır. Örneğin, Docker'a özgü yapılandırma seçeneklerini görmek için `docker.` yazmayı deneyin.

### Çıkarım

Sözdizimini ezberlemeden mevcut kanal işlemlerini, süreç yönergelerini ve yapılandırma seçeneklerini keşfetmek için VS Code'un akıllı otomatik tamamlamasını kullanabilirsiniz.

### Sırada ne var?

Gerçek zamanlı hata algılamanın, iş akışınızı çalıştırmadan önce sorunları sadece kodu okuyarak yakalamanıza nasıl yardımcı olduğunu öğrenin.

## 3. Hata Algılama ve Tanılama

VS Code'un gerçek zamanlı hata algılaması, iş akışınızı çalıştırmadan önce sorunları yakalamanıza yardımcı olur.

### 3.1. Sözdizimi Hatası Algılama

Algılamayı çalışırken görmek için kasıtlı bir hata oluşturalım. `basic_workflow.nf`'i açın ve süreç adını `FASTQC`'den `FASTQ`'ya (veya başka geçersiz bir ada) değiştirin. VS Code hatayı workflow bloğunda kırmızı dalgalı alt çizgiyle hemen vurgulayacaktır:

![Hata alt çizgisi](img/error_underline.png)

### 3.2. Sorunlar Paneli

Bireysel hata vurgulamanın ötesinde, VS Code çalışma alanınızdaki tüm hataları, uyarıları ve bilgi mesajlarını toplayan merkezi bir Sorunlar paneli sağlar. `Ctrl/Cmd+Shift+M` ile açın ve yalnızca geçerli dosyayla ilgili hataları göstermek için filtre simgesini kullanın:

![Sorunlar panelini filtrele](img/active_file.png)

Sorunlu satıra doğrudan atlamak için herhangi bir soruna tıklayın

![Sorunlar Paneli](img/problems_panel.png)

Süreç adını tekrar `FASTQC` olarak değiştirerek hatayı düzeltin.

### 3.3. Yaygın Hata Kalıpları

Nextflow sözdizimindeki yaygın hatalar şunlardır:

- **Eksik parantezler**: Eşleşmeyen `{` veya `}`
- **Eksik bloklar**: Süreçlerde gerekli bölümlerin eksik olması
- **Geçersiz sözdizimi**: Hatalı biçimlendirilmiş Nextflow DSL
- **Anahtar kelimelerde yazım hataları**: Yanlış yazılmış süreç yönergeleri
- **Kanal uyumsuzlukları**: Tür uyumsuzlukları

Nextflow dil sunucusu bu sorunları Sorunlar panelinde vurgular. Bir boru hattı çalıştırırken sözdizimi hatalarından kaçınmak için bunları erken kontrol edebilirsiniz.

### Çıkarım

İş akışınızı çalıştırmadan önce sözdizimi hatalarını ve sorunları yakalamak için VS Code'un hata algılamasını ve Sorunlar panelini kullanabilir, zamandan tasarruf edebilir ve hayal kırıklığını önleyebilirsiniz.

### Sırada ne var?

Karmaşık iş akışlarında süreçler, modüller ve tanımlar arasında nasıl verimli gezineceğinizi öğrenin.

---

## 4. Kod Gezinmesi ve Sembol Yönetimi

Birden fazla dosyaya yayılan karmaşık iş akışlarıyla çalışırken verimli gezinme çok önemlidir. Bunu anlamak için, `basic_workflow.nf`'deki süreç tanımını size sağladığımız modül için bir import ile değiştirin:

=== "Sonra"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Önce"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Tanıma Git

`FASTQC` gibi bir süreç adının üzerine fareyle gelirseniz, modül arayüzünü (girdiler ve çıktılar) içeren bir açılır pencere göreceksiniz:

![Tanıma git](img/syntax.png)

Bu özellik, modül dosyasını doğrudan açmadan modül arayüzünü anlamanıza olanak tanıdığı için iş akışları yazarken özellikle değerlidir.

**Ctrl/Cmd-tıklama** kullanarak herhangi bir süreç, modül veya değişken tanımına hızlıca gidebilirsiniz. Betiğin üst kısmındaki modül dosyasına bağlantı üzerine fareyle gelin ve önerildiği gibi bağlantıyı takip edin:

![Bağlantıyı takip et](img/follow_link.png)

Aynı şey süreç adları için de çalışır. `basic_workflow.nf`'e geri dönün ve bunu workflow bloğundaki `FASTQC` süreç adında deneyin. Bu sizi doğrudan süreç adına bağlar (bu örnekte modül dosyasıyla aynıdır, ancak çok daha büyük bir dosyanın ortasında olabilir).

Bulunduğunuz yere geri dönmek için **Alt+←** (veya Mac'te **Ctrl+-**) kullanın. Bu, yerinizi kaybetmeden kodu keşfetmenin güçlü bir yoludur.

Şimdi `complex_workflow.nf`'i (daha önce bahsedilen sadece gösterim amaçlı dosya) kullanarak daha karmaşık bir iş akışında gezinmeyi keşfedelim. Bu iş akışı, ayrı modül dosyalarında tanımlanan birden fazla süreç ve bazı satır içi süreçler içerir. Karmaşık çok dosyalı yapılar manuel olarak gezinmek zor olabilse de, tanımlara atlama yeteneği keşfetmeyi çok daha yönetilebilir hale getirir.

1. `complex_workflow.nf`'i açın
2. Modül tanımlarına gidin
3. Geri gitmek için **Alt+←** (veya **Ctrl+-**) kullanın
4. Workflow bloğundaki `FASTQC` süreç adına gidin. Bu sizi doğrudan süreç adına bağlar (bu örnekte modül dosyasıyla aynıdır, ancak çok daha büyük bir dosyanın ortasında olabilir).
5. Tekrar geri gidin
6. Workflow bloğundaki `TRIM_GALORE` sürecine gidin. Bu satır içi tanımlanmıştır, bu nedenle sizi ayrı bir dosyaya götürmez, ancak yine de süreç tanımını gösterir ve bulunduğunuz yere geri gidebilirsiniz.

### 4.2. Sembol Gezinmesi

`complex_workflow.nf` hala açıkken, VSCode'un üst kısmındaki arama çubuğuna `@` yazarak dosyadaki tüm sembollerin genel görünümünü elde edebilirsiniz (klavye kısayolu `Ctrl/Cmd+Shift+O`'dur, ancak Codespaces'te çalışmayabilir). Bu, geçerli dosyadaki tüm sembolleri listeleyen sembol gezinme panelini açar:

![Sembol gezinmesi](img/symbols.png)

Bu şunları gösterir:

- Tüm süreç tanımları
- Workflow tanımları (bu dosyada iki workflow tanımlanmıştır)
- Fonksiyon tanımları

Sonuçları filtrelemek için yazmaya başlayın.

### 4.3. Tüm Referansları Bul

Bir sürecin veya değişkenin kod tabanınız boyunca nerede kullanıldığını anlamak çok yardımcı olabilir. Örneğin, `FASTQC` sürecine yapılan tüm referansları bulmak istiyorsanız, önce tanımına giderek başlayın. Bunu doğrudan `modules/fastqc.nf`'i açarak veya yukarıda yaptığımız gibi VS Code'un `Ctrl/Cmd-tıklama` ile hızlı gezinme özelliğini kullanarak yapabilirsiniz. Süreç tanımına ulaştıktan sonra, `FASTQC` süreç adına sağ tıklayın ve kullanıldığı tüm örnekleri görmek için bağlam menüsünden "Find All References"ı seçin.

![Referansları bul](img/references.png)

Bu özellik, `FASTQC`'nin çalışma alanınızda referans verildiği tüm örnekleri, iki farklı iş akışındaki kullanımı dahil olmak üzere gösterir. Bu içgörü, `FASTQC` sürecindeki değişikliklerin potansiyel etkisini değerlendirmek için çok önemlidir.

### 4.4. Anahat Paneli

Explorer kenar çubuğunda bulunan (![Explorer simgesi](img/files_icon.png) tıklayın) Outline paneli, geçerli dosyanızdaki tüm sembollerin kullanışlı bir genel görünümünü sağlar. Bu özellik, fonksiyonları, değişkenleri ve diğer anahtar öğeleri hiyerarşik bir görünümde göstererek kodunuzun yapısında hızlıca gezinmenize ve yönetmenize olanak tanır.

![Anahat paneli](img/outline.png)

Dosya tarayıcısını kullanmadan kodunuzun farklı bölümlerine hızlıca gitmek için Outline panelini kullanın.

### 4.5. DAG görselleştirme

VS Code'un Nextflow eklentisi iş akışınızı Yönlendirilmiş Döngüsüz Grafik (DAG) olarak görselleştirebilir. Bu, veri akışını ve süreçler arasındaki bağımlılıkları anlamanıza yardımcı olur. `complex_workflow.nf`'i açın ve `workflow {`'ın üzerindeki "Preview DAG" düğmesine tıklayın (bu dosyadaki ikinci `workflow` bloğu):

![DAG önizlemesi](img/dag_preview.png)

Bu sadece 'giriş' iş akışıdır, ancak yukarıdaki `workflow RNASEQ_PIPELINE {`'nin üzerindeki "Preview DAG" düğmesine tıklayarak iç iş akışları için DAG'yi de önizleyebilirsiniz:

![DAG önizlemesi iç iş akışı](img/dag_preview_inner.png)

Bu iş akışı için, koddaki ilgili süreç tanımlarına gitmek için DAG'deki düğümleri kullanabilirsiniz. Bir düğüme tıklayın ve sizi editörde ilgili süreç tanımına götürecektir. Özellikle bir iş akışı büyük bir boyuta ulaştığında, bu gerçekten kod etrafında gezinmenize ve süreçlerin nasıl bağlandığını anlamanıza yardımcı olabilir.

### Çıkarım

Kod yapısını ve bağımlılıkları anlamak için tanıma git, sembol arama, referansları bul ve DAG görselleştirme özelliklerini kullanarak karmaşık iş akışlarında verimli gezinebilirsiniz.

### Sırada ne var?

Daha büyük Nextflow projelerinde birbirine bağlı birden fazla dosyayla nasıl etkili çalışacağınızı öğrenin.

## 5. Birden Fazla Dosyada Çalışma

Gerçek Nextflow geliştirme, birbirine bağlı birden fazla dosyayla çalışmayı içerir. VS Code'un karmaşık projeleri verimli bir şekilde yönetmenize nasıl yardımcı olduğunu keşfedelim.

### 5.1. Hızlı Dosya Gezinmesi

`complex_workflow.nf` açıkken, birkaç modülü içe aktardığını fark edeceksiniz. Aralarında hızlı gezinme pratiği yapalım.

**Ctrl+P** (veya **Cmd+P**) tuşuna basın ve "fast" yazmaya başlayın:

VS Code size eşleşen dosyaları gösterecektir. Anında atlamak için `modules/fastqc.nf`'i seçin. Bu, hangi dosyayı aradığınızı kabaca bildiğinizde dosya gezgini arasında tıklamaktan çok daha hızlıdır.

Bunu diğer kalıplarla deneyin:

- STAR hizalama modül dosyasını (`star.nf`) bulmak için "star" yazın
- Yardımcı fonksiyonlar dosyasını (`utils.nf`) bulmak için "utils" yazın
- Yapılandırma dosyalarına (`nextflow.config`) atlamak için "config" yazın

### 5.2. Çok Dosyalı Geliştirme için Bölünmüş Editör

Modüllerle çalışırken, genellikle hem ana iş akışını hem de modül tanımlarını aynı anda görmeniz gerekir. Bunu kuralım:

1. `complex_workflow.nf`'i açın
2. Yeni bir sekmede `modules/fastqc.nf`'i açın
3. `modules/fastqc.nf` sekmesine sağ tıklayın ve "Split Right" seçin
4. Şimdi her iki dosyayı yan yana görebilirsiniz

![Bölünmüş editör](img/split_editor.png)

Bu şu durumlarda paha biçilmezdir:

- Workflow çağrılarını yazarken modül arayüzlerini kontrol etme ve önizleme yeterli değil
- Farklı modüller arasında benzer süreçleri karşılaştırma
- Workflow ve modüller arasındaki veri akışında hata ayıklama

### 5.3. Proje Genelinde Arama

Bazen belirli kalıpların tüm projenizde nerede kullanıldığını bulmanız gerekir. Arama panelini açmak için `Ctrl/Cmd+Shift+F` tuşlarına basın.

Çalışma alanı genelinde `publishDir` aramayı deneyin:

![Proje araması](img/project_search.png)

Bu size yayınlama dizinlerini kullanan her dosyayı gösterir ve şunlarda yardımcı olur:

- Çıktı organizasyon kalıplarını anlama
- Belirli yönergelerin örneklerini bulma
- Modüller arasında tutarlılık sağlama

### Çıkarım

Hızlı dosya gezinmesi, bölünmüş editörler ve proje genelinde arama kullanarak karmaşık çok dosyalı projeleri yönetebilir ve iş akışları ve modüller arasında verimli çalışabilirsiniz.

### Sırada ne var?

Kod biçimlendirme ve bakım özelliklerinin iş akışlarınızı nasıl organize ve okunabilir tuttuğunu öğrenin.

---

## 6. Kod Biçimlendirme ve Bakım

Uygun kod biçimlendirmesi sadece estetik için değil, aynı zamanda okunabilirliği, anlaşılabilirliği ve karmaşık iş akışlarını güncelleme kolaylığını artırmak için de gereklidir.

### 6.1. Otomatik Biçimlendirme Çalışırken

`basic_workflow.nf`'i açın ve kasıtlı olarak biçimlendirmeyi bozun:

- Bazı girintileri kaldırın: Tüm belgeyi vurgulayın ve mümkün olduğunca çok girintiyi kaldırmak için `shift+tab` tuşlarına birçok kez basın.
- Rastgele yerlere ekstra boşluklar ekleyin: `channel.fromPath` ifadesinde, `(` sonrasına 30 boşluk ekleyin.
- Bazı satırları garip şekilde bölün: `.view {` operatörü ile `Processing sample:` dizesi arasına yeni bir satır ekleyin ancak kapanış parantezi `}` öncesine karşılık gelen bir yeni satır eklemeyin.

Şimdi otomatik biçimlendirme için `Shift+Alt+F` (veya MacOS'ta `Shift+Option+F`) tuşlarına basın:

VS Code hemen:

- Süreç yapısını net göstermek için girintiyi düzeltir
- Benzer öğeleri tutarlı bir şekilde hizalar
- Gereksiz boşlukları kaldırır
- Okunabilir satır sonlarını korur

Otomatik biçimlendirmenin her kod stili sorununu çözemeyebileceğini unutmayın. Nextflow dil sunucusu kodunuzu düzenli tutmayı amaçlar, ancak belirli alanlarda kişisel tercihlerinize de saygı gösterir. Örneğin, bir sürecin `script` bloğu içindeki girintiyi kaldırırsanız, biçimlendirici olduğu gibi bırakacaktır, çünkü kasıtlı olarak o stili tercih ediyor olabilirsiniz.

Şu anda Nextflow için katı bir stil zorlaması yoktur, bu nedenle dil sunucusu biraz esneklik sunar. Ancak, netliği korumak için metod ve fonksiyon tanımları etrafında biçimlendirme kurallarını tutarlı bir şekilde uygulayacaktır.

### 6.2. Kod Organizasyonu Özellikleri

#### Hızlı Yorum Yapma

İş akışınızda bir kod bloğu seçin ve yorum satırı haline getirmek için **Ctrl+/** (veya **Cmd+/**) tuşlarına basın:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Bu şu durumlarda mükemmeldir:

- Geliştirme sırasında iş akışlarının bölümlerini geçici olarak devre dışı bırakma
- Karmaşık kanal işlemlerine açıklayıcı yorumlar ekleme
- İş akışı bölümlerini belgeleme

Kodu yorumsuz hale getirmek için tekrar **Ctrl+/** (veya **Cmd+/**) kullanın.

#### Genel Bakış için Kod Katlama

`complex_workflow.nf`'de, süreç tanımlarının yanındaki küçük oklara dikkat edin. Süreçleri katlamak (daraltmak) için onlara tıklayın:

![Kod katlama](img/code_folding.png)

Bu, uygulama detaylarında kaybolmadan iş akışı yapınızın üst düzey bir genel görünümünü sağlar.

#### Parantez Eşleştirme

İmlecinizi herhangi bir `{` veya `}` parantezinin yanına yerleştirin ve VS Code eşleşen parantezi vurgular. Eşleşen parantezler arasında atlamak için **Ctrl+Shift+\\** (veya **Cmd+Shift+\\**) kullanın.

Bu şunlar için çok önemlidir:

- Süreç sınırlarını anlama
- Eksik veya fazla parantezleri bulma
- İç içe geçmiş iş akışı yapılarında gezinme

#### Çok Satırlı Seçim ve Düzenleme

Birden fazla satırı aynı anda düzenlemek için VS Code güçlü çok imleç yetenekleri sunar:

- **Çok satırlı seçim**: **Ctrl+Alt** (veya MacOS için **Cmd+Option**) basılı tutun ve birden fazla satır seçmek için ok tuşlarını kullanın
- **Çok satırlı girinti**: Birden fazla satır seçin ve tüm blokları girintili hale getirmek için **Tab** veya girintiyi kaldırmak için **Shift+Tab** kullanın

Bu özellikle şunlar için kullanışlıdır:

- Tüm süreç bloklarını tutarlı bir şekilde girintileme
- Aynı anda birden fazla satıra yorum ekleme
- Birden fazla süreçte benzer parametre tanımlarını düzenleme

### Çıkarım

Karmaşık iş akışlarını verimli bir şekilde organize etmek için otomatik biçimlendirme, yorum yapma özellikleri, kod katlama, parantez eşleştirme ve çok satırlı düzenleme kullanarak temiz, okunabilir kod koruyabilirsiniz.

### Sırada ne var?

VS Code'un sadece kod düzenlemesinin ötesinde daha geniş geliştirme iş akışınızla nasıl entegre olduğunu öğrenin.

---

## 7. Geliştirme İş Akışı Entegrasyonu

VS Code, sadece kod düzenlemenin ötesinde geliştirme iş akışınızla iyi entegre olur.

### 7.1. Versiyon Kontrol Entegrasyonu

!!! note "Codespaces ve Git Entegrasyonu"

    **GitHub Codespaces**'te çalışıyorsanız, bazı Git entegrasyon özellikleri, özellikle Source Control için klavye kısayolları beklendiği gibi çalışmayabilir. Ayrıca ilk kurulum sırasında dizini Git deposu olarak açmayı reddetmiş olabilirsiniz, bu eğitim amaçları için sorun değil.

Projeniz bir git deposu ise (bunun gibi), VS Code şunları gösterir:

- Renkli göstergelerle değiştirilmiş dosyalar
- Durum çubuğunda Git durumu
- Satır içi diff görünümleri
- Commit ve push yetenekleri

Git değişikliklerini görmek ve commit'leri doğrudan editörde hazırlamak için source control düğmesini (![Source control simgesi](img/source_control_icon.png)) kullanarak Source Control panelini açın (`Ctrl+Shift+G` veya VSCode'u yerel olarak çalıştırıyorsanız `Cmd+Shift+G`).

![Source Control Paneli](img/source_control.png)

### 7.2. İş Akışlarını Çalıştırma ve İnceleme

Bir iş akışı çalıştıralım ve ardından sonuçları inceleyelim. Entegre terminalde (hem Windows hem de MacOS için `Ctrl+Shift+` ters tırnak), temel iş akışını çalıştırın:

```bash title="Temel iş akışını çalıştır"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

İş akışı çalışırken terminalde gerçek zamanlı çıktı göreceksiniz. Tamamlandıktan sonra, editörünüzden ayrılmadan sonuçları incelemek için VS Code'u kullanabilirsiniz:

1. **Work dizinlerine gidin**: `.nextflow/work`'e göz atmak için dosya gezgini veya terminali kullanın
2. **Log dosyalarını açın**: Onları doğrudan VS Code'da açmak için terminal çıktısındaki log dosyası yollarına tıklayın
3. **Çıktıları inceleyin**: Dosya gezgininde yayınlanmış sonuç dizinlerine göz atın
4. **Yürütme raporlarını görüntüleyin**: HTML raporlarını doğrudan VS Code'da veya tarayıcınızda açın

Bu, birden fazla uygulama arasında geçiş yapmak yerine her şeyi tek bir yerde tutar.

### Çıkarım

Tüm geliştirme sürecinizi tek bir arayüzden yönetmek için VS Code'u versiyon kontrol ve iş akışı yürütme ile entegre edebilirsiniz.

### Sırada ne var?

Tüm bu IDE özelliklerinin günlük geliştirme iş akışınızda nasıl birlikte çalıştığını görün.

---

## 8. Özet ve Hızlı Notlar

Yukarıda tartışılan her bir IDE özelliği hakkında bazı hızlı notlar:

### 8.1. Yeni Bir Özelliğe Başlama

1. İlgili mevcut modülleri bulmak için **Hızlı dosya açma** (`Ctrl+P` veya `Cmd+P`)
2. Benzer süreçleri yan yana görüntülemek için **Bölünmüş editör**
3. Dosya yapısını anlamak için **Sembol gezinmesi** (`Ctrl+Shift+O` veya `Cmd+Shift+O`)
4. Yeni kodu hızlıca yazmak için **Otomatik tamamlama**

### 8.2. Sorunları Giderme

1. Tüm hataları bir kerede görmek için **Sorunlar paneli** (`Ctrl+Shift+M` veya `Cmd+Shift+M`)
2. Süreç arayüzlerini anlamak için **Tanıma git** (`Ctrl-tıklama` veya `Cmd-tıklama`)
3. Süreçlerin nasıl kullanıldığını görmek için **Tüm referansları bul**
4. Benzer kalıpları veya sorunları bulmak için **Proje genelinde arama**

### 8.3. Yeniden Düzenleme ve İyileştirme

1. Kalıpları bulmak için **Proje genelinde arama** (`Ctrl+Shift+F` veya `Cmd+Shift+F`)
2. Tutarlılığı korumak için **Otomatik biçimlendirme** (`Shift+Alt+F` veya `Shift+Option+F`)
3. Yapıya odaklanmak için **Kod katlama**
4. Değişiklikleri izlemek için **Git entegrasyonu**

---

## Özet

Şimdi Nextflow geliştirme için VS Code'un IDE özelliklerinin hızlı bir turunu tamamladınız. Bu araçlar sizi şu şekillerde önemli ölçüde daha üretken hale getirecektir:

- Gerçek zamanlı sözdizimi kontrolü ile **hataları azaltma**
- Akıllı otomatik tamamlama ile **geliştirmeyi hızlandırma**
- Karmaşık çok dosyalı iş akışlarında **gezinmeyi iyileştirme**
- Tutarlı biçimlendirme ile **kaliteyi koruma**
- Gelişmiş vurgulama ve yapı görselleştirme ile **anlamayı geliştirme**

Her şeyi hatırlamanızı beklemiyoruz, ancak şimdi bu özelliklerin var olduğunu biliyorsunuz ve ihtiyacınız olduğunda onları bulabileceksiniz. Nextflow iş akışları geliştirmeye devam ettikçe, bu IDE özellikleri doğal hale gelecek ve sözdizimi ve yapıyla boğuşmak yerine yüksek kaliteli kod yazmaya odaklanmanızı sağlayacaktır.

### Sırada ne var?

Bu IDE becerilerini diğer eğitim modülleri üzerinde çalışırken uygulayın, örneğin:

- **[nf-test](nf-test.md)**: İş akışlarınız için kapsamlı test paketleri oluşturun
- **[Hello nf-core](../../hello_nf-core/)**: Topluluk standartlarıyla üretim kalitesinde boru hatları oluşturun

Bu IDE özelliklerinin gerçek gücü, daha büyük, daha karmaşık projeler üzerinde çalışırken ortaya çıkar. Onları iş akışınıza kademeli olarak dahil etmeye başlayın—birkaç oturum içinde doğal hale gelecekler ve Nextflow geliştirmeye yaklaşım şeklinizi dönüştüreceklerdir.

Sizi yavaşlatmadan önce hataları yakalamaktan karmaşık kod tabanlarında kolayca gezinmeye kadar, bu araçlar sizi daha özgüvenli ve verimli bir geliştirici yapacaktır.

Mutlu kodlamalar!
