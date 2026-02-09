# Geliştirme Ortamı

Modern Entegre Geliştirme Ortamları (IDE'ler), Nextflow geliştirme deneyiminizi önemli ölçüde dönüştürebilir. Bu yan görev, özellikle VS Code ve Nextflow eklentisini kullanarak daha hızlı kod yazma, hataları erken yakalama ve karmaşık iş akışlarında verimli gezinme konularına odaklanır.

!!! note "Bu geleneksel bir eğitim değildir"

    Diğer eğitim modüllerinin aksine, bu kılavuz adım adım bir eğitimden ziyade hızlı ipuçları, öneriler ve pratik örnekler koleksiyonu olarak düzenlenmiştir. Her bölüm, ilgi alanlarınıza ve mevcut geliştirme ihtiyaçlarınıza göre bağımsız olarak keşfedilebilir. İstediğiniz gibi gezinin ve iş akışı geliştirmeniz için en yararlı olacak özelliklere odaklanın.

## Önce bilmeniz gerekenler

Bu kılavuz, [Hello Nextflow](../hello_nextflow/) eğitim kursunu tamamladığınızı ve aşağıdakiler dahil temel Nextflow kavramlarında rahat olduğunuzu varsayar:

- **Temel iş akışı yapısı**: Süreçleri, iş akışlarını ve bunların nasıl birbirine bağlandığını anlama
- **Kanal işlemleri**: Kanal oluşturma, süreçler arasında veri aktarma ve temel operatörleri kullanma
- **Modüller ve organizasyon**: Yeniden kullanılabilir modüller oluşturma ve include ifadelerini kullanma
- **Yapılandırma temelleri**: Parametreler, süreç yönergeleri ve profiller için `nextflow.config` kullanma

## Burada neler öğreneceksiniz

Bu kılavuz, sizi daha verimli bir Nextflow geliştiricisi yapacak **IDE üretkenlik özelliklerine** odaklanır:

- **Gelişmiş sözdizimi vurgulama**: VS Code'un kod yapınız hakkında size ne gösterdiğini anlama
- **Akıllı otomatik tamamlama**: Daha hızlı kod yazımı için bağlama duyarlı önerileri kullanma
- **Hata algılama ve tanılama**: İş akışınızı çalıştırmadan önce sözdizimi hatalarını yakalama
- **Kod gezinme**: Süreçler, modüller ve tanımlar arasında hızlıca hareket etme
- **Biçimlendirme ve organizasyon**: Tutarlı, okunabilir kod stili sürdürme
- **Yapay zeka destekli geliştirme** (isteğe bağlı): IDE'nizle entegre modern yapay zeka araçlarını kullanma

!!! info "Neden şimdi IDE özellikleri?"

    [Hello Nextflow](../hello_nextflow/) kursu boyunca muhtemelen zaten VS Code kullanıyordunuz, ancak odağı IDE özelliklerinden ziyade Nextflow temellerini öğrenmeye yönelttik. Artık süreçler, iş akışları, kanallar ve modüller gibi temel Nextflow kavramlarında rahat olduğunuza göre, sizi daha verimli bir geliştirici yapacak gelişmiş IDE özelliklerinden yararlanmaya hazırsınız.

    Bunu geliştirme ortamınızı "seviye atlama" olarak düşünün - kullandığınız aynı editör, neye yardımcı olduklarını anladığınızda gerçekten değerli hale gelen çok daha güçlü yeteneklere sahiptir.

---

## 0. Kurulum ve Isınma

IDE özelliklerini keşfetmek için özel olarak bir çalışma alanı kuralım:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Bu dizini VS Code'da açın:

```bash title="Open VS Code in current directory"
code .
```

`ide_features` dizini, çeşitli IDE özelliklerini gösteren örnek iş akışları içerir:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
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

    - `basic_workflow.nf`, çalıştırabileceğiniz ve değiştirebileceğiniz çalışan temel bir iş akışıdır
    - `complex_workflow.nf`, yalnızca gösterim için tasarlanmıştır ve gezinme özelliklerini gösterir - başarıyla çalışmayabilir ancak gerçekçi çok dosyalı iş akışı yapısını gösterir

### Klavye Kısayolları

Bu kılavuzdaki bazı özellikler isteğe bağlı klavye kısayolları kullanacaktır. Bu materyale tarayıcıda GitHub Codespaces üzerinden erişiyor olabilirsiniz ve bu durumda bazen kısayollar beklendiği gibi çalışmayacaktır çünkü sisteminizde başka şeyler için kullanılırlar.

VS Code'u yerel olarak çalıştırıyorsanız, muhtemelen gerçekten iş akışları yazarken olacağınız gibi, kısayollar açıklandığı gibi çalışacaktır.

Mac kullanıyorsanız, bazı (hepsi değil) klavye kısayolları "ctrl" yerine "cmd" kullanacaktır ve bunu metinde `Ctrl/Cmd` şeklinde belirteceğiz.

### 0.1. Nextflow Eklentisini Yükleme

!!! note "Zaten Devcontainer Kullanıyor musunuz?"

    **GitHub Codespaces**'te çalışıyorsanız veya **yerel bir devcontainer** kullanıyorsanız, Nextflow eklentisi muhtemelen zaten yüklenmiş ve sizin için yapılandırılmıştır. Aşağıdaki manuel yükleme adımlarını atlayabilir ve doğrudan eklenti özelliklerini keşfetmeye geçebilirsiniz.

Eklentiyi manuel olarak yüklemek için:

1. VS Code'u açın
2. Soldaki eklentiler simgesine tıklayarak Eklentiler görünümüne gidin: ![eklentiler simgesi](img/extensions_icon.png) (VSCode'u yerel olarak çalıştırıyorsanız kısayol `Ctrl/Cmd+Shift+X`)
3. "Nextflow" araması yapın
4. Resmi Nextflow eklentisini yükleyin

![Nextflow Eklentisini Yükle](img/install_extension.png)

### 0.2. Çalışma Alanı Düzeni

Hello Nextflow boyunca VS Code kullandığınız için, temellere zaten aşinasınız. Bu oturum için çalışma alanınızı verimli bir şekilde nasıl düzenleyeceğiniz aşağıda açıklanmıştır:

- **Editör Alanı**: Dosyaları görüntülemek ve düzenlemek için. Dosyaları yan yana karşılaştırmak için bunu birden fazla bölmeye ayırabilirsiniz.
- **Dosya Gezgini** tıklayın (![dosya gezgini simgesi](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Sisteminizdeki yerel dosyalar ve klasörler. Dosyalar arasında gezinmek için bunu solda açık tutun
- **Entegre Terminal** (`Ctrl+Shift+` ters tırnak hem Windows hem de MacOS için): Altta bilgisayarla etkileşim kurmak için bir terminal. Nextflow veya diğer komutları çalıştırmak için bunu kullanın.
- **Sorunlar Paneli** (`Ctrl+Shift+M`): VS Code, algıladığı hataları ve sorunları burada gösterecektir. Bu, sorunları bir bakışta vurgulamak için kullanışlıdır.

Örnekler üzerinde çalışırken çalışma alanınızı özelleştirmek için panelleri sürükleyebilir veya gizleyebilirsiniz (kenar çubuğunu açıp kapatmak için `Ctrl/Cmd+B`).

### Özet

VS Code'u Nextflow eklentisiyle kurdunuz ve verimli geliştirme için çalışma alanı düzenini anlıyorsunuz.

### Sırada ne var?

Sözdizimi vurgulamanın Nextflow kod yapısını bir bakışta anlamanıza nasıl yardımcı olduğunu öğrenin.

---

## 1. Sözdizimi Vurgulama ve Kod Yapısı

Artık çalışma alanınız kurulduğuna göre, VS Code'un sözdizimi vurgulamasının Nextflow kodunu daha etkili bir şekilde okumanıza ve yazmanıza nasıl yardımcı olduğunu keşfedelim.

### 1.1. Nextflow Sözdizimi Öğeleri

Sözdizimi vurgulamasını çalışırken görmek için `basic_workflow.nf` dosyasını açın:

![Sözdizimi Vitrini](img/syntax_showcase.png)

VS Code'un nasıl vurguladığına dikkat edin:

- **Anahtar kelimeler** (`process`, `workflow`, `input`, `output`, `script`) farklı renklerde
- **String değişmezleri** ve **parametreler** farklı stillerle
- **Yorumlar** soluk bir renkte
- **Değişkenler** ve **fonksiyon çağrıları** uygun vurgulamayla
- **Kod blokları** uygun girinti kılavuzlarıyla

!!! note "Temaya Bağlı Renkler"

    Gördüğünüz belirli renkler, VS Code temanıza (koyu/açık mod), renk ayarlarınıza ve yaptığınız özelleştirmelere bağlı olacaktır. Önemli olan, farklı sözdizimi öğelerinin birbirinden görsel olarak ayırt edilmesidir, bu da seçtiğiniz renk şemasından bağımsız olarak kod yapısını anlamayı kolaylaştırır.

### 1.2. Kod Yapısını Anlama

Sözdizimi vurgulaması, aşağıdakileri hızlıca belirlemenize yardımcı olur:

- **Süreç sınırları**: Farklı süreçler arasında net ayrım
- **Girdi/çıktı blokları**: Veri akışı tanımlarını kolayca belirleme
- **Script blokları**: Yürütülen gerçek komutlar
- **Kanal işlemleri**: Veri dönüştürme adımları
- **Yapılandırma yönergeleri**: Sürece özgü ayarlar

Bu görsel organizasyon, birden fazla süreç ve karmaşık veri akışları içeren karmaşık iş akışlarıyla çalışırken paha biçilmez hale gelir.

### Özet

VS Code'un sözdizimi vurgulamasının Nextflow kod yapısını okumanıza ve daha hızlı geliştirme için farklı dil öğelerini belirlemenize nasıl yardımcı olduğunu anlıyorsunuz.

### Sırada ne var?

Akıllı otomatik tamamlamanın bağlama duyarlı önerilerle kod yazmayı nasıl hızlandırdığını öğrenin.

---

## 2. Akıllı Otomatik Tamamlama

VS Code'un otomatik tamamlama özellikleri, bağlama göre uygun seçenekler önererek daha hızlı ve daha az hatayla kod yazmanıza yardımcı olur.

### 2.1. Bağlama Duyarlı Öneriler

Otomatik tamamlama seçenekleri, kodunuzda nerede olduğunuza bağlı olarak değişir:

#### Kanal İşlemleri

`basic_workflow.nf` dosyasını tekrar açın ve workflow bloğunda `channel.` yazmayı deneyin:

![Kanal otomatik tamamlama](img/autocomplete_channel.png)

Şunlar için öneriler göreceksiniz:

- `fromPath()` - Dosya yollarından kanal oluştur
- `fromFilePairs()` - Eşleştirilmiş dosyalardan kanal oluştur
- `of()` - Değerlerden kanal oluştur
- `fromSRA()` - SRA erişim numaralarından kanal oluştur
- Ve çok daha fazlası...

Bu, tam metod isimlerini hatırlamaya gerek kalmadan kullanılacak doğru kanal fabrikasını hızlıca bulmanıza yardımcı olur.

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

Bu, süreçleri yapılandırırken zaman kazandırır ve farklı yapılandırma kapsamlarında çalışır. Örneğin, Docker'a özgü yapılandırma seçeneklerini görmek için `docker.` yazmayı deneyin.

### Özet

Sözdizimini ezberlemeden mevcut kanal işlemlerini, süreç yönergelerini ve yapılandırma seçeneklerini keşfetmek için VS Code'un akıllı otomatik tamamlamasını kullanabilirsiniz.

### Sırada ne var?

Gerçek zamanlı hata algılamanın, iş akışınızı çalıştırmadan önce sorunları yakalamanıza nasıl yardımcı olduğunu, sadece kodu okuyarak öğrenin.

## 3. Hata Algılama ve Tanılama

VS Code'un gerçek zamanlı hata algılaması, iş akışınızı çalıştırmadan önce sorunları yakalamanıza yardımcı olur.

### 3.1. Sözdizimi Hatası Algılama

Algılamayı çalışırken görmek için kasıtlı bir hata oluşturalım. `basic_workflow.nf` dosyasını açın ve süreç adını `FASTQC`'den `FASTQ`'ya (veya başka bir geçersiz isme) değiştirin. VS Code, workflow bloğundaki hatayı hemen kırmızı dalgalı alt çizgiyle vurgulayacaktır:

![Hata alt çizgisi](img/error_underline.png)

### 3.2. Sorunlar Paneli

Bireysel hata vurgulamanın ötesinde, VS Code çalışma alanınızdaki tüm hataları, uyarıları ve bilgi mesajlarını toplayan merkezi bir Sorunlar paneli sağlar. `Ctrl/Cmd+Shift+M` ile açın ve yalnızca mevcut dosyayla ilgili hataları göstermek için filtre simgesini kullanın:

![Sorunlar panelini filtrele](img/active_file.png)

Sorunlu satıra doğrudan atlamak için herhangi bir soruna tıklayın

![Sorunlar Paneli](img/problems_panel.png)

Süreç adını tekrar `FASTQC` olarak değiştirerek hatayı düzeltin.

### 3.3. Yaygın Hata Kalıpları

Nextflow sözdizimindeki yaygın hatalar şunları içerir:

- **Eksik parantezler**: Eşleşmeyen `{` veya `}`
- **Eksik bloklar**: Süreçlerde gerekli bölümlerin eksik olması
- **Geçersiz sözdizimi**: Hatalı biçimlendirilmiş Nextflow DSL
- **Anahtar kelimelerdeki yazım hataları**: Yanlış yazılmış süreç yönergeleri
- **Kanal uyumsuzlukları**: Tür uyumsuzlukları

Nextflow dil sunucusu bu sorunları Sorunlar panelinde vurgular. Bir pipeline çalıştırırken sözdizimi hatalarından kaçınmak için bunları erken kontrol edebilirsiniz.

### Özet

İş akışınızı çalıştırmadan önce sözdizimi hatalarını ve sorunları yakalamak için VS Code'un hata algılamasını ve Sorunlar panelini kullanabilir, zaman kazanabilir ve hayal kırıklığını önleyebilirsiniz.

### Sırada ne var?

Karmaşık iş akışlarında süreçler, modüller ve tanımlar arasında verimli bir şekilde nasıl gezineceğinizi öğrenin.

---

## 4. Kod Gezinme ve Sembol Yönetimi

Birden fazla dosyaya yayılan karmaşık iş akışlarıyla çalışırken verimli gezinme çok önemlidir. Bunu anlamak için, `basic_workflow.nf` dosyasındaki süreç tanımını size sağladığımız modül için bir içe aktarma ile değiştirin:

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

Bu özellik, iş akışları yazarken özellikle değerlidir, çünkü modül dosyasını doğrudan açmadan modül arayüzünü anlamanıza olanak tanır.

**Ctrl/Cmd-tıklama** kullanarak herhangi bir süreç, modül veya değişken tanımına hızlıca gidebilirsiniz. Betiğin üst kısmındaki modül dosyasına bağlantının üzerine fareyle gelin ve önerildiği gibi bağlantıyı takip edin:

![Bağlantıyı takip et](img/follow_link.png)

Aynı şey süreç adları için de çalışır. `basic_workflow.nf` dosyasına geri dönün ve bunu workflow bloğundaki `FASTQC` süreç adında deneyin. Bu sizi doğrudan süreç adına bağlar (bu örnekte modül dosyasıyla aynıdır, ancak çok daha büyük bir dosyanın ortasında olabilir).

Bulunduğunuz yere geri dönmek için **Alt+←** (veya Mac'te **Ctrl+-**) kullanın. Bu, yerinizi kaybetmeden kodu keşfetmenin güçlü bir yoludur.

Şimdi `complex_workflow.nf` (daha önce bahsedilen yalnızca gösterim amaçlı dosya) kullanarak daha karmaşık bir iş akışında gezinmeyi keşfedelim. Bu iş akışı, ayrı modül dosyalarında tanımlanmış birden fazla sürecin yanı sıra bazı satır içi süreçler içerir. Karmaşık çok dosyalı yapılar manuel olarak gezinmek için zorlayıcı olabilirken, tanımlara atlama yeteneği keşfi çok daha yönetilebilir hale getirir.

1. `complex_workflow.nf` dosyasını açın
2. Modül tanımlarına gidin
3. Geri gitmek için **Alt+←** (veya **Ctrl+-**) kullanın
4. Workflow bloğundaki `FASTQC` süreç adına gidin. Bu sizi doğrudan süreç adına bağlar (bu örnekte modül dosyasıyla aynıdır, ancak çok daha büyük bir dosyanın ortasında olabilir).
5. Tekrar geri gidin
6. Workflow bloğundaki `TRIM_GALORE` sürecine gidin. Bu satır içi tanımlanmıştır, bu nedenle sizi ayrı bir dosyaya götürmeyecektir, ancak yine de süreç tanımını gösterecek ve bulunduğunuz yere geri gidebileceksiniz.

### 4.2. Sembol Gezinme

`complex_workflow.nf` hala açıkken, VSCode'un üst kısmındaki arama çubuğuna `@` yazarak dosyadaki tüm sembollere genel bir bakış elde edebilirsiniz (klavye kısayolu `Ctrl/Cmd+Shift+O`'dur, ancak Codespaces'te çalışmayabilir). Bu, mevcut dosyadaki tüm sembolleri listeleyen sembol gezinme panelini açar:

![Sembol gezinme](img/symbols.png)

Bu şunları gösterir:

- Tüm süreç tanımları
- İş akışı tanımları (bu dosyada iki iş akışı tanımlanmıştır)
- Fonksiyon tanımları

Sonuçları filtrelemek için yazmaya başlayın.

### 4.3. Tüm Referansları Bul

Bir sürecin veya değişkenin kod tabanınız boyunca nerede kullanıldığını anlamak çok yararlı olabilir. Örneğin, `FASTQC` sürecine yapılan tüm referansları bulmak istiyorsanız, tanımına giderek başlayın. Bunu `modules/fastqc.nf` dosyasını doğrudan açarak veya yukarıda yaptığımız gibi VS Code'un hızlı gezinme özelliğini `Ctrl/Cmd-tıklama` ile kullanarak yapabilirsiniz. Süreç tanımına ulaştığınızda, `FASTQC` süreç adına sağ tıklayın ve kullanıldığı tüm örnekleri görmek için bağlam menüsünden "Tüm Referansları Bul"ı seçin.

![Referansları bul](img/references.png)

Bu özellik, `FASTQC`'nin çalışma alanınızda referans verildiği tüm örnekleri, iki farklı iş akışındaki kullanımı dahil olmak üzere görüntüler. Bu içgörü, `FASTQC` sürecinde yapılacak değişikliklerin potansiyel etkisini değerlendirmek için çok önemlidir.

### 4.4. Anahat Paneli

Gezgin kenar çubuğunda (![Gezgin simgesi](img/files_icon.png) tıklayın) bulunan Anahat paneli, mevcut dosyanızdaki tüm sembollere uygun bir genel bakış sağlar. Bu özellik, fonksiyonları, değişkenleri ve diğer anahtar öğeleri hiyerarşik bir görünümde görüntüleyerek kodunuzun yapısında hızlıca gezinmenize ve yönetmenize olanak tanır.

![Anahat paneli](img/outline.png)

Dosya tarayıcısını kullanmadan kodunuzun farklı bölümlerine hızlıca gitmek için Anahat panelini kullanın.

### 4.5. DAG görselleştirme

VS Code'un Nextflow eklentisi, iş akışınızı Yönlendirilmiş Döngüsüz Grafik (DAG) olarak görselleştirebilir. Bu, süreçler arasındaki veri akışını ve bağımlılıkları anlamanıza yardımcı olur. `complex_workflow.nf` dosyasını açın ve `workflow {` üzerindeki "Preview DAG" düğmesine tıklayın (bu dosyadaki ikinci `workflow` bloğu):

![DAG önizleme](img/dag_preview.png)

Bu sadece 'giriş' iş akışıdır, ancak daha yukarıdaki `RNASEQ_PIPELINE {` iş akışının üzerindeki "Preview DAG" düğmesine tıklayarak iç iş akışları için de DAG'yi önizleyebilirsiniz:

![DAG önizleme iç iş akışı](img/dag_preview_inner.png)

Bu iş akışı için, editördeki ilgili süreç tanımlarına gitmek için DAG'deki düğümleri kullanabilirsiniz. Bir düğüme tıklayın ve sizi editördeki ilgili süreç tanımına götürecektir. Özellikle bir iş akışı büyük bir boyuta ulaştığında, bu gerçekten kod etrafında gezinmenize ve süreçlerin nasıl bağlandığını anlamanıza yardımcı olabilir.

### Özet

Kod yapısını ve bağımlılıkları anlamak için tanıma git, sembol arama, referansları bul ve DAG görselleştirme kullanarak karmaşık iş akışlarında verimli bir şekilde gezinebilirsiniz.

### Sırada ne var?

Daha büyük Nextflow projelerinde birbirine bağlı birden fazla dosyayla etkili bir şekilde nasıl çalışacağınızı öğrenin.

## 5. Birden Fazla Dosyada Çalışma

Gerçek Nextflow geliştirme, birbirine bağlı birden fazla dosyayla çalışmayı içerir. VS Code'un karmaşık projeleri verimli bir şekilde yönetmenize nasıl yardımcı olduğunu keşfedelim.

### 5.1. Hızlı Dosya Gezinme

`complex_workflow.nf` açıkken, birkaç modül içe aktardığını fark edeceksiniz. Aralarında hızlı gezinme pratiği yapalım.

**Ctrl+P** (veya **Cmd+P**) tuşlarına basın ve "fast" yazmaya başlayın:

VS Code size eşleşen dosyaları gösterecektir. Oraya anında atlamak için `modules/fastqc.nf` dosyasını seçin. Bu, hangi dosyayı aradığınızı kabaca bildiğinizde dosya gezgininde tıklamaktan çok daha hızlıdır.

Bunu diğer kalıplarla deneyin:

- STAR hizalama modül dosyasını (`star.nf`) bulmak için "star" yazın
- Yardımcı fonksiyonlar dosyasını (`utils.nf`) bulmak için "utils" yazın
- Yapılandırma dosyalarına (`nextflow.config`) atlamak için "config" yazın

### 5.2. Çok Dosyalı Geliştirme için Bölünmüş Editör

Modüllerle çalışırken, genellikle hem ana iş akışını hem de modül tanımlarını aynı anda görmeniz gerekir. Bunu kuralım:

1. `complex_workflow.nf` dosyasını açın
2. `modules/fastqc.nf` dosyasını yeni bir sekmede açın
3. `modules/fastqc.nf` sekmesine sağ tıklayın ve "Split Right"ı seçin
4. Şimdi her iki dosyayı yan yana görebilirsiniz

![Bölünmüş editör](img/split_editor.png)

Bu şu durumlarda paha biçilmezdir:

- İş akışı çağrıları yazarken modül arayüzlerini kontrol etme ve önizleme yeterli değilse
- Farklı modüller arasında benzer süreçleri karşılaştırma
- İş akışı ve modüller arasındaki veri akışında hata ayıklama

### 5.3. Proje Genelinde Arama

Bazen belirli kalıpların tüm projenizde nerede kullanıldığını bulmanız gerekir. Arama panelini açmak için `Ctrl/Cmd+Shift+F` tuşlarına basın.

Çalışma alanı genelinde `publishDir` aramayı deneyin:

![Proje araması](img/project_search.png)

Bu size yayınlama dizinlerini kullanan her dosyayı gösterir ve şunlara yardımcı olur:

- Çıktı organizasyon kalıplarını anlama
- Belirli yönergelerin örneklerini bulma
- Modüller arasında tutarlılığı sağlama

### Özet

İş akışları ve modüller arasında verimli çalışmak için hızlı dosya gezinme, bölünmüş editörler ve proje genelinde arama kullanarak karmaşık çok dosyalı projeleri yönetebilirsiniz.

### Sırada ne var?

Kod biçimlendirme ve bakım özelliklerinin iş akışlarınızı düzenli ve okunabilir tutmasını öğrenin.

---

## 6. Kod Biçimlendirme ve Bakım

Uygun kod biçimlendirme, yalnızca estetik için değil, aynı zamanda okunabilirliği, anlaşılabilirliği ve karmaşık iş akışlarını güncelleme kolaylığını artırmak için de önemlidir.

### 6.1. Otomatik Biçimlendirme Çalışırken

`basic_workflow.nf` dosyasını açın ve biçimlendirmeyi kasıtlı olarak bozun:

- Bazı girintileri kaldırın: Tüm belgeyi vurgulayın ve mümkün olduğunca çok girintiyi kaldırmak için `shift+tab` tuşlarına çok kez basın.
- Rastgele yerlere ekstra boşluklar ekleyin: `channel.fromPath` ifadesinde, `(` sonrasına 30 boşluk ekleyin.
- Bazı satırları garip bir şekilde bölün: `.view {` operatörü ile `Processing sample:` dizesi arasına yeni bir satır ekleyin ancak kapanış parantezi `}` öncesine karşılık gelen bir yeni satır eklemeyin.

Şimdi otomatik biçimlendirme için `Shift+Alt+F` (veya MacOS'ta `Shift+Option+F`) tuşlarına basın:

VS Code hemen:

- Süreç yapısını açıkça göstermek için girintiyi düzeltir
- Benzer öğeleri tutarlı bir şekilde hizalar
- Gereksiz boşlukları kaldırır
- Okunabilir satır sonlarını korur

Otomatik biçimlendirmenin her kod stili sorununu çözemeyebileceğini unutmayın. Nextflow dil sunucusu kodunuzu düzenli tutmayı amaçlar, ancak belirli alanlarda kişisel tercihlerinize de saygı gösterir. Örneğin, bir sürecin `script` bloğunun içindeki girintiyi kaldırırsanız, biçimlendirici onu olduğu gibi bırakacaktır, çünkü kasıtlı olarak bu stili tercih ediyor olabilirsiniz.

Şu anda Nextflow için katı bir stil zorlaması yoktur, bu nedenle dil sunucusu bir miktar esneklik sunar. Ancak, netliği korumak için metod ve fonksiyon tanımları etrafında tutarlı bir şekilde biçimlendirme kuralları uygulayacaktır.

### 6.2. Kod Organizasyon Özellikleri

#### Hızlı Yorum Yapma

İş akışınızda bir kod bloğu seçin ve yorum yapmak için **Ctrl+/** (veya **Cmd+/**) tuşlarına basın:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Bu şunlar için mükemmeldir:

- Geliştirme sırasında iş akışlarının bölümlerini geçici olarak devre dışı bırakma
- Karmaşık kanal işlemlerine açıklayıcı yorumlar ekleme
- İş akışı bölümlerini belgeleme

Kodu yorum satırından çıkarmak için tekrar **Ctrl+/** (veya **Cmd+/**) kullanın.

#### Genel Bakış için Kod Katlama

`complex_workflow.nf` dosyasında, süreç tanımlarının yanındaki küçük oklara dikkat edin. Süreçleri katlamak (daraltmak) için onlara tıklayın:

![Kod katlama](img/code_folding.png)

Bu, uygulama ayrıntılarında kaybolmadan iş akışı yapınızın üst düzey bir görünümünü verir.

#### Parantez Eşleştirme

İmlecinizi herhangi bir `{` veya `}` parantezinin yanına yerleştirin ve VS Code eşleşen parantezi vurgular. Eşleşen parantezler arasında atlamak için **Ctrl+Shift+\\** (veya **Cmd+Shift+\\**) kullanın.

Bu şunlar için çok önemlidir:

- Süreç sınırlarını anlama
- Eksik veya fazla parantezleri bulma
- İç içe iş akışı yapılarında gezinme

#### Çok Satırlı Seçim ve Düzenleme

Birden fazla satırı aynı anda düzenlemek için, VS Code güçlü çoklu imleç yetenekleri sunar:

- **Çok satırlı seçim**: **Ctrl+Alt** (veya MacOS için **Cmd+Option**) tuşunu basılı tutun ve birden fazla satır seçmek için ok tuşlarını kullanın
- **Çok satırlı girinti**: Birden fazla satır seçin ve tüm blokları girintilendirmek için **Tab** veya girintiyi azaltmak için **Shift+Tab** kullanın

Bu özellikle şunlar için kullanışlıdır:

- Tüm süreç bloklarını tutarlı bir şekilde girintilendirme
- Birden fazla satıra aynı anda yorum ekleme
- Birden fazla süreçte benzer parametre tanımlarını düzenleme

### Özet

Karmaşık iş akışlarını verimli bir şekilde düzenlemek için otomatik biçimlendirme, yorum yapma özellikleri, kod katlama, parantez eşleştirme ve çok satırlı düzenleme kullanarak temiz, okunabilir kod sürdürebilirsiniz.

### Sırada ne var?

VS Code'un sadece kod düzenlemesinin ötesinde daha geniş geliştirme iş akışınızla nasıl entegre olduğunu öğrenin.

---

## 7. Geliştirme İş Akışı Entegrasyonu

VS Code, sadece kod düzenlemenin ötesinde geliştirme iş akışınızla iyi entegre olur.

### 7.1. Sürüm Kontrolü Entegrasyonu

!!! note "Codespaces ve Git Entegrasyonu"

    **GitHub Codespaces**'te çalışıyorsanız, bazı Git entegrasyon özellikleri, özellikle Kaynak Kontrolü için klavye kısayolları beklendiği gibi çalışmayabilir. Ayrıca ilk kurulum sırasında dizini bir Git deposu olarak açmayı reddetmiş olabilirsiniz, bu eğitim amaçları için sorun değildir.

Projeniz bir git deposuysa (bu öyle), VS Code şunları gösterir:

- Renkli göstergelerle değiştirilmiş dosyalar
- Durum çubuğunda Git durumu
- Satır içi fark görünümleri
- Commit ve push yetenekleri

Git değişikliklerini görmek ve commit'leri doğrudan editörde aşamalandırmak için kaynak kontrol düğmesini (![Kaynak kontrol simgesi](img/source_control_icon.png)) kullanarak Kaynak Kontrolü panelini açın (`Ctrl+Shift+G` veya VSCode'u yerel olarak çalıştırıyorsanız `Cmd+Shift+G`).

![Kaynak Kontrolü Paneli](img/source_control.png)

### 7.2. İş Akışlarını Çalıştırma ve İnceleme

Bir iş akışı çalıştıralım ve ardından sonuçları inceleyelim. Entegre terminalde (`Ctrl+Shift+` ters tırnak hem Windows hem de MacOS'ta), temel iş akışını çalıştırın:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

İş akışı çalışırken, terminalde gerçek zamanlı çıktı göreceksiniz. Tamamlandıktan sonra, editörünüzden ayrılmadan sonuçları incelemek için VS Code'u kullanabilirsiniz:

1. **Work dizinlerine gidin**: `.nextflow/work` dizinine göz atmak için dosya gezginini veya terminali kullanın
2. **Log dosyalarını açın**: Doğrudan VS Code'da açmak için terminal çıktısındaki log dosyası yollarına tıklayın
3. **Çıktıları inceleyin**: Dosya gezgininde yayınlanmış sonuç dizinlerine göz atın
4. **Yürütme raporlarını görüntüleyin**: HTML raporlarını doğrudan VS Code'da veya tarayıcınızda açın

Bu, birden fazla uygulama arasında geçiş yapmak yerine her şeyi tek bir yerde tutar.

### Özet

Tüm geliştirme sürecinizi tek bir arayüzden yönetmek için VS Code'u sürüm kontrolü ve iş akışı yürütmesiyle entegre edebilirsiniz.

### Sırada ne var?

Tüm bu IDE özelliklerinin günlük geliştirme iş akışınızda nasıl birlikte çalıştığını görün.

---

## 8. Özet ve hızlı notlar

Yukarıda tartışılan IDE özelliklerinin her biri hakkında bazı hızlı notlar:

### 8.1. Yeni Bir Özellik Başlatma

1. İlgili mevcut modülleri bulmak için **Hızlı dosya açma** (`Ctrl+P` veya `Cmd+P`)
2. Benzer süreçleri yan yana görüntülemek için **Bölünmüş editör**
3. Dosya yapısını anlamak için **Sembol gezinme** (`Ctrl+Shift+O` veya `Cmd+Shift+O`)
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

Artık Nextflow geliştirme için VS Code'un IDE özelliklerinin hızlı bir turunu tamamladınız. Bu araçlar sizi şunlar sayesinde önemli ölçüde daha üretken yapacaktır:

- Gerçek zamanlı sözdizimi kontrolü ile **hataları azaltma**
- Akıllı otomatik tamamlama ile **geliştirmeyi hızlandırma**
- Karmaşık çok dosyalı iş akışlarında **gezinmeyi iyileştirme**
- Tutarlı biçimlendirme ile **kaliteyi sürdürme**
- Gelişmiş vurgulama ve yapı görselleştirme ile **anlamayı geliştirme**

Her şeyi hatırlamanızı beklemiyoruz, ancak artık bu özelliklerin var olduğunu bildiğinize göre, ihtiyacınız olduğunda bunları bulabileceksiniz. Nextflow iş akışları geliştirmeye devam ettikçe, bu IDE özellikleri doğal hale gelecek ve sözdizimi ve yapıyla boğuşmak yerine yüksek kaliteli kod yazmaya odaklanmanıza olanak tanıyacaktır.

### Sırada ne var?

Bu IDE becerilerini diğer eğitim modülleri üzerinde çalışırken uygulayın, örneğin:

- **[nf-test](nf-test.md)**: İş akışlarınız için kapsamlı test paketleri oluşturun
- **[Hello nf-core](../../hello_nf-core/)**: Topluluk standartlarıyla üretim kalitesinde pipeline'lar oluşturun

Bu IDE özelliklerinin gerçek gücü, daha büyük, daha karmaşık projeler üzerinde çalışırken ortaya çıkar. Bunları iş akışınıza kademeli olarak dahil etmeye başlayın—birkaç oturum içinde doğal hale gelecekler ve Nextflow geliştirmeye yaklaşımınızı dönüştüreceklerdir.

Sizi yavaşlatmadan önce hataları yakalamaktan karmaşık kod tabanlarında kolaylıkla gezinmeye kadar, bu araçlar sizi daha kendinden emin ve verimli bir geliştirici yapacaktır.

Mutlu kodlamalar!
