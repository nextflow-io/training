# Geliştirme Ortamı

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Modern Entegre Geliştirme Ortamları (IDE'ler), Nextflow geliştirme deneyiminizi köklü biçimde dönüştürebilir. Bu yan görev, özellikle VS Code ve Nextflow eklentisinden yararlanarak daha hızlı kod yazmaya, hataları erken yakalamaya ve karmaşık iş akışlarında verimli biçimde gezinmeye odaklanmaktadır.

!!! note "Bu geleneksel bir eğitim değil"

    Diğer eğitim modüllerinin aksine, bu rehber adım adım ilerleyen bir eğitim yerine hızlı ipuçları, pratik öneriler ve örneklerden oluşan bir koleksiyon olarak düzenlenmiştir. Her bölüm, ilgi alanlarınıza ve mevcut geliştirme ihtiyaçlarınıza göre bağımsız olarak incelenebilir. İş akışı geliştirmenize en hızlı katkı sağlayacak özelliklere odaklanarak istediğiniz bölüme atlayabilirsiniz.

## Önce bilmeniz gerekenler

Bu rehber, [Hello Nextflow](../hello_nextflow/) eğitim kursunu tamamladığınızı ve aşağıdakiler dahil temel Nextflow kavramlarına hâkim olduğunuzu varsaymaktadır:

- **Temel iş akışı yapısı**: Süreçleri, iş akışlarını ve bunların birbirine nasıl bağlandığını anlamak
- **Kanal işlemleri**: Kanal oluşturma, süreçler arasında veri aktarma ve temel operatörleri kullanma
- **Modüller ve organizasyon**: Yeniden kullanılabilir modüller oluşturma ve include ifadelerini kullanma
- **Yapılandırma temelleri**: Parametreler, süreç yönergeleri ve profiller için `nextflow.config` kullanımı

## Burada neler öğreneceksiniz

Bu rehber, sizi daha verimli bir Nextflow geliştiricisi yapacak **IDE üretkenlik özelliklerine** odaklanmaktadır:

- **Gelişmiş sözdizimi vurgulama**: VS Code'un kod yapınız hakkında size ne gösterdiğini anlamak
- **Akıllı otomatik tamamlama**: Daha hızlı kod yazımı için bağlama duyarlı önerilerden yararlanmak
- **Hata tespiti ve tanılama**: İş akışınızı çalıştırmadan önce sözdizimi hatalarını yakalamak
- **Kod gezintisi**: Süreçler, modüller ve tanımlar arasında hızlıca geçiş yapmak
- **Biçimlendirme ve organizasyon**: Tutarlı ve okunabilir kod stilini korumak
- **Yapay zeka destekli geliştirme** (isteğe bağlı): IDE'nize entegre modern yapay zeka araçlarını kullanmak

!!! info "IDE özellikleri neden şimdi?"

    [Hello Nextflow](../hello_nextflow/) kursu boyunca büyük olasılıkla VS Code kullanıyordunuz; ancak odağı IDE özellikleri yerine Nextflow temellerini öğrenmeye yönelttik. Artık süreçler, iş akışları, kanallar ve modüller gibi temel Nextflow kavramlarına hâkim olduğunuza göre, sizi daha verimli bir geliştirici yapacak gelişmiş IDE özelliklerinden yararlanmaya hazırsınız.

    Bunu geliştirme ortamınızı "bir üst seviyeye taşımak" olarak düşünebilirsiniz; kullandığınız editör, ne yaptıklarını anladığınızda gerçek değerini ortaya koyan çok daha güçlü yeteneklere sahiptir.

---

## 0. Kurulum ve Isınma

IDE özelliklerini keşfetmek için özel bir çalışma alanı oluşturalım:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Bu dizini VS Code'da açın:

```bash title="Open VS Code in current directory"
code .
```

`ide_features` dizini, çeşitli IDE özelliklerini gösteren örnek iş akışları içermektedir:

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

    - `basic_workflow.nf`, çalıştırıp değiştirebileceğiniz, işlevsel bir temel iş akışıdır
    - `complex_workflow.nf`, yalnızca gezinti özelliklerini göstermek amacıyla tasarlanmıştır; başarıyla çalışmayabilir, ancak gerçekçi bir çok dosyalı iş akışı yapısını göstermektedir

### Klavye Kısayolları

Bu rehberdeki bazı özellikler isteğe bağlı klavye kısayolları kullanmaktadır. Bu materyale tarayıcı üzerinden GitHub Codespaces aracılığıyla erişiyorsanız, kısayollar sisteminizde başka amaçlarla kullanıldığından zaman zaman beklendiği gibi çalışmayabilir.

VS Code'u yerel olarak çalıştırıyorsanız (iş akışı yazarken büyük olasılıkla böyle yapacaksınız), kısayollar açıklandığı gibi çalışacaktır.

Mac kullanıyorsanız, bazı (tümü değil) klavye kısayolları "ctrl" yerine "cmd" tuşunu kullanır; bunu metinde `Ctrl/Cmd` şeklinde belirteceğiz.

### 0.1. Nextflow Eklentisini Yükleme

!!! note "Devcontainer Kullanıyor musunuz?"

    **GitHub Codespaces**'te çalışıyorsanız veya **yerel bir devcontainer** kullanıyorsanız, Nextflow eklentisi büyük olasılıkla zaten yüklenmiş ve yapılandırılmıştır. Aşağıdaki manuel yükleme adımlarını atlayarak doğrudan eklenti özelliklerini keşfetmeye geçebilirsiniz.

Eklentiyi manuel olarak yüklemek için:

1. VS Code'u açın
2. Soldaki eklentiler simgesine tıklayarak Eklentiler görünümüne gidin: ![extensions icon](img/extensions_icon.png) (VS Code'u yerel olarak çalıştırıyorsanız kısayol: `Ctrl/Cmd+Shift+X`)
3. "Nextflow" araması yapın
4. Resmi Nextflow eklentisini yükleyin

![Nextflow Eklentisini Yükle](img/install_extension.png)

### 0.2. Çalışma Alanı Düzeni

Hello Nextflow boyunca VS Code kullandığınız için temel bilgilere zaten hâkimsiniz. Bu oturum için çalışma alanınızı verimli biçimde düzenlemenin yolları şunlardır:

- **Editör Alanı**: Dosyaları görüntülemek ve düzenlemek için. Dosyaları yan yana karşılaştırmak amacıyla birden fazla bölmeye bölebilirsiniz.
- **Dosya Gezgini** (![file explorer icon](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Sisteminizdeki yerel dosya ve klasörler. Dosyalar arasında gezinmek için bunu solda açık tutun.
- **Entegre Terminal** (hem Windows hem MacOS için `Ctrl+Shift+` backtick): Altta bilgisayarla etkileşim kurmak için bir terminal. Nextflow veya diğer komutları çalıştırmak için kullanın.
- **Sorunlar Paneli** (`Ctrl+Shift+M`): VS Code'un tespit ettiği hata ve sorunları burada gösterir. Sorunlara tek bakışta göz atmak için kullanışlıdır.

Örnekler üzerinde çalışırken düzeninizi özelleştirmek için panelleri sürükleyebilir veya gizleyebilirsiniz (`Ctrl/Cmd+B` ile kenar çubuğunu açıp kapatabilirsiniz).

### Özetle

VS Code'u Nextflow eklentisiyle kurdunuz ve verimli geliştirme için çalışma alanı düzenini öğrendiniz.

### Sırada ne var?

Sözdizimi vurgulamanın Nextflow kod yapısını tek bakışta anlamanıza nasıl yardımcı olduğunu öğrenin.

---

## 1. Sözdizimi Vurgulama ve Kod Yapısı

Çalışma alanınız hazır olduğuna göre, VS Code'un sözdizimi vurgulamasının Nextflow kodunu daha etkili biçimde okumanıza ve yazmanıza nasıl yardımcı olduğunu keşfedelim.

### 1.1. Nextflow Sözdizimi Öğeleri

Sözdizimi vurgulamasını görmek için `basic_workflow.nf` dosyasını açın:

![Sözdizimi Vitrini](img/syntax_showcase.png)

VS Code'un şunları nasıl vurguladığına dikkat edin:

- **Anahtar kelimeler** (`process`, `workflow`, `input`, `output`, `script`) farklı renklerle
- **String literaller** ve **parametreler** farklı stillerle
- **Yorumlar** soluk bir renkle
- **Değişkenler** ve **fonksiyon çağrıları** uygun vurguyla
- **Kod blokları** düzgün girinti kılavuzlarıyla

!!! note "Temaya Bağlı Renkler"

    Göreceğiniz renkler VS Code temanıza (koyu/açık mod), renk ayarlarınıza ve yaptığınız özelleştirmelere bağlıdır. Önemli olan, farklı sözdizimi öğelerinin birbirinden görsel olarak ayırt edilmesidir; bu sayede seçtiğiniz renk şemasından bağımsız olarak kod yapısını anlamak kolaylaşır.

### 1.2. Kod Yapısını Anlamak

Sözdizimi vurgulama şunları hızlıca tanımlamanıza yardımcı olur:

- **Süreç sınırları**: Farklı süreçler arasındaki net ayrım
- **Girdi/çıktı blokları**: Veri akışı tanımlarını kolayca fark etme
- **Script blokları**: Çalıştırılan gerçek komutlar
- **Kanal işlemleri**: Veri dönüşüm adımları
- **Yapılandırma yönergeleri**: Sürece özgü ayarlar

Bu görsel organizasyon, birden fazla süreç ve karmaşık veri akışları içeren iş akışlarıyla çalışırken son derece değerli hale gelir.

### Özetle

VS Code'un sözdizimi vurgulamasının Nextflow kod yapısını okumanıza ve farklı dil öğelerini tanımlamanıza nasıl yardımcı olduğunu öğrendiniz; bu sayede daha hızlı geliştirme yapabilirsiniz.

### Sırada ne var?

Akıllı otomatik tamamlamanın bağlama duyarlı önerilerle kod yazmayı nasıl hızlandırdığını öğrenin.

---

## 2. Akıllı Otomatik Tamamlama

VS Code'un otomatik tamamlama özellikleri, bağlama göre uygun seçenekler önererek daha hızlı ve daha az hatayla kod yazmanıza yardımcı olur.

### 2.1. Bağlama Duyarlı Öneriler

Otomatik tamamlama seçenekleri, kodunuzda bulunduğunuz konuma göre değişir:

#### Kanal İşlemleri

`basic_workflow.nf` dosyasını tekrar açın ve workflow bloğuna `channel.` yazmayı deneyin:

![Kanal otomatik tamamlama](img/autocomplete_channel.png)

Şu öneriler görünecektir:

- `fromPath()` - Dosya yollarından kanal oluşturur
- `fromFilePairs()` - Eşleştirilmiş dosyalardan kanal oluşturur
- `of()` - Değerlerden kanal oluşturur
- `fromSRA()` - SRA erişim numaralarından kanal oluşturur
- Ve daha fazlası...

Bu özellik, tam metot adlarını ezberlemek zorunda kalmadan doğru kanal fabrikasını hızlıca bulmanıza yardımcı olur.

Kanallara uygulanabilecek operatörleri de keşfedebilirsiniz. Örneğin, mevcut işlemleri görmek için `FASTQC.out.html.` yazın:

![Kanal işlemleri otomatik tamamlama](img/autocomplete_operators.png)

#### Süreç Yönergeleri

Bir sürecin script bloğunun içinde `task.` yazarak mevcut çalışma zamanı özelliklerini görün:

![Görev özellikleri otomatik tamamlama](img/autocomplete_task.png)

#### Yapılandırma

`nextflow.config` dosyasını açın ve herhangi bir yere `process.` yazarak mevcut süreç yönergelerini görün:

![Yapılandırma otomatik tamamlama](img/autocomplete_config.png)

Şu öneriler görünecektir:

- `executor`
- `memory`
- `cpus`

Bu özellik, süreçleri yapılandırırken zaman kazandırır ve farklı yapılandırma kapsamlarında çalışır. Örneğin, Docker'a özgü yapılandırma seçeneklerini görmek için `docker.` yazmayı deneyin.

### Özetle

VS Code'un akıllı otomatik tamamlamasını kullanarak sözdizimini ezberlemek zorunda kalmadan mevcut kanal işlemlerini, süreç yönergelerini ve yapılandırma seçeneklerini keşfedebilirsiniz.

### Sırada ne var?

Gerçek zamanlı hata tespitinin, yalnızca kodu okuyarak iş akışınızı çalıştırmadan önce sorunları yakalamanıza nasıl yardımcı olduğunu öğrenin.

## 3. Hata Tespiti ve Tanılama

VS Code'un gerçek zamanlı hata tespiti, iş akışınızı çalıştırmadan önce sorunları yakalamanıza yardımcı olur.

### 3.1. Sözdizimi Hatası Tespiti

Tespiti görmek için kasıtlı bir hata oluşturalım. `basic_workflow.nf` dosyasını açın ve süreç adını `FASTQC`'den `FASTQ`'ya (veya başka geçersiz bir ada) değiştirin. VS Code, workflow bloğundaki hatayı kırmızı dalgalı alt çizgiyle hemen vurgulayacaktır:

![Hata alt çizgisi](img/error_underline.png)

### 3.2. Sorunlar Paneli

Bireysel hata vurgulamanın ötesinde, VS Code çalışma alanınızdaki tüm hataları, uyarıları ve bilgi mesajlarını toplayan merkezi bir Sorunlar paneli sunar. `Ctrl/Cmd+Shift+M` ile açın ve yalnızca mevcut dosyayla ilgili hataları göstermek için filtre simgesini kullanın:

![Sorunlar panelini filtrele](img/active_file.png)

Doğrudan sorunlu satıra atlamak için herhangi bir soruna tıklayın.

![Sorunlar Paneli](img/problems_panel.png)

Süreç adını `FASTQC` olarak geri değiştirerek hatayı düzeltin.

### 3.3. Yaygın Hata Kalıpları

Nextflow sözdizimindeki yaygın hatalar şunlardır:

- **Eksik parantezler**: Eşleşmeyen `{` veya `}`
- **Tamamlanmamış bloklar**: Süreçlerde eksik zorunlu bölümler
- **Geçersiz sözdizimi**: Hatalı biçimlendirilmiş Nextflow DSL
- **Anahtar kelime yazım hataları**: Yanlış yazılmış süreç yönergeleri
- **Kanal uyumsuzlukları**: Tür uyumsuzlukları

Nextflow dil sunucusu bu sorunları Sorunlar panelinde vurgular. Bir pipeline çalıştırırken sözdizimi hatalarından kaçınmak için bunları erken kontrol edebilirsiniz.

### Özetle

VS Code'un hata tespiti ve Sorunlar panelini kullanarak iş akışınızı çalıştırmadan önce sözdizimi hatalarını ve sorunları yakalayabilir, zaman kazanabilir ve hayal kırıklığını önleyebilirsiniz.

### Sırada ne var?

Karmaşık iş akışlarında süreçler, modüller ve tanımlar arasında verimli biçimde nasıl gezineceğinizi öğrenin.

---

## 4. Kod Gezintisi ve Sembol Yönetimi

Birden fazla dosyaya yayılan karmaşık iş akışlarıyla çalışırken verimli gezinti son derece önemlidir. Bunu anlamak için `basic_workflow.nf` dosyasındaki süreç tanımını, size sağladığımız modül için bir import ifadesiyle değiştirin:

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

`FASTQC` gibi bir süreç adının üzerine fare ile geldiğinizde, modül arayüzünü (girdiler ve çıktılar) gösteren bir açılır pencere görürsünüz:

![Tanıma git](img/syntax.png)

Bu özellik, modül dosyasını doğrudan açmak zorunda kalmadan modül arayüzünü anlamanıza olanak tanıdığından iş akışı yazarken özellikle değerlidir.

**Ctrl/Cmd-tıklama** kullanarak herhangi bir süreç, modül veya değişken tanımına hızlıca gidebilirsiniz. Scriptin üst kısmındaki modül dosyasına olan bağlantının üzerine gelin ve önerilen bağlantıyı takip edin:

![Bağlantıyı takip et](img/follow_link.png)

Aynı şey süreç adları için de geçerlidir. `basic_workflow.nf` dosyasına geri dönün ve workflow bloğundaki `FASTQC` süreç adı üzerinde bunu deneyin. Bu sizi doğrudan süreç adına götürür (bu örnekte modül dosyasıyla aynıdır, ancak çok daha büyük bir dosyanın ortasında da olabilir).

Bulunduğunuz yere geri dönmek için **Alt+←** (veya Mac'te **Ctrl+-**) tuşlarını kullanın. Bu, yerinizi kaybetmeden kodu keşfetmenin güçlü bir yoludur.

Şimdi `complex_workflow.nf` (daha önce bahsedilen yalnızca gösterim amaçlı dosya) kullanarak daha karmaşık bir iş akışında gezinmeyi keşfedelim. Bu iş akışı, ayrı modül dosyalarında tanımlanmış birden fazla süreç ile bazı satır içi süreçler içermektedir. Karmaşık çok dosyalı yapılar manuel olarak gezinmesi zor olsa da tanımlara atlama yeteneği keşfetmeyi çok daha yönetilebilir kılar.

1. `complex_workflow.nf` dosyasını açın
2. Modül tanımlarına gidin
3. Geri gitmek için **Alt+←** (veya **Ctrl+-**) tuşlarını kullanın
4. Workflow bloğundaki `FASTQC` süreç adına gidin. Bu sizi doğrudan süreç adına götürür (bu örnekte modül dosyasıyla aynıdır, ancak çok daha büyük bir dosyanın ortasında da olabilir).
5. Tekrar geri gidin
6. Workflow bloğundaki `TRIM_GALORE` sürecine gidin. Bu süreç satır içi tanımlandığından sizi ayrı bir dosyaya götürmez; ancak yine de süreç tanımını gösterir ve bulunduğunuz yere geri dönebilirsiniz.

### 4.2. Sembol Gezintisi

`complex_workflow.nf` hâlâ açıkken, VSCode'un üst kısmındaki arama çubuğuna `@` yazarak dosyadaki tüm sembollere genel bir bakış elde edebilirsiniz (klavye kısayolu `Ctrl/Cmd+Shift+O`'dur, ancak Codespaces'te çalışmayabilir). Bu, mevcut dosyadaki tüm sembolleri listeleyen sembol gezinti panelini açar:

![Sembol gezintisi](img/symbols.png)

Bu panel şunları gösterir:

- Tüm süreç tanımları
- İş akışı tanımları (bu dosyada iki iş akışı tanımlanmıştır)
- Fonksiyon tanımları

Sonuçları filtrelemek için yazmaya başlayın.

### 4.3. Tüm Referansları Bul

Bir sürecin veya değişkenin kod tabanınızın neresinde kullanıldığını anlamak çok yardımcı olabilir. Örneğin, `FASTQC` sürecine yapılan tüm referansları bulmak istiyorsanız önce tanımına gidin. Bunu `modules/fastqc.nf` dosyasını doğrudan açarak veya yukarıda yaptığımız gibi `Ctrl/Cmd-tıklama` ile VS Code'un hızlı gezinti özelliğini kullanarak yapabilirsiniz. Süreç tanımına ulaştıktan sonra `FASTQC` süreç adına sağ tıklayın ve bağlam menüsünden "Find All References" seçeneğini seçerek kullanıldığı tüm yerleri görün.

![Referansları bul](img/references.png)

Bu özellik, `FASTQC`'nin çalışma alanınızda referans alındığı tüm yerleri (iki farklı iş akışındaki kullanımları dahil) gösterir. Bu bilgi, `FASTQC` sürecinde yapılacak değişikliklerin olası etkisini değerlendirmek için son derece önemlidir.

### 4.4. Taslak Paneli

Gezgin kenar çubuğunda (![Explorer icon](img/files_icon.png) simgesine tıklayın) bulunan Taslak paneli, mevcut dosyanızdaki tüm sembollere uygun bir genel bakış sunar. Bu özellik, fonksiyonları, değişkenleri ve diğer temel öğeleri hiyerarşik bir görünümde göstererek kodunuzun yapısında hızlıca gezinmenizi ve yönetmenizi sağlar.

![Taslak paneli](img/outline.png)

Dosya tarayıcısını kullanmadan kodunuzun farklı bölümlerine hızlıca gitmek için Taslak panelini kullanın.

### 4.5. DAG Görselleştirme

VS Code'un Nextflow eklentisi, iş akışınızı Yönlendirilmiş Döngüsüz Graf (DAG) olarak görselleştirebilir. Bu, süreçler arasındaki veri akışını ve bağımlılıkları anlamanıza yardımcı olur. `complex_workflow.nf` dosyasını açın ve `workflow {` ifadesinin üzerindeki "Preview DAG" düğmesine tıklayın (bu dosyadaki ikinci `workflow` bloğu):

![DAG önizleme](img/dag_preview.png)

Bu yalnızca 'giriş' iş akışıdır; daha yukarıdaki `RNASEQ_PIPELINE {` iş akışı için "Preview DAG" düğmesine tıklayarak iç iş akışlarının DAG'ını da önizleyebilirsiniz:

![DAG önizleme iç iş akışı](img/dag_preview_inner.png)

Bu iş akışı için DAG'daki düğümleri kullanarak koddaki ilgili süreç tanımlarına gidebilirsiniz. Bir düğüme tıkladığınızda editörde ilgili süreç tanımına götürülürsünüz. Özellikle bir iş akışı büyük boyutlara ulaştığında bu özellik, kod içinde gezinmenize ve süreçlerin nasıl bağlandığını anlamanıza gerçekten yardımcı olabilir.

### Özetle

Tanıma git, sembol arama, referansları bul ve DAG görselleştirme özelliklerini kullanarak karmaşık iş akışlarında verimli biçimde gezinebilir, kod yapısını ve bağımlılıkları anlayabilirsiniz.

### Sırada ne var?

Daha büyük Nextflow projelerinde birbirine bağlı birden fazla dosyayla nasıl etkili çalışacağınızı öğrenin.

## 5. Birden Fazla Dosyayla Çalışmak

Gerçek Nextflow geliştirmesi, birbirine bağlı birden fazla dosyayla çalışmayı gerektirir. VS Code'un karmaşık projeleri verimli biçimde yönetmenize nasıl yardımcı olduğunu keşfedelim.

### 5.1. Hızlı Dosya Gezintisi

`complex_workflow.nf` açıkken, birkaç modül içe aktardığını fark edeceksiniz. Aralarında hızlı gezinti yapmayı pratik edelim.

**Ctrl+P** (veya **Cmd+P**) tuşlarına basın ve "fast" yazmaya başlayın:

VS Code eşleşen dosyaları gösterecektir. Anında oraya atlamak için `modules/fastqc.nf` dosyasını seçin. Aradığınız dosyayı kabaca bildiğinizde bu yöntem, dosya gezgininde tıklamaktan çok daha hızlıdır.

Bunu diğer kalıplarla da deneyin:

- STAR hizalama modül dosyasını (`star.nf`) bulmak için "star" yazın
- Yardımcı fonksiyonlar dosyasını (`utils.nf`) bulmak için "utils" yazın
- Yapılandırma dosyalarına (`nextflow.config`) atlamak için "config" yazın

### 5.2. Çok Dosyalı Geliştirme için Bölünmüş Editör

Modüllerle çalışırken genellikle hem ana iş akışını hem de modül tanımlarını aynı anda görmeniz gerekir. Bunu ayarlayalım:

1. `complex_workflow.nf` dosyasını açın
2. `modules/fastqc.nf` dosyasını yeni bir sekmede açın
3. `modules/fastqc.nf` sekmesine sağ tıklayın ve "Split Right" seçeneğini seçin
4. Artık her iki dosyayı yan yana görebilirsiniz

![Bölünmüş editör](img/split_editor.png)

Bu özellik şu durumlarda son derece değerlidir:

- İş akışı çağrıları yazarken modül arayüzlerini kontrol etmek ve önizleme yeterli olmadığında
- Farklı modüllerdeki benzer süreçleri karşılaştırmak
- İş akışı ile modüller arasındaki veri akışında hata ayıklamak

### 5.3. Proje Genelinde Arama

Bazen belirli kalıpların tüm projenizde nerede kullanıldığını bulmanız gerekir. Arama panelini açmak için `Ctrl/Cmd+Shift+F` tuşlarına basın.

Çalışma alanında `publishDir` araması yapmayı deneyin:

![Proje araması](img/project_search.png)

Bu, yayımlama dizinlerini kullanan her dosyayı gösterir ve şunlara yardımcı olur:

- Çıktı organizasyon kalıplarını anlamak
- Belirli yönergelerin örneklerini bulmak
- Modüller arasında tutarlılığı sağlamak

### Özetle

Hızlı dosya gezintisi, bölünmüş editörler ve proje genelinde arama kullanarak karmaşık çok dosyalı projeleri yönetebilir ve iş akışları ile modüller arasında verimli biçimde çalışabilirsiniz.

### Sırada ne var?

Kod biçimlendirme ve bakım özelliklerinin iş akışlarınızı nasıl düzenli ve okunabilir tuttuğunu öğrenin.

---

## 6. Kod Biçimlendirme ve Bakım

Doğru kod biçimlendirmesi yalnızca estetik açıdan değil, okunabilirliği, anlaşılırlığı ve karmaşık iş akışlarını güncelleme kolaylığını artırmak açısından da son derece önemlidir.

### 6.1. Otomatik Biçimlendirme Uygulamalı

`basic_workflow.nf` dosyasını açın ve biçimlendirmeyi kasıtlı olarak bozun:

- Bazı girintileri kaldırın: Tüm belgeyi seçin ve mümkün olduğunca çok girinti kaldırmak için `shift+tab` tuşlarına birçok kez basın.
- Rastgele yerlere fazladan boşluk ekleyin: `channel.fromPath` ifadesinde `(` işaretinden sonra 30 boşluk ekleyin.
- Bazı satırları garip biçimde bölün: `.view {` operatörü ile `Processing sample:` dizesi arasına yeni bir satır ekleyin, ancak kapanış parantezi `}` öncesine karşılık gelen bir satır sonu eklemeyin.

Şimdi otomatik biçimlendirme için `Shift+Alt+F` (veya MacOS'ta `Shift+Option+F`) tuşlarına basın:

VS Code hemen şunları yapar:

- Süreç yapısını açıkça göstermek için girintileri düzeltir
- Benzer öğeleri tutarlı biçimde hizalar
- Gereksiz boşlukları kaldırır
- Okunabilir satır sonlarını korur

Otomatik biçimlendirmenin her kod stili sorununu çözmeyebileceğini unutmayın. Nextflow dil sunucusu kodunuzu düzenli tutmayı amaçlar; ancak belirli alanlarda kişisel tercihlerinize de saygı gösterir. Örneğin, bir sürecin `script` bloğu içindeki girintileri kaldırırsanız, biçimlendirici bunu olduğu gibi bırakır; çünkü bu stili kasıtlı olarak tercih ediyor olabilirsiniz.

Şu anda Nextflow için katı bir stil zorunluluğu bulunmamaktadır; bu nedenle dil sunucusu belirli esneklikler sunar. Ancak netliği korumak amacıyla metot ve fonksiyon tanımları etrafında biçimlendirme kurallarını tutarlı biçimde uygular.

### 6.2. Kod Organizasyon Özellikleri

#### Hızlı Yorum Ekleme

İş akışınızda bir kod bloğu seçin ve yorum satırına dönüştürmek için **Ctrl+/** (veya **Cmd+/**) tuşlarına basın:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Bu özellik şu durumlarda mükemmeldir:

- Geliştirme sırasında iş akışlarının bölümlerini geçici olarak devre dışı bırakmak
- Karmaşık kanal işlemlerine açıklayıcı yorumlar eklemek
- İş akışı bölümlerini belgelemek

Kodu yorum satırından çıkarmak için **Ctrl+/** (veya **Cmd+/**) tuşlarına tekrar basın.

#### Genel Bakış için Kod Katlama

`complex_workflow.nf` dosyasında süreç tanımlarının yanındaki küçük oklara dikkat edin. Süreçleri katlamak (daraltmak) için bunlara tıklayın:

![Kod katlama](img/code_folding.png)

Bu, uygulama ayrıntılarında kaybolmadan iş akışı yapınıza üst düzey bir bakış sağlar.

#### Parantez Eşleştirme

İmlecinizi herhangi bir `{` veya `}` parantezinin yanına getirdiğinizde VS Code eşleşen parantezi vurgular. Eşleşen parantezler arasında atlamak için **Ctrl+Shift+\\** (veya **Cmd+Shift+\\**) tuşlarını kullanın.

Bu özellik şunlar için son derece önemlidir:

- Süreç sınırlarını anlamak
- Eksik veya fazla parantezleri bulmak
- İç içe geçmiş iş akışı yapılarında gezinmek

#### Çok Satırlı Seçim ve Düzenleme

Birden fazla satırı aynı anda düzenlemek için VS Code güçlü çok imleçli özellikler sunar:

- **Çok satırlı seçim**: **Ctrl+Alt** (veya MacOS için **Cmd+Option**) tuşlarını basılı tutun ve birden fazla satır seçmek için ok tuşlarını kullanın
- **Çok satırlı girinti**: Birden fazla satır seçin ve tüm blokları girintili yapmak için **Tab**, girintisiz yapmak için **Shift+Tab** tuşlarını kullanın

Bu özellik özellikle şu durumlarda kullanışlıdır:

- Tüm süreç bloklarını tutarlı biçimde girintili yapmak
- Birden fazla satıra aynı anda yorum eklemek
- Birden fazla süreçteki benzer parametre tanımlarını düzenlemek

### Özetle

Otomatik biçimlendirme, yorum ekleme özellikleri, kod katlama, parantez eşleştirme ve çok satırlı düzenleme kullanarak karmaşık iş akışlarını verimli biçimde organize etmek için temiz ve okunabilir kod koruyabilirsiniz.

### Sırada ne var?

VS Code'un yalnızca kod düzenlemenin ötesinde geliştirme iş akışınızla nasıl entegre olduğunu öğrenin.

---

## 7. Geliştirme İş Akışı Entegrasyonu

VS Code, yalnızca kod düzenlemenin ötesinde geliştirme iş akışınızla iyi biçimde entegre olur.

### 7.1. Sürüm Kontrolü Entegrasyonu

!!! note "Codespaces ve Git Entegrasyonu"

    **GitHub Codespaces**'te çalışıyorsanız, bazı Git entegrasyon özellikleri beklendiği gibi çalışmayabilir; özellikle Kaynak Kontrolü için klavye kısayolları. Ayrıca ilk kurulum sırasında dizini bir Git deposu olarak açmayı reddettiyseniz, bu eğitim amaçları için sorun değildir.

Projeniz bir git deposuysa (bu proje gibi), VS Code şunları gösterir:

- Renkli göstergelerle değiştirilmiş dosyalar
- Durum çubuğunda Git durumu
- Satır içi fark görünümleri
- Commit ve push özellikleri

Git değişikliklerini görmek ve commit'leri doğrudan editörde hazırlamak için kaynak kontrolü düğmesini (![Source control icon](img/source_control_icon.png)) kullanarak Kaynak Kontrolü panelini açın (`Ctrl+Shift+G` veya VS Code'u yerel olarak çalıştırıyorsanız `Cmd+Shift+G`).

![Kaynak Kontrolü Paneli](img/source_control.png)

### 7.2. İş Akışlarını Çalıştırma ve İnceleme

Bir iş akışı çalıştıralım ve ardından sonuçları inceleyelim. Entegre terminalde (hem Windows hem MacOS için `Ctrl+Shift+` backtick), temel iş akışını çalıştırın:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

İş akışı çalışırken terminalde gerçek zamanlı çıktı göreceksiniz. Tamamlandıktan sonra editörden çıkmadan VS Code'u kullanarak sonuçları inceleyebilirsiniz:

1. **Work dizinlerine gidin**: `.nextflow/work` dizininde gezinmek için dosya gezginini veya terminali kullanın
2. **Log dosyalarını açın**: Terminal çıktısındaki log dosyası yollarına tıklayarak bunları doğrudan VS Code'da açın
3. **Çıktıları inceleyin**: Dosya gezgininde yayımlanan sonuç dizinlerine göz atın
4. **Çalıştırma raporlarını görüntüleyin**: HTML raporlarını doğrudan VS Code'da veya tarayıcınızda açın

Bu yaklaşım, birden fazla uygulama arasında geçiş yapmak yerine her şeyi tek bir yerde tutar.

### Özetle

VS Code'u sürüm kontrolü ve iş akışı çalıştırmayla entegre ederek tüm geliştirme sürecinizi tek bir arayüzden yönetebilirsiniz.

### Sırada ne var?

Tüm bu IDE özelliklerinin günlük geliştirme iş akışınızda birlikte nasıl çalıştığını görün.

---

## 8. Özet ve Hızlı Notlar

Yukarıda ele alınan IDE özelliklerinin her biri için bazı hızlı notlar:

### 8.1. Yeni Bir Özellik Başlatmak

1. İlgili mevcut modülleri bulmak için **Hızlı dosya açma** (`Ctrl+P` veya `Cmd+P`)
2. Benzer süreçleri yan yana görüntülemek için **Bölünmüş editör**
3. Dosya yapısını anlamak için **Sembol gezintisi** (`Ctrl+Shift+O` veya `Cmd+Shift+O`)
4. Yeni kod hızlıca yazmak için **Otomatik tamamlama**

### 8.2. Sorunları Ayıklamak

1. Tüm hataları bir arada görmek için **Sorunlar paneli** (`Ctrl+Shift+M` veya `Cmd+Shift+M`)
2. Süreç arayüzlerini anlamak için **Tanıma git** (`Ctrl-tıklama` veya `Cmd-tıklama`)
3. Süreçlerin nasıl kullanıldığını görmek için **Tüm referansları bul**
4. Benzer kalıpları veya sorunları bulmak için **Proje genelinde arama**

### 8.3. Yeniden Düzenleme ve İyileştirme

1. Kalıpları bulmak için **Proje genelinde arama** (`Ctrl+Shift+F` veya `Cmd+Shift+F`)
2. Tutarlılığı korumak için **Otomatik biçimlendirme** (`Shift+Alt+F` veya `Shift+Option+F`)
3. Yapıya odaklanmak için **Kod katlama**
4. Değişiklikleri takip etmek için **Git entegrasyonu**

---

## Özet

VS Code'un Nextflow geliştirme için IDE özelliklerine hızlı bir tur attınız. Bu araçlar şunları yaparak sizi önemli ölçüde daha üretken kılacaktır:

- **Hataları azaltmak**: Gerçek zamanlı sözdizimi denetimi sayesinde
- **Geliştirmeyi hızlandırmak**: Akıllı otomatik tamamlama ile
- **Gezinmeyi iyileştirmek**: Karmaşık çok dosyalı iş akışlarında
- **Kaliteyi korumak**: Tutarlı biçimlendirme aracılığıyla
- **Anlamayı geliştirmek**: Gelişmiş vurgulama ve yapı görselleştirmesi ile

Her şeyi hatırlamanızı beklemiyoruz; ancak artık bu özelliklerin var olduğunu bildiğinize göre, ihtiyaç duyduğunuzda onları bulabileceksiniz. Nextflow iş akışları geliştirmeye devam ettikçe bu IDE özellikleri ikinci doğanız haline gelecek ve sözdizimi ile yapıyla boğuşmak yerine yüksek kaliteli kod yazmaya odaklanmanızı sağlayacaktır.

### Sırada ne var?

Bu IDE becerilerini diğer eğitim modülleri üzerinde çalışırken uygulayın; örneğin:

- **[nf-test](nf-test.md)**: İş akışlarınız için kapsamlı test paketleri oluşturun
- **[Hello nf-core](../../hello_nf-core/)**: Topluluk standartlarıyla üretim kalitesinde pipeline'lar oluşturun

Bu IDE özelliklerinin gerçek gücü, daha büyük ve daha karmaşık projeler üzerinde çalıştıkça ortaya çıkar. Bunları iş akışınıza kademeli olarak dahil etmeye başlayın; birkaç oturum içinde ikinci doğanız haline gelecek ve Nextflow geliştirmeye yaklaşımınızı dönüştüreceklerdir.

Hataları sizi yavaşlatmadan yakalamaktan karmaşık kod tabanlarında kolaylıkla gezinmeye kadar bu araçlar, sizi daha özgüvenli ve verimli bir geliştirici yapacaktır.

İyi kodlamalar!
