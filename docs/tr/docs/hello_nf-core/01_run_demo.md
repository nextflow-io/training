# Bölüm 1: Bir demo pipeline çalıştırın

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirme önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core eğitim kursunun bu ilk bölümünde, bir nf-core pipeline'ını nasıl bulacağınızı ve deneyeceğinizi, kodun nasıl organize edildiğini anlamanızı ve [Hello Nextflow](../hello_nextflow/index.md)'da gösterilen sade Nextflow kodundan nasıl farklılaştığını tanımanızı göstereceğiz.

Kod yapısını ve araç işlemlerini göstermek için nf-core projesi tarafından pipeline envanterinin bir parçası olarak sürdürülen nf-core/demo adlı bir pipeline kullanacağız.

[Başlarken](./00_orientation.md) sayfasında belirtildiği gibi çalışma dizininizin `hello-nf-core/` olarak ayarlandığından emin olun.

---

## 1. nf-core/demo pipeline'ını bulun ve alın

[nf-co.re](https://nf-co.re) adresindeki proje web sitesinde nf-core/demo pipeline'ını bularak başlayalım. Bu site, genel dokümantasyon ve yardım makaleleri, her bir pipeline için dokümantasyon, blog yazıları, etkinlik duyuruları vb. gibi tüm bilgileri merkezileştirir.

### 1.1. Pipeline'ı web sitesinde bulun

Web tarayıcınızda [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) adresine gidin ve arama çubuğuna `demo` yazın.

![arama sonuçları](./img/search-results.png)

Pipeline dokümantasyon sayfasına erişmek için pipeline adı olan `demo`'ya tıklayın.

Yayınlanan her pipeline, aşağıdaki dokümantasyon bölümlerini içeren özel bir sayfaya sahiptir:

- **Introduction:** Pipeline'a giriş ve genel bakış
- **Usage:** Pipeline'ın nasıl çalıştırılacağına dair açıklamalar
- **Parameters:** Açıklamalarla birlikte gruplandırılmış pipeline parametreleri
- **Output:** Beklenen çıktı dosyalarının açıklamaları ve örnekleri
- **Results:** Tam test veri setinden üretilen örnek çıktı dosyaları
- **Releases & Statistics:** Pipeline sürüm geçmişi ve istatistikleri

Yeni bir pipeline'ı benimsemeyi düşündüğünüzde, çalıştırmayı denemeden önce ne yaptığını ve nasıl yapılandırılması gerektiğini anlamak için pipeline dokümantasyonunu dikkatlice okumalısınız.

Şimdi bir göz atın ve şunları bulup bulamayacağınızı görün:

- Pipeline hangi araçları çalıştıracak (Sekmeyi kontrol edin: `Introduction`)
- Pipeline hangi girdileri ve parametreleri kabul ediyor veya gerektiriyor (Sekmeyi kontrol edin: `Parameters`)
- Pipeline tarafından üretilen çıktılar neler (Sekmeyi kontrol edin: `Output`)

#### 1.1.1. Pipeline'a genel bakış

`Introduction` sekmesi, görsel bir temsil (metro haritası olarak adlandırılır) ve pipeline'ın bir parçası olarak çalıştırılan araçların bir listesi dahil olmak üzere pipeline'a genel bir bakış sağlar.

![pipeline metro haritası](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Örnek komut satırı

Dokümantasyon ayrıca bir örnek girdi dosyası (aşağıda daha ayrıntılı olarak tartışılacaktır) ve bir örnek komut satırı sağlar.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Örnek komutun bir iş akışı dosyası belirtmediğini, sadece pipeline deposuna referans olan `nf-core/demo`'yu belirttiğini fark edeceksiniz.

Bu şekilde çağrıldığında, Nextflow kodun belirli bir şekilde organize edildiğini varsayacaktır.
Bu yapıyı inceleyebilmemiz için kodu alalım.

### 1.2. Pipeline kodunu alın

Pipeline'ın amacımıza uygun göründüğünü belirledikten sonra, hadi deneyelim.
Neyse ki Nextflow, doğru biçimlendirilmiş depolardan pipeline'ları manuel olarak bir şey indirmek zorunda kalmadan almayı kolaylaştırır.

Terminale dönelim ve aşağıdakini çalıştıralım:

```bash
nextflow pull nf-core/demo
```

??? success "Komut çıktısı"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow, pipeline kodunun bir `pull` işlemini yapar, yani tüm depoyu yerel sürücünüze indirir.

Açık olmak gerekirse, bunu sadece nf-core pipeline'ları için değil, GitHub'da uygun şekilde kurulmuş herhangi bir Nextflow pipeline'ı için yapabilirsiniz.
Ancak nf-core, Nextflow pipeline'larının en büyük açık kaynak koleksiyonudur.

Nextflow'dan bu şekilde aldığınız pipeline'ların bir listesini almanızı sağlayabilirsiniz:

```bash
nextflow list
```

??? success "Komut çıktısı"

    ```console
    nf-core/demo
    ```

Dosyaların mevcut çalışma dizininizde olmadığını fark edeceksiniz.
Varsayılan olarak, Nextflow bunları `$NXF_HOME/assets` dizinine kaydeder.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Dizin içeriği"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note

    Eğitim ortamımızı kullanmıyorsanız, tam yol sisteminizde farklı olabilir.

Nextflow, indirilen kaynak kodunu kasıtlı olarak 'yolun dışında' tutar, çünkü bu pipeline'ların doğrudan etkileşimde bulunacağınız koddan ziyade kütüphaneler gibi kullanılması gerektiği ilkesine dayanır.

Ancak, bu eğitimin amaçları doğrultusunda, içeride ne olduğunu görebilmek istiyoruz.
Bu yüzden bunu kolaylaştırmak için, mevcut çalışma dizinimizden o konuma sembolik bir bağlantı oluşturalım.

```bash
ln -s $NXF_HOME/assets pipelines
```

Bu, az önce indirdiğimiz kodu keşfetmeyi kolaylaştıran bir kısayol oluşturur.

```bash
tree -L 2 pipelines
```

```console title="Dizin içeriği"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Artık gerektiğinde kaynak koduna daha kolay göz atabiliriz.

Ama önce, ilk nf-core pipeline'ımızı çalıştırmayı deneyelim!

### Özet

Artık nf-core web sitesi aracılığıyla bir pipeline'ı nasıl bulacağınızı ve kaynak kodun yerel bir kopyasını nasıl alacağınızı biliyorsunuz.

### Sırada ne var?

Bir nf-core pipeline'ını minimum çabayla nasıl deneyeceğinizi öğrenin.

---

## 2. Pipeline'ı test profiliyle deneyin

Kolaylık sağlayan bir şekilde, her nf-core pipeline'ı bir test profiliyle birlikte gelir.
Bu, [nf-core/test-datasets](https://github.com/nf-core/test-datasets) deposunda barındırılan küçük bir test veri seti kullanarak pipeline'ın çalışması için minimum yapılandırma ayarları setidir.
Küçük ölçekte bir pipeline'ı hızlıca denemek için harika bir yoldur.

!!! note

    Nextflow'un yapılandırma profili sistemi, farklı konteyner motorları veya yürütme ortamları arasında kolayca geçiş yapmanızı sağlar.
    Daha fazla ayrıntı için [Hello Nextflow Bölüm 6: Yapılandırma](../hello_nextflow/06_hello_config.md)'ya bakın.

### 2.1. Test profilini inceleyin

Bir pipeline'ın test profilinin çalıştırmadan önce ne belirttiğini kontrol etmek iyi bir uygulamadır.
`nf-core/demo` için `test` profili `conf/test.config` yapılandırma dosyasında bulunur ve aşağıda gösterilmiştir.

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Üstteki yorum bloğunun, bu test profiliyle pipeline'ın nasıl çalıştırılacağını gösteren bir kullanım örneği içerdiğini hemen fark edeceksiniz.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Sağlamamız gereken tek şeyler örnek komutta açılı parantezler arasında gösterilenlerdir: `<docker/singularity>` ve `<OUTDIR>`.

Hatırlatma olarak, `<docker/singularity>` konteyner sisteminin seçimini ifade eder. Tüm nf-core pipeline'ları, tekrarlanabilirliği sağlamak ve yazılım kurulum sorunlarını ortadan kaldırmak için konteynerlerle (Docker, Singularity, vb.) kullanılabilir olacak şekilde tasarlanmıştır.
Bu yüzden pipeline'ı test etmek için Docker mu yoksa Singularity mi kullanmak istediğimizi belirtmemiz gerekecek.

`--outdir <OUTDIR>` kısmı, Nextflow'un pipeline'ın çıktılarını yazacağı dizini ifade eder.
Bunun için bir isim vermemiz gerekiyor, bunu kendimiz oluşturabiliriz.
Eğer zaten mevcut değilse, Nextflow çalışma zamanında bizim için oluşturacaktır.

Yorum bloğundan sonraki bölüme geçersek, test profili bize test için önceden yapılandırılmış olanları gösterir: en önemlisi, `input` parametresi zaten bir test veri setine işaret edecek şekilde ayarlanmıştır, bu yüzden kendi verilerimizi sağlamamıza gerek yoktur.
Önceden yapılandırılmış girdiye giden bağlantıyı takip ederseniz, bunun birkaç deneysel örnek için örnek tanımlayıcıları ve dosya yolları içeren bir csv dosyası olduğunu göreceksiniz.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Buna samplesheet (örnek tablosu) denir ve nf-core pipeline'larına en yaygın girdi biçimidir.

!!! note

    Veri formatlarına ve türlerine aşina değilseniz endişelenmeyin, takip edenler için önemli değil.

Bu, pipeline'ı denemek için ihtiyacımız olan her şeye sahip olduğumuzu doğrular.

### 2.2. Pipeline'ı çalıştırın

Konteyner sistemi için Docker'ı ve çıktı dizini olarak `demo-results`'ı kullanmaya karar verelim ve test komutunu çalıştırmaya hazırız:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Çıktınız bununla eşleşiyorsa, tebrikler! İlk nf-core pipeline'ınızı çalıştırdınız.

Temel bir Nextflow pipeline'ı çalıştırdığınızda olduğundan çok daha fazla konsol çıktısı olduğunu fark edeceksiniz.
Pipeline'ın sürümünün, girdilerinin ve çıktılarının bir özetini ve birkaç yapılandırma öğesini içeren bir başlık var.

!!! note

    Çıktınız farklı zaman damgaları, yürütme adları ve dosya yolları gösterecektir, ancak genel yapı ve süreç yürütmesi benzer olmalıdır.

Yürütme çıktısına geçersek, hangi süreçlerin çalıştırıldığını bize söyleyen satırlara bir göz atalım:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Bu bize üç sürecin çalıştırıldığını söyler, nf-core web sitesindeki pipeline dokümantasyon sayfasında gösterilen üç araca karşılık gelir: FASTQC, SEQTK_TRIM ve MULTIQC.

Burada gösterildiği gibi `NFCORE_DEMO:DEMO:MULTIQC` gibi tam süreç adları, giriş niteliğindeki Hello Nextflow materyalinde görmüş olabileceğinizden daha uzundur.
Bunlar üst iş akışlarının adlarını içerir ve pipeline kodunun modülerliğini yansıtır.
Bunun hakkında biraz sonra daha fazla ayrıntıya gireceğiz.

### 2.3. Pipeline'ın çıktılarını inceleyin

Son olarak, pipeline tarafından üretilen `demo-results` dizinine bir göz atalım.

```bash
tree -L 2 demo-results
```

??? abstract "Dizin içeriği"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

Bu çok fazla görünebilir.
`nf-core/demo` pipeline'ının çıktıları hakkında daha fazla bilgi edinmek için [dokümantasyon sayfasına](https://nf-co.re/demo/1.0.2/docs/output/) göz atın.

Bu aşamada gözlemlenmesi gereken önemli şey, sonuçların modüle göre organize edilmiş olması ve ek olarak pipeline yürütmesi hakkında çeşitli zaman damgalı raporlar içeren `pipeline_info` adlı bir dizinin bulunmasıdır.

Örneğin, `execution_timeline_*` dosyası size hangi süreçlerin çalıştırıldığını, hangi sırayla ve ne kadar sürdüklerini gösterir:

![yürütme zaman çizelgesi raporu](./img/execution_timeline.png)

!!! note

    Burada görevler paralel olarak çalıştırılmadı çünkü Github Codespaces'te minimalist bir makinede çalışıyoruz.
    Bunların paralel çalıştığını görmek için, codespace'inizin CPU tahsisini ve test yapılandırmasındaki kaynak sınırlarını artırmayı deneyin.

Bu raporlar tüm nf-core pipeline'ları için otomatik olarak oluşturulur.

### Özet

Bir nf-core pipeline'ını yerleşik test profiliyle nasıl çalıştıracağınızı ve çıktılarını nerede bulacağınızı biliyorsunuz.

### Sırada ne var?

Pipeline kodunun nasıl organize edildiğini öğrenin.

---

## 3. Pipeline kod yapısını inceleyin

Artık pipeline'ı kullanıcı olarak başarıyla çalıştırdığımıza göre, nf-core pipeline'larının dahili olarak nasıl yapılandırıldığına bakmak için bakış açımızı değiştirelim.

nf-core projesi, pipeline'ların nasıl yapılandırılacağı ve kodun nasıl organize edileceği, yapılandırılacağı ve belgeleneceği konusunda güçlü yönergeler uygular.
Bunların hepsinin nasıl organize edildiğini anlamak, kursun 2. Bölümünde ele alacağımız kendi nf-core uyumlu pipeline'larınızı geliştirmeye doğru ilk adımdır.

Daha önce oluşturduğumuz `pipelines` sembolik bağlantısını kullanarak `nf-core/demo` deposunda pipeline kodunun nasıl organize edildiğine bir göz atalım.

`nf-core/demo` dizinini bulmak ve açmak için `tree` kullanabilir veya dosya gezginini kullanabilirsiniz.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Dizin içeriği"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

Orada çok şey oluyor, bu yüzden bunu adım adım ele alacağız.

İlk olarak, en üst seviyede özet bilgiler içeren bir README dosyası ve ayrıca lisanslama, katkı yönergeleri, alıntı ve davranış kuralları gibi proje bilgilerini özetleyen yardımcı dosyalar bulabileceğinizi not edelim.
Ayrıntılı pipeline dokümantasyonu `docs` dizininde bulunur.
Tüm bu içerik, nf-core web sitesindeki web sayfalarını programatik olarak oluşturmak için kullanılır, böylece her zaman kodla günceldir.

Şimdi, geri kalanı için keşfimizi üç aşamada bölüyoruz:

1. Pipeline kod bileşenleri (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline yapılandırması
3. Girdiler ve doğrulama

Pipeline kod bileşenleriyle başlayalım.
Bireysel dosyaların içindeki koda dalmak yerine, dosya hiyerarşisine ve yapısal organizasyona odaklanacağız.

### 3.1. Pipeline kod bileşenleri

Standart nf-core pipeline kod organizasyonu, [Hello Nextflow](../hello_nextflow/index.md) kursunun 4. Bölümü olan [Hello Modules](../hello_nextflow/04_hello_modules.md)'de tanıtıldığı gibi kod yeniden kullanımını maksimize etmek için tasarlanmış modüler bir yapıyı takip eder, ancak gerçek nf-core tarzında bu biraz ek karmaşıklıkla uygulanır.
Özellikle, nf-core pipeline'ları, yani bir üst iş akışı tarafından içe aktarılan iş akışı betiklerini bol miktarda kullanır.

Bu biraz soyut gelebilir, bu yüzden bunun `nf-core/demo` pipeline'ında pratikte nasıl kullanıldığına bir göz atalım.

!!! note

    Bu modüler bileşenlerin _nasıl_ bağlandığına dair gerçek kodu gözden geçirmeyeceğiz, çünkü kafa karıştırıcı olabilecek subworkflow'ların kullanımıyla ilişkili bazı ek karmaşıklıklar var ve bunu anlamak eğitimin bu aşamasında gerekli değil.
    Şimdilik, genel organizasyona ve mantığa odaklanacağız.

#### 3.1.1. Genel bakış

`nf-core/demo` pipeline'ı için ilgili kod bileşenleri arasındaki ilişkiler şöyle görünür:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

`main.nf` adlı bir _giriş noktası_ betiği vardır, bu iki tür iç içe iş akışı için bir sarmalayıcı görevi görür: `workflows/` altında bulunan ve `demo.nf` adlı gerçek analiz mantığını içeren iş akışı ve `subworkflows/` altında bulunan bir dizi bakım iş akışı.
`demo.nf` iş akışı, `modules/` altında bulunan **modülleri** çağırır; bunlar gerçek analiz adımlarını gerçekleştirecek **süreçleri** içerir.

!!! note

    Subworkflow'lar bakım işlevleriyle sınırlı değildir ve süreç modüllerini kullanabilirler.

    Burada gösterilen `nf-core/demo` pipeline'ı spektrumun daha basit tarafında olmaktadır, ancak diğer nf-core pipeline'ları (örneğin `nf-core/rnaseq`) gerçek analizde yer alan subworkflow'ları kullanır.

Şimdi, bu bileşenleri sırayla gözden geçirelim.

#### 3.1.2. Giriş noktası betiği: `main.nf`

`main.nf` betiği, `nextflow run nf-core/demo` komutunu çalıştırdığımızda Nextflow'un başladığı giriş noktasıdır.
Bu, `nextflow run nf-core/demo` komutunu pipeline'ı çalıştırmak için çalıştırdığınızda, Nextflow'un otomatik olarak `main.nf` betiğini bulup çalıştırdığı anlamına gelir.
Bu, sadece nf-core pipeline'ları için değil, bu geleneksel adlandırma ve yapıyı takip eden herhangi bir Nextflow pipeline'ı için çalışır.

Bir giriş noktası betiği kullanmak, gerçek analiz betiği çalıştırılmadan önce ve sonra standartlaştırılmış 'bakım' subworkflow'larını çalıştırmayı kolaylaştırır.
Gerçek analiz iş akışını ve modüllerini gözden geçirdikten sonra bunları ele alacağız.

#### 3.1.3. Analiz betiği: `workflows/demo.nf`

`workflows/demo.nf` iş akışı, pipeline'ın merkezi mantığının saklandığı yerdir.
Normal bir Nextflow iş akışı gibi yapılandırılmıştır, ancak bir üst iş akışından çağrılmak üzere tasarlanmıştır, bu da birkaç ekstra özellik gerektirir.
İlgili farklılıkları, Hello Nextflow'dan basit Hello pipeline'ını nf-core uyumlu bir forma dönüştürmeyi ele aldığımızda bu kursun bir sonraki bölümünde ele alacağız.

`demo.nf` iş akışı, bir sonraki gözden geçireceğimiz `modules/` altında bulunan **modülleri** çağırır.

!!! note

    Bazı nf-core analiz iş akışları, alt seviye subworkflow'ları çağırarak ek iç içe geçme seviyeleri gösterir.
    Bu çoğunlukla, yaygın olarak birlikte kullanılan iki veya daha fazla modülü kolayca yeniden kullanılabilir pipeline segmentlerine sarmak için kullanılır.
    nf-core web sitesinde mevcut [nf-core subworkflow'larına](https://nf-co.re/subworkflows/) göz atarak bazı örnekler görebilirsiniz.

    Analiz betiği subworkflow'lar kullandığında, bunlar `subworkflows/` dizini altında saklanır.

#### 3.1.4. Modüller

Modüller, [Hello Nextflow eğitim kursunun 4. Bölümünde](../hello_nextflow/04_hello_modules.md) açıklandığı gibi süreç kodunun bulunduğu yerdir.

nf-core projesinde, modüller hem kökenlerini hem de içeriklerini yansıtan çok seviyeli iç içe bir yapı kullanılarak organize edilir.
En üst seviyede, modüller `nf-core` veya `local` (nf-core projesinin parçası değil) olarak ayrılır ve ardından sardıkları araç(lar)ın adını taşıyan bir dizine yerleştirilir.
Araç bir araç setine aitse (yani birden fazla araç içeren bir paket) o zaman araç setinin adını taşıyan bir ara dizin seviyesi vardır.

Bunu `nf-core/demo` pipeline modüllerine uygulandığını görebilirsiniz:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Dizin içeriği"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Burada `fastqc` ve `multiqc` modüllerinin `nf-core` modülleri içinde en üst seviyede oturduğunu görüyorsunuz, oysa `trim` modülü ait olduğu araç seti olan `seqtk` altında oturuyor.
Bu durumda `local` modül yok.

Süreci tanımlayan modül kod dosyası her zaman `main.nf` olarak adlandırılır ve şimdilik görmezden geleceğimiz testler ve `.yml` dosyalarıyla birlikte gelir.

Birlikte ele alındığında, giriş noktası iş akışı, analiz iş akışı ve modüller pipeline'ın 'ilginç' kısımlarını çalıştırmak için yeterlidir.
Ancak, orada bakım subworkflow'ları da olduğunu biliyoruz, bu yüzden şimdi onlara bakalım.

#### 3.1.5. Bakım subworkflow'ları

Modüller gibi, subworkflow'lar da `local` ve `nf-core` dizinlerine ayrılır ve her subworkflow kendi iç içe dizin yapısına, kendi `main.nf` betiğine, testlere ve `.yml` dosyasına sahiptir.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Dizin içeriği"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Yukarıda belirtildiği gibi, `nf-core/demo` pipeline'ı analize özgü herhangi bir subworkflow içermez, bu yüzden burada gördüğümüz tüm subworkflow'lar, adlarındaki `utils_` öneki ile belirtildiği gibi, sözde 'bakım' veya 'yardımcı program' iş akışlarıdır.
Bu subworkflow'lar, diğer yardımcı işlevlerin yanı sıra konsol çıktısında süslü nf-core başlığını üreten şeydir.

!!! tip

    Adlandırma desenlerinin yanı sıra, bu subworkflow'ların gerçekten analizle ilgili herhangi bir işlev gerçekleştirmediğinin bir başka göstergesi de hiç süreç çağırmamasıdır.

Bu, `nf-core/demo` pipeline'ını oluşturan temel kod bileşenlerinin özetini tamamlar.
Şimdi geliştirmeye dalmadan önce hakkında biraz bilmeniz gereken kalan unsurlara bir göz atalım: pipeline yapılandırması ve girdi doğrulama.

### 3.2. Pipeline yapılandırması

Daha önce Nextflow'un, girdiler ve parametreler, bilgi işlem kaynakları ve orkestrasyon ile ilgili diğer yönler açısından pipeline yürütmesini yapılandırmak için birçok seçenek sunduğunu öğrendiniz.
nf-core projesi, Nextflow'un esnek özelleştirme seçenekleri üzerine inşa etmeyi amaçlayan, pipeline'lar arasında daha fazla tutarlılık ve sürdürülebilirlik sağlayan yüksek düzeyde standartlaştırılmış pipeline yapılandırma yönergeleri uygular.

Merkezi yapılandırma dosyası `nextflow.config`, parametreler ve diğer yapılandırma seçenekleri için varsayılan değerleri ayarlamak için kullanılır.
Bu yapılandırma seçeneklerinin çoğu varsayılan olarak uygulanırken diğerleri (örneğin, yazılım bağımlılığı profilleri) isteğe bağlı profiller olarak dahil edilir.

`conf` klasöründe saklanan ve varsayılan olarak veya isteğe bağlı olarak profiller olarak yapılandırmaya eklenebilen birkaç ek yapılandırma dosyası vardır:

- `base.config`: Çoğu yüksek performanslı bilgi işlem ortamında genel kullanım için uygun bir 'boş sayfa' yapılandırma dosyası. Bu, örneğin modüllere uygulanması kolay olan geniş kaynak kullanım gruplarını tanımlar.
- `modules.config`: Ek modül yönergeleri ve argümanları.
- `test.config`: Demo pipeline'ını çalıştırdığımızda kullandığımız, minimum test verileriyle pipeline'ı çalıştırmak için bir profil.
- `test_full.config`: Pipeline'ı tam boyutlu bir test veri setiyle çalıştırmak için bir profil.

Kursta bu dosyalardan birkaçına değineceğiz.

### 3.3. Girdiler ve doğrulama

Daha önce `nf-core/demo` pipeline'ının test profilini incelediğimizde belirttiğimiz gibi, girdi olarak dosya yolları ve örnek tanımlayıcıları içeren bir samplesheet almak üzere tasarlanmıştır.
Dosya yolları `nf-core/test-datasets` deposunda bulunan gerçek verilere bağlıdır.

`assets` dizini altında da bir örnek samplesheet sağlanmıştır, ancak bundaki yollar gerçek değildir.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Bu özel samplesheet oldukça basittir, ancak bazı pipeline'lar birincil girdilerle ilişkili çok daha fazla meta veri içeren daha karmaşık samplesheet'ler üzerinde çalışır.

Ne yazık ki, bu dosyaları gözle kontrol etmek zor olabileceğinden, girdi verilerinin uygunsuz biçimlendirilmesi çok yaygın bir pipeline başarısızlığı kaynağıdır.
İlgili bir sorun, parametrelerin yanlış sağlanmasıdır.

Bu sorunların çözümü, tüm girdi dosyalarında beklenen bilgi türlerini içerdiklerinden ve doğru biçimlendirildiğinden emin olmak için otomatik doğrulama kontrolleri yapmak ve parametrelerin beklenen türde olduğundan emin olmaktır.
Buna girdi doğrulama denir ve ideal olarak bir pipeline'ı çalıştırmayı denemeden _önce_ yapılmalıdır, girdilerle ilgili bir sorun olduğunu öğrenmek için pipeline'ın başarısız olmasını beklemek yerine.

Yapılandırma gibi, nf-core projesi girdi doğrulama konusunda çok kararlıdır ve Nextflow pipeline'ları için kapsamlı doğrulama yetenekleri sağlayan bir Nextflow eklentisi olan [nf-schema eklentisinin](https://nextflow-io.github.io/nf-schema/latest/) kullanılmasını önerir.

Bu konuyu bu kursun 5. Bölümünde daha ayrıntılı olarak ele alacağız.
Şimdilik, bu amaç için sağlanan iki JSON dosyası olduğunu bilin: `nextflow_schema.json` ve `assets/schema_input.json`.

`nextflow_schema.json`, tür, açıklama ve yardım metni dahil olmak üzere pipeline parametreleri hakkında bilgileri makine tarafından okunabilir bir formatta saklamak için kullanılan bir dosyadır.
Bu, otomatik parametre doğrulama, yardım metni oluşturma ve UI arayüzlerinde etkileşimli parametre formu oluşturma dahil olmak üzere çeşitli amaçlar için kullanılır.

`schema_input.json`, girdi samplesheet yapısını tanımlamak için kullanılan bir dosyadır.
Her sütun, makine tarafından okunabilir bir formatta bir tür, desen, açıklama ve yardım metnine sahip olabilir.
Şema, otomatik doğrulama ve yararlı hata mesajları sağlama dahil olmak üzere çeşitli amaçlar için kullanılır.

### Özet

Bir nf-core pipeline'ının ana bileşenlerinin neler olduğunu ve kodun nasıl organize edildiğini; yapılandırmanın ana unsurlarının nerede bulunduğunu biliyorsunuz; ve girdi doğrulamanın ne için olduğunun farkındasınız.

### Sırada ne var?

Bir mola verin! Bu çok fazlaydı. Kendinizi tazelenmiş ve hazır hissettiğinizde, öğrendiklerinizi nf-core uyumlu bir pipeline yazmak için uygulamak üzere bir sonraki bölüme geçin.

!!! tip

    Bir sonraki bölüme geçmeden önce subworkflow'larla iş akışlarını nasıl oluşturacağınızı öğrenmek isterseniz, [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Yan Görevine göz atın.
