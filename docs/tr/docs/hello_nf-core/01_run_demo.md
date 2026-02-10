# Bölüm 1: Demo pipeline çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core eğitim kursunun bu ilk bölümünde, bir nf-core pipeline'ının nasıl bulunacağını ve deneneceğini, kodun nasıl düzenlendiğini ve [Hello Nextflow](../hello_nextflow/index.md)'da gösterilen düz Nextflow kodundan nasıl farklılaştığını öğreneceksiniz.

nf-core projesi tarafından kod yapısını ve araç işlemlerini göstermek amacıyla pipeline envanterinin bir parçası olarak sürdürülen nf-core/demo adlı bir pipeline kullanacağız.

Çalışma dizininizin [Başlarken](./00_orientation.md) sayfasında belirtildiği şekilde `hello-nf-core/` olarak ayarlandığından emin olun.

---

## 1. nf-core/demo pipeline'ını bulma ve edinme

[nf-co.re](https://nf-co.re) adresindeki proje web sitesinde nf-core/demo pipeline'ını bularak başlayalım. Bu site, genel dokümantasyon ve yardım makaleleri, her pipeline için dokümantasyon, blog yazıları, etkinlik duyuruları vb. gibi tüm bilgileri merkezileştirir.

### 1.1. Pipeline'ı web sitesinde bulma

Web tarayıcınızda [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) adresine gidin ve arama çubuğuna `demo` yazın.

![arama sonuçları](./img/search-results.png)

Pipeline dokümantasyon sayfasına erişmek için pipeline adına, `demo`, tıklayın.

Yayınlanan her pipeline, aşağıdaki dokümantasyon bölümlerini içeren özel bir sayfaya sahiptir:

- **Introduction:** Pipeline'ın tanıtımı ve genel bakışı
- **Usage:** Pipeline'ın nasıl çalıştırılacağına dair açıklamalar
- **Parameters:** Açıklamalarla birlikte gruplandırılmış pipeline parametreleri
- **Output:** Beklenen çıktı dosyalarının açıklamaları ve örnekleri
- **Results:** Tam test veri setinden oluşturulan örnek çıktı dosyaları
- **Releases & Statistics:** Pipeline sürüm geçmişi ve istatistikleri

Yeni bir pipeline'ı kullanmayı düşündüğünüzde, ne yaptığını ve çalıştırmayı denemeden önce nasıl yapılandırılması gerektiğini anlamak için önce pipeline dokümantasyonunu dikkatlice okumalısınız.

Şimdi göz atın ve şunları bulup bulamayacağınızı görün:

- Pipeline hangi araçları çalıştıracak (Sekmeyi kontrol edin: `Introduction`)
- Pipeline hangi girdileri ve parametreleri kabul ediyor veya gerektiriyor (Sekmeyi kontrol edin: `Parameters`)
- Pipeline tarafından üretilen çıktılar nelerdir (Sekmeyi kontrol edin: `Output`)

#### 1.1.1. Pipeline'a genel bakış

`Introduction` sekmesi, görsel bir temsil (metro haritası olarak adlandırılır) ve pipeline'ın bir parçası olarak çalıştırılan araçların listesi dahil olmak üzere pipeline'a genel bir bakış sağlar.

![pipeline metro haritası](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Örnek komut satırı

Dokümantasyon ayrıca örnek bir girdi dosyası (aşağıda daha ayrıntılı olarak tartışılacak) ve örnek bir komut satırı sağlar.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Örnek komutun bir workflow dosyası belirtmediğini, sadece pipeline deposuna referans olan `nf-core/demo`'yu içerdiğini fark edeceksiniz.

Bu şekilde çağrıldığında, Nextflow kodun belirli bir şekilde düzenlendiğini varsayacaktır.
Bu yapıyı inceleyebilmemiz için kodu edinelim.

### 1.2. Pipeline kodunu edinme

Pipeline'ın amacımıza uygun göründüğünü belirledikten sonra, deneyelim.
Neyse ki Nextflow, doğru biçimlendirilmiş depolardan pipeline'ları manuel olarak bir şey indirmeye gerek kalmadan kolayca edinmeyi sağlar.

Terminale dönelim ve şunu çalıştıralım:

```bash
nextflow pull nf-core/demo
```

??? success "Komut çıktısı"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow, pipeline kodunun bir `pull` işlemini yapar, yani tüm depoyu yerel diskinize indirir.

Açık olmak gerekirse, bunu sadece nf-core pipeline'larıyla değil, GitHub'da uygun şekilde kurulmuş herhangi bir Nextflow pipeline'ı ile yapabilirsiniz.
Ancak nf-core, Nextflow pipeline'larının en büyük açık kaynak koleksiyonudur.

Nextflow'dan bu şekilde edindiğiniz pipeline'ların listesini almanızı sağlayabilir:

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

Nextflow, indirilen kaynak kodunu, bu pipeline'ların doğrudan etkileşimde bulunacağınız kod yerine daha çok kütüphaneler gibi kullanılması gerektiği ilkesiyle kasıtlı olarak 'yolun dışında' tutar.

Ancak, bu eğitimin amaçları doğrultusunda, oraya göz atabilmek ve içinde ne olduğunu görmek istiyoruz.
Bunu kolaylaştırmak için, mevcut çalışma dizinimizden o konuma sembolik bir bağlantı oluşturalım.

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

### Çıkarım

Artık nf-core web sitesi üzerinden bir pipeline'ı nasıl bulacağınızı ve kaynak kodunun yerel bir kopyasını nasıl edineceğinizi biliyorsunuz.

### Sırada ne var?

Minimum çabayla bir nf-core pipeline'ını nasıl deneyeceğinizi öğrenin.

---

## 2. Pipeline'ı test profiliyle deneme

Kolaylık sağlamak için, her nf-core pipeline'ı bir test profiliyle birlikte gelir.
Bu, [nf-core/test-datasets](https://github.com/nf-core/test-datasets) deposunda barındırılan küçük bir test veri setini kullanarak pipeline'ın çalıştırılması için minimum yapılandırma ayarları kümesidir.
Küçük ölçekte bir pipeline'ı hızlıca denemenin harika bir yoludur.

!!! note

    Nextflow'un yapılandırma profil sistemi, farklı konteyner motorları veya çalıştırma ortamları arasında kolayca geçiş yapmanızı sağlar.
    Daha fazla ayrıntı için [Hello Nextflow Bölüm 6: Yapılandırma](../hello_nextflow/06_hello_config.md) bölümüne bakın.

### 2.1. Test profilini inceleme

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

    // Girdi verileri
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Üstteki yorum bloğunun, bu test profiliyle pipeline'ın nasıl çalıştırılacağını gösteren bir kullanım örneği içerdiğini hemen fark edeceksiniz.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Sağlamamız gereken tek şeyler örnek komutta açılı parantezler içinde gösterilenlerdir: `<docker/singularity>` ve `<OUTDIR>`.

Hatırlatmak gerekirse, `<docker/singularity>` konteyner sistemi seçimini ifade eder. Tüm nf-core pipeline'ları, tekrarlanabilirliği sağlamak ve yazılım kurulum sorunlarını ortadan kaldırmak için konteynerler (Docker, Singularity, vb.) ile kullanılabilir olacak şekilde tasarlanmıştır.
Bu nedenle pipeline'ı test etmek için Docker veya Singularity kullanmak isteyip istemediğimizi belirtmemiz gerekecek.

`--outdir <OUTDIR>` kısmı, Nextflow'un pipeline'ın çıktılarını yazacağı dizini ifade eder.
Bunun için uydurulabilecek bir ad vermemiz gerekiyor.
Zaten mevcut değilse, Nextflow çalışma zamanında bizim için oluşturacaktır.

Yorum bloğundan sonraki bölüme geçersek, test profili bize test için önceden yapılandırılmış olanları gösterir: en önemlisi, `input` parametresi zaten bir test veri setine işaret edecek şekilde ayarlanmıştır, bu nedenle kendi verilerimizi sağlamamıza gerek yoktur.
Önceden yapılandırılmış girdinin bağlantısını takip ederseniz, bunun birkaç deneysel örnek için örnek tanımlayıcıları ve dosya yolları içeren bir csv dosyası olduğunu göreceksiniz.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Buna samplesheet denir ve nf-core pipeline'larına en yaygın girdi biçimidir.

!!! note

    Veri formatlarına ve türlerine aşina değilseniz endişelenmeyin, takip edenler için önemli değil.

Bu, pipeline'ı denemek için ihtiyacımız olan her şeye sahip olduğumuzu doğrular.

### 2.2. Pipeline'ı çalıştırma

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
Pipeline'ın sürümünün, girdilerinin ve çıktılarının bir özetini ve birkaç yapılandırma öğesini içeren bir başlık vardır.

!!! note

    Çıktınız farklı zaman damgaları, çalıştırma adları ve dosya yolları gösterecektir, ancak genel yapı ve process çalıştırması benzer olmalıdır.

Çalıştırma çıktısına geçerek, hangi process'lerin çalıştırıldığını bize söyleyen satırlara bakalım:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Bu bize üç process'in çalıştırıldığını söyler; bunlar nf-core web sitesindeki pipeline dokümantasyon sayfasında gösterilen üç araca karşılık gelir: FASTQC, SEQTK_TRIM ve MULTIQC.

Burada gösterildiği gibi `NFCORE_DEMO:DEMO:MULTIQC` gibi tam process adları, tanıtıcı Hello Nextflow materyalinde görmüş olabileceğinizden daha uzundur.
Bunlar üst workflow'larının adlarını içerir ve pipeline kodunun modülerliğini yansıtır.
Buna birazdan daha ayrıntılı gireceğiz.

### 2.3. Pipeline çıktılarını inceleme

Son olarak, pipeline tarafından üretilen `demo-results` dizinine bakalım.

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
`nf-core/demo` pipeline'ının çıktıları hakkında daha fazla bilgi edinmek için [dokümantasyon sayfasına](https://nf-co.re/demo/1.0.2/docs/output/) bakın.

Bu aşamada, gözlemlenmesi gereken önemli şey, sonuçların modüle göre düzenlenmiş olması ve ek olarak pipeline çalıştırması hakkında çeşitli zaman damgalı raporlar içeren `pipeline_info` adlı bir dizinin bulunmasıdır.

Örneğin, `execution_timeline_*` dosyası hangi process'lerin çalıştırıldığını, hangi sırayla ve çalışmalarının ne kadar sürdüğünü gösterir:

![çalıştırma zaman çizelgesi raporu](./img/execution_timeline.png)

!!! note

    Burada görevler paralel olarak çalıştırılmadı çünkü Github Codespaces'te minimalist bir makine üzerinde çalışıyoruz.
    Bunların paralel olarak çalıştığını görmek için, codespace'inizin CPU tahsisini ve test yapılandırmasındaki kaynak sınırlarını artırmayı deneyin.

Bu raporlar tüm nf-core pipeline'ları için otomatik olarak oluşturulur.

### Çıkarım

Yerleşik test profili kullanarak bir nf-core pipeline'ını nasıl çalıştıracağınızı ve çıktılarını nerede bulacağınızı biliyorsunuz.

### Sırada ne var?

Pipeline kodunun nasıl düzenlendiğini öğrenin.

---

## 3. Pipeline kod yapısını inceleme

Pipeline'ı kullanıcılar olarak başarıyla çalıştırdığımıza göre, şimdi nf-core pipeline'larının dahili olarak nasıl yapılandırıldığına bakmak için bakış açımızı değiştirelim.

nf-core projesi, pipeline'ların nasıl yapılandırılacağı ve kod, yapılandırma ve dokümantasyonun nasıl düzenleneceği konusunda güçlü yönergeler uygular.
Bunun nasıl düzenlendiğini anlamak, bu kursun 2. Bölümünde ele alacağımız kendi nf-core uyumlu pipeline'larınızı geliştirmeye yönelik ilk adımdır.

Daha önce oluşturduğumuz `pipelines` sembolik bağlantısını kullanarak pipeline kodunun `nf-core/demo` deposunda nasıl düzenlendiğine bakalım.

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

İlk olarak, en üst düzeyde, özet bilgileri içeren bir README dosyası ve lisanslama, katkı yönergeleri, alıntı ve davranış kuralları gibi proje bilgilerini özetleyen yardımcı dosyalar bulabilirsiniz.
Ayrıntılı pipeline dokümantasyonu `docs` dizininde bulunur.
Tüm bu içerik, nf-core web sitesindeki web sayfalarını programatik olarak oluşturmak için kullanılır, bu nedenle her zaman kodla günceldir.

Şimdi, geri kalanlar için keşfimizi üç aşamaya böleceğiz:

1. Pipeline kod bileşenleri (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline yapılandırması
3. Girdiler ve doğrulama

Pipeline kod bileşenleriyle başlayalım.
Bireysel dosyalar içindeki koda dalmak yerine, dosya hiyerarşisine ve yapısal organizasyona odaklanacağız.

### 3.1. Pipeline kod bileşenleri

Standart nf-core pipeline kod organizasyonu, [Hello Nextflow](../hello_nextflow/index.md) kursunun 4. Bölümü olan [Hello Modüller](../hello_nextflow/04_hello_modules.md)'de tanıtıldığı gibi, kod yeniden kullanımını en üst düzeye çıkarmak için tasarlanmış modüler bir yapıyı takip eder, ancak gerçek nf-core tarzında bu, biraz ek karmaşıklıkla uygulanır.
Özellikle, nf-core pipeline'ları, yani bir üst workflow tarafından içe aktarılan workflow betiklerini bolca kullanır.

Bu biraz soyut gelebilir, bu yüzden bunun pratikte `nf-core/demo` pipeline'ında nasıl kullanıldığına bakalım.

!!! note

    Bu modüler bileşenlerin _nasıl_ bağlandığına dair gerçek kodu incelemeyeceğiz, çünkü subworkflow'ların kullanımıyla ilişkili, kafa karıştırıcı olabilecek bazı ek karmaşıklıklar vardır ve bunu anlamak eğitimin bu aşamasında gerekli değildir.
    Şimdilik, genel organizasyona ve mantığa odaklanacağız.

#### 3.1.1. Genel bakış

`nf-core/demo` pipeline'ı için ilgili kod bileşenleri arasındaki ilişkiler şöyle görünür:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

`main.nf` adlı, iki tür iç içe workflow için sarmalayıcı görevi gören bir _giriş noktası_ betiği vardır: `workflows/` altında bulunan ve `demo.nf` olarak adlandırılan gerçek analiz mantığını içeren workflow ve `subworkflows/` altında bulunan bir dizi bakım workflow'u.
`demo.nf` workflow'u `modules/` altında bulunan **modülleri** çağırır; bunlar gerçek analiz adımlarını gerçekleştirecek **process'leri** içerir.

!!! note

    Subworkflow'lar bakım işlevleriyle sınırlı değildir ve process modülleri kullanabilirler.

    Burada gösterilen `nf-core/demo` pipeline'ı spektrumun daha basit tarafında olmaktadır, ancak diğer nf-core pipeline'ları (`nf-core/rnaseq` gibi) gerçek analizde yer alan subworkflow'lar kullanır.

Şimdi, bu bileşenleri sırayla gözden geçirelim.

#### 3.1.2. Giriş noktası betiği: `main.nf`

`main.nf` betiği, `nextflow run nf-core/demo` çalıştırdığımızda Nextflow'un başladığı giriş noktasıdır.
Bu, pipeline'ı çalıştırmak için `nextflow run nf-core/demo` komutunu çalıştırdığınızda, Nextflow'un otomatik olarak `main.nf` betiğini bulup çalıştırdığı anlamına gelir.
Bu, sadece nf-core pipeline'ları için değil, bu geleneksel adlandırma ve yapıyı takip eden herhangi bir Nextflow pipeline'ı için çalışır.

Bir giriş noktası betiği kullanmak, gerçek analiz betiği çalıştırılmadan önce ve sonra standartlaştırılmış 'bakım' subworkflow'larını çalıştırmayı kolaylaştırır.
Bunları, gerçek analiz workflow'unu ve modüllerini gözden geçirdikten sonra ele alacağız.

#### 3.1.3. Analiz betiği: `workflows/demo.nf`

`workflows/demo.nf` workflow'u, pipeline'ın merkezi mantığının saklandığı yerdir.
Normal bir Nextflow workflow'u gibi yapılandırılmıştır, ancak bir üst workflow'dan çağrılmak üzere tasarlanmıştır, bu da birkaç ekstra özellik gerektirir.
Kursun bir sonraki bölümünde, Hello Nextflow'daki basit Hello pipeline'ının nf-core uyumlu bir forma dönüştürülmesiyle uğraştığımızda ilgili farklılıkları ele alacağız.

`demo.nf` workflow'u, daha sonra gözden geçireceğimiz `modules/` altında bulunan **modülleri** çağırır.

!!! note

    Bazı nf-core analiz workflow'ları, alt düzey subworkflow'ları çağırarak ek iç içe geçme düzeyleri sergiler.
    Bu çoğunlukla, yaygın olarak birlikte kullanılan iki veya daha fazla modülü kolayca yeniden kullanılabilir pipeline segmentlerine sarmak için kullanılır.
    nf-core web sitesinde mevcut [nf-core subworkflow'larına](https://nf-co.re/subworkflows/) göz atarak bazı örnekler görebilirsiniz.

    Analiz betiği subworkflow'lar kullandığında, bunlar `subworkflows/` dizini altında saklanır.

#### 3.1.4. Modüller

Modüller, [Hello Nextflow eğitim kursunun 4. Bölümünde](../hello_nextflow/04_hello_modules.md) açıklandığı gibi process kodunun bulunduğu yerdir.

nf-core projesinde modüller, hem kökenlerini hem de içeriklerini yansıtan çok düzeyli iç içe geçmiş bir yapı kullanılarak düzenlenir.
En üst düzeyde, modüller `nf-core` veya `local` (nf-core projesinin parçası değil) olarak ayırt edilir ve daha sonra sardıkları araç(lar)dan sonra adlandırılan bir dizine yerleştirilir.
Araç bir araç kitine aitse (yani birden fazla araç içeren bir paket) o zaman araç kitinden sonra adlandırılan bir ara dizin düzeyi vardır.

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

Burada `fastqc` ve `multiqc` modüllerinin `nf-core` modüllerinin üst düzeyinde oturduğunu görüyorsunuz, oysa `trim` modülü ait olduğu araç kiti olan `seqtk` altında oturuyor.
Bu durumda `local` modül yoktur.

Process'i tanımlayan modül kod dosyası her zaman `main.nf` olarak adlandırılır ve şimdilik göz ardı edeceğimiz testler ve `.yml` dosyalarıyla birlikte gelir.

Birlikte ele alındığında, giriş noktası workflow'u, analiz workflow'u ve modüller, pipeline'ın 'ilginç' kısımlarını çalıştırmak için yeterlidir.
Ancak, orada bakım subworkflow'ları da olduğunu biliyoruz, o yüzden şimdi onlara bakalım.

#### 3.1.5. Bakım subworkflow'ları

Modüller gibi, subworkflow'lar da `local` ve `nf-core` dizinlerine ayrılır ve her subworkflow'un kendi `main.nf` betiği, testleri ve `.yml` dosyası olan kendi iç içe dizin yapısı vardır.

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

Yukarıda belirtildiği gibi, `nf-core/demo` pipeline'ı herhangi bir analize özgü subworkflow içermez, bu nedenle burada gördüğümüz tüm subworkflow'lar, adlarındaki `utils_` öneki ile gösterildiği gibi, 'bakım' veya 'yardımcı' workflow'lar olarak adlandırılır.
Bu subworkflow'lar, diğer yardımcı işlevlerin yanı sıra konsol çıktısında süslü nf-core başlığını üreten şeydir.

!!! tip

    Adlandırma desenlerinin yanı sıra, bu subworkflow'ların gerçekten analizle ilgili herhangi bir işlev gerçekleştirmediğinin bir başka göstergesi de hiçbir process çağırmamasıdır.

Bu, `nf-core/demo` pipeline'ını oluşturan temel kod bileşenlerinin bir özetini tamamlar.
Şimdi geliştirmeye dalmadan önce hakkında biraz bilgi sahibi olmanız gereken kalan öğelere bakalım: pipeline yapılandırması ve girdi doğrulama.

### 3.2. Pipeline yapılandırması

Daha önce Nextflow'un, girdiler ve parametreler, hesaplama kaynakları ve orkestrasyon ile ilgili diğer yönler açısından pipeline çalıştırmasını yapılandırmak için birçok seçenek sunduğunu öğrendiniz.
nf-core projesi, Nextflow'un esnek özelleştirme seçenekleri üzerine inşa ederek pipeline'lar arasında daha fazla tutarlılık ve sürdürülebilirlik sağlamayı amaçlayan, pipeline yapılandırması için son derece standartlaştırılmış yönergeler uygular.

Merkezi yapılandırma dosyası `nextflow.config`, parametreler ve diğer yapılandırma seçenekleri için varsayılan değerleri ayarlamak için kullanılır.
Bu yapılandırma seçeneklerinin çoğu varsayılan olarak uygulanırken diğerleri (örneğin, yazılım bağımlılık profilleri) isteğe bağlı profiller olarak dahil edilir.

`conf` klasöründe saklanan ve varsayılan olarak veya isteğe bağlı olarak profil olarak yapılandırmaya eklenebilen birkaç ek yapılandırma dosyası vardır:

- `base.config`: Çoğu yüksek performanslı hesaplama ortamında genel kullanım için uygun bir 'boş sayfa' yapılandırma dosyası. Bu, örneğin modüllere uygulanması kolay olan geniş kaynak kullanımı gruplarını tanımlar.
- `modules.config`: Ek modül yönergeleri ve argümanları.
- `test.config`: Demo pipeline'ı çalıştırdığımızda kullandığımız, minimum test verileriyle pipeline'ı çalıştırmak için bir profil.
- `test_full.config`: Pipeline'ı tam boyutlu bir test veri setiyle çalıştırmak için bir profil.

Kursta bu dosyalardan birkaçına değineceğiz.

### 3.3. Girdiler ve doğrulama

Daha önce `nf-core/demo` pipeline'ının test profilini incelerken belirttiğimiz gibi, girdi olarak `nf-core/test-datasets` deposunda bulunan gerçek verilere bağlanan dosya yolları ve örnek tanımlayıcıları içeren bir samplesheet alacak şekilde tasarlanmıştır.

`assets` dizini altında da örnek bir samplesheet sağlanmıştır, ancak bunda bulunan yollar gerçek değildir.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Bu belirli samplesheet oldukça basittir, ancak bazı pipeline'lar, birincil girdilerle ilişkili çok daha fazla meta veriye sahip, daha karmaşık samplesheet'ler üzerinde çalışır.

Ne yazık ki, bu dosyaların göz ile kontrol edilmesi zor olabileceğinden, girdi verilerinin uygunsuz biçimlendirilmesi çok yaygın bir pipeline başarısızlık kaynağıdır.
İlgili bir sorun, parametrelerin yanlış sağlanmasıdır.

Bu sorunların çözümü, tüm girdi dosyalarında beklenen bilgi türlerini içerdiklerinden ve doğru biçimlendirildiğinden emin olmak için otomatik doğrulama kontrolleri çalıştırmak ve parametrelerin beklenen türde olduğundan emin olmak içindir.
Buna girdi doğrulama denir ve ideal olarak pipeline'ı çalıştırmayı denemeden _önce_ yapılmalıdır, girdilerde bir sorun olduğunu öğrenmek için pipeline'ın başarısız olmasını beklemek yerine.

Yapılandırma için olduğu gibi, nf-core projesi girdi doğrulama konusunda çok kararlıdır ve Nextflow pipeline'ları için kapsamlı doğrulama yetenekleri sağlayan bir Nextflow eklentisi olan [nf-schema eklentisinin](https://nextflow-io.github.io/nf-schema/latest/) kullanılmasını önerir.

Bu konuyu bu kursun 5. Bölümünde daha ayrıntılı olarak ele alacağız.
Şimdilik, bu amaç için sağlanan `nextflow_schema.json` ve `assets/schema_input.json` olmak üzere iki JSON dosyası olduğunun farkında olun.

`nextflow_schema.json`, tip, açıklama ve yardım metni dahil olmak üzere pipeline parametreleri hakkında bilgileri makine tarafından okunabilir bir formatta saklamak için kullanılan bir dosyadır.
Bu, otomatik parametre doğrulama, yardım metni oluşturma ve UI arayüzlerinde etkileşimli parametre formu oluşturma dahil olmak üzere çeşitli amaçlar için kullanılır.

`schema_input.json`, girdi samplesheet yapısını tanımlamak için kullanılan bir dosyadır.
Her sütun, makine tarafından okunabilir bir formatta bir tip, desen, açıklama ve yardım metnine sahip olabilir.
Şema, otomatik doğrulama ve yararlı hata mesajları sağlama dahil olmak üzere çeşitli amaçlar için kullanılır.

### Çıkarım

Bir nf-core pipeline'ının ana bileşenlerinin neler olduğunu ve kodun nasıl düzenlendiğini; yapılandırmanın ana öğelerinin nerede bulunduğunu biliyorsunuz; ve girdi doğrulamanın ne için olduğunun farkındasınız.

### Sırada ne var?

Bir mola verin! Bu çoktu. Kendinizi tazelenmiş ve hazır hissettiğinizde, öğrendiklerinizi bir nf-core uyumlu pipeline yazmak için uygulamak üzere bir sonraki bölüme geçin.

!!! tip

    Bir sonraki bölüme geçmeden önce subworkflow'larla workflow'ların nasıl oluşturulacağını öğrenmek isterseniz, [Workflow'ların Workflow'ları](../side_quests/workflows_of_workflows.md) Yan Görevi'ne göz atın.
