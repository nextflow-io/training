# Bölüm 1: Demo pipeline çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Hello nf-core eğitim kursunun bu ilk bölümünde, bir nf-core pipeline'ının nasıl bulunacağını ve deneneceğini, ihtiyaçlarınıza göre çalıştırmanın nasıl yapılandırılıp özelleştirileceğini ve girdi doğrulamanın yaygın hatalara karşı nasıl koruma sağladığını öğreneceksiniz.

nf-core projesi tarafından kod yapısını ve araç işlemlerini göstermek amacıyla pipeline envanterinin bir parçası olarak sürdürülen nf-core/demo adlı bir pipeline kullanacağız.

Çalışma dizininizin [Başlarken](./00_orientation.md) sayfasında belirtildiği şekilde `hello-nf-core/` olarak ayarlandığından emin olun.

---

## 1. nf-core/demo pipeline'ını bulma ve edinme

[nf-co.re](https://nf-co.re) adresindeki proje web sitesinde nf-core/demo pipeline'ını bularak başlayalım. Bu site; genel dokümantasyon ve yardım makaleleri, her pipeline için dokümantasyon, blog yazıları, etkinlik duyuruları ve benzerleri gibi tüm bilgileri merkezileştirir.

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

`Introduction` sekmesi, metro haritası (subway map) adı verilen görsel bir temsil ve pipeline'ın bir parçası olarak çalıştırılan araçların listesi dahil olmak üzere pipeline'a genel bir bakış sağlar.

![pipeline metro haritası](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Örnek komut satırı

Dokümantasyon ayrıca örnek bir girdi dosyası (aşağıda daha ayrıntılı olarak ele alınacak) ve örnek bir komut satırı sağlar.

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
Neyse ki Nextflow, doğru biçimlendirilmiş depolardan pipeline'ları manuel olarak herhangi bir şey indirmeye gerek kalmadan kolayca edinmeyi sağlar.

#### 1.2.1. `nextflow pull` kullanımı

Terminale dönelim ve şunu çalıştıralım:

```bash
nextflow pull nf-core/demo
```

??? success "Komut çıktısı"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

Nextflow, pipeline kodunun bir `pull` işlemini yapar; yani tüm depoyu yerel diskinize indirir.

Açık olmak gerekirse, bunu sadece nf-core pipeline'larıyla değil, GitHub'da uygun şekilde kurulmuş herhangi bir Nextflow pipeline'ı ile yapabilirsiniz.
Ancak nf-core, Nextflow pipeline'larının en büyük açık kaynak koleksiyonudur.

#### 1.2.2. `nextflow list` kullanımı

Nextflow'dan bu şekilde edindiğiniz pipeline'ların listesini alabilirsiniz:

```bash
nextflow list
```

??? success "Komut çıktısı"

    ```console
    nf-core/demo
    ```

Birden fazla pipeline listelendiğinde nasıl göründüğünü görmek için birkaç pipeline daha indirmeyi deneyebilirsiniz.

#### 1.2.3. Pipeline'larınızı `$NXF_HOME/assets/` dizininde bulma

Dosyaların mevcut çalışma dizininizde olmadığını fark edeceksiniz.
Varsayılan olarak, Nextflow bunları `$NXF_HOME/assets` dizinine kaydeder.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Directory contents"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note "Not"

    Eğitim ortamımızı kullanmıyorsanız, tam yol sisteminizde farklı olabilir.

Nextflow, indirilen kaynak kodunu kasıtlı olarak erişimi doğrudan olmayan bir konumda tutar; bu pipeline'ların doğrudan etkileşimde bulunacağınız kod yerine daha çok kütüphaneler gibi kullanılması gerektiği ilkesine dayanır.

#### 1.2.4. Kaynak koda kolay erişim için sembolik bağlantı oluşturma

Koda ayrıntılı olarak bakmayacağız; ancak genel organizasyonun nasıl göründüğüne dair bir fikir edinmek için hızlıca göz atalım.

Pipeline kaynak koduna göz atmayı kolaylaştırmak için assets dizinine sembolik bir bağlantı oluşturun:

```bash
ln -s $NXF_HOME/assets pipelines
```

Bu, `tree -L 2 pipelines` komutuyla kodu keşfetmenizi veya dosyaları doğrudan açmanızı sağlayan bir kısayol oluşturur.

#### 1.2.5. Kod organizasyonuna genel bakış

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

Orada çok şey oluyor; ancak bunların büyük çoğunluğu hakkında endişelenmenize gerek yok.

Kısaca belirtmek gerekirse, en üst düzeyde özet bilgileri içeren bir README dosyası ve lisanslama, katkı yönergeleri, alıntı ve davranış kuralları gibi proje bilgilerini özetleyen yardımcı dosyalar bulabilirsiniz.
Ayrıntılı pipeline dokümantasyonu `docs` dizininde bulunur.
Tüm bu içerik, nf-core web sitesindeki web sayfalarını programatik olarak oluşturmak için kullanılır; bu nedenle her zaman kodla günceldir.

Geri kalanlar için üç işlevsel kod dosyası grubunu ayırt edebiliriz:

1. Pipeline kod bileşenleri (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline yapılandırması
3. Pipeline parametreleri / girdiler ve doğrulama

Bu kursun bu bölümünde pipeline kod bileşenlerini ele almayacağız; ancak nf-core pipeline'larının son kullanıcısı olarak sizinle ilgili olabilecek yapılandırma ve doğrulama öğelerine değineceğiz.

!!! tip "İpucu"

    Herhangi bir nf-core pipeline'ının kaynak koduna GitHub üzerinden de göz atabilirsiniz; örneğin [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Her nf-core pipeline'ı aynı dizin yapısını takip eder; bu nedenle yapıyı bir kez öğrendikten sonra herhangi bir pipeline için yapılandırma dosyalarını, modülleri ve workflow'ları aynı şekilde bulabilirsiniz.

Ama şimdi, pipeline'ı çalıştırmaya geçelim!

### Özetle

Artık nf-core web sitesi üzerinden bir pipeline'ı nasıl bulacağınızı ve kaynak kodunun yerel bir kopyasını nasıl edineceğinizi biliyorsunuz.

### Sırada ne var?

Minimum çabayla bir nf-core pipeline'ını nasıl deneyeceğinizi öğrenin.

---

## 2. Pipeline'ı test profiliyle deneme

Kolaylık sağlamak için, her nf-core pipeline'ı bir test profiliyle birlikte gelir.
Bu, [nf-core/test-datasets](https://github.com/nf-core/test-datasets) deposunda barındırılan küçük bir test veri setini kullanarak pipeline'ın çalıştırılması için minimum yapılandırma ayarları kümesidir.
Küçük ölçekte bir pipeline'ı hızlıca denemenin harika bir yoludur.

!!! note "Not"

    Nextflow'un yapılandırma profil sistemi, farklı konteyner motorları veya çalıştırma ortamları arasında kolayca geçiş yapmanızı sağlar.
    Daha fazla ayrıntı için [Hello Nextflow Bölüm 6: Yapılandırma](../hello_nextflow/06_hello_config.md) bölümüne bakın.

### 2.1. Test profilini inceleme

Bir pipeline'ın test profilinin çalıştırmadan önce ne belirttiğini kontrol etmek iyi bir uygulamadır.
`nf-core/demo` için `test` profili `conf/test.config` yapılandırma dosyasında bulunur.
`nextflow pull` ile indirilen pipeline kaynağının içinde yerel olarak bulabilirsiniz:

```bash
code $NXF_HOME/assets/nf-core/demo/conf/test.config
```

Bu dosyanın içeriği aşağıda gösterilmiştir:

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
        cpus: 2,
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

Sağlamamız gereken tek şeyler, örnek komutta açılı parantezler (yer tutucular) içinde gösterilenlerdir: `<docker/singularity>` ve `<OUTDIR>`.

Hatırlatmak gerekirse, `<docker/singularity>` konteyner sistemi seçimini ifade eder. Tüm nf-core pipeline'ları, tekrarlanabilirliği sağlamak ve yazılım kurulum sorunlarını ortadan kaldırmak için konteynerler (Docker, Singularity vb.) ile kullanılabilir olacak şekilde tasarlanmıştır.
Bu nedenle pipeline'ı test etmek için Docker veya Singularity kullanmak isteyip istemediğimizi belirtmemiz gerekecek.

`--outdir <OUTDIR>` kısmı, Nextflow'un pipeline'ın çıktılarını yazacağı dizini ifade eder.
Sizin belirleyeceğiniz bir isim vermemiz gerekiyor.
Dizin henüz mevcut değilse, Nextflow çalışma zamanında bizim için oluşturacaktır.

Yorum bloğundan sonraki bölüme geçersek, test profili bize test için önceden yapılandırılmış olanları gösterir: en önemlisi, `input` parametresi zaten bir test veri setine işaret edecek şekilde ayarlanmıştır; bu nedenle kendi verilerimizi sağlamamıza gerek yoktur.
Önceden yapılandırılmış girdinin bağlantısını takip ederseniz, bunun birkaç deneysel örnek için örnek tanımlayıcıları ve dosya yolları içeren bir CSV dosyası olduğunu göreceksiniz.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Buna samplesheet denir ve nf-core pipeline'larına en yaygın girdi biçimidir.

!!! note "Not"

    Veri formatlarına ve türlerine aşina değilseniz endişelenmeyin; takip edenler için önemli değil.

Bu, pipeline'ı denemek için ihtiyacımız olan her şeye sahip olduğumuzu doğrular.

### 2.2. Pipeline'ı çalıştırma

Konteyner sistemi için Docker'ı ve çıktı dizini olarak `demo-results`'ı kullanmaya karar verelim; test komutunu çalıştırmaya hazırız:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
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

!!! note "Not"

    Çıktınız farklı zaman damgaları, çalıştırma adları ve dosya yolları gösterecektir; ancak genel yapı ve süreç çalıştırması benzer olmalıdır.

Çıktının üst kısmındaki şu satıra dikkat edin:

```console
Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: 45904cb9d1 [master]
```

Bu satır, pipeline'ın hangi revizyonunun kullanıldığını gösterir.
Bir sürüm belirtmediğimiz için Nextflow, `master` dalındaki en son commit'i kullandı.
Tekrarlanabilir çalıştırmalar için `-r` bayrağını kullanarak belirli bir sürümü sabitlemelisiniz:

```bash
nextflow run nf-core/demo -r 1.1.0 -profile docker,test --outdir demo-results
```

Bu, yeni commit'ler veya sürümler yayınlansa da her seferinde aynı pipeline kodunun kullanılmasını sağlar.
Bu eğitimde basitlik adına `-r` bayrağını atlıyoruz; ancak üretim ortamında her zaman belirtmelisiniz.

Çalıştırma çıktısına geçerek, hangi süreçlerin çalıştırıldığını bize söyleyen satırlara bakalım:

```console
executor >  local (7)
[ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Bu bize üç sürecin çalıştırıldığını söyler; bunlar nf-core web sitesindeki pipeline dokümantasyon sayfasında gösterilen üç araca karşılık gelir: FASTQC, SEQTK_TRIM ve MULTIQC.

Burada gösterildiği gibi `NFCORE_DEMO:DEMO:MULTIQC` şeklindeki tam süreç adları, tanıtıcı Hello Nextflow materyalinde görmüş olabileceğinizden daha uzundur.
Bunlar üst iş akışlarının adlarını içerir ve pipeline kodunun modülerliğini yansıtır.
Buna bu kursun 2. Bölümünde daha ayrıntılı gireceğiz.

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
`nf-core/demo` pipeline'ının çıktıları hakkında daha fazla bilgi edinmek için [dokümantasyon sayfasına](https://nf-co.re/demo/1.1.0/docs/output/) bakın.

Bu aşamada, gözlemlenmesi gereken önemli şey, sonuçların modüle göre düzenlenmiş olması ve ek olarak pipeline çalıştırması hakkında çeşitli zaman damgalı raporlar içeren `pipeline_info` adlı bir dizinin bulunmasıdır.

Örneğin, `execution_timeline_*` dosyası hangi süreçlerin çalıştırıldığını, hangi sırayla ve çalışmalarının ne kadar sürdüğünü gösterir:

![çalıştırma zaman çizelgesi raporu](./img/execution_timeline.png)

!!! note "Not"

    Burada görevler paralel olarak çalıştırılmadı; çünkü Github Codespaces'te minimalist bir makine üzerinde çalışıyoruz.
    Bunların paralel olarak çalıştığını görmek için, codespace'inizin CPU tahsisini ve test yapılandırmasındaki kaynak sınırlarını artırmayı deneyin.

Bu raporlar tüm nf-core pipeline'ları için otomatik olarak oluşturulur.

### Özetle

Yerleşik test profili kullanarak bir nf-core pipeline'ını nasıl çalıştıracağınızı ve çıktılarını nerede bulacağınızı biliyorsunuz.

### Sırada ne var?

Pipeline'ı çalıştırmayı özelleştirmek için nasıl yapılandıracağınızı öğrenin.

---

## 3. Pipeline çalıştırmasını yapılandırma

[Hello Config](../hello_nextflow/06_hello_config.md) bölümünde açıklandığı gibi, pipeline kodunu değiştirmeden pipeline'ın hangi veriler üzerinde ve nasıl çalışacağını değiştirebilmek istiyoruz.
Bu amaçla Nextflow, pipeline yapılandırmasını kontrol etmenin birden fazla yolunu destekler; bu durum başlangıçta bunaltıcı gelebilir.

nf-core projesi, yapılandırma öğelerini düzenlemek için kurallar belirler ve en üst düzeyde iki tür yapılandırmayı birbirinden ayırt eder: **pipeline parametreleri** ve dar anlamda **yapılandırma**.

- **Pipeline parametreleri** (`params` sistemi aracılığıyla ayarlanır): Genellikle girdi dosyaları, araç davranış bayrakları ve analiz parametrelerini içerir.
- Dar anlamda **yapılandırma**: Pipeline'ın nasıl çalıştırıldığına ilişkin lojistiği ifade eder; yani yürütücü, hesaplama kaynağı tahsisleri ve benzerleri.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/params_vs_config.excalidraw.svg"
</figure>

Pipeline parametrelerini ele alarak başlayalım, ardından dar anlamda yapılandırmaya bakacağız.

### 3.1. Pipeline parametreleri

Tüm nf-core pipeline'larında, `--help` bayrağını kullanarak komut satırından doğrudan pipeline parametrelerinin tam listesini alabilirsiniz; bu bayrağın kendisi de bir pipeline parametresidir.

#### 3.1.1. `--help` ile parametre listesini alma

Demo pipeline için yardım komutunu çalıştırın:

```bash
nextflow run nf-core/demo --help
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `https://github.com/nf-core/demo` [run_name] DSL2 - revision: 45904cb9d1 [master]

    ----------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.1.0
    ----------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>

    Input/output options
      --input                       [string]           Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string]           The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string]           Email address for completion summary.
      --multiqc_title               [string]           MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string]           Name of iGenomes reference.
      --fasta                       [string]           Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean]          Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]           Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string]  Display the help message.
      --help_full                   [boolean]          Display the full detailed help message.
      --show_hidden                 [boolean]          Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 20 param(s), use the `--show_hidden` parameter to show them !!
    ----------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Görüldüğü gibi, çıktı parametreleri kategorilere göre gruplandırır (Input/output options, Reference genome options vb.) ve her biri için tür ve açıklama sağlar.

Bu kategorilendirme, aşağıda daha ayrıntılı ele alınan bir şema dosyası tarafından belirlenir.
Yalın Nextflow pipeline'larında `--help`, yalnızca geliştirici bunu manuel olarak uyguladıysa çalışır.

!!! tip "İpucu"

    `--publish_dir_mode` veya `--monochrome_logs` gibi varsayılan olarak gizlenen ek parametreleri görmek için `--help --show_hidden` kullanın.

#### 3.1.2. Parametre değerlerini ayarlama

[Hello Config](../hello_nextflow/06_hello_config.md) bölümünde ele alındığı gibi, parametre değerlerini komut satırında `--param_name` ile ayarlayabilir veya bir dizi parametreyi YAML dosyasında toplayıp `-params-file` ile geçirebilirsiniz.
Her iki yaklaşım da nf-core pipeline'larıyla aynı şekilde çalışır.

Örneğin, kırpma adımını atlamak için:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-notrim --skip_trim
```

??? success "Komut çıktısı"

    ```console
    executor >  local (4)
    [3f/a82c91] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [7d/c5e014] NFCORE_DEMO:DEMO:MULTIQC             | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`SEQTK_TRIM` süreci artık çıktıda görünmüyor.

!!! info "Bilgi"

    `-c` ile geçirilen özel bir yapılandırma dosyasında pipeline parametrelerini ayarlamak teknik olarak mümkün olsa da, Nextflow'un yapılandırma öncelik kurallarına bağlı olarak bu, pipeline'ın kendi `nextflow.config` dosyasında zaten ayarlanmış varsayılanları geçersiz kılmayabilir.
    Komut satırında `--param_name` veya `-params-file` kullanmak daha güvenilirdir; çünkü bunlar her zaman öncelik taşır.

    **Pratik kural olarak:** `--help` çıktısında görünüyorsa, bir yapılandırma dosyası yerine komut satırı veya params dosyası aracılığıyla ayarlayın.

#### 3.1.3. Parametre doğrulama

İlginç bir bilgi: `--help` komutu tüm nf-core pipeline'larında çalışır; çünkü nf-core projesi, geliştiricilerin tüm pipeline parametrelerini bir JSON şema dosyasında (`nextflow_schema.json`) resmi olarak tanımlamasını zorunlu kılar.
Bu şema, her parametrenin türünü, açıklamasını, varsayılan değerini ve gruplandırmasını kaydeder.

`--help` çıktısını desteklemenin yanı sıra, şema dosyası başlatma sırasında otomatik doğrulamayı da etkinleştirir.
Bu, Nextflow'un geçirdiğiniz her parametrenin var olup olmadığını ve uygun bir değer verilip verilmediğini (uygun türde, izin verilen değer aralığında vb.) kontrol edebildiği anlamına gelir.

Bunu [Bölüm 5: Girdi Doğrulama](05_input_validation.md) bölümünde daha ayrıntılı ele alıyoruz; ancak demo pipeline'a geçersiz parametre girdisi vererek bunu şimdiden uygulamada görebilirsiniz.

##### 3.1.3.1. Tanınmayan parametreler

Var olmayan bir parametre geçirmeyi deneyin:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --foobar "invalid"
```

Konsol çıktısı bir uyarı içerir:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

Pipeline çalışmaya devam eder; ancak uyarı, `--foobar`'ın tanınan bir parametre olmadığını hemen bildirir.
Bu, `--outdir` yerine `--outDir` gibi yazım hatalarını, çıktının neden yanlış yere gittiğini merak ederek hesaplama zamanı harcamadan önce yakalar.

##### 3.1.3.2. Geçersiz parametre değerleri

Doğrulama, parametre **değerlerini** de kontrol eder.
`--skip_trim` parametresi bir boolean bayraktır; bu nedenle string bir değer geçirilmesi pipeline'ın hemen başarısız olmasına neden olur:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]
```

Herhangi bir süreç çalışmadan önce pipeline durur; bu sizi başarısız veya hatalı bir çalıştırmadan korur.
Boolean parametreler değer olmadan bayrak olarak (`--skip_trim`) geçirilmeli ya da params dosyasında `true`/`false` olarak ayarlanmalıdır.

#### 3.1.4. Girdi doğrulama

Aynı doğrulama mantığı, girdi dosyalarının geçerliliğini kontrol etmek için de kullanılabilir.
Örneğin, bir pipeline ana veri girdisi olarak samplesheet bekliyorsa (birçok nf-core pipeline'ında bu durum geçerlidir), geliştirici girdi dosyasının nasıl yapılandırılması gerektiğini açıklayan bir girdi şeması (parametre şemasından ayrı) sağlayabilir.

Ardından çalışma zamanında Nextflow, sağlanan girdi dosyasının geçerli olup olmadığını kontrol edebilir.

Bunu da [Bölüm 5: Girdi Doğrulama](05_input_validation.md) bölümünde daha ayrıntılı ele alıyoruz; ancak demo pipeline'a geçersiz bir girdi samplesheet'i vererek bunu şimdiden uygulamada görebilirsiniz.

`nf-core/demo` pipeline'ı `sample`, `fastq_1` ve `fastq_2` sütunlarına sahip bir CSV dosyası bekler.
Bu, beklenen yapıyı, sütun türlerini ve kısıtlamaları belirten bir şema dosyasında (`assets/schema_input.json`) tanımlanmıştır.

??? abstract "assets/schema_input.json"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

Şema, `sample` ve `fastq_1`'in zorunlu olduğunu; `fastq_2`'nin ise isteğe bağlı olduğunu (hem çift uçlu hem de tek uçlu verileri destekler) belirtir.
Dosya yolları, varlık ve uzantı deseni açısından doğrulanır.

##### 3.1.4.1. Geçersiz bir samplesheet oluşturma

Eksik bir sütun ve var olmayan bir dosya yolu içeren bir samplesheet oluşturun:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Bu samplesheet, zorunlu `fastq_1` sütununu içermiyor ve `fastq_2`'de var olmayan bir dosya yolu barındırıyor.
Her iki sorun da bir sonraki adımda doğrulama hatası üretecektir.

##### 3.1.4.2. Demo pipeline'ı geçersiz samplesheet ile çalıştırma

Demo pipeline'ı `malformed_samplesheet.csv` dosyasını girdi olarak kullanarak çalıştırın.

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory
       '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces
       and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1
```

Görüldüğü gibi, pipeline hemen başarısız olur ve **tüm** doğrulama hatalarını aynı anda raporlar.
nf-schema ilk hatada durmaz; tüm sorunları toplar ve birlikte listeler; böylece sorunları tek tek keşfetmek yerine hepsini bir seferde düzeltebilirsiniz.

Her hata, soruna neden olan tam girişi ve alanı tanımlar; böylece samplesheet'inizi düzeltip pipeline'ı yeniden başlatabilirsiniz. Nextflow'un dosya yoluna gerçekten erişmeye çalıştığı ilerleyen bir aşamada başarısız olacağından endişe etmeden.

Geliştiriciler için tüm bunlar bu kursun [Bölüm 5](./05_input_validation.md)'inde daha ayrıntılı ele alınmaktadır.

### 3.2. Yapılandırma

Dar anlamda yapılandırma, pipeline'ın **nasıl** çalıştığını kontrol eder: kaynak tahsisi, araca özgü argümanlar, görevlerin nerede çalıştırıldığı ve hangi yazılım paketleme sisteminin kullanılacağı.

nf-core pipeline'ları, `nextflow.config` ve `conf/` dizininde varsayılan yapılandırmayı içerir.
Herhangi bir şeyi geçersiz kılmadan önce, varsayılanların nerede bulunduğunu bilmek faydalıdır.

2.1. bölümünde pipeline kaynak kodunun `$NXF_HOME/assets` dizininde bulunduğunu gördünüz.
Mevcut yapılandırma dosyalarını listelemek için:

```bash
ls $NXF_HOME/assets/nf-core/demo/conf/
```

```console
base.config  igenomes.config  igenomes_ignored.config  modules.config  test.config  test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/nfcore_config_files.excalidraw.svg"
</figure>

En önemli yapılandırma dosyaları şunlardır:

- **`conf/base.config`**: Süreçlere CPU, bellek ve zaman atayan kaynak etiketlerini (`process_low`, `process_medium`, `process_high`) tanımlar. Bir sürecin beklenenden fazla kaynak kullandığını gördüğünüzde, bu varsayılanlar buradan gelir.
- **`conf/modules.config`**: Süreç başına araç argümanlarını (`ext.args`) ve çıktı yayımlama ayarlarını (`publishDir`) belirler. Her aracın varsayılan olarak hangi argümanları aldığını görmek için bu dosyayı açın.
- **`conf/test.config`**: 2.1. bölümünde kullandığınız test profili; `resourceLimits` aracılığıyla kaynakları sınırlar ve bir test samplesheet'i ayarlar. `-profile test` ile etkinleştirilir.
  Tam boyutlu bir test veri setiyle çalıştırmak için, kıyaslama açısından kullanışlı olan `conf/test_full.config` de mevcuttur.

Merkezi `nextflow.config`, yukarıdakilerin tümünü yükler ve her şey için uygun varsayılan değerleri ayarlar.

Bu dosyalarda belirtilen ayarlardan herhangi birini değiştirmek isterseniz, bu dosyaları doğrudan değiştirmeyin.
Bunun yerine kendi yapılandırma dosyanızı oluşturun ve `-c` ile geçirin.
Belirttiğiniz değerler, diğer dosyalarda ayarlanan varsayılan değerleri geçersiz kılar.

Bunu pratikte yapmak için birkaç alıştırma üzerinden geçelim.

#### 3.2.1. Bir süreç için kaynak tahsisini değiştirme

Demo pipeline, `base.config` dosyasında tanımlanan etiketleri kullanarak kaynakları atar.
Örneğin, `FASTQC` süreci 6 CPU ve 36 GB bellek tahsis eden `process_medium` etiketini kullanır.

Test profili kaynakları `resourceLimits` aracılığıyla sınırlar; ancak belirli süreçler için kaynakları da geçersiz kılabilirsiniz.

`custom.config` adlı bir dosya oluşturun:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
}
```

Pipeline'ı özel yapılandırmanızla çalıştırın:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-custom -c custom.config
```

??? success "Komut çıktısı"

    ```console
    executor >  local (7)
    [2a/f17b3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [9c/e4d028] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [5b/a93c71] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

`-c` bayrağı, yapılandırmanızı pipeline'ın yerleşik yapılandırmasının üzerine ekler.

#### 3.2.2. `ext.args` ile araç argüman değerlerini ayarlama

Birçok komut satırı aracının, çok yaygın kullanılmadıkça pipeline parametresi olarak ayarlanmayan isteğe bağlı argümanları vardır.
Bu araç argümanları için nf-core modülleri, argümanları bir yapılandırma dosyası aracılığıyla temel araca geçirmek amacıyla `ext.args` adlı bir Nextflow kuralını kullanır.

Örneğin, `ext.args` kullanarak `SEQTK_TRIM` modülüne bir kırpma argümanı ekleyelim.

##### 3.2.2.1. Özel yapılandırmayı güncelleme

`custom.config` dosyanızı güncelleyin:

```groovy title="custom.config" linenums="1" hl_lines="6 7 8"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

Bu, `seqtk trimfq`'ya kalite kırpmasına ek olarak her okumanın başından 5 baz kırpmasını söyler.

##### 3.2.2.2. Pipeline'ı çalıştırma

Etkisini görmek için pipeline'ı bu yapılandırmayla tekrar çalıştırın:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results-extargs -c custom.config
```

??? success "Komut çıktısı"

    ```console
    executor >  local (7)
    [1e/b7a392] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [ab/cd1234] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [4f/c8d105] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Argümanın uygulandığını doğrulamak için, çalıştırma çıktısından `SEQTK_TRIM` work dizini hash'ini bulun (örn. `work/ab/cd1234...`) ve içindeki `.command.sh` dosyasını kontrol edin:

```bash
cat work/ab/cd1234/.command.sh
```

??? success "Komut çıktısı"

    ```console
    #!/usr/bin/env bash
    ...
    seqtk trimfq -b 5 SAMPLE3_SE.fastq.gz | gzip -c > SAMPLE3_SE.trimmed.fastq.gz
    ```

`seqtk trimfq` komutunda `-b 5`'i görmelisiniz; bu, `ext.args` geçersiz kılmanızın etkili olduğunu doğrular.

##### 3.2.2.3. Varsayılan değerleri geçersiz kılma

Bazı modüllerin `ext.args` değerleri varsayılan olarak zaten ayarlanmıştır.
Örneğin, `FASTQC` modülü varsayılan olarak `ext.args = '--quiet'` ile yapılandırılmıştır (`conf/modules.config` dosyasında tanımlanmıştır).

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}"
        ]
    }
```

Özel bir yapılandırma dosyası aracılığıyla `ext.args` için bir değer sağlarsanız, bu değer söz konusu süreç için ayarlanan varsayılanın tamamen yerini alır.

Örneğin, varsayılan `'--quiet'` iken `ext.args = '--kmers 8'` ayarlarsanız, `--quiet` bayrağı artık uygulanmayacaktır.
Her ikisini de korumak için `ext.args = '--quiet --kmers 8'` olarak ayarlayın.

Bu, `ext.args` ile argüman değerleri sağlamak istediğiniz araçların varsayılan yapılandırmasını kontrol etmekten sorumlu olduğunuz anlamına gelir.

### Özetle

Bir nf-core pipeline'ından nasıl yardım alacağınızı, parametreleri nasıl ayarlayacağınızı ve bunların nasıl doğrulandığını anladığınızı; yapılandırma dosyaları aracılığıyla yapılandırmayı nasıl özelleştireceğinizi biliyorsunuz.

### Sırada ne var?

Bir mola verin! Hazır hissettiğinizde, kendi nf-core uyumlu pipeline'ınızı sıfırdan oluşturacağınız 2. Bölüme geçin.
