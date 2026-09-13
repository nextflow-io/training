# Bölüm 1: Demo Pipeline'ı Çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Use nf-core eğitim kursunun bu ilk bölümünde, bir nf-core pipeline'ının nasıl bulunacağını ve yerleşik test profili kullanılarak nasıl deneneceğini gösteriyoruz.

nf-core projesi tarafından, demonstrasyon ve eğitim amaçlı pipeline envanterinin bir parçası olarak sürdürülen nf-core/demo adlı bir pipeline kullanacağız.

Çalışma dizininizin [Başlarken](./00_orientation.md) sayfasında belirtildiği şekilde `nfcore-use/` olarak ayarlandığından emin olun.

---

## 1. nf-core/demo Pipeline'ını Bulma ve İndirme

Tüm belgeleri, yardım makalelerini, her pipeline'a ait dokümantasyonu, blog yazılarını, etkinlik duyurularını ve benzeri içerikleri merkezileştiren proje web sitesi [nf-co.re](https://nf-co.re) üzerinden nf-core/demo pipeline'ını bulmaya başlayalım.

### 1.1. Pipeline'ı Web Sitesinde Bulma

Web tarayıcınızda [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) adresine gidin ve arama çubuğuna `demo` yazın.

![arama sonuçları](./img/search-results.png)

Pipeline dokümantasyon sayfasına erişmek için pipeline adına, yani `demo`'ya tıklayın.

Yayımlanan her pipeline'ın aşağıdaki dokümantasyon bölümlerini içeren özel bir sayfası bulunur:

- **Introduction:** Pipeline'a giriş ve genel bakış
- **Usage:** Pipeline'ın nasıl çalıştırılacağına dair açıklamalar
- **Parameters:** Açıklamalarıyla birlikte gruplandırılmış pipeline parametreleri
- **Output:** Beklenen çıktı dosyalarının açıklamaları ve örnekleri
- **Results:** Tam test veri setinden üretilen örnek çıktı dosyaları
- **Releases & Statistics:** Pipeline sürüm geçmişi ve istatistikler

Yeni bir pipeline benimsemeyi düşündüğünüzde, çalıştırmayı denemeden önce pipeline dokümantasyonunu dikkatlice okuyarak ne yaptığını ve nasıl yapılandırılması gerektiğini anlamalısınız.

Şimdi bir göz atın ve şunları bulmaya çalışın:

- Pipeline'ın hangi araçları çalıştıracağı (`Introduction` sekmesini kontrol edin)
- Pipeline'ın hangi girdileri ve parametreleri kabul ettiği veya gerektirdiği (`Parameters` sekmesini kontrol edin)
- Pipeline tarafından üretilen çıktılar (`Output` sekmesini kontrol edin)

#### 1.1.1. Pipeline'a Genel Bakış

`Introduction` sekmesi, pipeline'ın metro haritası (subway map) adı verilen görsel bir temsili ve pipeline'ın parçası olarak çalıştırılan araçların listesi dahil olmak üzere genel bir bakış sunar.

![pipeline metro haritası](./img/nf-core-demo-subway-cropped.png)

1. Okuma kalite kontrolü ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adaptör ve kalite kırpma ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Ham okumalar için kalite kontrolü sunumu ([MULTIQC](http://multiqc.info/))
4. Bir inek aracılığıyla eğlenceli bir metin mesajı oluşturma ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Örnek Komut Satırı

Dokümantasyon ayrıca bir örnek girdi dosyası (aşağıda daha ayrıntılı ele alınmaktadır) ve bir örnek komut satırı sağlar.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Örnek komutun bir iş akışı dosyası belirtmediğini, yalnızca pipeline deposuna atıfta bulunduğunu fark edeceksiniz: `nf-core/demo`.

Bu şekilde çağrıldığında Nextflow, kodun belirli bir şekilde düzenlendiğini varsayar.
Bu yapıyı inceleyebilmek için kodu indirelim.

### 1.2. Pipeline Kodunu İndirme

Pipeline'ın amaçlarımıza uygun göründüğünü belirledikten sonra deneyelim.
Neyse ki Nextflow, doğru biçimde yapılandırılmış depolardan pipeline'ları manuel olarak herhangi bir şey indirmenize gerek kalmadan kolayca almanızı sağlar.

#### 1.2.1. `nextflow pull` Kullanımı

Terminale dönelim ve aşağıdaki komutu çalıştıralım:

```bash
nextflow pull nf-core/demo
```

??? success "Komut çıktısı"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

Nextflow, pipeline kodunu `pull` ederek tam depoyu yerel sürücünüze indirir.

Bunu, yalnızca nf-core pipeline'larıyla değil, GitHub'da uygun şekilde kurulmuş herhangi bir Nextflow pipeline'ıyla yapabileceğinizi belirtmek gerekir.
Ancak nf-core, en büyük açık kaynaklı Nextflow pipeline koleksiyonudur.

#### 1.2.2. `nextflow list` Kullanımı

Nextflow'un bu şekilde indirdiğiniz pipeline'ların listesini size vermesini sağlayabilirsiniz:

```bash
nextflow list
```

??? success "Komut çıktısı"

    ```console
    nf-core/demo
    ```

Birden fazla pipeline'ınız olduğunda nasıl listelendiğini görmek için birkaç pipeline daha indirmeyi deneyebilirsiniz.

#### 1.2.3. Pipeline'ın İndirildiği Yeri Bulma

Dosyaların mevcut çalışma dizininizde olmadığını fark edeceksiniz.
Nextflow, varsayılan olarak indirilen pipeline'ları `$NXF_HOME/assets` altına kaydeder.

Belirli bir pipeline'ın nerede bulunduğunu öğrenmek için doğrudan Nextflow'a sorun:

```bash
nextflow info nf-core/demo
```

??? success "Komut çıktısı"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Bilgi"

    Eğitim ortamamızı kullanmıyorsanız tam yol sisteminizde farklı olabilir.

Nextflow, indirilen kaynak kodu kasıtlı olarak erişimi doğrudan olmayan bir yerde tutar; çünkü bu pipeline'ların doğrudan etkileşime gireceğiniz koddan ziyade kütüphaneler gibi kullanılması gerektiği ilkesini benimser.

Arka planda Nextflow, indirilen her pipeline'ı `$NXF_HOME/assets/.repos/` altında bir git deposu olarak saklar ve her revizyon için kodu `clones/<commit>/` alt dizinine çıkarır.
`.repos` gizli bir dizin olduğundan, düz bir `tree -L 2 $NXF_HOME/assets/` komutu boş görünecektir.

#### 1.2.4. Kaynak Koda Kolayca Erişmek için Sembolik Bağlantı Oluşturma

Koda ayrıntılı olarak bakmayacağız; ancak genel organizasyonun nasıl göründüğüne dair bir fikir edinmek için hızlıca bir göz atalım.

Pipeline kaynak koduna göz atmayı kolaylaştırmak için pipeline'ın çıkarılmış kopyasına işaret eden bir sembolik bağlantı oluşturun:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Bu, `tree -L 2 pipelines/nf-core/demo` komutuyla kodu incelemenizi veya dosyaları doğrudan açmanızı sağlayan bir kısayol oluşturur.

#### 1.2.5. Kod Organizasyonuna Genel Bakış

`nf-core/demo` dizinini bulmak ve açmak için `tree` komutunu ya da dosya gezginini kullanabilirsiniz.

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

    7 directories, 12 files
    ```

Gördüğünüz gibi, orada pek çok şey var; ancak bunların büyük çoğunluğu hakkında endişelenmenize gerek yok.

Kısaca belirtmek gerekirse, en üst düzeyde lisanslama, katkı yönergeleri, atıf ve davranış kuralları gibi proje bilgilerini özetleyen yardımcı dosyaların yanı sıra özet bilgiler içeren bir README dosyası bulabilirsiniz.
Ayrıntılı pipeline dokümantasyonu `docs` dizininde yer almaktadır.
Bu içeriklerin tamamı, nf-core web sitesindeki sayfaları programatik olarak oluşturmak için kullanılır; bu sayede her zaman kodla güncel kalırlar.

Geri kalanı için üç işlevsel kod dosyası grubunu ayırt edebiliriz:

1. Pipeline kod bileşenleri (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Pipeline yapılandırması
3. Pipeline parametreleri / girdiler ve doğrulama

Bu kursun bu bölümünde pipeline kod bileşenlerini ele almayacağız; ancak nf-core pipeline'larının son kullanıcısı olarak sizin için muhtemelen ilgili olacak yapılandırma ve doğrulama öğelerine değineceğiz.

!!! tip "İpucu"

    Herhangi bir nf-core pipeline'ının kaynak koduna GitHub üzerinden de göz atabilirsiniz; örneğin [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Her nf-core pipeline'ı aynı dizin düzenini izler; bu nedenle yapıyı bir kez öğrendikten sonra herhangi bir pipeline için yapılandırma dosyalarını, modülleri ve iş akışlarını aynı şekilde bulabilirsiniz.

Şimdi pipeline'ı çalıştırmaya geçelim!

### Özetle

nf-core web sitesi aracılığıyla bir pipeline'ı nasıl bulacağınızı ve kaynak kodun yerel bir kopyasını nasıl indireceğinizi artık biliyorsunuz.

### Sırada ne var?

Minimum çabayla bir nf-core pipeline'ını nasıl deneyeceğinizi öğrenin.

---

## 2. Pipeline'ı Test Profiliyle Deneme

Her nf-core pipeline'ı, kullanışlı bir şekilde bir test profiliyle birlikte gelir.
Bu, [nf-core/test-datasets](https://github.com/nf-core/test-datasets) deposunda barındırılan küçük bir test veri seti kullanılarak pipeline'ın çalışması için gereken minimum yapılandırma ayarları kümesidir.
Bir pipeline'ı küçük ölçekte hızlıca denemek için harika bir yöntemdir.

!!! tip "İpucu"

    Nextflow'un yapılandırma profili sistemi, farklı konteyner motorları veya yürütme ortamları arasında kolayca geçiş yapmanızı sağlar.
    Daha fazla ayrıntı için bkz. [Hello Nextflow Bölüm 6: Yapılandırma](../hello_nextflow/06_hello_config.md).

### 2.1. Test Profilini İnceleme

Bir pipeline'ın test profilinin ne belirttiğini çalıştırmadan önce kontrol etmek iyi bir uygulamadır.
`nf-core/demo` için `test` profili, `conf/test.config` yapılandırma dosyasında yer almaktadır.
Bunu, `nextflow pull`'un indirdiği pipeline kaynağının içinde, 1.2.4. bölümünde oluşturulan `pipelines` sembolik bağlantısı aracılığıyla yerel olarak bulabilirsiniz:

```bash
code pipelines/nf-core/demo/conf/test.config
```

Bu dosyanın içeriği şöyledir:

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
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Girdi verisi
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Üstteki yorum bloğunun, bu test profiliyle pipeline'ın nasıl çalıştırılacağını gösteren bir kullanım örneği içerdiğini hemen fark edeceksiniz.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

Sağlamamız gereken tek şey, örnek komuttaki köşeli parantezler arasında gösterilenlerdir: `<docker/singularity>` ve `<OUTDIR>`.

Hatırlatmak gerekirse, `<docker/singularity>` konteyner sisteminin seçimini ifade eder. Tüm nf-core pipeline'ları, tekrarlanabilirliği sağlamak ve yazılım kurulum sorunlarını ortadan kaldırmak için konteynerlerle (Docker, Singularity vb.) kullanılabilecek şekilde tasarlanmıştır.
Bu nedenle pipeline'ı test etmek için Docker mu yoksa Singularity mi kullanmak istediğimizi belirtmemiz gerekecek.

`--outdir <OUTDIR>` kısmı ise Nextflow'un pipeline çıktılarını yazacağı dizini ifade eder.
Bunun için bir isim belirlememiz gerekiyor; istediğimiz bir isim verebiliriz.
Henüz mevcut değilse Nextflow çalışma zamanında bizim için oluşturacaktır.

Yorum bloğundan sonraki bölüme geçersek, test profili test için önceden yapılandırılmış olanları gösterir: en önemlisi, `input` parametresi zaten bir test veri setine işaret edecek şekilde ayarlanmıştır; dolayısıyla kendi verilerimizi sağlamamıza gerek yoktur.
Önceden yapılandırılmış girdinin bağlantısını takip ederseniz, bunun çeşitli deneysel örnekler için örnek tanımlayıcıları ve dosya yolları içeren bir CSV dosyası olduğunu göreceksiniz.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Buna örnek sayfası (samplesheet) denir ve nf-core pipeline'larına en yaygın girdi biçimidir.
Veri formatlarına ve türlerine aşina değilseniz endişelenmeyin; bundan sonrası için önemli değil.

Artık pipeline'ı denemek için ihtiyacımız olan her şeye sahibiz.

### 2.2. Pipeline'ı Çalıştırma

Yukarıda belirtildiği gibi, örnek test komutunu neredeyse olduğu gibi kullanabiliriz; yalnızca hangi yazılım paketlemenin kullanılacağını ve çıktı dizinine ne ad verileceğini belirtmemiz gerekiyor.
Burada konteyner sistemi olarak Docker'ı ve çıktı dizini adı olarak `demo-results`'ı kullanacağız.

Bununla birlikte test komutunu çalıştırabiliriz:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Çıktınız bununla eşleşiyorsa tebrikler! İlk nf-core pipeline'ınızı çalıştırdınız.

Temel bir Nextflow pipeline'ı çalıştırdığınızdakinden çok daha fazla konsol çıktısı olduğunu fark edeceksiniz.
Pipeline'ın sürümünü, girdileri ve çıktıları ile bazı yapılandırma öğelerini özetleyen bir başlık bulunmaktadır.

!!! info "Bilgi"

    Çıktınız farklı zaman damgaları, yürütme adları ve dosya yolları gösterecektir; ancak genel yapı ve süreç yürütmesi benzer olmalıdır.

Çıktının üst kısmına yakın şu satıra dikkat edin:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Bu, pipeline'ın hangi revizyonunun kullanıldığını gösterir.
Bir sürüm belirtmediğimiz için Nextflow, `master` üzerindeki en son commit'i kullandı.
Tekrarlanabilir çalıştırmalar için `-r` bayrağını kullanarak belirli bir sürümü sabitlemelisiniz:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

Bu, yeni commit'ler veya sürümlerden bağımsız olarak her seferinde aynı pipeline kodunun kullanılmasını sağlar.
Bu eğitimde basitlik adına `-r` bayrağını atlıyoruz; ancak üretim ortamında her zaman belirtmelisiniz.

Yürütme çıktısına geçersek, hangi süreçlerin çalıştırıldığını bize söyleyen satırlara bakalım:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Bu bize, nf-core web sitesindeki pipeline dokümantasyon sayfasında gösterilen dört araca karşılık gelen dört sürecin çalıştırıldığını söyler: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` ve `COWPY`.

Burada gösterilen tam süreç adları, örneğin `NFCORE_DEMO:DEMO:MULTIQC`, giriş niteliğindeki Hello Nextflow materyalinde görmüş olabileceğinizden daha uzundur.
Bunlar, üst iş akışlarının adlarını içerir ve pipeline kodunun modülerliğini yansıtır.
nf-core tarzı pipeline'lar geliştirmeyi kendiniz öğrenmek istiyorsanız [Build with nf-core](../nfcore_build/index.md) kursuna bakın.

### 2.3. Pipeline Çıktılarını İnceleme

Son olarak, pipeline tarafından üretilen `demo-results` dizinine bir göz atalım.

```bash
tree -L 2 demo-results
```

??? abstract "Dizin içeriği"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
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
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Bu çok fazla gibi görünebilir.
`nf-core/demo` pipeline'ının çıktıları hakkında daha fazla bilgi edinmek için [dokümantasyon sayfasına](https://nf-co.re/demo/1.2.0/docs/output/) bakın.

Bu aşamada dikkat edilmesi gereken önemli nokta, sonuçların modüle göre düzenlendiği ve ayrıca pipeline yürütmesiyle ilgili çeşitli zaman damgalı raporlar içeren `pipeline_info` adlı bir dizinin bulunduğudur.

Örneğin, `execution_timeline_*` dosyası hangi süreçlerin çalıştırıldığını, hangi sırayla ve ne kadar sürdüğünü gösterir:

![yürütme zaman çizelgesi raporu](./img/execution_timeline.png)

!!! info "Bilgi"

    Burada görevler paralel olarak çalıştırılmadı; çünkü Github Codespaces'te minimalist bir makinede çalışıyoruz.
    Bu görevlerin paralel çalıştığını görmek için codespace'inizin CPU tahsisini ve test yapılandırmasındaki kaynak sınırlarını artırmayı deneyin.

Bu raporlar tüm nf-core pipeline'ları için otomatik olarak oluşturulur.

### Özetle

Bir nf-core pipeline'ını yerleşik test profiliyle nasıl çalıştıracağınızı ve çıktılarını nerede bulacağınızı artık biliyorsunuz.

### Sırada ne var?

Pipeline yürütmesini nasıl yapılandıracağınızı öğreneceğiniz [Bölüm 2](./02_configure_execution.md)'ye geçin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Bir nf-core pipeline'ını bulma ve indirme ile kod yapısını inceleme
- Yerleşik test profiliyle bir pipeline'ı çalıştırma
