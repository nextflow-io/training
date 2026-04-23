# Dosya Girdi İşleme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bilimsel analiz iş akışları çoğunlukla çok sayıda dosyanın işlenmesini gerektirir.
Nextflow, dosyaları verimli şekilde yönetmek için güçlü araçlar sunar; bu araçlar sayesinde verilerinizi minimum kodla düzenleyip işleyebilirsiniz.

### Öğrenme hedefleri

Bu yan görevde, Nextflow'un dosyaları nasıl ele aldığını keşfedeceğiz: temel dosya işlemlerinden dosya koleksiyonlarıyla çalışmaya yönelik daha gelişmiş tekniklere kadar.
Dosya adlarından üst veri (metadata) çıkarmayı öğreneceksiniz; bu, bilimsel analiz pipeline'larında sıkça karşılaşılan bir gerekliliktir.

Bu yan görevin sonunda şunları yapabileceksiniz:

- Nextflow'un `file()` metodunu kullanarak dosya yolu dizelerinden Path nesneleri oluşturmak
- Ad, uzantı ve üst dizin gibi dosya özelliklerine erişmek
- URI'lar aracılığıyla hem yerel hem de uzak dosyaları şeffaf biçimde yönetmek
- `channel.fromPath()` ve `channel.fromFilePairs()` ile dosya yönetimini otomatikleştirmek için kanallar kullanmak
- Dize işleme yöntemleriyle dosya adlarından üst veri çıkarmak ve yapılandırmak
- Desen eşleştirme ve glob ifadeleri kullanarak ilgili dosyaları gruplamak
- Dosya işlemlerini uygun girdi yönetimiyle Nextflow süreçlerine entegre etmek
- Üst veriye dayalı dizin yapıları kullanarak süreç çıktılarını düzenlemek

Bu beceriler, farklı türde dosya girdilerini büyük bir esneklikle işleyebilen iş akışları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmış olmanız gerekir:

- [Hello Nextflow](../../hello_nextflow/) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmaları (süreçler, kanallar, operatörler) konusunda rahat olmak.

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../envsetup/index.md) sayfasında açıklandığı şekilde açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Şimdi bu eğitime ait dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/working_with_files
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

`main.nf` adlı basit bir iş akışı dosyası, iki modül dosyası içeren bir `modules` dizini ve bazı örnek veri dosyaları içeren bir `data` dizini bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Bu dizin, üç hastaya (A, B, C) ait çift uçlu dizileme verilerini içermektedir.

Her hasta için `tumor` (genellikle tümör biyopsisinden elde edilen) veya `normal` (sağlıklı doku ya da kandan alınan) türünde örnekler mevcuttur.
Kanser analiziyle tanışık değilseniz, bunun karşılaştırmalı analizler yapmak amacıyla eşleştirilmiş tümör/normal örnekler kullanan deneysel bir modele karşılık geldiğini bilmeniz yeterlidir.

Özellikle A hastası için iki teknik replikat (tekrar) setimiz bulunmaktadır.

Dizileme verisi dosyaları, 'ileri okumalar' ve 'ters okumalar' olarak bilinen okumalar için tipik `_R1_` ve `_R2_` kuralına göre adlandırılmıştır.

_Bu deneysel tasarıma aşina değilseniz endişelenmeyin; bu eğitimi anlamak için kritik değildir._

#### Görevi inceleyin

Göreviniz, aşağıdakileri yapacak bir Nextflow iş akışı yazmaktır:

1. Nextflow'un dosya yönetimi metodlarını kullanarak girdi dosyalarını **yüklemek**
2. Dosya adı yapısından üst veriyi (hasta kimliği, replikat, örnek türü) **çıkarmak**
3. `channel.fromFilePairs()` kullanarak eşleştirilmiş dosyaları (R1/R2) **gruplamak**
4. Sağlanan analiz modülüyle dosyaları **işlemek**
5. Çıktıları, çıkarılan üst veriye dayalı bir dizin yapısına göre **düzenlemek**

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışır durumda
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebildiyseniz, başlayabilirsiniz.

---

## 1. Temel dosya işlemleri

### 1.1. `.class` ile bir nesnenin türünü belirleme

`main.nf` iş akışı dosyasına bir göz atın:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Bir dize yolundan Path nesnesi oluştur
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Bu, iş akışında tek bir dosya yoluna başvuran (herhangi bir süreç içermeyen) küçük bir iş akışıdır; dosyayı sınıfıyla birlikte konsola yazdırır.

??? info "`.class` nedir?"

    Nextflow'da `.class`, üzerinde çalıştığımız nesnenin türünü bize söyler. "Bu ne tür bir şey?" diye sormak gibidir; bir dize mi, sayı mı, dosya mı yoksa başka bir şey mi olduğunu öğrenmemizi sağlar.
    Bu, sonraki bölümlerde düz bir dize ile Path nesnesi arasındaki farkı göstermemize yardımcı olacaktır.

İş akışını çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Gördüğünüz gibi, Nextflow dize yolunu tam olarak yazdığımız şekilde yazdırdı.

Bu yalnızca metin çıktısıdır; Nextflow henüz bu konuda özel bir işlem yapmamıştır.
Ayrıca Nextflow açısından bunun yalnızca bir dize (`java.lang.String` sınıfından) olduğunu doğruladık.
Nextflow'a bunun bir dosyaya karşılık geldiğini henüz söylemediğimiz için bu mantıklıdır.

### 1.2. `file()` ile Path nesnesi oluşturma

Nextflow'a dosyaları nasıl yöneteceğini, yol dizelerinden [Path nesneleri](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) oluşturarak söyleyebiliriz.

İş akışımızda, `data/patientA_rep1_normal_R1_001.fastq.gz` dize yolunu, dosya özelliklerine ve işlemlerine erişim sağlayan `file()` metodunu kullanarak bir Path nesnesine dönüştürebiliriz.

`main.nf` dosyasını düzenleyerek dizeyi aşağıdaki gibi `file()` ile sarın:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Şimdi iş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Bu sefer, girdi olarak verdiğimiz göreli yol yerine tam mutlak yolu görüyorsunuz.

Nextflow, dizimizi bir Path nesnesine dönüştürerek sistemdeki gerçek dosya konumuna çözdü.
Dosya yolu artık mutlak olacak; örneğin `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz` gibi.

Path nesnesinin sınıfının `sun.nio.fs.UnixPath` olduğuna da dikkat edin: bu, Nextflow'un yerel dosyaları temsil etme biçimidir.
İleride göreceğimiz gibi, uzak dosyalar farklı sınıf adlarına sahip olacaktır (örneğin HTTP dosyaları için `nextflow.file.http.XPath`), ancak hepsi tamamen aynı şekilde çalışır ve iş akışlarınızda aynı biçimde kullanılabilir.

!!! tip "İpucu"

    **Temel fark:**

    - **Path dizesi**: Nextflow'un karakter dizisi olarak ele aldığı yalnızca metin
    - **Path nesnesi**: Nextflow'un üzerinde işlem yapabildiği akıllı bir dosya referansı

    Şöyle düşünebilirsiniz: bir yol dizesi kağıda yazılmış bir adres gibidir; Path nesnesi ise o adresin bir GPS cihazına yüklenmiş hali gibidir; cihaz oraya nasıl gidileceğini bilir ve yolculuk hakkında ayrıntılar sunabilir.

### 1.3. Dosya özelliklerine erişme

Bu neden yararlıdır? Artık Nextflow, `myFile`'ın yalnızca bir dize değil, bir Path nesnesi olduğunu anladığına göre, Path nesnesinin çeşitli özelliklerine erişebiliriz.

İş akışımızı yerleşik dosya özelliklerini yazdıracak şekilde güncelleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Yukarıda konsola yazdırılan çeşitli dosya özelliklerini görebilirsiniz.

### 1.4. Dosyayı bir sürece beslemek

Dizeler ile Path nesneleri arasındaki fark, süreçlerle gerçek iş akışları oluşturmaya başladığınızda kritik hale gelir.
Şimdiye kadar Nextflow'un girdi dosyamızı artık bir dosya olarak ele aldığını doğruladık; şimdi de bu dosya üzerinde bir süreçte gerçekten bir şeyler çalıştırıp çalıştıramayacağımıza bakalım.

#### 1.4.1. Süreci içe aktarın ve kodu inceleyin

Size, bir dosya girdisi alıp içindeki satır sayısını sayan `COUNT_LINES` adlı önceden yazılmış bir süreç modülü sunuyoruz.

Süreci iş akışında kullanmak için, iş akışı bloğundan önce bir include ifadesi eklemeniz yeterlidir:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Modül dosyasını açarak kodunu inceleyebilirsiniz:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Gördüğünüz gibi, dosyayı açıp içindeki satır sayısını sayan oldukça basit bir betiktir.

??? info "`debug true` ne işe yarar?"

    Süreç tanımındaki `debug true` yönergesi, Nextflow'un betiğinizden gelen çıktıyı (satır sayısı "40" gibi) doğrudan yürütme günlüğüne yazdırmasını sağlar.
    Bu olmadan yalnızca süreç yürütme durumunu görürsünüz; betiğinizden gelen gerçek çıktıyı göremezsiniz.

    Nextflow süreçlerinde hata ayıklama hakkında daha fazla bilgi için [Nextflow İş Akışlarında Hata Ayıklama](debugging.md) yan görevine bakın.

#### 1.4.2. `COUNT_LINES` çağrısı ekleme

Süreç artık iş akışı için kullanılabilir olduğuna göre, girdi dosyası üzerinde çalıştırmak için `COUNT_LINES` sürecine bir çağrı ekleyebiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Şimdi iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Bu, dosya üzerinde bir süreç içinde uygun şekilde işlem yapabildiğimizi göstermektedir.

Özellikle Nextflow şu işlemleri başarıyla gerçekleştirdi:

- Dosyayı çalışma dizinine taşıdı (staged)
- .gz dosyasını açtı
- Satırları saydı (bu durumda 40 satır)
- Hatasız tamamladı

Bu sorunsuz işlemin anahtarı, Nextflow'a girdimizin bir dosya olduğunu ve bu şekilde ele alınması gerektiğini açıkça söylememizdir.

### 1.5. Temel dosya girdi hatalarını gidermek

Bu durum Nextflow'a yeni başlayanların sık takıldığı bir noktadır; bu nedenle yanlış yapıldığında ne olduğuna birkaç dakika ayıralım.

Dosya yönetimini yanlış yapabileceğiniz iki ana yer vardır: iş akışı düzeyi ve süreç düzeyi.

#### 1.5.1. İş akışı düzeyinde hata

İş akışı bloğunda girdiyi belirtirken dosyayı yeniden bir dize olarak ele alırsak ne olduğunu görelim.

İş akışında aşağıdaki düzenlemeleri yapın; yola özgü yazdırma ifadelerini yorum satırına aldığınızdan emin olun:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

Şimdi iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Önemli olan kısım şudur:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Bir `path` girdisi belirttiğinizde, Nextflow gerçek dosya referansları ilettiğinizi doğrular; yalnızca dizeler değil.
Bu hata, `'data/patientA_rep1_normal_R1_001.fastq.gz'`'nin bir dize olduğunu, Path nesnesi olmadığını söylemektedir; dolayısıyla geçerli bir yol değeri değildir.

Nextflow sorunu hemen tespit etti ve süreci başlatmadan önce durdu.

#### 1.5.2. Süreç düzeyinde hata

Nextflow'a girdiyi dosya olarak ele almasını söylemeyi unutabileceğimiz diğer yer, süreç tanımıdır.

!!! warning "1.5.1'deki iş akışı hatasını koruyun"

    Bu testin doğru çalışması için iş akışını bozuk halde bırakın (`file()` yerine düz dize kullanarak).
    Süreçteki `val` ile birleştirildiğinde, aşağıda gösterilen hata üretilir.

Modülde aşağıdaki düzenlemeyi yapın:

=== "Sonra"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Önce"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

Şimdi iş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Süreç yukarıda belirtildiği gibi hata ayıklama bilgisi çıktısı verecek şekilde ayarlandığından, hata hakkında çok sayıda ayrıntı gösterilmektedir.

En ilgili bölümler şunlardır:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Bu, sistemin dosyayı bulamadığını söylemektedir; ancak yola bakarsanız, o konumda o adda bir dosya mevcuttur.

Bu işlemi çalıştırdığımızda Nextflow, dize değerini betiğe iletti; ancak gerçek dosyayı çalışma dizinine _taşımadı_ (stage etmedi).
Dolayısıyla süreç, göreli dizeyi (`data/patientA_rep1_normal_R1_001.fastq.gz`) kullanmaya çalıştı; ancak bu dosya sürecin çalışma dizininde mevcut değildir.

Bu iki örnek birlikte ele alındığında, bir girdinin dosya olarak yönetilmesi gerektiğini Nextflow'a bildirmenin ne kadar önemli olduğunu göstermektedir.

!!! note "Not"

    Bir sonraki bölüme geçmeden önce her iki kasıtlı hatayı da düzelttiğinizden emin olun.

### Özetle

- Path dizeleri ve Path nesneleri: Dizeler yalnızca metindir; Path nesneleri ise akıllı dosya referanslarıdır.
- `file()` metodu, bir dize yolunu Nextflow'un üzerinde çalışabileceği bir Path nesnesine dönüştürür.
- `name`, `simpleName`, `extension` ve `parent` gibi dosya özelliklerine [dosya özelliklerini kullanarak](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) erişebilirsiniz.
- Dizeler yerine Path nesneleri kullanmak, Nextflow'un iş akışınızdaki dosyaları düzgün şekilde yönetmesini sağlar.
- Süreç Girdi Sonuçları: Dosyaların süreçler tarafından doğru şekilde taşınması (staged) ve kullanılabilir olması için uygun dosya yönetimi, dizeler değil Path nesneleri gerektirir.

---

## 2. Uzak dosyaları kullanmak

Nextflow'un temel özelliklerinden biri, yerel dosyalar (aynı makinedeki) ile internet üzerinden erişilebilen uzak dosyalar arasında sorunsuz geçiş yapabilmektir.

Doğru yapıyorsanız, farklı konumlardan gelen dosyalara uyum sağlamak için iş akışınızın mantığını hiçbir zaman değiştirmeniz gerekmemelidir.
Uzak bir dosya kullanmak için tek yapmanız gereken, iş akışına sağlarken dosya yolunda uygun öneki belirtmektir.

Örneğin, `/path/to/data` öneksizdir ve 'normal' bir yerel dosya yolunu gösterirken, `s3://path/to/data` ifadesi `s3://` önekini içerir ve Amazon'un S3 nesne depolama alanında bulunduğunu gösterir.

Pek çok farklı protokol desteklenmektedir:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Bunlardan herhangi birini kullanmak için, dizedeki ilgili öneki belirtmeniz yeterlidir; bu durumda teknik olarak dosya yolu yerine Tekdüzen Kaynak Tanımlayıcısı (URI) adı verilir.
Nextflow, kimlik doğrulamayı ve dosyaları doğru yere taşımayı, indirme veya yüklemeyi ve beklediğiniz diğer tüm dosya işlemlerini otomatik olarak gerçekleştirir.

Bu sistemin temel gücü, herhangi bir pipeline mantığını değiştirmeden ortamlar arasında geçiş yapabilmemizi sağlamasıdır.
Örneğin, yalnızca URI'yi değiştirerek küçük bir yerel test setiyle geliştirme yapabilir, ardından uzak depolama alanında bulunan tam ölçekli bir test setine geçebilirsiniz.

### 2.1. İnternetten bir dosya kullanmak

Bunu test etmek için, iş akışımıza sağladığımız yerel yolu, Github'da depolanan aynı verinin bir kopyasına işaret eden bir HTTPS yoluyla değiştirelim.

!!! warning "Uyarı"

    Bu yalnızca aktif bir internet bağlantınız varsa çalışacaktır.

`main.nf` dosyasını tekrar açın ve girdi yolunu aşağıdaki gibi değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // İnternetten uzak bir dosya kullanma
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

İş akışını çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Çalışıyor! Çok az şeyin değiştiğini görebilirsiniz.

Konsol çıktısındaki tek fark, yol nesnesi sınıfının artık `nextflow.file.http.XPath` olmasıdır; yerel yol için sınıf `sun.nio.fs.UnixPath`'ti.
Bu sınıfları hatırlamanıza gerek yoktur; bunu yalnızca Nextflow'un farklı konumları uygun şekilde tanımlayıp yönettiğini göstermek için belirtiyoruz.

Arka planda Nextflow, dosyayı çalışma dizini içindeki bir hazırlama dizinine indirdi.
Bu hazırlanmış dosya daha sonra yerel bir dosya olarak ele alınabilir ve ilgili süreç dizinine sembolik bağlantıyla bağlanabilir.

Bunu, sürecin karma değerinde bulunan çalışma dizininin içeriğine bakarak doğrulayabilirsiniz.

??? abstract "Çalışma dizini içeriği"

    Süreç karması `8a/2ab7ca` ise çalışma dizinini şu şekilde inceleyebilirsiniz:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Sembolik bağlantı, Nextflow'un otomatik olarak indirdiği uzak dosyanın hazırlanmış bir kopyasına işaret etmektedir.

Daha büyük dosyalar için indirme adımı, yerel dosyalar üzerinde çalışmaya kıyasla ek süre alacaktır.
Ancak Nextflow, gereksiz indirmeleri önlemek için hazırlanmış bir kopyasının zaten mevcut olup olmadığını kontrol eder.
Dolayısıyla aynı dosya üzerinde tekrar çalıştırırsanız ve hazırlanmış dosyayı silmediyseniz, Nextflow hazırlanmış kopyayı kullanacaktır.

Bu, Nextflow'un temel özelliklerinden biri olan yerel ve uzak veriler arasında geçişin ne kadar kolay olduğunu göstermektedir.

!!! note "Not"

    Bu ilkenin önemli bir istisnası, HTTPS ile glob desenleri veya dizin yolları kullanamayacağınızdır; çünkü HTTPS birden fazla dosyayı listeleyemez; bu nedenle tam dosya URL'lerini belirtmeniz gerekir.
    Ancak blob depolama (`s3://`, `az://`, `gs://`) gibi diğer depolama protokolleri hem glob hem de dizin yollarını kullanabilir.

    Bulut depolama ile glob desenlerini şu şekilde kullanabilirsiniz:

    ```groovy title="Cloud storage examples (not runnable in this environment)"
    // Glob desenleriyle S3 - birden fazla dosyayla eşleşir
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Glob desenleriyle Azure Blob Storage
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Glob desenleriyle Google Cloud Storage
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Bir sonraki bölümde glob'larla pratikte nasıl çalışılacağını göstereceğiz.

### 2.2. Yerel dosyaya geri dönmek

Bu yan görevin geri kalanında yerel örnek dosyalarımızı kullanmaya devam edeceğiz; bu nedenle iş akışı girdisini orijinal dosyaya geri döndürelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Özetle

- Uzak verilere URI kullanılarak erişilir (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow, bu yollar süreçlere beslendiği sürece verileri otomatik olarak indirir ve doğru yere taşır
- Uzak dosyaları indirmek veya yüklemek için mantık yazmayın!
- Yerel ve uzak dosyalar farklı nesne türleri üretir ancak aynı şekilde çalışır
- **Önemli**: HTTP/HTTPS yalnızca tek dosyalarla çalışır (glob desenleri desteklenmez)
- Bulut depolama (S3, Azure, GCS) hem tek dosyaları hem de glob desenlerini destekler
- Kod mantığını değiştirmeden yerel ve uzak veri kaynakları arasında sorunsuz geçiş yapabilirsiniz (protokol gerekli işlemleri desteklediği sürece)

---

## 3. `fromPath()` kanal fabrikasını kullanmak

Şimdiye kadar tek bir dosyayla çalıştık; ancak Nextflow'da genellikle işlemek için birden fazla girdi dosyası içeren bir girdi kanalı oluşturmak isteyeceğiz.

Bunu yapmanın naif bir yolu, `file()` metodunu [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) ile şu şekilde birleştirmek olurdu:

```groovy title="Syntax example"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Bu çalışır; ancak hantaldır.

!!! tip "`file()` ile `channel.fromPath()` ne zaman kullanılır?"

    - Doğrudan işleme için tek bir Path nesnesine ihtiyaç duyduğunuzda `file()` kullanın (bir dosyanın var olup olmadığını kontrol etmek, özelliklerine erişmek veya tek bir süreç çağrısına iletmek için)
    - Özellikle glob desenleriyle birden fazla dosya tutabilen bir kanala ihtiyaç duyduğunuzda ya da dosyalar birden fazla süreçten geçecekse `channel.fromPath()` kullanın

İşte burada [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) devreye girer: bir veya daha fazla statik dosya dizesinden ve glob desenlerinden kanal oluşturmak için ihtiyaç duyduğumuz tüm işlevselliği bir araya getiren kullanışlı bir kanal fabrikasıdır.

### 3.1. Kanal fabrikasını ekleme

İş akışımızı `channel.fromPath` kullanacak şekilde güncelleyelim.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Dosya özelliklerini yazdır
        /* Şimdilik bunları yorum satırına alalım, sonra geri döneceğiz!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Dosyadaki satırları say
        // COUNT_LINES(myFile)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

Özellikleri yazdıran kodu şimdilik yorum satırına aldık ve bunun yerine yalnızca dosya adını yazdırmak için bir `.view` ifadesi ekledik.

İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Gördüğünüz gibi, dosya yolu kanalda `Path` türünde bir nesne olarak yüklenmektedir.
Bu, `file()`'ın yapacağına benzer; ancak şimdi istediğimizde daha fazla dosya yükleyebileceğimiz bir kanalımız var.

`channel.fromPath()` kullanmak, bir dosya listesiyle doldurulmuş yeni bir kanal oluşturmanın kullanışlı bir yoludur.

### 3.2. Kanaldaki dosyaların özelliklerini görüntüleme

Kanal fabrikasını ilk kullandığımızda kodu basitleştirip yalnızca dosya adını yazdırdık.

Şimdi tam dosya özelliklerini yazdırmaya geri dönelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Dosyadaki satırları say
        COUNT_LINES(ch_files)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Dosyadaki satırları say
        // COUNT_LINES(ch_files)
    ```

Kanal tabanlı yaklaşımımızla dosya işlemenin hâlâ doğru çalıştığını doğrulamak için `COUNT_LINES` süreç çağrısını da yeniden etkinleştiriyoruz.

İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

İşte bu kadar; öncekiyle aynı sonuçlar, ancak artık dosya bir kanalda olduğu için daha fazlasını ekleyebiliriz.

### 3.3. Birden fazla dosyayı eşleştirmek için glob kullanmak

Kanala daha fazla dosya yüklemenin birkaç yolu vardır.
Burada size glob desenlerini nasıl kullanacağınızı göstereceğiz; bunlar, joker karakterlere dayalı olarak dosya ve dizin adlarını eşleştirip almak için kullanışlı bir yöntemdir.
Bu desenleri eşleştirme işlemine "globbing" veya "dosya adı genişletme" denir.

!!! note "Not"

    Daha önce belirtildiği gibi, Nextflow; HTTPS dosya yolları dışında çoğu durumda girdi ve çıktı dosyalarını yönetmek için globbing'i destekler; çünkü HTTPS birden fazla dosyayı listeleyemez.

Belirli bir hastayla (`patientA`) ilişkili bir dosya çiftindeki her iki dosyayı almak istediğimizi varsayalım:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Dosya adları arasındaki tek fark replikat numarası, yani `R`'den sonraki sayı olduğundan, joker karakter `*`'ı şu şekilde sayının yerine kullanabiliriz:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

İhtiyacımız olan glob deseni budur.

Şimdi tek yapmamız gereken, kanal fabrikasındaki dosya yolunu bu glob desenini kullanacak şekilde güncellemektir:

=== "Sonra"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7"
      // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow, bunun bir glob deseni olduğunu otomatik olarak tanıyacak ve uygun şekilde işleyecektir.

Bunu test etmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Gördüğünüz gibi, kanalımızda artık iki Path nesnesi var; bu, Nextflow'un dosya adı genişletmeyi doğru şekilde yaptığını ve her iki dosyayı da beklendiği gibi yükleyip işlediğini göstermektedir.

Bu yöntemi kullanarak, glob desenini değiştirerek istediğimiz kadar çok veya az dosya alabiliriz. Deseni daha geniş yaparsak, örneğin dosya adlarının değişken tüm bölümlerini `*` ile değiştirerek (_örn._ `data/patient*_rep*_*_R*_001.fastq.gz`) `data` dizinindeki tüm örnek dosyaları alabiliriz.

### Özetle

- `channel.fromPath()`, bir desenle eşleşen dosyalardan oluşan bir kanal oluşturur
- Her dosya kanalda ayrı bir öğe olarak yayınlanır
- Birden fazla dosyayı eşleştirmek için glob deseni kullanabiliriz
- Dosyalar otomatik olarak tam özelliklerle Path nesnelerine dönüştürülür
- `.view()` metodu kanal içeriğinin incelenmesine olanak tanır

---

## 4. Dosya adlarından temel üst veri çıkarmak

Çoğu bilimsel alanda, veriyi içeren dosyaların adlarına üst verinin kodlanması çok yaygındır.
Örneğin biyoinformatikte, dizileme verisi içeren dosyalar genellikle örnek, koşul, replikat ve okuma numarası hakkında bilgi kodlayacak şekilde adlandırılır.

Dosya adları tutarlı bir kurala göre oluşturulmuşsa, bu üst veriyi standart bir şekilde çıkarabilir ve analiziniz boyunca kullanabilirsiniz.
Elbette bu büyük bir 'eğer'dir ve dosya adı yapısına güvendiğinizde çok dikkatli olmalısınız; ancak gerçek şu ki bu yaklaşım çok yaygın kullanılmaktadır; bu nedenle Nextflow'da nasıl yapıldığına bakalım.

Örnek verilerimizde, dosya adlarının tutarlı biçimde yapılandırılmış üst veri içerdiğini biliyoruz.
Örneğin, `patientA_rep1_normal_R2_001` dosya adı şunları kodlamaktadır:

- hasta kimliği: `patientA`
- replikat kimliği: `rep1`
- örnek türü: `normal` (`tumor`'ın karşıtı)
- okuma seti: `R1` (`R2`'nin karşıtı)

İş akışımızı bu bilgileri üç adımda almak üzere değiştireceğiz:

1. Üst veriyi içeren dosyanın `simpleName`'ini almak
2. `tokenize()` adlı bir metot kullanarak üst veriyi ayırmak
3. Üst veriyi düzenlemek için bir map kullanmak

!!! warning "Uyarı"

    Hasta adları veya diğer tanımlayıcı özellikler gibi hassas bilgileri dosya adlarına asla kodlamamalısınız; bu, hasta gizliliğini veya diğer ilgili güvenlik kısıtlamalarını tehlikeye atabilir.

### 4.1. `simpleName`'i almak

`simpleName`, yolu ve uzantısı çıkarılmış dosya adına karşılık gelen bir dosya özelliğidir.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Bu, `simpleName`'i alır ve bir `map()` işlemi kullanarak tam dosya nesnesiyle ilişkilendirir.

Çalıştığını test etmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Kanaldaki her öğe artık `simpleName`'i ve orijinal dosya nesnesini içeren bir demettir.

### 4.2. `simpleName`'den üst veriyi çıkarmak

Bu noktada istediğimiz üst veri `simpleName`'e gömülüdür; ancak öğelere doğrudan erişemeyiz.
Bu nedenle `simpleName`'i bileşenlerine ayırmamız gerekir.
Neyse ki bu bileşenler orijinal dosya adında yalnızca alt çizgilerle ayrılmıştır; bu nedenle bu görev için mükemmel olan `tokenize()` adlı yaygın bir Nextflow metodunu uygulayabiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

`tokenize()` metodu, `simpleName` dizesini alt çizgi bulduğu her yerde böler ve alt dizeleri içeren bir liste döndürür.

İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Artık kanalımızdaki her öğenin demeti, üst veri listesini (_örn._ `[patientA, rep1, normal, R1, 001]`) ve orijinal dosya nesnesini içermektedir.

Harika!
Hasta bilgilerini tek bir dizeden bir dize listesine ayırdık.
Artık hasta bilgisinin her bölümünü ayrı ayrı ele alabiliriz.

### 4.3. Üst veriyi düzenlemek için map kullanmak

Üst verimiz şu an yalnızca düz bir listedir.
Kullanımı yeterince kolaydır; ancak okunması güçtür.

```console
[patientA, rep1, normal, R1, 001]
```

3. indeksteki öğe nedir? Üst veri yapısının orijinal açıklamasına başvurmadan söyleyebilir misiniz?

Bu, her öğenin bir anahtar kümesine ve bunlarla ilişkili değerlere sahip olduğu, karşılık gelen değeri almak için her anahtara kolayca başvurabildiğiniz bir anahtar-değer deposu kullanmak için harika bir fırsattır.

Örneğimizde bu, şu düzenden:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Şu düzene geçmek anlamına gelir:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Nextflow'da buna [map](https://nextflow.io/docs/latest/script.html#maps) denir.

Şimdi düz listemizi bir map'e dönüştürelim.
İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Buradaki temel değişiklikler şunlardır:

- **Yıkımlama ataması**: `def (patient, replicate, type, readNum) = ...` ifadesi, tokenize edilmiş değerleri tek satırda adlandırılmış değişkenlere çıkarır
- **Map değişmez sözdizimi**: `[id: patient, replicate: ...]` ifadesi, her anahtarın (`id` gibi) bir değerle (`patient` gibi) ilişkilendirildiği bir map oluşturur
- **İç içe yapı**: Dıştaki liste `[..., myFile]`, üst veri map'ini orijinal dosya nesnesiyle eşleştirir

Ayrıca gereksiz bazı karakterleri kaldırmak için `replace()` adlı bir dize değiştirme metodunu kullanarak üst veri dizelerinin bir kısmını basitleştirdik (_örn._ replikat kimliklerinden yalnızca sayıyı korumak için `replicate.replace('rep', '')`).

İş akışını tekrar çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Artık üst veri düzgün biçimde etiketlenmiştir (_örn._ `[id:patientA, replicate:1, type:normal, readNum:2]`); bu sayede neyin ne olduğunu anlamak çok daha kolaydır.

Üst verinin öğelerini iş akışında kullanmak da çok daha kolay olacak; kodumuz daha okunabilir ve bakımı daha kolay hale gelecektir.

### Özetle

- Nextflow'da dosya adlarını tam bir programlama dilinin gücüyle işleyebiliriz
- İlgili bilgileri çıkarmak için dosya adlarını dize olarak ele alabiliriz
- `tokenize()` ve `replace()` gibi metodların kullanımı, dosya adındaki dizeleri işlememize olanak tanır
- `.map()` işlemi, yapıyı koruyarak kanal öğelerini dönüştürür
- Yapılandırılmış üst veri (map'ler), kodu konumsal listelerden daha okunabilir ve bakımı daha kolay hale getirir

Sırada eşleştirilmiş veri dosyalarının nasıl yönetileceğine bakacağız.

---

## 5. Eşleştirilmiş veri dosyalarını yönetmek

Pek çok deneysel tasarım, açıkça eşleştirilmiş bir şekilde yönetilmesinden fayda sağlayan eşleştirilmiş veri dosyaları üretir.
Örneğin biyoinformatikte, dizileme verisi genellikle eşleştirilmiş okumalar biçiminde üretilir; yani aynı DNA parçasından kaynaklanan dizi dizeleri (genellikle zıt uçlardan okunduğu için 'ileri' ve 'ters' olarak adlandırılır).

Örnek verilerimizde de durum böyledir; R1 ve R2, iki okuma setini ifade eder.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow, bu tür eşleştirilmiş dosyalarla çalışmak için `channel.fromFilePairs()` adlı özel bir kanal fabrikası sunar; bu fabrika, paylaşılan bir adlandırma desenine göre dosyaları otomatik olarak gruplar. Bu sayede eşleştirilmiş dosyaları daha az çabayla daha sıkı bir şekilde ilişkilendirebilirsiniz.

İş akışımızı bundan yararlanacak şekilde değiştireceğiz.
Bu iki adım gerektirecektir:

1. Kanal fabrikasını `channel.fromFilePairs()` olarak değiştirmek
2. Üst veriyi çıkarmak ve eşleştirmek

### 5.1. Kanal fabrikasını `channel.fromFilePairs()` olarak değiştirmek

`channel.fromFilePairs` kullanmak için, Nextflow'un bir çiftteki iki üyeyi tanımlamak amacıyla kullanması gereken deseni belirtmemiz gerekir.

Örnek verilerimize geri dönerek, adlandırma desenini şu şekilde formalize edebiliriz:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Bu, daha önce kullandığımız glob desenine benzer; ancak bu, çiftin iki üyesini tanımlayan alt dizeleri (R'den hemen sonra gelen `1` veya `2`) özellikle sıralar.

`main.nf` iş akışını buna göre güncelleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Şimdilik eşleştirmeyi yorum satırına alalım, sonra geri döneceğiz!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Kanal fabrikasını değiştirdik ve dosya eşleştirme desenini uyarladık; bu arada map işlemini de yorum satırına aldık.
Bunu daha sonra birkaç değişiklikle geri ekleyeceğiz.

Test etmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="7-8"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_koch] DSL2 - revision: 44fdf66105

    [-        ] COUNT_LINES -
    [-        ] COUNT_LINES -
    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ERROR ~ Error executing process > 'COUNT_LINES (1)'

    Caused by:
      Not a valid path value: 'patientA_rep1_normal_R'



    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Eyvah, bu sefer çalıştırma başarısız oldu!

Hata mesajının ilgili kısmı şurasıdır:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Bunun nedeni kanal fabrikasını değiştirmiş olmamızdır.
Şimdiye kadar orijinal girdi kanalı yalnızca dosya yollarını içeriyordu.
Yaptığımız tüm üst veri işlemleri aslında kanal içeriğini etkilemiyordu.

Artık `.fromFilePairs` kanal fabrikasını kullandığımıza göre, elde edilen kanalın içeriği farklıdır.
Yalnızca bir kanal öğesi görüyoruz; bu öğe iki öğeden oluşan bir demet içermektedir: iki dosya tarafından paylaşılan `simpleName` bölümü (tanımlayıcı olarak kullanılır) ve `id, [ file1, file2 ]` biçiminde iki dosya nesnesini içeren bir demet.

Bu harika; çünkü Nextflow, paylaşılan öneki inceleyerek hasta adını çıkarma ve bunu hasta tanımlayıcısı olarak kullanma gibi zor işi yapmıştır.

Ancak bu, mevcut iş akışımızı bozmaktadır.
Süreci değiştirmeden `COUNT_LINES`'ı aynı şekilde çalıştırmak isteseydik, dosya yollarını çıkarmak için bir eşleştirme işlemi uygulamamız gerekirdi.
Ancak bunu yapmayacağız; çünkü nihai hedefimiz, dosya çiftlerini uygun şekilde işleyen farklı bir süreç olan `ANALYZE_READS`'i kullanmaktır.

Bu nedenle `COUNT_LINES` çağrısını yorum satırına alalım (veya silelim) ve devam edelim.

=== "Sonra"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Dosyadaki satırları say
        // COUNT_LINES(ch_files)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="26" hl_lines="2"
        // Dosyadaki satırları say
        COUNT_LINES(ch_files)
    ```

`COUNT_LINES` include ifadesini de yorum satırına alabilir veya silebilirsiniz; ancak bunun işlevsel bir etkisi olmayacaktır.

Şimdi iş akışını tekrar çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="5"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [fabulous_davinci] DSL2 - revision: 22b53268dc

    [patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Yaşasın, bu sefer iş akışı başarılı oldu!

Ancak `id` alanından geri kalan üst veriyi hâlâ çıkarmamız gerekiyor.

### 5.2. Dosya çiftlerinden üst veriyi çıkarmak ve düzenlemek

Önceki `map` işlemimiz veri yapısıyla eşleşmediği için çalışmayacaktır; ancak çalışacak şekilde değiştirebiliriz.

`fromFilePairs()` tarafından tanımlayıcı olarak kullanılan dizede gerçek hasta tanımlayıcısına zaten erişimimiz var; bu nedenle daha önce yaptığımız gibi Path nesnesinden `simpleName` almak yerine bunu kullanarak üst veriyi çıkarabiliriz.

İş akışındaki map işleminin yorum satırını kaldırın ve aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3-4 9 11 13"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id, files ->
            def (sample, replicate, type) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                files
            ]
        }
        .view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3-5 11 13"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Şimdilik eşleştirmeyi yorum satırına alalım, sonra geri döneceğiz!
        ch_files.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        */
        .view()
    ```

Bu sefer map, yalnızca `myFile` yerine `id, files` ile başlıyor ve `tokenize()`, `myFile.simpleName` yerine `id`'ye uygulanıyor.

Ayrıca `tokenize()` satırından `readNum`'u çıkardığımıza dikkat edin; özellikle adlandırmadığımız alt dizeler (soldan başlayarak) sessizce atılacaktır.
Bunu yapabiliriz; çünkü eşleştirilmiş dosyalar artık sıkıca ilişkilendirilmiştir; dolayısıyla üst veri map'inde `readNum`'a artık ihtiyacımız yoktur.

İş akışını çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console

    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

İşte bu: çıktı demetinin ilk konumunda üst veri map'i (`[id:patientA, replicate:1, type:normal]`), ardından tasarlandığı gibi eşleştirilmiş dosyaların demeti yer almaktadır.

Elbette bu yalnızca o belirli dosya çiftini alıp işleyecektir.
Birden fazla çifti işlemeyi denemek istiyorsanız, girdi desenine joker karakterler ekleyip ne olduğunu görebilirsiniz.
Örneğin `data/patientA_rep1_*_R{1,2}_001.fastq.gz` kullanmayı deneyin.

### Özetle

- [`channel.fromFilePairs()` ilgili dosyaları otomatik olarak bulur ve eşleştirir](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Bu, pipeline'ınızda çift uçlu okumaların yönetimini basitleştirir
- Eşleştirilmiş dosyalar `[id, [file1, file2]]` demetleri olarak gruplandırılabilir
- Üst veri çıkarma, tek tek dosyalar yerine eşleştirilmiş dosya kimliğinden yapılabilir

---

## 6. Süreçlerde dosya işlemlerini kullanmak

Şimdi tüm bunları, Nextflow süreci içinde dosya işlemlerinin nasıl kullanılacağını pekiştirmek için basit bir süreçte bir araya getirelim.

Size, üst veri demeti ve bir çift girdi dosyası alan ve bunları analiz eden `ANALYZE_READS` adlı önceden yazılmış bir süreç modülü sunuyoruz.
Bu işlemin dizi hizalaması, varyant çağırma veya bu veri türü için mantıklı olan başka bir adım olduğunu hayal edebiliriz.

Başlayalım.

### 6.1. Süreci içe aktarın ve kodu inceleyin

Bu süreci iş akışında kullanmak için, iş akışı bloğundan önce bir modül include ifadesi eklememiz yeterlidir.

İş akışında aşağıdaki düzenlemeyi yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { ANALYZE_READS } from './modules/analyze_reads.nf'

    workflow {
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Modül dosyasını açarak kodunu inceleyebilirsiniz:

```groovy title="modules/analyze_reads.nf - process example" linenums="1"
#!/usr/bin/env nextflow

process ANALYZE_READS {
    tag { meta.id }

    publishDir { "results/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta.id), path("${meta.id}_stats.txt")

    script:
    """
    echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Read 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Read 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Read 1 size: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    echo "Read 2 size: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
    """
}
```

!!! note "Not"

    `tag` ve `publishDir` yönergeleri, dize enterpolasyonu (`"${...}"`) yerine closure sözdizimi (`{ ... }`) kullanır.
    Bunun nedeni, bu yönergelerin çalışma zamanına kadar kullanılamayan girdi değişkenlerine (`meta`) başvurmasıdır.
    Closure sözdizimi, değerlendirmeyi süreç gerçekten çalışana kadar erteler.

!!! note "Not"

    Üst veri map'imizi kurala uygun olarak `meta` olarak adlandırıyoruz.
    Meta map'lere daha ayrıntılı bir bakış için [Üst Veri ve Meta Map'ler](../metadata/) yan görevine bakın.

### 6.2. Süreci iş akışında çağırmak

Süreç artık iş akışı için kullanılabilir olduğuna göre, çalıştırmak için `ANALYZE_READS` sürecine bir çağrı ekleyebiliriz.

Örnek verilerimiz üzerinde çalıştırmak için iki şey yapmamız gerekecektir:

1. Yeniden eşleştirilmiş kanala bir ad vermek
2. Sürece bir çağrı eklemek

#### 6.2.1. Yeniden eşleştirilmiş girdi kanalını adlandırmak

Daha önce eşleştirme işlemlerini doğrudan girdi kanalına uyguladık.
Yeniden eşleştirilmiş içeriği `ANALYZE_READS` sürecine beslemek (ve bunu açık ve okunması kolay bir şekilde yapmak) için `ch_samples` adlı yeni bir kanal oluşturmak istiyoruz.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana iş akışında `.view()` operatörünü `.set { ch_samples }` ile değiştirin ve kanala adıyla başvurabildiğimizi test eden bir satır ekleyin.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="14 16-17"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
            .set { ch_samples }

        // Geçici: ch_samples'a göz at
        ch_samples.view()
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="14"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_files.map { id,  files ->
           def (sample, replicate, type, readNum) = id.tokenize('_')
           [
               [
                   id: sample,
                   replicate: replicate.replace('rep', ''),
                   type: type
               ],
               files
           ]
        }
        .view()
    }
    ```

Bunu çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [goofy_kirch] DSL2 - revision: 3313283e42

    [[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
    ```

Bu, artık kanala adıyla başvurabildiğimizi doğrulamaktadır.

#### 6.2.2. Veriler üzerinde süreci çağırmak

Şimdi `ch_samples` kanalı üzerinde `ANALYZE_READS` sürecini gerçekten çağıralım.

Ana iş akışında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="23"
        // Analizi çalıştır
        ANALYZE_READS(ch_samples)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="23"
        // Geçici: ch_samples'a göz at
        ch_samples.view()
    ```

Bunu çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [shrivelled_cori] DSL2 - revision: b546a31769

    executor >  local (1)
    [b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
    ```

Bu süreç çıktılarını bir `results` dizinine yayımlamak üzere ayarlanmıştır; bu nedenle oraya bir göz atın.

??? abstract "Dizin ve dosya içeriği"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Sample metadata: patientA
    Replicate: 1
    Type: normal
    Read 1: patientA_rep1_normal_R1_001.fastq.gz
    Read 2: patientA_rep1_normal_R2_001.fastq.gz
    Read 1 size: 10 reads
    Read 2 size: 10 reads
    ```

Süreç girdilerimizi alarak tasarlandığı gibi hasta üst verisini içeren yeni bir dosya oluşturdu.
Mükemmel!

### 6.3. Çok daha fazla hasta eklemek

Elbette bu yalnızca tek bir hasta için tek bir dosya çiftini işlemektedir; bu, Nextflow ile elde etmeyi umduğunuz yüksek verimlilik türü değildir.
Muhtemelen aynı anda çok daha fazla veri işlemek isteyeceksiniz.

`channel.fromPath()`'in girdi olarak bir _glob_ kabul ettiğini hatırlayın; bu, desenle eşleşen herhangi bir sayıda dosyayı kabul edebileceği anlamına gelir.
Bu nedenle tüm hastaları dahil etmek istiyorsak, daha önce kısaca belirtildiği gibi girdi dizesini daha fazla hastayı içerecek şekilde değiştirmemiz yeterlidir.

Mümkün olduğunca geniş kapsamlı olmak istediğimizi varsayalım.
İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="2"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Pipeline'ı tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

    executor >  local (8)
    [d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Sonuçlar dizini artık mevcut tüm veriler için sonuçlar içermelidir.

??? abstract "Dizin içeriği"

    ```console
    results
    ├── patientA
    │   └── patientA_stats.txt
    ├── patientB
    │   └── patientB_stats.txt
    └── patientC
        └── patientC_stats.txt
    ```

Başarı! Tüm hastaları tek seferde analiz ettik! Değil mi?

Belki de hayır.
Daha yakından bakarsanız bir sorunumuz var: patientA için iki replikatımız var; ancak yalnızca bir çıktı dosyası var!
Her seferinde çıktı dosyasının üzerine yazıyoruz.

### 6.4. Yayımlanan dosyaları benzersiz kılmak

Hasta üst verisine erişimimiz olduğundan, bunu yayımlanan dosyaları benzersiz kılmak için kullanabiliriz; ya dizin yapısına ya da dosya adlarının kendisine farklılaştırıcı üst veri ekleyerek.

İş akışında aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Önce"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Burada örnek türleri ve replikatlar için ek dizin düzeyleri kullanma seçeneğini gösteriyoruz; ancak bunu dosya adı düzeyinde de deneyebilirsiniz.

Şimdi pipeline'ı bir kez daha çalıştırın; ancak kendinize temiz bir çalışma alanı sağlamak için önce results dizinini sildiğinizden emin olun:

```bash
rm -r results
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    executor >  local (8)
    [e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
    ```

Şimdi results dizinini kontrol edin:

??? abstract "Dizin içeriği"

    ```console
    results/
    ├── normal
    │   ├── patientA
    │   │   ├── 1
    │   │   │   └── patientA_stats.txt
    │   │   └── 2
    │   │       └── patientA_stats.txt
    │   ├── patientB
    │   │   └── 1
    │   │       └── patientB_stats.txt
    │   └── patientC
    │       └── 1
    │           └── patientC_stats.txt
    └── tumor
        ├── patientA
        │   ├── 1
        │   │   └── patientA_stats.txt
        │   └── 2
        │       └── patientA_stats.txt
        ├── patientB
        │   └── 1
        │       └── patientB_stats.txt
        └── patientC
            └── 1
                └── patientC_stats.txt
    ```

İşte bu: tüm üst verilerimiz düzgün biçimde organize edilmiş. Başarı!

Bu şekilde bir map'e yüklenmiş üst veriniz olduğunda yapabileceğiniz çok daha fazla şey vardır:

1. Hasta özelliklerine göre organize edilmiş çıktı dizinleri oluşturmak
2. Süreçlerde hasta özelliklerine göre kararlar almak
3. Üst veri değerlerine göre verileri bölmek, birleştirmek ve yeniden birleştirmek

Üst veriyi açık tutma ve veriye bağlı (dosya adlarına kodlanmış yerine) bu deseni, sağlam ve bakımı kolay analiz iş akışları oluşturmayı mümkün kılan Nextflow'daki temel en iyi uygulamalardan biridir.
Bu konuda daha fazla bilgi edinmek için [Üst Veri ve Meta Map'ler](../metadata/) yan görevine bakabilirsiniz.

### Özetle

- `publishDir` yönergesi, çıktıları üst veri değerlerine göre düzenleyebilir
- Demetlerdeki üst veri, sonuçların yapılandırılmış organizasyonunu sağlar
- Bu yaklaşım, net veri kaynağı izlenebilirliğiyle bakımı kolay iş akışları oluşturur
- Süreçler, girdi olarak üst veri ve dosya demetleri alabilir
- `tag` yönergesi, yürütme günlüklerinde süreç tanımlaması sağlar
- İş akışı yapısı, kanal oluşturmayı süreç yürütmesinden ayırır

---

## Özet

Bu yan görevde, Nextflow'da dosyalarla nasıl çalışılacağını öğrendiniz: temel işlemlerden dosya koleksiyonlarını yönetmeye yönelik daha gelişmiş tekniklere kadar.

Bu teknikleri kendi çalışmalarınızda uygulamak, özellikle karmaşık adlandırma kurallarına sahip çok sayıda dosyayla çalışırken daha verimli ve bakımı kolay iş akışları oluşturmanızı sağlayacaktır.

### Temel desenler

1.  **Temel Dosya İşlemleri:** `file()` ile Path nesneleri oluşturduk ve ad, uzantı ve üst dizin gibi dosya özelliklerine eriştik; dizeler ile Path nesneleri arasındaki farkı öğrendik.

    - `file()` ile Path nesnesi oluşturma

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Dosya özelliklerini alma

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Uzak Dosyaları Kullanmak**: URI'lar kullanarak yerel ve uzak dosyalar arasında şeffaf biçimde geçiş yapmayı öğrendik; bu, Nextflow'un iş akışı mantığını değiştirmeden çeşitli kaynaklardan dosyaları yönetme yeteneğini göstermektedir.

    - Yerel dosya

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - FTP

    ```groovy
    myFile = file('ftp://path/to/file.txt')
    ```

    - HTTPS

    ```groovy
    myFile = file('https://path/to/file.txt')
    ```

    - Amazon S3

    ```groovy
    myFile = file('s3://path/to/file.txt')
    ```

    - Azure Blob Storage

    ```groovy
    myFile = file('az://path/to/file.txt')
    ```

    - Google Cloud Storage

    ```groovy
    myFile = file('gs://path/to/file.txt')
    ```

3.  **`fromPath()` kanal fabrikasını kullanarak dosya yükleme:** `channel.fromPath()` ile dosya desenlerinden kanallar oluşturduk ve nesne türleri dahil dosya özelliklerini görüntüledik.

    - Dosya deseninden kanal oluşturma

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Dosya özelliklerini alma

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Dosya Adlarından Hasta Üst Verisini Çıkarmak:** Dosya adlarından üst veriyi çıkarmak ve yapılandırmak için `tokenize()` ve `replace()` kullandık; bunları organize map'lere dönüştürdük.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **`channel.fromFilePairs` ile basitleştirme:** İlgili dosyaları otomatik olarak eşleştirmek ve eşleştirilmiş dosya kimliklerinden üst veri çıkarmak için `channel.fromFilePairs()` kullandık.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Süreçlerde Dosya İşlemlerini Kullanmak:** Dosya işlemlerini uygun girdi yönetimiyle Nextflow süreçlerine entegre ettik; çıktıları üst veriye göre düzenlemek için `publishDir` kullandık.

    - Süreç girdileriyle meta map ilişkilendirme

    ```groovy
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
             files
        ]
    }
        .set { ch_samples }

    ANALYZE_READS(ch_samples)
    ```

    - Çıktıları üst veriye göre düzenleme

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Ek kaynaklar

- [Nextflow Belgeleri: Dosyalarla Çalışmak](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Sırada ne var?

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
