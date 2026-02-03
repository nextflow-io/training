# Dosya Girdi İşleme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bilimsel analiz iş akışları genellikle çok sayıda dosyanın işlenmesini içerir.
Nextflow, dosyaları verimli bir şekilde yönetmek için güçlü araçlar sağlar ve verilerinizi minimum kodla düzenlemenize ve işlemenize yardımcı olur.

### Öğrenme hedefleri

Bu yan görevde, Nextflow'un dosyaları nasıl işlediğini, temel dosya işlemlerinden dosya koleksiyonlarıyla çalışmak için daha gelişmiş tekniklere kadar keşfedeceğiz.
Bilimsel analiz pipeline'larında yaygın bir gereksinim olan dosya adlarından meta veri çıkarmayı öğreneceksiniz.

Bu yan görevin sonunda şunları yapabileceksiniz:

- Nextflow'un `file()` metodunu kullanarak dosya yolu dizelerinden Path nesneleri oluşturma
- Ad, uzantı ve üst dizin gibi dosya özelliklerine erişme
- URI'ler kullanarak hem yerel hem de uzak dosyaları şeffaf bir şekilde işleme
- `channel.fromPath()` ve `channel.fromFilePairs()` ile dosya işlemeyi otomatikleştirmek için channel'ları kullanma
- Dize manipülasyonu kullanarak dosya adlarından meta veri çıkarma ve yapılandırma
- Desen eşleştirme ve glob ifadeleri kullanarak ilişkili dosyaları gruplama
- Doğru girdi işleme ile dosya işlemlerini Nextflow process'lerine entegre etme
- Meta veri odaklı dizin yapıları kullanarak process çıktılarını düzenleme

Bu beceriler, farklı dosya girdilerini büyük esneklikle işleyebilen iş akışları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../../hello_nextflow/) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (process'ler, channel'lar, operatörler) rahatça kullanabiliyor olmalısınız

<!-- Meta map'ler yan görevini önce yapmayı öneren kısmı kaldırdım çünkü bu daha sonra doğal olarak çalışıyor -->

---

## 0. Başlangıç

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md)'nda açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/working_with_files
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

`main.nf` adlı basit bir workflow dosyası, iki modül dosyası içeren bir `modules` dizini ve bazı örnek veri dosyaları içeren bir `data` dizini bulacaksınız.

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

Bu dizin, üç hastadan (A, B, C) paired-end dizileme verilerini içerir.

Her hasta için `tumor` (genellikle tümör biyopsilerinden kaynaklanan) veya `normal` (sağlıklı doku veya kandan alınan) türünde örneklerimiz var.
Kanser analiziyle tanışık değilseniz, bunun karşılaştırmalı analizler gerçekleştirmek için eşleştirilmiş tümör/normal örnekleri kullanan deneysel bir modele karşılık geldiğini bilin.

Özellikle A hastası için iki set teknik replik (tekrar) var.

Dizileme veri dosyaları, 'forward reads' ve 'reverse reads' olarak bilinen şeyler için tipik `_R1_` ve `_R2_` adlandırma kuralıyla adlandırılmıştır.

_Bu deneysel tasarıma aşina değilseniz endişelenmeyin, bu eğitimi anlamak için kritik değil._

#### Görevi inceleyin

Göreviniz şunları yapacak bir Nextflow iş akışı yazmaktır:

1. Nextflow'un dosya işleme yöntemlerini kullanarak girdi dosyalarını **yükleme**
2. Dosya adı yapısından meta verileri (hasta kimliği, replik, örnek türü) **çıkarma**
3. `channel.fromFilePairs()` kullanarak eşleştirilmiş dosyaları (R1/R2) **gruplama**
4. Sağlanan bir analiz modülü ile dosyaları **işleme**
5. Çıkarılan meta verilere dayalı bir dizin yapısına çıktıları **düzenleme**

#### Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Temel dosya işlemleri

### 1.1. `.class` ile bir nesnenin türünü belirleme

`main.nf` workflow dosyasına bakın:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Bir dize yolundan Path nesnesi oluştur
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} sınıfı ${myFile.class}"
}
```

Bu, workflow'unda tek bir dosya yoluna başvuran, ardından onu sınıfıyla birlikte konsola yazdıran (herhangi bir process olmadan) mini bir workflow'dur.

??? info "`.class` nedir?"

    Nextflow'da `.class` bize hangi tür nesneyle çalıştığımızı söyler. "Bu ne tür bir şey?" diye sormak gibidir; bir dize mi, sayı mı, dosya mı yoksa başka bir şey mi olduğunu öğrenmek için.
    Bu, sonraki bölümlerde düz bir dize ile Path nesnesi arasındaki farkı göstermemize yardımcı olacaktır.

İş akışını çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz sınıfı java.lang.String
    ```

Gördüğünüz gibi, Nextflow yazdığımız dize yolunu aynen yazdırdı.

Bu sadece metin çıktısıdır; Nextflow henüz bununla özel bir şey yapmadı.
Ayrıca Nextflow'un bakış açısına göre bunun sadece bir dize (`java.lang.String` sınıfında) olduğunu doğruladık.
Bu mantıklı, çünkü henüz Nextflow'a bunun bir dosyaya karşılık geldiğini söylemedik.

### 1.2. `file()` ile Path nesnesi oluşturma

Nextflow'a dosyaları nasıl işleyeceğini, yol dizelerinden [Path nesneleri](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) oluşturarak söyleyebiliriz.

Workflow'umuzda, `file()` metodunu kullanarak `data/patientA_rep1_normal_R1_001.fastq.gz` dize yolunu dosya özelliklerine ve işlemlerine erişim sağlayan bir Path nesnesine dönüştürebiliriz.

`main.nf`'i düzenleyerek dizeyi aşağıdaki gibi `file()` ile sarın:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} sınıfı ${myFile.class}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} sınıfı ${myFile.class}"
    ```

Şimdi iş akışını tekrar çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz sınıfı class sun.nio.fs.UnixPath
    ```

Bu sefer, girdi olarak sağladığımız göreli yol yerine tam mutlak yolu görüyorsunuz.

Nextflow, dizemizi bir Path nesnesine dönüştürdü ve sistemdeki gerçek dosya konumuna çözümledi.
Dosya yolu artık `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz` gibi mutlak olacaktır.

Ayrıca Path nesne sınıfının `sun.nio.fs.UnixPath` olduğuna dikkat edin: bu Nextflow'un yerel dosyaları temsil etme şeklidir.
Daha sonra göreceğimiz gibi, uzak dosyalar farklı sınıf adlarına sahip olacak (HTTP dosyaları için `nextflow.file.http.XPath` gibi), ancak hepsi aynı şekilde çalışır ve iş akışlarınızda aynı şekilde kullanılabilir.

!!! tip "İpucu"

    **Temel fark:**

    - **Yol dizesi**: Nextflow'un karakter olarak işlediği sadece metin
    - **Path nesnesi**: Nextflow'un çalışabileceği akıllı bir dosya referansı

    Şöyle düşünün: yol dizesi kağıda bir adres yazmak gibidir, Path nesnesi ise oraya nasıl gidileceğini bilen ve yolculuk hakkında ayrıntılar söyleyebilen bir GPS cihazına yüklenmiş adresi gibidir.

### 1.3. Dosya özelliklerine erişme

Bu neden yararlı? Şimdi Nextflow, `myFile`'ın sadece bir dize değil, bir Path nesnesi olduğunu anladığından, Path nesnesinin çeşitli özelliklerine erişebiliriz.

Yerleşik dosya özelliklerini yazdırmak için iş akışımızı güncelleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} sınıfı ${myFile.class}"
    ```

İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    Dosya nesne sınıfı: sun.nio.fs.UnixPath
    Dosya adı: patientA_rep1_normal_R1_001.fastq.gz
    Basit ad: patientA_rep1_normal_R1_001
    Uzantı: gz
    Üst dizin: /workspaces/training/side-quests/working_with_files/data
    ```

Çeşitli dosya özelliklerinin yukarıda konsola yazdırıldığını görüyorsunuz.

### 1.4. Dosyayı bir process'e besleme

Dizeler ve Path nesneleri arasındaki fark, process'lerle gerçek iş akışları oluşturmaya başladığınızda kritik hale gelir.
Şimdiye kadar Nextflow'un girdi dosyamızı bir dosya olarak işlediğini doğruladık, ancak bu dosyayı bir process'te gerçekten çalıştırıp çalıştıramayacağımıza bakalım.

#### 1.4.1. Process'i içe aktarma ve kodu inceleme

Size, bir dosya girdisi alan ve içinde kaç satır olduğunu sayan `COUNT_LINES` adlı önceden yazılmış bir process modülü sağlıyoruz.

Process'i iş akışında kullanmak için, workflow bloğundan önce bir include ifadesi eklemeniz yeterlidir:

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
    echo "Dosya işleniyor: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Gördüğünüz gibi, dosyayı açan ve kaç satır içerdiğini sayan oldukça basit bir script.

??? info "`debug true` ne yapar?"

    Process tanımındaki `debug true` yönergesi, Nextflow'un script'inizden çıktıyı (satır sayısı "40" gibi) doğrudan yürütme günlüğünde yazdırmasına neden olur.
    Bu olmadan, yalnızca process yürütme durumunu görürsünüz, ancak script'inizden gerçek çıktıyı görmezsiniz.

    Nextflow process'lerini hata ayıklama hakkında daha fazla bilgi için [Nextflow İş Akışlarında Hata Ayıklama](debugging.md) yan görevine bakın.

#### 1.4.2. `COUNT_LINES`'a bir çağrı ekleme

Artık process iş akışında kullanılabilir olduğundan, girdi dosyası üzerinde çalıştırmak için `COUNT_LINES` process'ine bir çağrı ekleyebiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    ```

Şimdi iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    Dosya nesne sınıfı: class sun.nio.fs.UnixPath
    Dosya adı: patientA_rep1_normal_R1_001.fastq.gz
    Basit ad: patientA_rep1_normal_R1_001
    Uzantı: gz
    Üst dizin: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Dosya işleniyor: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Bu, bir process içinde dosya üzerinde uygun şekilde işlem yapabildiğimizi gösterir.

Özellikle, Nextflow aşağıdaki işlemleri başarıyla gerçekleştirdi:

- Dosyayı çalışma dizinine hazırladı
- .gz dosyasını açtı
- Satırları saydı (bu durumda 40 satır)
- Hatasız tamamlandı

Bu sorunsuz işlemin anahtarı, girdimizin bir dosya olduğunu ve bu şekilde işlenmesi gerektiğini Nextflow'a açıkça söylememizdir.

### 1.5. Temel dosya girdi hatalarını giderme

Bu genellikle Nextflow'a yeni başlayanları yanıltır, bu yüzden yanlış yaptığınızda ne olacağına bakmak için birkaç dakika ayıralım.

Dosya işlemeyi yanlış yapabileceğiniz iki ana yer vardır: iş akışı düzeyinde ve process düzeyinde.

#### 1.5.1. İş akışı düzeyinde hata

Workflow bloğunda girdiyi belirtirken dosyayı bir dize olarak işlemeye geri dönersek ne olacağını görelim.

Path'e özgü print ifadelerini yorumlamayı unutmadan iş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        /*
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
        */

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"

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

Önemli olan kısım şu:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Bir `path` girdisi belirttiğinizde, Nextflow sadece dizeler değil, gerçek dosya referansları geçirdiğinizi doğrular.
Bu hata, `'data/patientA_rep1_normal_R1_001.fastq.gz'`'nin geçerli bir yol değeri olmadığını söylüyor çünkü bu bir dize, Path nesnesi değil.

Nextflow sorunu hemen algıladı ve process'i başlatmadan önce durdurdu.

#### 1.5.2. Process düzeyinde hata

Nextflow'a girdiyi bir dosya olarak işlemesini söylemeyi unutabileceğimiz diğer yer process tanımındadır.

!!! warning "1.5.1'deki iş akışı hatasını koruyun"

    Bu testin doğru çalışması için, iş akışını bozuk durumunda tutun (`file()` yerine düz dize kullanarak).
    Process'te `val` ile birleştirildiğinde, bu aşağıda gösterilen hatayı üretir.

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
      echo "Dosya işleniyor: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Dosya işleniyor: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Dosya işleniyor: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Bu, process hata ayıklama bilgisi çıkarmak üzere ayarlandığı için hata hakkında çok fazla ayrıntı gösterir.

En alakalı bölümler şunlardır:

```console
Command executed:

  set -o pipefail
  echo "Dosya işleniyor: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Dosya işleniyor: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Bu, sistemin dosyayı bulamadığını söylüyor; ancak yola bakarsanız, o konumda o isimde bir dosya var.

Bunu çalıştırdığımızda, Nextflow dize değerini script'e geçirdi, ancak gerçek dosyayı çalışma dizinine _hazırlamadı_.
Bu nedenle process göreli dizeyi, `data/patientA_rep1_normal_R1_001.fastq.gz`'yi kullanmaya çalıştı, ancak bu dosya process çalışma dizini içinde mevcut değil.

Bu iki örnek birlikte, Nextflow'a bir girdinin dosya olarak işlenmesi gerektiğini söylemenin ne kadar önemli olduğunu gösterir.

!!! note "Not"

    Bir sonraki bölüme geçmeden önce her iki kasıtlı hatayı da düzelttiğinizden emin olun.

### Özet

- Yol dizeleri vs Path nesneleri: Dizeler sadece metindir, Path nesneleri akıllı dosya referanslarıdır
- `file()` metodu bir dize yolunu Nextflow'un çalışabileceği bir Path nesnesine dönüştürür
- [Dosya özelliklerini kullanarak](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) `name`, `simpleName`, `extension` ve `parent` gibi dosya özelliklerine erişebilirsiniz
- Dizeler yerine Path nesneleri kullanmak, Nextflow'un iş akışınızdaki dosyaları düzgün bir şekilde yönetmesine olanak tanır
- Process Girdi Sonuçları: Doğru dosya işleme, dosyaların process'ler tarafından kullanılmak üzere doğru şekilde hazırlanmasını ve erişilebilir olmasını sağlamak için dizeler yerine Path nesneleri gerektirir.

---

## 2. Uzak dosyaları kullanma

Nextflow'un temel özelliklerinden biri, yerel dosyalar (aynı makinede) ile internet üzerinden erişilebilen uzak dosyalar arasında sorunsuz geçiş yapabilme yeteneğidir.

Doğru yapıyorsanız, farklı konumlardan gelen dosyaları barındırmak için iş akışınızın mantığını asla değiştirmenize gerek kalmamalıdır.
Uzak bir dosya kullanmak için tek yapmanız gereken, iş akışına sağlarken dosya yolunda uygun öneki belirtmektir.

Örneğin, `/path/to/data` önek içermez, bu da 'normal' bir yerel dosya yolu olduğunu gösterir, oysa `s3://path/to/data` `s3://` önekini içerir, bu da Amazon'un S3 nesne deposunda bulunduğunu gösterir.

Birçok farklı protokol desteklenir:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Bunlardan herhangi birini kullanmak için, dizedeki ilgili öneki belirtmeniz yeterlidir; bu teknik olarak dosya yolu yerine Tekdüzen Kaynak Tanımlayıcısı (URI) olarak adlandırılır.
Nextflow, kimlik doğrulama ve dosyaları doğru yere hazırlama, indirme veya yükleme ve beklediğiniz tüm diğer dosya işlemlerini yönetir.

Bu sistemin temel gücü, herhangi bir pipeline mantığını değiştirmeden ortamlar arasında geçiş yapmamızı sağlamasıdır.
Örneğin, sadece URI'yi değiştirerek uzak depoda bulunan tam ölçekli bir test setine geçmeden önce küçük, yerel bir test setiyle geliştirebilirsiniz.

### 2.1. İnternetten bir dosya kullanma

İş akışımıza sağladığımız yerel yolu, GitHub'da depolanan aynı verilerin bir kopyasına işaret eden bir HTTPS yoluyla değiştirerek bunu test edelim.

!!! warning "Uyarı"

    Bu yalnızca aktif bir internet bağlantınız varsa çalışır.

`main.nf`'i tekrar açın ve girdi yolunu aşağıdaki gibi değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // İnternetten uzak bir dosya kullanma
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    ```

İş akışını çalıştıralım:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    Dosya nesne sınıfı: class nextflow.file.http.XPath
    Dosya adı: patientA_rep1_normal_R1_001.fastq.gz
    Basit ad: patientA_rep1_normal_R1_001
    Uzantı: gz
    Üst dizin: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Dosya işleniyor: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Çalışıyor! Çok az şeyin değiştiğini görebilirsiniz.

Konsol çıktısındaki tek fark, path nesne sınıfının artık `nextflow.file.http.XPath` olmasıdır, yerel yol için sınıf `sun.nio.fs.UnixPath` idi.
Bu sınıfları hatırlamanıza gerek yok; bunu sadece Nextflow'un farklı konumları tanımladığını ve uygun şekilde işlediğini göstermek için belirtiyoruz.

Perde arkasında, Nextflow dosyayı çalışma dizini içinde bulunan bir hazırlama dizinine indirdi.
Bu hazırlanan dosya daha sonra yerel bir dosya olarak işlenebilir ve ilgili process dizinine sembolik olarak bağlanabilir.

Bunun burada gerçekleştiğini, process'in hash değerinde bulunan çalışma dizininin içeriğine bakarak doğrulayabilirsiniz.

??? abstract "Çalışma dizini içeriği"

    Process hash'i `8a/2ab7ca` ise, çalışma dizinini keşfedebilirsiniz:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Sembolik bağlantı, Nextflow'un otomatik olarak indirdiği uzak dosyanın hazırlanmış bir kopyasına işaret eder.

Daha büyük dosyalar için, indirme adımı yerel dosyalarda çalışmaya kıyasla biraz ekstra zaman alacaktır.
Ancak, Nextflow gereksiz indirmeleri önlemek için hazırlanmış bir kopyasının zaten olup olmadığını kontrol eder.
Bu nedenle, aynı dosyayı tekrar çalıştırır ve hazırlanmış dosyayı silmediyseniz, Nextflow hazırlanmış kopyayı kullanır.

Bu, Nextflow kullanarak yerel ve uzak veriler arasında geçiş yapmanın ne kadar kolay olduğunu gösterir; bu Nextflow'un temel bir özelliğidir.

!!! note "Not"

    Bu ilkenin tek önemli istisnası, HTTPS birden fazla dosyayı listeleyemediği için HTTPS ile glob desenleri veya dizin yolları kullanamamanızdır.
    Ancak, blob depolama (`s3://`, `az://`, `gs://`) gibi diğer depolama protokolleri hem glob'ları hem de dizin yollarını kullanabilir.

    Bulut depolamayla glob desenlerini şu şekilde kullanabilirsiniz:

    ```groovy title="Bulut depolama örnekleri (bu ortamda çalıştırılamaz)"
    // Glob desenleriyle S3 - birden fazla dosyayla eşleşir
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Glob desenleriyle Azure Blob Storage
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Glob desenleriyle Google Cloud Storage
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Bir sonraki bölümde glob'larla pratikte nasıl çalışacağınızı göstereceğiz.

### 2.2. Yerel dosyaya geri dönme

Bu yan görevin geri kalanında yerel örnek dosyalarımızı kullanmaya devam edeceğiz, bu yüzden iş akışı girdisini orijinal dosyaya geri çevirelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    ```

### Özet

- Uzak verilere bir URI (HTTP, FTP, S3, Azure, Google Cloud) kullanılarak erişilir
- Nextflow, bu yollar process'lere beslendiği sürece verileri otomatik olarak indirir ve doğru yere hazırlar
- Uzak dosyaları indirmek veya yüklemek için mantık yazmayın!
- Yerel ve uzak dosyalar farklı nesne türleri üretir ancak aynı şekilde çalışır
- **Önemli**: HTTP/HTTPS yalnızca tekil dosyalarla çalışır (glob desenleri yok)
- Bulut depolama (S3, Azure, GCS) hem tekil dosyaları hem de glob desenlerini destekler
- Protokol gerekli işlemlerinizi desteklediği sürece kod mantığını değiştirmeden yerel ve uzak veri kaynakları arasında sorunsuz geçiş yapabilirsiniz

---

## 3. `fromPath()` channel factory kullanma

Şimdiye kadar tek seferde tek bir dosyayla çalışıyorduk, ancak Nextflow'da genellikle işlenecek birden fazla girdi dosyasıyla bir girdi channel'ı oluşturmak isteyeceğiz.

Bunu yapmanın naif bir yolu, `file()` metodunu [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) ile şu şekilde birleştirmek olurdu:

```groovy title="Sözdizimi örneği"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Bu çalışır, ancak hantal.

!!! tip "`file()` vs `channel.fromPath()` ne zaman kullanılır"

    - Doğrudan manipülasyon için tek bir Path nesnesine ihtiyacınız olduğunda (bir dosyanın var olup olmadığını kontrol etme, özelliklerini okuma veya tek bir process çağrısına geçirme) `file()` kullanın
    - Özellikle glob desenleriyle birden fazla dosya tutabilen veya dosyalar birden fazla process'ten geçecekse bir channel'a ihtiyacınız olduğunda `channel.fromPath()` kullanın

[`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) burada devreye girer: bir veya daha fazla statik dosya dizesinden ve glob desenlerinden bir channel oluşturmak için ihtiyacımız olan tüm işlevselliği bir araya getiren kullanışlı bir channel factory.

### 3.1. Channel factory'yi ekleme

İş akışımızı `channel.fromPath` kullanacak şekilde güncelleyelim.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Dosya bulundu: $myFile" }

        // Dosya özelliklerini yazdır
        /* Şimdilik bunları yorumlayın, daha sonra geri döneceğiz!
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
        */

        // Dosyadaki satırları say
        // COUNT_LINES(myFile)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Bir dize yolundan Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özelliklerini yazdır
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

Ayrıca özellikleri yazdıran kodu şimdilik yorumladık ve bunun yerine sadece dosya adını yazdırmak için bir `.view` ifadesi ekledik.

İş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Dosya bulundu: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Gördüğünüz gibi, dosya yolu channel'da `Path` türü nesne olarak yükleniyor.
Bu, `file()`'ın yapacağına benzer, ancak şimdi istersek daha fazla dosya yükleyebileceğimiz bir channel'ımız var.

`channel.fromPath()` kullanmak, bir dosya listesiyle doldurulmuş yeni bir channel oluşturmanın uygun bir yoludur.

### 3.2. Channel'daki dosyaların özelliklerini görüntüleme

Channel factory kullanımımızın ilk geçişinde, kodu basitleştirdik ve sadece dosya adını yazdırdık.

Tam dosya özelliklerini yazdırmaya geri dönelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "Dosya nesne sınıfı: ${myFile.class}"
            println "Dosya adı: ${myFile.name}"
            println "Basit ad: ${myFile.simpleName}"
            println "Uzantı: ${myFile.extension}"
            println "Üst dizin: ${myFile.parent}"
        }

        // Dosyadaki satırları say
        COUNT_LINES(ch_files)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Dosya bulundu: $myFile" }

        // Dosyadaki satırları say
        // COUNT_LINES(ch_files)
    ```

Ayrıca dosya işlemenin channel tabanlı yaklaşımımızla hala doğru çalıştığını doğrulamak için `COUNT_LINES` process çağrısını yeniden etkinleştiriyoruz.

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
    Dosya nesne sınıfı: sun.nio.fs.UnixPath
    Dosya adı: patientA_rep1_normal_R1_001.fastq.gz
    Basit ad: patientA_rep1_normal_R1_001
    Uzantı: gz
    Üst dizin: /workspaces/training/side-quests/working_with_files/data
    Dosya işleniyor: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

İşte, öncekiyle aynı sonuçlar ancak şimdi dosya bir channel'da, bu yüzden daha fazla ekleyebiliriz.

### 3.3. Birden fazla dosyayla eşleşmek için glob kullanma

Channel'a daha fazla dosya yüklemenin birkaç yolu var.
Burada, joker karakterlere dayalı olarak dosya ve dizin adlarını eşleştirmenin ve almanın uygun bir yolu olan glob desenlerini nasıl kullanacağınızı göstereceğiz.
Bu desenleri eşleştirme işlemine "globbing" veya "dosya adı genişletme" denir.

!!! note "Not"

    Daha önce belirtildiği gibi, Nextflow girdi ve çıktı dosyalarını yönetmek için çoğu durumda globbing'i destekler, ancak HTTPS birden fazla dosyayı listeleyemediği için HTTPS dosya yollarıyla istisna olarak desteklemez.

Diyelim ki belirli bir hastayla, `patientA` ile ilişkili bir çift dosyadaki her iki dosyayı da almak istiyoruz:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Dosya adları arasındaki tek fark replik numarası olduğundan, yani `R`'den sonraki sayı, aşağıdaki gibi sayı yerine joker karakter `*`'yi kullanabiliriz:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Bu, ihtiyacımız olan glob desenidir.

Şimdi tek yapmamız gereken channel factory'deki dosya yolunu bu glob desenini kullanacak şekilde güncellemektir:

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
    Dosya nesne sınıfı: class sun.nio.fs.UnixPath
    Dosya adı: patientA_rep1_normal_R1_001.fastq.gz
    Basit ad: patientA_rep1_normal_R1_001
    Uzantı: gz
    Üst dizin: /workspaces/training/side-quests/working_with_files/data
    Dosya nesne sınıfı: class sun.nio.fs.UnixPath
    Dosya adı: patientA_rep1_normal_R2_001.fastq.gz
    Basit ad: patientA_rep1_normal_R2_001
    Uzantı: gz
    Üst dizin: /workspaces/training/side-quests/working_with_files/data
    Dosya işleniyor: patientA_rep1_normal_R1_001.fastq.gz
    40

    Dosya işleniyor: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Gördüğünüz gibi, şimdi channel'ımızda iki Path nesnemiz var, bu da Nextflow'un dosya adı genişletmesini doğru şekilde yaptığını ve her iki dosyayı da beklendiği gibi yükleyip işlediğini gösteriyor.

Bu yöntemi kullanarak, sadece glob desenini değiştirerek istediğimiz kadar çok veya az dosyayı alabiliriz. Daha cömert yaparsak, örneğin dosya adlarının tüm değişken kısımlarını `*` ile değiştirerek (_örneğin_ `data/patient*_rep*_*_R*_001.fastq.gz`) `data` dizinindeki tüm örnek dosyaları yakalayabiliriz.

### Özet

- `channel.fromPath()` bir desenle eşleşen dosyalarla bir channel oluşturur
- Her dosya channel'da ayrı bir öğe olarak yayılır
- Birden fazla dosyayla eşleşmek için bir glob deseni kullanabiliriz
- Dosyalar otomatik olarak tam özelliklere sahip Path nesnelerine dönüştürülür
- `.view()` metodu channel içeriklerinin incelenmesine olanak tanır

---

## 4. Dosya adlarından temel meta veri çıkarma

Çoğu bilimsel alanda, meta verilerin verileri içeren dosyaların adlarına kodlanması çok yaygındır.
Örneğin, biyoinformatikte, dizileme verilerini içeren dosyalar genellikle örnek, koşul, replik ve okuma numarası hakkında bilgi kodlayan bir şekilde adlandırılır.

Dosya adları tutarlı bir kurala göre oluşturulmuşsa, bu meta veriyi standart bir şekilde çıkarabilir ve analiziniz sırasında kullanabilirsiniz.
Bu elbette büyük bir "eğer"dir ve dosya adı yapısına güvenirken çok dikkatli olmalısınız; ancak gerçek şu ki bu yaklaşım çok yaygın olarak kullanılmaktadır, bu yüzden Nextflow'da nasıl yapıldığına bir göz atalım.

Örnek verilerimiz söz konusu olduğunda, dosya adlarının tutarlı bir şekilde yapılandırılmış meta veri içerdiğini biliyoruz.
Örneğin, `patientA_rep1_normal_R2_001` dosya adı aşağıdakileri kodlar:

- hasta kimliği: `patientA`
- replik kimliği: `rep1`
- örnek türü: `normal` (`tumor`'un aksine)
- okuma seti: `R1` (`R2`'nin aksine)

Bu bilgiyi üç adımda almak için iş akışımızı değiştireceğiz:

1. Meta veriyi içeren dosyanın `simpleName`'ini alma
2. `tokenize()` adlı bir metot kullanarak meta veriyi ayırma
3. Meta veriyi düzenlemek için bir map kullanma

!!! warning "Uyarı"

    Hasta adları veya diğer tanımlayıcı özellikler gibi hassas bilgileri asla dosya adlarına kodlamamalısınız, çünkü bu hasta gizliliğini veya diğer ilgili güvenlik kısıtlamalarını tehlikeye atabilir.

### 4.1. `simpleName`'i alma

`simpleName`, yolu ve uzantısından soyulmuş dosya adına karşılık gelen bir dosya özelliğidir.

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
            println "Dosya nesne sınıfı: ${myFile.class}"
            println "Dosya adı: ${myFile.name}"
            println "Basit ad: ${myFile.simpleName}"
            println "Uzantı: ${myFile.extension}"
            println "Üst dizin: ${myFile.parent}"
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
    Dosya işleniyor: patientA_rep1_normal_R1_001.fastq.gz
    40

    Dosya işleniyor: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Channel'daki her öğe artık `simpleName`'i ve orijinal dosya nesnesini içeren bir tuple'dır.

### 4.2. `simplename`'den meta veriyi çıkarma

Bu noktada, istediğimiz meta veri `simplename`'e gömülü, ancak tek tek öğelere doğrudan erişemiyoruz.
Bu yüzden `simplename`'i bileşenlerine ayırmamız gerekiyor.
Neyse ki, bu bileşenler orijinal dosya adında sadece alt çizgilerle ayrılmış, bu yüzden bu görev için mükemmel olan `tokenize()` adlı yaygın bir Nextflow metodunu uygulayabiliriz.

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

`tokenize()` metodu, alt çizgi bulduğu her yerde `simpleName` dizesini böler ve alt dizeleri içeren bir liste döndürür.

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
    Dosya işleniyor: patientA_rep1_normal_R2_001.fastq.gz
    40

    Dosya işleniyor: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Şimdi channel'ımızdaki her öğenin tuple'ı meta veri listesini (_örneğin_ `[patientA, rep1, normal, R1, 001]`) ve orijinal dosya nesnesini içeriyor.

Harika!
Hasta bilgilerimizi tek bir dizeden bir dize listesine ayırdık.
Artık hasta bilgilerinin her bir parçasını ayrı ayrı işleyebiliriz.

### 4.3. Meta veriyi düzenlemek için map kullanma

Meta verimiz şu anda sadece düz bir liste.
Kullanması kolay ama okumak zor.

```console
[patientA, rep1, normal, R1, 001]
```

Dizin 3'teki öğe nedir? Meta veri yapısının orijinal açıklamasına başvurmadan söyleyebilir misiniz?

Bu, her öğenin bir dizi anahtar ve ilişkili değerlere sahip olduğu bir anahtar-değer deposu kullanmak için harika bir fırsattır, böylece karşılık gelen değeri almak için her anahtara kolayca başvurabilirsiniz.

Örneğimizde, bu şu organizasyondan:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

Şu organizasyona geçmek anlamına gelir:

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

- **Yapı bozma ataması**: `def (patient, replicate, type, readNum) = ...` tokenize edilmiş değerleri tek satırda adlandırılmış değişkenlere çıkarır
- **Map literal sözdizimi**: `[id: patient, replicate: ...]` her anahtarın (örneğin `id`) bir değerle (örneğin `patient`) ilişkilendirildiği bir map oluşturur
- **İç içe yapı**: Dıştaki liste `[..., myFile]` meta veri map'ini orijinal dosya nesnesiyle eşleştirir

Ayrıca gereksiz olan bazı karakterleri kaldırmak için `replace()` adlı bir dize değiştirme metodu kullanarak meta veri dizelerinin birkaçını basitleştirdik (_örneğin_ replik kimliklerinden sadece sayıyı tutmak için `replicate.replace('rep', '')`).

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
    Dosya işleniyor: patientA_rep1_normal_R2_001.fastq.gz
    40

    Dosya işleniyor: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Şimdi meta veri düzgün bir şekilde etiketlenmiş (_örneğin_ `[id:patientA, replicate:1, type:normal, readNum:2]`) bu yüzden neyin ne olduğunu söylemek çok daha kolay.

İş akışında meta veri öğelerini gerçekten kullanmak da çok daha kolay olacak ve kodumuzu okumayı ve bakımını kolaylaştıracak.

### Özet

- Nextflow'da dosya adlarını tam bir programlama dilinin gücüyle işleyebiliriz
- İlgili bilgileri çıkarmak için dosya adlarını dize olarak işleyebiliriz
- `tokenize()` ve `replace()` gibi metotların kullanımı, dosya adındaki dizeleri manipüle etmemize olanak tanır
- `.map()` işlemi yapıyı korurken channel öğelerini dönüştürür
- Yapılandırılmış meta veri (map'ler) konumsal listelerden daha okunabilir ve bakımı kolay kod sağlar

Şimdi, eşleştirilmiş veri dosyalarını nasıl işleyeceğimize bakacağız.

---

## 5. Eşleştirilmiş veri dosyalarını işleme

Birçok deneysel tasarım, açıkça eşleştirilmiş bir şekilde işlenmekten fayda sağlayan eşleştirilmiş veri dosyaları üretir.
Örneğin, biyoinformatikte, dizileme verileri genellikle eşleştirilmiş okumalar biçiminde üretilir, yani aynı DNA parçasından kaynaklanan dizi dizileri (genellikle zıt uçlardan okudukları için 'forward' ve 'reverse' olarak adlandırılır).

Bu, R1 ve R2'nin iki okuma setine atıfta bulunduğu örnek verilerimiz için de geçerlidir.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow, paylaşılan bir adlandırma desenine göre dosyaları otomatik olarak gruplayan `channel.fromFilePairs()` adlı eşleştirilmiş dosyalarla çalışmak için özel bir channel factory sağlar. Bu, eşleştirilmiş dosyaları daha az çabayla daha sıkı bir şekilde ilişkilendirmenize olanak tanır.

Bundan yararlanmak için iş akışımızı değiştireceğiz.
İki adım alacak:

1. Channel factory'yi `channel.fromFilePairs()`'a geçirme
2. Meta veriyi çıkarma ve eşleme

### 5.1. Channel factory'yi `channel.fromFilePairs()`'a geçirme

`channel.fromFilePairs` kullanmak için, Nextflow'un bir çiftteki iki üyeyi tanımlamak için kullanması gereken deseni belirtmemiz gerekir.

Örnek verilerimize geri dönersek, adlandırma desenini aşağıdaki gibi resmileştirebiliriz:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Bu, daha önce kullandığımız glob desenine benzer, ancak bu özellikle çiftin iki üyesini tanımlayan alt dizeleri (R'den hemen sonra `1` veya `2` gelen) numaralandırır.

İş akışı `main.nf`'i buna göre güncelleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Eşlemeyi şimdilik yorumlayın, daha sonra geri döneceğiz!
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

Channel factory'yi değiştirdik ve dosya eşleştirme desenini uyarladık, ve bunu yaparken map işlemini yorumladık.
Daha sonra birkaç değişiklikle bunu geri ekleyeceğiz.

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

Hata mesajının ilgili kısmı burada:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Bunun nedeni channel factory'yi değiştirmiş olmamızdır.
Şimdiye kadar, orijinal girdi channel'ı yalnızca dosya yollarını içeriyordu.
Yaptığımız tüm meta veri manipülasyonu aslında channel içeriklerini etkilemedi.

Artık `.fromFilePairs` channel factory'yi kullandığımıza göre, ortaya çıkan channel'ın içerikleri farklı.
Sadece bir channel öğesi görüyoruz, iki öğe içeren bir tuple'dan oluşuyor: iki dosya tarafından paylaşılan `simpleName` kısmı, bir tanımlayıcı görevi görüyor, ve iki dosya nesnesini içeren bir tuple, `id, [ file1, file2 ]` formatında.

Bu harika, çünkü Nextflow paylaşılan öneki inceleyerek ve bunu hasta tanımlayıcısı olarak kullanarak hasta adını çıkarmanın zor işini yaptı.

Ancak, mevcut iş akışımızı bozuyor.
`COUNT_LINES`'ı process'i değiştirmeden aynı şekilde çalıştırmak istersek, dosya yollarını çıkarmak için bir eşleme işlemi uygulamamız gerekir.
Ancak bunu yapmayacağız, çünkü nihai hedefimiz dosya çiftlerini uygun şekilde işleyen farklı bir process, `ANALYZE_READS` kullanmaktır.

Bu yüzden sadece `COUNT_LINES`'a yapılan çağrıyı yorumlayalım (veya silelim) ve devam edelim.

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

`COUNT_LINES` include ifadesini de yorumlayabilir veya silebilirsiniz, ancak bunun işlevsel bir etkisi olmayacaktır.

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

Yaşasın, bu sefer iş akışı başarılı!

Ancak, yine de `id` alanından meta verinin geri kalanını almamız gerekiyor.

### 5.2. Dosya çiftlerinden meta veriyi çıkarma ve düzenleme

Önceki `map` işlemimiz veri yapısıyla eşleşmediği için çalışmayacak, ancak çalışması için değiştirebiliriz.

Zaten `fromFilePairs()`'ın tanımlayıcı olarak kullandığı dizede gerçek hasta tanımlayıcısına erişimimiz var, bu yüzden daha önce yaptığımız gibi Path nesnesinden `simpleName`'i almadan meta veriyi çıkarmak için bunu kullanabiliriz.

İş akışındaki map işlemini yorumdan çıkarın ve aşağıdaki düzenlemeleri yapın:

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
        /* Eşlemeyi şimdilik yorumlayın, daha sonra geri döneceğiz!
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

Bu sefer map sadece `myFile` yerine `id, files`'dan başlıyor ve `tokenize()` `myFile.simpleName` yerine `id`'ye uygulanıyor.

Ayrıca `readNum`'u `tokenize()` satırından çıkardığımıza dikkat edin; özellikle adlandırmadığımız alt dizeler (soldan başlayarak) sessizce atılacaktır.
Bunu yapabiliriz çünkü eşleştirilmiş dosyalar artık sıkıca ilişkili, bu yüzden meta veri map'inde artık `readNum`'a ihtiyacımız yok.

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

İşte burada: çıktı tuple'ının ilk pozisyonunda meta veri map'imiz (`[id:patientA, replicate:1, type:normal]`) var, ardından amaçlandığı gibi eşleştirilmiş dosyalar tuple'ı geliyor.

Elbette, bu yalnızca belirli dosya çiftini alacak ve işleyecektir.
Birden fazla çifti işlemeyi denemek isterseniz, girdi desenine joker karakterler eklemeyi ve ne olacağını görmeyi deneyebilirsiniz.
Örneğin, `data/patientA_rep1_*_R{1,2}_001.fastq.gz` kullanmayı deneyin

### Özet

- [`channel.fromFilePairs()` ilişkili dosyaları otomatik olarak bulur ve eşleştirir](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Bu, pipeline'ınızda paired-end okumaları işlemeyi basitleştirir
- Eşleştirilmiş dosyalar `[id, [file1, file2]]` tuple'ları olarak gruplanabilir
- Meta veri çıkarma, tek tek dosyalar yerine eşleştirilmiş dosya kimliğinden yapılabilir

---

## 6. Process'lerde dosya işlemlerini kullanma

Şimdi, bir Nextflow process içinde dosya işlemlerini nasıl kullanacağımızı pekiştirmek için tüm bunları basit bir process'te bir araya getirelim.

Size, bir meta veri tuple'ı ve bir çift girdi dosyası alan ve bunları analiz eden `ANALYZE_READS` adlı önceden yazılmış bir process modülü sağlıyoruz.
Bu işlemin dizilim hizalaması, varyant çağırma veya bu veri türü için anlamlı olan başka herhangi bir adımı yaptığını hayal edebiliriz.

Başlayalım.

### 6.1. Process'i içe aktarma ve kodu inceleme

Bu process'i iş akışında kullanmak için, workflow bloğundan önce bir modül include ifadesi eklememiz yeterlidir.

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

```groovy title="modules/analyze_reads.nf - process örneği" linenums="1"
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
    echo "Örnek meta verisi: ${meta.id}" > ${meta.id}_stats.txt
    echo "Replik: ${meta.replicate}" >> ${meta.id}_stats.txt
    echo "Tür: ${meta.type}" >> ${meta.id}_stats.txt
    echo "Okuma 1: ${files[0]}" >> ${meta.id}_stats.txt
    echo "Okuma 2: ${files[1]}" >> ${meta.id}_stats.txt
    echo "Okuma 1 boyutu: \$(gunzip -dc ${files[0]} | wc -l | awk '{print \$1/4}') okuma" >> ${meta.id}_stats.txt
    echo "Okuma 2 boyutu: \$(gunzip -dc ${files[1]} | wc -l | awk '{print \$1/4}') okuma" >> ${meta.id}_stats.txt
    """
}
```

!!! note "Not"

    `tag` ve `publishDir` yönergeleri, dize interpolasyonu (`"${...}"`) yerine closure sözdizimi (`{ ... }`) kullanır.
    Bunun nedeni, bu yönergelerin çalışma zamanına kadar mevcut olmayan girdi değişkenlerine (`meta`) başvurmasıdır.
    Closure sözdizimi, değerlendirmeyi process gerçekten çalışana kadar erteler.

!!! note "Not"

    Meta veri map'imizi kural olarak `meta` olarak adlandırıyoruz.
    Meta map'lere daha derin bir bakış için [Meta veri ve meta map'ler](./metadata.md) yan görevine bakın.

### 6.2. İş akışında process'i çağırma

Artık process iş akışında kullanılabilir olduğundan, çalıştırmak için `ANALYZE_READS` process'ine bir çağrı ekleyebiliriz.

Örnek verilerimizde çalıştırmak için iki şey yapmamız gerekecek:

1. Yeniden eşlenen channel'a bir ad verme
2. Process'e bir çağrı ekleme

#### 6.2.1. Yeniden eşlenen girdi channel'ını adlandırma

Daha önce eşleme manipülasyonlarını doğrudan girdi channel'ına uygulamıştık.
Yeniden eşlenen içerikleri `ANALYZE_READS` process'ine beslemek için (ve bunu açık ve okunması kolay bir şekilde yapmak için) `ch_samples` adlı yeni bir channel oluşturmak istiyoruz.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana iş akışında, `.view()` operatörünü `.set { ch_samples }` ile değiştirin ve channel'a adla başvurabileceğimizi test eden bir satır ekleyin.

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

        // Geçici: ch_samples'a bakış
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

Bu, artık channel'a adla başvurabileceğimizi doğrular.

#### 6.2.2. Veri üzerinde process'i çağırma

Şimdi `ch_samples` channel'ı üzerinde `ANALYZE_READS` process'ini gerçekten çağıralım.

Ana iş akışında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="23"
        // Analizi çalıştır
        ANALYZE_READS(ch_samples)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="23"
        // Geçici: ch_samples'a bakış
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

Bu process çıktılarını bir `results` dizinine yayınlamak üzere ayarlanmış, bu yüzden oraya bakın.

??? abstract "Dizin ve dosya içeriği"

    ```console
    results
    └── patientA
        └── patientA_stats.txt
    ```

    ```txt title="patientA_stats.txt"
    Örnek meta verisi: patientA
    Replik: 1
    Tür: normal
    Okuma 1: patientA_rep1_normal_R1_001.fastq.gz
    Okuma 2: patientA_rep1_normal_R2_001.fastq.gz
    Okuma 1 boyutu: 10 okuma
    Okuma 2 boyutu: 10 okuma
    ```

Process girdilerimizi aldı ve tasarlandığı gibi hasta meta verilerini içeren yeni bir dosya oluşturdu.
Mükemmel!

### 6.3. Çok daha fazla hasta dahil etme

Elbette, bu sadece tek bir hasta için tek bir dosya çiftini işliyor, bu da Nextflow ile elde etmeyi umduğunuz yüksek verimlilik türü değil.
Muhtemelen aynı anda çok daha fazla veri işlemek isteyeceksiniz.

`channel.fromPath()`'in girdi olarak bir _glob_ kabul ettiğini unutmayın, bu da desenle eşleşen herhangi sayıda dosyayı kabul edebileceği anlamına gelir.
Bu nedenle, tüm hastaları dahil etmek istiyorsak, daha önce geçerken belirtildiği gibi daha fazla hasta içerecek şekilde girdi dizesini değiştirebiliriz.

Mümkün olduğunca açgözlü olmak istediğimizi varsayalım.
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

Results dizini artık tüm mevcut veriler için sonuçlar içermelidir.

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

Başarılı! Tüm hastaları tek seferde analiz ettik! Değil mi?

Belki de hayır.
Daha yakından bakarsanız, bir sorunumuz var: patientA için iki replikimiz var, ama sadece bir çıktı dosyası var!
Her seferinde çıktı dosyasını üzerine yazıyoruz.

### 6.4. Yayınlanan dosyaları benzersiz yapma

Hasta meta verisine erişimimiz olduğundan, bunu dizin yapısında veya dosya adlarının kendisinde farklılaştırıcı meta veri ekleyerek yayınlanan dosyaları benzersiz yapmak için kullanabiliriz.

İş akışında aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Önce"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Burada örnek türlerini ve replikleri hesaba katmak için ek dizin seviyeleri kullanma seçeneğini gösteriyoruz, ancak bunu dosya adı seviyesinde yapmayı da deneyebilirsiniz.

Şimdi pipeline'ı bir kez daha çalıştırın, ancak kendinize temiz bir çalışma alanı sağlamak için önce results dizinini sildiğinizden emin olun:

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
    │           └── patientA_stats.txt
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
                └── patientA_stats.txt
    ```

İşte burada, tüm meta verilerimiz düzgünce organize edilmiş. Bu başarı!

Meta verinizi böyle bir map'e yükledikten sonra yapabileceğiniz çok daha fazla şey var:

1. Hasta özelliklerine göre organize çıktı dizinleri oluşturma
2. Hasta özelliklerine göre process'lerde kararlar verme
3. Meta veri değerlerine göre verileri bölme, birleştirme ve yeniden birleştirme

Meta veriyi açık ve veriye bağlı tutma deseni (dosya adlarına kodlanmak yerine), sağlam, bakımı kolay analiz iş akışları oluşturmayı sağlayan Nextflow'un temel en iyi uygulamasıdır.
Bu hakkında daha fazla bilgiyi [Meta veri ve meta map'ler](./metadata.md) yan görevinde öğrenebilirsiniz.

### Özet

- `publishDir` yönergesi, çıktıları meta veri değerlerine göre organize edebilir
- Tuple'lardaki meta veri, sonuçların yapılandırılmış organizasyonuna olanak tanır
- Bu yaklaşım, açık veri kaynağı ile bakımı kolay iş akışları oluşturur
- Process'ler meta veri ve dosya tuple'larını girdi olarak alabilir
- `tag` yönergesi, yürütme günlüklerinde process tanımlama sağlar
- İş akışı yapısı, channel oluşturmayı process yürütmesinden ayırır

---

## Özet

Bu yan görevde, Nextflow'da dosyalarla nasıl çalışacağınızı, temel işlemlerden dosya koleksiyonlarıyla çalışmak için daha gelişmiş tekniklere kadar öğrendiniz.

Bu teknikleri kendi çalışmanıza uygulamak, özellikle karmaşık adlandırma kurallarına sahip çok sayıda dosyayla çalışırken daha verimli ve bakımı kolay iş akışları oluşturmanıza olanak tanır.

### Temel desenler

1.  **Temel Dosya İşlemleri:** `file()` ile Path nesneleri oluşturduk ve ad, uzantı ve üst dizin gibi dosya özelliklerine eriştik, dizeler ve Path nesneleri arasındaki farkı öğrendik.

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

2.  **Uzak Dosyaları Kullanma**: URI'ler kullanarak yerel ve uzak dosyalar arasında şeffaf geçiş yapmayı öğrendik, Nextflow'un iş akışı mantığını değiştirmeden çeşitli kaynaklardan dosyaları işleme yeteneğini gösterdik.

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

3.  **`fromPath()` channel factory ile dosyaları yükleme:** `channel.fromPath()` ile dosya desenlerinden channel'lar oluşturduk ve nesne türleri dahil dosya özelliklerini görüntüledik.

    - Dosya deseninden channel oluşturma

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Dosya özelliklerini alma

    ```groovy
     ch_files.view { myFile ->
        println "Dosya nesne sınıfı: ${myFile.class}"
        println "Dosya adı: ${myFile.name}"
        println "Basit ad: ${myFile.simpleName}"
        println "Uzantı: ${myFile.extension}"
        println "Üst dizin: ${myFile.parent}"
    }
    ```

4.  **Dosya Adlarından Hasta Meta Verisi Çıkarma:** Dosya adlarından meta veriyi çıkarmak ve yapılandırmak için `tokenize()` ve `replace()` kullandık, bunları organize map'lere dönüştürdük.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **channel.fromFilePairs ile Basitleştirme:** İlişkili dosyaları otomatik olarak eşleştirmek ve eşleştirilmiş dosya kimliklerinden meta veri çıkarmak için `channel.fromFilePairs()` kullandık.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Process'lerde Dosya İşlemlerini Kullanma:** Doğru girdi işleme ile dosya işlemlerini Nextflow process'lerine entegre ettik, çıktıları meta veriye göre organize etmek için `publishDir` kullandık.

    - Process girdileriyle bir meta map ilişkilendirme

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

    - Çıktıları meta veriye göre düzenleme

    ```groovy
    publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

### Ek kaynaklar

- [Nextflow Belgeleri: Dosyalarla Çalışma](https://www.nextflow.io/docs/latest/working-with-files.html)
- [channel.fromPath](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)
- [channel.fromFilePairs](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
