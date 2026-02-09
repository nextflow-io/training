# Dosya girdi işleme

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bilimsel analiz iş akışları genellikle çok sayıda dosyanın işlenmesini içerir.
Nextflow, dosyaları verimli bir şekilde işlemek için güçlü araçlar sağlayarak verilerinizi minimum kodla düzenlemenize ve işlemenize yardımcı olur.

### Öğrenme hedefleri

Bu yan görevde, Nextflow'un dosyaları nasıl işlediğini, temel dosya işlemlerinden dosya koleksiyonlarıyla çalışmak için daha gelişmiş tekniklere kadar keşfedeceğiz.
Bilimsel analiz boru hatlarında yaygın bir gereklilik olan dosya adlarından metadata çıkarmayı öğreneceksiniz.

Bu yan görevin sonunda şunları yapabileceksiniz:

- Nextflow'un `file()` metodunu kullanarak dosya yolu dizelerinden Path nesneleri oluşturmak
- Ad, uzantı ve üst dizin gibi dosya özniteliklerine erişmek
- URI'ler kullanarak hem yerel hem de uzak dosyaları şeffaf bir şekilde işlemek
- `channel.fromPath()` ve `channel.fromFilePairs()` ile dosya işlemeyi otomatikleştirmek için kanalları kullanmak
- Dize manipülasyonu kullanarak dosya adlarından metadata çıkarmak ve yapılandırmak
- Desen eşleştirme ve glob ifadeleri kullanarak ilgili dosyaları gruplamak
- Dosya işlemlerini uygun girdi işleme ile Nextflow süreçlerine entegre etmek
- Metadata odaklı dizin yapıları kullanarak süreç çıktılarını düzenlemek

Bu beceriler, farklı türdeki dosya girdilerini büyük esneklikle işleyebilen iş akışları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan görevi üstlenmeden önce şunları yapmalısınız:

- [Hello Nextflow](../../hello_nextflow/) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler) rahatça kullanabiliyor olmak

---

## 0. Başlarken

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![GitHub Codespaces'te Aç](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine gidin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/working_with_files
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri gözden geçirin

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

Bu dizin, üç hastadan (A, B, C) çift uçlu dizileme verisi içermektedir.

Her hasta için, `tumor` (genellikle tümör biyopsilerinden kaynaklanan) veya `normal` (sağlıklı doku veya kandan alınan) türünde örneklerimiz var.
Kanser analiziyle tanışık değilseniz, bunun karşılaştırmalı analizler yapmak için eşleştirilmiş tümör/normal örnekler kullanan bir deneysel modele karşılık geldiğini bilin.

Özellikle hasta A için, iki teknik tekrar (repeat) setimiz var.

Dizileme veri dosyaları, 'ileri okumalar' ve 'ters okumalar' olarak bilinen şeyler için tipik bir `_R1_` ve `_R2_` kuralıyla adlandırılmıştır.

_Bu deneysel tasarıma aşina değilseniz endişelenmeyin, bu eğitimi anlamak için kritik değil._

#### Görevi gözden geçirin

Göreviniz şunları yapacak bir Nextflow iş akışı yazmaktır:

1. Nextflow'un dosya işleme metodlarını kullanarak girdi dosyalarını **yüklemek**
2. Dosya adı yapısından metadata (hasta kimliği, tekrar, örnek türü) **çıkarmak**
3. `channel.fromFilePairs()` kullanarak eşleştirilmiş dosyaları (R1/R2) **gruplamak**
4. Sağlanan analiz modülü ile dosyaları **işlemek**
5. Çıktıları çıkarılan metadata'ya dayalı bir dizin yapısına **düzenlemek**

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun hedefini ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Temel dosya işlemleri

### 1.1. `.class` ile bir nesnenin türünü belirleyin

`main.nf` iş akışı dosyasına bir göz atın:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Bir dize yolundan bir Path nesnesi oluştur
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Bu, iş akışında tek bir dosya yoluna atıfta bulunan ve ardından onu sınıfıyla birlikte konsola yazdıran bir mini iş akışıdır (herhangi bir süreç olmadan).

??? info "`.class` nedir?"

    Nextflow'da `.class`, hangi tür nesneyle çalıştığımızı söyler. Bu, "bu ne tür bir şey?" diye sormak gibidir; bir dize mi, bir sayı mı, bir dosya mı yoksa başka bir şey mi olduğunu öğrenmek için.
    Bu, sonraki bölümlerde düz bir dize ile bir Path nesnesi arasındaki farkı göstermemize yardımcı olacaktır.

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

Gördüğünüz gibi, Nextflow dize yolunu tam olarak yazdığımız gibi yazdırdı.

Bu sadece metin çıktısıdır; Nextflow henüz onunla özel bir şey yapmadı.
Ayrıca Nextflow açısından bunun sadece bir dize (`java.lang.String` sınıfından) olduğunu doğruladık.
Bu mantıklı, çünkü henüz Nextflow'a bunun bir dosyaya karşılık geldiğini söylemedik.

### 1.2. file() ile bir Path nesnesi oluşturun

Nextflow'a dosyaları nasıl işleyeceğini yol dizelerinden [Path nesneleri](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) oluşturarak söyleyebiliriz.

İş akışımızda, `data/patientA_rep1_normal_R1_001.fastq.gz` dize yolunu `file()` metodunu kullanarak bir Path nesnesine dönüştürebiliriz; bu, dosya özelliklerine ve işlemlerine erişim sağlar.

`main.nf` dosyasını aşağıdaki gibi dizeyi `file()` ile sarmalayacak şekilde düzenleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Bir dize yolundan bir Path nesnesi oluştur
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

Bu sefer, girdi olarak sağladığımız göreceli yol yerine tam mutlak yolu görüyorsunuz.

Nextflow dizimizi bir Path nesnesine dönüştürdü ve onu sistemdeki gerçek dosya konumuna çözdü.
Dosya yolu artık `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz` gibi mutlak olacaktır.

Ayrıca Path nesnesi sınıfının `sun.nio.fs.UnixPath` olduğuna dikkat edin: bu, Nextflow'un yerel dosyaları temsil etme yöntemidir.
Daha sonra göreceğimiz gibi, uzak dosyalar farklı sınıf adlarına sahip olacaktır (HTTP dosyaları için `nextflow.file.http.XPath` gibi), ancak hepsi tamamen aynı şekilde çalışır ve iş akışlarınızda aynı şekilde kullanılabilir.

!!! tip

    **Temel fark:**

    - **Yol dizesi**: Nextflow'un karakter olarak ele aldığı sadece metin
    - **Path nesnesi**: Nextflow'un çalışabileceği akıllı bir dosya referansı

    Bunu şöyle düşünün: bir yol dizesi, kağıda bir adres yazmak gibidir, bir Path nesnesi ise adresi oraya nasıl gidileceğini bilen ve yolculuk hakkında size ayrıntılar verebilen bir GPS cihazına yüklenmiş olmak gibidir.

### 1.3. Dosya özniteliklerine erişin

Bu neden yararlıdır? Şimdi Nextflow `myFile`'ın bir Path nesnesi olduğunu ve sadece bir dize olmadığını anladığına göre, Path nesnesinin çeşitli özniteliklerine erişebiliriz.

İş akışımızı yerleşik dosya özniteliklerini yazdıracak şekilde güncelleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Bir dize yolundan bir Path nesnesi oluştur
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

Yukarıda konsola yazdırılan çeşitli dosya özniteliklerini görüyorsunuz.

### 1.4. Dosyayı bir sürece besleyin

Dizeler ve Path nesneleri arasındaki fark, süreçlerle gerçek iş akışları oluşturmaya başladığınızda kritik hale gelir.
Şimdiye kadar Nextflow'un artık girdi dosyamızı bir dosya olarak ele aldığını doğruladık, ancak o dosya üzerinde bir süreçte gerçekten bir şey çalıştırıp çalıştıramayacağımızı görelim.

#### 1.4.1. Süreci içe aktarın ve kodu inceleyin

Size, bir dosya girdisi alan ve içinde kaç satır olduğunu sayan `COUNT_LINES` adlı önceden yazılmış bir süreç modülü sağlıyoruz.

Süreci iş akışında kullanmak için, workflow bloğundan önce bir include ifadesi eklemeniz yeterlidir:

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

Kodunu incelemek için modül dosyasını açabilirsiniz:

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

Gördüğünüz gibi, dosyayı açan ve kaç satır içerdiğini sayan oldukça basit bir betiktir.

??? info "`debug true` ne yapar?"

    Süreç tanımındaki `debug true` yönergesi, Nextflow'un betiğinizden gelen çıktıyı (satır sayısı "40" gibi) doğrudan yürütme günlüğünde yazdırmasına neden olur.
    Bu olmadan, yalnızca süreç yürütme durumunu görürsünüz ancak betiğinizden gelen gerçek çıktıyı görmezsiniz.

    Nextflow süreçlerinde hata ayıklama hakkında daha fazla bilgi için [Nextflow İş Akışlarında Hata Ayıklama](debugging.md) yan görevine bakın.

#### 1.4.2. `COUNT_LINES`'a bir çağrı ekleyin

Artık süreç iş akışı için kullanılabilir olduğuna göre, girdi dosyası üzerinde çalıştırmak için `COUNT_LINES` sürecine bir çağrı ekleyebiliriz.

İş akışında aşağıdaki düzenlemeleri yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
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
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Ve şimdi iş akışını çalıştırın:

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

Bu, dosya üzerinde bir süreç içinde uygun şekilde işlem yapabildiğimizi gösterir.

Özellikle, Nextflow aşağıdaki işlemleri başarıyla gerçekleştirdi:

- Dosyayı çalışma dizinine hazırladı
- .gz dosyasını açtı
- Satırları saydı (bu durumda 40 satır)
- Hatasız tamamlandı

Bu sorunsuz işlemin anahtarı, Nextflow'a girdimizin bir dosya olduğunu ve bu şekilde ele alınması gerektiğini açıkça söylememizdir.

### 1.5. Temel dosya girdi hatalarını giderin

Bu genellikle Nextflow'a yeni başlayanları yanıltır, bu yüzden yanlış yaptığınızda ne olduğuna bakmak için birkaç dakika ayıralım.

Dosya işlemeyi yanlış yapabileceğiniz iki ana yer vardır: iş akışı düzeyinde ve süreç düzeyinde.

#### 1.5.1. İş akışı düzeyinde hata

İş akışı bloğunda girdiyi belirtirken dosyayı bir dize olarak ele almaya geri dönersek ne olduğunu görelim.

İş akışında aşağıdaki düzenlemeleri yapın, yola özgü yazdırma ifadelerini yorumlamayı unutmayın:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Dosya özniteliklerini yazdır
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
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

Ve şimdi iş akışını çalıştırın:

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

Önemli kısım bu:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Bir `path` girdisi belirttiğinizde, Nextflow sadece dizeler değil, gerçek dosya referansları geçirdiğinizi doğrular.
Bu hata, `'data/patientA_rep1_normal_R1_001.fastq.gz'`'nin geçerli bir yol değeri olmadığını söylüyor çünkü bu bir dize, bir Path nesnesi değil.

Nextflow sorunu hemen tespit etti ve süreci başlatmadan önce durdu.

#### 1.5.2. Süreç düzeyinde hata

Nextflow'a girdinin bir dosya olarak ele alınmasını istediğimizi belirtmeyi unutabileceğimiz diğer yer, süreç tanımındadır.

!!! warning "1.5.1'deki iş akışı hatasını koruyun"

    Bu testin doğru çalışması için, iş akışını bozuk durumunda tutun (`file()` yerine düz dize kullanarak).
    Süreçteki `val` ile birleştirildiğinde, bu aşağıda gösterilen hatayı üretir.

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

Ve şimdi iş akışını tekrar çalıştırın:

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

Bu, yukarıda belirtildiği gibi süreç hata ayıklama bilgisi çıktısı verecek şekilde ayarlandığı için hata hakkında birçok ayrıntı gösterir.

Bunlar en ilgili bölümlerdir:

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

Bu, sistemin dosyayı bulamadığını söylüyor; ancak yolu kontrol ederseniz, o konumda o adda bir dosya var.

Bunu çalıştırdığımızda, Nextflow dize değerini betiğe geçirdi, ancak gerçek dosyayı çalışma dizinine _hazırlamadı_.
Bu nedenle süreç göreceli dizeyi, `data/patientA_rep1_normal_R1_001.fastq.gz`'yi kullanmaya çalıştı, ancak bu dosya süreç çalışma dizininde mevcut değil.

Bu iki örnek birlikte, Nextflow'a bir girdinin dosya olarak işlenmesi gerektiğini söylemenin ne kadar önemli olduğunu gösterir.

!!! note

    Bir sonraki bölüme geçmeden önce her iki kasıtlı hatayı da geri alıp düzelttiğinizden emin olun.

### Özet

- Yol dizeleri vs Path nesneleri: Dizeler sadece metindir, Path nesneleri akıllı dosya referanslarıdır
- `file()` metodu bir yol dizesini Nextflow'un çalışabileceği bir Path nesnesine dönüştürür
- `name`, `simpleName`, `extension` ve `parent` gibi dosya özelliklerine [dosya özniteliklerini kullanarak](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes) erişebilirsiniz
- Dizeler yerine Path nesneleri kullanmak, Nextflow'un iş akışınızdaki dosyaları düzgün bir şekilde yönetmesine olanak tanır
- Süreç Girdi Sonuçları: Uygun dosya işleme, dosyaların süreçler tarafından kullanılmak üzere doğru şekilde hazırlanmasını ve erişilebilir olmasını sağlamak için dizeler değil, Path nesneleri gerektirir.

---

## 2. Uzak dosyaları kullanma

Nextflow'un temel özelliklerinden biri, yerel dosyalar (aynı makinede) ile internet üzerinden erişilebilen uzak dosyalar arasında sorunsuz bir şekilde geçiş yapabilme yeteneğidir.

Doğru yapıyorsanız, farklı konumlardan gelen dosyaları barındırmak için iş akışınızın mantığını asla değiştirmeniz gerekmemelidir.
Tek yapmanız gereken, iş akışına sağlarken dosya yolunda uygun öneki belirtmektir.

Örneğin, `/path/to/data` öneki yoktur, bu da 'normal' bir yerel dosya yolu olduğunu gösterir, oysa `s3://path/to/data`, Amazon'un S3 nesne depolamasında bulunduğunu gösteren `s3://` önekini içerir.

Birçok farklı protokol desteklenir:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Bunlardan herhangi birini kullanmak için, dizede ilgili öneki belirtmeniz yeterlidir; bu durumda teknik olarak dosya yolu yerine Uniform Resource Identifier (URI) olarak adlandırılır.
Nextflow kimlik doğrulamayı ve dosyaları doğru yere hazırlamayı, indirmeyi veya yüklemeyi ve beklediğiniz diğer tüm dosya işlemlerini yönetecektir.

Bu sistemin temel gücü, ortamlar arasında geçiş yapmamızı herhangi bir boru hattı mantığını değiştirmeden sağlamasıdır.
Örneğin, URI'yi değiştirerek uzak depolamada bulunan tam ölçekli bir test setine geçmeden önce küçük, yerel bir test setiyle geliştirebilirsiniz.

### 2.1. İnternetten bir dosya kullanın

Bunu, iş akışımıza sağladığımız yerel yolu, Github'da depolanan aynı verinin bir kopyasına işaret eden bir HTTPS yolu ile değiştirerek test edelim.

!!! warning

    Bu yalnızca aktif bir internet bağlantınız varsa çalışacaktır.

`main.nf`'yi tekrar açın ve girdi yolunu aşağıdaki gibi değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // İnternetten uzak bir dosya kullanma
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
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

Çalışıyor! Konsol çıktısında çok az şeyin değiştiğini görebilirsiniz.

Konsol çıktısındaki tek fark, yol nesnesi sınıfının artık `nextflow.file.http.XPath` olmasıdır, oysa yerel yol için sınıf `sun.nio.fs.UnixPath` idi.
Bu sınıfları hatırlamanıza gerek yok; bunu sadece Nextflow'un farklı konumları uygun şekilde tanımladığını ve işlediğini göstermek için belirtiyoruz.

Perde arkasında, Nextflow dosyayı çalışma dizini içinde bulunan bir hazırlama dizinine indirdi.
Bu hazırlanmış dosya daha sonra yerel bir dosya olarak ele alınabilir ve ilgili süreç dizinine sembolik bağlantı yapılabilir.

Bunun burada gerçekleştiğini, sürecin hash değerinde bulunan çalışma dizininin içeriğine bakarak doğrulayabilirsiniz.

??? abstract "Çalışma dizini içeriği"

    Süreç hash'i `8a/2ab7ca` ise, çalışma dizinini keşfedebilirsiniz:

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

Daha büyük dosyalar için, indirme adımının yerel dosyalarda çalıştırmaya kıyasla ekstra zaman alacağını unutmayın.
Ancak, Nextflow gereksiz indirmeleri önlemek için zaten hazırlanmış bir kopyaya sahip olup olmadığını kontrol eder.
Bu nedenle aynı dosya üzerinde tekrar çalıştırırsanız ve hazırlanmış dosyayı silmediyseniz, Nextflow hazırlanmış kopyayı kullanacaktır.

Bu, Nextflow kullanarak yerel ve uzak veriler arasında geçiş yapmanın ne kadar kolay olduğunu gösterir; bu, Nextflow'un temel bir özelliğidir.

!!! note

    Bu ilkenin önemli bir istisnası, HTTPS ile glob desenlerini veya dizin yollarını kullanamamanızdır çünkü HTTPS birden fazla dosyayı listeleyemez, bu nedenle tam dosya URL'lerini belirtmelisiniz.
    Ancak, blob depolama (`s3://`, `az://`, `gs://`) gibi diğer depolama protokolleri hem glob'ları hem de dizin yollarını kullanabilir.

    Bulut depolama ile glob desenlerini nasıl kullanabileceğiniz aşağıda açıklanmıştır:

    ```groovy title="Bulut depolama örnekleri (bu ortamda çalıştırılamaz)"
    // Glob desenleriyle S3 - birden fazla dosyayla eşleşir
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Glob desenleriyle Azure Blob Storage
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Glob desenleriyle Google Cloud Storage
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Glob'larla pratikte nasıl çalışılacağını bir sonraki bölümde göstereceğiz.

### 2.2. Yerel dosyaya geri dönün

Bu yan görevin geri kalanında yerel örnek dosyalarımızı kullanmaya devam edeceğiz, bu yüzden iş akışı girdisini orijinal dosyaya geri çevirelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Özet

- Uzak verilere bir URI (HTTP, FTP, S3, Azure, Google Cloud) kullanılarak erişilir
- Nextflow, bu yollar süreçlere beslendiği sürece verileri otomatik olarak doğru yere indirecek ve hazırlayacaktır
- Uzak dosyaları indirmek veya yüklemek için mantık yazmayın!
- Yerel ve uzak dosyalar farklı nesne türleri üretir ancak aynı şekilde çalışır
- **Önemli**: HTTP/HTTPS yalnızca tek dosyalarla çalışır (glob desenleri yok)
- Bulut depolama (S3, Azure, GCS) hem tek dosyaları hem de glob desenlerini destekler
- Kod mantığını değiştirmeden yerel ve uzak veri kaynakları arasında sorunsuz bir şekilde geçiş yapabilirsiniz (protokol gerekli işlemlerinizi desteklediği sürece)

---

## 3. `fromPath()` kanal fabrikasını kullanma

Şimdiye kadar bir seferde tek bir dosyayla çalışıyorduk, ancak Nextflow'da genellikle işlemek için birden fazla girdi dosyası içeren bir girdi kanalı oluşturmak isteyeceğiz.

Bunu yapmanın naif bir yolu, `file()` metodunu [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) ile şu şekilde birleştirmek olabilir:

```groovy title="Sözdizimi örneği"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Bu çalışır, ancak hantal.

!!! tip "`file()` vs `channel.fromPath()` ne zaman kullanılır"

    - Doğrudan manipülasyon için tek bir Path nesnesine ihtiyacınız olduğunda `file()` kullanın (bir dosyanın var olup olmadığını kontrol etme, özniteliklerini okuma veya tek bir süreç çağrısına geçirme)
    - Birden fazla dosya tutabilen bir kanala ihtiyacınız olduğunda, özellikle glob desenleriyle veya dosyalar birden fazla süreçten geçecekse `channel.fromPath()` kullanın

Burası [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath)'in devreye girdiği yerdir: bir veya daha fazla statik dosya dizesinden ve glob desenlerinden bir kanal oluşturmak için ihtiyacımız olan tüm işlevselliği bir araya getiren kullanışlı bir kanal fabrikası.

### 3.1. Kanal fabrikasını ekleyin

İş akışımızı `channel.fromPath` kullanacak şekilde güncelleyelim.

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // channel.fromPath ile dosyaları yükle
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Dosya özniteliklerini yazdır
        /* Şimdilik bunları yorumlayın, geri döneceğiz!
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
        // Bir dize yolundan bir Path nesnesi oluştur
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Dosya özniteliklerini yazdır
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Dosyadaki satırları say
        COUNT_LINES(myFile)
    ```

Ayrıca şimdilik öznitelikleri yazdıran kodu yorumladık ve bunun yerine sadece dosya adını yazdırmak için bir `.view` ifadesi ekledik.

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

Gördüğünüz gibi, dosya yolu kanalda bir `Path` türü nesnesi olarak yükleniyor.
Bu, `file()`'ın yapacağına benzer, ancak şimdi istersek daha fazla dosya yükleyebileceğimiz bir kanalımız var.

`channel.fromPath()` kullanmak, bir dosya listesiyle doldurulmuş yeni bir kanal oluşturmanın kullanışlı bir yoludur.

### 3.2. Kanaldaki dosyaların özniteliklerini görüntüleyin

Kanal fabrikasını kullanmadaki ilk geçişimizde, kodu basitleştirdik ve sadece dosya adını yazdırdık.

Tam dosya özniteliklerini yazdırmaya geri dönelim:

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

Ayrıca dosya işlemenin kanal tabanlı yaklaşımımızla hala doğru çalıştığını doğrulamak için `COUNT_LINES` süreç çağrısını yeniden etkinleştiriyoruz.

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

Ve işte, öncekiyle aynı sonuçlar ama şimdi dosyayı bir kanalda tutuyoruz, böylece daha fazla ekleyebiliriz.

### 3.3. Birden fazla dosyayı eşleştirmek için glob kullanma

Kanala daha fazla dosya yüklemenin birkaç yolu var.
Burada size glob desenlerini nasıl kullanacağınızı göstereceğiz; bunlar, joker karakterlere dayalı olarak dosya ve dizin adlarını eşleştirmenin ve almanın kullanışlı bir yoludur.
Bu desenleri eşleştirme işlemine "globbing" veya "dosya adı genişletme" denir.

!!! note

    Daha önce belirtildiği gibi, Nextflow çoğu durumda girdi ve çıktı dosyalarını yönetmek için globbing'i destekler, HTTPS dosya yolları hariç çünkü HTTPS birden fazla dosyayı listeleyemez.

Diyelim ki belirli bir hasta, `patientA` ile ilişkili bir dosya çiftindeki her iki dosyayı da almak istiyoruz:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Dosya adları arasındaki tek fark tekrar numarası olduğundan, _yani_ `R`'den sonraki sayı, joker karakter `*`'ı aşağıdaki gibi sayının yerine kullanabiliriz:

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

Nextflow bunun bir glob deseni olduğunu otomatik olarak tanıyacak ve uygun şekilde işleyecektir.

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

Gördüğünüz gibi, artık kanalımızda iki Path nesnesi var; bu, Nextflow'un dosya adı genişletmesini doğru yaptığını ve her iki dosyayı da beklendiği gibi yüklediğini ve işlediğini gösterir.

Bu yöntemi kullanarak, glob desenini değiştirerek istediğimiz kadar çok veya az dosya alabiliriz. Daha cömert hale getirirsek, örneğin dosya adlarının tüm değişken kısımlarını `*` ile değiştirerek (_örn._ `data/patient*_rep*_*_R*_001.fastq.gz`) `data` dizinindeki tüm örnek dosyaları alabiliriz.

### Özet

- `channel.fromPath()` bir desenle eşleşen dosyalardan bir kanal oluşturur
- Her dosya kanalda ayrı bir öğe olarak yayınlanır
- Birden fazla dosyayı eşleştirmek için bir glob deseni kullanabiliriz
- Dosyalar otomatik olarak tam özniteliklere sahip Path nesnelerine dönüştürülür
- `.view()` metodu kanal içeriğinin incelenmesine olanak tanır

---

## 4. Dosya adlarından temel metadata çıkarma

Çoğu bilimsel alanda, verileri içeren dosyaların adlarında kodlanmış metadata'ya sahip olmak çok yaygındır.
Örneğin, biyoinformatikte, dizileme verisi içeren dosyalar genellikle örnek, koşul, tekrar ve okuma numarası hakkında bilgi kodlayan bir şekilde adlandırılır.

Dosya adları tutarlı bir kurala göre oluşturulmuşsa, bu metadata'yı standart bir şekilde çıkarabilir ve analiziniz sırasında kullanabilirsiniz.
Bu büyük bir 'eğer' tabii ki ve dosya adı yapısına güvendiğinizde çok dikkatli olmalısınız; ancak gerçek şu ki bu yaklaşım çok yaygın olarak kullanılıyor, bu yüzden Nextflow'da nasıl yapıldığına bir göz atalım.

Örnek verilerimiz söz konusu olduğunda, dosya adlarının tutarlı bir şekilde yapılandırılmış metadata içerdiğini biliyoruz.
Örneğin, `patientA_rep1_normal_R2_001` dosya adı şunları kodlar:

- hasta kimliği: `patientA`
- tekrar kimliği: `rep1`
- örnek türü: `normal` (`tumor`'un aksine)
- okuma seti: `R1` (`R2`'nin aksine)

İş akışımızı bu bilgiyi üç adımda alacak şekilde değiştireceğiz:

1. Metadata içeren dosyanın `simpleName`'ini alın
2. `tokenize()` adlı bir metot kullanarak metadata'yı ayırın
3. Metadata'yı düzenlemek için bir map kullanın

!!! warning

    Hasta adları veya diğer tanımlayıcı özellikler gibi hassas bilgileri asla dosya adlarına kodlamamalısınız, çünkü bu hasta gizliliğini veya diğer ilgili güvenlik kısıtlamalarını tehlikeye atabilir.

### 4.1. `simpleName`'i alın

`simpleName`, yolu ve uzantısı çıkarılmış dosya adına karşılık gelen bir dosya özniteliğidir.

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

Kanaldaki her öğe artık `simpleName` ve orijinal dosya nesnesini içeren bir demet.

### 4.2. `simplename`'den metadata'yı çıkarın

Bu noktada, istediğimiz metadata `simplename`'e gömülü, ancak tek tek öğelere doğrudan erişemiyoruz.
Bu yüzden `simplename`'i bileşenlerine ayırmamız gerekiyor.
Neyse ki, bu bileşenler orijinal dosya adında basitçe alt çizgilerle ayrılmış, bu yüzden bu görev için mükemmel olan `tokenize()` adlı yaygın bir Nextflow metodunu uygulayabiliriz.

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

Şimdi kanalımızdaki her öğe için demet, metadata listesini (_örn._ `[patientA, rep1, normal, R1, 001]`) ve orijinal dosya nesnesini içeriyor.

Bu harika!
Hasta bilgilerimizi tek bir dizeden bir dize listesine ayırdık.
Artık hasta bilgisinin her bir parçasını ayrı ayrı işleyebiliriz.

### 4.3. Metadata'yı düzenlemek için bir map kullanın

Metadata'mız şu anda sadece düz bir liste.
Kullanımı yeterince kolay ama okunması zor.

```console
[patientA, rep1, normal, R1, 001]
```

İndeks 3'teki öğe nedir? Metadata yapısının orijinal açıklamasına geri dönmeden söyleyebilir misiniz?

Bu, her öğenin bir dizi anahtar ve bunlarla ilişkili değerlere sahip olduğu bir anahtar-değer deposu kullanmak için harika bir fırsattır, böylece karşılık gelen değeri almak için her anahtara kolayca başvurabilirsiniz.

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

Düz listemizi şimdi bir map'e dönüştürelim.
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

- **Destructuring ataması**: `def (patient, replicate, type, readNum) = ...` tokenize edilmiş değerleri tek satırda adlandırılmış değişkenlere çıkarır
- **Map literal sözdizimi**: `[id: patient, replicate: ...]` her anahtarın (örneğin `id`) bir değerle (örneğin `patient`) ilişkilendirildiği bir map oluşturur
- **İç içe yapı**: Dış liste `[..., myFile]` metadata map'ini orijinal dosya nesnesiyle eşleştirir

Ayrıca gereksiz olan bazı karakterleri kaldırmak için `replace()` adlı bir dize değiştirme metodu kullanarak birkaç metadata dizesini basitleştirdik (_örn._ tekrar kimliklerinden sadece numarayı tutmak için `replicate.replace('rep', '')`).

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

Şimdi metadata düzgün bir şekilde etiketlenmiş (_örn._ `[id:patientA, replicate:1, type:normal, readNum:2]`) bu yüzden neyin ne olduğunu söylemek çok daha kolay.

Ayrıca iş akışında metadata öğelerini gerçekten kullanmak çok daha kolay olacak ve kodumuzun okunmasını ve bakımını kolaylaştıracak.

### Özet

- Nextflow'da dosya adlarını tam bir programlama dilinin gücüyle işleyebiliriz
- İlgili bilgileri çıkarmak için dosya adlarını dizeler olarak ele alabiliriz
- `tokenize()` ve `replace()` gibi metodların kullanımı, dosya adındaki dizeleri manipüle etmemize olanak tanır
- `.map()` işlemi, yapıyı korurken kanal öğelerini dönüştürür
- Yapılandırılmış metadata (map'ler) kodu konumsal listelerden daha okunabilir ve bakımı kolay hale getirir

Sırada, eşleştirilmiş veri dosyalarını nasıl işleyeceğimize bakacağız.

---

## 5. Eşleştirilmiş veri dosyalarını işleme

Birçok deneysel tasarım, açıkça eşleştirilmiş bir şekilde işlenmekten fayda sağlayan eşleştirilmiş veri dosyaları üretir.
Örneğin, biyoinformatikte, dizileme verisi genellikle eşleştirilmiş okumalar şeklinde üretilir; bu, aynı DNA parçasından kaynaklanan dizi dizelerini ifade eder (genellikle zıt uçlardan okundukları için 'ileri' ve 'ters' olarak adlandırılır).

Bu, R1 ve R2'nin iki okuma setine atıfta bulunduğu örnek verilerimizin durumudur.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow, bunun gibi eşleştirilmiş dosyalarla çalışmak için `channel.fromFilePairs()` adlı özel bir kanal fabrikası sağlar; bu, dosyaları paylaşılan bir adlandırma desenine göre otomatik olarak gruplar. Bu, eşleştirilmiş dosyaları daha az çabayla daha sıkı bir şekilde ilişkilendirmenize olanak tanır.

İş akışımızı bundan yararlanacak şekilde değiştireceğiz.
İki adım alacak:

1. Kanal fabrikasını `channel.fromFilePairs()`'e geçirin
2. Metadata'yı çıkarın ve eşleyin

### 5.1. Kanal fabrikasını `channel.fromFilePairs()`'e geçirin

`channel.fromFilePairs` kullanmak için, Nextflow'un bir çiftteki iki üyeyi tanımlamak için kullanması gereken deseni belirtmemiz gerekir.

Örnek verilerimize geri dönersek, adlandırma desenini şu şekilde resmileştirebiliriz:

```console
data/patientA_rep1_normal_R{1,2}_001.fastq.gz
```

Bu, daha önce kullandığımız glob desenine benzer, ancak bu özellikle çiftin iki üyesini tanımlayan alt dizeleri (R'den hemen sonra gelen `1` veya `2`) numaralandırır.

İş akışı `main.nf`'yi buna göre güncelleyelim:

=== "Sonra"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // channel.fromFilePairs ile dosyaları yükle
        ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        /* Şimdilik eşlemeyi yorumlayın, geri döneceğiz!
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

Kanal fabrikasını değiştirdik ve dosya eşleştirme desenini uyarladık ve bunu yaparken map işlemini yorumladık.
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

Hata mesajının ilgili kısmı burada:

```console
Not a valid path value: 'patientA_rep1_normal_R'
```

Bunun nedeni kanal fabrikasını değiştirmiş olmamızdır.
Şimdiye kadar, orijinal girdi kanalı yalnızca dosya yollarını içeriyordu.
Yaptığımız tüm metadata manipülasyonu aslında kanal içeriğini etkilemedi.

Artık `.fromFilePairs` kanal fabrikasını kullandığımıza göre, ortaya çıkan kanalın içeriği farklı.
Yalnızca bir kanal öğesi görüyoruz; bu, iki öğe içeren bir demetten oluşuyor: iki dosya tarafından paylaşılan `simpleName`'in bir kısmı (tanımlayıcı olarak hizmet eder) ve iki dosya nesnesini içeren bir demet, `id, [ file1, file2 ]` formatında.

Bu harika, çünkü Nextflow paylaşılan öneki inceleyerek ve bunu hasta tanımlayıcısı olarak kullanarak hasta adını çıkarma konusunda zor işi yaptı.

Ancak, mevcut iş akışımızı bozuyor.
Süreci değiştirmeden `COUNT_LINES`'ı hala aynı şekilde çalıştırmak isteseydik, dosya yollarını çıkarmak için bir eşleme işlemi uygulamamız gerekirdi.
Ancak bunu yapmayacağız, çünkü nihai hedefimiz dosya çiftlerini uygun şekilde işleyen farklı bir süreç olan `ANALYZE_READS`'i kullanmak.

Bu yüzden `COUNT_LINES` çağrısını yorumlayalım (veya silelim) ve devam edelim.

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

Yaşasın, bu sefer iş akışı başarılı oldu!

Ancak, `id` alanından metadata'nın geri kalanını çıkarmamız gerekiyor.

### 5.2. Dosya çiftlerinden metadata çıkarın ve düzenleyin

Önceki `map` işlemimiz veri yapısıyla eşleşmediği için çalışmayacak, ancak çalışması için değiştirebiliriz.

Tanımlayıcı olarak `fromFilePairs()`'in kullandığı dizede gerçek hasta tanımlayıcısına zaten erişimimiz var, bu yüzden daha önce yaptığımız gibi Path nesnesinden `simpleName`'i almadan metadata'yı çıkarmak için bunu kullanabiliriz.

İş akışındaki map işleminin yorumunu kaldırın ve aşağıdaki düzenlemeleri yapın:

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
        /* Şimdilik eşlemeyi yorumlayın, geri döneceğiz!
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

Bu sefer map, sadece `myFile` yerine `id, files` ile başlıyor ve `tokenize()`, `myFile.simpleName` yerine `id`'ye uygulanıyor.

Ayrıca `tokenize()` satırından `readNum`'u çıkardığımıza dikkat edin; özellikle adlandırmadığımız (soldan başlayarak) tüm alt dizeler sessizce bırakılacaktır.
Bunu yapabiliriz çünkü eşleştirilmiş dosyalar artık sıkı bir şekilde ilişkilendirilmiş, bu yüzden artık metadata map'inde `readNum`'a ihtiyacımız yok.

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

Ve işte burada: çıktı demetinin ilk konumunda metadata map'i (`[id:patientA, replicate:1, type:normal]`), ardından amaçlandığı gibi eşleştirilmiş dosyaların demeti var.

Tabii ki, bu yalnızca o belirli dosya çiftini alacak ve işleyecektir.
Birden fazla çifti işlemeyi denemek isterseniz, girdi desenine joker karakterler eklemeyi deneyebilir ve ne olduğunu görebilirsiniz.
Örneğin, `data/patientA_rep1_*_R{1,2}_001.fastq.gz` kullanmayı deneyin

### Özet

- [`channel.fromFilePairs()` otomatik olarak ilgili dosyaları bulur ve eşleştirir](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- Bu, boru hattınızda çift uçlu okumaları işlemeyi basitleştirir
- Eşleştirilmiş dosyalar `[id, [file1, file2]]` demetleri olarak gruplandırılabilir
- Metadata çıkarma, tek tek dosyalar yerine eşleştirilmiş dosya kimliğinden yapılabilir

---

## 6. Süreçlerde dosya işlemlerini kullanma

Şimdi tüm bunları, Nextflow sürecinde dosya işlemlerinin nasıl kullanılacağını pekiştirmek için basit bir süreçte bir araya getirelim.

Size, bir metadata demeti ve bir çift girdi dosyası alan ve bunları analiz eden `ANALYZE_READS` adlı önceden yazılmış bir süreç modülü sağlıyoruz.
Bunun dizi hizalama veya varyant çağırma veya bu veri türü için mantıklı olan başka herhangi bir adım yaptığını hayal edebiliriz.

Başlayalım.

### 6.1. Süreci içe aktarın ve kodu inceleyin

Bu süreci iş akışında kullanmak için, workflow bloğundan önce bir modül include ifadesi eklememiz yeterlidir.

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

Kodunu incelemek için modül dosyasını açabilirsiniz:

```groovy title="modules/analyze_reads.nf - süreç örneği" linenums="1"
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

!!! note

    `tag` ve `publishDir` yönergeleri, dize enterpolasyonu (`"${...}"`) yerine closure sözdizimi (`{ ... }`) kullanır.
    Bunun nedeni, bu yönergelerin çalışma zamanına kadar kullanılamayan girdi değişkenlerine (`meta`) başvurmasıdır.
    Closure sözdizimi, süreç gerçekten çalışana kadar değerlendirmeyi erteler.

!!! note

    Metadata map'imizi geleneksel olarak `meta` olarak adlandırıyoruz.
    Meta map'lere daha derin bir dalış için, [Metadata ve meta map'ler](./metadata.md) yan görevine bakın.

### 6.2. İş akışında süreci çağırın

Artık süreç iş akışı için kullanılabilir olduğuna göre, çalıştırmak için `ANALYZE_READS` sürecine bir çağrı ekleyebiliriz.

Örnek verilerimiz üzerinde çalıştırmak için iki şey yapmamız gerekecek:

1. Yeniden eşlenmiş kanala bir ad verin
2. Sürece bir çağrı ekleyin

#### 6.2.1. Yeniden eşlenmiş girdi kanalını adlandırın

Daha önce eşleme manipülasyonlarını doğrudan girdi kanalına uyguladık.
Yeniden eşlenmiş içeriği `ANALYZE_READS` sürecine beslemek için (ve bunu açık ve okunması kolay bir şekilde yapmak için) `ch_samples` adlı yeni bir kanal oluşturmak istiyoruz.

Bunu [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set) operatörünü kullanarak yapabiliriz.

Ana iş akışında, `.view()` operatörünü `.set { ch_samples }` ile değiştirin ve kanala adıyla başvurabileceğimizi test eden bir satır ekleyin.

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

        // Geçici: ch_samples'a göz atın
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

Bu, artık kanala adıyla başvurabileceğimizi doğrular.

#### 6.2.2. Veri üzerinde süreci çağırın

Şimdi `ch_samples` kanalı üzerinde `ANALYZE_READS` sürecini gerçekten çağıralım.

Ana iş akışında, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="main.nf" linenums="23"
        // Analizi çalıştır
        ANALYZE_READS(ch_samples)
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="23"
        // Geçici: ch_samples'a göz atın
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

Bu süreç, çıktılarını bir `results` dizinine yayınlayacak şekilde ayarlanmış, bu yüzden oraya bir göz atın.

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

Süreç girdilerimizi aldı ve tasarlandığı gibi hasta metadata'sını içeren yeni bir dosya oluşturdu.
Muhteşem!

### 6.3. Çok daha fazla hasta ekleyin

Tabii ki, bu sadece tek bir hasta için tek bir dosya çiftini işliyor; bu, Nextflow ile elde etmeyi umduğunuz yüksek verim türü değil.
Muhtemelen bir seferde çok daha fazla veri işlemek isteyeceksiniz.

`channel.fromPath()`'in girdi olarak bir _glob_ kabul ettiğini unutmayın; bu, desenle eşleşen herhangi bir sayıda dosyayı kabul edebileceği anlamına gelir.
Bu nedenle tüm hastaları dahil etmek istersek, daha önce geçerken belirtildiği gibi, daha fazla hasta içerecek şekilde girdi dizesini değiştirebiliriz.

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

Boru hattını tekrar çalıştırın:

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

Belki de değil.
Daha yakından bakarsanız, bir sorunumuz var: patientA için iki tekrarımız var, ancak sadece bir çıktı dosyası!
Her seferinde çıktı dosyasının üzerine yazıyoruz.

### 6.4. Yayınlanan dosyaları benzersiz yapın

Hasta metadata'sına erişimimiz olduğundan, farklılaştırıcı metadata'yı dizin yapısına veya dosya adlarının kendisine dahil ederek yayınlanan dosyaları benzersiz yapmak için kullanabiliriz.

İş akışında aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.type}/${meta.id}/${meta.replicate}" }, mode: 'copy'
    ```

=== "Önce"

    ```groovy title="modules/analyze_reads.nf" linenums="6"
        publishDir { "results/${meta.id}" }, mode: 'copy'
    ```

Burada örnek türleri ve tekrarları hesaba katmak için ek dizin seviyeleri kullanma seçeneğini gösteriyoruz, ancak bunu dosya adı seviyesinde yapmayı da deneyebilirsiniz.

Şimdi boru hattını bir kez daha çalıştırın, ancak kendinize temiz bir çalışma alanı vermek için önce sonuçlar dizinini kaldırdığınızdan emin olun:

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

Şimdi sonuçlar dizinini kontrol edin:

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

Ve işte burada, tüm metadata'mız düzgün bir şekilde düzenlenmiş. Bu başarı!

Metadata'nızı bunun gibi bir map'e yükledikten sonra yapabileceğiniz çok daha fazla şey var:

1. Hasta özniteliklerine dayalı düzenli çıktı dizinleri oluşturun
2. Süreçlerde hasta özelliklerine dayalı kararlar alın
3. Metadata değerlerine dayalı verileri bölün, birleştirin ve yeniden birleştirin

Metadata'yı açık tutma ve veriye ekleme (dosya adlarında kodlanmış yerine) bu modeli, sağlam, bakımı kolay analiz iş akışları oluşturmayı sağlayan Nextflow'da temel bir en iyi uygulamadır.
Bu konuda daha fazla bilgiyi [Metadata ve meta map'ler](./metadata.md) yan görevinde öğrenebilirsiniz.

### Özet

- `publishDir` yönergesi, çıktıları metadata değerlerine göre düzenleyebilir
- Demetlerdeki metadata, sonuçların yapılandırılmış organizasyonunu sağlar
- Bu yaklaşım, net veri kökeni ile bakımı kolay iş akışları oluşturur
- Süreçler, metadata ve dosyaların demetlerini girdi olarak alabilir
- `tag` yönergesi, yürütme günlüklerinde süreç tanımlaması sağlar
- İş akışı yapısı, kanal oluşturmayı süreç yürütmesinden ayırır

---

## Özet

Bu yan görevde, Nextflow'da dosyalarla nasıl çalışılacağını, temel işlemlerden dosya koleksiyonlarını işlemek için daha gelişmiş tekniklere kadar öğrendiniz.

Bu teknikleri kendi çalışmanızda uygulamak, özellikle karmaşık adlandırma kurallarına sahip çok sayıda dosyayla çalışırken daha verimli ve bakımı kolay iş akışları oluşturmanızı sağlayacaktır.

### Temel desenler

1.  **Temel Dosya İşlemleri:** `file()` ile Path nesneleri oluşturduk ve ad, uzantı ve üst dizin gibi dosya özniteliklerine eriştik; dizeler ve Path nesneleri arasındaki farkı öğrendik.

    - `file()` ile bir Path nesnesi oluşturun

    ```groovy
    myFile = file('path/to/file.txt')
    ```

    - Dosya özniteliklerini alın

    ```groovy
    println myFile.name       // file.txt
    println myFile.baseName   // file
    println myFile.extension  // txt
    println myFile.parent     // path/to
    ```

2.  **Uzak Dosyaları Kullanma**: URI'ler kullanarak yerel ve uzak dosyalar arasında şeffaf bir şekilde nasıl geçiş yapılacağını öğrendik; Nextflow'un iş akışı mantığını değiştirmeden çeşitli kaynaklardan dosyaları işleme yeteneğini gösterdik.

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

3.  **`fromPath()` kanal fabrikasını kullanarak dosyaları yükleme:** `channel.fromPath()` ile dosya desenlerinden kanallar oluşturduk ve nesne türleri dahil dosya özniteliklerini görüntüledik.

    - Bir dosya deseninden kanal oluşturun

    ```groovy
     ch_files = channel.fromPath('data/*.fastq.gz')
    ```

    - Dosya özniteliklerini alın

    ```groovy
     ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }
    ```

4.  **Dosya Adlarından Hasta Metadata'sını Çıkarma:** Dosya adlarından metadata çıkarmak ve yapılandırmak için `tokenize()` ve `replace()` kullandık, bunları düzenli map'lere dönüştürdük.

    ```groovy
    def name = file.name.tokenize('_')
    def patientId = name[0]
    def replicate = name[1].replace('rep', '')
    def type = name[2]
    def readNum = name[3].replace('R', '')
    ```

5.  **channel.fromFilePairs ile basitleştirme:** İlgili dosyaları otomatik olarak eşleştirmek ve eşleştirilmiş dosya kimliklerinden metadata çıkarmak için `channel.fromFilePairs()` kullandık.

    ```groovy
    ch_pairs = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

6.  **Süreçlerde Dosya İşlemlerini Kullanma:** Dosya işlemlerini uygun girdi işleme ile Nextflow süreçlerine entegre ettik, çıktıları metadata'ya göre düzenlemek için `publishDir` kullandık.

    - Süreç girdileriyle bir meta map ilişkilendirin

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

    - Çıktıları metadata'ya göre düzenleyin

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
