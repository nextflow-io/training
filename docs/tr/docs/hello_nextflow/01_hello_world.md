# Bölüm 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=tr" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube kanalında [tüm oynatma listesine](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) bakın.

:green_book: Video transkripti [burada](./transcripts/01_hello_world.md) mevcuttur.
///

Hello Nextflow eğitim kursunun bu ilk bölümünde, temel Nextflow mantığını ve bileşenlerini göstermek için kademeli olarak oluşturacağımız çok basit, alana bağımlı olmayan bir Hello World örneğiyle konuya giriyoruz.

??? info "Hello World örneği nedir?"

    "Hello World!", bir programlama dilinin veya yazılım çerçevesinin temel söz dizimini ve yapısını göstermek için tasarlanmış minimalist bir örnektir.
    Örnek genellikle "Hello, World!" ifadesini konsol veya terminal gibi çıktı cihazına yazdırmaktan veya bir dosyaya yazmaktan oluşur.

---

## 0. Isınma: Doğrudan Hello World örneği çalıştırın

Bunu, Nextflow'a sarmadan önce ne yaptığını göstermek için terminalde doğrudan çalıştırdığımız basit bir komutla gösterelim.

!!! tip "İpucu"

    [Başlarken](00_orientation.md) sayfasında açıklandığı gibi artık `hello-nextflow/` dizininin içinde olmalısınız.

### 0.1. Terminalin merhaba demesini sağlayın

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
echo 'Hello World!'
```

??? success "Komut çıktısı"

    ```console
    Hello World!
    ```

Bu, terminalde 'Hello World' metnini çıktı olarak verir.

### 0.2. Çıktıyı bir dosyaya yazın

İş akışlarını çalıştırmak çoğunlukla dosyalardan veri okumayı ve sonuçları diğer dosyalara yazmayı içerir, bu nedenle örneği biraz daha ilgili hale getirmek için metin çıktısını bir dosyaya yazacak şekilde komutu değiştirelim.

```bash
echo 'Hello World!' > output.txt
```

??? success "Komut çıktısı"

    ```console

    ```

Bu, terminale hiçbir şey çıktı vermez.

### 0.3. Çıktıyı bulun

'Hello World' metni artık belirttiğimiz `output.txt` adlı çıktı dosyasında olmalıdır.
Dosya gezgininde veya örneğin `cat` yardımcı programını kullanarak komut satırından açabilirsiniz.

??? abstract "Dosya içerikleri"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

İşte ilk Nextflow iş akışımızla çoğaltmaya çalışacağımız şey bu.

### Özet

Artık terminalde metin çıktısı veren basit bir komutu nasıl çalıştıracağınızı ve isteğe bağlı olarak çıktıyı bir dosyaya nasıl yazdıracağınızı biliyorsunuz.

### Sırada ne var?

Bunun bir Nextflow iş akışı olarak nasıl yazıldığını öğrenin.

---

## 1. Betiği inceleyin ve çalıştırın

Size daha önce yaptığımız şeyi (Hello World!' yazmak) Nextflow ile yapan `hello-world.nf` adında tam işlevsel ama minimalist bir iş akışı betiği sağlıyoruz.

Başlamanız için, nasıl yapılandırıldığını anlamanız için iş akışı betiğini açalım.
Sonra çalıştırıp çıktılarını arayacağız.

### 1.1. Kodu inceleyin

`hello-world.nf` betiğini mevcut dizininizde bulacaksınız, bu `hello-nextflow` olmalıdır. Düzenleyici bölmesinde açın.

??? full-code "Tam kod dosyası"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' ifadesini bir dosyaya yazdırmak için echo kullan
    */
    process sayHello {

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }

    workflow {

        main:
        // bir selamlama yayınla
        sayHello()
    }
    ```

Bir Nextflow iş akışı betiği genellikle bir veya daha fazla [**process**](https://nextflow.io/docs/latest/process.html) tanımı ve [**workflow**](https://nextflow.io/docs/latest/workflow.html)'un kendisini, ayrıca daha sonra tanıtacağımız birkaç isteğe bağlı blok (burada mevcut değil) içerir.

Her **process**, iş akışındaki karşılık gelen adımın hangi işlem(ler)i gerçekleştirmesi gerektiğini açıklarken, **workflow** çeşitli adımları birbirine bağlayan veri akışı mantığını tanımlar.

Önce **process** bloğuna daha yakından bakacağız, sonra **workflow** bloğuna bakacağız.

#### 1.1.1. `process` tanımı

İlk kod bloğu bir **process**'i tanımlar.

Süreç tanımı `process` anahtar kelimesiyle başlar, ardından süreç adı ve son olarak süslü parantezlerle sınırlandırılmış süreç gövdesi gelir.
Süreç gövdesi, çalıştırılacak komutu belirten bir script bloğu içermelidir; bu, bir komut satırı terminalinde çalıştırabileceğiniz herhangi bir şey olabilir.

```groovy title="hello-world.nf" linenums="3"
/*
* 'Hello World!' ifadesini bir dosyaya yazdırmak için echo kullan
*/
process sayHello {

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Burada **çıktısını** `output.txt` adlı bir dosyaya yazan `sayHello` adında bir **process** var.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Bu, yalnızca bir `output` tanımı ve yürütülecek `script`'i içeren çok minimal bir süreç tanımıdır.

`output` tanımı, bunun bir yol olarak ele alınması gerektiğini Nextflow'a söyleyen `path` niteleyicisini içerir (hem dizin yollarını hem de dosyaları içerir).
Yaygın bir diğer niteleyici `val`'dır.

Önemli olarak, çıktı tanımı hangi çıktının oluşturulacağını _belirlemez_.
Yalnızca beklenen çıktının ne olduğunu _bildirir_, böylece Nextflow yürütme tamamlandıktan sonra onu arayabilir.
Bu, komutun başarıyla yürütüldüğünü doğrulamak ve gerekirse çıktıyı aşağı akış süreçlerine iletmek için gereklidir. Çıktı bloğunda bildirilen şeyle eşleşmeyen üretilen çıktı, aşağı akış süreçlerine iletilmez.

!!! warning "Uyarı"

    Bu örnek kırılgandır çünkü çıktı dosya adını iki ayrı yerde (script ve output blokları) sabit kodladık.
    Birini değiştirip diğerini değiştirmezsek, betik bozulur.
    Daha sonra, bu sorunu azaltmak için değişkenleri kullanmanın yollarını öğreneceksiniz.

Gerçek dünya iş akışlarında, bir süreç genellikle biraz sonra tanıtacağımız yönergeler ve girdiler gibi ek bloklar içerir.

#### 1.1.2. `workflow` tanımı

İkinci kod bloğu **workflow**'un kendisini tanımlar.
İş akışı tanımı `workflow` anahtar kelimesiyle başlar, ardından isteğe bağlı bir ad, sonra süslü parantezlerle sınırlandırılmış iş akışı gövdesi gelir.

Burada `sayHello` sürecine bir çağrı içeren `main:` bloğundan (iş akışının ana gövdesi olduğunu söyler) oluşan bir **workflow** var.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // bir selamlama yayınla
    sayHello()
}
```

Bu çok minimal bir **workflow** tanımıdır.
Gerçek dünya iş akışlarında, iş akışı genellikle **channel**'larla bağlanmış **process**'lere birden fazla çağrı içerir ve süreçler bir veya daha fazla değişken **girdi** bekler.

Değişken girdilerin nasıl ekleneceğini bu eğitim modülünün ilerleyen bölümlerinde öğreneceksiniz; ve daha fazla süreç eklemeyi ve bunları kanallarla bağlamayı bu kursun 3. Bölümünde öğreneceksiniz.

!!! tip "İpucu"

    Teknik olarak `main:` satırı bunun gibi basit iş akışları için gerekli değildir, bu nedenle sahip olmayan iş akışlarıyla karşılaşabilirsiniz.
    Ancak iş akışı düzeyinde çıktılardan yararlanmak için buna ihtiyacımız olacak, bu yüzden en başından dahil edebiliriz.

### 1.2. İş akışını çalıştırın

Koda bakmak, çalıştırmak kadar eğlenceli değil, o halde bunu pratikte deneyelim.

#### 1.2.1. İş akışını başlatın ve yürütmeyi izleyin

Terminalde aşağıdaki komutu çalıştırın:

```bash
nextflow run hello-world.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [65/7be2fa] sayHello | 1 of 1 ✔
    ```

Konsol çıktınız buna benziyorsa, tebrikler, ilk Nextflow iş akışınızı çalıştırdınız!

Buradaki en önemli çıktı, yukarıdaki çıktıda vurgulanan son satırdır:

```console
[65/7be2fa] sayHello | 1 of 1 ✔
```

Bu bize `sayHello` sürecinin bir kez başarıyla yürütüldüğünü (`1 of 1 ✔`) söyler.

Önemli olarak, bu satır aynı zamanda `sayHello` süreç çağrısının çıktısını nerede bulacağınızı da söyler.
Şimdi buna bakalım.

#### 1.2.2. `work` dizininde çıktıyı ve günlükleri bulun

Nextflow'u belirli bir dizinde ilk kez çalıştırdığınızda, yürütme sırasında oluşturulan tüm dosyaları (ve sembolik bağlantıları) yazacağı `work` adlı bir dizin oluşturur.

`work` dizini içinde, Nextflow çıktıları ve günlükleri süreç çağrısı başına düzenler.
Her süreç çağrısı için Nextflow, benzersiz olması için hash ile adlandırılmış iç içe bir alt dizin oluşturur ve gerekli tüm girdileri hazırlar (varsayılan olarak sembolik bağlantılar kullanarak), yardımcı dosyalar yazar ve sürecin herhangi bir çıktısını ve günlüklerini yazar.

Bu alt dizinin yolu, konsol çıktısında köşeli parantez içinde kısaltılmış biçimde gösterilir.
Yukarıda gösterilen çalıştırma için elde ettiğimize bakıldığında, sayHello süreci için konsol günlük satırı `[65/7be2fa]` ile başlar. Bu şu dizin yoluna karşılık gelir: `work/65/7be2fa7be2fad5e71e5f49998f795677fd68`

Orada ne olduğuna bir bakalım.

??? abstract "Dizin içerikleri"

    ```console
    work
    └── 65
        └── 7be2fad5e71e5f49998f795677fd68
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Aynı şeyi görmüyor musunuz?"

    Tam alt dizin adları sisteminizde farklı olacaktır.

    VSCode dosya gezgininde görev alt dizininin içeriğine göz atarsanız, tüm dosyaları hemen göreceksiniz.
    Ancak günlük dosyaları terminalde görünmez olarak ayarlanmıştır, bu nedenle bunları görüntülemek için `ls` veya `tree` kullanmak istiyorsanız, görünmez dosyaları görüntülemek için ilgili seçeneği ayarlamanız gerekir.

    ```bash
    tree -a work
    ```

Bakmak istediğiniz ilk şey, iş akışının gerçek çıktısı, yani `sayHello` süreci tarafından üretilen `output.txt` dosyasıdır.
Açın ve `Hello World!` selamlamasını bulacaksınız, ki minimalist iş akışımızın amacı buydu.

??? abstract "Dosya içerikleri"

    ```console title="output.txt"
    Hello World!
    ```

İşe yaradı!

Kuşkusuz, bu kadar küçük bir sonuç için çok fazla sarmalayıcı kod gibi görünebilir, ancak tüm bu sarmalayıcı kodun değeri, girdi dosyalarını okumaya ve birden fazla adımı bir araya getirmeye başladığımızda daha belirgin hale gelecektir.

Bununla birlikte, o dizindeki diğer dosyalara da bakalım. Bunlar, görev yürütmenin bir parçası olarak Nextflow tarafından üretilen yardımcı ve günlük dosyalarıdır.

- **`.command.begin`**: Süreç çağrısının yürütülmesinin başlangıcıyla ilgili meta veriler
- **`.command.err`**: Süreç çağrısı tarafından yayılan hata mesajları (`stderr`)
- **`.command.log`**: Süreç çağrısı tarafından yayılan tam günlük çıktısı
- **`.command.out`**: Süreç çağrısı tarafından yayılan normal çıktı (`stdout`)
- **`.command.run`**: Süreç çağrısını yürütmek için Nextflow tarafından çalıştırılan tam betik
- **`.command.sh`**: Süreç çağrısı tarafından gerçekte çalıştırılan komut
- **`.exitcode`**: Komuttan kaynaklanan çıkış kodu

`.command.sh` dosyası özellikle kullanışlıdır çünkü tüm defter tutma ve görev/ortam kurulumu dahil değil, Nextflow'un yürüttüğü ana komutu söyler.

??? abstract "Dosya içerikleri"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Bu, daha önce manuel olarak çalıştırdığımız şeyle eşleşir.

Bu durumda süreç komutu sabit kodlanmış olduğu için çok basittir, ancak kursun ilerleyen bölümlerinde bazı değişken enterpolasyonu içeren süreç komutları göreceksiniz.
Bu, başarısız bir çalıştırmayı sorun giderirken Nextflow'un kodu nasıl yorumladığını ve hangi komutun üretildiğini tam olarak görebilmeyi özellikle değerli kılar.

### 1.3. İş akışını tekrar çalıştırın

İş akışını birkaç kez daha çalıştırmayı deneyin, ardından `work/` altındaki görev dizinlerine bakın.

??? abstract "Dizin içerikleri"

    ```console
    work
    ├── 0f
    │   └── 52b7e07b0e274a80843fca48ed21b8
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 65
        └── 7be2fad5e71e5f49998f795677fd68
    │   │   ├── .command.begin
    │   │   ├── .command.err
    │   │   ├── .command.log
    │   │   ├── .command.out
    │   │   ├── .command.run
    │   │   ├── .command.sh
    │   │   ├── .exitcode
    │   │   └── output.txt
    │   └── e029f2e75305874a9ab263d21ebc2c
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 6c
    │   └── d4fd787e0b01b3c82e85696c297500
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── e8
        └── ab99fad46ade52905ec973ff39bb80
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Her çalıştırma için eksiksiz çıktı ve günlük dosyaları içeren yeni bir alt dizin oluşturulduğunu görüyorsunuz.
Bu, aynı iş akışını birkaç kez çalıştırmanın önceki çalıştırmaların sonuçlarını üzerine yazmayacağını gösterir.

### Özet

Basit bir Nextflow betiğini nasıl çözeceğinizi, çalıştıracağınızı ve work dizininde çıktı ve ilgili günlük dosyalarını nasıl bulacağınızı biliyorsunuz.

### Sırada ne var?

İş akışı çıktılarını daha uygun bir konuma nasıl yayınlayacağınızı öğrenin.

---

## 2. Çıktıları yayınlayın

Az önce öğrendiğiniz gibi, iş akışımız tarafından üretilen çıktı, birkaç seviye derinliğindeki bir çalışma dizininde gömülüdür.
Bu kasıtlı olarak yapılır; Nextflow bu dizini kontrol eder ve biz onunla etkileşime girmememiz gerekir.
Ancak bu, önemsediğimiz çıktıları almayı zorlaştırır.

Neyse ki Nextflow, [iş akışı düzeyinde çıktı tanımlarını](https://www.nextflow.io/docs/latest/workflow.html#workflow-outputs) kullanarak çıktıları belirlenmiş bir dizine yayınlamanın bir yolunu sağlar.

### 2.1. Temel kullanım

Bu, iki yeni kod parçası içerecektir:

1. `workflow` gövdesi içinde süreç çıktılarını bildiren bir `publish:` bloğu.
2. Mod ve konum gibi çıktı seçeneklerini belirten betiğe bir `output` bloğu.

#### 2.1.1. `sayHello` sürecinin çıktısını bildirin

İş akışı gövdesine (`main:` bloğuyla aynı türde kod öğesi) bir `publish:` bloğu eklememiz ve `sayHello()` sürecinin çıktısını listelememiz gerekir.

İş akışı betik dosyası `hello-world.nf`'de aşağıdaki kod satırlarını ekleyin:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="7-8"
    workflow {

        main:
        // bir selamlama yayınla
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // bir selamlama yayınla
        sayHello()
    }
    ```

Sürecin çıktısına basitçe `sayHello().out` yaparak başvurabileceğimizi ve ona rastgele bir ad, `first_output` atayabileceğimizi görüyorsunuz.

#### 2.1.2. Betiğe bir `output:` bloğu ekleyin

Şimdi çıktı dizini yolunun belirtileceği `output:` bloğunu eklememiz yeterli. Bu yeni bloğun betik içindeki `workflow` bloğunun **dışında** ve **altında** yer aldığını unutmayın.

İş akışı betik dosyası `hello-world.nf`'de aşağıdaki kod satırlarını ekleyin:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="17" hl_lines="11-15"
    workflow {

        main:
        // bir selamlama yayınla
        sayHello()

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '.'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="17"
    workflow {

        main:
        // bir selamlama yayınla
        sayHello()

        publish:
        first_output = sayHello.out
    }
    ```

Bunu, `workflow` bloğunda bildirilen herhangi bir süreç çıktısına belirli yollar atamak için kullanabiliriz.
Daha sonra, sofistike çıktı dizin yapıları oluşturmanın yollarını öğreneceksiniz, ancak şimdilik basitlik için minimal bir yolu sabit kodluyoruz.

#### 2.1.3. İş akışını çalıştırın

Şimdi değiştirilmiş iş akışı betiğini çalıştırın:

```bash
nextflow run hello-world.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [jovial_mayer] DSL2 - revision: 35bd3425e5

    executor >  local (1)
    [9f/48ef97] sayHello | 1 of 1 ✔
    ```

Terminal çıktısı tanıdık görünmeli. Dışarıdan hiçbir şey değişmedi.

Ancak dosya gezgininizi kontrol edin: bu sefer Nextflow `results/` adında yeni bir dizin oluşturdu.

??? abstract "Dizin içerikleri"

    ```console hl_lines="10-11 22"
    .
    ├── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── results
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/9f/48ef97f110b0dbd83635d7cbe288d2/output.txt
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── work
        ├── 65
        └── 9f
    ```

`results` dizininin içinde, az önce çalıştırdığımız komut tarafından work dizininde üretilen `output.txt`'ye sembolik bir bağlantı buluyoruz.

Bu, work alt dizinini araştırmak zorunda kalmadan çıktı dosyalarını kolayca almamızı sağlar.

### 2.2. Özel bir konum ayarlayın

Varsayılan bir konuma sahip olmak harikadır, ancak sonuçların nereye kaydedildiğini ve nasıl düzenlendiğini özelleştirmek isteyebilirsiniz.

Örneğin, çıktılarınızı alt dizinler halinde düzenlemek isteyebilirsiniz.
Bunu yapmanın en basit yolu, çıktı başına belirli çıktı yolu atamaktır.

#### 2.2.1. Çıktı yolunu değiştirin

Bir kez daha, belirli bir çıktı için yayınlama davranışını değiştirmek gerçekten basittir.
Özel bir konum ayarlamak için `path`'i buna göre düzenlemeniz yeterlidir:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="3"
    output {
        first_output {
            path '.'
        }
    }
    ```

Bu, bireysel çıktı seviyesinde ayarlandığından, ihtiyaçlarınıza uygun farklı konumlar ve alt dizinler belirleyebilirsiniz.

#### 2.2.2. İş akışını tekrar çalıştırın

Deneyelim.

```bash
nextflow run hello-world.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [8c/79499c] process > sayHello [100%] 1 of 1 ✔
    ```

Bu sefer sonuç belirtilen alt dizine yazılır.

??? abstract "Dizin içerikleri"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Önceki yürütmenin sonucunun hâlâ orada olduğunu görüyorsunuz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

İstediğiniz kadar iç içe geçme seviyesi kullanabilirsiniz.
Sonuçları düzenlemek için kullanılan dizinleri adlandırmak için süreç adını veya diğer değişkenleri kullanmak ve üst düzey çıktı dizininin varsayılan adını değiştirmek (`-o` CLI bayrağı veya config değişkeni `outputDir` tarafından kontrol edilir) de mümkündür.
Bu seçenekleri eğitimin ilerleyen bölümlerinde ele alacağız.

### 2.3. Yayınlama modunu kopyalamaya ayarlayın

Varsayılan olarak, çıktılar `work` dizininden sembolik bağlantılar olarak yayınlanır.
Bu, dosya sisteminde yalnızca tek bir dosya olduğu anlamına gelir.

Bu, birden fazla kopya saklamak istemediğiniz çok büyük dosyalarla uğraşırken harikadır.
Ancak work dizinini herhangi bir noktada silerseniz (kısa süre sonra temizleme işlemlerini ele alacağız), dosyaya erişimi kaybedersiniz.
Bu nedenle, önemli dosyaların kopyalarını güvenli bir yere kaydetmek için bir planınız olması gerekir.

Kolay bir seçenek, önemsediğiniz çıktılar için yayınlama modunu kopyalamaya geçirmektir.

#### 2.3.1. Mod yönergesini ekleyin

Bu kısım gerçekten basittir.
İlgili iş akışı düzeyinde çıktı tanımına `mode 'copy'` eklemeniz yeterlidir:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="27" hl_lines="4"
    output {
        first_output {
            path 'hello_world'
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="27"
    output {
        first_output {
            path 'hello_world'
        }
    }
    ```

Bu, söz konusu çıktı için yayınlama modunu ayarlar.

#### 2.3.2. İş akışını tekrar çalıştırın

Deneyelim.

```bash
nextflow run hello-world.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [tiny_shaw] DSL2 - revision: 757723adc1

    executor >  local (1)
    [df/521638] process > sayHello [100%] 1 of 1 ✔
    ```

Bu sefer, sonuçlara bakarsanız, dosya yalnızca bir sembolik bağlantı yerine düzgün bir kopyadır.

??? abstract "Dizin içerikleri"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Bu da bireysel çıktı seviyesinde ayarlandığından, yayınlama modunu ayrıntılı bir şekilde ayarlamanıza olanak tanır.
Bu, daha sonra çok adımlı iş akışlarına geçtiğimizde özellikle işe yarayacaktır; örneğin yalnızca son çıktıları kopyalamak ve ara çıktıları sembolik bağlantı olarak bırakmak isteyebilirsiniz.

Daha önce belirtildiği gibi, çıktıların nasıl yayınlandığını kontrol etmek için başka, daha sofistike seçenekler de vardır.
Nextflow yolculuğunuzda bunları nasıl kullanacağınızı zamanı geldiğinde göstereceğiz.

### 2.4. Süreç düzeyinde `publishDir` yönergeleri hakkında not

Çok yakın zamana kadar, çıktıları yayınlamanın yerleşik yolu, her bir süreç düzeyinde bir `publishDir` yönergesi kullanmaktı.

`sayHello` sürecinin çıktıları için az önce yaptığımız şeyi başarmak için, bunun yerine süreç tanımına aşağıdaki satırı eklemiş olurduk:

```groovy title="hello-world.nf" linenums="6" hl_lines="3"
process sayHello {

    publishDir 'results/hello_world', mode: 'copy'

    output:
    path 'output.txt'

    script:
    """
    echo 'Hello World!' > output.txt
    """
}
```

Bu kod kalıbını eski Nextflow iş akışlarında ve süreç modüllerinde hâlâ her yerde bulacaksınız, bu nedenle bunun farkında olmak önemlidir.
Ancak, Nextflow dilinin gelecek sürümlerinde sonunda izin verilmeyeceğinden, yeni çalışmalarda kullanmanızı önermiyoruz.

### Özet

İş akışı çıktılarını daha uygun bir konuma nasıl yayınlayacağınızı biliyorsunuz.

### Sırada ne var?

Komut satırı parametresi aracılığıyla değişken girdi sağlamayı ve varsayılan değerleri etkili bir şekilde kullanmayı öğrenin.

---

## 3. Komut satırında iletilen değişken girdi kullanın

Mevcut durumunda, iş akışımız süreç komutuna sabit kodlanmış bir selamlama kullanır.
Çalışma zamanında selamlamayı daha kolay değiştirebilmemiz için bir girdi değişkeni kullanarak biraz esneklik eklemek istiyoruz.

Bu, betiğimizde üç dizi değişiklik yapmamızı gerektirir:

1. Sürecin değişken girdi beklemesini sağlayın
2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın
3. Girdiyi iş akışı gövdesindeki sürece iletin

Bu değişiklikleri tek tek yapalım.

### 3.1. `sayHello` sürecinin değişken girdi beklemesini sağlayın

Süreç tanımını (1) bir girdi değişkenini kabul edecek ve (2) bu değişkeni komut satırında kullanacak şekilde düzenlememiz gerekir.

#### 3.1.1. Süreç tanımına bir girdi bloğu ekleyin

İlk olarak, süreç tanımını `greeting` adlı bir girdiyi kabul edecek şekilde uyarlayalım.

Süreç bloğunda aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="6" hl_lines="3-4"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="6"
    process sayHello {

        output:
        path 'output.txt'
    ```

`greeting` değişkeni, bunun bir değer olduğunu (yol değil) Nextflow'a söylemek için `val` ile ön eklenmiştir.

#### 3.1.2. Süreç komutunu girdi değişkenini kullanacak şekilde düzenleyin

Şimdi orijinal sabit kodlanmış değeri, almayı beklediğimiz girdi değişkeninin değeriyle değiştiriyoruz.

Süreç bloğunda aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo '${greeting}' > output.txt
    """
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="14" hl_lines="3"
    script:
    """
    echo 'Hello World!' > output.txt
    """
    ```

`$` simgesi ve süslü parantezler (`{ }`) Nextflow'a bunun gerçek girdi değeriyle değiştirilmesi gereken (=enterpolasyon) bir değişken adı olduğunu söyler.

!!! tip "İpucu"

    Süslü parantezler (`{ }`) teknik olarak Nextflow'un önceki sürümlerinde isteğe bağlıydı, bu nedenle bunun `echo '$greeting' > output.txt` şeklinde yazıldığı eski iş akışları görebilirsiniz.

Artık `sayHello()` süreci değişken girdi kabul etmeye hazır olduğuna göre, iş akışı düzeyinde süreç çağrısına bir girdi değeri sağlamanın bir yoluna ihtiyacımız var.

### 3.2. Kullanıcı girdisini yakalamak için komut satırı parametresi ayarlayın

Girdiyi doğrudan süreç çağrısını `sayHello('Hello World!')` yaparak sabit kodlayabiliriz.
Ancak, iş akışımızla gerçek iş yaptığımızda, girdilerini komut satırından kontrol edebilmek isteyeceğiz, böylece şöyle bir şey yapabiliriz:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Neyse ki, Nextflow'un [yerleşik](https://nextflow.io/docs/latest/config.html#params) iş akışı parametre sistemi `params` vardır ve bu, CLI parametrelerini bildirmeyi ve kullanmayı kolaylaştırır.

Genel söz dizimi, komut satırında bir `--<parametre_adı>` parametresi beklemek için `params.<parametre_adı>` bildirmektir.

Burada `--input` adlı bir parametre oluşturmak istiyoruz, bu nedenle iş akışında bir yerde `params.input` bildirmemiz gerekiyor.
Prensipte bunu herhangi bir yere yazabiliriz; ancak bunu `sayHello()` süreç çağrısına vermek isteyeceğimiz için, doğrudan `sayHello(params.input)` yazarak oraya ekleyebiliriz.

İş akışı bloğunda aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // bir selamlama yayınla
    sayHello(params.input)
    ```

=== "Önce"

    ```groovy title="hello-world.nf" linenums="23" hl_lines="2"
    // bir selamlama yayınla
    sayHello()
    ```

Bu, Nextflow'a `sayHello` sürecini `--input` parametresi aracılığıyla sağlanan değer üzerinde çalıştırmasını söyler.

Aslında, bölümün başında belirtilen (2) ve (3) adımlarını tek seferde başardık.

### 3.3. İş akışı komutunu çalıştırın

Çalıştıralım!

```bash
nextflow run hello-world.nf --input 'Bonjour le monde!'
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elated_lavoisier] DSL2 - revision: 7c031b42ea

    executor >  local (1)
    [4b/654319] sayHello | 1 of 1 ✔
    ```

Tüm bu düzenlemeleri doğru yaptıysanız, başka bir başarılı yürütme elde etmelisiniz.

Artık selamlamanın yeni sürümüne sahip olduğunuzdan emin olmak için çıktı dosyasını açtığınızdan emin olun.

??? abstract "Dosya içerikleri"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Voilà!

Yeni yürütmenin `results` dizinine yayınlanan çıktı dosyasını nasıl üzerine yazdığına dikkat edin.
Ancak, önceki çalıştırmaların sonuçları `work` altındaki görev dizinlerinde hâlâ korunmaktadır.

!!! tip "İpucu"

    Nextflow düzeyindeki parametreleri iş akışı düzeyindeki parametrelerden kolayca ayırt edebilirsiniz.

    - Bir iş akışına uygulanan parametreler her zaman çift tire (`--`) alır.
    - Bir Nextflow ayarını değiştiren parametreler, örneğin daha önce kullandığımız `-resume` özelliği, tek tire (`-`) alır.

### 3.4. Komut satırı parametreleri için varsayılan değerler kullanın

Tamam, bu uygundu, ancak birçok durumda, her çalıştırma için belirtmek zorunda kalmamanız için belirli bir parametre için varsayılan bir değer sağlamak mantıklıdır.

#### 3.4.1. CLI parametresi için varsayılan değer ayarlayın

İş akışı tanımından önce bildirerek `input` parametresine varsayılan bir değer verelim.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline parametreleri
 */
params {
    input: String = 'Holà mundo!'
}
```

Gördüğünüz gibi, iş akışının beklediği girdi türünü belirtebiliriz (Nextflow 25.10.2 ve sonrası).
Söz dizimi `ad: Tür = varsayılan_değer`'dir.
Desteklenen türler arasında `String`, `Integer`, `Float`, `Boolean` ve `Path` bulunur.

!!! info "Bilgi"

    Eski iş akışlarında, tüm bu `params` bloğunun sadece `input = 'Holà mundo!'` olarak yazıldığını görebilirsiniz.

İş akışınıza daha fazla parametre ekledikçe, varsayılan değer vermeniz gerekip gerekmediğine bakılmaksızın hepsini bu bloğa eklemelisiniz.
Bu, yapılandırılabilir tüm parametreleri bir bakışta bulmayı kolaylaştıracaktır.

#### 3.4.2. İş akışını parametreyi belirtmeden tekrar çalıştırın

Artık varsayılan bir değer ayarladığınıza göre, komut satırında bir değer belirtmek zorunda kalmadan iş akışını tekrar çalıştırabilirsiniz.

```bash
nextflow run hello-world.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [determined_edison] DSL2 - revision: 3539118582

    executor >  local (1)
    [72/394147] sayHello | 1 of 1 ✔
    ```

Çıktı daha önce olduğu gibi aynı yerde olacak, ancak içerikler yeni metinle güncellenmiş olmalıdır.

??? abstract "Dosya içerikleri"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow, çıktıyı oluşturmak için selamlama parametresinin varsayılan değerini kullandı.

#### 3.4.3. Varsayılan değeri geçersiz kılın

Parametreyi komut satırında sağlarsanız, CLI değeri varsayılan değeri geçersiz kılar.

Deneyin:

```bash
nextflow run hello-world.nf --input 'Konnichiwa!'
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [elegant_faraday] DSL2 - revision: 3539118582

    executor >  local (1)
    [6f/a12a91] sayHello | 1 of 1 ✔
    ```

Bir kez daha, sonuçlar dizininizde karşılık gelen güncellenmiş çıktıyı bulmalısınız.

??? abstract "Dosya içerikleri"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note "Not"

    Nextflow'ta, parametreler için değerleri belirleyebileceğiniz birden fazla yer vardır.
    Aynı parametre birden fazla yerde farklı değerlere ayarlanırsa, Nextflow hangi değeri kullanacağını [burada](https://www.nextflow.io/docs/latest/config.html) açıklanan öncelik sırasına göre belirler.

    Bunu Bölüm 6'da (Yapılandırma) daha ayrıntılı olarak ele alacağız.

### Özet

Komut satırı parametresi aracılığıyla çalışma zamanında sağlanan basit değişken girdiyi, ayrıca varsayılan değerleri ayarlamayı, kullanmayı ve geçersiz kılmayı nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

İş akışı yürütmelerini daha rahat nasıl yöneteceğinizi öğrenin.

---

## 4. İş akışı yürütmelerini yönetin

İş akışlarını nasıl başlatacağınızı ve çıktıları nasıl alacağınızı bilmek harikadır, ancak özellikle kendi iş akışlarınızı geliştiriyorsanız, hayatınızı kolaylaştıracak iş akışı yönetiminin birkaç diğer yönü olduğunu çabucak göreceksiniz.

Burada, aynı iş akışını yeniden başlatmanız gerektiğinde [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) özelliğini nasıl kullanacağınızı, [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) ile geçmiş yürütmelerin günlüğünü nasıl inceleyeceğinizi ve [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) ile eski work dizinlerini nasıl sileceğinizi gösteriyoruz.

### 4.1. `-resume` ile bir iş akışını yeniden başlatın

Bazen, daha önce başlattığınız bir iş akışını, zaten başarıyla tamamlanmış adımları tekrarlamadan yeniden çalıştırmak isteyeceksiniz.

Nextflow'un bunu yapmanızı sağlayan [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) adlı bir seçeneği vardır.
Özellikle, bu modda, zaten tam olarak aynı kod, ayarlar ve girdilerle çalıştırılmış olan tüm süreçler atlanacaktır.
Bu, Nextflow'un yalnızca son çalıştırmadan bu yana eklediğiniz veya değiştirdiğiniz süreçleri veya yeni ayarlar veya girdiler sağladığınız süreçleri çalıştıracağı anlamına gelir.

Bunu yapmanın iki önemli avantajı vardır:

- İş akışınızı geliştirme aşamasındaysanız, değişikliklerinizi test etmek için yalnızca üzerinde aktif olarak çalıştığınız süreci(leri) çalıştırmanız gerektiğinden daha hızlı iterasyon yapabilirsiniz.
- Bir iş akışını üretimde çalıştırıyorsanız ve bir şeyler ters giderse, birçok durumda sorunu düzeltip iş akışını yeniden başlatabilirsiniz ve arıza noktasından çalışmaya devam eder, bu da size çok zaman ve hesaplama tasarrufu sağlayabilir.

Kullanmak için komutunuza `-resume` eklemeniz ve çalıştırmanız yeterlidir:

```bash
nextflow run hello-world.nf -resume
```

??? success "Komut çıktısı"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

Konsol çıktısı tanıdık görünmeli, ancak daha öncesine kıyasla biraz farklı bir şey var.

Süreç durumu satırına (satır 5) eklenen `cached:` kısmına bakın, bu Nextflow'un bu işi zaten yaptığını ve önceki başarılı çalıştırmanın sonucunu yeniden kullandığını gösterir.

Work alt dizini hash'inin önceki çalıştırmayla aynı olduğunu da görebilirsiniz.
Nextflow kelimenin tam anlamıyla sizi önceki yürütmeye yönlendiriyor ve "Bunu zaten orada yaptım" diyor.

!!! tip "İpucu"

    `resume` ile bir iş akışını yeniden çalıştırdığınızda, Nextflow daha önce başarıyla çalıştırılmış yürütmeler tarafından work dizini dışına yayınlanan dosyaların üzerine yazmaz.

### 4.2. Geçmiş yürütmelerin günlüğünü inceleyin

İster yeni bir iş akışı geliştiriyor olun, ister iş akışlarını üretimde çalıştırıyor olun, bir noktada muhtemelen geçmiş çalıştırmalar hakkında bilgi aramanız gerekecektir.
İşte bunu nasıl yapacağınız.

Bir nextflow iş akışı başlattığınızda, mevcut çalışma dizininde `.nextflow` adlı gizli bir dizin altında `history` adlı bir günlük dosyasına bir satır yazılır.

??? abstract "Dosya içerikleri"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Bu dosya, mevcut çalışma dizininden başlatılan her Nextflow çalıştırması için zaman damgası, çalıştırma adı, durum, revizyon kimliği, oturum kimliği ve tam komut satırını verir.

Bu bilgilere erişmenin daha uygun bir yolu `nextflow log` komutunu kullanmaktır.

```bash
nextflow log
```

??? success "Komut çıktısı"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Bu, günlük dosyasının içeriğini terminale çıkarır ve bir başlık satırıyla zenginleştirir.

Yeni bir `nextflow run` komutu çalıştırdığınızda oturum kimliğinin değiştiğini, ANCAK `-resume` seçeneğini kullanıyorsanız değişmediğini fark edeceksiniz.
Bu durumda oturum kimliği aynı kalır.

Nextflow, çalıştırma önbellek bilgilerini `.nextflow` altında bulunan `cache` dizini altında gruplandırmak için oturum kimliğini kullanır.

### 4.3. Eski work dizinlerini silin

Geliştirme sürecinde, taslak iş akışınızı genellikle çok sayıda çalıştırırsınız, bu da birçok alt dizinde birçok dosyanın birikmesine yol açabilir.

Neyse ki Nextflow, artık umursamadığınız geçmiş çalıştırmaların work alt dizinlerini otomatik olarak silebilen yararlı bir `clean` alt komutu içerir.

#### 4.3.1. Silme kriterlerini belirleyin

Neyin silineceğini belirlemek için birden fazla [seçenek](https://www.nextflow.io/docs/latest/reference/cli.html#clean) vardır.

Burada, çalıştırma adı kullanılarak belirtilen belirli bir çalıştırmadan önce tüm alt dizinleri silen bir örnek gösteriyoruz.

`-resume` kullanmadığınız en son başarılı çalıştırmaya bakın; bizim durumumuzda çalıştırma adı `golden_cantor` idi.

Çalıştırma adı, `Launching (...)` konsol çıktı satırında köşeli parantez içinde gösterilen makine tarafından oluşturulan iki parçalı dizedir.
Zaman damgası ve/veya komut satırına göre bir çalıştırmayı aramak için Nextflow günlüğünü de kullanabilirsiniz.

#### 4.3.2. Kuru çalıştırma yapın

Önce komut verildiğinde nelerin silineceğini kontrol etmek için kuru çalıştırma bayrağı `-n`'yi kullanıyoruz:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Komut çıktısı"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Çıktınız farklı görev dizini adlarına sahip olacak ve farklı sayıda satır olabilir, ancak örneğe benzer görünmelidir.

Herhangi bir satır çıktısı görmüyorsanız, ya geçerli bir çalıştırma adı sağlamadınız ya da silinecek geçmiş çalıştırma yok. Örnek komuttaki `golden_cantor`'u günlüğünüzdeki karşılık gelen en son çalıştırma adına değiştirdiğinizden emin olun.

#### 4.3.3. Silme işlemine devam edin

Çıktı beklendiği gibi görünüyorsa ve silme işlemine devam etmek istiyorsanız, `-n` yerine `-f` bayrağıyla komutu yeniden çalıştırın:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Komut çıktısı"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Çıktı daha öncekine benzer olmalı, ancak şimdi 'Would remove' yerine 'Removed' diyor.
Bunun iki karakterlik alt dizinleri (yukarıdaki `a3/` gibi) kaldırmadığını, ancak içeriklerini boşalttığını unutmayın.

!!! Warning "Uyarı"

    Geçmiş çalıştırmalardan work alt dizinlerini silmek onları Nextflow'un önbelleğinden kaldırır ve bu dizinlerde depolanan tüm çıktıları siler.
    Bu, Nextflow'un karşılık gelen süreçleri yeniden çalıştırmadan yürütmeyi sürdürme yeteneğini bozar.

    Umursadığınız veya güvenmeyi planladığınız çıktıları kaydetmekten siz sorumlusunuz! Bu, `publish` yönergesi için `symlink` modu yerine `copy` modunu tercih etmemizin ana nedenidir.

### Özet

Çıktıları belirli bir dizine nasıl yayınlayacağınızı, zaten aynı şekilde çalıştırılmış adımları tekrarlamadan bir iş akışını nasıl yeniden başlatacağınızı ve eski work dizinlerini temizlemek için `nextflow clean` komutunu nasıl kullanacağınızı biliyorsunuz.

Daha genel olarak, basit bir Nextflow iş akışını nasıl yorumlayacağınızı, yürütmesini nasıl yöneteceğinizi ve çıktıları nasıl alacağınızı biliyorsunuz.

### Sırada ne var?

Küçük bir mola verin, hak ettiniz!

Hazır olduğunuzda, iş akışınıza girdileri beslemek için kanalları nasıl kullanacağınızı öğrenmek için [**Bölüm 2: Hello Channels**](./02_hello_channels.md)'a geçin; bu, Nextflow'un yerleşik veri akışı paralelliğinden ve diğer güçlü özelliklerden yararlanmanızı sağlayacaktır.

---

## Quiz

<quiz>
Bir Nextflow sürecinin minimum gerekli bileşenleri nelerdir?
- [ ] Yalnızca input ve output blokları
- [x] Output ve script blokları
- [ ] Input, output ve script blokları
- [ ] Yalnızca script bloğu

Daha fazla bilgi: [1.1.1. Süreç tanımı](#111-process-tanımı)
</quiz>

<quiz>
Bir süreçteki output bloğunun amacı nedir?
- [ ] Sonuçları konsola yazdırmak
- [ ] Dosyaları work dizinine kaydetmek
- [x] Süreçten beklenen çıktıları bildirmek
- [ ] Ortam değişkenlerini tanımlamak

Daha fazla bilgi: [1.1.1. Süreç tanımı](#111-process-tanımı)
</quiz>

<quiz>
Bir Nextflow iş akışını çalıştırmak için hangi komut kullanılır?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Bir görevin work dizinine bakıldığında, hangi dosya gerçekte yürütülen komutu içerir?

```
work/a3/7be2fa.../
├── .command.begin
├── .command.err
├── .command.log
├── .command.out
├── .command.run
├── .command.sh
├── .exitcode
└── output.txt
```

- [ ] `.command.run`
- [x] `.command.sh`
- [ ] `.command.log`
- [ ] `.command.out`

Daha fazla bilgi: [1.2.2. `work` dizininde çıktıyı ve günlükleri bulun](#122-work-dizininde-çıktıyı-ve-günlükleri-bulun)
</quiz>

<quiz>
`-resume` bayrağı ne yapar?
- [ ] İş akışını baştan başlatır
- [ ] İş akışını duraklatır
- [x] Zaten başarıyla tamamlanmış süreçleri atlar
- [ ] İş akışının yedeğini oluşturur

Daha fazla bilgi: [4.1. `-resume` ile bir iş akışını yeniden başlatın](#41--resume-ile-bir-iş-akışını-yeniden-başlatın)
</quiz>

<quiz>
İş akışı çıktılarını yayınlamak için varsayılan mod nedir?
- [ ] Dosyaları çıktı dizinine kopyala
- [x] Çıktı dizininde sembolik bağlantılar oluştur
- [ ] Dosyaları çıktı dizinine taşı
- [ ] Dosyaları çıktı dizininde sıkıştır

Daha fazla bilgi: [2.3. Yayınlama modunu kopyalamaya ayarlayın](#23-yayınlama-modunu-kopyalamaya-ayarlayın)
</quiz>

<quiz>
Komut satırından bir Nextflow iş akışına parametre değeri nasıl iletilir?
- [ ] `-parameter değer`
- [ ] `--parameter:değer`
- [x] `--parameter değer`
- [ ] `-p parameter=değer`

Daha fazla bilgi: [3.2. Kullanıcı girdisini yakalamak için komut satırı parametresi ayarlayın](#32-kullanıcı-girdisini-yakalamak-için-komut-satırı-parametresi-ayarlayın)
</quiz>

<quiz>
Nextflow script bloğu içinde bir değişkene nasıl başvurulur?
- [ ] `%değişken%` söz dizimini kullan
- [x] `#!groovy ${değişken}` söz dizimini kullan
- [ ] `{{değişken}}` söz dizimini kullan
- [ ] `[değişken]` söz dizimini kullan
</quiz>
