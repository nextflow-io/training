# Bölüm 1: Hello World

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube kanalında görün.

:green_book: Video metni [burada](./transcripts/01_hello_world.md) mevcuttur.
///

Hello Nextflow eğitim kursunun bu ilk bölümünde, temel Nextflow mantığı ve bileşenlerinin kullanımını göstermek için aşamalı olarak geliştireceğimiz çok basit, alana özgü olmayan bir Hello World örneğiyle konuya giriş yapıyoruz.

??? info "Hello World örneği nedir?"

    "Hello World!" bir programlama dilinin veya yazılım çerçevesinin temel sözdizimini ve yapısını göstermeyi amaçlayan minimalist bir örnektir.
    Örnek tipik olarak "Hello, World!" ifadesini konsol veya terminal gibi çıktı cihazına yazdırmaktan veya bir dosyaya yazmaktan oluşur.

---

## 0. Isınma: Bir Hello World örneğini doğrudan çalıştırın

Bunu, Nextflow'a sarmadan önce ne yaptığını göstermek için doğrudan terminalde çalıştırdığımız basit bir komutla gösterelim.

!!! tip

    [Başlarken](00_orientation.md) sayfasında açıklandığı gibi şu anda `hello-nextflow/` dizininin içinde olmanız gerektiğini unutmayın.

### 0.1. Terminalin merhaba demesini sağlayın

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
echo 'Hello World!'
```

??? success "Komut çıktısı"

    ```console
    Hello World!
    ```

Bu, 'Hello World' metnini doğrudan terminalde çıktı olarak verir.

### 0.2. Çıktıyı bir dosyaya yazın

Pipeline'ları çalıştırmak çoğunlukla dosyalardan veri okumayı ve sonuçları diğer dosyalara yazmayı içerir, bu yüzden örneği biraz daha alakalı hale getirmek için komutu metin çıktısını bir dosyaya yazacak şekilde değiştirelim.

```bash
echo 'Hello World!' > output.txt
```

??? success "Komut çıktısı"

    ```console

    ```

Bu terminale herhangi bir şey çıktı vermez.

### 0.3. Çıktıyı bulun

'Hello World' metni şimdi belirttiğimiz `output.txt` adlı çıktı dosyasında olmalıdır.
Dosya gezgininde veya örneğin `cat` yardımcı programını kullanarak komut satırından açabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

İlk Nextflow iş akışımızla çoğaltmaya çalışacağımız şey budur.

### Özet

Artık terminalde biraz metin çıktısı veren basit bir komutun nasıl çalıştırılacağını ve isteğe bağlı olarak çıktıyı bir dosyaya nasıl yazdıracağınızı biliyorsunuz.

### Sırada ne var?

Bunun Nextflow iş akışı olarak yazıldığında nasıl görüneceğini öğrenin.

---

## 1. Betiği inceleyin ve çalıştırın

Size daha önce yaptığımız şeyi (yani 'Hello World!' yazmayı) yapan ancak Nextflow ile yapan, tam işlevsel ama minimalist bir iş akışı betiği sağlıyoruz: `hello-world.nf`.

Başlamanız için, iş akışı betiğini açalım, böylece nasıl yapılandırıldığına dair bir fikir edinebilirsiniz.
Ardından onu çalıştıracağız ve çıktılarını arayacağız.

### 1.1. Kodu inceleyin

`hello-world.nf` betiğini mevcut dizininizde bulacaksınız, bu dizin `hello-nextflow` olmalıdır. Editör panelinde açın.

??? full-code "Tam kod dosyası"

    ```groovy title="hello-world.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' yazdırmak için echo kullan
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

Bir Nextflow iş akışı betiği tipik olarak bir veya daha fazla [**process**](https://nextflow.io/docs/latest/process.html) tanımı ve [**workflow**](https://nextflow.io/docs/latest/workflow.html)'un kendisini içerir, ayrıca birkaç isteğe bağlı blok (burada mevcut değil) içerir ki bunları daha sonra tanıtacağız.

Her **process**, pipeline'daki ilgili adımın hangi işlem(ler)i gerçekleştirmesi gerektiğini tanımlarken, **workflow** çeşitli adımları birbirine bağlayan veri akışı mantığını tanımlar.

Önce **process** bloğuna daha yakından bakacağız, ardından **workflow** bloğuna bakacağız.

#### 1.1.1. `process` tanımı

İlk kod bloğu bir **process**'i tanımlar.

Process tanımı `process` anahtar kelimesiyle başlar, ardından süreç adı ve son olarak süslü parantezlerle sınırlandırılmış süreç gövdesi gelir.
Süreç gövdesi, çalıştırılacak komutu belirten bir script bloğu içermelidir; bu, komut satırı terminalinde çalıştırabileceğiniz herhangi bir şey olabilir.

```groovy title="hello-world.nf" linenums="3"
/*
* 'Hello World!' yazdırmak için echo kullan
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

Burada `output.txt` adlı bir dosyaya **output** yazan `sayHello` adlı bir **process**'imiz var.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

Bu, sadece bir `output` tanımı ve çalıştırılacak `script`'i içeren çok minimal bir süreç tanımıdır.

`output` tanımı `path` niteleyicisini içerir; bu, Nextflow'a bunun bir yol olarak ele alınması gerektiğini söyler (hem dizin yollarını hem de dosyaları içerir).
Bir diğer yaygın niteleyici `val`'dir.

Önemli olarak, output tanımı hangi çıktının oluşturulacağını _belirlemez_.
Sadece beklenen çıktının ne olduğunu _bildirir_, böylece Nextflow yürütme tamamlandığında onu arayabilir.
Bu, komutun başarıyla yürütüldüğünü doğrulamak ve gerekirse çıktıyı aşağı akış süreçlerine iletmek için gereklidir. Output bloğunda bildirilenle eşleşmeyen üretilen çıktı, aşağı akış süreçlerine iletilmeyecektir.

!!! warning

    Bu örnek kırılgandır çünkü çıktı dosya adını iki ayrı yerde (script ve output blokları) sabit kodladık.
    Birini değiştirip diğerini değiştirmezsek, betik bozulacaktır.
    Daha sonra, bu sorunu azaltmak için değişkenleri kullanmanın yollarını öğreneceksiniz.

Gerçek dünya pipeline'ında, bir süreç genellikle yönergeler ve girdiler gibi ek bloklar içerir; bunları birazdan tanıtacağız.

#### 1.1.2. `workflow` tanımı

İkinci kod bloğu **workflow**'un kendisini tanımlar.
Workflow tanımı `workflow` anahtar kelimesiyle başlar, ardından isteğe bağlı bir ad gelir, sonra süslü parantezlerle sınırlandırılmış workflow gövdesi gelir.

Burada `main:` bloğu ('bu workflow'un ana gövdesidir' diyen) içeren ve `sayHello` sürecine bir çağrı içeren bir **workflow**'umuz var.

```groovy title="hello-world.nf" linenums="17"
workflow {

    main:
    // bir selamlama yayınla
    sayHello()
}
```

Bu çok minimal bir **workflow** tanımıdır.
Gerçek dünya pipeline'ında, workflow tipik olarak **channel**'lar tarafından bağlanan **process**'lere birden fazla çağrı içerir ve süreçler bir veya daha fazla değişken **input** bekler.

Daha sonra bu eğitim modülünde değişken girdilerin nasıl ekleneceğini öğreneceksiniz; ve bu kursun 3. Bölümünde daha fazla süreç eklemeyi ve bunları kanallarla bağlamayı öğreneceksiniz.

!!! tip

    Teknik olarak `main:` satırı bunun gibi basit iş akışları için gerekli değildir, bu nedenle buna sahip olmayan iş akışlarıyla karşılaşabilirsiniz.
    Ancak workflow düzeyinde çıktılardan yararlanmak için buna ihtiyacımız olacak, bu yüzden baştan dahil etsek iyi olur.

### 1.2. İş akışını çalıştırın

Koda bakmak onu çalıştırmak kadar eğlenceli değil, o yüzden bunu pratikte deneyelim.

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

Bu bize `sayHello` sürecinin başarıyla bir kez yürütüldüğünü söyler (`1 of 1 ✔`).

Önemli olarak, bu satır ayrıca `sayHello` süreç çağrısının çıktısını nerede bulacağınızı da söyler.
Şimdi buna bakalım.

#### 1.2.2. `work` dizininde çıktıyı ve günlükleri bulun

Nextflow'u belirli bir dizinde ilk kez çalıştırdığınızda, yürütme sırasında oluşturulan tüm dosyaları (ve sembolik bağlantıları) yazacağı `work` adlı bir dizin oluşturur.

`work` dizini içinde, Nextflow çıktıları ve günlükleri süreç çağrısı başına düzenler.
Her süreç çağrısı için Nextflow, benzersiz olması için bir hash ile adlandırılmış iç içe bir alt dizin oluşturur; burada tüm gerekli girdileri hazırlar (varsayılan olarak sembolik bağlantılar kullanarak), yardımcı dosyalar yazar ve günlükleri ve sürecin herhangi bir çıktısını yazar.

Bu alt dizinin yolu, konsol çıktısında kısaltılmış biçimde köşeli parantez içinde gösterilir.
Yukarıda gösterilen çalıştırma için aldığımıza bakarsak, sayHello süreci için konsol günlük satırı `[65/7be2fa]` ile başlar. Bu, şu dizin yoluna karşılık gelir: `work/65/7be2fad5e71e5f49998f795677fd68`

Orada ne olduğuna bir göz atalım.

??? abstract "Dizin içeriği"

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

Bakmak istediğiniz ilk şey, iş akışının gerçek çıktısıdır, yani `sayHello` süreci tarafından üretilen `output.txt` dosyası.
Açın ve minimalist iş akışımızın amacı olan `Hello World!` selamlamasını bulacaksınız.

??? abstract "Dosya içeriği"

    ```console title="output.txt"
    Hello World!
    ```

İşe yaradı!

Kabul etmeliyiz ki, bu kadar küçük bir sonuç için çok fazla sarmalayıcı kod gibi görünebilir, ancak tüm bu sarmalayıcı kodun değeri, girdi dosyalarını okumaya ve birden fazla adımı birbirine bağlamaya başladığımızda daha belirgin hale gelecektir.

Bununla birlikte, o dizindeki diğer dosyalara da bakalım. Bunlar, Nextflow tarafından görev yürütmesinin bir parçası olarak üretilen yardımcı ve günlük dosyalarıdır.

- **`.command.begin`**: Süreç çağrısının yürütmesinin başlangıcıyla ilgili meta veriler
- **`.command.err`**: Süreç çağrısı tarafından yayılan hata mesajları (`stderr`)
- **`.command.log`**: Süreç çağrısı tarafından yayılan tam günlük çıktısı
- **`.command.out`**: Süreç çağrısı tarafından yayılan normal çıktı (`stdout`)
- **`.command.run`**: Süreç çağrısını yürütmek için Nextflow tarafından çalıştırılan tam betik
- **`.command.sh`**: Süreç çağrısı tarafından gerçekten çalıştırılan komut
- **`.exitcode`**: Komuttan kaynaklanan çıkış kodu

`.command.sh` dosyası özellikle yararlıdır çünkü Nextflow'un yürüttüğü ana komutu size söyler, tüm defter tutma ve görev/ortam kurulumunu içermez.

??? abstract "Dosya içeriği"

    ```console title=".command.sh"
    #!/bin/bash -ue
    echo 'Hello World!' > output.txt
    ```

Bu, daha önce manuel olarak çalıştırdığımızla eşleşir.

Bu durumda çok basittir çünkü süreç komutu sabit kodlanmıştı, ancak kursta ilerledikçe bazı değişken interpolasyonu içeren süreç komutları göreceksiniz.
Bu, başarısız bir çalıştırmada sorun giderirken Nextflow'un kodu tam olarak nasıl yorumladığını ve hangi komutun üretildiğini görebilmenizi özellikle değerli kılar.

### 1.3. İş akışını tekrar çalıştırın

İş akışını birkaç kez yeniden çalıştırmayı deneyin, ardından `work/` altındaki görev dizinlerine bakın.

??? abstract "Dizin içeriği"

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

Her çalıştırma için tam bir çıktı ve günlük dosyası seti ile yeni bir alt dizinin oluşturulduğunu görüyorsunuz.
Bu size aynı iş akışını birkaç kez çalıştırmanın önceki çalıştırmaların sonuçlarının üzerine yazmayacağını gösterir.

### Özet

Basit bir Nextflow betiğini nasıl çözeceğinizi, çalıştıracağınızı ve work dizininde çıktıyı ve ilgili günlük dosyalarını nasıl bulacağınızı biliyorsunuz.

### Sırada ne var?

İş akışı çıktılarını daha uygun bir konuma nasıl yayınlayacağınızı öğrenin.

---

## 2. Çıktıları yayınlayın

Az önce öğrendiğiniz gibi, pipeline'ımız tarafından üretilen çıktı, birkaç katman derinlikte bir çalışma dizinine gömülüdür.
Bu kasıtlı olarak yapılır; Nextflow bu dizinin kontrolündedir ve onunla etkileşime girmememiz gerekir.
Ancak bu, önemsediğimiz çıktıları almayı zorlaştırır.

Neyse ki, Nextflow [workflow output definitions](https://nextflow.io/docs/latest/workflow.html#workflow-outputs) kullanarak çıktıları belirlenmiş bir dizine yayınlamanın bir yolunu sağlar.

### 2.1. Temel kullanım

Bu, iki yeni kod parçası içerecektir:

1. `workflow` gövdesi içinde süreç çıktılarını bildiren bir `publish:` bloğu.
2. Mod ve konum gibi çıktı seçeneklerini belirten betiğe bir `output` bloğu.

#### 2.1.1. `sayHello` sürecinin çıktısını bildirin

Workflow gövdesine bir `publish:` bloğu eklememiz (`main:` bloğuyla aynı tür kod öğesi) ve `sayHello()` sürecinin çıktısını listelememiz gerekir.

`hello-world.nf` iş akışı betik dosyasında aşağıdaki kod satırlarını ekleyin:

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

Sürecin çıktısına basitçe `sayHello().out` yaparak başvurabileceğimizi ve ona rastgele bir ad, `first_output`, atayabileceğimizi görüyorsunuz.

#### 2.1.2. Betiğe bir `output:` bloğu ekleyin

Şimdi sadece çıktı dizin yolunun belirtileceği `output:` bloğunu eklememiz gerekiyor. Bu yeni bloğun betik içinde `workflow` bloğunun **dışında** ve **altında** olduğunu unutmayın.

`hello-world.nf` iş akışı betik dosyasında aşağıdaki kod satırlarını ekleyin:

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
Daha sonra, sofistike çıktı dizin yapıları oluşturmanın yollarını öğreneceksiniz, ancak şimdilik basitlik için sadece minimal bir yol sabit kodluyoruz.

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

Terminal çıktısı tanıdık görünmelidir. Dışarıdan hiçbir şey değişmedi.

Ancak dosya gezgininizi kontrol edin: bu sefer Nextflow, `results/` adlı yeni bir dizin oluşturdu.

??? abstract "Dizin içeriği"

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

`results` dizininin içinde, az önce çalıştırdığımız komut tarafından work dizininde üretilen `output.txt`'ye bir sembolik bağlantı buluyoruz.

Bu, work alt dizinini kazmak zorunda kalmadan çıktı dosyalarını kolayca almamızı sağlar.

### 2.2. Özel bir konum ayarlayın

Varsayılan bir konuma sahip olmak harikadır, ancak sonuçların nereye kaydedileceğini ve nasıl düzenleneceğini özelleştirmek isteyebilirsiniz.

Örneğin, çıktılarınızı alt dizinlere düzenlemek isteyebilirsiniz.
Bunu yapmanın en basit yolu, çıktı başına belirli çıktı yolu atamaktır.

#### 2.2.1. Çıktı yolunu değiştirin

Bir kez daha, belirli bir çıktı için yayınlama davranışını değiştirmek gerçekten basittir.
Özel bir konum ayarlamak için `path`'i buna göre düzenleyin:

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

Bu, bireysel çıktı düzeyinde ayarlandığından, ihtiyaçlarınıza uygun farklı konumlar ve alt dizinler belirtebilirsiniz.

#### 2.2.2. İş akışını tekrar çalıştırın

Hadi deneyelim.

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

Bu sefer sonuç belirtilen alt dizin altına yazılır.

??? abstract "Dizin içeriği"

    ```console hl_lines="2-3"
    results/
    ├── hello_world
    │   └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c2e506b79e2e01acb808d9d12/output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Önceki yürütmeden gelen sonucun hala orada olduğunu görüyorsunuz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_output.svg"
</figure>

İstediğiniz kadar iç içe geçme seviyesi kullanabilirsiniz.
Sonuçları düzenlemek için kullanılan dizinleri adlandırmak için süreç adını veya diğer değişkenleri kullanmak da mümkündür ve üst düzey çıktı dizininin varsayılan adını değiştirmek mümkündür (bu, `-o` CLI bayrağı veya config değişkeni `outputDir` tarafından kontrol edilir).
Bu seçenekleri eğitimde daha sonra ele alacağız.

### 2.3. Yayınlama modunu kopyalama olarak ayarlayın

Varsayılan olarak, çıktılar `work` dizininden sembolik bağlantılar olarak yayınlanır.
Bu, dosya sisteminde yalnızca tek bir dosya olduğu anlamına gelir.

Bu, birden fazla kopyasını saklamak istemediğiniz çok büyük dosyalarla uğraşırken harikadır.
Ancak, herhangi bir noktada work dizinini silerseniz (temizleme işlemlerini kısaca ele alacağız), dosyaya erişiminizi kaybedersiniz.
Bu nedenle, önemli dosyaların kopyalarını güvenli bir yere kaydetmek için bir planınız olması gerekir.

Kolay bir seçenek, önemsediğiniz çıktılar için yayınlama modunu kopyalama olarak değiştirmektir.

#### 2.3.1. Mod yönergesini ekleyin

Bu kısım gerçekten basittir.
İlgili workflow düzeyinde çıktı tanımına sadece `mode 'copy'` ekleyin:

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

Bu, o belirli çıktı için yayınlama modunu ayarlar.

#### 2.3.2. İş akışını tekrar çalıştırın

Hadi deneyelim.

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

Bu sefer sonuçlara bakarsanız, dosya sadece bir sembolik bağlantı yerine uygun bir kopyadır.

??? abstract "Dizin içeriği"

    ```console hl_lines="3"
    results/
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/65/f56f2cd75df1352e106fcdd084b97b/output.txt
    ```

Bu da bireysel çıktı düzeyinde ayarlandığından, yayınlama modunu ayrıntılı bir şekilde ayarlamanıza olanak tanır.
Bu, özellikle daha sonra çok adımlı pipeline'lara geçtiğimizde kullanışlı olacaktır; örneğin, yalnızca nihai çıktıları kopyalamak ve ara çıktıları sembolik bağlantılar olarak bırakmak isteyebilirsiniz.

Daha önce belirtildiği gibi, çıktıların nasıl yayınlanacağını kontrol etmek için başka, daha sofistike seçenekler vardır.
Nextflow yolculuğunuzda zamanı geldiğinde bunları nasıl kullanacağınızı size göstereceğiz.

### 2.4. Süreç düzeyinde `publishDir` yönergeleri hakkında not

Çok yakın zamana kadar, çıktıları yayınlamanın yerleşik yolu, her bir süreç düzeyinde bir `publishDir` yönergesi kullanarak yapmaktı.

`sayHello` sürecinin çıktıları için az önce yaptığımızı başarmak için, bunun yerine süreç tanımına aşağıdaki satırı eklerdik:

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

Bu kod desenini eski Nextflow pipeline'larında ve süreç modüllerinde her yerde bulacaksınız, bu nedenle bunun farkında olmak önemlidir.
Ancak, Nextflow dilinin gelecekteki sürümlerinde sonunda izin verilmeyeceği için bunu yeni herhangi bir çalışmada kullanmanızı önermiyoruz.

### Özet

İş akışı çıktılarını daha uygun bir konuma nasıl yayınlayacağınızı biliyorsunuz.

### Sırada ne var?

Komut satırı parametresi aracılığıyla değişken bir girdi sağlamayı ve varsayılan değerleri etkili bir şekilde kullanmayı öğrenin.

---

## 3. Komut satırında iletilen değişken bir girdi kullanın

Mevcut durumunda, iş akışımız süreç komutuna sabit kodlanmış bir selamlama kullanır.
Çalışma zamanında selamlamayı daha kolay değiştirebilmemiz için bir girdi değişkeni kullanarak biraz esneklik eklemek istiyoruz.

Bu, betiğimizde üç değişiklik yapmamızı gerektirir:

1. Süreci değişken bir girdi bekleyecek şekilde değiştirin
2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın
3. Girdiyi workflow gövdesindeki sürece iletin

Bu değişiklikleri birer birer yapalım.

### 3.1. `sayHello` sürecini değişken bir girdi bekleyecek şekilde değiştirin

Süreç tanımını (1) bir girdi değişkenini kabul edecek ve (2) bu değişkeni komut satırında kullanacak şekilde düzenlememiz gerekiyor.

#### 3.1.1. Süreç tanımına bir input bloğu ekleyin

İlk olarak, süreç tanımını `greeting` adlı bir girdiyi kabul edecek şekilde uyarlayalım.

Process bloğunda aşağıdaki kod değişikliğini yapın:

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

`greeting` değişkeni, Nextflow'a bunun bir değer (yol değil) olduğunu söylemek için `val` ile öneklenir.

#### 3.1.2. Girdi değişkenini kullanmak için süreç komutunu düzenleyin

Şimdi orijinal sabit kodlanmış değeri almayı beklediğimiz girdi değişkeninin değeriyle değiştiriyoruz.

Process bloğunda aşağıdaki kod değişikliğini yapın:

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

`$` sembolü ve süslü parantezler (`{ }`) Nextflow'a bunun gerçek girdi değeriyle değiştirilmesi gereken (=interpolasyon) bir değişken adı olduğunu söyler.

!!! tip

    Süslü parantezler (`{ }`) Nextflow'un önceki sürümlerinde teknik olarak isteğe bağlıydı, bu nedenle bunun `echo '$greeting' > output.txt` olarak yazıldığı eski iş akışları görebilirsiniz.

Artık `sayHello()` süreci değişken bir girdi kabul etmeye hazır olduğuna göre, workflow düzeyinde süreç çağrısına bir girdi değeri sağlamanın bir yoluna ihtiyacımız var.

### 3.2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın

Süreç çağrısını `sayHello('Hello World!')` yaparak doğrudan bir girdiyi sabit kodlayabiliriz.
Ancak, iş akışımızla gerçek iş yaptığımızda, girdilerini komut satırından kontrol edebilmek isteyeceğiz, böylece şöyle bir şey yapabiliriz:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world_input.svg"
</figure>

Neyse ki, Nextflow'un CLI parametrelerini bildirmeyi ve kullanmayı kolaylaştıran [`params`](https://nextflow.io/docs/latest/config.html#params) adlı yerleşik bir iş akışı parametre sistemi vardır.

Genel sözdizimi, Nextflow'a komut satırında bir `--<parameter_name>` parametresi beklemesini söylemek için `params.<parameter_name>` bildirmektir.

Burada `--input` adlı bir parametre oluşturmak istiyoruz, bu nedenle iş akışında bir yerde `params.input` bildirmemiz gerekiyor.
Prensipte bunu herhangi bir yere yazabiliriz; ancak bunu `sayHello()` süreç çağrısına vermek isteyeceğimiz için, `sayHello(params.input)` yazarak doğrudan oraya takabiliriz.

Workflow bloğunda aşağıdaki kod değişikliğini yapın:

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

Aslında, bölümün başında özetlenen (2) ve (3) adımlarını tek seferde gerçekleştirdik.

### 3.3. İş akışı komutunu çalıştırın

Hadi çalıştıralım!

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

Tüm bu düzenlemeleri doğru yaptıysanız, başka bir başarılı yürütme almalısınız.

Artık selamlamanın yeni sürümüne sahip olduğunuzu kontrol etmek için çıktı dosyasını açtığınızdan emin olun.

??? abstract "Dosya içeriği"

    ```console title="results/hello_world/output.txt"
    Bonjour le monde!
    ```

Et voilà!

Yeni yürütmenin `results` dizinine yayınlanan çıktı dosyasının üzerine nasıl yazdığına dikkat edin.
Ancak, önceki çalıştırmaların sonuçları hala `work` altındaki görev dizinlerinde korunur.

!!! tip

    Nextflow düzeyinde parametreleri pipeline düzeyinde parametrelerden kolayca ayırt edebilirsiniz.

    - Bir pipeline'a uygulanan parametreler her zaman çift tire (`--`) alır.
    - Daha önce kullandığımız `-resume` özelliği gibi bir Nextflow ayarını değiştiren parametreler tek tire (`-`) alır.

### 3.4. Komut satırı parametreleri için varsayılan değerler kullanın

Tamam, bu kullanışlıydı, ancak birçok durumda, her çalıştırma için belirtmek zorunda kalmamak için belirli bir parametre için varsayılan bir değer sağlamak mantıklıdır.

#### 3.4.1. CLI parametresi için varsayılan bir değer ayarlayın

`input` parametresine workflow tanımından önce bildirerek varsayılan bir değer verelim.

```groovy title="hello-world.nf" linenums="20"
/*
 * Pipeline parametreleri
 */
params {
    input: String = 'Holà mundo!'
}
```

Gördüğünüz gibi, iş akışının beklediği girdi türünü belirtebiliriz (Nextflow 25.10.2 ve sonrası).
Sözdizimi `name: Type = default_value`'dur.
Desteklenen türler arasında `String`, `Integer`, `Float`, `Boolean` ve `Path` bulunur.

!!! info

    Eski iş akışlarında, tüm `params` bloğunun sadece `input = 'Holà mundo!'` olarak yazıldığını görebilirsiniz.

Pipeline'ınıza daha fazla parametre ekledikçe, varsayılan bir değer vermeniz gerekip gerekmediğine bakılmaksızın hepsini bu bloğa eklemelisiniz.
Bu, tüm yapılandırılabilir parametreleri bir bakışta bulmayı kolaylaştıracaktır.

#### 3.4.2. Parametreyi belirtmeden iş akışını tekrar çalıştırın

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

Çıktı daha önce olduğu gibi aynı yerde olacak, ancak içerik yeni metinle güncellenmiş olmalıdır.

??? abstract "Dosya içeriği"

    ```console title="results/hello_world/output.txt"
    Holà mundo!
    ```

Nextflow, çıktıyı oluşturmak için greeting parametresinin varsayılan değerini kullandı.

#### 3.4.3. Varsayılan değeri geçersiz kılın

Parametreyi komut satırında sağlarsanız, CLI değeri varsayılan değeri geçersiz kılacaktır.

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

Bir kez daha, results dizininizde ilgili güncellenmiş çıktıyı bulmalısınız.

??? abstract "Dosya içeriği"

    ```console title="results/hello_world/output.txt"
    Konnichiwa!
    ```

!!! note

    Nextflow'da, parametreler için değerleri belirtebileceğiniz birden fazla yer vardır.
    Aynı parametre birden fazla yerde farklı değerlere ayarlanırsa, Nextflow hangi değerin kullanılacağını [burada](https://www.nextflow.io/docs/latest/config.html) açıklanan öncelik sırasına göre belirleyecektir.

    Bunu Bölüm 6'da (Yapılandırma) daha ayrıntılı olarak ele alacağız.

### Özet

Komut satırı parametresi aracılığıyla çalışma zamanında sağlanan basit bir değişken girdinin nasıl kullanılacağını ve varsayılan değerlerin nasıl ayarlanacağını, kullanılacağını ve geçersiz kılınacağını biliyorsunuz.

### Sırada ne var?

İş akışı yürütmelerini daha rahat nasıl yöneteceğinizi öğrenin.

---

## 4. İş akışı yürütmelerini yönetin

İş akışlarını nasıl başlatacağınızı ve çıktıları nasıl alacağınızı bilmek harikadır, ancak özellikle kendi iş akışlarınızı geliştiriyorsanız, hayatınızı kolaylaştıracak iş akışı yönetiminin birkaç başka yönü olduğunu hızla göreceksiniz.

Burada size aynı iş akışını yeniden başlatmanız gerektiğinde [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) özelliğinin nasıl kullanılacağını, geçmiş yürütmelerin günlüğünü [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) ile nasıl inceleyeceğinizi ve eski work dizinlerini [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) ile nasıl sileceğinizi gösteriyoruz.

### 4.1. `-resume` ile bir iş akışını yeniden başlatın

Bazen, daha önce başlattığınız bir pipeline'ı, zaten başarıyla tamamlanmış adımları yeniden yapmadan yeniden çalıştırmak isteyeceksiniz.

Nextflow'un bunu yapmanıza izin veren [`-resume`](https://nextflow.io/docs/latest/cache-and-resume.html) adlı bir seçeneği vardır.
Özellikle, bu modda, tam olarak aynı kod, ayarlar ve girdilerle zaten çalıştırılmış olan tüm süreçler atlanacaktır.
Bu, Nextflow'un yalnızca son çalıştırmadan bu yana eklediğiniz veya değiştirdiğiniz süreçleri veya yeni ayarlar veya girdiler sağladığınız süreçleri çalıştıracağı anlamına gelir.

Bunu yapmanın iki temel avantajı vardır:

- Pipeline'ınızı geliştirme aşamasındaysanız, değişikliklerinizi test etmek için yalnızca aktif olarak üzerinde çalıştığınız süreci çalıştırmanız gerektiğinden daha hızlı yineleme yapabilirsiniz.
- Bir pipeline'ı üretimde çalıştırıyorsanız ve bir şeyler ters giderse, birçok durumda sorunu düzeltebilir ve pipeline'ı yeniden başlatabilirsiniz ve başarısızlık noktasından çalışmaya devam edecektir, bu da size çok fazla zaman ve hesaplama tasarrufu sağlayabilir.

Kullanmak için komutunuza `-resume` ekleyin ve çalıştırın:

```bash
nextflow run hello-world.nf -resume
```

??? success "Komut çıktısı"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-world.nf` [golden_cantor] DSL2 - revision: 35bd3425e5

    [62/49a1f8] sayHello | 1 of 1, cached: 1 ✔
    ```

Konsol çıktısı tanıdık görünmelidir, ancak daha öncekine göre biraz farklı bir şey var.

Süreç durum satırına (satır 5) eklenmiş olan `cached:` kısmını arayın; bu, Nextflow'un bu işi zaten yaptığını tanıdığı ve önceki başarılı çalıştırmadan sonucu yeniden kullandığı anlamına gelir.

Work alt dizin hash'inin önceki çalıştırmadakiyle aynı olduğunu da görebilirsiniz.
Nextflow kelimenin tam anlamıyla size önceki yürütmeyi gösteriyor ve "Bunu zaten orada yaptım" diyor.

!!! tip

    Bir pipeline'ı `resume` ile yeniden çalıştırdığınızda, Nextflow daha önce başarıyla çalıştırılan yürütmeler tarafından work dizini dışında yayınlanan hiçbir dosyanın üzerine yazmaz.

### 4.2. Geçmiş yürütmelerin günlüğünü inceleyin

Yeni bir pipeline geliştiriyor veya pipeline'ları üretimde çalıştırıyor olun, bir noktada muhtemelen geçmiş çalıştırmalar hakkında bilgi aramanız gerekecektir.
İşte bunu nasıl yapacağınız.

Bir nextflow iş akışını her başlattığınızda, mevcut çalışma dizinindeki `.nextflow` adlı gizli bir dizin altında `history` adlı bir günlük dosyasına bir satır yazılır.

??? abstract "Dosya içeriği"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Bu dosya size mevcut çalışma dizininden başlatılan her Nextflow çalıştırması için zaman damgası, çalıştırma adı, durum, revizyon kimliği, oturum kimliği ve tam komut satırını verir.

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

Bu, günlük dosyasının içeriğini bir başlık satırıyla birlikte terminale çıktı olarak verecektir.

Oturum kimliğinin, `-resume` seçeneğini kullanmadığınız SÜRECE yeni bir `nextflow run` komutu çalıştırdığınızda her değiştiğini fark edeceksiniz.
Bu durumda, oturum kimliği aynı kalır.

Nextflow, çalıştırma önbellekleme bilgilerini `.nextflow` altında bulunan `cache` dizini altında gruplamak için oturum kimliğini kullanır.

### 4.3. Eski work dizinlerini silin

Geliştirme sürecinde, tipik olarak taslak pipeline'ınızı çok sayıda kez çalıştıracaksınız, bu da birçok alt dizinde birçok dosyanın birikmesine yol açabilir.

Neyse ki Nextflow, artık önemsemediğiniz geçmiş çalıştırmalar için work alt dizinlerini otomatik olarak silebilen yararlı bir `clean` alt komutu içerir.

#### 4.3.1. Silme kriterlerini belirleyin

Neyin silineceğini belirlemek için birden fazla [seçenek](https://www.nextflow.io/docs/latest/reference/cli.html#clean) vardır.

Burada size, çalıştırma adı kullanılarak belirtilen belirli bir çalıştırmadan önceki çalıştırmalardan tüm alt dizinleri silen bir örnek gösteriyoruz.

`-resume` kullanmadığınız en son başarılı çalıştırmayı arayın; bizim durumumuzda çalıştırma adı `golden_cantor` idi.

Çalıştırma adı, `Launching (...)` konsol çıktı satırında köşeli parantez içinde gösterilen makine tarafından oluşturulan iki parçalı dizedir.
Zaman damgasına ve/veya komut satırına göre bir çalıştırmayı aramak için Nextflow günlüğünü de kullanabilirsiniz.

#### 4.3.2. Kuru çalıştırma yapın

İlk olarak, komut verildiğinde neyin silineceğini kontrol etmek için kuru çalıştırma bayrağı `-n` kullanırız:

```bash
nextflow clean -before golden_cantor -n
```

??? success "Komut çıktısı"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Çıktınız farklı görev dizin adlarına sahip olacak ve farklı sayıda satıra sahip olabilir, ancak örneğe benzer görünmelidir.

Herhangi bir satır çıktısı görmüyorsanız, geçerli bir çalıştırma adı sağlamadınız veya silinecek geçmiş çalıştırma yok. Örnek komuttaki `golden_cantor`'u günlüğünüzdeki ilgili en son çalıştırma adıyla değiştirdiğinizden emin olun.

#### 4.3.3. Silme işlemine devam edin

Çıktı beklendiği gibi görünüyorsa ve silme işlemine devam etmek istiyorsanız, komutu `-n` yerine `-f` bayrağıyla yeniden çalıştırın:

```bash
nextflow clean -before golden_cantor -f
```

??? success "Komut çıktısı"

    ```console
    Removed /workspaces/training/hello-nextflow/work/a3/7be2fad5e71e5f49998f795677fd68
    ```

Çıktı daha öncekine benzer olmalıdır, ancak şimdi 'Would remove' yerine 'Removed' diyor.
Bunun iki karakterli alt dizinleri (yukarıdaki `a3/` gibi) kaldırmadığını ancak içeriklerini boşalttığını unutmayın.

!!! Warning

    Geçmiş çalıştırmalardan work alt dizinlerini silmek, bunları Nextflow'un önbelleğinden kaldırır ve bu dizinlerde depolanan tüm çıktıları siler.
    Bu, Nextflow'un ilgili süreçleri yeniden çalıştırmadan yürütmeyi sürdürme yeteneğini bozduğu anlamına gelir.

    Önemsediğiniz veya güvenmeyi planladığınız çıktıları kaydetmekten siz sorumlusunuz! Bu, `publish` yönergesi için `symlink` modu yerine `copy` modunu kullanmayı tercih etmemizin ana nedenidir.

### Özet

Çıktıları belirli bir dizine nasıl yayınlayacağınızı, zaten aynı şekilde çalıştırılmış adımları tekrarlamadan bir pipeline'ı nasıl yeniden başlatacağınızı ve eski work dizinlerini temizlemek için `nextflow clean` komutunu nasıl kullanacağınızı biliyorsunuz.

Daha genel olarak, basit bir Nextflow iş akışını nasıl yorumlayacağınızı, yürütmesini nasıl yöneteceğinizi ve çıktıları nasıl alacağınızı biliyorsunuz.

### Sırada ne var?

Küçük bir mola verin, hak ettiniz!

Hazır olduğunuzda, Nextflow'un yerleşik veri akışı paralelliğinden ve diğer güçlü özelliklerden yararlanmanıza olanak tanıyacak girdileri iş akışınıza beslemek için kanalların nasıl kullanılacağını öğrenmek için [**Bölüm 2: Hello Channels**](./02_hello_channels.md)'a geçin.

---

## Quiz

<quiz>
Bir Nextflow sürecinin minimum gerekli bileşenleri nelerdir?
- [ ] Yalnızca input ve output blokları
- [x] Output ve script blokları
- [ ] Input, output ve script blokları
- [ ] Yalnızca bir script bloğu

Daha fazla bilgi: [1.1.1. Process tanımı](#111-process-tanımı)
</quiz>

<quiz>
Bir süreçteki output bloğunun amacı nedir?
- [ ] Sonuçları konsola yazdırmak
- [ ] Dosyaları work dizinine kaydetmek
- [x] Süreçten beklenen çıktıları bildirmek
- [ ] Ortam değişkenlerini tanımlamak

Daha fazla bilgi: [1.1.1. Process tanımı](#111-process-tanımı)
</quiz>

<quiz>
Bir Nextflow iş akışını çalıştırmak için hangi komut kullanılır?
- [ ] `nextflow start`
- [ ] `nextflow execute`
- [x] `nextflow run`
- [ ] `nextflow launch`
</quiz>

<quiz>
Bir görevin work dizinine bakıldığında, yürütülen gerçek komutu hangi dosya içerir?

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
- [ ] İş akışını baştan yeniden başlatır
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

Daha fazla bilgi: [2.3. Yayınlama modunu kopyalama olarak ayarlayın](#23-yayınlama-modunu-kopyalama-olarak-ayarlayın)
</quiz>

<quiz>
Komut satırından bir Nextflow iş akışına bir parametre değeri nasıl iletilir?
- [ ] `-parameter value`
- [ ] `--parameter:value`
- [x] `--parameter value`
- [ ] `-p parameter=value`

Daha fazla bilgi: [3.2. Kullanıcı girdisini yakalamak için bir komut satırı parametresi ayarlayın](#32-kullanıcı-girdisini-yakalamak-için-bir-komut-satırı-parametresi-ayarlayın)
</quiz>

<quiz>
Bir Nextflow script bloğu içinde bir değişkene nasıl başvurulur?
- [ ] `%variable%` sözdizimini kullan
- [x] `#!groovy ${variable}` sözdizimini kullan
- [ ] `{{variable}}` sözdizimini kullan
- [ ] `[variable]` sözdizimini kullan
</quiz>
