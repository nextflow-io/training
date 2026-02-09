# Bölüm 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube kanalında izleyin.

:green_book: Video transkripti [burada](./transcripts/02_hello_channels.md) mevcuttur.
///

Bu kursun Bölüm 1'inde (Hello World), bir sürece değişken girdi sağlamanın nasıl yapılacağını, girdiyi doğrudan süreç çağrısında vererek gösterdik: `sayHello(params.input)`.
Bu kasıtlı olarak basitleştirilmiş bir yaklaşımdı.
Pratikte, bu yaklaşımın önemli sınırlamaları vardır; yalnızca süreci bir kez, tek bir değer üzerinde çalıştırmak istediğimiz çok basit durumlar için işe yarar.
Çoğu gerçekçi iş akışı kullanım durumunda, birden fazla değeri (örneğin, birden fazla örnek için deneysel veriler) işlemek isteriz, bu nedenle girdileri işlemek için daha sofistike bir yola ihtiyacımız vardır.

Nextflow [**kanalları**](https://nextflow.io/docs/latest/channel.html) bunun için vardır.
Kanallar, girdileri verimli bir şekilde işlemek ve çok adımlı iş akışlarında bir adımdan diğerine aktarmak için tasarlanmış kuyruklardır; yerleşik paralellik ve birçok ek fayda sağlarlar.

Kursun bu bölümünde, çeşitli farklı kaynaklardan gelen birden fazla girdiyi işlemek için bir kanalın nasıl kullanılacağını öğreneceksiniz.
Ayrıca kanal içeriklerini gerektiği gibi dönüştürmek için [**operatörleri**](https://nextflow.io/docs/latest/reference/operator.html) kullanmayı öğreneceksiniz.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Hello Nextflow](./index.md) kursunun Bölüm 1'ini tamamladığınızı varsayar, ancak o bölümde ele alınan temel konularda rahat hissediyorsanız, özel bir şey yapmadan buradan başlayabilirsiniz.

---

## 0. Isınma: `hello-channels.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-channels.nf` iş akışı betiğini kullanacağız.
Bu, bu eğitim kursunun Bölüm 1'ini tamamlayarak üretilen betiğe eşdeğerdir, ancak çıktı hedefini değiştirdik:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Her şeyin çalıştığından emin olmak için, herhangi bir değişiklik yapmadan önce betiği bir kez çalıştırın:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Daha önce olduğu gibi, `output.txt` adlı çıktı dosyasını `results/hello_channels` dizininde bulacaksınız (yukarıda gösterilen iş akışı betiğinin `output` bloğunda belirtildiği gibi).

??? abstract "Dizin içeriği"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Dosya içeriği"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Bu sizin için çalıştıysa, kanallar hakkında bilgi edinmeye hazırsınız.

---

## 1. Bir kanal aracılığıyla açıkça değişken girdiler sağlayın

`sayHello()` sürecine değişken girdiyi aktarmak için örtük işlemeye güvenmek yerine bir **kanal** oluşturacağız; bu yaklaşımın belirli sınırlamaları vardır.

### 1.1. Bir girdi kanalı oluşturun

Bir kanal kurmak için kullanabileceğimiz çeşitli [**kanal fabrikaları**](https://nextflow.io/docs/latest/reference/channel.html) vardır.
Şimdilik işleri basit tutmak için, [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of) adlı en temel kanal fabrikasını kullanacağız; bu, tek bir değer içeren bir kanal oluşturacaktır.
İşlevsel olarak bu, daha önce kurduğumuz şekle benzer olacak, ancak Nextflow'un örtük olarak bir kanal oluşturmasını sağlamak yerine, bunu şimdi açıkça yapıyoruz.

Kullanacağımız kod satırı şudur:

```console title="Sözdizimi"
greeting_ch = channel.of('Hello Channels!')
```

Bu, `channel.of()` kanal fabrikasını kullanarak `greeting_ch` adlı bir kanal oluşturur; bu, basit bir kuyruk kanalı kurar ve selamlama değeri olarak kullanılacak `'Hello Channels!'` dizesini yükler.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note

    Okunabilirlik adına geçici olarak CLI parametresi kullanmak yerine sabit kodlanmış dizelere geri dönüyoruz. Kanal düzeyinde neler olduğunu ele aldıktan sonra CLI parametrelerini kullanmaya geri döneceğiz.

Workflow bloğunda, kanal fabrikası kodunu ekleyin:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello Channels!')
        // bir selamlama yayınla
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // bir selamlama yayınla
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Bu henüz işlevsel değil çünkü girdiyi süreç çağrısına henüz değiştirmedik.

### 1.2. Kanalı süreç çağrısına girdi olarak ekleyin

Şimdi yeni oluşturduğumuz kanalı `sayHello()` süreç çağrısına gerçekten bağlamamız gerekiyor; daha önce doğrudan sağladığımız CLI parametresini değiştirerek.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello Channels!')
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello Channels!')
        // bir selamlama yayınla
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Bu, Nextflow'a `sayHello` sürecini `greeting_ch` kanalının içeriği üzerinde çalıştırmasını söyler.

Artık iş akışımız düzgün bir şekilde işlevseldir; `sayHello('Hello Channels!')` yazmanın açık eşdeğeridir.

### 1.3. İş akışını çalıştırın

Hadi çalıştıralım!

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Her iki düzenlemeyi de doğru yaptıysanız, başarılı bir yürütme elde etmelisiniz.
Sonucun hala daha önce olduğu gibi aynı olduğundan emin olmak için sonuçlar dizinini kontrol edebilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Böylece iş akışımızın esnekliğini artırdık ve aynı sonuca ulaştık.
Bu, somut bir fayda olmadan daha fazla kod yazmak gibi görünebilir, ancak değer, daha fazla girdiyi işlemeye başlar başlamaz netleşecektir.

Bunun bir önizlemesi olarak, devam etmeden önce bir şeye daha bakalım: veri girdisini yönetmek için açık bir kanal kullanmanın küçük ama kullanışlı bir faydası.

### 1.4. Kanal içeriğini incelemek için `view()` kullanın

Nextflow kanalları, içeriklerini operatörler kullanarak işlememize olanak tanıyacak şekilde oluşturulmuştur; bunları bu bölümde daha sonra ayrıntılı olarak ele alacağız.

Şimdilik, bir kanalın içeriğini incelemek için [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) adlı süper basit bir operatörün nasıl kullanılacağını göstereceğiz.
`view()`'ı Python'daki `print()` ifadesi veya diğer dillerdeki eşdeğeri gibi bir hata ayıklama aracı olarak düşünebilirsiniz.

Workflow bloğuna bu küçük satırı ekleyin:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello Channels!')
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Tam boşluk miktarı, 4'ün katı olduğu sürece önemli değildir; sadece `.view()` ifadesinin başlangıcını kanal yapısının `.of()` kısmına hizalamayı amaçlıyoruz.

Şimdi iş akışını tekrar çalıştırın:

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Gördüğünüz gibi, bu kanal içeriğini konsola çıktı olarak verir.
Burada sadece bir öğemiz var, ancak bir sonraki bölümde kanala birden fazla değer yüklemeye başladığımızda, bunun satır başına bir öğe çıktısı verecek şekilde ayarlandığını göreceksiniz.

### Özet

Bir sürece girdi sağlamak için temel bir kanal fabrikasının nasıl kullanılacağını biliyorsunuz.

### Sırada ne var?

İş akışını birden fazla girdi değeri üzerinde yineleme yapacak şekilde nasıl değiştireceğinizi öğrenin.

---

## 2. İş akışını birden fazla girdi değeri üzerinde çalışacak şekilde değiştirin

İş akışları genellikle toplu olarak işlenmesi amaçlanan girdi grupları üzerinde çalışır, bu nedenle iş akışını birden fazla girdi değerini kabul edecek şekilde yükseltmek istiyoruz.

### 2.1. Girdi kanalına birden fazla selamlama yükleyin

Kullanışlı bir şekilde, kullandığımız `channel.of()` kanal fabrikası birden fazla değeri kabul etmekten oldukça memnundur, bu nedenle bunu hiç değiştirmemize gerek yoktur.
Sadece kanala birden fazla değer yükleyebiliriz.

Bunları `'Hello'`, `'Bonjour'` ve `'Holà'` yapalım.

#### 2.1.1. Daha fazla selamlama ekleyin

Workflow bloğundan önce, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // girdiler için bir kanal oluştur
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // girdiler için bir kanal oluştur
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

Dokümantasyon bize bunun işe yaraması gerektiğini söylüyor. Gerçekten bu kadar basit olabilir mi?

#### 2.1.2. Komutu çalıştırın ve günlük çıktısına bakın

Hadi deneyelim.

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Kesinlikle sorunsuz çalışmış gibi görünüyor.
Yürütme monitörü, `sayHello` süreci için `3 of 3` çağrı yapıldığını gösteriyor ve `view()` ifadesi tarafından numaralandırılan üç selamlamayı görüyoruz, vaat edildiği gibi satır başına bir tane.

Ancak, sonuçlar dizininde hala sadece bir çıktı var:

??? abstract "Dizin içeriği"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Dosya içeriği"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Orada üç selamlamadan birini görmelisiniz, ancak aldığınız selamlama burada gösterilenden farklı olabilir.
Bunun neden olabileceğini düşünebiliyor musunuz?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Diyagramda, kanal yeşil renkte temsil edilir ve öğelerin sırası bir borudaki bilye gibi temsil edilir: yüklenen ilk öğe sağda, sonra ikincisi ortada, sonra üçüncüsü soldadır._

Yürütme monitörüne geri baktığımızda, bize sadece bir alt dizin yolu verdi (`f4/c9962c`).
Oraya bir göz atalım.

??? abstract "Dizin içeriği"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Dosya içeriği"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Bu, sonuçlar dizininde aldığımız selamlama bile değil! Neler oluyor?

Bu noktada, size varsayılan olarak ANSI günlük sisteminin aynı sürece yapılan birden fazla çağrıdan gelen günlükleri aynı satıra yazdığını söylememiz gerekiyor.
Yani sayHello() sürecine yapılan üç çağrının durumu da aynı yere iniyor.

Neyse ki, tam süreç çağrıları listesini görmek için bu davranışı devre dışı bırakabiliriz.

#### 2.1.3. Komutu `-ansi-log false` seçeneği ile tekrar çalıştırın

Günlüğü süreç çağrısı başına bir satır görüntüleyecek şekilde genişletmek için komuta `-ansi-log false` ekleyin.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Bu sefer üç süreç çalıştırmasının tümünü ve ilişkili çalışma alt dizinlerini çıktıda listelenen şekilde görüyoruz.

Bu çok daha iyi, en azından basit bir iş akışı için.
Karmaşık bir iş akışı veya çok sayıda girdi için, terminale çıktı olarak verilen tam listenin olması biraz bunaltıcı olurdu.
Bu yüzden `-ansi-log false` varsayılan davranış değildir.

!!! tip

    Durumun raporlanma şekli iki günlük modu arasında biraz farklıdır.
    Yoğunlaştırılmış modda, Nextflow çağrıların başarıyla tamamlanıp tamamlanmadığını bildirir.
    Bu genişletilmiş modda, yalnızca gönderildiklerini bildirir.

Her neyse, artık her süreç çağrısının alt dizinlerine sahip olduğumuza göre, günlüklerini ve çıktılarını arayabiliriz.

??? abstract "Dizin içeriği"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Dosya içeriği"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Bu, üç sürecin de başarıyla çalıştığını gösterir (yaşasın).

Bununla birlikte, sonuçlar dizininde hala sadece bir çıktı dosyası olması sorunu var.

`sayHello` süreci için çıktı dosyası adını sabit kodladığımızı hatırlayabilirsiniz, bu nedenle üç çağrının tümü `output.txt` adlı bir dosya üretti.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Çıktı dosyaları çalışma alt dizinlerinde kaldığı ve diğer süreçlerden izole edildiği sürece, bu sorun değildir.
Ancak aynı sonuçlar dizinine yayınlandıklarında, oraya ilk kopyalanan hangisiyse bir sonraki tarafından üzerine yazılır ve bu böyle devam eder.

### 2.2. Çıktı dosyası adlarının benzersiz olacağından emin olun

Tüm çıktıları aynı sonuçlar dizinine yayınlamaya devam edebiliriz, ancak benzersiz adlara sahip olacaklarından emin olmamız gerekir.
Özellikle, son dosya adlarının benzersiz olması için ilk süreci dinamik olarak bir dosya adı oluşturacak şekilde değiştirmemiz gerekir.

Peki dosya adlarını nasıl benzersiz yaparız?
Bunu yapmanın yaygın bir yolu, (girdi kanalından alınan) girdilerden benzersiz bir meta veri parçasını çıktı dosyası adının bir parçası olarak kullanmaktır.
Burada, kolaylık sağlamak için, sadece kısa bir dize olduğu için selamlamanın kendisini kullanacağız ve bunu temel çıktı dosya adının önüne ekleyeceğiz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Dinamik bir çıktı dosyası adı oluşturun

Süreç bloğunda, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

Hem çıktı tanımında hem de `script:` komut bloğunda `output.txt`'yi değiştirdiğinizden emin olun.

!!! tip

    Çıktı tanımında, çıktı dosya adı ifadesinin etrafında çift tırnak kullanMALIsınız (tek tırnak DEĞİL), aksi takdirde başarısız olur.

Bu, süreç her çağrıldığında benzersiz bir çıktı dosyası adı üretmeli, böylece çıktı dizininde aynı sürece yapılan diğer çağrılardan gelen çıktılardan ayırt edilebilir.

#### 2.2.2. İş akışını çalıştırın

Hadi çalıştıralım. Varsayılan ANSI günlük ayarlarıyla çalıştırmaya geri döndüğümüzü unutmayın.

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Özet görünümüne geri dönerek, çıktı tekrar bir satırda özetleniyor.
Tüm çıktı selamlamalarının orada olup olmadığını görmek için `results` dizinine bir göz atın.

??? abstract "Dizin içeriği"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Evet! Ve her birinin beklenen içeriği var.

??? abstract "Dosya içeriği"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Başarı! Artık çıktı dosyalarının üzerine yazılması konusunda endişelenmeden istediğimiz kadar selamlama ekleyebiliriz.

!!! tip

    Pratikte, dosyaları girdi verilerinin kendisine göre adlandırmak neredeyse her zaman pratik değildir.
    Dinamik dosya adları oluşturmanın daha iyi yolu, girdi dosyalarıyla birlikte bir sürece meta veri aktarmaktır.
    Meta veriler genellikle bir 'örnek sayfası' veya eşdeğerleri aracılığıyla sağlanır.
    Bunu daha sonra Nextflow eğitiminizde nasıl yapacağınızı öğreneceksiniz (bkz. [Meta veri yan görevi](../side_quests/metadata.md)).

### Özet

Bir kanal aracılığıyla birden fazla girdi öğesini nasıl besleyeceğinizi biliyorsunuz.

### Sırada ne var?

Bir kanalın içeriğini dönüştürmek için bir operatör kullanmayı öğrenin.

---

## 3. Bir dizi aracılığıyla birden fazla girdi sağlayın

Size az önce doğrudan kanal fabrikasında sabit kodlanmış birden fazla girdi öğesinin nasıl işleneceğini gösterdik.
Ya bu birden fazla girdiyi farklı bir şekilde sağlamak isteseydik?

Örneğin, şöyle bir öğe dizisi içeren bir girdi değişkeni kurduğumuzu hayal edin:

`greetings_array = ['Hello','Bonjour','Holà']`

Bunu çıktı kanalımıza yükleyebilir ve çalışmasını bekleyebilir miyiz?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Hadi öğrenelim.

### 3.1. Kanala girdi olarak bir değer dizisi sağlayın

Sağduyu, tek bir değer yerine basitçe bir değer dizisi geçirebilmemiz gerektiğini öne sürüyor.
Hadi deneyelim; girdi değişkenini kurmamız ve kanal fabrikasına yüklememiz gerekecek.

#### 3.1.1. Girdi değişkenini kurun

Az önce hayal ettiğimiz `greetings_array` değişkenini alıp workflow bloğuna ekleyerek gerçeğe dönüştürelim:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Bu henüz işlevsel değil, sadece dizi için bir bildirim ekledik.

#### 3.1.2. Selamlamalar dizisini kanal fabrikasına girdi olarak ayarlayın

Şimdi kanal fabrikasında şu anda sabit kodlanmış olan `'Hello','Bonjour','Holà'` değerlerini az önce oluşturduğumuz `greetings_array` ile değiştireceğiz.

Workflow bloğunda, aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Bu artık işlevsel olmalı.

#### 3.1.3. İş akışını çalıştırın

Hadi çalıştırmayı deneyelim:

```bash
nextflow run hello-channels.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Ah hayır! Bir hata var!

`view()` çıktısına ve hata mesajlarına bakın.

Nextflow'un dizideki üç dizeyi ayrı değerler olarak kullanmak yerine, `[Hello, Bonjour, Holà]`'yı tek bir dize değeri olarak kullanarak tek bir süreç çağrısı çalıştırmaya çalıştığı görülüyor.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Yani soruna neden olan 'paketleme'dir.
Nextflow'un diziyi açmasını ve bireysel dizeleri kanala yüklemesini nasıl sağlarız?

### 3.2. Kanal içeriğini dönüştürmek için bir operatör kullanın

[**Operatörlerin**](https://nextflow.io/docs/latest/reference/operator.html) devreye girdiği yer burasıdır.
Zaten `.view()` operatörünü kullandınız, bu sadece orada ne olduğuna bakar.
Şimdi bir kanalın içeriği üzerinde işlem yapmamıza izin veren operatörlere bakacağız.

Nextflow dokümantasyonundaki [operatörler listesine](https://nextflow.io/docs/latest/reference/operator.html) göz atarsanız, tam olarak ihtiyacımız olanı yapan [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten)'i bulacaksınız: bir dizinin içeriğini açmak ve bunları bireysel öğeler olarak yaymak.

#### 3.2.1. `flatten()` operatörünü ekleyin

`flatten()` operatörünü girdi kanalımıza uygulamak için, onu kanal fabrikası bildirimine ekliyoruz.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Burada okunabilirlik için operatörü bir sonraki satıra ekledik, ancak isterseniz operatörleri kanal fabrikasıyla aynı satıra ekleyebilirsiniz, şöyle:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. `view()` ifadelerini iyileştirin

Bunu hemen çalışıp çalışmadığını test etmek için çalıştırabiliriz, ancak bunu yaparken, kanal içeriğini nasıl incelediğimizi iyileştireceğiz.

`flatten()` operatörü uygulanmadan önce ve sonra içeriklerin nasıl göründüğünü karşılaştırabilmek istiyoruz, bu nedenle ikinci bir tane ekleyeceğiz VE çıktıda daha net etiketlenmelerini sağlamak için biraz kod ekleyeceğiz.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

İkinci bir `.view` ifadesi eklediğimizi ve her biri için boş parantezleri (`()`) `{ greeting -> "Before flatten: $greeting" }` gibi bazı kodlar içeren süslü parantezlerle değiştirdiğimizi görüyorsunuz.

Bunlara _closure_ denir. İçerdikleri kod, kanaldaki her öğe için yürütülecektir.
İç değer için geçici bir değişken tanımlıyoruz, burada `greeting` olarak adlandırılmış (ancak herhangi bir rastgele ad olabilir), bu yalnızca o closure'ın kapsamı içinde kullanılır.

Bu örnekte, `$greeting` kanala yüklenen her bir öğeyi temsil eder.
Bu, düzgün etiketlenmiş konsol çıktısı ile sonuçlanacaktır.

!!! info

    Bazı pipeline'larda operatör closure'ları içinde kullanılan `$it` adlı özel bir değişken görebilirsiniz.
    Bu, `->` ile tanımlamaya gerek kalmadan iç değişkene kısa yoldan erişim sağlayan _örtük_ bir değişkendir.

    Kod netliğine yardımcı olmak için açık olmayı tercih ediyoruz, bu nedenle `$it` sözdizimi önerilmez ve Nextflow dilinden yavaş yavaş kaldırılacaktır.

#### 3.2.3. İş akışını çalıştırın

Son olarak, iş akışını tekrar çalıştırmayı deneyebilirsiniz!

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

Bu sefer çalışıyor VE bize `flatten()` operatörünü çalıştırmadan önce ve sonra kanalın içeriğinin nasıl göründüğüne dair ek içgörü veriyor.

- Tek bir `Before flatten:` ifadesi çünkü o noktada kanal bir öğe içeriyor, orijinal dizi.
- Üç ayrı `After flatten:` ifadesi, her selamlama için bir tane, bunlar artık kanaldaki bireysel öğelerdir.

Önemlisi, bu her öğenin artık iş akışı tarafından ayrı ayrı işlenebileceği anlamına gelir.

!!! tip

    Teknik olarak, işleminde örtük bir eşleme adımı içeren farklı bir kanal fabrikası olan [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist) kullanarak aynı sonuçları elde etmek mümkündür.
    Burada basit bir kullanım durumunda bir operatörün kullanımını göstermek için bunu kullanmamayı seçtik.

### Özet

Bir kanalın içeriğini dönüştürmek için `flatten()` gibi bir operatörün nasıl kullanılacağını ve bir operatör uygulamadan önce ve sonra kanal içeriğini incelemek için `view()` operatörünün nasıl kullanılacağını biliyorsunuz.

### Sırada ne var?

İş akışının girdi değerlerinin kaynağı olarak bir dosya almasını nasıl sağlayacağınızı öğrenin.

---

## 4. Bir CSV dosyasından girdi değerlerini okuyun

Gerçekçi olarak, bir değer dizisinden başlayacağımız durumlar nadiren olur.
Büyük olasılıkla, işlenmesi gereken verileri içeren bir veya daha fazla dosyamız olacak, bir tür yapılandırılmış formatta.

Birkaç girdi selamlaması içeren `greetings.csv` adlı bir CSV dosyası hazırladık; bu, gerçek bir veri analizinde işlemek isteyebileceğiniz sütunlu veri türünü taklit ediyor, `data/` altında saklanıyor.
(Sayılar anlamlı değildir, sadece açıklayıcı amaçlar için oradadırlar.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Bir sonraki görevimiz, iş akışımızı bu dosyadaki değerleri okuyacak şekilde uyarlamaktır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Bunun nasıl gerçekleşebileceğini görelim.

### 4.1. Betiği selamlamaların kaynağı olarak bir CSV dosyası bekleyecek şekilde değiştirin

Başlamak için, betikte iki temel değişiklik yapmamız gerekecek:

- Girdi parametresini CSV dosyasına işaret edecek şekilde değiştirin
- Kanal fabrikasını bir dosyayı işlemek için tasarlanmış bir fabrikaya değiştirin

#### 4.1.1. Girdi parametresini CSV dosyasına işaret edecek şekilde değiştirin

Bölüm 1'de kurduğumuz `params.input` parametresini hatırlıyor musunuz?
Bunu selamlamalarımızı içeren CSV dosyasına işaret edecek şekilde güncelleyeceğiz.

Parametre bildirimine aşağıdaki düzenlemeyi yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline parametreleri
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parametreleri
     */
    input: String = 'Holà mundo!'
    ```

Bu, dosyanın iş akışı koduyla birlikte bulunduğunu varsayar.
Diğer veri konumlarıyla nasıl başa çıkacağınızı daha sonra Nextflow yolculuğunuzda öğreneceksiniz.

#### 4.1.2. Bir dosyayı işlemek için tasarlanmış bir kanal fabrikasına geçin

Artık girdi olarak basit dizeler yerine bir dosya kullanmak istediğimize göre, daha önceki `channel.of()` kanal fabrikasını kullanamayız.
Dosya yollarını işlemek için bazı yerleşik işlevlere sahip yeni bir kanal fabrikası olan [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath)'e geçmemiz gerekiyor.

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // girdi selamlamalarının bir dizisini bildir
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Kanal girdisini `param.input`'a geri değiştirdiğimizi ve artık ihtiyacımız olmayacağı için `greetings_array` bildirimini sildiğimizi fark edeceksiniz.
Ayrıca `flatten()` ve ikinci `view()` ifadesini yorum satırı haline getirdik.

#### 4.1.3. İş akışını çalıştırın

Yeni kanal fabrikası ve girdi dosyasıyla iş akışını çalıştırmayı deneyelim.

```bash
nextflow run hello-channels.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Ah hayır, çalışmıyor. Konsol çıktısının başına ve hata mesajına bakın.
`Command executed:` kısmı burada özellikle yararlıdır.

Bu biraz tanıdık gelebilir.
Nextflow'un dosya yolunun kendisini bir dize değeri olarak kullanarak tek bir süreç çağrısı çalıştırmaya çalıştığı görülüyor.
Yani dosya yolunu doğru bir şekilde çözdü, ancak aslında içeriğini ayrıştırmadı, ki biz bunu istiyorduk.

Nextflow'un dosyayı açmasını ve içeriğini kanala yüklemesini nasıl sağlarız?

Başka bir [operatöre](https://nextflow.io/docs/latest/reference/operator.html) ihtiyacımız var gibi görünüyor!

### 4.2. Dosyayı ayrıştırmak için `splitCsv()` operatörünü kullanın

Operatörler listesine tekrar baktığımızda, CSV biçimli metni ayrıştırmak ve bölmek için tasarlanmış [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv)'yi buluyoruz.

#### 4.2.1. Kanala `splitCsv()` uygulayın

Operatörü uygulamak için, onu daha önce olduğu gibi kanal fabrikası satırına ekliyoruz.

Workflow bloğunda, `flatten()`'i `splitcsv()` ile değiştirmek için aşağıdaki kod değişikliğini yapın (yorum satırı olmadan):

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Gördüğünüz gibi, önce/sonra `view()` ifadelerini de güncelledik.
Teknik olarak aynı değişken adını (`greeting`) kullanabilirdik, ancak kodu başkaları tarafından daha okunabilir hale getirmek için daha uygun bir şeye (`csv`) güncelledik.

#### 4.2.2. İş akışını tekrar çalıştırın

Eklenen CSV ayrıştırma mantığıyla iş akışını çalıştırmayı deneyelim.

```bash
nextflow run hello-channels.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

İlginç bir şekilde, bu da başarısız oluyor, ancak farklı bir hatayla.
Bu sefer Nextflow dosyanın içeriğini ayrıştırdı (yaşasın!) ancak her satırı bir dizi olarak yükledi ve her dizi kanaldaki bir öğedir.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Ona her satırdaki yalnızca ilk sütunu almasını söylememiz gerekiyor.
Peki bunu nasıl açarız?

Daha önce bir kanalın içeriğini açmak için `flatten()` kullanmıştık, ancak bu burada işe yaramaz çünkü flatten _her şeyi_ açar (isterseniz kendiniz deneyebilirsiniz).

Bunun yerine, Nextflow pipeline'larında gerçekten kullanışlı olan ve çok sık karşılaşılan `map()` adlı başka bir operatör kullanacağız.

### 4.3. Selamlamaları çıkarmak için `map()` operatörünü kullanın

[`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) operatörü, bir kanalın içeriğine her türlü eşleme yapmamıza olanak tanıyan çok kullanışlı küçük bir araçtır.

Bu durumda, veri dosyamızdaki her satırdan istediğimiz o bir öğeyi çıkarmak için kullanacağız.
Sözdizimi şöyle görünür:

```groovy title="Sözdizimi"
.map { row -> row[0] }
```

Bu, 'kanaldaki her satır için, içerdiği 0. (ilk) öğeyi al' anlamına gelir.

Öyleyse bunu CSV ayrıştırmamıza uygulayalım.

#### 4.3.1. Kanala `map()` uygulayın

Workflow bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Önce"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Operatörün beklediğimiz şeyi yaptığını doğrulamak için başka bir `view()` çağrısı eklediğimizi görüyorsunuz.

#### 4.3.2. İş akışını çalıştırın

Bunu bir kez daha çalıştıralım:

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

Bu sefer hatasız çalışmalı.

`view()` ifadelerinin çıktısına baktığınızda, şunları görürsünüz:

- Tek bir `Before splitCsv:` ifadesi: o noktada kanal bir öğe içeriyor, orijinal dosya yolu.
- Üç ayrı `After splitCsv:` ifadesi: her selamlama için bir tane, ancak her biri dosyadaki o satıra karşılık gelen bir dizi içinde yer alıyor.
- Üç ayrı `After map:` ifadesi: her selamlama için bir tane, bunlar artık kanaldaki bireysel öğelerdir.

_Satırların çıktınızda farklı bir sırada görünebileceğini unutmayın._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

Her selamlamanın doğru bir şekilde çıkarıldığını ve iş akışı boyunca işlendiğini doğrulamak için çıktı dosyalarına da bakabilirsiniz.

Daha önce olduğu gibi aynı sonuca ulaştık, ancak şimdi herhangi bir kodu değiştirmeden bir girdi dosyasını değiştirerek işlemek istediğimiz selamlamalar kanalına daha fazla öğe ekleme konusunda çok daha fazla esnekliğe sahibiz.
Karmaşık girdileri işlemek için daha sofistike yaklaşımları daha sonraki bir eğitimde öğreneceksiniz.

### Özet

Bir girdi değerleri dosyasını okumak ve bunları uygun şekilde işlemek için `.fromPath()` kanal yapıcısını ve `splitCsv()` ve `map()` operatörlerini nasıl kullanacağınızı biliyorsunuz.

Daha genel olarak, Nextflow'un süreçlere girdileri yönetmek için **kanalları** ve içeriklerini dönüştürmek için **operatörleri** nasıl kullandığına dair temel bir anlayışa sahipsiniz.
Ayrıca kanalların paralel yürütmeyi örtük olarak nasıl işlediğini gördünüz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### Sırada ne var?

Büyük bir mola verin, bu bölümde çok çalıştınız!

Hazır olduğunuzda, daha fazla adım eklemeyi ve bunları uygun bir iş akışına bağlamayı öğrenmek için [**Bölüm 3: Hello Workflow**](./03_hello_workflow.md)'a geçin.

---

## Test

<quiz>
Nextflow'da kanal nedir?
- [ ] Bir dosya yolu belirtimi
- [ ] Bir süreç tanımı
- [x] Süreçler arasında veri aktarmak için kuyruk benzeri bir yapı
- [ ] Bir yapılandırma ayarı

Daha fazla bilgi: [1.1. Bir girdi kanalı oluşturun](#11-bir-girdi-kanali-olusturun)
</quiz>

<quiz>
Bu kod ne çıktı verir?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (tek bir liste)
- [x] Her öğe ayrı bir satırda: `Hello`, `Bonjour`, `Hola`
- [ ] Hiçbir şey (kanallar varsayılan olarak yazdırılmaz)
- [ ] Bir hata (geçersiz sözdizimi)

Daha fazla bilgi: [1.1. Bir girdi kanalı oluşturun](#11-bir-girdi-kanali-olusturun)
</quiz>

<quiz>
Bir kanal birden fazla değer içerdiğinde, Nextflow süreç yürütmesini nasıl işler?
- [ ] Süreç tüm değerlerle bir kez çalışır
- [x] Süreç kanaldaki her değer için bir kez çalışır
- [ ] Süreç yalnızca ilk değerle çalışır
- [ ] Süreç yalnızca son değerle çalışır

Daha fazla bilgi: [2. İş akışını birden fazla girdi değeri üzerinde çalışacak şekilde değiştirin](#2-is-akisini-birden-fazla-girdi-degeri-uzerinde-calisacak-sekilde-degistirin)
</quiz>

<quiz>
`flatten()` operatörü ne yapar?
- [ ] Birden fazla kanalı birleştirir
- [ ] Kanal öğelerini sıralar
- [x] Dizileri bireysel öğelere açar
- [ ] Yinelenen öğeleri kaldırır

Daha fazla bilgi: [3.2.1. `flatten()` operatörünü ekleyin](#321-flatten-operatorunu-ekleyin)
</quiz>

<quiz>
`view()` operatörünün amacı nedir?
- [ ] Kanal içeriğini filtrelemek
- [ ] Kanal öğelerini dönüştürmek
- [x] Kanal içeriğini incelemek ve hata ayıklamak
- [ ] Kanal içeriğini bir dosyaya kaydetmek

Daha fazla bilgi: [1.4. Kanal içeriğini incelemek için `view()` kullanın](#14-kanal-icerigini-incelemek-icin-view-kullanin)
</quiz>

<quiz>
`splitCsv()` ne yapar?
- [ ] Kanal içeriğinden bir CSV dosyası oluşturur
- [ ] Bir dizeyi virgüllerle böler
- [x] Bir CSV dosyasını her satırı temsil eden dizilere ayrıştırır
- [ ] Birden fazla CSV dosyasını birleştirir

Daha fazla bilgi: [4.2. Dosyayı ayrıştırmak için `splitCsv()` operatörünü kullanın](#42-dosyayi-ayristirmak-icin-splitcsv-operatorunu-kullanin)
</quiz>

<quiz>
`map()` operatörünün amacı nedir?
- [ ] Bir kanaldan öğeleri filtrelemek
- [ ] Birden fazla kanalı birleştirmek
- [x] Bir kanaldaki her öğeyi dönüştürmek
- [ ] Bir kanaldaki öğeleri saymak

Daha fazla bilgi: [4.3. Selamlamaları çıkarmak için `map()` operatörünü kullanın](#43-selamlamalari-cikarmak-icin-map-operatorunu-kullanin)
</quiz>

<quiz>
Birden fazla girdiyi işlerken dinamik çıktı dosya adları kullanmak neden önemlidir?
- [ ] Performansı artırmak için
- [ ] Disk alanını azaltmak için
- [x] Çıktı dosyalarının birbirinin üzerine yazılmasını önlemek için
- [ ] Resume işlevselliğini etkinleştirmek için

Daha fazla bilgi: [2.2. Çıktı dosyası adlarının benzersiz olacağından emin olun](#22-cikti-dosyasi-adlarinin-benzersiz-olacagindan-emin-olun)
</quiz>
