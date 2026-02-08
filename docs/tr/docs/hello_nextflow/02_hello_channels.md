# Bölüm 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=tr" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube kanalında [tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) izleyin.

:green_book: Video transkripti [burada](./transcripts/02_hello_channels.md) mevcuttur.
///

Bu kursun 1. Bölümünde (Hello World), bir sürece değişken girdi sağlamak için girdiyi doğrudan süreç çağrısında nasıl ileteceğinizi gösterdik: `sayHello(params.input)`.
Bu kasıtlı olarak basitleştirilmiş bir yaklaşımdı.
Pratikte, bu yaklaşımın önemli sınırlamaları vardır; yani yalnızca süreci tek bir değer üzerinde yalnızca bir kez çalıştırmak istediğimiz çok basit durumlar için çalışır.
Çoğu gerçekçi iş akışı kullanım durumunda, birden fazla değeri işlemek istiyoruz (örneğin, birden fazla örnek için deneysel veriler), bu nedenle girdileri işlemek için daha sofistike bir yönteme ihtiyacımız var.

Nextflow [**kanalları**](https://nextflow.io/docs/latest/channel.html) tam da bunun için var.
Kanallar, girdileri verimli bir şekilde işlemek ve çok adımlı iş akışlarında bir adımdan diğerine taşımak için tasarlanmış kuyruklardır; yerleşik paralellik ve birçok ek avantaj sağlarlar.

Bu kursun bu bölümünde, çeşitli farklı kaynaklardan gelen birden fazla girdiyi işlemek için bir kanalı nasıl kullanacağınızı öğreneceksiniz.
Ayrıca kanal içeriklerini gerektiği gibi dönüştürmek için [**operatörleri**](https://nextflow.io/docs/latest/reference/operator.html) nasıl kullanacağınızı da öğreneceksiniz.

??? info "Bu bölümden nasıl başlanır"

    Bu kursun bu bölümü, [Hello Nextflow](./index.md) kursunun 1. Bölümünü tamamladığınızı varsayar, ancak o bölümde ele alınan temel konulara hakimseniz, özel bir şey yapmadan buradan başlayabilirsiniz.

---

## 0. Isınma: `hello-channels.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-channels.nf` iş akışı betiğini kullanacağız.
Bu betik, bu eğitim kursunun 1. Bölümünde üretilen betiğe eşdeğerdir; ancak çıktı hedefini değiştirdik:

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

Daha önce olduğu gibi, çıktı dosyasını `output.txt` adıyla `results/hello_channels` dizininde bulacaksınız (yukarıda gösterilen iş akışı betiğinin `output` bloğunda belirtildiği gibi).

??? abstract "Dizin içerikleri"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Dosya içerikleri"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Bu sizin için çalıştıysa, kanallar hakkında öğrenmeye hazırsınız.

---

## 1. Değişken girdileri bir kanal aracılığıyla açıkça sağlayın

Örtük işlemeye güvenmek yerine, değişken girdiyi `sayHello()` sürecine iletmek için bir **kanal** oluşturacağız; örtük işlemenin belirli sınırlamaları vardır.

### 1.1. Bir girdi kanalı oluşturun

Bir kanal kurmak için kullanabileceğimiz çeşitli [**kanal fabrikaları**](https://nextflow.io/docs/latest/reference/channel.html) vardır.
Şimdilik işleri basit tutmak için, tek bir değer içeren bir kanal oluşturacak olan [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of) adlı en temel kanal fabrikasını kullanacağız.
İşlevsel olarak bu, daha önce kurduğumuz yönteme benzer olacak, ancak Nextflow'un örtük olarak bir kanal oluşturmasına izin vermek yerine, bunu artık açıkça yapıyoruz.

Kullanacağımız kod satırı şu:

```console title="Sözdizimi"
greeting_ch = channel.of('Hello Channels!')
```

Bu, `channel.of()` kanal fabrikasını kullanarak `greeting_ch` adında bir kanal oluşturur; bu fabrika basit bir kuyruk kanalı kurar ve selamlama değeri olarak kullanılmak üzere `'Hello Channels!'` dizesini yükler.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note "Not"

    Okunabilirlik adına geçici olarak sabit kodlanmış dizelere dönüyoruz ve CLI parametresi kullanmıyoruz. Kanal düzeyinde neler olduğunu ele aldıktan sonra CLI parametrelerini kullanmaya geri döneceğiz.

İş akışı bloğunda, kanal fabrikası kodunu ekleyin:

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

Süreç çağrısına girdiyi henüz değiştirmediğimiz için bu henüz işlevsel değil.

### 1.2. Kanalı süreç çağrısına girdi olarak ekleyin

Şimdi yeni oluşturduğumuz kanalı `sayHello()` süreç çağrısına bağlamamız gerekiyor; daha önce doğrudan sağladığımız CLI parametresinin yerine geçecek.

İş akışı bloğunda, aşağıdaki kod değişikliğini yapın:

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

Bu, Nextflow'a `sayHello` sürecini `greeting_ch` kanalının içerikleri üzerinde çalıştırmasını söyler.

Artık iş akışımız düzgün bir şekilde işlevsel; `sayHello('Hello Channels!')` yazmanın açık eşdeğeridir.

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
Sonucun hâlâ daha önce olduğu gibi olduğundan emin olmak için results dizinini kontrol edebilirsiniz.

??? abstract "Dosya içerikleri"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Böylece aynı sonuca ulaşırken iş akışımızın esnekliğini artırdık.
Bu, somut bir fayda olmadan daha fazla kod yazmak gibi görünebilir, ancak daha fazla girdiyi işlemeye başladığımızda değeri netleşecek.

Bunun bir önizlemesi olarak, devam etmeden önce bir şeye daha bakalım: veri girişini yönetmek için açık bir kanal kullanmanın küçük ama kullanışlı bir faydası.

### 1.4. Kanal içeriklerini incelemek için `view()` kullanın

Nextflow kanalları, içerikleri üzerinde operatörler kullanarak işlem yapmamıza olanak tanıyan bir şekilde oluşturulmuştur; operatörleri bu bölümde daha ayrıntılı ele alacağız.

Şimdilik, size bir kanalın içeriklerini incelemek için [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) adlı çok basit bir operatörü nasıl kullanacağınızı göstereceğiz.
`view()`'ı Python'daki `print()` ifadesi veya diğer dillerdeki eşdeğeri gibi bir hata ayıklama aracı olarak düşünebilirsiniz.

Bu küçük satırı iş akışı bloğuna ekleyin:

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

Tam boşluk miktarı 4'ün katı olduğu sürece önemli değil; amacımız `.view()` ifadesinin başlangıcını kanal yapımının `.of()` kısmına hizalamak.

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

Gördüğünüz gibi, bu kanal içeriklerini konsola çıktı olarak verir.
Burada yalnızca bir öğemiz var, ancak bir sonraki bölümde kanala birden fazla değer yüklemeye başladığımızda, bunun satır başına bir öğe çıktısı verecek şekilde ayarlandığını göreceksiniz.

### Özet

Bir sürece girdi sağlamak için temel bir kanal fabrikasını nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

İş akışının birden fazla girdi değeri üzerinde yineleme yapmasını sağlamak için kanalları nasıl kullanacağınızı öğrenin.

---

## 2. İş akışını birden fazla girdi değeri üzerinde çalışacak şekilde değiştirin

İş akışları genellikle toplu olarak işlenmesi gereken girdi grupları üzerinde çalışır, bu nedenle iş akışını birden fazla girdi değerini kabul edecek şekilde yükseltmek istiyoruz.

### 2.1. Girdi kanalına birden fazla selamlama yükleyin

Kullanışlı bir şekilde, kullandığımız `channel.of()` kanal fabrikası birden fazla değeri kabul etmekten memnuniyet duyar, bu nedenle bunu hiç değiştirmemize gerek yok.
Sadece kanala birden fazla değer yükleyebiliriz.

Bunları `'Hello'`, `'Bonjour'` ve `'Holà'` yapalım.

#### 2.1.1. Daha fazla selamlama ekleyin

İş akışı bloğundan önce, aşağıdaki kod değişikliğini yapın:

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

Dokümantasyon bunun çalışması gerektiğini söylüyor. Gerçekten bu kadar basit olabilir mi?

#### 2.1.2. Komutu çalıştırın ve log çıktısına bakın

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

Kesinlikle sorunsuz çalışmış görünüyor.
Yürütme monitörü `sayHello` süreci için `3 of 3` çağrı yapıldığını gösteriyor ve söz verildiği gibi `view()` ifadesi tarafından sıralanan üç selamlamayı satır başına bir tane olarak görüyoruz.

Ancak, results dizininde hâlâ yalnızca bir çıktı var:

??? abstract "Dizin içerikleri"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Dosya içerikleri"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Orada üç selamlamadan birini görmelisiniz, ancak aldığınız burada gösterilenden farklı olabilir.
Bunun neden olabileceğini düşünebilir misiniz?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_Diyagramda, kanal yeşil renkte temsil edilmektedir ve öğelerin sırası bir borudaki misketler gibi temsil edilmektedir: ilk yüklenen sağda, sonra ikincisi ortada, sonra üçüncüsü soldadır._

Yürütme monitörüne geri baktığımızda, bize yalnızca bir alt dizin yolu verdi (`f4/c9962c`).
Hadi oraya bir göz atalım.

??? abstract "Dizin içerikleri"

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

??? abstract "Dosya içerikleri"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Bu results dizininde aldığımız selamlama bile değil! Neler oluyor?

Bu noktada, varsayılan olarak ANSI loglama sisteminin aynı sürece yapılan birden fazla çağrının loglamasını aynı satıra yazdığını söylememiz gerekiyor.
Dolayısıyla sayHello() sürecine yapılan üç çağrının tamamının durumu aynı yere düşüyor.

Neyse ki, süreç çağrılarının tam listesini görmek için bu davranışı devre dışı bırakabiliriz.

#### 2.1.3. Komutu `-ansi-log false` seçeneğiyle tekrar çalıştırın

Loglamayı süreç çağrısı başına bir satır görüntüleyecek şekilde genişletmek için komuta `-ansi-log false` ekleyin.

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

Bu sefer çıktıda listelenen üç süreç çalışmasını ve bunlarla ilişkili work alt dizinlerini görüyoruz.

Bu çok daha iyi, en azından basit bir iş akışı için.
Karmaşık bir iş akışı veya çok sayıda girdi için, tam listenin terminale çıktı olarak verilmesi biraz bunaltıcı olurdu.
Bu nedenle `-ansi-log false` varsayılan davranış değil.

!!! tip "İpucu"

    Durumun raporlanma şekli iki loglama modu arasında biraz farklıdır.
    Yoğunlaştırılmış modda, Nextflow çağrıların başarıyla tamamlanıp tamamlanmadığını raporlar.
    Bu genişletilmiş modda, yalnızca gönderildiklerini raporlar.

Her neyse, artık her süreç çağrısının alt dizinlerine sahip olduğumuza göre, loglarına ve çıktılarına bakabiliriz.

??? abstract "Dizin içerikleri"

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

??? abstract "Dosya içerikleri"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Bu, üç sürecin de başarıyla çalıştığını gösteriyor (yaşasın).

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

Bununla birlikte, results dizininde hâlâ yalnızca bir çıktı dosyası olması sorunu var.

`sayHello` süreci için çıktı dosya adını sabit kodladığımızı hatırlayabilirsiniz, bu nedenle üç çağrının tümü `output.txt` adlı bir dosya üretti.

Çıktı dosyaları, diğer süreçlerden izole edilmiş work alt dizinlerinde kaldığı sürece, bu sorun değil.
Ancak aynı results dizinine yayınlandıklarında, ilk kopyalanan bir sonraki tarafından üzerine yazılır ve bu böyle devam eder.

### 2.2. Çıktı dosya adlarının benzersiz olmasını sağlayın

Tüm çıktıları aynı results dizinine yayınlamaya devam edebiliriz, ancak benzersiz adlara sahip olmalarını sağlamamız gerekiyor.
Özellikle, son dosya adlarının benzersiz olması için ilk süreci dinamik olarak bir dosya adı oluşturacak şekilde değiştirmemiz gerekiyor.

Peki dosya adlarını nasıl benzersiz yaparız?
Bunu yapmanın yaygın bir yolu, girdilerden (girdi kanalından alınan) bazı benzersiz meta verileri çıktı dosya adının bir parçası olarak kullanmaktır.
Burada, kolaylık olsun diye, yalnızca kısa bir dize olduğu için selamlamanın kendisini kullanacağız ve bunu temel çıktı dosya adının önüne ekleyeceğiz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. Dinamik bir çıktı dosya adı oluşturun

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

`output.txt` ifadesini hem çıktı tanımında hem de `script:` komut bloğunda değiştirdiğinizden emin olun.

!!! tip "İpucu"

    Çıktı tanımında, çıktı dosya adı ifadesinin etrafında çift tırnak kullanmalısınız (tek tırnak DEĞİL), aksi takdirde başarısız olacaktır.

Bu, süreç her çağrıldığında benzersiz bir çıktı dosya adı üretmelidir, böylece çıktı dizinindeki aynı sürece yapılan diğer çağrıların çıktılarından ayırt edilebilir.

#### 2.2.2. İş akışını çalıştırın

Hadi çalıştıralım. Varsayılan ANSI log ayarlarıyla çalıştırmaya geri döndüğümüzü unutmayın.

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

Özet görünümüne geri döndüğümüzde, çıktı tekrar bir satırda özetleniyor.
Tüm çıktı selamlamalarının orada olup olmadığını görmek için `results` dizinine bakın.

??? abstract "Dizin içerikleri"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Evet! Ve her birinin beklenen içerikleri var.

??? abstract "Dosya içerikleri"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Başarılı! Artık çıktı dosyalarının üzerine yazılması konusunda endişelenmeden istediğimiz kadar selamlama ekleyebiliriz.

!!! tip "İpucu"

    Pratikte, dosyaları girdi verilerine dayalı olarak adlandırmak neredeyse her zaman pratik değildir.
    Dinamik dosya adları oluşturmanın daha iyi yolu, meta verileri girdi dosyalarıyla birlikte bir sürece iletmektir.
    Meta veriler genellikle bir 'örnek sayfası' veya eşdeğerleri aracılığıyla sağlanır.
    Bunu daha sonra Nextflow eğitiminizde öğreneceksiniz (bkz. [Meta veri yan görevi](../side_quests/metadata.md)).

### Özet

Bir kanal aracılığıyla birden fazla girdi öğesini nasıl besleyeceğinizi biliyorsunuz.

### Sırada ne var?

Bir kanalın içeriklerini dönüştürmek için bir operatörü nasıl kullanacağınızı öğrenin.

---

## 3. Bir dizi aracılığıyla birden fazla girdi sağlayın

Az önce, doğrudan kanal fabrikasında sabit kodlanmış birden fazla girdi öğesini nasıl işleyeceğinizi gösterdik.
Ya bu birden fazla girdiyi farklı bir şekilde sağlamak istersek?

Örneğin, şöyle bir öğe dizisi içeren bir girdi değişkeni kurduğumuzu düşünün:

`greetings_array = ['Hello','Bonjour','Holà']`

Bunu çıktı kanalımıza yükleyip çalışmasını bekleyebilir miyiz?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Hadi öğrenelim.

### 3.1. Kanala girdi olarak bir değer dizisi sağlayın

Sağduyu, tek bir değer yerine bir değer dizisini basitçe iletebilmemiz gerektiğini öne sürüyor.
Hadi deneyelim; girdi değişkenini kurmamız ve kanal fabrikasına yüklememiz gerekecek.

#### 3.1.1. Girdi değişkenini kurun

Az önce hayal ettiğimiz `greetings_array` değişkenini iş akışı bloğuna ekleyerek gerçeğe dönüştürelim:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // bir girdi selamlama dizisi tanımla
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

Bu henüz işlevsel değil, sadece dizi için bir tanımlama ekledik.

#### 3.1.2. Selamlama dizisini kanal fabrikasına girdi olarak ayarlayın

Şimdi kanal fabrikasında şu anda sabit kodlanmış `'Hello','Bonjour','Holà'` değerlerini az önce oluşturduğumuz `greetings_array` ile değiştireceğiz.

İş akışı bloğunda, aşağıdaki değişikliği yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // bir girdi selamlama dizisi tanımla
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
        // bir girdi selamlama dizisi tanımla
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

Hay aksi! Bir hata var!

`view()` çıktısına ve hata mesajlarına bakın.

Görünüşe göre Nextflow, dizideki üç dizeyi ayrı değerler olarak kullanmak yerine, `[Hello, Bonjour, Holà]`'yı bir dize değeri olarak kullanarak tek bir süreç çağrısı çalıştırmaya çalıştı.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

Yani soruna neden olan "paketleme".
Nextflow'un diziyi açmasını ve bireysel dizeleri kanala yüklemesini nasıl sağlarız?

### 3.2. Kanal içeriklerini dönüştürmek için bir operatör kullanın

İşte [**operatörler**](https://nextflow.io/docs/latest/reference/operator.html) burada devreye giriyor.
Zaten `.view()` operatörünü kullandınız, bu sadece içinde ne olduğuna bakar.
Şimdi bir kanalın içerikleri üzerinde işlem yapmamıza olanak tanıyan operatörlere bakacağız.

Nextflow belgelerindeki [operatör listesine](https://www.nextflow.io/docs/latest/reference/operator.html) göz atarsanız, tam olarak ihtiyacımız olanı yapan [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten) operatörünü bulacaksınız: bir dizinin içeriğini açar ve bunları ayrı öğeler olarak yayınlar.

#### 3.2.1. `flatten()` operatörünü ekleyin

`flatten()` operatörünü girdi kanalımıza uygulamak için, kanal fabrikası tanımlamasına ekliyoruz.

İş akışı bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // bir girdi selamlama dizisi tanımla
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
        // bir girdi selamlama dizisi tanımla
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

Burada okunabilirlik için operatörü bir sonraki satıra ekledik, ancak tercih ederseniz operatörleri kanal fabrikasıyla aynı satıra ekleyebilirsiniz, şöyle:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. `view()` ifadesini/ifadelerini iyileştirin

Bunu hemen test etmek için çalıştırabiliriz, ancak bu sırada kanal içeriklerini nasıl incelediğimizi iyileştireceğiz.

`flatten()` operatörü uygulanmadan önce ve sonra içeriklerin nasıl göründüğünü karşılaştırabilmek istiyoruz, bu yüzden ikinci bir tane ekleyeceğiz VE çıktıda daha net etiketlenmelerini sağlamak için biraz kod ekleyeceğiz.

İş akışı bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // bir girdi selamlama dizisi tanımla
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "flatten öncesi: $greeting" }
                             .flatten()
                             .view { greeting -> "flatten sonrası: $greeting" }
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
        // bir girdi selamlama dizisi tanımla
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

İkinci bir `.view` ifadesi eklediğimizi ve her biri için boş parantezleri (`()`) biraz kod içeren süslü parantezlerle değiştirdiğimizi görüyorsunuz, örneğin `{ greeting -> "flatten öncesi: $greeting" }`.

Bunlara _closure_ denir. İçerdikleri kod, kanaldaki her öğe için yürütülecektir.
İç değer için burada `greeting` adında (ancak herhangi bir rastgele ad olabilir) geçici bir değişken tanımlıyoruz, bu yalnızca o closure kapsamında kullanılır.

Bu örnekte, `$greeting` kanalda yüklenen her bir öğeyi temsil eder.
Bu, düzgün etiketlenmiş konsol çıktısı ile sonuçlanacaktır.

!!! info "Bilgi"

    Bazı pipeline'larda operatör closure'ları içinde `$it` adlı özel bir değişkenin kullanıldığını görebilirsiniz.
    Bu, `->` ile tanımlamaya gerek kalmadan iç değişkene kısa yoldan erişim sağlayan bir _örtük_ değişkendir.

    Kod netliğine yardımcı olmak için açık olmayı tercih ediyoruz, bu nedenle `$it` sözdizimi önerilmemektedir ve Nextflow dilinden yavaş yavaş kaldırılacaktır.

#### 3.2.3. İş akışını çalıştırın

Son olarak, iş akışını tekrar çalıştırabilirsiniz!

```bash
nextflow run hello-channels.nf
```

??? success "Komut çıktısı"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    flatten öncesi: [Hello, Bonjour, Holà]
    flatten sonrası: Hello
    flatten sonrası: Bonjour
    flatten sonrası: Holà
    ```

Bu sefer çalışıyor VE `flatten()` operatörünü çalıştırmadan önce ve sonra kanal içeriklerinin nasıl göründüğüne dair ek bilgi veriyor.

- Tek bir `flatten öncesi:` ifadesi görüyorsunuz çünkü o noktada kanal bir öğe içeriyor, orijinal dizi.
- Sonra, artık kanalda bireysel öğeler olan her selamlama için bir tane olmak üzere üç ayrı `flatten sonrası:` ifadesi görüyorsunuz.

Önemli olarak, bu her öğenin artık iş akışı tarafından ayrı ayrı işlenebileceği anlamına gelir.

!!! tip "İpucu"

    Aynı sonuçları, işleminde örtük bir eşleme adımı içeren farklı bir kanal fabrikası [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist) kullanarak elde etmek teknik olarak mümkündür.
    Burada, basit bir kullanım durumunda bir operatörün kullanımını göstermek için bunu kullanmamayı tercih ettik.

### Özet

Bir kanalın içeriklerini dönüştürmek için `flatten()` gibi bir operatörü nasıl kullanacağınızı ve bir operatör uygulamadan önce ve sonra kanal içeriklerini incelemek için `view()` operatörünü nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

İş akışının girdi değerlerinin kaynağı olarak bir dosya almasını nasıl sağlayacağınızı öğrenin.

---

## 4. Bir CSV dosyasından girdi değerlerini okuyun

Gerçekçi olarak, nadiren bir değer dizisinden başlayacağız.
Büyük olasılıkla, işlenmesi gereken verileri içeren bir veya daha fazla dosyamız olacak, bir tür yapılandırılmış formatta.

`greetings.csv` adlı, gerçek bir veri analizinde işlemek isteyebileceğiniz türden sütunlu verileri taklit eden birkaç girdi selamlaması içeren bir CSV dosyası hazırladık; bu dosya `data/` altında depolanıyor.
(Sayılar anlamlı değil, sadece gösterim amaçlı oradalar.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Bir sonraki görevimiz, bu dosyadan değerleri okumak için iş akışımızı uyarlamak.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Bunu nasıl gerçekleştirebileceğimizi görelim.

### 4.1. Betiği selamlamaların kaynağı olarak bir CSV dosyası bekleyecek şekilde değiştirin

Başlamak için, betikte iki önemli değişiklik yapmamız gerekecek:

- Girdi parametresini CSV dosyasına işaret edecek şekilde değiştirin
- Kanal fabrikasını bir dosyayı işlemek için tasarlanmış birine değiştirin

#### 4.1.1. Girdi parametresini CSV dosyasına işaret edecek şekilde değiştirin

1. Bölümde kurduğumuz `params.input` parametresini hatırlıyor musunuz?
   Selamlamalarımızı içeren CSV dosyasına işaret edecek şekilde güncelleyeceğiz.

Parametre tanımlamasında aşağıdaki düzenlemeyi yapın:

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

Bu, dosyanın iş akışı koduyla aynı konumda olduğunu varsayar.
Diğer veri konumlarıyla nasıl başa çıkacağınızı Nextflow yolculuğunuzda daha sonra öğreneceksiniz.

#### 4.1.2. Bir dosyayı işlemek için tasarlanmış bir kanal fabrikasına geçin

Artık girdi olarak basit dizeler yerine bir dosya kullanmak istediğimiz için, önceki `channel.of()` kanal fabrikasını kullanamayız.
Dosya yollarını işlemek için yerleşik işlevselliğe sahip yeni bir kanal fabrikası olan [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) kullanmaya geçmemiz gerekiyor.

İş akışı bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "flatten öncesi: $greeting" }
                             // .flatten()
                             // .view { greeting -> "flatten sonrası: $greeting" }
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
        // bir girdi selamlama dizisi tanımla
        greetings_array = ['Hello','Bonjour','Holà']
        // girdiler için bir kanal oluştur
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "flatten öncesi: $greeting" }
                             .flatten()
                             .view { greeting -> "flatten sonrası: $greeting" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Kanal girdisini `param.input`'a geri döndürdüğümüzü ve artık ihtiyacımız olmayacağı için `greetings_array` tanımlamasını sildiğimizi fark edeceksiniz.
Ayrıca `flatten()` ve ikinci `view()` ifadesini yorum satırı haline getirdik.

#### 4.1.3. İş akışını çalıştırın

Yeni kanal fabrikası ve girdi dosyası ile iş akışını çalıştırmayı deneyelim.

```bash
nextflow run hello-channels.nf
```

??? failure "Komut çıktısı"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    flatten öncesi: /workspaces/training/hello-nextflow/data/greetings.csv
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

Hay aksi, çalışmıyor. Konsol çıktısının ve hata mesajının başlangıcına bakın.
`Command executed:` kısmı burada özellikle faydalı.

Bu biraz tanıdık görünebilir.
Görünüşe göre Nextflow, dosya yolunun kendisini bir dize değeri olarak kullanarak tek bir süreç çağrısı çalıştırmaya çalıştı.
Dosya yolunu doğru şekilde çözümledi, ancak istediğimiz şey olan içeriğini ayrıştırmadı.

Nextflow'un dosyayı açmasını ve içeriğini kanala yüklemesini nasıl sağlarız?

Görünüşe göre başka bir [operatöre](https://nextflow.io/docs/latest/reference/operator.html) ihtiyacımız var!

### 4.2. Dosyayı ayrıştırmak için `splitCsv()` operatörünü kullanın

Operatör listesine tekrar baktığımızda, CSV formatındaki metni ayrıştırmak ve bölmek için tasarlanmış [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) operatörünü buluyoruz.

#### 4.2.1. `splitCsv()`'yi kanala uygulayın

Operatörü uygulamak için, daha önce olduğu gibi kanal fabrikası satırına ekliyoruz.

İş akışı bloğunda, `flatten()` yerine `splitcsv()` (yorum satırından çıkarılmış) koymak için aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "splitCsv öncesi: $csv" }
                             .splitCsv()
                             .view { csv -> "splitCsv sonrası: $csv" }
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
                             .view { greeting -> "flatten öncesi: $greeting" }
                             // .flatten()
                             // .view { greeting -> "flatten sonrası: $greeting" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Gördüğünüz gibi, önce/sonra `view()` ifadelerini de güncelledik.
Teknik olarak aynı değişken adını (`greeting`) kullanabilirdik ama kodu başkaları için daha okunabilir kılmak için daha uygun bir şeye (`csv`) güncelledik.

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
    splitCsv öncesi: /workspaces/training/hello-nextflow/data/greetings.csv
    splitCsv sonrası: [Hello, English, 123]
    splitCsv sonrası: [Bonjour, French, 456]
    splitCsv sonrası: [Holà, Spanish, 789]
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
Bu sefer Nextflow dosyanın içeriğini ayrıştırdı (yaşasın!) ama her satırı bir dizi olarak yükledi ve her dizi kanalda bir öğe.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

Ona her satırdaki yalnızca ilk sütunu almasını söylememiz gerekiyor.
Peki bunu nasıl açarız?

Daha önce bir kanalın içeriğini açmak için `flatten()` kullandık, ama bu burada işe yaramaz çünkü flatten _her şeyi_ açar (kendiniz görmek isterseniz denemekten çekinmeyin).

Bunun yerine, Nextflow pipeline'larında gerçekten faydalı olan ve sık sık karşımıza çıkan `map()` adlı başka bir operatör kullanacağız.

### 4.3. Selamlamaları çıkarmak için `map()` operatörünü kullanın

[`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) operatörü, bir kanalın içeriklerine her türlü eşlemeyi yapmamıza olanak tanıyan çok kullanışlı küçük bir araçtır.

Bu durumda, veri dosyamızdaki her satırdan istediğimiz o tek öğeyi çıkarmak için kullanacağız.
Sözdizimi şöyle görünüyor:

```groovy title="Sözdizimi"
.map { row -> row[0] }
```

Bu, 'kanaldaki her satır için, içerdiği 0. (birinci) öğeyi al' anlamına gelir.

Şimdi bunu CSV ayrıştırmamıza uygulayalım.

#### 4.3.1. `map()`'i kanala uygulayın

İş akışı bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "splitCsv öncesi: $csv" }
                             .splitCsv()
                             .view { csv -> "splitCsv sonrası: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "map sonrası: $csv" }
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
                             .view { csv -> "splitCsv öncesi: $csv" }
                             .splitCsv()
                             .view { csv -> "splitCsv sonrası: $csv" }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Operatörün beklediğimizi yaptığını doğrulamak için başka bir `view()` çağrısı eklediğimizi görüyorsunuz.

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
    splitCsv öncesi: /workspaces/training/hello-nextflow/data/greetings.csv
    splitCsv sonrası: [Hello, English, 123]
    splitCsv sonrası: [Bonjour, French, 456]
    splitCsv sonrası: [Holà, Spanish, 789]
    map sonrası: Hello
    map sonrası: Bonjour
    map sonrası: Holà
    ```

Bu sefer hatasız çalışmalı.

`view()` ifadelerinin çıktısına bakarak şunları görüyorsunuz:

- Tek bir `splitCsv öncesi:` ifadesi: o noktada kanal bir öğe içeriyor, orijinal dosya yolu.
- Üç ayrı `splitCsv sonrası:` ifadesi: her selamlama için bir tane, ama her biri dosyadaki o satıra karşılık gelen bir dizi içinde bulunuyor.
- Üç ayrı `map sonrası:` ifadesi: her selamlama için bir tane, bunlar artık kanalda bireysel öğeler.

_Satırların çıktınızda farklı bir sırada görünebileceğini unutmayın._

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

Her selamlamanın doğru şekilde çıkarıldığını ve iş akışı boyunca işlendiğini doğrulamak için çıktı dosyalarına da bakabilirsiniz.

Daha önce olduğu gibi aynı sonuca ulaştık, ancak şimdi herhangi bir kodu değiştirmeden bir girdi dosyasını değiştirerek işlemek istediğimiz selamlama kanalına daha fazla öğe ekleme konusunda çok daha fazla esnekliğe sahibiz.
Karmaşık girdileri işlemek için daha sofistike yaklaşımları daha sonraki bir eğitimde öğreneceksiniz.

### Özet

Bir girdi değerleri dosyasını okumak ve bunları uygun şekilde işlemek için `.fromPath()` kanal yapıcısını ve `splitCsv()` ile `map()` operatörlerini nasıl kullanacağınızı biliyorsunuz.

Daha genel olarak, Nextflow'un süreçlere girdileri yönetmek için **kanalları** ve içeriklerini dönüştürmek için **operatörleri** nasıl kullandığına dair temel bir anlayışa sahipsiniz.
Ayrıca kanalların paralel yürütmeyi örtük olarak nasıl işlediğini de gördünüz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### Sırada ne var?

Büyük bir mola verin, bu bölümde çok çalıştınız!

Hazır olduğunuzda, daha fazla adım eklemeyi ve bunları düzgün bir iş akışına bağlamayı öğrenmek için [**Bölüm 3: Hello Workflow**](./03_hello_workflow.md)'na geçin.

---

## Quiz

<quiz>
Nextflow'da kanal nedir?
- [ ] Bir dosya yolu belirtimi
- [ ] Bir süreç tanımı
- [x] Süreçler arasında veri iletmek için kuyruk benzeri bir yapı
- [ ] Bir yapılandırma ayarı

Daha fazla bilgi: [1.1. Bir girdi kanalı oluşturun](#11-bir-girdi-kanali-olusturun)
</quiz>

<quiz>
Bu kod ne çıktı verecek?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (tek bir liste)
- [x] Her öğe ayrı bir satırda: `Hello`, `Bonjour`, `Hola`
- [ ] Hiçbir şey (kanallar varsayılan olarak yazdırmaz)
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
- [ ] Birden fazla kanalı bir araya getirir
- [ ] Kanal öğelerini sıralar
- [x] Dizileri bireysel öğelere açar
- [ ] Yinelenen öğeleri kaldırır

Daha fazla bilgi: [3.2.1. `flatten()` operatörünü ekleyin](#321-flatten-operatorunu-ekleyin)
</quiz>

<quiz>
`view()` operatörünün amacı nedir?
- [ ] Kanal içeriklerini filtrelemek
- [ ] Kanal öğelerini dönüştürmek
- [x] Kanal içeriklerini incelemek ve hata ayıklamak
- [ ] Kanal içeriklerini bir dosyaya kaydetmek

Daha fazla bilgi: [1.4. Kanal içeriklerini incelemek için `view()` kullanın](#14-kanal-iceriklerini-incelemek-icin-view-kullanin)
</quiz>

<quiz>
`splitCsv()` ne yapar?
- [ ] Kanal içeriklerinden bir CSV dosyası oluşturur
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
- [ ] Performansı iyileştirmek için
- [ ] Disk alanını azaltmak için
- [x] Çıktı dosyalarının birbirinin üzerine yazmasını önlemek için
- [ ] Devam etme işlevselliğini etkinleştirmek için

Daha fazla bilgi: [2.2. Çıktı dosya adlarının benzersiz olmasını sağlayın](#22-cikti-dosya-adlarinin-benzersiz-olmasini-saglayin)
</quiz>
