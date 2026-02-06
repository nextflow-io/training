# Bölüm 1: Temel işlemleri çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow Run eğitim kursunun bu ilk bölümünde, temel işlemleri göstermek ve ilgili Nextflow kod bileşenlerini işaret etmek için kullanacağımız çok temel, alandan bağımsız bir Hello World örneğiyle konuya yumuşak bir giriş yapıyoruz.

??? info "Hello World örneği nedir?"

    "Hello World!", bir programlama dilinin veya yazılım çerçevesinin temel sözdizimini ve yapısını göstermek için tasarlanmış minimalist bir örnektir.
    Örnek genellikle "Hello, World!" ifadesini konsol veya terminal gibi çıktı aygıtına yazdırmak veya bir dosyaya yazmaktan oluşur.

---

## 1. Doğrudan bir Hello World çalıştırın

Nextflow'a sarmadan önce ne yaptığını göstermek için doğrudan terminalde çalıştırdığımız basit bir komutla bu kavramı gösterelim.

!!! tip "İpucu"

    [Başlangıç](00_orientation.md) sayfasında açıklandığı gibi artık `nextflow-run/` dizininde olmanız gerektiğini unutmayın.

### 1.1. Terminalin merhaba demesini sağlayın

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
echo 'Hello World!'
```

??? success "Komut çıktısı"

    ```console
    Hello World!
    ```

Bu, 'Hello World' metnini doğrudan terminalde çıktı olarak verir.

### 1.2. Çıktıyı bir dosyaya yazın

Pipeline'ları çalıştırmak çoğunlukla dosyalardan veri okumayı ve sonuçları başka dosyalara yazmayı içerir, bu yüzden örneği biraz daha ilgili hale getirmek için metin çıktısını bir dosyaya yazacak şekilde komutu değiştirelim.

```bash
echo 'Hello World!' > output.txt
```

??? success "Komut çıktısı"

    ```console

    ```

Bu, terminale hiçbir şey çıktı vermez.

### 1.3. Çıktıyı bulun

'Hello World' metni artık belirttiğimiz çıktı dosyasında, `output.txt` adıyla olmalıdır.
Dosya gezgininde açabilir veya örneğin `cat` yardımcı programını kullanarak komut satırından açabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

İlk Nextflow workflow'umuzla çoğaltmaya çalışacağımız şey budur.

### Özet

Artık terminalde bazı metinler çıktı veren basit bir komutu nasıl çalıştıracağınızı ve isteğe bağlı olarak çıktıyı bir dosyaya nasıl yazdıracağınızı biliyorsunuz.

### Sırada ne var?

Aynı sonucu elde eden bir Nextflow workflow'u çalıştırmanın ne gerektirdiğini öğrenin.

---

## 2. Workflow'u çalıştırın

Size `--input` adlı bir komut satırı argümanı aracılığıyla bir girdi selamlama alan ve bu selamlamayı içeren bir metin dosyası üreten `1-hello.nf` adlı bir workflow betiği sağlıyoruz.

Henüz koda bakmayacağız; önce çalıştırmanın nasıl göründüğünü görelim.

### 2.1. Workflow'u başlatın ve çalışmayı izleyin

Terminalde aşağıdaki komutu çalıştırın:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Konsol çıktınız buna benzer görünüyorsa, tebrikler, ilk Nextflow workflow'unuzu çalıştırdınız!

Buradaki en önemli çıktı, yukarıdaki çıktıda vurgulanan son satırdır:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Bu bize `sayHello` sürecinin bir kez başarıyla çalıştırıldığını söyler (`1 of 1 ✔`).

Bu harika, ama merak ediyor olabilirsiniz: çıktı nerede?

### 2.2. Çıktı dosyasını `results` dizininde bulun

Bu workflow, çıktısını bir results dizinine yayınlamak üzere yapılandırılmıştır.
Geçerli dizininize bakarsanız, workflow'u çalıştırdığınızda Nextflow'un `results` adında yeni bir dizin ve bunun altında `output.txt` adlı bir dosya içeren `1-hello` adlı bir alt dizin oluşturduğunu göreceksiniz.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Dosyayı açın; içerik, komut satırında belirttiğiniz dizeyle eşleşmelidir.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Harika, workflow'umuz yapması gerekeni yaptı!

### 2.3. Sonuçları farklı bir dizine kaydedin

Varsayılan olarak, Nextflow pipeline çıktılarını geçerli yolunuzda `results` adlı bir dizine kaydeder.
Dosyalarınızın nereye yayınlanacağını değiştirmek için `-output-dir` CLI bayrağını (veya kısaca `-o`) kullanın.

!!! danger "Dikkat"

    `--input`'un iki tire, `-output-dir`'in bir tire aldığını unutmayın!
    Bunun nedeni `--input`'un bir pipeline _parametresi_ ve `-output-dir`'in çekirdek bir Nextflow CLI bayrağı olmasıdır.
    Bunlar hakkında daha sonra daha fazla bilgi.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Artık çıktılarınızın `results` yerine `hello_results` adlı bir dizine yayınlandığını görmelisiniz:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Bu dizindeki dosyalar öncekiyle aynıdır, sadece üst düzey dizin farklıdır.
Ancak, her iki durumda da 'yayınlanan' sonucun, Nextflow workflow'u çalıştırdığında ürettiği gerçek çıktının bir kopyası (veya bazı durumlarda sembolik bir bağlantı) olduğunu unutmayın.

Şimdi, Nextflow'un işi gerçekte nerede yaptığını görmek için kaputun altına bir göz atacağız.

!!! Warning "Uyarı"

    Tüm workflow'lar çıktıları bir results dizinine yayınlamak üzere ayarlanmayacaktır ve/veya dizin adları ve yapısı farklı olabilir.
    Bu bölümde biraz ileride, bu davranışın nerede belirtildiğini nasıl bulacağınızı göstereceğiz.

### 2.4. Orijinal çıktıyı ve logları `work/` dizininde bulun

Bir workflow çalıştırdığınızda, Nextflow workflow'daki her bir sürecin her çağrısı için (=pipeline'daki her adım) ayrı bir 'görev dizini' oluşturur.
Her biri için gerekli girdileri hazırlar, ilgili talimatları yürütür ve çıktıları ve log dosyalarını, benzersiz hale getirmek için otomatik olarak bir hash kullanılarak adlandırılan bu tek dizine yazar.

Bu görev dizinlerinin tümü, geçerli dizininizdeki (komutu çalıştırdığınız yer) `work` adlı bir dizin altında yaşayacaktır.

Bu karmaşık gelebilir, bu yüzden pratikte nasıl göründüğüne bakalım.

Daha önce çalıştırdığımız workflow için konsol çıktısına geri dönersek, şu satır vardı:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Satırın `[a3/1e1535]` ile nasıl başladığını görüyor musunuz?
Bu, o süreç çağrısı için görev dizini yolunun kısaltılmış halidir ve `work/` dizin yolu içinde `sayHello` süreç çağrısının çıktısını nerede bulacağınızı söyler.

Aşağıdaki komutu yazarak (kendi terminalinizde gördüğünüzle `a3/1e1535`'i değiştirerek) ve tab tuşuna basarak yolu otomatik tamamlayarak veya bir yıldız işareti ekleyerek tam yolu bulabilirsiniz:

```bash
ls work/a3/1e1535*
```

Bu, tam dizin yolunu vermeli: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

İçinde ne olduğuna bir göz atalım.

??? abstract "Dizin içeriği"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
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
    Ancak, log dosyaları terminalde görünmez olarak ayarlanmıştır, bu nedenle `ls` veya `tree` kullanarak görüntülemek istiyorsanız, görünmez dosyaları görüntülemek için ilgili seçeneği ayarlamanız gerekir.

    ```bash
    tree -a work
    ```

`work/`'de yaptığımız iki farklı pipeline çalıştırmasından iki dizin seti vardır.
Her görev yürütmesi, üzerinde çalışmak için kendi izole dizinini alır.
Bu durumda pipeline her iki seferde de aynı şeyi yaptı, bu nedenle her görev dizininin içeriği özdeştir.

`output.txt` dosyasını hemen tanımalısınız, bu aslında `results` dizinine yayınlanan `sayHello` sürecinin orijinal çıktısıdır.
Açarsanız, `Hello World!` selamlamasını tekrar bulacaksınız.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

Peki ya diğer tüm dosyalar?

Bunlar, Nextflow'un görev yürütmesinin bir parçası olarak yazdığı yardımcı ve log dosyalarıdır:

- **`.command.begin`**: Görev başlatıldığında oluşturulan sentinel dosyası.
- **`.command.err`**: Süreç çağrısı tarafından yayılan hata mesajları (`stderr`)
- **`.command.log`**: Süreç çağrısı tarafından yayılan tam log çıktısı
- **`.command.out`**: Süreç çağrısı tarafından normal çıktı (`stdout`)
- **`.command.run`**: Nextflow tarafından süreç çağrısını yürütmek için çalıştırılan tam betik
- **`.command.sh`**: Süreç çağrısı tarafından gerçekte çalıştırılan komut
- **`.exitcode`**: Komuttan kaynaklanan çıkış kodu

`.command.sh` dosyası özellikle yararlıdır çünkü size Nextflow'un yürüttüğü ana komutu gösterir, tüm defter tutma ve görev/ortam kurulumu dahil değildir.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Bu, workflow'un daha önce doğrudan komut satırında çalıştırdığımız aynı komutu oluşturduğunu doğrular.

Bir şeyler ters gittiğinde ve ne olduğunu gidermeniz gerektiğinde, Nextflow'un workflow talimatlarına, değişken enterpolasyonuna vb. göre tam olarak hangi komutu oluşturduğunu kontrol etmek için `command.sh` betiğine bakmak yararlı olabilir.

### 2.5. Workflow'u farklı selamlamalarla yeniden çalıştırın

Workflow'u `--input` argümanı için farklı değerlerle birkaç kez yeniden çalıştırmayı deneyin, ardından görev dizinlerine bakın.

??? abstract "Dizin içeriği"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Her çalıştırma için tam bir çıktı ve log dosyası seti içeren yeni bir alt dizin oluşturulduğunu görüyorsunuz.

Buna karşılık, `results` dizinine bakarsanız, hala yalnızca bir sonuç seti vardır ve çıktı dosyasının içeriği en son çalıştırdığınız şeye karşılık gelir.

??? abstract "Dizin içeriği"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Bu size yayınlanan sonuçların sonraki çalıştırmalar tarafından üzerine yazılacağını, oysa `work/` altındaki görev dizinlerinin korunduğunu gösterir.

### Özet

Basit bir Nextflow betiğini nasıl çalıştıracağınızı, yürütülmesini nasıl izleyeceğinizi ve çıktılarını nasıl bulacağınızı biliyorsunuz.

### Sırada ne var?

Temel bir Nextflow betiğini nasıl okuyacağınızı ve bileşenlerinin işlevselliğiyle nasıl ilişkili olduğunu öğrenin.

---

## 3. Hello World workflow başlangıç betiğini inceleyin

Orada yaptığımız şey temelde workflow betiğini bir kara kutu gibi ele almaktı.
Şimdi ne yaptığını gördüğümüze göre, kutuyu açalım ve içine bakalım.

Buradaki amacımız Nextflow kodunun sözdizimini ezberlemek değil, ana bileşenlerin ne olduğu ve nasıl organize edildikleri hakkında temel bir sezgi oluşturmaktır.

### 3.1. Genel kod yapısını inceleyin

`1-hello.nf` betiğini `nextflow-run` olması gereken geçerli dizininizde bulacaksınız. Editör bölmesinde açın.

??? full-code "Tam kod dosyası"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' ifadesini bir dosyaya yazdırmak için echo kullan
    */
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

    /*
    * Pipeline parametreleri
    */
    params {
        input: String
    }

    workflow {

        main:
        // bir selamlama yayınla
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

Bir Nextflow workflow betiği tipik olarak bir veya daha fazla **process** tanımı, **workflow**'un kendisi ve **params** ve **output** gibi birkaç isteğe bağlı blok içerir.

Her **process**, pipeline'daki ilgili adımın hangi işlem(ler)i gerçekleştirmesi gerektiğini açıklarken, **workflow** çeşitli adımları birbirine bağlayan veri akışı mantığını tanımlar.

Önce **process** bloğuna daha yakından bakalım, ardından **workflow** bloğuna bakacağız.

### 3.2. `process` tanımı

İlk kod bloğu bir [**process**](https://nextflow.io/docs/latest/process.html)'i tanımlar.
Process tanımı `process` anahtar kelimesiyle başlar, ardından process adı ve son olarak süslü parantezlerle sınırlandırılmış process gövdesi gelir.
Process gövdesi, çalıştırılacak komutu belirten bir script bloğu içermelidir; bu, komut satırı terminalinde çalıştırabileceğiniz herhangi bir şey olabilir.

```groovy title="1-hello.nf" linenums="3"
/*
* Bir selamlamayı bir dosyaya yazdırmak için echo kullan
*/
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

Burada `greeting` adlı bir **input** değişkeni alan ve **output**'unu `output.txt` adlı bir dosyaya yazan `sayHello` adlı bir **process**'imiz var.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Bu, yalnızca bir `input` tanımı, bir `output` tanımı ve yürütülecek `script`'i içeren çok minimal bir process tanımıdır.

`input` tanımı, Nextflow'a bir tür değer beklediğini söyleyen `val` niteleyicisini içerir (bir dize, bir sayı veya herhangi bir şey olabilir).

`output` tanımı, bunun bir yol olarak işlenmesi gerektiğini söyleyen `path` niteleyicisini içerir (hem dizin yollarını hem de dosyaları içerir).

### 3.3. `workflow` tanımı

İkinci kod bloğu [**workflow**](https://nextflow.io/docs/latest/workflow.html)'un kendisini tanımlar.
Workflow tanımı `workflow` anahtar kelimesiyle başlar, ardından isteğe bağlı bir ad ve süslü parantezlerle sınırlandırılmış workflow gövdesi gelir.

Burada `main:` bloğu ve `publish:` bloğundan oluşan bir **workflow**'umuz var.
`main:` bloğu workflow'un ana gövdesidir ve `publish:` bloğu `results` dizinine yayınlanması gereken çıktıları listeler.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // bir selamlama yayınla
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

Bu durumda `main:` bloğu `sayHello` sürecine bir çağrı içerir ve selamlama olarak kullanması için `params.input` adlı bir girdi verir.

Birazdan daha ayrıntılı tartışacağımız gibi, `params.input` komut satırımızda `--input` parametresine verdiğimiz değeri tutar.

`publish:` bloğu, `sayHello()` süreç çağrısının çıktısını listeler, buna `sayHello.out` olarak atıfta bulunur ve `first_output` adını verir (bu, workflow yazarının istediği herhangi bir şey olabilir).

Bu çok minimal bir **workflow** tanımıdır.
Gerçek dünya pipeline'larında, workflow tipik olarak **channel**'larla bağlanan birden fazla **process** çağrısı içerir ve değişken girdiler için varsayılan değerler ayarlanmış olabilir.

Kursun 2. Bölümünde buna gireceğiz.
Şimdilik, workflow'umuzun girdileri ve çıktıları nasıl işlediğine daha yakından bakalım.

### 3.4. Komut satırı parametreleri için `params` sistemi

`sayHello()` süreç çağrısına sağladığımız `params.input`, şık bir Nextflow kod parçasıdır ve üzerinde ekstra bir dakika harcamaya değer.

Yukarıda belirtildiği gibi, `--input` komut satırı parametresinin değerini `sayHello()` süreç çağrısına bu şekilde geçiriyoruz.
Aslında, sadece `params.someParameterName` bildirmek, workflow'a komut satırından `--someParameterName` adlı bir parametre vermek için yeterlidir.

Burada bu parametre bildirimini, workflow'un beklediği girdi türünü belirten bir `params` bloğu kurarak resmileştirdik (Nextflow 25.10.2 ve sonrası).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parametreleri
 */
params {
    input: String
}
```

Desteklenen türler arasında `String`, `Integer`, `Float`, `Boolean` ve `Path` bulunur.
Daha fazla bilgi için, Nextflow referans dokümantasyonundaki [Workflow parametreleri](https://nextflow.io/docs/latest/config.html#workflow-parameters) bölümüne bakın.

!!! tip "İpucu"

    `params` sistemi kullanılarak bildirilen _workflow_ parametreleri komut satırında her zaman iki tire alır (`--`).
    Bu, onları yalnızca bir tire alan (`-`) Nextflow düzeyindeki CLI bayraklarından ayırır.

### 3.5. `publish` direktifi

Workflow'un diğer ucunda, `publish:` bloğuna zaten göz attık.
Bu, çıktı işleme sisteminin yarısıdır; diğer yarısı aşağıda bulunan `output` bloğudur.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Bu, `publish:` bloğunda listelenen `first_output` çıktısının varsayılan `results` çıktı dizini altında `1-hello` adlı bir alt dizine kopyalanması gerektiğini belirtir.

`mode 'copy'` satırı, sistemin varsayılan davranışını geçersiz kılar, bu davranış düzgün bir kopya yerine `work/` dizinindeki orijinal dosyaya sembolik bir bağlantı (veya symlink) yapmaktır.

Yayınlama davranışını kontrol etmek için burada gösterilenden daha fazla seçenek vardır; daha sonra birkaçını ele alacağız.
Bir workflow birden fazla çıktı ürettiğinde, her birinin `output` bloğunda bu şekilde listelendiğini de göreceksiniz.

Daha fazla bilgi için, Nextflow referans dokümantasyonundaki [Çıktıları yayınlama](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) bölümüne bakın.

??? info "`publishDir` kullanarak çıktıları yayınlamak için eski sözdizimi"

    Çok yakın zamana kadar, çıktıları yayınlamanın yerleşik yolu, `publishDir` direktifi kullanarak her bir süreç düzeyinde yapmaktı.

    Bu kod kalıbını eski Nextflow pipeline'larında ve süreç modüllerinde hala her yerde bulacaksınız, bu yüzden bunun farkında olmak önemlidir.

    Workflow'da `publish:` bloğu ve üst düzeyde `output` bloğu yerine, `sayHello` süreç tanımında bir `publishDir` satırı görürdünüz:

    ```groovy title="Sözdizimi örneği" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    Ancak, gelecekteki Nextflow dil sürümlerinde sonunda yasaklanacağı için bunu herhangi bir yeni çalışmada kullanmanızı önermiyoruz.

### Özet

Artık basit bir Nextflow workflow'unun nasıl yapılandırıldığını ve temel bileşenlerin işlevselliğiyle nasıl ilişkili olduğunu biliyorsunuz.

### Sırada ne var?

Workflow çalıştırmalarınızı rahatça yönetmeyi öğrenin.

---

## 4. Workflow çalıştırmalarını yönetin

Workflow'ları başlatmayı ve çıktıları almayı bilmek harikadır, ancak hayatınızı kolaylaştıracak workflow yönetiminin birkaç başka yönü olduğunu hızla keşfedeceksiniz.

Burada size aynı workflow'u yeniden başlatmanız gerektiğinde `resume` özelliğinden nasıl yararlanacağınızı, çalıştırma loglarını `nextflow log` ile nasıl inceleyeceğinizi ve eski çalışma dizinlerini `nextflow clean` ile nasıl sileceğinizi gösteriyoruz.

### 4.1. `-resume` ile bir workflow'u yeniden başlatın

Bazen, daha önce başlattığınız bir pipeline'ı, önceden başarıyla tamamlanmış herhangi bir işi yeniden yapmadan çalıştırmak isteyeceksiniz.

Nextflow'un bunu yapmanıza olanak tanıyan `-resume` adlı bir seçeneği vardır.
Özellikle, bu modda, tam olarak aynı kod, ayarlar ve girdilerle zaten çalıştırılmış olan tüm süreçler atlanacaktır.
Bu, Nextflow'un yalnızca son çalıştırmadan bu yana eklediğiniz veya değiştirdiğiniz ya da yeni ayarlar veya girdiler sağladığınız süreçleri çalıştıracağı anlamına gelir.

Bunu yapmanın iki önemli avantajı vardır:

- Bir pipeline geliştirmenin ortasındaysanız, değişikliklerinizi test etmek için yalnızca aktif olarak üzerinde çalıştığınız süreç(ler)i çalıştırmanız gerektiğinden daha hızlı iterasyon yapabilirsiniz.
- Üretimde bir pipeline çalıştırıyorsanız ve bir şeyler ters giderse, birçok durumda sorunu düzeltebilir ve pipeline'ı yeniden başlatabilirsiniz ve başarısızlık noktasından itibaren çalışmaya devam eder, bu da size çok zaman ve hesaplama tasarrufu sağlayabilir.

Kullanmak için, komutunuza `-resume` ekleyin ve çalıştırın:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Komut çıktısı"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

Konsol çıktısı tanıdık görünmeli, ancak öncekine kıyasla biraz farklı olan bir şey var.

Süreç durum satırında (satır 5) eklenen `cached:` kısmına bakın, bu Nextflow'un bu işi zaten yaptığını tanıdığı ve önceki başarılı çalıştırmanın sonucunu yeniden kullandığı anlamına gelir.

Ayrıca çalışma alt dizini hash'inin önceki çalıştırmayla aynı olduğunu görebilirsiniz.
Nextflow kelimenin tam anlamıyla size önceki çalıştırmayı gösteriyor ve "Bunu zaten orada yaptım" diyor.

!!! tip "İpucu"

    Bir pipeline'ı `resume` ile yeniden çalıştırdığınızda, Nextflow daha önce başarıyla çalıştırılan çalıştırmalar tarafından çalışma dizininin dışında yayınlanan hiçbir dosyanın üzerine yazmaz.

    Daha fazla bilgi için, Nextflow referans dokümantasyonundaki [Cache ve resume](https://nextflow.io/docs/latest/cache-and-resume.html) bölümüne bakın.

### 4.2. Geçmiş çalıştırmaların logunu inceleyin

Bir nextflow workflow'u başlattığınızda, geçerli çalışma dizininde `.nextflow` adlı gizli bir dizin altında `history` adlı bir log dosyasına bir satır yazılır.

??? abstract "Dosya içeriği"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Bu dosya size geçerli çalışma dizininden başlatılan her Nextflow çalıştırması için zaman damgası, çalıştırma adı, durum, revizyon kimliği, oturum kimliği ve tam komut satırını verir.

Bu bilgilere erişmenin daha uygun bir yolu [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) komutunu kullanmaktır.

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

Bu, log dosyasının içeriğini bir başlık satırıyla zenginleştirerek terminale çıktı verecektir.

`-resume` seçeneğini kullanmıyorsanız, yeni bir `nextflow run` komutu çalıştırdığınızda oturum kimliğinin değiştiğini fark edeceksiniz.
Bu durumda, oturum kimliği aynı kalır.

Nextflow, çalıştırma önbellekleme bilgilerini gruplamak için oturum kimliğini, ayrıca `.nextflow` altında bulunan `cache` dizini altında kullanır.

### 4.3. Eski çalışma dizinlerini silin

Çok sayıda pipeline çalıştırırsanız, birçok alt dizin boyunca çok sayıda dosya biriktirebilirsiniz.
Alt dizinler rastgele adlandırıldığından, adlarından hangilerinin daha eski veya daha yeni çalıştırmalar olduğunu söylemek zordur.

Neyse ki Nextflow, artık umursamadığınız geçmiş çalıştırmalar için çalışma alt dizinlerini otomatik olarak silebilen yararlı bir [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) komutu içerir.

#### 4.3.1. Silme kriterlerini belirleyin

Neyin silineceğini belirlemek için birden fazla seçenek vardır, yukarıda bağlantısı verilen dokümantasyonda bunları keşfedebilirsiniz.
Burada size belirli bir çalıştırmadan önceki tüm alt dizinleri silen bir örnek gösteriyoruz, çalıştırma adı kullanılarak belirtilir.

`-resume` kullanmadığınız en son başarılı çalıştırmayı bulun; bizim durumumuzda çalıştırma adı `backstabbing_swartz` idi.

Çalıştırma adı, `Launching (...)` konsol çıktı satırında köşeli parantezler içinde gösterilen makine tarafından oluşturulan iki parçalı dizedir.
Ayrıca zaman damgası ve/veya komut satırına göre bir çalıştırmayı aramak için Nextflow logunu kullanabilirsiniz.

#### 4.3.2. Kuru çalıştırma yapın

Önce, komut verildiğinde neyin silineceğini kontrol etmek için kuru çalıştırma bayrağı `-n` kullanıyoruz:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Komut çıktısı"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Çıktınız farklı görev dizini adlarına sahip olacak ve farklı sayıda satır olabilir, ancak örneğe benzer görünmelidir.

Herhangi bir satır çıktı görmüyorsanız, ya geçerli bir çalıştırma adı sağlamadınız ya da silecek geçmiş çalıştırma yok. Örnek komuttaki `backstabbing_swartz`'ı loginizdeki ilgili en son çalıştırma adıyla değiştirdiğinizden emin olun.

#### 4.3.3. Silme işlemine devam edin

Çıktı beklendiği gibi görünüyorsa ve silme işlemine devam etmek istiyorsanız, `-n` yerine `-f` bayrağıyla komutu yeniden çalıştırın:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Komut çıktısı"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Çıktı daha önceye benzer olmalı, ancak şimdi 'Would remove' yerine 'Removed' diyor.
Bunun iki karakterli alt dizinleri (yukarıdaki `eb/` gibi) kaldırmadığını, ancak içeriklerini boşalttığını unutmayın.

!!! Warning "Uyarı"

    Geçmiş çalıştırmalardan çalışma alt dizinlerini silmek, onları Nextflow'un önbelleğinden kaldırır ve bu dizinlerde depolanan tüm çıktıları siler.
    Bu, Nextflow'un ilgili süreçleri yeniden çalıştırmadan çalışmaya devam etme yeteneğini bozar.

    Önemsediğiniz çıktıları kaydetmekten siz sorumlusunuz! Bu, `publish` direktifi için `symlink` modu yerine `copy` modunu tercih etmemizin ana nedenidir.

### Özet

Aynı şekilde zaten çalıştırılmış adımları tekrarlamadan bir pipeline'ı nasıl yeniden başlatacağınızı, çalıştırma logunu nasıl inceleyeceğinizi ve eski çalışma dizinlerini temizlemek için `nextflow clean` komutunu nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Biraz ara verin! Nextflow sözdiziminin ve temel kullanım talimatlarının yapı taşlarını yeni absorbe ettiniz.

Bu eğitimin bir sonraki bölümünde, Nextflow'un birden fazla girdiyi verimli bir şekilde işlemenize, birbirine bağlı birden fazla adımdan oluşan workflow'ları çalıştırmanıza, modüler kod bileşenlerinden yararlanmanıza ve daha fazla tekrar üretilebilirlik ve taşınabilirlik için konteynerları kullanmanıza nasıl izin verdiğini gösteren Hello World pipeline'ının art arda daha gerçekçi dört versiyonuna bakacağız.

---

## Quiz

<quiz>
`[a3/7be2fa] SAYHELLO | 1 of 1 ✔` konsol çıktı satırında `[a3/7be2fa]` neyi temsil eder?
- [ ] Süreç versiyon numarası
- [ ] Benzersiz bir çalıştırma tanımlayıcısı
- [x] Görevin çalışma dizinine kısaltılmış yol
- [ ] Çıktı dosyasının sağlama toplamı

Daha fazla bilgi: [2.4. Orijinal çıktıyı ve logları `work/` dizininde bulun](#24-orijinal-ciktiyi-ve-loglari-work-dizininde-bulun)
</quiz>

<quiz>
Görev dizinindeki `.command.sh` dosyasının amacı nedir?
- [ ] Görevin yapılandırma ayarlarını saklar
- [x] Süreç tarafından yürütülen gerçek komutu gösterir
- [ ] Başarısız görevlerden hata mesajlarını içerir
- [ ] Görev için hazırlanan girdi dosyalarını listeler

Daha fazla bilgi: [2.4. Orijinal çıktıyı ve logları `work/` dizininde bulun](#24-orijinal-ciktiyi-ve-loglari-work-dizininde-bulun)
</quiz>

<quiz>
`-resume` olmadan bir workflow'u yeniden çalıştırdığınızda yayınlanan sonuçlara ne olur?
- [ ] Ayrı zaman damgalı dizinlerde korunurlar
- [x] Yeni çalıştırma tarafından üzerine yazılırlar
- [ ] Nextflow üzerine yazmayı engeller ve başarısız olur
- [ ] Otomatik olarak yedeklenirler

Daha fazla bilgi: [2.5. Workflow'u farklı selamlamalarla yeniden çalıştırın](#25-workflowu-farkli-selamlamalarla-yeniden-calistirin)
</quiz>

<quiz>
Bu konsol çıktısı neyi gösterir?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Görev başarısız oldu ve atlandı
- [ ] Görev bir kuyrukta bekliyor
- [x] Nextflow önceki özdeş bir çalıştırmadan sonuçları yeniden kullandı
- [ ] Görev manuel olarak iptal edildi

Daha fazla bilgi: [4.1. `-resume` ile bir workflow'u yeniden başlatın](#41--resume-ile-bir-workflowu-yeniden-baslatin)
</quiz>

<quiz>
`nextflow log` komutunun görüntülediği çalıştırma geçmişini Nextflow nerede saklar?
- [ ] results dizininde
- [ ] work dizininde
- [x] `.nextflow/history` dosyasında
- [ ] `nextflow.config` dosyasında

Daha fazla bilgi: [4.2. Geçmiş çalıştırmaların logunu inceleyin](#42-gecmis-calistirmalarin-logunu-inceleyin)
</quiz>

<quiz>
Bir workflow dosyasındaki `params` bloğunun amacı nedir?
- [ ] Süreç kaynak gereksinimlerini tanımlamak
- [ ] Executor'ı yapılandırmak
- [x] Workflow girdi parametrelerini bildirmek ve türünü belirlemek
- [ ] Çıktı yayınlama seçeneklerini belirtmek

Daha fazla bilgi: [3.4. Komut satırı parametreleri için params sistemi](#34-komut-satiri-parametreleri-icin-params-sistemi)
</quiz>

<quiz>
Workflow'un `output` bloğunda `mode 'copy'` ne yapar?
- [ ] work dizininin bir yedeğini oluşturur
- [x] Sembolik bağlantılar yerine dosyaların tam kopyasını yapar
- [ ] Workflow betiğini sonuçlara kopyalar
- [ ] Artımlı dosya kopyalamayı etkinleştirir

Daha fazla bilgi: [3.5. publish direktifi](#35-publish-direktifi)
</quiz>

<quiz>
Dosyaları gerçekten silmeden önce `nextflow clean` komutuyla kullanılması önerilen bayrak nedir?
- [x] `-n` (kuru çalıştırma) neyin silineceğini önizlemek için
- [ ] `-v` (ayrıntılı) ayrıntılı çıktı görmek için
- [ ] `-a` (tümü) tüm dizinleri seçmek için
- [ ] `-q` (sessiz) uyarıları bastırmak için

Daha fazla bilgi: [4.3. Eski çalışma dizinlerini silin](#43-eski-calisma-dizinlerini-silin)
</quiz>
