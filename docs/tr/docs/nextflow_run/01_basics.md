# Bölüm 1: Temel işlemleri çalıştırma

Nextflow Run eğitim kursunun bu ilk bölümünde, konuya çok temel, alana özgü olmayan bir Hello World örneğiyle başlıyoruz. Bu örneği temel işlemleri göstermek ve ilgili Nextflow kod bileşenlerine işaret etmek için kullanacağız.

??? info "Hello World örneği nedir?"

    "Hello World!" minimalist bir örnektir ve bir programlama dilinin veya yazılım çerçevesinin temel sözdizimini ve yapısını göstermeyi amaçlar.
    Örnek tipik olarak "Hello, World!" ifadesini konsol veya terminal gibi bir çıktı cihazına yazdırmayı veya bir dosyaya yazmayı içerir.

---

## 1. Doğrudan bir Hello World çalıştırma

Bu konsepti, terminalde doğrudan çalıştırdığımız basit bir komutla gösterelim; Nextflow'a sarmadan önce ne yaptığını görmek için.

!!! tip

    [Başlarken](00_orientation.md) sayfasında açıklandığı gibi şu anda `nextflow-run/` dizininde olmanız gerektiğini unutmayın.

### 1.1. Terminalin merhaba demesini sağlama

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
echo 'Hello World!'
```

??? success "Komut çıktısı"

    ```console
    Hello World!
    ```

Bu, 'Hello World' metnini doğrudan terminalde çıktı olarak verir.

### 1.2. Çıktıyı bir dosyaya yazma

Boru hatlarını çalıştırmak çoğunlukla dosyalardan veri okumayı ve sonuçları diğer dosyalara yazmayı içerir, bu yüzden örneği biraz daha alakalı hale getirmek için komutu metin çıktısını bir dosyaya yazacak şekilde değiştirelim.

```bash
echo 'Hello World!' > output.txt
```

??? success "Komut çıktısı"

    ```console

    ```

Bu terminale herhangi bir şey çıktı vermez.

### 1.3. Çıktıyı bulma

'Hello World' metni şimdi belirttiğimiz `output.txt` adlı çıktı dosyasında olmalıdır.
Dosyayı dosya gezgininde veya komut satırından örneğin `cat` yardımcı programını kullanarak açabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

İlk Nextflow iş akışımızla çoğaltmaya çalışacağımız şey budur.

### Özet

Artık terminalde bazı metinleri çıktı olarak veren basit bir komutu nasıl çalıştıracağınızı ve isteğe bağlı olarak çıktıyı bir dosyaya nasıl yazdıracağınızı biliyorsunuz.

### Sırada ne var?

Aynı sonucu elde eden bir Nextflow iş akışını çalıştırmanın ne gerektirdiğini öğrenin.

---

## 2. İş akışını çalıştırma

Size `--input` adlı bir komut satırı argümanı aracılığıyla bir giriş selamlaması alan ve bu selamlamayı içeren bir metin dosyası üreten `1-hello.nf` adlı bir iş akışı betiği sağlıyoruz.

Henüz koda bakmayacağız; önce çalıştırmanın nasıl göründüğünü görelim.

### 2.1. İş akışını başlatma ve yürütmeyi izleme

Terminalde aşağıdaki komutu çalıştırın.

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

Konsol çıktınız buna benziyorsa, tebrikler, ilk Nextflow iş akışınızı çalıştırdınız!

Buradaki en önemli çıktı, yukarıdaki çıktıda vurgulanan son satırdır:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Bu bize `sayHello` sürecinin başarıyla bir kez yürütüldüğünü söyler (`1 of 1 ✔`).

Bu harika, ancak merak ediyor olabilirsiniz: çıktı nerede?

### 2.2. Çıktı dosyasını `results` dizininde bulma

Bu iş akışı, çıktısını bir sonuçlar dizinine yayınlayacak şekilde yapılandırılmıştır.
Mevcut dizininize bakarsanız, iş akışını çalıştırdığınızda Nextflow'un `results` adlı yeni bir dizin ve bunun altında `output.txt` adlı bir dosya içeren `1-hello` adlı bir alt dizin oluşturduğunu göreceksiniz.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Dosyayı açın; içerik komut satırında belirttiğiniz dizeyle eşleşmelidir.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Harika, iş akışımız yapması gerekeni yaptı!

### 2.3. Sonuçları farklı bir dizine kaydetme

Varsayılan olarak, Nextflow boru hattı çıktılarını mevcut yolunuzda `results` adlı bir dizine kaydeder.
Dosyalarınızın nereye yayınlanacağını değiştirmek için `-output-dir` CLI bayrağını (veya kısaca `-o`) kullanın

!!! danger

    `--input`'un iki tire ve `-output-dir`'in bir tire aldığına dikkat edin!
    Bunun nedeni `--input`'un bir boru hattı _parametresi_ ve `-output-dir`'in bir çekirdek Nextflow CLI bayrağı olmasıdır.
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

Çıktılarınızın artık `results` yerine `hello_results` adlı bir dizine yayınlandığını görmelisiniz:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Bu dizindeki dosyalar öncekiyle aynıdır, sadece üst düzey dizin farklıdır.
Ancak, her iki durumda da 'yayınlanan' sonucun, Nextflow'un iş akışını yürüttüğünde ürettiği gerçek çıktının bir kopyası (veya bazı durumlarda sembolik bir bağlantı) olduğunu unutmayın.

Şimdi, Nextflow'un işi gerçekte nerede yürüttüğünü görmek için kaputun altına bakacağız.

!!! Warning

    Tüm iş akışları çıktıları bir sonuçlar dizinine yayınlayacak şekilde ayarlanmayacaktır ve/veya dizin adları ve yapısı farklı olabilir.
    Bu bölümde biraz daha ileride, bu davranışın nerede belirtildiğini nasıl öğreneceğinizi göstereceğiz.

### 2.4. Orijinal çıktıyı ve günlükleri `work/` dizininde bulma

Bir iş akışı çalıştırdığınızda, Nextflow iş akışındaki her sürecin her çağrısı için ayrı bir 'görev dizini' oluşturur (=boru hattındaki her adım).
Her biri için gerekli girdileri hazırlar, ilgili talimat(lar)ı yürütür ve çıktıları ve günlük dosyalarını o dizin içinde yazar; bu dizin benzersiz olması için bir hash kullanılarak otomatik olarak adlandırılır.

Tüm bu görev dizinleri, mevcut dizininizde (komutu çalıştırdığınız yerde) `work` adlı bir dizin altında bulunur.

Bu kafa karıştırıcı gelebilir, o yüzden pratikte bunun nasıl göründüğünü görelim.

Daha önce çalıştırdığımız iş akışının konsol çıktısına geri dönersek, şu satıra sahiptik:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Satırın `[a3/1e1535]` ile nasıl başladığını görüyor musunuz?
Bu, o süreç çağrısı için görev dizini yolunun kısaltılmış bir biçimidir ve `sayHello` süreç çağrısının çıktısını `work/` dizin yolu içinde nerede bulacağınızı söyler.

Aşağıdaki komutu yazarak (kendi terminalinizde gördüğünüz ile `a3/1e1535`'i değiştirerek) ve yolu otomatik tamamlamak için tab tuşuna basarak veya bir yıldız işareti ekleyerek tam yolu bulabilirsiniz:

```bash
ls work/a3/1e1535*
```

Bu, tam dizin yolunu vermelidir: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

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

    Görev alt dizininin içeriğine VSCode dosya gezgininde göz atarsanız, tüm dosyaları hemen göreceksiniz.
    Ancak, günlük dosyaları terminalde görünmez olarak ayarlanmıştır, bu nedenle bunları görüntülemek için `ls` veya `tree` kullanmak istiyorsanız, görünmez dosyaları görüntülemek için ilgili seçeneği ayarlamanız gerekir.

    ```bash
    tree -a work
    ```

`work/` içinde yaptığımız iki farklı boru hattı çalıştırmasından iki dizin seti vardır.
Her görev yürütmesi, çalışmak için kendi izole dizinini alır.
Bu durumda boru hattı her iki seferde de aynı şeyi yaptı, bu nedenle her görev dizininin içeriği aynıdır

`results` dizinine yayınlanan `sayHello` sürecinin orijinal çıktısı olan `output.txt` dosyasını hemen tanımalısınız.
Açarsanız, `Hello World!` selamlamasını tekrar bulacaksınız.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

Peki ya diğer tüm dosyalar?

Bunlar Nextflow'un görev yürütmesinin bir parçası olarak yazdığı yardımcı ve günlük dosyalarıdır:

- **`.command.begin`**: Görev başlatılır başlatılmaz oluşturulan gözcü dosyası.
- **`.command.err`**: Süreç çağrısı tarafından yayınlanan hata mesajları (`stderr`)
- **`.command.log`**: Süreç çağrısı tarafından yayınlanan tam günlük çıktısı
- **`.command.out`**: Süreç çağrısı tarafından normal çıktı (`stdout`)
- **`.command.run`**: Süreç çağrısını yürütmek için Nextflow tarafından çalıştırılan tam betik
- **`.command.sh`**: Süreç çağrısı tarafından gerçekte çalıştırılan komut
- **`.exitcode`**: Komuttan kaynaklanan çıkış kodu

`.command.sh` dosyası özellikle yararlıdır çünkü Nextflow'un yürüttüğü ana komutu gösterir, tüm defter tutma ve görev/ortam kurulumunu içermez.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Bu, iş akışının daha önce komut satırında doğrudan çalıştırdığımız komutun aynısını oluşturduğunu doğrular.

Bir şeyler ters gittiğinde ve ne olduğunu araştırmanız gerektiğinde, Nextflow'un iş akışı talimatlarına, değişken enterpolasyonuna vb. dayanarak tam olarak hangi komutu oluşturduğunu kontrol etmek için `command.sh` betiğine bakmak yararlı olabilir.

### 2.5. İş akışını farklı selamlamalarla yeniden çalıştırma

`--input` argümanı için farklı değerlerle iş akışını birkaç kez yeniden çalıştırmayı deneyin, ardından görev dizinlerine bakın.

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

Her çalıştırma için tam bir çıktı ve günlük dosyası seti içeren yeni bir alt dizinin oluşturulduğunu görüyorsunuz.

Buna karşılık, `results` dizinine bakarsanız, hala yalnızca bir sonuç seti vardır ve çıktı dosyasının içeriği en son çalıştırdığınız şeye karşılık gelir.

??? abstract "Dizin içeriği"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Bu size yayınlanan sonuçların sonraki yürütmeler tarafından üzerine yazılacağını, oysa `work/` altındaki görev dizinlerinin korunduğunu gösterir.

### Özet

Basit bir Nextflow betiğini nasıl çalıştıracağınızı, yürütmesini nasıl izleyeceğinizi ve çıktılarını nasıl bulacağınızı biliyorsunuz.

### Sırada ne var?

Temel bir Nextflow betiğini nasıl okuyacağınızı ve bileşenlerinin işlevselliğiyle nasıl ilişkili olduğunu öğrenin.

---

## 3. Hello World iş akışı başlangıç betiğini inceleme

Orada yaptığımız şey temelde iş akışı betiğini kara kutu gibi ele almaktı.
Şimdi ne yaptığını gördüğümüze göre, kutuyu açalım ve içine bakalım.

Buradaki amacımız Nextflow kodunun sözdizimini ezberlemek değil, ana bileşenlerin neler olduğu ve nasıl organize edildikleri hakkında bazı temel sezgiler oluşturmaktır.

### 3.1. Genel kod yapısını inceleme

`1-hello.nf` betiğini `nextflow-run` olması gereken mevcut dizininizde bulacaksınız. Editör bölmesinde açın.

??? full-code "Tam kod dosyası"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' yazdırmak için echo kullan
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
    * Boru hattı parametreleri
    */
    params {
        input: String
    }

    workflow {

        main:
        // emit a greeting
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

Bir Nextflow iş akışı betiği tipik olarak bir veya daha fazla **process** tanımı, **workflow**'un kendisi ve **params** ve **output** gibi birkaç isteğe bağlı blok içerir.

Her **process**, boru hattındaki ilgili adımın hangi işlem(ler)i gerçekleştirmesi gerektiğini tanımlarken, **workflow** çeşitli adımları birbirine bağlayan veri akışı mantığını tanımlar.

Önce **process** bloğuna daha yakından bakalım, ardından **workflow** bloğuna bakacağız.

### 3.2. `process` tanımı

İlk kod bloğu bir [**process**](https://nextflow.io/docs/latest/process.html) tanımlar.
Süreç tanımı `process` anahtar kelimesiyle başlar, ardından süreç adı ve son olarak süslü parantezlerle sınırlandırılmış süreç gövdesi gelir.
Süreç gövdesi, çalıştırılacak komutu belirten bir betik bloğu içermelidir; bu, komut satırı terminalinde çalıştırabileceğiniz herhangi bir şey olabilir.

```groovy title="1-hello.nf" linenums="3"
/*
* Bir dosyaya selamlama yazdırmak için echo kullan
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

Bu, sadece bir `input` tanımı, bir `output` tanımı ve yürütülecek `script`'i içeren çok minimal bir süreç tanımıdır.

`input` tanımı, Nextflow'a bir tür değer beklemesini söyleyen `val` niteleyicisini içerir (bir dize, bir sayı, her neyse olabilir).

`output` tanımı, Nextflow'a bunun bir yol olarak ele alınması gerektiğini söyleyen `path` niteleyicisini içerir (hem dizin yollarını hem de dosyaları içerir).

### 3.3. `workflow` tanımı

İkinci kod bloğu [**workflow**](https://nextflow.io/docs/latest/workflow.html)'un kendisini tanımlar.
İş akışı tanımı `workflow` anahtar kelimesiyle başlar, ardından isteğe bağlı bir ad, ardından süslü parantezlerle sınırlandırılmış iş akışı gövdesi gelir.

Burada bir `main:` bloğu ve bir `publish:` bloğundan oluşan bir **workflow**'umuz var.
`main:` bloğu iş akışının ana gövdesidir ve `publish:` bloğu `results` dizinine yayınlanması gereken çıktıları listeler.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

Bu durumda `main:` bloğu, `sayHello` sürecine bir çağrı içerir ve ona selamlama olarak kullanması için `params.input` adlı bir girdi verir.

Birazdan daha ayrıntılı olarak tartışacağımız gibi, `params.input` komut satırımızda `--input` parametresine verdiğimiz değeri tutar.

`publish:` bloğu, `sayHello()` süreç çağrısının çıktısını listeler, buna `sayHello.out` olarak atıfta bulunur ve `first_output` adını verir (bu, iş akışı yazarının istediği herhangi bir şey olabilir).

Bu çok minimal bir **workflow** tanımıdır.
Gerçek dünya boru hattında, iş akışı tipik olarak **channel**'lar tarafından bağlanan **process**'lere birden fazla çağrı içerir ve değişken girdiler için varsayılan değerler ayarlanmış olabilir.

Kursun 2. Bölümünde buna gireceğiz.
Şimdilik, iş akışımızın girdileri ve çıktıları nasıl ele aldığına daha yakından bakalım.

### 3.4. Komut satırı parametrelerinin `params` sistemi

`sayHello()` süreç çağrısına sağladığımız `params.input`, düzgün bir Nextflow kodu parçasıdır ve üzerinde fazladan bir dakika harcamaya değer.

Yukarıda belirtildiği gibi, `--input` komut satırı parametresinin değerini `sayHello()` süreç çağrısına bu şekilde aktarıyoruz.
Aslında, sadece `params.someParameterName` bildirmek, iş akışına komut satırından `--someParameterName` adlı bir parametre vermek için yeterlidir.

Burada, iş akışının beklediği girdi türünü belirten bir `params` bloğu kurarak bu parametre bildirimini resmileştirdik (Nextflow 25.10.2 ve sonrası).

```groovy title="1-hello.nf" linenums="20"
/*
 * Boru hattı parametreleri
 */
params {
    input: String
}
```

Desteklenen türler arasında `String`, `Integer`, `Float`, `Boolean` ve `Path` bulunur.
Daha fazla bilgi için Nextflow referans belgelerindeki [İş akışı parametreleri](https://nextflow.io/docs/latest/config.html#workflow-parameters) bölümüne bakın.

!!! tip

    `params` sistemi kullanılarak bildirilen _iş akışı_ parametrelerinin komut satırında her zaman iki tire aldığını (`--`) unutmayın.
    Bu, onları yalnızca bir tire (`-`) alan _Nextflow düzeyindeki_ CLI bayraklarından ayırır.

### 3.5. `publish` yönergesi

İş akışının diğer ucunda, `publish:` bloğuna zaten göz atmıştık.
Bu, çıktı işleme sisteminin bir yarısıdır; diğer yarısı aşağıda bulunan `output` bloğudur.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Bu, `publish:` bloğunda listelenen `first_output` çıktısının varsayılan `results` çıktı dizini altında `1-hello` adlı bir alt dizine kopyalanması gerektiğini belirtir.

`mode 'copy'` satırı, sistemin varsayılan davranışını geçersiz kılar; bu davranış, düzgün bir kopya yerine `work/` dizinindeki orijinal dosyaya sembolik bir bağlantı (veya symlink) yapmaktır.

Yayınlama davranışını kontrol etmek için burada gösterilenden daha fazla seçenek vardır; birkaçını daha sonra ele alacağız.
Ayrıca bir iş akışı birden fazla çıktı ürettiğinde, her birinin `output` bloğunda bu şekilde listelendiğini göreceksiniz.

Daha fazla bilgi için Nextflow referans belgelerindeki [Çıktıları yayınlama](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) bölümüne bakın.

??? info "`publishDir` kullanarak çıktıları yayınlamak için eski sözdizimi"

    Çok yakın zamana kadar, çıktıları yayınlamanın yerleşik yolu, bunu her bir süreç düzeyinde bir `publishDir` yönergesi kullanarak yapmaktı.

    Bu kod desenini eski Nextflow boru hatlarında ve süreç modüllerinde her yerde bulacaksınız, bu nedenle bunun farkında olmak önemlidir.

    İş akışında bir `publish:` bloğu ve üst düzeyde bir `output` bloğu yerine, `sayHello` süreç tanımında bir `publishDir` satırı görürdünüz:

    ```groovy title="Syntax example" linenums="1" hl_lines="3"
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

    Ancak, Nextflow dilinin gelecekteki sürümlerinde sonunda izin verilmeyeceği için bunu yeni herhangi bir çalışmada kullanmanızı önermiyoruz.

### Özet

Artık basit bir Nextflow iş akışının nasıl yapılandırıldığını ve temel bileşenlerin işlevselliğiyle nasıl ilişkili olduğunu biliyorsunuz.

### Sırada ne var?

İş akışı yürütmelerinizi rahatça yönetmeyi öğrenin.

---

## 4. İş akışı yürütmelerini yönetme

İş akışlarını nasıl başlatacağınızı ve çıktıları nasıl alacağınızı bilmek harika, ancak hayatınızı kolaylaştıracak iş akışı yönetiminin birkaç başka yönü olduğunu hızla göreceksiniz.

Burada size aynı iş akışını yeniden başlatmanız gerektiğinde `resume` özelliğinden nasıl yararlanacağınızı, yürütme günlüklerini `nextflow log` ile nasıl inceleyeceğinizi ve eski çalışma dizinlerini `nextflow clean` ile nasıl sileceğinizi gösteriyoruz.

### 4.1. Bir iş akışını `-resume` ile yeniden başlatma

Bazen, daha önce başlattığınız bir boru hattını, zaten başarıyla tamamlanmış herhangi bir işi yeniden yapmadan yeniden çalıştırmak isteyeceksiniz.

Nextflow'un bunu yapmanıza izin veren `-resume` adlı bir seçeneği vardır.
Özellikle, bu modda, tam olarak aynı kod, ayarlar ve girdilerle zaten çalıştırılmış olan tüm süreçler atlanacaktır.
Bu, Nextflow'un yalnızca son çalıştırmadan bu yana eklediğiniz veya değiştirdiğiniz süreçleri veya yeni ayarlar veya girdiler sağladığınız süreçleri çalıştıracağı anlamına gelir.

Bunu yapmanın iki temel avantajı vardır:

- Bir boru hattı geliştiriyorsanız, değişikliklerinizi test etmek için yalnızca aktif olarak üzerinde çalıştığınız süreci çalıştırmanız gerektiğinden daha hızlı yineleme yapabilirsiniz.
- Üretimde bir boru hattı çalıştırıyorsanız ve bir şeyler ters giderse, çoğu durumda sorunu düzeltebilir ve boru hattını yeniden başlatabilirsiniz ve hata noktasından çalışmaya devam edecektir, bu da size çok fazla zaman ve hesaplama tasarrufu sağlayabilir.

Kullanmak için komutunuza `-resume` ekleyin ve çalıştırın:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Komut çıktısı"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

Konsol çıktısı tanıdık görünmeli, ancak öncekine göre biraz farklı olan bir şey var.

Süreç durum satırına (satır 5) eklenmiş olan `cached:` kısmına bakın; bu, Nextflow'un bu işi zaten yaptığını tanıdığı ve önceki başarılı çalıştırmadan sonucu yeniden kullandığı anlamına gelir.

Ayrıca çalışma alt dizini hash'inin önceki çalıştırmadakiyle aynı olduğunu görebilirsiniz.
Nextflow kelimenin tam anlamıyla size önceki yürütmeyi gösteriyor ve "Bunu zaten orada yaptım" diyor.

!!! tip

    Bir boru hattını `resume` ile yeniden çalıştırdığınızda, Nextflow daha önce başarıyla çalıştırılan yürütmeler tarafından çalışma dizini dışında yayınlanan hiçbir dosyanın üzerine yazmaz.

    Daha fazla bilgi için Nextflow referans belgelerindeki [Önbellek ve devam ettirme](https://nextflow.io/docs/latest/cache-and-resume.html) bölümüne bakın.

### 4.2. Geçmiş yürütmelerin günlüğünü inceleme

Bir nextflow iş akışı başlattığınızda, mevcut çalışma dizinindeki `.nextflow` adlı gizli bir dizin altında `history` adlı bir günlük dosyasına bir satır yazılır.

??? abstract "Dosya içeriği"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Bu dosya size mevcut çalışma dizini içinden başlatılan her Nextflow çalıştırması için zaman damgası, çalıştırma adı, durum, revizyon kimliği, oturum kimliği ve tam komut satırını verir.

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

Bu, günlük dosyasının içeriğini bir başlık satırıyla genişletilmiş olarak terminale çıktı verecektir.

Oturum kimliğinin, `-resume` seçeneğini kullanıyorsanız HARIÇ, yeni bir `nextflow run` komutu çalıştırdığınızda her değiştiğini fark edeceksiniz.
Bu durumda, oturum kimliği aynı kalır.

Nextflow, çalıştırma önbelleğe alma bilgilerini `.nextflow` altında bulunan `cache` dizini altında gruplamak için oturum kimliğini kullanır.

### 4.3. Eski çalışma dizinlerini silme

Çok sayıda boru hattı çalıştırırsanız, birçok alt dizinde çok sayıda dosya biriktirmeye başlayabilirsiniz.
Alt dizinler rastgele adlandırıldığından, adlarından eski ve daha yeni çalıştırmaları ayırt etmek zordur.

Neyse ki Nextflow, artık umursamadığınız geçmiş çalıştırmaların çalışma alt dizinlerini otomatik olarak silebilen [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) adlı yararlı bir komut içerir.

#### 4.3.1. Silme kriterlerini belirleme

Neyin silineceğini belirlemek için yukarıda bağlantısı verilen belgelerde keşfedebileceğiniz birden fazla seçenek vardır.
Burada size, çalıştırma adı kullanılarak belirtilen belirli bir çalıştırmadan önceki çalıştırmaların tüm alt dizinlerini silen bir örnek gösteriyoruz.

`-resume` kullanmadığınız en son başarılı çalıştırmayı arayın; bizim durumumuzda çalıştırma adı `backstabbing_swartz` idi.

Çalıştırma adı, `Launching (...)` konsol çıktı satırında köşeli parantez içinde gösterilen makine tarafından oluşturulan iki parçalı dizedir.
Ayrıca zaman damgasına ve/veya komut satırına göre bir çalıştırmayı aramak için Nextflow günlüğünü kullanabilirsiniz.

#### 4.3.2. Kuru çalıştırma yapma

Önce komut verildiğinde neyin silineceğini kontrol etmek için kuru çalıştırma bayrağı `-n` kullanırız:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Komut çıktısı"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Çıktınız farklı görev dizini adlarına sahip olacak ve farklı sayıda satıra sahip olabilir, ancak örneğe benzer görünmelidir.

Herhangi bir satır çıktısı görmüyorsanız, geçerli bir çalıştırma adı sağlamadınız veya silinecek geçmiş çalıştırma yok. Örnek komuttaki `backstabbing_swartz`'ı günlüğünüzdeki karşılık gelen en son çalıştırma adıyla değiştirdiğinizden emin olun.

#### 4.3.3. Silme işlemine devam etme

Çıktı beklendiği gibi görünüyorsa ve silme işlemine devam etmek istiyorsanız, komutu `-n` yerine `-f` bayrağıyla yeniden çalıştırın:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Komut çıktısı"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Çıktı öncekine benzer olmalıdır, ancak şimdi 'Would remove' yerine 'Removed' diyor.
Bunun iki karakterli alt dizinleri (yukarıdaki `eb/` gibi) kaldırmadığını ancak içeriklerini boşalttığını unutmayın.

!!! Warning

    Geçmiş çalıştırmalardan çalışma alt dizinlerini silmek, onları Nextflow'un önbelleğinden kaldırır ve bu dizinlerde depolanan tüm çıktıları siler.
    Bu, Nextflow'un ilgili süreçleri yeniden çalıştırmadan yürütmeye devam etme yeteneğini bozar.

    Önem verdiğiniz çıktıları kaydetmek sizin sorumluluğunuzdadır! `publish` yönergesi için `symlink` modu yerine `copy` modunu kullanmayı tercih etmemizin ana nedeni budur.

### Özet

Zaten aynı şekilde çalıştırılmış adımları tekrarlamadan bir boru hattını nasıl yeniden başlatacağınızı, yürütme günlüğünü nasıl inceleyeceğinizi ve eski çalışma dizinlerini temizlemek için `nextflow clean` komutunu nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Kısa bir mola verin! Nextflow sözdiziminin yapı taşlarını ve temel kullanım talimatlarını yeni emdiniz.

Bu eğitimin bir sonraki bölümünde, Nextflow'un birden fazla girdiyi verimli bir şekilde işlemenize, birbirine bağlı birden fazla adımdan oluşan iş akışlarını çalıştırmanıza, modüler kod bileşenlerinden yararlanmanıza ve daha fazla tekrarlanabilirlik ve taşınabilirlik için konteynerlerden yararlanmanıza nasıl izin verdiğini gösterecek olan Hello World boru hattının art arda daha gerçekçi dört versiyonuna bakacağız.

---

## Quiz

<quiz>
`[a3/7be2fa] SAYHELLO | 1 of 1 ✔` konsol çıktı satırında, `[a3/7be2fa]` neyi temsil eder?
- [ ] Süreç sürüm numarası
- [ ] Benzersiz bir çalıştırma tanımlayıcısı
- [x] Görevin çalışma dizinine kısaltılmış yol
- [ ] Çıktı dosyasının sağlama toplamı

Daha fazla bilgi: [2.3. Orijinal çıktıyı ve günlükleri `work/` dizininde bulma](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
Bir görev dizinindeki `.command.sh` dosyasının amacı nedir?
- [ ] Görevin yapılandırma ayarlarını saklar
- [x] Süreç tarafından yürütülen gerçek komutu gösterir
- [ ] Başarısız görevlerden hata mesajları içerir
- [ ] Görev için hazırlanan girdi dosyalarını listeler

Daha fazla bilgi: [2.3. Orijinal çıktıyı ve günlükleri `work/` dizininde bulma](#23-find-the-original-output-and-logs-in-the-work-directory)
</quiz>

<quiz>
`-resume` olmadan bir iş akışını yeniden çalıştırdığınızda yayınlanan sonuçlara ne olur?
- [ ] Ayrı zaman damgalı dizinlerde korunurlar
- [x] Yeni yürütme tarafından üzerine yazılırlar
- [ ] Nextflow üzerine yazmayı önler ve başarısız olur
- [ ] Otomatik olarak yedeklenir

Daha fazla bilgi: [2.4. İş akışını farklı selamlamalarla yeniden çalıştırma](#24-re-run-the-workflow-with-different-greetings)
</quiz>

<quiz>
Bu konsol çıktısı neyi gösterir?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] Görev başarısız oldu ve atlandı
- [ ] Görev bir kuyrukta bekliyor
- [x] Nextflow önceki aynı yürütmeden sonuçları yeniden kullandı
- [ ] Görev manuel olarak iptal edildi

Daha fazla bilgi: [4.1. Bir iş akışını `-resume` ile yeniden başlatma](#41-re-launch-a-workflow-with--resume)
</quiz>

<quiz>
Nextflow, `nextflow log` komutunun görüntülediği yürütme geçmişini nerede saklar?
- [ ] Sonuçlar dizininde
- [ ] Çalışma dizininde
- [x] `.nextflow/history` dosyasında
- [ ] `nextflow.config` dosyasında

Daha fazla bilgi: [4.2. Geçmiş yürütmelerin günlüğünü inceleme](#42-inspect-the-log-of-past-executions)
</quiz>

<quiz>
Bir iş akışı dosyasındaki `params` bloğunun amacı nedir?
- [ ] Süreç kaynak gereksinimlerini tanımlamak
- [ ] Yürütücüyü yapılandırmak
- [x] İş akışı girdi parametrelerini bildirmek ve türlerini belirtmek
- [ ] Çıktı yayınlama seçeneklerini belirtmek

Daha fazla bilgi: [3.4. Komut satırı parametrelerinin `params` sistemi](#34-the-params-system-of-command-line-parameters)
</quiz>

<quiz>
İş akışının `output` bloğunda `mode 'copy'` ne yapar?
- [ ] Çalışma dizininin yedeğini oluşturur
- [x] Sembolik bağlantılar yerine dosyaların tam kopyasını yapar
- [ ] İş akışı betiğini sonuçlara kopyalar
- [ ] Artımlı dosya kopyalamayı etkinleştirir

Daha fazla bilgi: [3.5. `publish` yönergesi](#35-the-publish-directive)
</quiz>

<quiz>
Dosyaları gerçekten silmeden önce `nextflow clean` komutuyla kullanılması önerilen bayrak nedir?
- [x] `-n` (kuru çalıştırma) neyin silineceğini önizlemek için
- [ ] `-v` (ayrıntılı) detaylı çıktı görmek için
- [ ] `-a` (tümü) tüm dizinleri seçmek için
- [ ] `-q` (sessiz) uyarıları bastırmak için

Daha fazla bilgi: [4.3. Eski çalışma dizinlerini silme](#43-delete-older-work-directories)
</quiz>
