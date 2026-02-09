# Bölüm 2: Gerçek boru hatlarını çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kursun Bölüm 1'inde (Temel İşlemleri Çalıştırma), kod karmaşıklığını düşük tutmak için yalnızca minimal özelliklere sahip örnek bir iş akışı ile başladık.
Örneğin, `1-hello.nf` bir seferde tek bir değer sağlamak için bir komut satırı parametresi (`--input`) kullandı.

Ancak, gerçek dünyadaki boru hatlarının çoğu, büyük miktarda veriyi ölçekte verimli bir şekilde işlemeyi sağlamak ve bazen karmaşık mantıkla birbirine zincirlenmiş birden fazla işleme adımı uygulamak için daha gelişmiş özellikler kullanır.

Eğitimin bu bölümünde, orijinal Hello World boru hattının genişletilmiş sürümlerini deneyerek gerçek dünyadaki boru hatlarının temel özelliklerini gösteriyoruz.

## 1. Bir dosyadan girdi verilerini işleme

Gerçek dünyadaki bir boru hattında, genellikle bir veya daha fazla girdi dosyasında bulunan birden fazla veri noktasını (veya veri serisini) işlemek isteriz.
Ve mümkün olduğunda, analiz için bekleme süresini kısaltmak için bağımsız verilerin işlenmesini paralel olarak çalıştırmak isteriz.

Nextflow'un bunu nasıl yaptığını göstermek için, gerçek bir veri analizinde işlemek isteyebileceğiniz sütunlu veri türünü taklit eden birkaç girdi selamlaması içeren `greetings.csv` adlı bir CSV dosyası hazırladık.
Sayıların anlamlı olmadığını, sadece açıklayıcı amaçlarla orada olduklarını unutmayın.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Ayrıca, CSV dosyasını okuyacak, selamlamaları çıkaracak ve her birini ayrı bir dosyaya yazacak olan, şimdi `2a-inputs.nf` adlı orijinal iş akışının geliştirilmiş bir sürümünü yazdık.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Önce iş akışını çalıştıralım ve daha sonra ilgili Nextflow koduna bir göz atalım.

### 1.1. İş akışını çalıştırma

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Heyecan verici bir şekilde, bu '3 of 3' çağrısının süreç için yapıldığını gösteriyor gibi görünüyor, bu da umut verici çünkü girdi olarak sağladığımız CSV'de üç satır veri vardı.
Bu, `sayHello()` sürecinin her girdi satırında bir kez olmak üzere üç kez çağrıldığını gösteriyor.

### 1.2. Yayınlanan çıktıları `results` dizininde bulma

İş akışımızın hala çıktılarımızın bir kopyasını oraya yazıp yazmadığını görmek için 'results' dizinine bakalım.

??? abstract "Dizin içeriği"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Evet! Yeterince uygun bir şekilde, farklı adlara sahip üç çıktı dosyası içeren `2a-inputs` adlı yeni bir dizin görüyoruz.

Her birinin uygun selamlama dizesini içerdiğinden emin olmak için her birini açabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Bu, girdi dosyasındaki her selamlamanın uygun şekilde işlendiğini doğrular.

### 1.3. Orijinal çıktıları ve günlükleri bulma

Yukarıdaki konsol çıktısının yalnızca bir görev dizinine atıfta bulunduğunu fark etmiş olabilirsiniz.
Bu, `sayHello()` için yapılan üç çağrının da o tek görev dizininde yürütüldüğü anlamına mı geliyor?

#### 1.3.1. Terminalde verilen görev dizinini inceleme

O `8e/0eb066` görev dizininin içine bir göz atalım.

??? abstract "Dizin içeriği"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Yalnızca selamlamalardan birine karşılık gelen çıktıyı buluyoruz (gizli dosyaların görüntülenmesini etkinleştirirsek yardımcı dosyalarla birlikte).

Peki burada neler oluyor?

Varsayılan olarak, ANSI günlük sistemi aynı süreç için tüm çağrıların durum bilgilerini aynı satıra yazar.
Sonuç olarak, konsol çıktısında bize üç görev dizini yolundan (`8e/0eb066`) yalnızca birini gösterdi.
Orada listelenmeyen iki tane daha var.

#### 1.3.2. Terminalin daha fazla ayrıntı göstermesini sağlama

Komuta `-ansi-log false` ekleyerek günlük davranışını değiştirerek süreç çağrılarının tam listesini görebiliriz:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Komut çıktısı"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Bu sefer üç süreç çalıştırmasının tümünü ve ilişkili çalışma alt dizinlerini çıktıda listelenmiş olarak görüyoruz.
ANSI günlüğünü devre dışı bırakmak ayrıca Nextflow'un terminal çıktısında renk kullanmasını da engelledi.

Durumun iki günlük modu arasında biraz farklı şekilde raporlandığına dikkat edin.
Yoğunlaştırılmış modda, Nextflow çağrıların başarıyla tamamlanıp tamamlanmadığını bildirir.
Bu genişletilmiş modda, yalnızca gönderildiklerini bildirir.

Bu, `sayHello()` sürecinin üç kez çağrıldığını ve her biri için ayrı bir görev dizini oluşturulduğunu doğrular.

Orada listelenen görev dizinlerinin her birinin içine bakarsak, her birinin selamlamalardan birine karşılık geldiğini doğrulayabiliriz.

??? abstract "Dizin içeriği"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Bu, her süreç çağrısının diğerlerinden izole bir şekilde yürütüldüğünü doğrular.
Bunun, süreç benzersiz olmayan adlara sahip ara dosyalar üretirse çakışmaları önlemek de dahil olmak üzere birçok avantajı vardır.

!!! tip

    Karmaşık bir iş akışı veya çok sayıda girdi için, terminale çıktı olarak verilen tam listenin biraz bunaltıcı olabileceğini unutmayın, bu nedenle insanlar normalde rutin kullanımda `-ansi-log false` kullanmazlar.

### 1.4. İş akışı kodunu inceleme

Bu iş akışı sürümü, bir CSV girdi dosyasını okuyabilir, girdileri ayrı ayrı işleyebilir ve çıktıları benzersiz şekilde adlandırabilir.

İş akışı kodunda bunu mümkün kılan şeye bir göz atalım.

??? full-code "Tam kod dosyası"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' yazdırmak için echo kullan
    */
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

    /*
    * Boru hattı parametreleri
    */
    params {
        input: Path
    }

    workflow {

        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Bir kez daha, kod sözdizimini ezberlemek zorunda değilsiniz, ancak önemli işlevsellik sağlayan iş akışının temel bileşenlerini tanımayı öğrenmek iyidir.

#### 1.4.1. CSV'den girdi verilerini yükleme

Bu en ilginç kısım: komut satırından tek bir değer almaktan, bir CSV dosyası almaya, onu ayrıştırmaya ve içerdiği bireysel selamlamaları işlemeye nasıl geçtik?

Nextflow'da bunu bir [**kanal**](https://nextflow.io/docs/latest/channel.html) ile yapıyoruz: girdileri verimli bir şekilde işlemek ve çok adımlı iş akışlarında bir adımdan diğerine aktarmak için tasarlanmış, yerleşik paralellik ve birçok ek fayda sağlayan bir kuyruk yapısı.

Hadi parçalayalım.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // bir selamlama yayınla
    sayHello(greeting_ch)
```

Bu kod, CSV dosyasını okuyan, ayrıştıran ve her satırdan ilk sütunu çıkaran `greeting_ch` adlı bir kanal oluşturur.
Sonuç, `Hello`, `Bonjour` ve `Holà` içeren bir kanaldır.

??? tip "Bu nasıl çalışır?"

    İşte bu satırın sade Türkçe anlamı:

    - `channel.fromPath` dosya yolu(ları)ndan bir kanal oluşturan bir **kanal fabrikası**dır
    - `(params.input)` dosya yolunun komut satırında `--input` ile sağlandığını belirtir

    Başka bir deyişle, bu satır Nextflow'a şunu söyler: `--input` ile verilen dosya yolunu al ve içeriğini girdi verisi olarak işlemeye hazırlan.

    Ardından sonraki iki satır, dosyanın gerçek ayrıştırmasını ve verilerin uygun veri yapısına yüklenmesini yapan **operatörler** uygular:

    - `.splitCsv()` Nextflow'a CSV dosyasını satırları ve sütunları temsil eden bir diziye ayrıştırmasını söyler
    - `.map { line -> line[0] }` Nextflow'a her satırdan yalnızca ilk sütundaki öğeyi almasını söyler

    Yani pratikte, aşağıdaki CSV dosyasından başlayarak:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Bunu şuna benzeyen bir diziye dönüştürdük:

    ```txt title="Dizi içeriği"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    Ve sonra üç satırın her birinden ilk öğeyi aldık ve bunları şimdi şunları içeren bir Nextflow kanalına yükledik: `Hello`, `Bonjour` ve `Holà`.

    Kanalları ve operatörleri derinlemesine anlamak istiyorsanız, bunları kendiniz nasıl yazacağınız da dahil olmak üzere, [Hello Nextflow Bölüm 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) bölümüne bakın.

#### 1.4.2. Her selamlama için süreci çağırma

Ardından, iş akışının `main:` bloğunun son satırında, yüklenen `greeting_ch` kanalını `sayHello()` sürecine girdi olarak sağlıyoruz.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // bir selamlama yayınla
    sayHello(greeting_ch)
```

Bu, Nextflow'a süreci kanaldaki her öğe üzerinde, _yani_ her selamlama üzerinde ayrı ayrı çalıştırmasını söyler.
Ve Nextflow akıllı olduğu için, mevcut bilgi işlem altyapısına bağlı olarak mümkünse bu süreç çağrılarını paralel olarak çalıştıracaktır.

Nispeten çok az kodla çok miktarda verinin (birçok örnek veya veri noktası, araştırma biriminiz ne olursa olsun) verimli ve ölçeklenebilir işlenmesini bu şekilde başarabilirsiniz.

#### 1.4.3. Çıktıların nasıl adlandırıldığı

Son olarak, çıktı dosyalarının benzersiz şekilde adlandırılmasını nasıl sağladığımızı görmek için süreç koduna hızlıca bir göz atmaya değer.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

`1-hello.nf`'deki bu sürecin sürümüyle karşılaştırıldığında, çıktı bildiriminin ve komutun ilgili kısmının çıktı dosyası adına selamlama değerini içerecek şekilde değiştiğini görüyorsunuz.

Bu, çıktı dosyası adlarının ortak sonuçlar dizinine yayınlandıklarında çakışmayacağından emin olmanın bir yoludur.

Ve süreç bildirimi içinde yapmamız gereken tek değişiklik bu!

### Özet

Kanalların ve operatörlerin birden fazla girdiyi verimli bir şekilde işlememizi nasıl sağladığını temel düzeyde anlıyorsunuz.

### Sırada ne var?

Çok adımlı iş akışlarının nasıl oluşturulduğunu ve nasıl çalıştıklarını keşfedin.

---

## 2. Çok adımlı iş akışlarını çalıştırma

Gerçek dünyadaki iş akışlarının çoğu birden fazla adım içerir.
Kanallar hakkında öğrendiklerimizin üzerine inşa edelim ve Nextflow'un çok adımlı bir iş akışında süreçleri birbirine bağlamak için kanalları ve operatörleri nasıl kullandığına bakalım.

Bu amaçla, size üç ayrı adımı birbirine zincirleyen ve aşağıdakileri gösteren örnek bir iş akışı sunuyoruz:

1. Verilerin bir süreçten diğerine akmasını sağlama
2. Birden fazla süreç çağrısından çıktıları tek bir süreç çağrısında toplama

Özellikle, her girdi selamlamasını alan, büyük harfe dönüştüren ve ardından tüm büyük harfli selamlamaları tek bir çıktı dosyasında toplayan `2b-multistep.nf` adlı iş akışının genişletilmiş bir sürümünü yaptık.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Daha önce olduğu gibi, önce iş akışını çalıştıracağız, sonra yeni olan şeyi görmek için koda bakacağız.

### 2.1. İş akışını çalıştırma

Terminalinizde aşağıdaki komutu çalıştırın:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Komut çıktısı"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Söz verildiği gibi, iş akışının bir parçası olarak birden fazla adımın çalıştırıldığını görüyorsunuz; ilk ikisi (`sayHello` ve `convertToUpper`) muhtemelen her bir selamlama üzerinde çalıştırıldı ve üçüncüsü (`collectGreetings`) üç `convertToUpper` çağrısının çıktıları üzerinde yalnızca bir kez çalıştırılmış olacak.

### 2.2. Çıktıları bulma

Bunun gerçekten olup olmadığını doğrulamak için `results` dizinine bir göz atalım.

??? abstract "Dizin içeriği"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Gördüğünüz gibi, `2b-multistep` adlı yeni bir dizinimiz var ve öncekinden çok daha fazla dosya içeriyor.
Dosyalardan bazıları `intermediates` adlı bir alt dizinde gruplandırılmış, iki dosya ise üst düzeyde bulunuyor.

Bu ikisi, çok adımlı iş akışının nihai sonuçlarıdır.
Dosya adlarına bakmak ve beklediğiniz şey olduklarını doğrulamak için içeriklerini kontrol etmek için bir dakikanızı ayırın.

??? abstract "Dosya içeriği"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

İlki, söz verildiği gibi büyük harfe dönüştürülmüş ve tek bir dosyada toplanmış üç selamlamayı içerir.
İkincisi, çalıştırma hakkında bazı bilgileri özetleyen bir rapor dosyasıdır.

### 2.3. Kodu inceleme

Koda bakalım ve çok adımlı iş akışları için temel kalıpları belirleyelim.

??? full-code "Tam kod dosyası"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' yazdırmak için echo kullan
    */
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

    /*
    * Selamlamayı büyük harfe dönüştürmek için bir metin değiştirme aracı kullan
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Büyük harfli selamlamaları tek bir çıktı dosyasında topla
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Boru hattı parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Orada çok şey oluyor, ancak iş akışının önceki sürümüne kıyasla en belirgin fark, şimdi birden fazla süreç tanımı olması ve buna bağlı olarak iş akışı bloğunda birkaç süreç çağrısı olmasıdır.

Daha yakından bakalım ve en ilginç parçaları belirleyip belirleyemeyeceğimize bakalım.

#### 2.3.1. İş akışı yapısını görselleştirme

Nextflow uzantısı ile VSCode kullanıyorsanız, herhangi bir Nextflow betiğindeki iş akışı bloğunun hemen üstünde görüntülenen küçük `DAG preview` bağlantısına tıklayarak süreçlerin nasıl bağlandığına dair yararlı bir diyagram alabilirsiniz.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Bu size süreçlerin nasıl bağlandığına ve ne ürettiklerine dair güzel bir genel bakış sunar.

Orijinal `sayHello` sürecine ek olarak, şimdi konsol çıktısında gördüğümüz süreçlerin adlarıyla eşleşen `convertToUpper` ve `collectGreetings` süreçlerine de sahibiz.
İki yeni süreç tanımı, `collectGreetings`'in `batch` adlı ek bir girdi parametresi alması ve iki çıktı üretmesi dışında `sayHello` süreci ile aynı şekilde yapılandırılmıştır.

Her birinin koduna ayrıntılı olarak girmeyeceğiz, ancak merak ediyorsanız, ayrıntılara [Hello Nextflow'un Bölüm 2'sinde](../hello_nextflow/03_hello_workflow.md) bakabilirsiniz.

Şimdilik, süreçlerin birbirine nasıl bağlandığını inceleyelim.

#### 2.3.2. Süreçler nasıl bağlanır

Burada bakılması gerçekten ilginç olan şey, süreç çağrılarının iş akışının `main:` bloğunda nasıl birbirine zincirlendiğidir.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // bir selamlama yayınla
    sayHello(greeting_ch)
    // selamlamayı büyük harfe dönüştür
    convertToUpper(sayHello.out)
    // tüm selamlamaları bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

İlk süreç çağrısının, `sayHello(greeting_ch)`, değişmediğini görebilirsiniz.
Ardından bir sonraki süreç çağrısı, `convertToUpper`'a, `sayHello`'nun çıktısına `sayHello.out` olarak atıfta bulunur.

Kalıp basittir: `processName.out` bir sürecin çıktı kanalına atıfta bulunur ve bu doğrudan bir sonraki sürece aktarılabilir.
Nextflow'da verileri bir adımdan diğerine bu şekilde aktarırız.

#### 2.3.3. Bir süreç birden fazla girdi alabilir

Üçüncü süreç çağrısı, `collectGreetings`'e yapılan çağrı, biraz farklıdır.

```groovy title="2b-multistep.nf" linenums="77"
    // tüm selamlamaları bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Bu çağrının iki girdi aldığını görüyorsunuz, `convertToUpper.out.collect()` ve `params.batch`.
`.collect()` kısmını şimdilik göz ardı edersek, bunu `collectGreetings(input1, input2)` olarak genelleştirebiliriz.

Bu, süreç modülündeki iki girdi bildirimiyle eşleşir:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Nextflow bunu ayrıştırdığında, çağrıdaki ilk girdiyi `path input_files`'a ve ikincisini `val batch_name`'e atayacaktır.

Artık bir sürecin birden fazla girdi alabileceğini ve çağrının iş akışı bloğunda nasıl göründüğünü biliyorsunuz.

Şimdi o ilk girdiye, `convertToUpper.out.collect()`'e daha yakından bakalım.

#### 2.3.4. `collectGreetings` çağrısında `collect()` ne yapar

`sayHello`'nun çıktısını `convertToUpper`'a aktarmak için, `sayHello`'nun çıktı kanalına basitçe `sayHello.out` olarak atıfta bulunduk. Ancak bir sonraki adım için, `convertToUpper.out.collect()` referansını görüyoruz.

Bu `collect()` kısmı nedir ve ne yapar?

Elbette bir operatördür. Daha önce karşılaştığımız `splitCsv` ve `map` operatörleri gibi.
Bu sefer operatör `collect` olarak adlandırılır ve `convertToUpper` tarafından üretilen çıktı kanalına uygulanır.

`collect` operatörü, aynı sürecin birden fazla çağrısından çıktıları toplamak ve bunları tek bir kanal öğesine paketlemek için kullanılır.

Bu iş akışı bağlamında, `convertToUpper.out` kanalındaki üç büyük harfli selamlamayı (bunlar üç ayrı kanal öğesidir ve normalde bir sonraki süreç tarafından ayrı çağrılarda işlenirdi) alıyor ve bunları tek bir öğeye paketliyor.
Tüm selamlamaları aynı dosyaya geri almanın yolu budur.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

Buna karşılık, `convertToUpper()` çıktısını `collectGreetings()`'e beslemeden önce `collect()` uygulamazsak, Nextflow basitçe `collectGreetings()`'i her selamlama üzerinde bağımsız olarak çalıştırır, bu da hedefimize ulaşmaz.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Süreç çağrıları arasında kanalların içeriğine dönüşümler uygulamak için kullanılabilecek birçok başka [operatör](https://nextflow.io/docs/latest/reference/operator.html) vardır.

Bu, boru hattı geliştiricilerine boru hattının akış mantığını özelleştirmek için çok fazla esneklik sağlar.
Dezavantajı, bazen boru hattının ne yaptığını çözmeyi zorlaştırabilmesidir.

#### 2.3.5. Bir girdi parametresi varsayılan bir değere sahip olabilir

`collectGreetings`'in ikinci bir girdi aldığını, `params.batch`'i fark etmiş olabilirsiniz:

```groovy title="2b-multistep.nf" linenums="77"
    // tüm selamlamaları bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Bu, iş akışına `--batch` adlı bir CLI parametresi geçirir.
Ancak, iş akışını daha önce başlattığımızda, bir `--batch` parametresi belirtmedik.

Burada neler oluyor?
`params` bloğuna bir göz atın:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

İş akışında yapılandırılmış bir varsayılan değer var, bu yüzden onu sağlamak zorunda değiliz.
Ancak komut satırında bir tane sağlarsak, belirttiğimiz değer varsayılan yerine kullanılacaktır.

Deneyin:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Komut çıktısı"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Özel toplu adınızla adlandırılmış yeni nihai çıktılar görmelisiniz.

??? abstract "Dizin içeriği"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Bu, Bölüm 3'te daha ayrıntılı olarak ele alacağımız girdi yapılandırmasının bir yönüdür, ancak şimdilik önemli olan girdi parametrelerine varsayılan değerler verilebileceğini bilmektir.

#### 2.3.6. Bir süreç birden fazla çıktı üretebilir

`collectGreetings` süreç tanımında, aşağıdaki çıktı bildirimlerini görüyoruz:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Bunlar daha sonra `publish:` bloğunda `emit:` ile verilen adla anılır:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Bu, belirli çıktıları iş akışındaki diğer süreçlere çeşitli operatörlerle birlikte ayrı ayrı aktarmayı kolaylaştırır.

#### 2.3.7. Yayınlanan çıktılar düzenlenebilir

`output` bloğunda, iş akışının yalnızca nihai çıktılarını seçmeyi kolaylaştırmak için ara sonuçları gruplamak için özel yollar kullandık.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Yayınlanan çıktıları düzenlemenin daha gelişmiş yolları vardır; yapılandırma bölümünde birkaçına değineceğiz.

!!! tip "İş akışları oluşturma hakkında daha fazla bilgi edinmek ister misiniz?"

    Çok adımlı iş akışları oluşturmanın ayrıntılı kapsamı için, [Hello Nextflow Bölüm 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md) bölümüne bakın.

### Özet

Çok adımlı iş akışlarının kanallar ve operatörler kullanılarak nasıl oluşturulduğunu ve nasıl çalıştıklarını temel düzeyde anlıyorsunuz.
Ayrıca süreçlerin birden fazla girdi alabileceğini ve birden fazla çıktı üretebileceğini ve bunların yapılandırılmış bir şekilde yayınlanabileceğini gördünüz.

### Sırada ne var?

Nextflow boru hatlarının kod yeniden kullanımını ve sürdürülebilirliği teşvik etmek için nasıl modülerleştirilebileceğini öğrenin.

---

## 3. Modülerleştirilmiş boru hatlarını çalıştırma

Şimdiye kadar baktığımız tüm iş akışları, tüm ilgili kodu içeren tek bir iş akışı dosyasından oluşuyordu.

Ancak, gerçek dünyadaki boru hatları genellikle _modülerleştirilmekten_, yani kodun farklı dosyalara bölünmesinden fayda sağlar.
Bu, geliştirilmelerini ve bakımlarını daha verimli ve sürdürülebilir hale getirebilir.

Burada Nextflow'da en yaygın kod modülerliği biçimini göstereceğiz, bu da **modüllerin** kullanımıdır.

Nextflow'da bir [**modül**](https://nextflow.io/docs/latest/module.html), bağımsız bir kod dosyasında kendi başına kapsüllenmiş tek bir süreç tanımıdır.
Bir iş akışında bir modül kullanmak için, iş akışı kod dosyanıza tek satırlık bir içe aktarma ifadesi eklemeniz yeterlidir; ardından süreci normalde yaptığınız gibi iş akışına entegre edebilirsiniz.
Bu, kodun birden fazla kopyasını üretmeden süreç tanımlarını birden fazla iş akışında yeniden kullanmayı mümkün kılar.

Şimdiye kadar tüm süreçlerini monolitik bir kod dosyasına dahil edilmiş iş akışları çalıştırıyorduk.
Şimdi süreçler ayrı modüllerde saklandığında nasıl göründüğünü göreceğiz.

Elbette bir kez daha gösteri amaçları için `2c-modules.nf` adlı uygun bir iş akışı ve `modules/` dizininde bulunan bir dizi modül hazırladık.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Dizin içeriği"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Her biri süreçlerden birinin adını taşıyan dört Nextflow dosyası olduğunu görüyorsunuz.
`cowpy.nf` dosyasını şimdilik göz ardı edebilirsiniz; ona daha sonra geleceğiz.

### 3.1. Kodu inceleme

Bu sefer önce koda bakacağız.
`2c-modules.nf` iş akışı dosyasını açarak başlayın.

??? full-code "Tam kod dosyası"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Boru hattı parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

İş akışı mantığının iş akışının önceki sürümüyle tamamen aynı olduğunu görüyorsunuz.
Ancak, süreç kodu iş akışı dosyasından gitmiş ve bunun yerine `modules` altındaki ayrı dosyalara işaret eden `include` ifadeleri var.

```groovy title="hello-modules.nf" linenums="3"
// Modülleri dahil et
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Bu dosyalardan birini açın ve ilgili sürecin kodunu bulacaksınız.

??? full-code "Tam kod dosyası"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' yazdırmak için echo kullan
    */
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

Gördüğünüz gibi, süreç kodu değişmedi; sadece ana iş akışı dosyasında olmak yerine ayrı bir modül dosyasına kopyalandı.
Aynı şey diğer iki süreç için de geçerlidir.

Şimdi bu yeni sürümü çalıştırmanın nasıl göründüğüne bakalım.

### 3.2. İş akışını çalıştırma

Terminalinizde bu komutu `-resume` bayrağıyla çalıştırın:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Süreç yürütmelerinin tümünün başarıyla önbelleğe alındığını fark edeceksiniz, bu da Nextflow'un kod bölünmüş ve ana iş akışı dosyası yeniden adlandırılmış olsa bile istenen işi zaten yaptığını tanıdığı anlamına gelir.

Bunların hiçbiri Nextflow için önemli değil; önemli olan, tüm kod bir araya getirildikten ve değerlendirildikten sonra oluşturulan iş betiğidir.

!!! tip

    Bir iş akışının bir bölümünü daha büyük bir boru hattına aktarılabilen bir 'alt iş akışı' olarak kapsüllemek de mümkündür, ancak bu bu kursun kapsamı dışındadır.

    Birleştirilebilir iş akışları geliştirme hakkında daha fazla bilgiyi [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/) Yan Görevinde öğrenebilirsiniz.

### Özet

Süreçlerin kod yeniden kullanımını teşvik etmek ve sürdürülebilirliği artırmak için bağımsız modüllerde nasıl saklanabileceğini biliyorsunuz.

### Sırada ne var?

Yazılım bağımlılıklarını yönetmek için konteynırları kullanmayı öğrenin.

---

## 4. Konteynırlaştırılmış yazılım kullanma

Şimdiye kadar örnek olarak kullandığımız iş akışları, ortamımızda varsayılan olarak bulunan UNIX araçlarını kullanarak çok temel metin işleme işlemlerini çalıştırmaya ihtiyaç duyuyordu.

Ancak, gerçek dünyadaki boru hatları genellikle çoğu ortamda varsayılan olarak dahil edilmeyen özel araçlar ve paketler gerektirir.
Normalde, bu araçları kurmanız, bağımlılıklarını yönetmeniz ve çakışmaları çözmeniz gerekir.

Bunların hepsi çok sıkıcı ve can sıkıcıdır.
Bu sorunu çözmenin çok daha iyi bir yolu **konteynırları** kullanmaktır.

Bir **konteyner**, kod, sistem kütüphaneleri ve ayarlar dahil olmak üzere bir uygulamayı çalıştırmak için gereken her şeyi içeren bir konteyner **imajından** oluşturulan hafif, bağımsız, yürütülebilir bir yazılım birimidir.

!!! Tip

    Bunu [Docker](https://www.docker.com/get-started/) teknolojisini kullanarak öğretiyoruz, ancak Nextflow başka birkaç konteyner teknolojisini de destekler.
    Nextflow'un konteynırlar için desteği hakkında daha fazla bilgiyi [buradan](https://nextflow.io/docs/latest/container.html) öğrenebilirsiniz.

### 4.1. Bir konteynırı doğrudan kullanma

İlk olarak, bir konteynırla doğrudan etkileşim kurmayı deneyelim.
Bu, Nextflow'da kullanmaya başlamadan önce konteynırların ne olduğuna dair anlayışınızı sağlamlaştırmaya yardımcı olacaktır.

#### 4.1.1. Konteyner imajını çekme

Bir konteyner kullanmak için, genellikle bir konteyner kayıt defterinden bir konteyner imajı indirirsiniz veya "çekersiniz" ve ardından bir konteyner örneği oluşturmak için konteyner imajını çalıştırırsınız.

Genel sözdizimi şu şekildedir:

```bash title="Sözdizimi"
docker pull '<container>'
```

- `docker pull` konteyner sistemine bir depodan bir konteyner imajı çekmesi talimatıdır.
- `'<container>'` konteyner imajının URI adresidir.

Örnek olarak, rastgele metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII sanatı üreten `cowsay` adlı bir aracın python uygulaması olan [cowpy](https://github.com/jeffbuttars/cowpy) içeren bir konteyner imajı çekelim.

Yayınlanmış konteynırları bulabileceğiniz çeşitli depolar vardır.
Bu Docker konteyner imajını `cowpy` Conda paketinden oluşturmak için [Seqera Containers](https://seqera.io/containers/) hizmetini kullandık: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Tam çekme komutunu çalıştırın:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Komut çıktısı"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Bu, sisteme belirtilen imajı indirmesini söyler.
İndirme tamamlandığında, konteyner imajının yerel bir kopyasına sahip olursunuz.

#### 4.1.2. Konteynırı başlatma

Konteynırlar tek seferlik bir komut olarak çalıştırılabilir, ancak bunları etkileşimli olarak da kullanabilirsiniz, bu size konteyner içinde bir kabuk istemi verir ve komutla oynamanıza izin verir.

Genel sözdizimi şu şekildedir:

```bash title="Sözdizimi"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` konteyner sistemine bir konteyner imajından bir konteyner örneği başlatması ve içinde bir komut yürütmesi talimatıdır.
- `--rm` sisteme komut tamamlandıktan sonra konteyner örneğini kapatmasını söyler.

Tamamen monte edilmiş konteyner yürütme komutu şöyle görünür:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Bu komutu çalıştırın ve isteminizin `(base) root@b645838b3314:/tmp#` gibi bir şeye değiştiğini görmelisiniz, bu da artık konteynırın içinde olduğunuzu gösterir.

Dizin içeriğini listelemek için `ls` çalıştırarak bunu doğrulayabilirsiniz:

```bash
ls /
```

??? success "Komut çıktısı"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Konteyner içindeki dosya sisteminin ana sistem dosya sisteminizden farklı olduğunu görüyorsunuz.

!!! Tip

    Bir konteyner çalıştırdığınızda, varsayılan olarak ana sistemden izole edilir.
    Bu, açıkça izin vermediğiniz sürece konteynırın ana sistemdeki herhangi bir dosyaya erişemeyeceği anlamına gelir, bunu `docker run` komutunun bir parçası olarak aşağıdaki sözdizimini kullanarak bir birim bağlamak istediğinizi belirterek yaparsınız:

    ```bash title="Sözdizimi"
    -v <outside_path>:<inside_path>
    ```

    Bu, dosya sisteminizin o kısmına erişmek için kullanabileceğiniz konteyner duvarından etkili bir şekilde bir tünel oluşturur.

    Bu, [Hello Nextflow'un Bölüm 5'inde](../hello_nextflow/05_hello_containers.md) daha ayrıntılı olarak ele alınmıştır.

#### 4.1.3. `cowpy` aracını çalıştırma

Konteynırın içinden, `cowpy` komutunu doğrudan çalıştırabilirsiniz.

```bash
cowpy "Hello Containers"
```

??? success "Komut çıktısı"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Bu, belirttiğimiz metni içeren bir konuşma balonu ile varsayılan inek karakterinin (veya 'cowacter') ASCII sanatını üretir.

Artık temel kullanımı test ettiğinize göre, ona bazı parametreler vermeyi deneyebilirsiniz.
Örneğin, araç belgeleri karakteri `-c` ile ayarlayabileceğimizi söylüyor.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Komut çıktısı"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Bu sefer ASCII sanat çıktısı Linux pengueni Tux'u gösteriyor, çünkü `-c tux` parametresini belirttik.

Konteynırın içinde olduğunuz için, sisteminizin kendisine herhangi bir kütüphane kurmak konusunda endişelenmeden cowpy komutunu istediğiniz kadar çalıştırabilir, girdi parametrelerini değiştirebilirsiniz.

??? tip "Diğer mevcut karakterler"

    Farklı bir karakter seçmek için '-c' bayrağını kullanın, bunlar dahil:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Bununla oynamaktan çekinmeyin.
İşiniz bittiğinde, `exit` komutunu kullanarak konteynırdan çıkın:

```bash
exit
```

Kendinizi normal kabuğunuzda bulacaksınız.

### 4.2. Bir iş akışında konteyner kullanma

Bir boru hattı çalıştırdığımızda, Nextflow'a her adımda hangi konteynırı kullanacağını söyleyebilmek istiyoruz ve önemlisi, az önce yaptığımız tüm işi yapmasını istiyoruz: konteynırı çekmek, başlatmak, komutu çalıştırmak ve işi bittiğinde konteynırı kapatmak.

İyi haber: Nextflow'un bizim için yapacağı tam olarak bu.
Her süreç için bir konteyner belirtmemiz yeterli.

Bunun nasıl çalıştığını göstermek için, üçüncü adımda üretilen toplanan selamlamalar dosyası üzerinde `cowpy` çalıştıran iş akışımızın başka bir sürümünü yaptık.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Bu, konuşma balonunda üç selamlamayı içeren ASCII sanatını içeren bir dosya çıktısı vermelidir.

#### 4.2.1. Kodu inceleme

İş akışı, `cowpy` çalıştırmak için ekstra adım artı öncekine çok benzer.

??? full-code "Tam kod dosyası"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Boru hattı parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
        // tüm selamlamaları bir dosyada topla
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy ile selamlamaların ASCII sanatını oluştur
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Bu iş akışının bir modül dosyasından bir `cowpy` süreci içe aktardığını ve bunu `collectGreetings()` çağrısının çıktısı artı `params.character` adlı bir girdi parametresi üzerinde çağırdığını görüyorsunuz.

```groovy title="2d-container.nf" linenums="31"
// cowpy ile selamlamaların ASCII sanatını oluştur
cowpy(collectGreetings.out.outfile, params.character)
```

ASCII sanatı oluşturmak için cowpy komutunu sarmalayan `cowpy` süreci, `cowpy.nf` modülünde tanımlanmıştır.

??? full-code "Tam kod dosyası"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // cowpy ile ASCII sanatı oluştur (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

`cowpy` süreci iki girdi gerektirir: konuşma balonuna konulacak metni içeren bir girdi dosyasının yolu (`input_file`) ve karakter değişkeni için bir değer.

Önemlisi, daha önce kullandığımız konteyner URI'sine işaret eden `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'` satırını da içerir.

#### 4.2.2. Docker'ın yapılandırmada etkin olduğunu kontrol etme

Eğitim kursunun Bölüm 3'ünü biraz önceden alarak, iş akışı yürütmesini yapılandırmak için Nextflow'un sunduğu ana yollardan biri olan `nextflow.config` yapılandırma dosyasını tanıtacağız.
Geçerli dizinde `nextflow.config` adlı bir dosya bulunduğunda, Nextflow onu otomatik olarak yükleyecek ve içerdiği herhangi bir yapılandırmayı uygulayacaktır.

Bu amaçla, Docker'ı etkinleştiren tek satırlık bir kodla bir `nextflow.config` dosyası ekledik.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Bu yapılandırma, Nextflow'a uyumlu bir konteyner belirten herhangi bir süreç için Docker kullanmasını söyler.

!!! tip

    Teknik olarak Docker yürütmesini komut satırından, çalıştırma başına `-with-docker <container>` parametresini kullanarak etkinleştirmek mümkündür.
    Ancak, bu bize tüm iş akışı için yalnızca bir konteyner belirtmemize izin verirken, size gösterdiğimiz yaklaşım süreç başına farklı bir konteyner belirtmemize izin verir.
    İkincisi modülerlik, kod bakımı ve tekrarlanabilirlik için çok daha iyidir.

#### 4.2.3. İş akışını çalıştırma

Özetlemek gerekirse, çalıştırmak üzere olduğumuz şey bu:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Sizce işe yarayacak mı?

İş akışını `-resume` bayrağıyla çalıştıralım ve karakterin hindi olmasını istediğimizi belirtelim.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

İlk üç adım önbelleğe alındı çünkü bunları daha önce çalıştırmıştık, ancak `cowpy` süreci yeni olduğu için gerçekten çalıştırılıyor.

`cowpy` adımının çıktısını `results` dizininde bulabilirsiniz.

??? abstract "Dosya içeriği"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Karakterin tüm selamlamaları söylediğini görüyorsunuz, çünkü toplanan büyük harfli selamlamalar dosyası üzerinde çalıştı.

Daha da önemlisi, bunu cowpy'nin ve tüm bağımlılıklarının düzgün bir kurulumunu yapmak zorunda kalmadan boru hattımızın bir parçası olarak çalıştırabildik.
Ve artık boru hattını işbirlikçilerimizle paylaşabilir ve onların da Docker veya yukarıda belirtildiği gibi alternatiflerinden birini (Singularity/Apptainer gibi) kurmak dışında herhangi bir şey kurmalarına gerek kalmadan altyapılarında çalıştırmalarını sağlayabiliriz.

#### 4.2.4. Nextflow'un konteynırlaştırılmış görevi nasıl başlattığını inceleme

Bu bölümün son bir kodası olarak, Nextflow'un konteynırlarla nasıl çalıştığına dair biraz daha içgörü elde etmek için `cowpy` süreç çağrılarından birinin çalışma alt dizinine bir göz atalım.

`cowpy` süreci için çalışma alt dizininin yolunu bulmak için `nextflow run` komutunuzdan çıktıyı kontrol edin.
Yukarıda gösterilen çalıştırma için aldığımıza bakarsak, `cowpy` süreci için konsol günlük satırı `[7f/caf718]` ile başlıyor.
Bu, aşağıdaki kısaltılmış dizin yoluna karşılık gelir: `work/7f/caf718`.

O dizinde, Nextflow'un boru hattını yürütme sürecinde sizin adınıza çalıştırdığı tüm komutları içeren `.command.run` dosyasını bulacaksınız.

??? abstract "Dosya içeriği"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # girdi dosyalarını hazırla
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Bu dosyada `nxf_launch` araması yaparsanız, şuna benzer bir şey görmelisiniz:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Bu başlatma komutu, Nextflow'un manuel olarak çalıştırdığımızda yaptığımıza çok benzer bir `docker run` komutu kullanarak süreç çağrısını başlattığını gösterir.
Ayrıca ilgili çalışma alt dizinini konteynıra bağlar, konteyner içindeki çalışma dizinini buna göre ayarlar ve `.command.sh` dosyasındaki şablonlu bash betiğimizi çalıştırır.

Bu, önceki bölümde manuel olarak yapmak zorunda kaldığımız tüm zor işin artık Nextflow tarafından bizim için yapıldığını doğrular!

### Özet

Konteynırların yazılım aracı sürümlerini yönetmede ve tekrarlanabilirliği sağlamada ne rol oynadığını anlıyorsunuz.

Daha genel olarak, gerçek dünyadaki Nextflow boru hatlarının temel bileşenlerinin neler olduğuna ve nasıl organize edildiklerine dair temel bir anlayışa sahipsiniz.
Nextflow'un birden fazla girdiyi verimli bir şekilde nasıl işleyebileceğinin, birbirine bağlı birden fazla adımdan oluşan iş akışlarını nasıl çalıştırabileceğinin, modüler kod bileşenlerinden nasıl yararlanabileceğinin ve daha fazla tekrarlanabilirlik ve taşınabilirlik için konteynırları nasıl kullanabileceğinin temellerini biliyorsunuz.

### Sırada ne var?

Bir mola daha verin! Bu, Nextflow boru hatlarının nasıl çalıştığına dair büyük bir bilgi yığınıydı.

Bu eğitimin son bölümünde, yapılandırma konusunu daha derinlemesine inceleyeceğiz.
Boru hattınızın yürütülmesini altyapınıza uyacak şekilde nasıl yapılandıracağınızı ve girdilerin ve parametrelerin yapılandırmasını nasıl yöneteceğinizi öğreneceksiniz.

---

## Quiz

<quiz>
Nextflow neden her süreç çağrısı için ayrı bir görev dizini oluşturur?
- [ ] Yürütme hızını artırmak için
- [ ] Bellek kullanımını azaltmak için
- [x] Yürütmeleri izole etmek ve çıktılar arasındaki çakışmaları önlemek için
- [ ] Paralel dosya sıkıştırmayı etkinleştirmek için

Daha fazla bilgi: [1.3. Orijinal çıktıları ve günlükleri bulma](#13-orijinal-çıktıları-ve-günlükleri-bulma)
</quiz>

<quiz>
Bir iş akışı çalıştırırken `-ansi-log false` seçeneği ne yapar?
- [ ] Tüm konsol çıktısını devre dışı bırakır
- [x] Çıktıdan rengi kaldırır
- [x] Tüm görev dizini yollarını tek bir satırda yoğunlaştırmak yerine gösterir
- [ ] Ayrıntılı hata ayıklama modunu etkinleştirir

Daha fazla bilgi: [1.3.2. Terminalin daha fazla ayrıntı göstermesini sağlama](#132-terminalin-daha-fazla-ayrıntı-göstermesini-sağlama)

Bu stili tercih ediyorsanız aşağıdaki ortam değişkenlerinden birini de kullanabilirsiniz:

```bash
export NXF_ANSI_LOG=0
# veya
export NO_COLOR=1
```

</quiz>

<quiz>
`#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }` kodunda, `#!groovy .map { line -> line[0] }` ne yapar?
- [ ] Boş satırları filtreler
- [ ] Satırları alfabetik olarak sıralar
- [x] Her CSV satırından ilk sütunu çıkarır
- [ ] Satır sayısını sayar

Daha fazla bilgi: [1.4.1. CSV'den girdi verilerini yükleme](#141-csvden-girdi-verilerini-yükleme)
</quiz>

<quiz>
Çıktı dosya adlarına girdi değerini dahil etmek neden önemlidir (örn. `#!groovy "${greeting}-output.txt"`)?
- [ ] İşleme hızını artırmak için
- [ ] Devam etme işlevselliğini etkinleştirmek için
- [x] Birden fazla girdi işlenirken çıktı dosyalarının birbirinin üzerine yazılmasını önlemek için
- [ ] Dosyaları sıkıştırmayı kolaylaştırmak için

Daha fazla bilgi: [1.4.3. Çıktıların nasıl adlandırıldığı](#143-çıktıların-nasıl-adlandırıldığı)
</quiz>

<quiz>
Modülerleştirilmiş bir iş akışında `include` ifadesinin amacı nedir?
- [ ] Süreç kodunu iş akışı dosyasına kopyalamak
- [x] Harici bir modül dosyasından bir süreç tanımını içe aktarmak
- [ ] Yapılandırma ayarlarını dahil etmek
- [ ] Belge yorumları eklemek

Daha fazla bilgi: [3. Modülerleştirilmiş boru hatlarını çalıştırma](#3-modülerleştirilmiş-boru-hatlarını-çalıştırma)
</quiz>

<quiz>
Bir iş akışını modülerleştirdiğinizde ve `-resume` ile çalıştırdığınızda ne olur?
- [ ] Modüler süreçler için önbelleğe alma devre dışı bırakılır
- [ ] Tüm görevler yeniden yürütülmelidir
- [x] Önbelleğe alma, oluşturulan iş betiklerine göre normal şekilde çalışır
- [ ] Yalnızca ana iş akışı dosyası önbelleğe alınır

Daha fazla bilgi: [3.2. İş akışını çalıştırma](#32-iş-akışını-çalıştırma)
</quiz>

<quiz>
Bir süreç tanımındaki `container` yönergesi neyi belirtir?
- [ ] Süreç için çalışma dizini
- [ ] Maksimum bellek tahsisi
- [x] Süreci çalıştırmak için kullanılacak konteyner imajı URI'si
- [ ] Çıktı dosyası formatı

Daha fazla bilgi: [4.2. Bir iş akışında konteyner kullanma](#42-bir-iş-akışında-konteyner-kullanma)
</quiz>

<quiz>
`.command.run` dosyasında, `nxf_launch` işlevi ne içerir?
- [ ] Nextflow sürüm bilgisi
- [ ] İş akışı parametreleri
- [x] Birim bağlamaları ve konteyner ayarları ile `docker run` komutu
- [ ] Süreç girdi bildirimleri

Daha fazla bilgi: [4.2.4. Nextflow'un konteynırlaştırılmış görevi nasıl başlattığını inceleme](#424-nextflowun-konteynırlaştırılmış-görevi-nasıl-başlattığını-inceleme)
</quiz>

<quiz>
Nextflow konteynırlaştırılmış bir süreç çalıştırırken otomatik olarak ne yapar? (Tümünü seçin)
- [x] Gerekirse konteyner imajını çeker
- [x] Çalışma dizinini konteynıra bağlar
- [x] Süreç betiğini konteynırın içinde çalıştırır
- [x] Yürütmeden sonra konteyner örneğini temizler

Daha fazla bilgi: [4. Konteynırlaştırılmış yazılım kullanma](#4-konteynırlaştırılmış-yazılım-kullanma)
</quiz>
