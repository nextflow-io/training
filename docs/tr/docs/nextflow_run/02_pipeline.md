# Bölüm 2: Gerçek pipeline'ları çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu kursun 1. Bölümünde (Temel İşlemleri Çalıştırma), kod karmaşıklığını düşük tutmak için yalnızca minimal özelliklere sahip örnek bir workflow ile başladık.
Örneğin, `1-hello.nf` bir seferde tek bir değer sağlamak için komut satırı parametresi (`--input`) kullandı.

Ancak, gerçek dünya pipeline'larının çoğu, büyük miktarda veriyi ölçekte verimli bir şekilde işlemek ve bazen karmaşık mantıkla birbirine zincirlenen birden fazla işleme adımı uygulamak için daha sofistike özellikler kullanır.

Eğitimin bu bölümünde, orijinal Hello World pipeline'ının genişletilmiş versiyonlarını deneyerek gerçek dünya pipeline'larının temel özelliklerini gösteriyoruz.

## 1. Bir dosyadan girdi verilerini işleme

Gerçek dünya pipeline'larında, tipik olarak bir veya daha fazla girdi dosyasında bulunan birden fazla veri noktasını (veya veri serisini) işlemek istiyoruz.
Ve mümkün olduğunda, analiz için harcanan süreyi kısaltmak için bağımsız verilerin işlenmesini paralel olarak çalıştırmak istiyoruz.

Nextflow'un bunu nasıl yaptığını göstermek için, gerçek bir veri analizinde işlemek isteyebileceğiniz türden sütunlu verileri taklit eden birkaç girdi selamlama içeren `greetings.csv` adlı bir CSV dosyası hazırladık.
Sayıların anlamlı olmadığını, sadece açıklama amaçlı olduğunu unutmayın.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Ayrıca, CSV dosyasını okuyacak, selamlamaları çıkaracak ve her birini ayrı bir dosyaya yazacak olan `2a-inputs.nf` adlı orijinal workflow'un geliştirilmiş bir versiyonunu yazdık.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Önce workflow'u çalıştıralım ve ardından ilgili Nextflow koduna bakacağız.

### 1.1. Workflow'u çalıştırın

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

Heyecan verici bir şekilde, bu process için '3 of 3' çağrının yapıldığını gösteriyor gibi görünüyor, bu cesaret verici, çünkü girdi olarak sağladığımız CSV'de üç veri satırı vardı.
Bu, `sayHello()` process'inin her girdi satırında bir kez olmak üzere üç kez çağrıldığını gösteriyor.

### 1.2. Yayınlanan çıktıları `results` dizininde bulun

Workflow'umuzun hala çıktılarımızın bir kopyasını oraya yazıp yazmadığını görmek için 'results' dizinine bakalım.

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

Evet! Yeterince uygun bir şekilde farklı adlara sahip üç çıktı dosyası içeren `2a-inputs` adlı yeni bir dizin görüyoruz.

Uygun selamlama dizesini içerdiklerinden emin olmak için her birini açabilirsiniz.

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

### 1.3. Orijinal çıktıları ve logları bulun

Yukarıdaki konsol çıktısının yalnızca bir görev dizinine atıfta bulunduğunu fark etmiş olabilirsiniz.
Bu, `sayHello()` çağrısının üçünün de bu tek görev dizininde yürütüldüğü anlamına mı geliyor?

#### 1.3.1. Terminalde verilen görev dizinini inceleyin

`8e/0eb066` görev dizininin içine bakalım.

??? abstract "Dizin içeriği"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Yalnızca selamlamalardan birine karşılık gelen çıktıyı buluyoruz (gizli dosyaların görüntülenmesini etkinleştirirsek yardımcı dosyaları da).

Peki ne oluyor?

Varsayılan olarak, ANSI loglama sistemi aynı process'e yapılan tüm çağrılar için durum bilgisini aynı satıra yazar.
Sonuç olarak, konsol çıktısında bize üç görev dizini yolundan yalnızca birini (`8e/0eb066`) gösterdi.
Orada listelenmeyen iki tane daha var.

#### 1.3.2. Terminalin daha fazla ayrıntı göstermesini sağlayın

Process çağrılarının tam listesini görmek için loglama davranışını aşağıdaki gibi komuta `-ansi-log false` ekleyerek değiştirebiliriz:

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

Bu sefer çıktıda üç process çalıştırmasını ve ilişkili çalışma alt dizinlerini görüyoruz.
ANSI loglamasını devre dışı bırakmak, Nextflow'un terminal çıktısında renk kullanmasını da engelledi.

Durumun iki loglama modu arasında biraz farklı raporlandığına dikkat edin.
Yoğunlaştırılmış modda, Nextflow çağrıların başarıyla tamamlanıp tamamlanmadığını raporlar.
Bu genişletilmiş modda, yalnızca gönderildiklerini raporlar.

Bu, `sayHello()` process'inin üç kez çağrıldığını ve her biri için ayrı bir görev dizini oluşturulduğunu doğrular.

Orada listelenen her görev dizininin içine bakarsak, her birinin selamlamalardan birine karşılık geldiğini doğrulayabiliriz.

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

Bu, her process çağrısının diğerlerinden izole olarak yürütüldüğünü doğrular.
Bunun, process'in benzersiz olmayan adlara sahip ara dosyalar üretmesi durumunda çakışmaları önlemek de dahil olmak üzere birçok avantajı vardır.

!!! tip "İpucu"

    Karmaşık bir workflow veya çok sayıda girdi için, tam listenin terminale çıktı verilmesi biraz bunaltıcı olabilir, bu nedenle insanlar normalde rutin kullanımda `-ansi-log false` kullanmazlar.

### 1.4. Workflow kodunu inceleyin

Bu workflow versiyonu bir CSV girdi dosyasını okuyabilir, girdileri ayrı ayrı işleyebilir ve çıktıları benzersiz şekilde adlandırabilir.

Workflow kodunda bunu mümkün kılan şeylere bir göz atalım.

??? full-code "Tam kod dosyası"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' ifadesini bir dosyaya yazdırmak için echo kullan
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
    * Pipeline parametreleri
    */
    params {
        input: Path
    }

    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
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

Yine, kod sözdizimini ezberlememeniz gerekmez, ancak önemli işlevsellik sağlayan workflow'un temel bileşenlerini tanımayı öğrenmek iyidir.

#### 1.4.1. CSV'den girdi verilerini yükleme

En ilginç kısım şudur: komut satırından tek bir değer almaktan, bir CSV dosyası almaya, ayrıştırmaya ve içerdiği bireysel selamlamaları işlemeye nasıl geçtik?

Nextflow'da bunu bir [**channel**](https://nextflow.io/docs/latest/channel.html) ile yapıyoruz: girdileri verimli bir şekilde işlemek ve çok adımlı workflow'larda bir adımdan diğerine taşımak için tasarlanmış, yerleşik paralellik ve birçok ek avantaj sağlayan bir yapı.

Parçalayalım.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // bir CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // bir selamlama yayınla
    sayHello(greeting_ch)
```

Bu kod, CSV dosyasını okuyan, ayrıştıran ve her satırdan ilk sütunu çıkaran `greeting_ch` adlı bir channel oluşturur.
Sonuç, `Hello`, `Bonjour` ve `Holà` içeren bir channel'dır.

??? tip "Bu nasıl çalışır?"

    Bu satırın düz Türkçe'de anlamı şudur:

    - `channel.fromPath`, dosya yollarından bir channel oluşturan bir **channel factory**'dir
    - `(params.input)`, dosya yolunun komut satırında `--input` ile sağlandığını belirtir

    Başka bir deyişle, bu satır Nextflow'a şunu söyler: `--input` ile verilen dosya yolunu al ve içeriğini girdi verisi olarak işlemeye hazırlan.

    Sonra, sonraki iki satır dosyanın gerçek ayrıştırmasını yapan ve verileri uygun veri yapısına yükleyen **operatör**leri uygular:

    - `.splitCsv()` Nextflow'a CSV dosyasını satırları ve sütunları temsil eden bir diziye ayrıştırmasını söyler
    - `.map { line -> line[0] }` Nextflow'a her satırdan yalnızca ilk sütundaki elemanı almasını söyler

    Yani pratikte, aşağıdaki CSV dosyasından başlayarak:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Bunu şuna benzer bir diziye dönüştürdük:

    ```txt title="Dizi içeriği"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    Ve sonra üç satırın her birinden ilk elemanı aldık ve artık `Hello`, `Bonjour` ve `Holà` içeren bir Nextflow channel'ına yükledik.

    Channel'ları ve operatörleri derinlemesine anlamak istiyorsanız, bunları kendiniz nasıl yazacağınız da dahil, [Hello Nextflow Bölüm 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) bölümüne bakın.

#### 1.4.2. Her selamlamada process'i çağırma

Sonra, workflow'un `main:` bloğunun son satırında, yüklenen `greeting_ch` channel'ını `sayHello()` process'ine girdi olarak sağlıyoruz.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // bir CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // bir selamlama yayınla
    sayHello(greeting_ch)
```

Bu, Nextflow'a process'i channel'daki her eleman üzerinde, yani her selamlamada bireysel olarak çalıştırmasını söyler.
Ve Nextflow bu kadar akıllı olduğundan, mevcut hesaplama altyapısına bağlı olarak bu process çağrılarını mümkünse paralel olarak çalıştıracaktır.

Çok sayıda veriyi (birçok örnek veya veri noktası, araştırma biriminiz ne olursa olsun) karşılaştırmalı olarak çok az kodla verimli ve ölçeklenebilir bir şekilde bu şekilde işleyebilirsiniz.

#### 1.4.3. Çıktılar nasıl adlandırılır

Son olarak, çıktı dosyalarının benzersiz şekilde adlandırılmasını nasıl sağladığımızı görmek için process koduna hızlıca bakmaya değer.

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

Gördüğünüz gibi, `1-hello.nf`'deki bu process versiyonuna kıyasla, çıktı bildirimi ve komutun ilgili kısmı, selamlama değerini çıktı dosya adına dahil edecek şekilde değişti.

Bu, ortak results dizinine yayınlandıklarında çıktı dosya adlarının çakışmamasını sağlamanın bir yoludur.

Ve process bildiriminde yapmak zorunda olduğumuz tek değişiklik bu!

### Özet

Channel'ların ve operatörlerin birden fazla girdiyi verimli bir şekilde işlememizi nasıl sağladığını temel düzeyde anlıyorsunuz.

### Sırada ne var?

Çok adımlı workflow'ların nasıl oluşturulduğunu ve nasıl çalıştığını keşfedin.

---

## 2. Çok adımlı workflow'ları çalıştırma

Gerçek dünya workflow'larının çoğu birden fazla adım içerir.
Channel'lar hakkında öğrendiklerimizi temel alarak, Nextflow'un çok adımlı bir workflow'da process'leri birbirine bağlamak için channel'ları ve operatörleri nasıl kullandığına bakalım.

Bu amaçla, aşağıdakileri gösteren, üç ayrı adımı birbirine zincirleyen örnek bir workflow sağlıyoruz:

1. Bir process'ten diğerine veri akışı sağlama
2. Birden fazla process çağrısından gelen çıktıları tek bir process çağrısında toplama

Özellikle, her girdi selamlamasını alan, büyük harfe dönüştüren ve ardından tüm büyük harfli selamlamaları tek bir çıktı dosyasında toplayan `2b-multistep.nf` adlı genişletilmiş bir workflow versiyonu yaptık.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Daha önce olduğu gibi, önce workflow'u çalıştıracağız, ardından yeni olan şeyi görmek için koda bakacağız.

### 2.1. Workflow'u çalıştırın

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

Gördüğünüz gibi, vaat edildiği gibi, workflow'un bir parçası olarak birden fazla adım çalıştırıldı; ilk ikisi (`sayHello` ve `convertToUpper`) muhtemelen her bir selamlamada çalıştırıldı ve üçüncüsü (`collectGreetings`) yalnızca bir kez, üç `convertToUpper` çağrısının çıktılarında çalıştırıldı.

### 2.2. Çıktıları bulun

Gerçekten olan şeyin bu olduğunu `results` dizinine bakarak doğrulayalım.

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

Gördüğünüz gibi, `2b-multistep` adında yeni bir dizinimiz var ve öncekinden çok daha fazla dosya içeriyor.
Dosyaların bazıları `intermediates` adlı bir alt dizinde gruplandırılmış, iki dosya ise üst düzeyde bulunuyor.

Bu ikisi, çok adımlı workflow'un nihai sonuçlarıdır.
Dosya adlarına bakıp beklediğiniz gibi olduklarını doğrulamak için içeriklerini kontrol etmek üzere bir dakikanızı ayırın.

??? abstract "Dosya içeriği"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

İlki, vaat edildiği gibi büyük harfe dönüştürülmüş ve tek bir dosyada toplanmış üç selamlamamızı içeriyor.
İkincisi, çalıştırma hakkında bazı bilgileri özetleyen bir rapor dosyasıdır.

### 2.3. Kodu inceleyin

Koda bakalım ve çok adımlı workflow'lar için temel kalıpları tanımlayalım.

??? full-code "Tam kod dosyası"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' ifadesini bir dosyaya yazdırmak için echo kullan
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
    * Pipeline parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
        // tüm selamlamaları tek bir dosyada topla
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

Orada çok şey oluyor, ancak önceki workflow versiyonuna kıyasla en belirgin fark, şimdi birden fazla process tanımının olması ve buna bağlı olarak workflow bloğunda birkaç process çağrısının olmasıdır.

Daha yakından bakalım ve en ilginç parçaları tanımlayabilecek miyiz görelim.

#### 2.3.1. Workflow yapısını görselleştirme

Nextflow uzantısıyla VSCode kullanıyorsanız, herhangi bir Nextflow betiğinde workflow bloğunun hemen üzerinde görüntülenen küçük `DAG preview` bağlantısına tıklayarak process'lerin nasıl bağlandığına dair yararlı bir diyagram alabilirsiniz.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Bu size process'lerin nasıl bağlandığına ve ne ürettiklerine dair güzel bir genel bakış sunar.

Orijinal `sayHello` process'ine ek olarak, şimdi konsol çıktısında gördüğümüz process'lerin adlarıyla eşleşen `convertToUpper` ve `collectGreetings`'in de olduğunu görüyorsunuz.
İki yeni process tanımı, `sayHello` process'i ile aynı şekilde yapılandırılmıştır, ancak `collectGreetings` `batch` adlı ek bir girdi parametresi alır ve iki çıktı üretir.

Her birinin koduna ayrıntılı olarak girmeyeceğiz, ancak merak ediyorsanız, ayrıntıları [Hello Nextflow Bölüm 3](../hello_nextflow/03_hello_workflow.md)'te bulabilirsiniz.

Şimdilik, process'lerin birbirine nasıl bağlandığını inceleyelim.

#### 2.3.2. Process'ler nasıl bağlanır

Burada bakılması gereken gerçekten ilginç şey, process çağrılarının workflow'un `main:` bloğunda nasıl birbirine zincirlendiğidir.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // bir CSV dosyasından girdiler için bir kanal oluştur
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // bir selamlama yayınla
    sayHello(greeting_ch)
    // selamlamayı büyük harfe dönüştür
    convertToUpper(sayHello.out)
    // tüm selamlamaları tek bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

İlk process çağrısı olan `sayHello(greeting_ch)`'nin değişmediğini görebilirsiniz.
Sonra, `convertToUpper`'a yapılan sonraki process çağrısı, `sayHello`'nun çıktısına `sayHello.out` olarak atıfta bulunur.

Kalıp basit: `processName.out` bir process'in çıktı channel'ına atıfta bulunur ve bu doğrudan sonraki process'e geçirilebilir.
Nextflow'da verileri bir adımdan diğerine bu şekilde taşıyoruz.

#### 2.3.3. Bir process birden fazla girdi alabilir

Üçüncü process çağrısı olan `collectGreetings` biraz farklıdır.

```groovy title="2b-multistep.nf" linenums="77"
    // tüm selamlamaları tek bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Bu çağrıya `convertToUpper.out.collect()` ve `params.batch` olmak üzere iki girdi verildiğini görüyorsunuz.
Şimdilik `.collect()` kısmını göz ardı edersek, bunu `collectGreetings(input1, input2)` olarak genelleyebiliriz.

Bu, process modülündeki iki girdi bildirimiyle eşleşir:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Nextflow bunu ayrıştırdığında, çağrıdaki ilk girdiyi `path input_files`'a, ikinciyi `val batch_name`'e atayacaktır.

Böylece artık bir process'in birden fazla girdi alabileceğini ve workflow bloğundaki çağrının nasıl göründüğünü biliyorsunuz.

Şimdi ilk girdiye daha yakından bakalım, `convertToUpper.out.collect()`.

#### 2.3.4. `collectGreetings` çağrısında `collect()` ne yapar

`sayHello`'nun çıktısını `convertToUpper`'a geçirmek için, `sayHello`'nun çıktı channel'ına `sayHello.out` olarak atıfta bulunduk. Ancak bir sonraki adım için `convertToUpper.out.collect()` referansı görüyoruz.

Bu `collect()` kısmı nedir ve ne yapar?

Bu elbette bir operatördür. Daha önce karşılaştığımız `splitCsv` ve `map` operatörleri gibi.
Bu sefer operatör `collect` olarak adlandırılır ve `convertToUpper` tarafından üretilen çıktı channel'ına uygulanır.

`collect` operatörü, aynı process'e yapılan birden fazla çağrıdan gelen çıktıları toplamak ve tek bir channel elemanına paketlemek için kullanılır.

Bu workflow bağlamında, `convertToUpper.out` channel'ındaki üç büyük harfli selamlamayı --bunlar üç ayrı channel öğesidir ve normalde sonraki process tarafından ayrı çağrılarda işlenirler-- alıp tek bir öğeye paketliyor.

Daha pratik terimlerle: `collectGreetings()`'e beslemeden önce `convertToUpper()`'ın çıktısına `collect()` uygulamasaydık, Nextflow basitçe `collectGreetings()`'i her selamlamada bağımsız olarak çalıştırırdı, bu da amacımıza ulaşmazdı.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Buna karşılık, `collect()` kullanmak, workflow'un ikinci adımı tarafından üretilen tüm ayrı büyük harfli selamlamaları almamıza ve pipeline'ın üçüncü adımında tek bir çağrıya hepsini birlikte beslememize olanak tanır.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

Tüm selamlamaları aynı dosyaya bu şekilde geri alıyoruz.

Process çağrıları arasında channel içeriklerine dönüşümler uygulamak için birçok başka [operatör](https://www.nextflow.io/docs/latest/reference/operator.html#operator-page) mevcuttur.

Bu, pipeline geliştiricilerine pipeline'larının akış mantığını özelleştirmek için çok fazla esneklik sağlar.
Dezavantajı, bazen pipeline'ın ne yaptığını çözmeyi zorlaştırabilmesidir.

#### 2.3.5. Bir girdi parametresinin varsayılan değeri olabilir

`collectGreetings`'in ikinci bir girdi aldığını, `params.batch`'i fark etmiş olabilirsiniz:

```groovy title="2b-multistep.nf" linenums="77"
    // tüm selamlamaları tek bir dosyada topla
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Bu, workflow'a `--batch` adlı bir CLI parametresi geçirir.
Ancak, workflow'u daha önce başlattığımızda bir `--batch` parametresi belirtmedik.

Burada ne oluyor?
`params` bloğuna bir göz atın:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Workflow'da yapılandırılmış bir varsayılan değer var, bu yüzden sağlamamız gerekmiyor.
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

Özel batch adınızla adlandırılmış yeni nihai çıktılar görmelisiniz.

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

#### 2.3.6. Bir process birden fazla çıktı üretebilir

`collectGreetings` process tanımında aşağıdaki çıktı bildirimlerini görüyoruz:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Bunlara daha sonra `publish:` bloğunda `emit:` ile verilen adla atıfta bulunulur:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Bu, belirli çıktıları çeşitli operatörlerle birlikte workflow'daki diğer process'lere bireysel olarak geçirmeyi kolaylaştırır.

#### 2.3.7. Yayınlanan çıktılar organize edilebilir

`output` bloğunda, workflow'un yalnızca nihai çıktılarını seçmeyi kolaylaştırmak için ara sonuçları gruplamak üzere özel yollar kullandık.

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

Yayınlanan çıktıları organize etmenin daha sofistike yolları var; yapılandırma bölümünde birkaçına değineceğiz.

!!! tip "Workflow oluşturma hakkında daha fazla bilgi edinmek ister misiniz?"

    Çok adımlı workflow'lar oluşturmanın ayrıntılı ele alınması için [Hello Nextflow Bölüm 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md) bölümüne bakın.

### Özet

Çok adımlı workflow'ların channel'lar ve operatörler kullanılarak nasıl oluşturulduğunu ve nasıl çalıştığını temel düzeyde anlıyorsunuz.
Ayrıca process'lerin birden fazla girdi alabileceğini ve birden fazla çıktı üretebileceğini ve bunların yapılandırılmış bir şekilde yayınlanabileceğini gördünüz.

### Sırada ne var?

Nextflow pipeline'larının kod yeniden kullanımını ve sürdürülebilirliği teşvik etmek için nasıl modülerleştirilebileceğini öğrenin.

---

## 3. Modülerleştirilmiş pipeline'ları çalıştırma

Şu ana kadar baktığımız tüm workflow'lar, tüm ilgili kodu içeren tek bir workflow dosyasından oluşuyordu.

Ancak, gerçek dünya pipeline'ları tipik olarak _modülerleştirilmiş_ olmaktan fayda görür, yani kod farklı dosyalara bölünür.
Bu, geliştirmelerini ve bakımlarını daha verimli ve sürdürülebilir hale getirebilir.

Burada Nextflow'da en yaygın kod modülerlik biçimini göstereceğiz, bu da **modül** kullanımıdır.

Nextflow'da bir [**modül**](https://nextflow.io/docs/latest/module.html), bağımsız bir kod dosyasında kendi başına kapsüllenmiş tek bir process tanımıdır.
Bir workflow'da modül kullanmak için, workflow kod dosyanıza tek satırlık bir import ifadesi eklemeniz yeterlidir; ardından process'i normalde yapacağınız şekilde workflow'a entegre edebilirsiniz.
Bu, kodun birden fazla kopyasını üretmeden birden fazla workflow'da process tanımlarını yeniden kullanmayı mümkün kılar.

Şu ana kadar tüm process'leri monolitik bir kod dosyasına dahil edilmiş workflow'ları çalıştırıyorduk.
Şimdi process'ler ayrı modüllerde saklandığında nasıl göründüğünü göreceğiz.

Elbette yine gösterim amacıyla uygun bir workflow, `2c-modules.nf` adlı, `modules/` dizininde bulunan bir modül setiyle birlikte hazırladık.

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

Her biri process'lerden birinin adını taşıyan dört Nextflow dosyası olduğunu görüyorsunuz.
Şimdilik `cowpy.nf` dosyasını görmezden gelebilirsiniz; ona daha sonra geleceğiz.

### 3.1. Kodu inceleyin

Bu sefer önce koda bakacağız.
`2c-modules.nf` workflow dosyasını açarak başlayın.

??? full-code "Tam kod dosyası"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
        // tüm selamlamaları tek bir dosyada topla
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

Workflow mantığının önceki workflow versiyonuyla tamamen aynı olduğunu görüyorsunuz.
Ancak, process kodu workflow dosyasından gitmiş ve yerine `modules` altındaki ayrı dosyalara işaret eden `include` ifadeleri var.

```groovy title="hello-modules.nf" linenums="3"
// Modülleri dahil et
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Bu dosyalardan birini açın ve ilgili process'in kodunu bulacaksınız.

??? full-code "Tam kod dosyası"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * 'Hello World!' ifadesini bir dosyaya yazdırmak için echo kullan
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

Gördüğünüz gibi, process kodu değişmedi; sadece ana workflow dosyasında olmak yerine bireysel bir modül dosyasına kopyalandı.
Aynısı diğer iki process için de geçerlidir.

Şimdi bu yeni versiyonu çalıştırmanın nasıl göründüğüne bakalım.

### 3.2. Workflow'u çalıştırın

Bu komutu `-resume` bayrağıyla terminalinizde çalıştırın:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [j6/cdfa66] sayHello (1)       | 3 of 3, cached: ✔
    [95/79484f] convertToUpper (2) | 3 of 3, cached: ✔
    [5e/4358gc] collectGreetings   | 1 of 1, cached: ✔
    ```

Process çalıştırmalarının hepsinin başarıyla önbelleğe alındığını fark edeceksiniz, yani Nextflow istenen işi zaten yaptığını tanıdı, kod bölünmüş ve ana workflow dosyası yeniden adlandırılmış olsa bile.

Bunların hiçbiri Nextflow için önemli değil; önemli olan tüm kod bir araya getirilip değerlendirildikten sonra oluşturulan iş betiğidir.

!!! tip "İpucu"

    Bir workflow'un bir bölümünü daha büyük bir pipeline'a aktarılabilecek bir 'alt workflow' olarak kapsüllemek de mümkündür, ancak bu kursun kapsamı dışındadır.

    [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/) Yan Görevinde birleştirilebilir workflow'lar geliştirme hakkında daha fazla bilgi edinebilirsiniz.

### Özet

Process'lerin kod yeniden kullanımını teşvik etmek ve sürdürülebilirliği artırmak için bağımsız modüllerde nasıl saklanabileceğini biliyorsunuz.

### Sırada ne var?

Yazılım bağımlılıklarını yönetmek için konteynerleri kullanmayı öğrenin.

---

## 4. Konteynerleştirilmiş yazılım kullanma

Şu ana kadar örnek olarak kullandığımız workflow'lar, ortamımızda bulunan UNIX araçlarını kullanarak çok temel metin işleme işlemlerini çalıştırmak zorunda kaldı.

Ancak, gerçek dünya pipeline'ları tipik olarak çoğu ortamda varsayılan olarak dahil edilmeyen özel araçlar ve paketler gerektirir.
Genellikle bu araçları yüklemeniz, bağımlılıklarını yönetmeniz ve herhangi bir çakışmayı çözmeniz gerekirdi.

Bunların hepsi çok sıkıcı ve can sıkıcıdır.
Bu sorunu ele almanın çok daha iyi bir yolu **konteyner** kullanmaktır.

Bir **konteyner**, kod, sistem kütüphaneleri ve ayarlar dahil bir uygulamayı çalıştırmak için gereken her şeyi içeren bir konteyner **imajından** oluşturulan hafif, bağımsız, çalıştırılabilir bir yazılım birimidir.

!!! Tip "İpucu"

    Bunu [Docker](https://www.docker.com/get-started/) teknolojisini kullanarak öğretiyoruz, ancak Nextflow birkaç başka konteyner teknolojisini de destekler.
    Nextflow'un konteynerler için desteği hakkında daha fazla bilgiyi [buradan](https://nextflow.io/docs/latest/container.html) öğrenebilirsiniz.

### 4.1. Bir konteyneri doğrudan kullanın

Önce bir konteynerle doğrudan etkileşime geçelim.
Bu, Nextflow'da kullanmaya başlamadan önce konteynerlerin ne olduğunu anlamanıza yardımcı olacaktır.

#### 4.1.1. Konteyner imajını çekin

Bir konteyneri kullanmak için, genellikle bir konteyner kaydından bir konteyner imajı indirirsiniz veya "çekersiniz" ve ardından bir konteyner örneği oluşturmak için konteyner imajını çalıştırırsınız.

Genel sözdizimi şöyledir:

```bash title="Sözdizimi"
docker pull '<container>'
```

- `docker pull`, konteyner sistemine bir havuzdan konteyner imajı çekmesi için talimat verir.
- `'<container>'`, konteyner imajının URI adresidir.

Örnek olarak, rastgele metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII art üreten `cowsay` adlı bir aracın python uygulaması olan [cowpy](https://github.com/jeffbuttars/cowpy) içeren bir konteyner imajı çekelim.

Yayınlanmış konteynerleri bulabileceğiniz çeşitli havuzlar vardır.
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

#### 4.1.2. Konteyneri başlatın

Konteynerler tek seferlik bir komut olarak çalıştırılabilir, ancak bunları etkileşimli olarak da kullanabilirsiniz, bu size konteynerin içinde bir shell promptu verir ve komutla oynamanıza olanak tanır.

Genel sözdizimi şöyledir:

```bash title="Sözdizimi"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'`, konteyner sistemine bir konteyner imajından bir konteyner örneği başlatması ve içinde bir komut yürütmesi için talimat verir.
- `--rm`, sisteme komut tamamlandıktan sonra konteyner örneğini kapatmasını söyler.

Tam olarak birleştirilmiş konteyner yürütme komutu şöyle görünür:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Bu komutu çalıştırın ve promptunuzun `(base) root@b645838b3314:/tmp#` gibi bir şeye dönüştüğünü görmelisiniz, bu artık konteynerin içinde olduğunuzu gösterir.

Bunu, dizin içeriğini listelemek için `ls` çalıştırarak doğrulayabilirsiniz:

```bash
ls /
```

??? success "Komut çıktısı"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Konteynerin içindeki dosya sisteminin ana sisteminizden farklı olduğunu görüyorsunuz.

!!! Tip "İpucu"

    Bir konteyneri çalıştırdığınızda, varsayılan olarak ana sistemden izole edilir.
    Bu, konteynerin ana sistemdeki hiçbir dosyaya erişemeyeceği anlamına gelir, bunu aşağıdaki sözdizimini kullanarak `docker run` komutunun bir parçası olarak bir birim bağlamak istediğinizi belirterek açıkça izin vermediğiniz sürece:

    ```bash title="Sözdizimi"
    -v <outside_path>:<inside_path>
    ```

    Bu, dosya sisteminizin o bölümüne erişmek için kullanabileceğiniz, konteyner duvarından bir tünel oluşturur.

    Bu, [Hello Nextflow Bölüm 5](../hello_nextflow/05_hello_containers.md)'te daha ayrıntılı olarak ele alınmaktadır.

#### 4.1.3. `cowpy` aracını çalıştırın

Konteynerin içinden `cowpy` komutunu doğrudan çalıştırabilirsiniz.

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

Bu, belirttiğimiz metni içeren bir konuşma balonuyla varsayılan inek karakterinin (veya 'cowacter') ASCII art'ını üretir.

Temel kullanımı test ettiğinize göre, ona bazı parametreler vermeyi deneyebilirsiniz.
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

Bu sefer ASCII art çıktısı Linux pengueni Tux'u gösteriyor, çünkü `-c tux` parametresini belirttik.

Konteynerin içinde olduğunuz için, sisteminize herhangi bir kütüphane yüklemek zorunda kalmadan girdi parametrelerini değiştirerek cowpy komutunu istediğiniz kadar çalıştırabilirsiniz.

??? tip "Diğer mevcut karakterler"

    Farklı bir karakter seçmek için '-c' bayrağını kullanın, bunlar dahil:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Bununla oynamaktan çekinmeyin.
İşiniz bittiğinde, `exit` komutunu kullanarak konteynerden çıkın:

```bash
exit
```

Normal shell'inize geri döneceksiniz.

### 4.2. Bir workflow'da konteyner kullanın

Bir pipeline çalıştırdığımızda, Nextflow'a her adımda hangi konteynerin kullanılacağını söyleyebilmek istiyoruz ve önemlisi, az önce yaptığımız tüm işi onun halletmesini istiyoruz: konteyneri çekmek, başlatmak, komutu çalıştırmak ve bittiğinde konteyneri kapatmak.

İyi haber: tam olarak Nextflow'un bizim için yapacağı şey bu.
Sadece her process için bir konteyner belirtmemiz gerekiyor.

Bunun nasıl çalıştığını göstermek için, üçüncü adımda üretilen toplanan selamlamalar dosyasında `cowpy` çalıştıran workflow'umuzun başka bir versiyonunu yaptık.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Bu, konuşma balonunda üç selamlama içeren ASCII art içeren bir dosya çıktısı vermelidir.

#### 4.2.1. Kodu inceleyin

Workflow, öncekine çok benzer, artı `cowpy` çalıştırmak için ekstra adım.

??? full-code "Tam kod dosyası"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
        // selamlamayı büyük harfe dönüştür
        convertToUpper(sayHello.out)
        // tüm selamlamaları tek bir dosyada topla
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

Bu workflow'un bir modül dosyasından bir `cowpy` process'i içe aktardığını ve bunu `collectGreetings()` çağrısının çıktısında, artı `params.character` adlı bir girdi parametresinde çağırdığını görüyorsunuz.

```groovy title="2d-container.nf" linenums="25"
// cowpy ile ASCII sanatı oluştur
cowpy(collectGreetings.out, params.character)
```

ASCII art oluşturmak için cowpy komutunu saran `cowpy` process'i, `cowpy.nf` modülünde tanımlanmıştır.

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

`cowpy` process'i iki girdi gerektirir: konuşma balonuna konulacak metni içeren bir girdi dosyasının yolu (`input_file`) ve karakter değişkeni için bir değer.

Önemlisi, daha önce kullandığımız konteyner URI'sine işaret eden `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'` satırını da içerir.

#### 4.2.2. Yapılandırmada Docker'ın etkinleştirildiğini kontrol edin

Bu eğitim kursunun 3. Bölümünü biraz öne çekerek, Nextflow'un workflow çalıştırmasını yapılandırmanın ana yollarından biri olan `nextflow.config` yapılandırma dosyasını tanıtacağız.
Geçerli dizinde `nextflow.config` adlı bir dosya bulunduğunda, Nextflow otomatik olarak yükleyecek ve içerdiği yapılandırmayı uygulayacaktır.

Bu amaçla, Docker'ı etkinleştiren tek satırlık kod içeren bir `nextflow.config` dosyası ekledik.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Bu yapılandırma, Nextflow'a uyumlu bir konteyner belirten herhangi bir process için Docker kullanmasını söyler.

!!! tip "İpucu"

    Docker çalıştırmasını komut satırından, çalıştırma başına bazında `-with-docker <container>` parametresini kullanarak etkinleştirmek teknik olarak mümkündür.
    Ancak, bu yalnızca tüm workflow için bir konteyner belirtmemize izin verirken, az önce gösterdiğimiz yaklaşım process başına farklı bir konteyner belirtmemize olanak tanır.
    İkincisi modülerlik, kod bakımı ve tekrar üretilebilirlik için çok daha iyidir.

#### 4.2.3. Workflow'u çalıştırın

Özetlemek gerekirse, çalıştırmak üzere olduğumuz şey:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Çalışacak mı sizce?

Workflow'u `-resume` bayrağıyla çalıştıralım ve karakterin hindi olmasını istediğimizi belirtelim.

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

İlk üç adım daha önce çalıştırdığımız için önbelleğe alındı, ancak `cowpy` process'i yeni olduğu için gerçekten çalıştırılıyor.

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

Karakterin tüm selamlamaları söylediğini görüyorsunuz, çünkü toplanan büyük harfli selamlamalar dosyasında çalıştı.

Daha da önemlisi, cowpy ve tüm bağımlılıklarının düzgün bir kurulumunu yapmak zorunda kalmadan bunu pipeline'ımızın bir parçası olarak çalıştırabildik.
Ve şimdi pipeline'ı iş arkadaşlarımızla paylaşabilir ve Docker veya yukarıda belirtildiği gibi alternatiflerinden biri (Singularity/Apptainer gibi) dışında hiçbir şey yüklemelerine gerek kalmadan kendi altyapılarında çalıştırmalarını sağlayabiliriz.

#### 4.2.4. Nextflow'un konteynerleştirilmiş görevi nasıl başlattığını inceleyin

Bu bölüme son bir koda olarak, Nextflow'un kaputun altında konteynerlerle nasıl çalıştığı hakkında biraz daha fazla bilgi edinmek için `cowpy` process çağrılarından birinin çalışma alt dizinine bir göz atalım.

`cowpy` process'inin çalışma alt dizininin yolunu bulmak için `nextflow run` komutunuzun çıktısını kontrol edin.
Yukarıda gösterilen çalıştırma için aldığımıza bakarsak, `cowpy` process'i için konsol log satırı `[7f/caf718]` ile başlıyor.
Bu, şu kısaltılmış dizin yoluna karşılık gelir: `work/7f/caf718`.

Bu dizinde, Nextflow'un pipeline'ı yürütme sürecinde sizin adınıza çalıştırdığı tüm komutları içeren `.command.run` dosyasını bulacaksınız.

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

Bu dosyada `nxf_launch` için arama yaparsanız, şuna benzer bir şey görmelisiniz:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Bu başlatma komutu, Nextflow'un process çağrısını başlatmak için manuel olarak yaptığımıza çok benzer bir `docker run` komutu kullandığını gösterir.
Ayrıca ilgili çalışma alt dizinini konteynere bağlar, konteyner içindeki çalışma dizinini buna göre ayarlar ve `.command.sh` dosyasındaki şablonlu bash betiğimizi çalıştırır.

Bu, önceki bölümde manuel olarak yapmak zorunda kaldığımız tüm zor işlerin şimdi Nextflow tarafından bizim için yapıldığını doğrular!

### Özet

Konteynerlerin araç sürümlerini yönetmede ve tekrar üretilebilirliği sağlamada ne gibi bir rol oynadığını anlıyorsunuz.

Daha genel olarak, gerçek dünya Nextflow pipeline'larının temel bileşenlerinin neler olduğu ve nasıl organize edildikleri hakkında temel bir anlayışa sahipsiniz.
Nextflow'un birden fazla girdiyi verimli bir şekilde nasıl işleyebileceğinin, birbirine bağlı birden fazla adımdan oluşan workflow'ları nasıl çalıştırabileceğinin, modüler kod bileşenlerinden nasıl yararlanabileceğinin ve daha fazla tekrar üretilebilirlik ve taşınabilirlik için konteynerleri nasıl kullanabileceğinin temellerini biliyorsunuz.

### Sırada ne var?

Bir mola daha verin! Bu, Nextflow pipeline'larının nasıl çalıştığına dair büyük bir bilgi yığınıydı.

Bu eğitimin son bölümünde, yapılandırma konusunu daha derinlemesine inceleyeceğiz.
Pipeline'ınızın yürütülmesini altyapınıza uyacak şekilde nasıl yapılandıracağınızı ve girdilerin ve parametrelerin yapılandırmasını nasıl yöneteceğinizi öğreneceksiniz.

---

## Kısa sınav

<quiz>
Nextflow neden her process çağrısı için ayrı bir görev dizini oluşturur?
- [ ] Yürütme hızını artırmak için
- [ ] Bellek kullanımını azaltmak için
- [x] Yürütmeleri izole etmek ve çıktılar arasındaki çakışmaları önlemek için
- [ ] Paralel dosya sıkıştırmayı etkinleştirmek için

Daha fazla bilgi: [1.3. Orijinal çıktıları ve logları bulun](#13-orijinal-çıktıları-ve-logları-bulun)
</quiz>

<quiz>
Bir workflow çalıştırılırken `-ansi-log false` seçeneği ne yapar?
- [ ] Tüm konsol çıktısını devre dışı bırakır
- [x] Çıktıdan rengi kaldırır
- [x] Tüm görev dizini yollarını tek bir satırda yoğunlaştırmak yerine gösterir
- [ ] Ayrıntılı hata ayıklama modunu etkinleştirir

Daha fazla bilgi: [1.3.2. Terminalin daha fazla ayrıntı göstermesini sağlayın](#132-terminalin-daha-fazla-ayrıntı-göstermesini-sağlayın)

Alternatif olarak, bu stili tercih ediyorsanız aşağıdaki ortam değişkenlerinden birini de kullanabilirsiniz:

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
- [ ] Resume işlevselliğini etkinleştirmek için
- [x] Birden fazla girdi işlenirken çıktı dosyalarının birbirinin üzerine yazılmasını önlemek için
- [ ] Dosyaların sıkıştırılmasını kolaylaştırmak için

Daha fazla bilgi: [1.4.3. Çıktılar nasıl adlandırılır](#143-çıktılar-nasıl-adlandırılır)
</quiz>

<quiz>
Modülerleştirilmiş bir workflow'da `include` ifadesinin amacı nedir?
- [ ] Process kodunu workflow dosyasına kopyalamak için
- [x] Harici bir modül dosyasından bir process tanımını içe aktarmak için
- [ ] Yapılandırma ayarlarını dahil etmek için
- [ ] Belge açıklamaları eklemek için

Daha fazla bilgi: [3. Modülerleştirilmiş pipeline'ları çalıştırma](#3-modülerleştirilmiş-pipelineları-çalıştırma)
</quiz>

<quiz>
Bir workflow'u modülerleştirdiğinizde ve `-resume` ile çalıştırdığınızda ne olur?
- [ ] Modüler process'ler için önbellekleme devre dışı bırakılır
- [ ] Tüm görevler yeniden yürütülmelidir
- [x] Önbellekleme, oluşturulan iş betiklerine göre normal şekilde çalışır
- [ ] Yalnızca ana workflow dosyası önbelleğe alınır

Daha fazla bilgi: [3.2. Workflow'u çalıştırın](#32-workflowu-çalıştırın)
</quiz>

<quiz>
Bir process tanımındaki `container` yönergesi neyi belirtir?
- [ ] Process için çalışma dizini
- [ ] Maksimum bellek tahsisi
- [x] Process'i çalıştırmak için kullanılacak konteyner imajı URI'si
- [ ] Çıktı dosya formatı

Daha fazla bilgi: [4.2. Bir workflow'da konteyner kullanın](#42-bir-workflowda-konteyner-kullanın)
</quiz>

<quiz>
`.command.run` dosyasında, `nxf_launch` fonksiyonu ne içerir?
- [ ] Nextflow sürüm bilgisi
- [ ] Workflow parametreleri
- [x] Birim bağlantıları ve konteyner ayarları içeren `docker run` komutu
- [ ] Process girdi bildirimleri

Daha fazla bilgi: [4.2.4. Nextflow'un konteynerleştirilmiş görevi nasıl başlattığını inceleyin](#424-nextflowun-konteynerleştirilmiş-görevi-nasıl-başlattığını-inceleyin)
</quiz>

<quiz>
Nextflow, konteynerleştirilmiş bir process çalıştırırken otomatik olarak ne halleder? (Tümünü seçin)
- [x] Gerekirse konteyner imajını çekme
- [x] Çalışma dizinini konteynere bağlama
- [x] Process betiğini konteyner içinde çalıştırma
- [x] Yürütme sonrasında konteyner örneğini temizleme

Daha fazla bilgi: [4. Konteynerleştirilmiş yazılım kullanma](#4-konteynerleştirilmiş-yazılım-kullanma)
</quiz>
