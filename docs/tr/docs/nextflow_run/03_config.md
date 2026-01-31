# Bölüm 3: Çalıştırma yapılandırması

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu bölümde, bir Nextflow pipeline'ının yapılandırmasını yönetmeyi keşfedeceğiz; davranışını özelleştirmek, farklı ortamlara uyarlamak ve kaynak kullanımını optimize etmek için _workflow kodunun tek bir satırını değiştirmeden_.

Bunu yapmanın birden fazla yolu vardır; bunlar birlikte kullanılabilir ve [burada](https://www.nextflow.io/docs/latest/config.html) açıklanan öncelik sırasına göre yorumlanır.

Kursun bu bölümünde, Bölüm 2'de konteynerler bölümünde zaten karşılaştığınız `nextflow.config` dosyası olan en basit ve en yaygın yapılandırma dosyası mekanizmasını göstereceğiz.

Process direktifleri, executor'lar, profiller ve parametre dosyaları gibi Nextflow yapılandırmasının temel bileşenlerini inceleyeceğiz.
Bu yapılandırma seçeneklerini etkili bir şekilde kullanmayı öğrenerek, Nextflow pipeline'larının esnekliğinden, ölçeklenebilirliğinden ve performansından tam olarak yararlanabilirsiniz.

Bu yapılandırma öğelerini uygulamak için, bu eğitim kursunun 2. Bölümünün sonunda çalıştırdığımız workflow'un `3-main.nf` olarak yeniden adlandırılmış yeni bir kopyasını çalıştıracağız.

Hello pipeline'ına aşina değilseniz veya bir hatırlatmaya ihtiyacınız varsa, [bu bilgi sayfasına](../info/hello_pipeline.md) bakın.

---

## 1. Workflow girdi parametrelerini yönetme

??? example "Senaryo"

    Bir pipeline indirdiniz ve aynı girdi dosyaları ve ayarlarla tekrar tekrar çalıştırmak istiyorsunuz, ancak her seferinde tüm parametreleri yazmak istemiyorsunuz.
    Ya da belki komut satırı argümanlarıyla rahat olmayan bir meslektaşınız için pipeline'ı kuruyorsunuz.

Şu ana kadar üzerinde çalıştığımız şeyin bir uzantısı olan bir yapılandırma yönüyle başlayacağız: girdi parametrelerinin yönetimi.

Şu anda workflow'umuz, workflow betiğinde bir `params` bloğunda bildirilen komut satırı aracılığıyla birkaç parametre değeri kabul edecek şekilde ayarlanmıştır.
Bunlardan birinin bildiriminin bir parçası olarak varsayılan değeri ayarlanmıştır.

Ancak, hepsi için varsayılan değerler ayarlamak veya mevcut varsayılanı komut satırında parametre belirtmek ya da orijinal betik dosyasını değiştirmek zorunda kalmadan geçersiz kılmak isteyebilirsiniz.

Bunu yapmanın birden fazla yolu vardır; size çok yaygın olarak kullanılan üç temel yolu göstereceğiz.

### 1.1. `nextflow.config`'de değerler ayarlayın

Bu en basit yaklaşımdır, ancak muhtemelen en az esnektir çünkü ana `nextflow.config` dosyası her çalıştırma için düzenlemek isteyeceğiniz bir şey değildir.
Ancak, workflow'da parametreleri _bildirme_ (ki kesinlikle oraya aittir) ile bir yapılandırma dosyasında daha evde olan _varsayılan değerler_ sağlama endişelerini ayırma avantajına sahiptir.

Bunu iki adımda yapalım.

#### 1.1.1. Yapılandırma dosyasında bir `params` bloğu oluşturun

`nextflow.config` dosyasında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parametreleri
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

`params` bloğunu workflow'dan yapılandırma dosyasına sadece kopyalamadığımıza dikkat edin.
Zaten varsayılan değeri bildirilen `batch` parametresi için sözdizimi biraz farklıdır.
Workflow dosyasında bu türlenmiş bir bildirimdir.
Yapılandırmada bunlar değer atamalarıdır.

Teknik olarak, bu workflow dosyasında hala belirtilen varsayılan değerleri geçersiz kılmak için yeterlidir.
`batch` için varsayılan değeri değiştirebilir ve yapılandırma dosyasında ayarlanan değerin workflow dosyasında ayarlanan değeri geçersiz kıldığını kendiniz doğrulamak için workflow'u çalıştırabilirsiniz.

Ancak yapılandırmayı tamamen yapılandırma dosyasına taşıma ruhuyla, bu varsayılan değeri workflow dosyasından tamamen kaldıralım.

#### 1.1.2. Workflow dosyasındaki `batch` için varsayılan değeri kaldırın

`3-main.nf` workflow dosyasında aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parametreleri
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Önce"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parametreleri
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Artık workflow dosyasının kendisi bu parametreler için herhangi bir varsayılan değer ayarlamıyor.

#### 1.1.3. Pipeline'ı çalıştırın

Komut satırında herhangi bir parametre belirtmeden doğru çalışıp çalışmadığını test edelim.

```bash
nextflow run 3-main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretiyor.

Nihai ASCII art çıktısı `results/3-main/` dizininde, öncekiyle aynı şekilde `cowpy-COLLECTED-batch-output.txt` adı altında.

??? abstract "Dosya içeriği"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
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

İşlevsel olarak bu taşıma hiçbir şeyi değiştirmedi, ancak kavramsal olarak varsayılan değerlerin yapılandırma dosyasında ayarlanması biraz daha temizdir.

### 1.2. Çalıştırmaya özel yapılandırma dosyası kullanın

??? example "Senaryo"

    Ana yapılandırma dosyanızı değiştirmeden farklı ayarlarla denemeler yapmak istiyorsunuz.

Bunu, deneyler için çalışma dizini olarak kullanacağınız bir alt dizinde yeni bir `nextflow.config` dosyası oluşturarak yapabilirsiniz.

#### 1.2.1. Boş bir yapılandırmayla çalışma dizini oluşturun

Yeni bir dizin oluşturup içine girerek başlayalım:

```bash
mkdir -p tux-run
cd tux-run
```

Ardından, bu dizinde boş bir yapılandırma dosyası oluşturun:

```bash
touch nextflow.config
```

Bu boş bir dosya üretir.

#### 1.2.2. Deneysel yapılandırmayı kurun

Şimdi yeni dosyayı açın ve özelleştirmek istediğiniz parametreleri ekleyin:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Girdi dosyasının yolunun dizin yapısını yansıtması gerektiğini unutmayın.

#### 1.2.3. Pipeline'ı çalıştırın

Artık pipeline'ımızı yeni çalışma dizinimizin içinden çalıştırabiliriz.
Yolu buna göre uyarladığınızdan emin olun!

```bash
nextflow run ../3-main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Bu, `tux-run/` altında `tux-run/work/` ve `tux-run/results/` dahil yeni dizinler oluşturacaktır.

Bu çalıştırmada, Nextflow mevcut dizinimizdeki `nextflow.config`'i pipeline'ın kök dizinindeki `nextflow.config` ile birleştirir ve böylece varsayılan karakteri (turkey) tux karakteriyle geçersiz kılar.

Nihai çıktı dosyası, selamlamaları söyleyen tux karakterini içermelidir.

??? abstract "Dosya içeriği"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

İşte bu kadar; artık 'normal' yapılandırmanızı değiştirmeden denemeler yapmak için bir alanınız var.

!!! warning "Uyarı"

    Bir sonraki bölüme geçmeden önce önceki dizine geri dönmeyi unutmayın!

    ```bash
    cd ..
    ```

Şimdi parametre değerlerini ayarlamanın başka bir yararlı yoluna bakalım.

### 1.3. Parametre dosyası kullanın

??? example "Senaryo"

    Tam çalıştırma parametrelerini bir iş arkadaşınızla paylaşmanız veya bir yayın için kaydetmeniz gerekiyor.

Alt dizin yaklaşımı denemeler için harika çalışır, ancak biraz kurulum içerir ve yolları buna göre uyarlamanızı gerektirir.
Pipeline'ınızı belirli bir değer setiyle çalıştırmak veya başka birinin bunu minimum çabayla yapmasını sağlamak istediğinizde daha basit bir yaklaşım vardır.

Nextflow, parametreleri YAML veya JSON formatında bir parametre dosyası aracılığıyla belirtmemize olanak tanır; bu, örneğin alternatif varsayılan değer setlerini ve çalıştırmaya özel parametre değerlerini yönetmeyi ve dağıtmayı çok uygun hale getirir.

#### 1.3.1. Örnek parametre dosyasını inceleyin

Bunu göstermek için, geçerli dizinde `test-params.yaml` adlı örnek bir parametre dosyası sağlıyoruz:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Bu parametre dosyası, belirtmek istediğimiz girdilerin her biri için bir anahtar-değer çifti içerir.
Sözdizimini yapılandırma dosyasıyla karşılaştırırsanız eşittir işaretleri (`=`) yerine iki nokta üst üste (`:`) kullanımına dikkat edin.
Yapılandırma dosyası Groovy'de yazılırken, parametre dosyası YAML'de yazılır.

!!! info "Bilgi"

    Ayrıca örnek olarak parametre dosyasının JSON versiyonunu da sağlıyoruz ancak burada onunla çalıştırmayacağız.
    Onu kendi başınıza denemeye çekinmeyin.

#### 1.3.2. Pipeline'ı çalıştırın

Workflow'u bu parametre dosyasıyla çalıştırmak için, temel komuta `-params-file <filename>` eklemeniz yeterlidir.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Nihai çıktı dosyası, selamlamaları söyleyen stegosaurus karakterini içermelidir.

??? abstract "Dosya içeriği"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Yalnızca birkaç parametreniz olduğunda parametre dosyası kullanmak aşırı görünebilir, ancak bazı pipeline'lar düzinelerce parametre bekler.
Bu durumlarda, bir parametre dosyası kullanmak, devasa komut satırları yazmak zorunda kalmadan ve workflow betiğini değiştirmeden çalışma zamanında parametre değerleri sağlamamıza olanak tanır.

Ayrıca parametre setlerini iş arkadaşlarına veya örneğin bir yayın için destekleyici bilgi olarak dağıtmayı kolaylaştırır.
Bu, çalışmanızı başkaları tarafından daha tekrar üretilebilir hale getirir.

### Özet

Workflow girdilerini yönetmek için temel yapılandırma seçeneklerinden nasıl yararlanacağınızı biliyorsunuz.

### Sırada ne var?

Workflow çıktılarınızın nerede ve nasıl yayınlandığını yönetmeyi öğrenin.

---

## 2. Workflow çıktılarını yönetme

??? example "Senaryo"

    Pipeline'ınız çıktıları sabit kodlanmış bir dizine yayınlıyor, ancak her seferinde workflow kodunu düzenlemeden sonuçları proje veya deney adına göre organize etmek istiyorsunuz.

Miras aldığımız workflow, workflow düzeyinde çıktı bildirimleri için yollar kullanıyor, bu çok esnek değil ve çok tekrar içeriyor.

Bunu daha esnek yapılandırmak için birkaç yaygın yola bakalım.

### 2.1. `outputDir` dizin adını özelleştirme

Şu ana kadar çalıştırdığımız workflow'un her versiyonu çıktılarını çıktı tanımlarına sabit kodlanmış farklı bir alt dizine yayınladı.

Bunu kullanıcı tarafından yapılandırılabilir bir parametre kullanacak şekilde değiştirelim.
Bunun için tamamen yeni bir parametre oluşturabilirdik, ancak `batch` parametresi tam orada olduğu için onu kullanalım.

#### 2.1.1. Yapılandırma dosyasında `outputDir` için bir değer ayarlayın

Nextflow'un çıktıları yayınlamak için kullandığı yol `outputDir` seçeneğiyle kontrol edilir.
Tüm çıktılar için yolu değiştirmek için bu seçenek için `nextflow.config` yapılandırma dosyasında bir değer ayarlayabilirsiniz.

`nextflow.config` dosyasına aşağıdaki kodu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parametreleri
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parametreleri
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Bu, yerleşik varsayılan yol olan `results/`'ı `results/` artı alt dizin olarak `batch` parametresinin değeriyle değiştirecektir.
İsterseniz `results` kısmını da değiştirebilirsiniz.

Geçici bir değişiklik için, komutunuzda `-output-dir` parametresini kullanarak komut satırından bu seçeneği ayarlayabilirsiniz (ancak bu durumda `batch` parametre değerini kullanamazsınız).

#### 2.1.2. Sabit kodlanmış yolun tekrarlanan kısmını kaldırın

Çıktı seçeneklerinde hala sabit kodlanmış bir alt dizinimiz var, bu yüzden şimdi onu kaldıralım.

Workflow dosyasında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

`outputDir` varsayılanını değiştirmek yerine her yola `${params.batch}` eklemiş de olabilirdik, ancak bu daha özlüdür.

#### 2.1.3. Pipeline'ı çalıştırın

Batch adını komut satırından `outdir` olarak ayarlayarak doğru çalışıp çalışmadığını test edelim.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Bu hala öncekiyle aynı çıktıyı üretiyor, ancak bu sefer çıktılarımızı `results/outdir/` altında buluyoruz.

??? abstract "Dizin içeriği"

    ```console
    results/outdir/
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Bu yaklaşımı özel yol tanımlarıyla birleştirerek istediğiniz herhangi bir dizin hiyerarşisini oluşturabilirsiniz.

### 2.2. Çıktıları process'e göre organize etme

Çıktıları daha fazla organize etmenin popüler bir yolu, process'e göre yapmaktır, yani pipeline'da çalıştırılan her process için alt dizinler oluşturmak.

#### 2.2.1. Çıktı yollarını process adlarına referansla değiştirin

Tek yapmanız gereken, çıktı yolu bildiriminde process adına `<task>.name` olarak referans vermektir.

Workflow dosyasında aşağıdaki değişiklikleri yapın:

=== "Sonra"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Bu, çıktı yolu yapılandırmasından kalan sabit kodlanmış öğeleri kaldırır.

#### 2.2.2. Pipeline'ı çalıştırın

Batch adını komut satırından `pnames` olarak ayarlayarak doğru çalışıp çalışmadığını test edelim.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Bu hala öncekiyle aynı çıktıyı üretiyor, ancak bu sefer çıktılarımızı `results/pnames/` altında buluyoruz ve process'e göre gruplandırılmışlar.

??? abstract "Dizin içeriği"

    ```console
    results/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Burada `intermediates` ile nihai çıktıların üst düzeyde olması arasındaki ayrımı sildiğimizi unutmayın.
Elbette bu yaklaşımları karıştırıp eşleştirebilirsiniz, örneğin ilk çıktının yolunu `intermediates/${sayHello.name}` olarak ayarlayarak.

### 2.3. Workflow düzeyinde yayınlama modunu ayarlama

Son olarak, tekrarlayan kod miktarını azaltma ruhuyla, çıktı başına `mode` bildirimlerini yapılandırmada tek bir satırla değiştirebiliriz.

#### 2.3.1. Yapılandırma dosyasına `workflow.output.mode` ekleyin

`nextflow.config` dosyasına aşağıdaki kodu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

`outputDir` seçeneği gibi, yapılandırma dosyasında `workflow.output.mode`'a bir değer vermek, workflow dosyasında ayarlanmış olanı geçersiz kılmak için yeterli olurdu, ancak yine de gereksiz kodu kaldıralım.

#### 2.3.2. Workflow dosyasından çıktı modunu kaldırın

Workflow dosyasında aşağıdaki değişiklikleri yapın:

=== "Sonra"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Önce"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

Bu daha özlü, değil mi?

#### 2.3.3. Pipeline'ı çalıştırın

Batch adını komut satırından `outmode` olarak ayarlayarak doğru çalışıp çalışmadığını test edelim.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Bu hala öncekiyle aynı çıktıyı üretiyor, ancak bu sefer çıktılarımızı `results/outmode/` altında buluyoruz.
Hala hepsi düzgün kopyalar, symlink değil.

??? abstract "Dizin içeriği"

    ```console
    results/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Çıktı başına mod ayarlama yolunu hala kullanmak istemenizin ana nedeni, aynı workflow içinde karıştırıp eşleştirmek istemeniz durumudur, yani bazı çıktıların kopyalanmasını ve bazılarının symlink olmasını isterseniz.

Bu şekilde özelleştirebileceğiniz birçok başka seçenek var, ancak umarız bu size seçeneklerin kapsamı ve bunları tercihlerinize uygun şekilde etkili bir şekilde nasıl kullanacağınız hakkında bir fikir verir.

### Özet

Çıktılarınızın yayınlandığı dizinlerin adlandırmasını ve yapısını ve workflow çıktı yayınlama modunu nasıl kontrol edeceğinizi biliyorsunuz.

### Sırada ne var?

Yazılım paketleme teknolojisiyle başlayarak workflow yapılandırmanızı hesaplama ortamınıza nasıl uyarlayacağınızı öğrenin.

---

## 3. Yazılım paketleme teknolojisi seçme

Şu ana kadar girdilerin nasıl girdiğini ve çıktıların nereden çıktığını kontrol eden yapılandırma öğelerine baktık. Şimdi workflow yapılandırmanızı hesaplama ortamınıza uyarlamaya daha özel olarak odaklanma zamanı.

Bu yoldaki ilk adım, her adımda çalıştırılacak yazılım paketlerinin nereden geleceğini belirlemektir.
Zaten yerel hesaplama ortamında yüklü mü?
İmajları alıp bir konteyner sistemi aracılığıyla çalıştırmamız mı gerekiyor?
Yoksa Conda paketlerini alıp yerel bir Conda ortamı mı oluşturmamız gerekiyor?

Bu eğitim kursunun ilk bölümünde (Bölümler 1-4) workflow'umuzda yalnızca yerel olarak yüklenmiş yazılımları kullandık.
Ardından Bölüm 5'te Docker konteynerlarını ve Docker konteynerların kullanımını etkinleştirmek için kullandığımız `nextflow.config` dosyasını tanıttık.

Şimdi `nextflow.config` dosyası aracılığıyla alternatif bir yazılım paketleme seçeneğini nasıl yapılandırabileceğimizi görelim.

### 3.1. Docker'ı devre dışı bırakın ve yapılandırma dosyasında Conda'yı etkinleştirin

??? example "Senaryo"

    Pipeline'ınızı güvenlik nedenleriyle Docker'ın izin verilmediği bir HPC kümesine taşıyorsunuz.
    Küme Singularity ve Conda'yı destekliyor, bu yüzden yapılandırmanızı buna göre değiştirmeniz gerekiyor.

Nextflow, HPC'de daha yaygın olarak kullanılan Singularity ve Conda gibi yazılım paket yöneticileri dahil birden fazla konteyner teknolojisini destekler.

Yapılandırma dosyamızı Docker yerine Conda kullanacak şekilde değiştirebiliriz.
Bunu yapmak için `docker.enabled` değerini `false` olarak değiştirelim ve Conda kullanımını etkinleştiren bir direktif ekleyelim:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Bu, Nextflow'un Conda paketleri belirtilmiş process'ler için Conda ortamları oluşturmasına ve kullanmasına izin verecektir.
Bu da artık `cowpy` process'imize bunlardan birini eklememiz gerektiği anlamına gelir!

### 3.2. Process tanımında bir Conda paketi belirtin

`cowpy` aracını içeren bir Conda paketi için URI'yi zaten aldık: `conda-forge::cowpy==1.1.5`

Şimdi `conda` direktifini kullanarak URI'yi `cowpy` process tanımına ekliyoruz:

=== "Sonra"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Önce"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Açık olmak gerekirse, `docker` direktifini _değiştirmiyoruz_, alternatif bir seçenek _ekliyoruz_.

!!! tip "İpucu"

    Belirli bir conda paketi için URI almanın birkaç farklı yolu vardır.
    [Seqera Containers](https://seqera.io/containers/) arama sorgusunu kullanmanızı öneririz; bu size kopyalayıp yapıştırabileceğiniz bir URI verecektir, ondan bir konteyner oluşturmayı planlamasanız bile.

### 3.3. Conda'yı kullanabildiğini doğrulamak için workflow'u çalıştırın

Deneyelim.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Komut çıktısı"

    ```console title="Çıktı"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Bu sorunsuz çalışmalı ve `results/conda` altında öncekiyle aynı çıktıları üretmelidir.

Arka planda, Nextflow Conda paketlerini aldı ve ortamı oluşturdu, bu normalde biraz iş gerektirir; bu yüzden bunların hiçbirini kendimiz yapmak zorunda kalmamamız güzel!

!!! info "Bilgi"

    Bu hızlı çalışır çünkü `cowpy` paketi oldukça küçüktür, ancak büyük paketlerle çalışıyorsanız, ilk seferinde normalden biraz daha uzun sürebilir ve tamamlanmadan önce konsol çıktısının bir dakika kadar 'takılı' kaldığını görebilirsiniz.
    Bu normaldir ve Nextflow'un yeni bir paket kullandığınızda ilk seferinde yaptığı ekstra işten kaynaklanır.

Bizim bakış açımızdan, arka planda mekanikler biraz farklı olsa bile, Docker ile çalıştırmakla tamamen aynı görünüyor.

Bu, gerekirse Conda ortamlarıyla çalıştırmaya hazır olduğumuz anlamına gelir.

??? info "Docker ve Conda'yı karıştırma ve eşleştirme"

    Bu direktifler process başına atandığından, 'karıştırıp eşleştirmek' mümkündür, yani örneğin kullandığınız hesaplama altyapısı her ikisini de destekliyorsa, workflow'unuzdaki bazı process'leri Docker ile ve diğerlerini Conda ile çalışacak şekilde yapılandırmak mümkündür.
    Bu durumda, yapılandırma dosyanızda hem Docker'ı hem de Conda'yı etkinleştirirsiniz.
    Belirli bir process için her ikisi de mevcutsa, Nextflow konteynerleri önceliklendirecektir.

    Ve daha önce belirtildiği gibi, Nextflow birden fazla başka yazılım paketleme ve konteyner teknolojisini destekler, bu yüzden sadece bu ikisiyle sınırlı değilsiniz.

### Özet

Her process'in hangi yazılım paketini kullanması gerektiğini nasıl yapılandıracağınızı ve teknolojiler arasında nasıl geçiş yapacağınızı biliyorsunuz.

### Sırada ne var?

Nextflow'un işi gerçekten yapmak için kullandığı çalıştırma platformunu nasıl değiştireceğinizi öğrenin.

---

## 4. Çalıştırma platformu seçme

??? example "Senaryo"

    Pipeline'ınızı dizüstü bilgisayarınızda geliştiriyordunuz ve test ediyordunuz, ancak şimdi binlerce örnek üzerinde çalıştırmanız gerekiyor.
    Kurumunuzun bunun yerine kullanmak istediğiniz Slurm scheduler'lı bir HPC kümesi var.

Şu ana kadar pipeline'ımızı yerel executor ile çalıştırıyorduk.
Bu, her görevi Nextflow'un çalıştığı makinede yürütür.
Nextflow başladığında, mevcut CPU'lara ve belleğe bakar.
Çalıştırılmaya hazır görevlerin kaynakları mevcut kaynakları aşarsa, Nextflow önceki görevlerden bir veya daha fazlası tamamlanıp gerekli kaynakları serbest bırakana kadar son görevleri çalıştırmaktan geri tutar.

Yerel executor uygun ve verimlidir, ancak tek bir makineyle sınırlıdır. Çok büyük iş yükleri için, yerel makinenizin bir darboğaz olduğunu keşfedebilirsiniz; ya mevcut kaynaklardan daha fazlasını gerektiren tek bir göreviniz var ya da tek bir makinenin onları çalıştırmasını beklemenin çok uzun süreceği kadar çok göreviniz var.

Nextflow, HPC scheduler'ları (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor ve diğerleri) ve bulut çalıştırma backend'leri (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes ve daha fazlası) dahil [birçok farklı çalıştırma backend'ini](https://www.nextflow.io/docs/latest/executor.html) destekler.

### 4.1. Farklı bir backend'i hedefleme

Executor seçimi `executor` adlı bir process direktifi ile ayarlanır.
Varsayılan olarak `local` olarak ayarlanmıştır, bu yüzden aşağıdaki yapılandırma ima edilir:

```groovy title="Yerleşik yapılandırma"
process {
    executor = 'local'
}
```

Executor'ı farklı bir backend'i hedefleyecek şekilde ayarlamak için, kaynak tahsisleri için yukarıda açıklanan benzer sözdizimini kullanarak istediğiniz executor'ı belirtmeniz yeterlidir (tüm seçenekler için [belgelere](https://www.nextflow.io/docs/latest/executor.html) bakın).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Uyarı"

    Bunu eğitim ortamında gerçekten test edemiyoruz çünkü bir HPC'ye bağlanacak şekilde ayarlanmamış.

### 4.2. Çalıştırma parametreleri için backend'e özgü sözdizimi ile başa çıkma

Çoğu yüksek performanslı hesaplama platformu, kaynak tahsis istekleri ve sınırlamaları (örn. CPU sayısı ve bellek) ve kullanılacak iş kuyruğunun adı gibi belirli parametreleri belirtmenize izin verir (ve bazen gerektirir).

Ne yazık ki, bu sistemlerin her biri, bir işin nasıl tanımlanması ve ilgili scheduler'a nasıl gönderilmesi gerektiğini tanımlamak için farklı teknolojiler, sözdizimler ve yapılandırmalar kullanır.

??? abstract "Örnekler"

    Örneğin, 8 CPU ve 4GB RAM gerektiren ve "my-science-work" kuyruğunda çalıştırılacak aynı iş, backend'e bağlı olarak farklı şekillerde ifade edilmelidir.

    ```bash title="SLURM için yapılandırma / sbatch ile gönder"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="PBS için yapılandırma / qsub ile gönder"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="SGE için yapılandırma / qsub ile gönder"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Neyse ki, Nextflow tüm bunları basitleştirir.
İlgili özellikleri, `cpus`, `memory` ve `queue` gibi (diğer özellikler için belgelere bakın) yalnızca bir kez belirtebilmeniz için standartlaştırılmış bir sözdizimi sağlar.
Ardından, çalışma zamanında, Nextflow executor ayarına göre uygun backend'e özgü betikleri oluşturmak için bu ayarları kullanacaktır.

Bu standartlaştırılmış sözdizimini bir sonraki bölümde ele alacağız.

### Özet

Artık farklı türde hesaplama altyapısı kullanmak için executor'ı nasıl değiştireceğinizi biliyorsunuz.

### Sırada ne var?

Nextflow'da kaynak tahsislerini ve sınırlamalarını değerlendirmeyi ve ifade etmeyi öğrenin.

---

## 5. Hesaplama kaynak tahsislerini kontrol etme

??? example "Senaryo"

    Pipeline'ınız kümede sürekli başarısız oluyor çünkü görevler bellek sınırlarını aştığı için sonlandırılıyor.
    Ya da belki kullanmadığınız kaynaklar için ücretlendiriliyorsunuz ve maliyetleri optimize etmek istiyorsunuz.

Çoğu yüksek performanslı hesaplama platformu, CPU sayısı ve bellek gibi belirli kaynak tahsis parametrelerini belirtmenize izin verir (ve bazen gerektirir).

Varsayılan olarak, Nextflow her process için tek bir CPU ve 2GB bellek kullanacaktır.
İlgili process direktifleri `cpus` ve `memory` olarak adlandırılır, bu yüzden aşağıdaki yapılandırma ima edilir:

```groovy title="Yerleşik yapılandırma" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Bu değerleri, yapılandırma dosyanızda ek process direktifleri kullanarak tüm process'ler veya belirli adlandırılmış process'ler için değiştirebilirsiniz.
Nextflow bunları seçilen executor için uygun talimatlara çevirecektir.

Peki hangi değerleri kullanacağınızı nasıl biliyorsunuz?

### 5.1. Kaynak kullanım raporu oluşturmak için workflow'u çalıştırın

??? example "Senaryo"

    Process'lerinizin ne kadar bellek veya CPU'ya ihtiyaç duyduğunu bilmiyorsunuz ve kaynakları israf etmekten veya işlerin sonlandırılmasından kaçınmak istiyorsunuz.

Process'lerinizin ne kadar CPU ve belleğe ihtiyaç duyacağını önceden bilmiyorsanız, biraz kaynak profilleme yapabilirsiniz; yani workflow'u bazı varsayılan tahsislerle çalıştırırsınız, her process'in ne kadar kullandığını kaydedersiniz ve oradan temel tahsisleri nasıl ayarlayacağınızı tahmin edersiniz.

Uygun bir şekilde, Nextflow bunu yapmak için yerleşik araçlar içerir ve istek üzerine sizin için bir rapor oluşturmaktan mutluluk duyar.

Bunu yapmak için, komut satırınıza `-with-report <filename>.html` ekleyin.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Rapor, tarayıcınızda indirebileceğiniz ve açabileceğiniz bir html dosyasıdır. Ayrıca eğitim ortamında görüntülemek için soldaki dosya gezgininde sağ tıklayıp `Show preview`'a tıklayabilirsiniz.

Kaynakları ayarlama fırsatlarını tanımlayıp tanımlayamayacağınızı görmek için rapora bakıp birkaç dakikanızı ayırın.
Kullanım sonuçlarını tahsis edilenin yüzdesi olarak gösteren sekmelere tıkladığınızdan emin olun.
Mevcut tüm özellikleri açıklayan bazı [belgeler](https://www.nextflow.io/docs/latest/reports.html) var.

### 5.2. Tüm process'ler için kaynak tahsislerini ayarlama

Profilleme, eğitim workflow'umuzdaki process'lerin çok hafif olduğunu gösteriyor, bu yüzden process başına varsayılan bellek tahsisini 1GB'a düşürelim.

`nextflow.config` dosyanıza, pipeline parametreleri bölümünden önce aşağıdakini ekleyin:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

Bu, tükettiğimiz hesaplama miktarını azaltmaya yardımcı olacaktır.

### 5.3. Belirli bir process için kaynak tahsislerini ayarlama

Aynı zamanda, tek bir process için tahsislerin nasıl ayarlanacağını gösterebilmemiz için `cowpy` process'inin diğerlerinden daha fazla kaynak gerektirdiğini varsayacağız.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Bu yapılandırmayla, `cowpy` process'i hariç tüm process'ler 1GB bellek ve tek bir CPU (ima edilen varsayılan) isteyecektir; `cowpy` process'i 2GB ve 2 CPU isteyecektir.

!!! info "Bilgi"

    Birkaç CPU'ya sahip bir makineniz varsa ve process başına yüksek sayıda tahsis ederseniz, process çağrılarının birbiri ardına kuyruğa alındığını görebilirsiniz.
    Bunun nedeni, Nextflow'un mevcut olandan daha fazla CPU istememizi sağlamasıdır.

### 5.4. Güncellenmiş yapılandırmayla workflow'u çalıştırın

Bunu deneyelim, yapılandırma değişikliklerinden önce ve sonra performansı karşılaştırabilmemiz için profilleme raporu için farklı bir dosya adı sağlayalım.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Bu kadar küçük bir iş yükü olduğu için muhtemelen gerçek bir fark fark etmeyeceksiniz, ancak gerçek dünya workflow'unun performansını ve kaynak gereksinimlerini analiz etmek için kullanacağınız yaklaşım budur.

Process'leriniz farklı kaynak gereksinimlerine sahip olduğunda çok yararlıdır. Tahmin değil, gerçek verilere dayalı olarak her process için ayarladığınız kaynak tahsislerini doğru boyutlandırmanızı sağlar.

!!! tip "İpucu"

    Bu, kaynakları optimize etmek için yapabileceğiniz şeylerin sadece küçük bir tadımlığı.
    Nextflow'un kendisi, kaynak sınırlamaları nedeniyle başarısız olan işleri yeniden denemek için gerçekten şık bir [dinamik yeniden deneme mantığı](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) yerleşik olarak içerir.
    Ayrıca, Seqera Platform kaynak tahsislerinizi otomatik olarak optimize etmek için AI odaklı araçlar da sunmaktadır.

### 5.5. Kaynak sınırları ekleme

Hangi hesaplama executor'ını ve hesaplama altyapısını kullandığınıza bağlı olarak, tahsis edebileceğiniz (veya etmeniz gereken) konusunda bazı kısıtlamalar olabilir.
Örneğin, kümeniz belirli sınırlar içinde kalmanızı gerektirebilir.

İlgili sınırlamaları ayarlamak için `resourceLimits` direktifini kullanabilirsiniz. Sözdizimi, tek başına bir process bloğunda olduğunda şöyle görünür:

```groovy title="Sözdizimi örneği"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow bu değerleri belirttiğiniz executor'a bağlı olarak uygun talimatlara çevirecektir.

Bunu çalıştırmayacağız, çünkü eğitim ortamında ilgili altyapıya erişimimiz yok.
Ancak, bu sınırları aşan kaynak tahsisleriyle workflow'u çalıştırmaya çalışsaydınız, ardından `.command.run` betik dosyasındaki `sbatch` komutuna baksaydınız, executor'a gerçekten gönderilen isteklerin `resourceLimits` tarafından belirtilen değerlerle sınırlandırıldığını görürdünüz.

??? info "Kurumsal referans yapılandırmaları"

    nf-core projesi, çok çeşitli HPC ve bulut executor'larını kapsayan, dünya genelindeki çeşitli kurumlar tarafından paylaşılan bir [yapılandırma dosyaları koleksiyonu](https://nf-co.re/configs/) derlemiştir.

    Bu paylaşılan yapılandırmalar hem orada çalışan ve bu nedenle kurumlarının yapılandırmasını hazır olarak kullanabilen insanlar için hem de kendi altyapıları için bir yapılandırma geliştirmek isteyen insanlar için bir model olarak değerlidir.

### Özet

Kaynak kullanımını değerlendirmek için bir profilleme raporu oluşturmayı ve tüm process'ler ve/veya bireysel process'ler için kaynak tahsislerini nasıl değiştireceğinizi ve HPC'de çalıştırmak için kaynak sınırlamalarını nasıl ayarlayacağınızı biliyorsunuz.

### Sırada ne var?

Önceden ayarlanmış yapılandırma profillerini nasıl kuracağınızı ve çalışma zamanında bunlar arasında nasıl geçiş yapacağınızı öğrenin.

---

## 6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın

??? example "Senaryo"

    Geliştirme için dizüstü bilgisayarınızda ve üretim çalıştırmaları için kurumunuzun HPC'sinde pipeline'ları çalıştırmak arasında düzenli olarak geçiş yapıyorsunuz.
    Her ortam değiştirdiğinizde yapılandırma ayarlarını manuel olarak değiştirmekten yoruldunuz.

Size pipeline yapılandırmanızı üzerinde çalıştığınız projeye veya kullandığınız hesaplama ortamına bağlı olarak özelleştirebileceğiniz birkaç yol gösterdik.

Hangi hesaplama altyapısını kullandığınıza bağlı olarak alternatif ayarlar arasında geçiş yapmak isteyebilirsiniz. Örneğin, dizüstü bilgisayarınızda küçük ölçekli testler geliştirmek ve yerel olarak çalıştırmak, ardından HPC veya bulutta tam ölçekli iş yükleri çalıştırmak isteyebilirsiniz.

Nextflow, farklı yapılandırmaları tanımlayan istediğiniz sayıda profil kurmanıza olanak tanır; bunları daha sonra yapılandırma dosyasını değiştirmek yerine bir komut satırı argümanı kullanarak çalışma zamanında seçebilirsiniz.

### 6.1. Yerel geliştirme ve HPC'de çalıştırma arasında geçiş yapmak için profiller oluşturun

İki alternatif profil kuralım; biri Docker konteynerları kullanacağımız normal bir bilgisayarda küçük ölçekli yükler çalıştırmak için, biri de Conda paketleri kullanacağımız Slurm scheduler'lı bir üniversite HPC'sinde çalıştırmak için.

#### 6.1.1. Profilleri kurun

`nextflow.config` dosyanıza, pipeline parametreleri bölümünden sonra ancak çıktı ayarlarından önce aşağıdakini ekleyin:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

Üniversite HPC için ayrıca kaynak sınırlamalarını da belirttiğimizi görüyorsunuz.

#### 6.1.2. Workflow'u bir profille çalıştırın

Nextflow komut satırımızda bir profil belirtmek için `-profile` argümanını kullanıyoruz.

Workflow'u `my_laptop` yapılandırmasıyla çalıştırmayı deneyelim.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Gördüğünüz gibi, bu bize çalışma zamanında yapılandırmalar arasında çok rahatça geçiş yapma imkanı veriyor.

!!! warning "Uyarı"

    `univ_hpc` profili, bir Slurm scheduler'a erişimimiz olmadığı için eğitim ortamında düzgün çalışmayacaktır.

Gelecekte bunlarla her zaman birlikte oluşan başka yapılandırma öğeleri bulursak, bunları ilgili profil(ler)e ekleyebiliriz.
Ayrıca birlikte gruplamak istediğimiz başka yapılandırma öğeleri varsa ek profiller oluşturabiliriz.

### 6.2. Test parametreleri profili oluşturun

??? example "Senaryo"

    Başkalarının kendi girdi verilerini toplamadan pipeline'ınızı hızlıca denemesini istiyorsunuz.

Profiller yalnızca altyapı yapılandırması için değildir.
Başkalarının uygun girdi değerlerini kendileri toplamak zorunda kalmadan workflow'u denemelerini kolaylaştırmak için workflow parametreleri için varsayılan değerler ayarlamak için de kullanabiliriz.
Bunu bir parametre dosyası kullanmaya alternatif olarak düşünebilirsiniz.

#### 6.2.1. Profili kurun

Bu bağlamda varsayılan değerleri ifade etmek için sözdizimi, `test` olarak adlandırdığımız bir profil için şöyle görünür:

```groovy title="Sözdizimi örneği"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Workflow'umuz için bir test profili eklersek, `profiles` bloğu şöyle olur:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Teknik yapılandırma profilleri için olduğu gibi, istediğiniz herhangi bir ad altında parametreleri belirten birden fazla farklı profil kurabilirsiniz.

#### 6.2.2. Workflow'u yerel olarak test profiliyle çalıştırın

Uygun bir şekilde, profiller karşılıklı olarak dışlayıcı değildir, bu yüzden komut satırımızda `-profile <profile1>,<profile2>` sözdizimini kullanarak birden fazla profil belirtebiliriz (herhangi bir sayıda profil için).

Aynı yapılandırma öğeleri için değerler ayarlayan ve aynı yapılandırma dosyasında tanımlanan profilleri birleştirirseniz, Nextflow hangi değeri en son okuduğunu kullanarak çakışmayı çözecektir (yani dosyada daha sonra gelen).
Çakışan ayarlar farklı yapılandırma kaynaklarında ayarlanmışsa, varsayılan [öncelik sırası](https://www.nextflow.io/docs/latest/config.html) geçerlidir.

Önceki komutumza test profilini eklemeyi deneyelim:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Bu, mümkün olduğunda Docker kullanacak ve `results/test` altında çıktılar üretecek ve bu sefer karakter komik ikili `dragonandcow`.

??? abstract "Dosya içeriği"

    ```console title="results/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Bu, herhangi bir test veri dosyasını workflow koduyla birlikte dağıttığımız sürece, herkesin komut satırı veya parametre dosyası aracılığıyla kendi girdilerini sağlamak zorunda kalmadan workflow'u hızlıca deneyebileceği anlamına gelir.

!!! tip "İpucu"

    Harici olarak depolanan daha büyük dosyalar için URL'lere işaret edebiliriz.
    Nextflow, açık bir bağlantı olduğu sürece bunları otomatik olarak indirecektir.

    Daha fazla ayrıntı için [Working with Files](../side_quests/working_with_files.md) Yan Görevine bakın.

### 6.3. Çözülmüş yapılandırmayı görmek için `nextflow config` kullanın

Yukarıda belirtildiği gibi, bazen aynı parametre birleştirmek istediğiniz profillerde farklı değerlere ayarlanabilir.
Ve daha genel olarak, yapılandırma öğelerinin saklanabileceği çok sayıda yer vardır ve bazen aynı özellikler farklı yerlerde farklı değerlere ayarlanabilir.

Nextflow, herhangi bir çakışmayı çözmek için belirli bir [öncelik sırası](https://www.nextflow.io/docs/latest/config.html) uygular, ancak bunu kendiniz belirlemeniz zor olabilir.
Ve hiçbir şey çakışmasa bile, şeylerin yapılandırılabileceği tüm olası yerlere bakmak sıkıcı olabilir.

Neyse ki, Nextflow bu süreci sizin için otomatikleştirebilen `config` adlı uygun bir yardımcı araç içerir.

`config` aracı, mevcut çalışma dizininizdeki tüm içeriği keşfedecek, tüm yapılandırma dosyalarını toplayacak ve Nextflow'un workflow'u çalıştırmak için kullanacağı tamamen çözülmüş yapılandırmayı üretecektir.
Bu, herhangi bir şey başlatmak zorunda kalmadan hangi ayarların kullanılacağını öğrenmenizi sağlar.

#### 6.3.1. Varsayılan yapılandırmayı çözün

Varsayılan olarak uygulanacak yapılandırmayı çözmek için bu komutu çalıştırın.

```bash
nextflow config
```

??? success "Komut çıktısı"

    ```groovy
    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    params {
      input = 'greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }
    ```

Bu, komut satırında ekstra bir şey belirtmezseniz elde ettiğiniz temel yapılandırmayı gösterir.

#### 6.3.2. Belirli ayarlar etkinleştirilmiş yapılandırmayı çözün

Komut satırı parametreleri sağlarsanız, örneğin bir veya daha fazla profil etkinleştirmek veya bir parametre dosyası yüklemek, komut bunları da dikkate alacaktır.

```bash
nextflow config -profile my_laptop,test
```

??? success "Komut çıktısı"

    ```groovy
    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    params {
      input = 'greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }
    ```

Bu, birden fazla yapılandırma katmanı içeren karmaşık projeler için özellikle yararlı olur.

### Özet

Çalışma zamanında minimum güçlükle önceden ayarlanmış bir yapılandırmayı seçmek için profilleri nasıl kullanacağınızı biliyorsunuz.
Daha genel olarak, workflow çalıştırmalarınızı farklı hesaplama platformlarına uyacak şekilde ve analizlerinizin tekrar üretilebilirliğini artırmak için nasıl yapılandıracağınızı biliyorsunuz.

### Sırada ne var?

Pipeline'ları doğrudan GitHub gibi uzak havuzlardan nasıl çalıştıracağınızı öğrenin.

---

## 7. Uzak havuzlardan pipeline'ları çalıştırma

??? example "Senaryo"

    nf-core'dan olanlar gibi köklü bir pipeline'ı kodu kendiniz indirmeden ve yönetmeden çalıştırmak istiyorsunuz.

Şu ana kadar mevcut dizinde bulunan workflow betiklerini çalıştırıyorduk.
Pratikte, genellikle GitHub gibi uzak havuzlarda depolanan pipeline'ları çalıştırmak isteyeceksiniz.

Nextflow bunu kolaylaştırır: herhangi bir pipeline'ı önce manuel olarak indirmeden doğrudan bir Git havuzu URL'sinden çalıştırabilirsiniz.

### 7.1. GitHub'dan bir pipeline çalıştırma

Uzak bir pipeline çalıştırmak için temel sözdizimi `nextflow run <repository>`'dir; burada `<repository>` `nextflow-io/hello` gibi bir GitHub havuz yolu, tam bir URL veya GitLab, Bitbucket veya diğer Git barındırma hizmetlerine giden bir yol olabilir.

Resmi Nextflow "hello" demo pipeline'ını çalıştırmayı deneyin:

```bash
nextflow run nextflow-io/hello
```

Uzak bir pipeline'ı ilk kez çalıştırdığınızda, Nextflow onu indirir ve yerel olarak önbelleğe alır.
Sonraki çalıştırmalar, açıkça bir güncelleme istemediğiniz sürece önbelleğe alınmış versiyonu kullanır.

### 7.2. Tekrar üretilebilirlik için versiyon belirtin

Varsayılan olarak, Nextflow varsayılan daldan en son versiyonu çalıştırır.
`-r` bayrağını kullanarak belirli bir versiyon, dal veya commit belirtebilirsiniz:

```bash
nextflow run nextflow-io/hello -r v1.1
```

Tam versiyonları belirtmek tekrar üretilebilirlik için önemlidir.

### Özet

Pipeline'ları doğrudan GitHub ve diğer uzak havuzlardan nasıl çalıştıracağınızı ve tekrar üretilebilirlik için versiyonları nasıl belirteceğinizi biliyorsunuz.

### Sırada ne var?

Kendinize büyük bir tebrik verin!
Nextflow pipeline'larını çalıştırmaya ve yönetmeye başlamak için bilmeniz gereken her şeyi biliyorsunuz.

Bu, bu kursu sonlandırıyor, ancak öğrenmeye devam etmek istiyorsanız, iki ana önerimiz var:

- Kendi pipeline'larınızı geliştirmeyi daha derinlemesine incelemek istiyorsanız, bu kursla aynı genel ilerlemeyi kapsayan ancak channel'lar ve operatörler hakkında çok daha ayrıntılı giden yeni başlayanlar için bir kurs olan [Hello Nextflow](../hello_nextflow/index.md)'a bakın.
- Koda daha derinlemesine girmeden Nextflow pipeline'larını çalıştırmayı öğrenmeye devam etmek istiyorsanız, son derece popüler [nf-core](https://nf-co.re/) projesinden pipeline'ları bulmak ve çalıştırmak için araçları tanıtan [Hello nf-core](../hello_nf-core/index.md)'un ilk bölümüne bakın.

İyi eğlenceler!

---

## Quiz

<quiz>
Parametre değerleri hem workflow dosyasında hem de `nextflow.config`'de ayarlandığında, hangisi önceliklidir?
- [ ] Workflow dosyası değeri
- [x] Yapılandırma dosyası değeri
- [ ] Karşılaşılan ilk değer
- [ ] Bir hataya neden olur

Daha fazla bilgi: [1.1. `nextflow.config`'de değerler ayarlayın](#11-nextflowconfigde-degerler-ayarlayin)
</quiz>

<quiz>
Workflow dosyasında vs. yapılandırma dosyasında varsayılan parametre ayarlama arasındaki sözdizimi farkı nedir?
- [ ] Aynı sözdizimini kullanırlar
- [x] Workflow türlenmiş bildirim kullanır (`#!groovy param: Type = value`), yapılandırma atama kullanır (`#!groovy param = value`)
- [ ] Yapılandırma türlenmiş bildirim kullanır, workflow atama kullanır
- [ ] Yalnızca yapılandırma dosyaları varsayılan değerler ayarlayabilir

Daha fazla bilgi: [1.1. `nextflow.config`'de değerler ayarlayın](#11-nextflowconfigde-degerler-ayarlayin)
</quiz>

<quiz>
Bir workflow çalıştırırken parametre dosyası nasıl belirtilir?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Daha fazla bilgi: [1.3. Parametre dosyası kullanın](#13-parametre-dosyasi-kullanin)
</quiz>

<quiz>
`outputDir` yapılandırma seçeneği neyi kontrol eder?
- [ ] work dizininin konumu
- [x] Workflow çıktılarının yayınlandığı temel yol
- [ ] Log dosyaları için dizin
- [ ] Modül dosyalarının konumu

Daha fazla bilgi: [2.1. outputDir dizin adını özelleştirme](#21-outputdir-dizin-adini-ozellestirme)
</quiz>

<quiz>
Çıktı yolu yapılandırmasında bir process adına dinamik olarak nasıl referans verilir?
- [ ] `#!groovy ${processName}`
- [ ] `process.name`
- [x] `#!groovy { meta.id }`
- [ ] `@processName`

Daha fazla bilgi: [2.2. Çıktıları process'e göre organize etme](#22-ciktilari-processe-gore-organize-etme)
</quiz>

<quiz>
Hem Docker hem de Conda etkinleştirilmişse ve bir process'in her iki direktifi de varsa, hangisi önceliklendirilir?
- [x] Docker (konteynerlar)
- [ ] Conda
- [ ] Process'te tanımlanan ilki
- [ ] Bir hataya neden olur

Daha fazla bilgi: [3. Yazılım paketleme teknolojisi seçme](#3-yazilim-paketleme-teknolojisi-secme)
</quiz>

<quiz>
Nextflow'daki varsayılan executor nedir?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Daha fazla bilgi: [4. Çalıştırma platformu seçme](#4-calistirma-platformu-secme)
</quiz>

<quiz>
Kaynak kullanım raporu oluşturan komut nedir?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Daha fazla bilgi: [5.1. Kaynak kullanım raporu oluşturmak için workflow'u çalıştırın](#51-kaynak-kullanim-raporu-olusturmak-icin-workflowu-calistirin)
</quiz>

<quiz>
Yapılandırma dosyasında `cowpy` adlı belirli bir process için kaynak gereksinimleri nasıl ayarlanır?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Daha fazla bilgi: [5.3. Belirli bir process için kaynak tahsislerini ayarlama](#53-belirli-bir-process-icin-kaynak-tahsislerini-ayarlama)
</quiz>

<quiz>
`resourceLimits` direktifi ne yapar?
- [ ] Minimum kaynak gereksinimlerini ayarlar
- [ ] Process'lere kaynak tahsis eder
- [x] İstenebilecek maksimum kaynakları sınırlar
- [ ] Kaynak kullanımını gerçek zamanlı izler

Daha fazla bilgi: [5.5. Kaynak sınırları ekleme](#55-kaynak-sinirlari-ekleme)
</quiz>

<quiz>
Tek bir komutta birden fazla profil nasıl belirtilir?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Daha fazla bilgi: [6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın](#6-onceden-ayarlanmis-yapilandirmalar-arasinda-gecis-yapmak-icin-profilleri-kullanin)
</quiz>

<quiz>
Nextflow'un kullanacağı tamamen çözülmüş yapılandırmayı gösteren komut nedir?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Daha fazla bilgi: [6.3. Çözülmüş yapılandırmayı görmek için `nextflow config` kullanın](#63-cozulmus-yapilandirmayi-gormek-icin-nextflow-config-kullanin)
</quiz>

<quiz>
Profiller ne için kullanılabilir? (Geçerli olanların hepsini seçin)
- [x] Altyapıya özgü ayarları tanımlama (executor'lar, konteynerlar)
- [x] Farklı ortamlar için kaynak sınırlarını ayarlama
- [x] Kolay workflow testi için test parametreleri sağlama
- [ ] Yeni process'ler tanımlama

Daha fazla bilgi: [6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın](#6-onceden-ayarlanmis-yapilandirmalar-arasinda-gecis-yapmak-icin-profilleri-kullanin)
</quiz>
