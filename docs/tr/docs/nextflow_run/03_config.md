# Bölüm 3: Çalıştırma yapılandırması

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } YZ destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu bölüm, bir Nextflow iş akışının yapılandırmasını nasıl yöneteceğinizi keşfedecektir; böylece davranışını özelleştirebilir, farklı ortamlara uyarlayabilir ve kaynak kullanımını optimize edebilirsiniz - _iş akışı kodunun kendisinde tek bir satır bile değiştirmeden_.

Bunu yapmanın birden fazla yolu vardır ve bunlar birlikte kullanılabilir ve [Configuration](https://nextflow.io/docs/latest/config.html) belgelerinde açıklanan öncelik sırasına göre yorumlanır.

Kursun bu bölümünde, size en basit ve en yaygın yapılandırma dosyası mekanizmasını göstereceğiz: `nextflow.config` dosyası. Bu dosyayla Bölüm 2'deki konteynerler bölümünde zaten karşılaştınız.

Süreç yönergeleri, yürütücüler, profiller ve parametre dosyaları gibi Nextflow yapılandırmasının temel bileşenlerini ele alacağız.
Bu yapılandırma seçeneklerini etkili bir şekilde kullanmayı öğrenerek, Nextflow iş akışlarının esnekliğinden, ölçeklenebilirliğinden ve performansından tam olarak yararlanabilirsiniz.

Bu yapılandırma öğelerini uygulamak için, bu eğitim kursunun Bölüm 2'sinin sonunda çalıştırdığımız iş akışının yeni bir kopyasını çalıştıracağız; adı `3-main.nf` olarak değiştirilmiş.

Hello iş akışına aşina değilseniz veya bir hatırlatmaya ihtiyacınız varsa, [bu bilgi sayfasına](../info/hello_pipeline.md) bakın.

---

## 1. İş akışı girdi parametrelerini yönetme

??? example "Senaryo"

    Bir iş akışı indirdiniz ve aynı girdi dosyaları ve ayarlarla tekrar tekrar çalıştırmak istiyorsunuz, ancak her seferinde tüm parametreleri yazmak istemiyorsunuz.
    Ya da belki de komut satırı argümanlarıyla rahat olmayan bir meslektaşınız için iş akışını kuruyorsunuz.

Şimdiye kadar üzerinde çalıştığımız şeyin bir uzantısı olan bir yapılandırma yönüyle başlayacağız: girdi parametrelerinin yönetimi.

Şu anda, iş akışımız komut satırı aracılığıyla birkaç parametre değerini kabul edecek şekilde ayarlanmış durumda; bunlar iş akışı betiğinin kendisindeki bir `params` bloğunda bildirilmiş.
Birinin bildirimin bir parçası olarak ayarlanmış varsayılan bir değeri var.

Ancak, hepsi için varsayılan değerler ayarlamak veya mevcut varsayılanı komut satırında parametre belirtmek ya da orijinal betik dosyasını değiştirmek zorunda kalmadan geçersiz kılmak isteyebilirsiniz.

Bunu yapmanın birden fazla yolu vardır; size çok yaygın olarak kullanılan üç temel yolu göstereceğiz.

### 1.1. `nextflow.config` dosyasında değerleri ayarlama

Bu en basit yaklaşımdır, ancak ana `nextflow.config` dosyası her çalıştırma için düzenlemek isteyeceğiniz bir şey olmadığından muhtemelen en az esnektir.
Ancak parametreleri iş akışında _bildirme_ (ki bu kesinlikle oraya aittir) ile _varsayılan değerleri_ sağlama (bunlar bir yapılandırma dosyasında daha uygun) endişelerini ayırma avantajına sahiptir.

Bunu iki adımda yapalım.

#### 1.1.1. Yapılandırma dosyasında bir `params` bloğu oluşturma

`nextflow.config` dosyasında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
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

`params` bloğunu iş akışından yapılandırma dosyasına basitçe kopyalamadığımıza dikkat edin.
Zaten varsayılan bir değer bildirimi olan `batch` parametresi için sözdizimi biraz farklı.
İş akışı dosyasında, bu türü belirtilmiş bir bildirimdir.
Yapılandırmada ise bunlar değer atamalardır.

Teknik olarak, bu iş akışı dosyasında hala belirtilen varsayılan değerleri geçersiz kılmak için yeterlidir.
`batch` için varsayılan değeri değiştirebilir ve yapılandırma dosyasında ayarlanan değerin iş akışı dosyasında ayarlanan değeri geçersiz kıldığından emin olmak için iş akışını çalıştırabilirsiniz.

Ancak yapılandırmayı tamamen yapılandırma dosyasına taşıma ruhuyla, iş akışı dosyasından bu varsayılan değeri tamamen kaldıralım.

#### 1.1.2. İş akışı dosyasındaki `batch` için varsayılan değeri kaldırma

`3-main.nf` iş akışı dosyasında aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
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
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Artık iş akışı dosyasının kendisi bu parametreler için herhangi bir varsayılan değer ayarlamıyor.

#### 1.1.3. İş akışını çalıştırma

Komut satırında herhangi bir parametre belirtmeden doğru çalıştığını test edelim.

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

Bu hala daha önce olduğu gibi aynı çıktıyı üretir.

Son ASCII sanat çıktısı `results/3-main/` dizininde, `cowpy-COLLECTED-batch-output.txt` adı altında, daha önce olduğu gibi.

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

İşlevsel olarak, bu hareket hiçbir şeyi değiştirmedi, ancak kavramsal olarak varsayılan değerlerin yapılandırma dosyasında ayarlanması biraz daha temiz.

### 1.2. Çalıştırmaya özgü bir yapılandırma dosyası kullanma

??? example "Senaryo"

    Ana yapılandırma dosyanızı değiştirmeden farklı ayarlarla denemeler yapmak istiyorsunuz.

Bunu, deneyler için çalışma dizini olarak kullanacağınız bir alt dizinde yeni bir `nextflow.config` dosyası oluşturarak yapabilirsiniz.

#### 1.2.1. Boş bir yapılandırmayla çalışma dizini oluşturma

Yeni bir dizin oluşturarak ve içine girerek başlayalım:

```bash
mkdir -p tux-run
cd tux-run
```

Ardından, bu dizinde boş bir yapılandırma dosyası oluşturun:

```bash
touch nextflow.config
```

Bu boş bir dosya üretir.

#### 1.2.2. Deneysel yapılandırmayı ayarlama

Şimdi yeni dosyayı açın ve özelleştirmek istediğiniz parametreleri ekleyin:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Girdi dosyasının yolunun dizin yapısını yansıtması gerektiğine dikkat edin.

#### 1.2.3. İş akışını çalıştırma

Artık iş akışımızı yeni çalışma dizinimizin içinden çalıştırabiliriz.
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

Bu, `tux-run/work/` ve `tux-run/results/` dahil olmak üzere `tux-run/` altında yeni bir dizin seti oluşturacaktır.

Bu çalıştırmada, Nextflow mevcut dizinimizdeki `nextflow.config` dosyasını iş akışının kök dizinindeki `nextflow.config` ile birleştirir ve böylece varsayılan karakteri (turkey) tux karakteriyle geçersiz kılar.

Son çıktı dosyası, selamlamaları söyleyen tux karakterini içermelidir.

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

!!! warning

    Bir sonraki bölüme geçmeden önce önceki dizine geri döndüğünüzden emin olun!

    ```bash
    cd ..
    ```

Şimdi parametre değerlerini ayarlamanın başka bir kullanışlı yoluna bakalım.

### 1.3. Bir parametre dosyası kullanma

??? example "Senaryo"

    Tam çalıştırma parametrelerini bir meslektaşınızla paylaşmanız veya bir yayın için kaydetmeniz gerekiyor.

Alt dizin yaklaşımı denemeler için harika çalışır, ancak biraz kurulum gerektirir ve yolları buna göre uyarlamanızı gerektirir.
İş akışınızı belirli bir değer setiyle çalıştırmak istediğinizde veya başka birinin bunu minimum çabayla yapmasını sağlamak istediğinizde daha basit bir yaklaşım vardır.

Nextflow, YAML veya JSON formatında bir [parametre dosyası](https://nextflow.io/docs/latest/config.html#parameter-file) aracılığıyla parametreleri belirtmemize olanak tanır; bu da örneğin alternatif varsayılan değer setlerini ve çalıştırmaya özgü parametre değerlerini yönetmeyi ve dağıtmayı çok kolaylaştırır.

#### 1.3.1. Örnek parametre dosyasını inceleme

Bunu göstermek için, mevcut dizinde `test-params.yaml` adlı bir örnek parametre dosyası sağlıyoruz:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Bu parametre dosyası, belirtmek istediğimiz her girdi için bir anahtar-değer çifti içerir.
Yapılandırma dosyasıyla karşılaştırırsanız, eşittir işaretleri (`=`) yerine iki nokta üst üste (`:`) kullanımına dikkat edin.
Yapılandırma dosyası Groovy'de yazılırken, parametre dosyası YAML'de yazılır.

!!! info

    Ayrıca örnek olarak parametre dosyasının bir JSON sürümünü de sağlıyoruz ancak burada bununla çalıştırmayacağız.
    Bunu kendi başınıza denemekten çekinmeyin.

#### 1.3.2. İş akışını çalıştırma

İş akışını bu parametre dosyasıyla çalıştırmak için, temel komuta basitçe `-params-file <dosyaadı>` ekleyin.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

Son çıktı dosyası, selamlamaları söyleyen stegosaurus karakterini içermelidir.

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

Bir parametre dosyası kullanmak, belirtilecek yalnızca birkaç parametreniz olduğunda aşırıya kaçmak gibi görünebilir, ancak bazı iş akışları düzinelerce parametre bekler.
Bu durumlarda, bir parametre dosyası kullanmak, çalışma zamanında büyük komut satırları yazmak ve iş akışı betiğini değiştirmek zorunda kalmadan parametre değerleri sağlamamıza olanak tanır.

Ayrıca parametre setlerini meslektaşlara dağıtmayı veya örneğin bir yayın için destekleyici bilgi olarak sunmayı kolaylaştırır.
Bu, çalışmanızı başkaları tarafından daha tekrarlanabilir hale getirir.

### Özet

İş akışı girdilerini yönetmek için temel yapılandırma seçeneklerinden nasıl yararlanacağınızı biliyorsunuz.

### Sırada ne var?

İş akışı çıktılarınızın nerede ve nasıl yayınlanacağını nasıl yöneteceğinizi öğrenin.

---

## 2. İş akışı çıktılarını yönetme

??? example "Senaryo"

    İş akışınız çıktıları sabit kodlanmış bir dizine yayınlıyor, ancak her seferinde iş akışı kodunu düzenlemeden sonuçları proje veya deney adına göre düzenlemek istiyorsunuz.

Miras aldığımız iş akışı, iş akışı düzeyinde çıktı bildirimleri için yollar kullanıyor; bu pek esnek değil ve çok fazla tekrar içeriyor.

Bunu daha esnek hale getirmek için yapılandırabileceğiniz birkaç yaygın yola bakalım.

### 2.1. `outputDir` dizin adını özelleştirme

Şimdiye kadar çalıştırdığımız iş akışının her sürümü, çıktılarını çıktı tanımlarına sabit kodlanmış farklı bir alt dizine yayınladı.

Bölüm 1'de `-output-dir` CLI bayrağını kullanarak bu alt dizinin nerede olduğunu değiştirdik, ancak bu hala sadece statik bir dize.
Bunun yerine bunu bir yapılandırma dosyasında yapılandıralım; burada daha karmaşık dinamik yollar tanımlayabiliriz.
Bunun için tamamen yeni bir parametre oluşturabiliriz, ancak zaten orada olduğu için `batch` parametresini kullanalım.

#### 2.1.1. Yapılandırma dosyasında `outputDir` için bir değer ayarlama

Nextflow'un çıktıları yayınlamak için kullandığı yol `outputDir` seçeneği tarafından kontrol edilir.
Tüm çıktılar için yolu değiştirmek için, `nextflow.config` yapılandırma dosyasında bu seçenek için bir değer ayarlayabilirsiniz.

`nextflow.config` dosyasına aşağıdaki kodu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Bu, yerleşik varsayılan yol olan `results/` yerine `results_config/` artı `batch` parametresinin değerini alt dizin olarak kullanacaktır.

Bu seçeneği komut satırından komutunuzda `-output-dir` parametresini kullanarak da ayarlayabileceğinizi unutmayın (`-o` kısaca), ancak o zaman `batch` parametre değerini kullanamazsınız.
CLI bayrağını kullanmak, ayarlanmışsa yapılandırmadaki `outputDir`'i üzerine yazacaktır.

#### 2.1.2. Sabit kodlanmış yolun tekrarlanan kısmını kaldırma

Çıktı seçeneklerinde hala sabit kodlanmış bir alt dizinimiz var, o yüzden şimdi bundan kurtulalım.

İş akışı dosyasında aşağıdaki kod değişikliklerini yapın:

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

`outputDir` varsayılanını değiştirmek yerine her yola sadece `${params.batch}` ekleyebilirdik, ancak bu daha özlü.

#### 2.1.3. İş akışını çalıştırma

Doğru çalıştığını test edelim, batch adını komut satırından `outdir` olarak ayarlayalım.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretir, ancak bu sefer çıktılarımızı `results_config/outdir/` altında buluyoruz.

??? abstract "Dizin içeriği"

    ```console
    results_config/outdir
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

### 2.2. Çıktıları sürece göre düzenleme

Çıktıları daha fazla düzenlemenin popüler bir yolu, bunu sürece göre yapmaktır, _yani_ iş akışında çalıştırılan her süreç için alt dizinler oluşturmak.

#### 2.2.1. Çıktı yollarını süreç adlarına referansla değiştirme

Tek yapmanız gereken, çıktı yolu bildiriminde sürecin adını `<process>.name` olarak referans göstermektir.

İş akışı dosyasında aşağıdaki değişiklikleri yapın:

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

#### 2.2.2. İş akışını çalıştırma

Doğru çalıştığını test edelim, batch adını komut satırından `pnames` olarak ayarlayalım.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretir, ancak bu sefer çıktılarımızı `results_config/pnames/` altında buluyoruz ve sürece göre gruplandırılmışlar.

??? abstract "Dizin içeriği"

    ```console
    results_config/pnames/
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

!!! note

    Burada `intermediates` ile son çıktıların üst düzeyde olması arasındaki ayrımı sildiğimize dikkat edin.
    Bu yaklaşımları karıştırıp eşleştirebilir ve hatta birden fazla değişken ekleyebilirsiniz, örneğin ilk çıktının yolunu `#!groovy "${params.batch}/intermediates/${sayHello.name}"` olarak ayarlayarak

### 2.3. Yayınlama modunu iş akışı düzeyinde ayarlama

Son olarak, tekrarlayan kod miktarını azaltma ruhuyla, çıktı başına `mode` bildirimlerini yapılandırmada tek bir satırla değiştirebiliriz.

#### 2.3.1. Yapılandırma dosyasına `workflow.output.mode` ekleme

`nextflow.config` dosyasına aşağıdaki kodu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

Tıpkı `outputDir` seçeneği gibi, yapılandırma dosyasında `workflow.output.mode`'a bir değer vermek iş akışı dosyasında ayarlanmış olanı geçersiz kılmak için yeterli olurdu, ancak yine de gereksiz kodu kaldıralım.

#### 2.3.2. İş akışı dosyasından çıktı modunu kaldırma

İş akışı dosyasında aşağıdaki değişiklikleri yapın:

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

#### 2.3.3. İş akışını çalıştırma

Doğru çalıştığını test edelim, batch adını komut satırından `outmode` olarak ayarlayalım.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretir, ancak bu sefer çıktılarımızı `results_config/outmode/` altında buluyoruz.
Hepsi hala uygun kopyalar, sembolik bağlantılar değil.

??? abstract "Dizin içeriği"

    ```console
    results_config/outmode/
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

Çıktı başına mod ayarlama yöntemini kullanmak isteyebileceğiniz ana neden, aynı iş akışı içinde karıştırmak istemenizdir, _yani_ bazı çıktıların kopyalanmasını ve bazılarının sembolik bağlantı olmasını istemenizdir.

Bu şekilde özelleştirebileceğiniz başka birçok seçenek vardır, ancak umarım bu size seçenek yelpazesi ve tercihlerinize uygun şekilde bunları etkili bir şekilde nasıl kullanacağınız konusunda bir fikir verir.

### Özet

Çıktılarınızın yayınlandığı dizinlerin adlandırmasını ve yapısını, ayrıca iş akışı çıktı yayınlama modunu nasıl kontrol edeceğinizi biliyorsunuz.

### Sırada ne var?

İş akışı yapılandırmanızı bilgi işlem ortamınıza nasıl uyarlayacağınızı öğrenin, yazılım paketleme teknolojisiyle başlayarak.

---

## 3. Bir yazılım paketleme teknolojisi seçme

Şimdiye kadar girdilerin nasıl girdiğini ve girdilerin nereye çıktığını kontrol eden yapılandırma öğelerine bakıyorduk. Şimdi iş akışı yapılandırmanızı bilgi işlem ortamınıza uyarlamaya daha spesifik olarak odaklanmanın zamanı geldi.

Bu yoldaki ilk adım, her adımda çalıştırılacak yazılım paketlerinin nereden geleceğini belirtmektir.
Yerel bilgi işlem ortamında zaten yüklü mü?
Görüntüleri almamız ve bir konteyner sistemi aracılığıyla çalıştırmamız mı gerekiyor?
Yoksa Conda paketlerini almamız ve yerel bir Conda ortamı oluşturmamız mı gerekiyor?

Bu eğitim kursunun ilk bölümünde (Bölüm 1-4) iş akışımızda sadece yerel olarak yüklenmiş yazılımı kullandık.
Ardından Bölüm 5'te Docker konteynerlerini ve `nextflow.config` dosyasını tanıttık; bunu Docker konteynerlerinin kullanımını etkinleştirmek için kullandık.

Şimdi `nextflow.config` dosyası aracılığıyla alternatif bir yazılım paketleme seçeneğini nasıl yapılandırabileceğimizi görelim.

### 3.1. Yapılandırma dosyasında Docker'ı devre dışı bırakma ve Conda'yı etkinleştirme

??? example "Senaryo"

    İş akışınızı, güvenlik nedenleriyle Docker'a izin verilmeyen bir HPC kümesine taşıyorsunuz.
    Küme Singularity ve Conda'yı destekliyor, bu nedenle yapılandırmanızı buna göre değiştirmeniz gerekiyor.

Daha önce belirtildiği gibi, Nextflow, HPC'de daha yaygın olarak kullanılan Singularity dahil olmak üzere birden fazla konteyner teknolojisini ve Conda gibi yazılım paket yöneticilerini destekler.

Yapılandırma dosyamızı Docker yerine Conda kullanacak şekilde değiştirebiliriz.
Bunu yapmak için, `docker.enabled` değerini `false` olarak değiştirelim ve Conda kullanımını etkinleştiren bir yönerge ekleyelim:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Bu, Nextflow'un Conda paketleri belirtilmiş süreçler için Conda ortamları oluşturmasına ve kullanmasına olanak tanır.
Bu da şimdi `cowpy` sürecimize bunlardan birini eklememiz gerektiği anlamına gelir!

### 3.2. Süreç tanımında bir Conda paketi belirtme

`cowpy` aracını içeren bir Conda paketi için URI'yi zaten aldık: `conda-forge::cowpy==1.1.5`

Şimdi URI'yi `conda` yönergesini kullanarak `cowpy` süreç tanımına ekleyelim:

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

Açık olmak gerekirse, `docker` yönergesini _değiştirmiyoruz_, alternatif bir seçenek _ekliyoruz_.

!!! tip

    Belirli bir conda paketi için URI'yi almanın birkaç farklı yolu vardır.
    Bir konteyner oluşturmayı planlamasanız bile kopyalayıp yapıştırabileceğiniz bir URI verecek olan [Seqera Containers](https://seqera.io/containers/) arama sorgusunu kullanmanızı öneririz.

### 3.3. Conda kullanabileceğini doğrulamak için iş akışını çalıştırma

Hadi deneyelim.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Komut çıktısı"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Bu sorunsuz çalışmalı ve `results_config/conda` altında daha önce olduğu gibi aynı çıktıları üretmelidir.

Perde arkasında, Nextflow Conda paketlerini almış ve ortamı oluşturmuştur; bu normalde biraz iş gerektirir; bu yüzden bunların hiçbirini kendimiz yapmak zorunda kalmamamız güzel!

!!! info

    `cowpy` paketi oldukça küçük olduğu için bu hızlı çalışır, ancak büyük paketlerle çalışıyorsanız, ilk seferinde normalden biraz daha uzun sürebilir ve konsol çıktısının tamamlanmadan önce bir dakika kadar 'takılı' kaldığını görebilirsiniz.
    Bu normaldir ve Nextflow'un yeni bir paketi ilk kez kullandığınızda yaptığı ekstra işten kaynaklanır.

Bizim açımızdan, arka planda mekanikler biraz farklı olsa da Docker ile çalıştırmakla tamamen aynı şekilde çalışıyor gibi görünüyor.

Bu, gerekirse Conda ortamlarıyla çalıştırmaya hazır olduğumuz anlamına gelir.

??? info "Docker ve Conda'yı karıştırma ve eşleştirme"

    Bu yönergeler süreç başına atandığından, 'karıştırma ve eşleştirme' mümkündür, _yani_ iş akışınızdaki bazı süreçleri Docker ile, diğerlerini örneğin Conda ile çalışacak şekilde yapılandırmak, eğer kullandığınız bilgi işlem altyapısı her ikisini de destekliyorsa.
    Bu durumda, yapılandırma dosyanızda hem Docker'ı hem de Conda'yı etkinleştirirsiniz.
    Belirli bir süreç için her ikisi de mevcutsa, Nextflow konteynerlere öncelik verecektir.

    Ve daha önce belirtildiği gibi, Nextflow sadece bu ikisiyle sınırlı olmayan birden fazla başka yazılım paketleme ve konteyner teknolojisini destekler.

### Özet

Her sürecin hangi yazılım paketini kullanması gerektiğini nasıl yapılandıracağınızı ve teknolojiler arasında nasıl geçiş yapacağınızı biliyorsunuz.

### Sırada ne var?

Nextflow tarafından işi gerçekten yapmak için kullanılan yürütme platformunu nasıl değiştireceğinizi öğrenin.

---

## 4. Bir yürütme platformu seçme

??? example "Senaryo"

    İş akışınızı dizüstü bilgisayarınızda geliştirip test ediyordunuz, ancak şimdi binlerce örnek üzerinde çalıştırmanız gerekiyor.
    Kurumunuzun bunun yerine kullanmak istediğiniz Slurm zamanlayıcılı bir HPC kümesi var.

Şimdiye kadar, iş akışımızı yerel yürütücü ile çalıştırıyorduk.
Bu, her görevi Nextflow'un üzerinde çalıştığı makinede yürütür.
Nextflow başladığında, mevcut CPU'lara ve belleğe bakar.
Çalıştırılmaya hazır görevlerin kaynakları mevcut kaynakları aşarsa, Nextflow son görevleri, önceki görevlerden biri veya daha fazlası bitene ve gerekli kaynakları serbest bırakana kadar yürütmeden geri tutar.

Yerel yürütücü kullanışlı ve verimlidir, ancak tek bir makineyle sınırlıdır. Çok büyük iş yükleri için, yerel makinenizin bir darboğaz olduğunu keşfedebilirsiniz; ya sahip olduğunuzdan daha fazla kaynak gerektiren tek bir göreviniz olduğu için ya da tek bir makinenin bunları çalıştırmasını beklemenin çok uzun süreceği kadar çok göreviniz olduğu için.

Nextflow, HPC zamanlayıcıları (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor ve diğerleri) ve bulut yürütme arka uçları (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes ve daha fazlası) dahil olmak üzere [birçok farklı yürütme arka ucunu](https://nextflow.io/docs/latest/executor.html) destekler.

### 4.1. Farklı bir arka ucu hedefleme

Yürütücü seçimi, `executor` adlı bir süreç yönergesi tarafından ayarlanır.
Varsayılan olarak `local` olarak ayarlanmıştır, bu nedenle aşağıdaki yapılandırma ima edilir:

```groovy title="Yerleşik yapılandırma"
process {
    executor = 'local'
}
```

Yürütücüyü farklı bir arka ucu hedefleyecek şekilde ayarlamak için, kaynak tahsisleri için yukarıda açıklandığı gibi benzer sözdizimini kullanarak istediğiniz yürütücüyü belirtirsiniz (tüm seçenekler için [Executors](https://nextflow.io/docs/latest/executor.html)'a bakın).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Eğitim ortamı bir HPC'ye bağlanacak şekilde ayarlanmadığı için bunu gerçekten test edemeyiz.

### 4.2. Yürütme parametreleri için arka uca özgü sözdizimi ile başa çıkma

Çoğu yüksek performanslı bilgi işlem platformu, belirli parametreleri (örneğin CPU sayısı ve bellek gibi kaynak tahsis istekleri ve sınırlamaları) ve kullanılacak iş kuyruğunun adını belirtmenize izin verir (ve bazen bunu gerektirir).

Ne yazık ki, bu sistemlerin her biri bir işin nasıl tanımlanması ve ilgili zamanlayıcıya gönderilmesi gerektiğini tanımlamak için farklı teknolojiler, sözdizimi ve yapılandırmalar kullanır.

??? abstract "Örnekler"

    Örneğin, "my-science-work" kuyruğunda yürütülmek üzere 8 CPU ve 4GB RAM gerektiren aynı işin, arka uca bağlı olarak aşağıdaki farklı şekillerde ifade edilmesi gerekir.

    ```bash title="SLURM için yapılandırma / sbatch kullanarak gönderme"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="PBS için yapılandırma / qsub kullanarak gönderme"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="SGE için yapılandırma / qsub kullanarak gönderme"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Neyse ki, Nextflow tüm bunları basitleştirir.
`cpus`, `memory` ve `queue` gibi ilgili özellikleri sadece bir kez belirtebilmeniz için standartlaştırılmış bir sözdizimi sağlar (tüm mevcut seçenekler için [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives)'e bakın).
Ardından, çalışma zamanında, Nextflow bu ayarları yürütücü ayarına göre uygun arka uca özgü betikleri oluşturmak için kullanacaktır.

Bu standartlaştırılmış sözdizimini bir sonraki bölümde ele alacağız.

### Özet

Artık farklı türde bilgi işlem altyapılarını kullanmak için yürütücüyü nasıl değiştireceğinizi biliyorsunuz.

### Sırada ne var?

Nextflow'da kaynak tahsislerini ve sınırlamalarını nasıl değerlendireceğinizi ve ifade edeceğinizi öğrenin.

---

## 5. Bilgi işlem kaynak tahsislerini kontrol etme

??? example "Senaryo"

    İş akışınız kümede başarısız olmaya devam ediyor çünkü görevler bellek sınırlarını aştıkları için sonlandırılıyor.
    Ya da belki kullanmadığınız kaynaklar için ücretlendiriliyorsunuz ve maliyetleri optimize etmek istiyorsunuz.

Çoğu yüksek performanslı bilgi işlem platformu, CPU sayısı ve bellek gibi belirli kaynak tahsis parametrelerini belirtmenize izin verir (ve bazen bunu gerektirir).

Varsayılan olarak, Nextflow her süreç için tek bir CPU ve 2GB bellek kullanacaktır.
Karşılık gelen süreç yönergeleri `cpus` ve `memory` olarak adlandırılır, bu nedenle aşağıdaki yapılandırma ima edilir:

```groovy title="Yerleşik yapılandırma" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Bu değerleri, tüm süreçler için veya belirli adlandırılmış süreçler için, yapılandırma dosyanızda ek süreç yönergeleri kullanarak değiştirebilirsiniz.
Nextflow bunları seçilen yürütücü için uygun talimatlara çevirecektir.

Ancak hangi değerleri kullanacağınızı nasıl bilirsiniz?

### 5.1. Bir kaynak kullanım raporu oluşturmak için iş akışını çalıştırma

??? example "Senaryo"

    Süreçlerinizin ne kadar bellek veya CPU'ya ihtiyacı olduğunu bilmiyorsunuz ve kaynakları boşa harcamaktan veya işlerin sonlandırılmasından kaçınmak istiyorsunuz.

Süreçlerinizin muhtemelen ne kadar CPU ve belleğe ihtiyaç duyacağını önceden bilmiyorsanız, biraz kaynak profilleme yapabilirsiniz; yani iş akışını bazı varsayılan tahsislerle çalıştırır, her sürecin ne kadar kullandığını kaydeder ve oradan temel tahsisleri nasıl ayarlayacağınızı tahmin edersiniz.

Uygun bir şekilde, Nextflow bunu yapmak için yerleşik araçlar içerir ve talep üzerine sizin için bir rapor oluşturmaktan mutluluk duyar.

Bunu yapmak için, komut satırınıza `-with-report <dosyaadı>.html` ekleyin.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Rapor bir html dosyasıdır; bunu indirebilir ve tarayıcınızda açabilirsiniz. Ayrıca soldaki dosya gezgininde sağ tıklayıp `Show preview`'a tıklayarak eğitim ortamında görüntüleyebilirsiniz.

Raporu incelemek için birkaç dakikanızı ayırın ve kaynakları ayarlamak için bazı fırsatlar belirleyip belirleyemeyeceğinize bakın.
Kullanım sonuçlarını tahsis edilenin yüzdesi olarak gösteren sekmelere tıkladığınızdan emin olun.

Mevcut tüm özelliklerle ilgili belgeler için [Reports](https://nextflow.io/docs/latest/reports.html)'a bakın.

### 5.2. Tüm süreçler için kaynak tahsisleri ayarlama

Profilleme, eğitim iş akışımızdaki süreçlerin çok hafif olduğunu gösteriyor, bu yüzden varsayılan bellek tahsisini süreç başına 1GB'a düşürelim.

`nextflow.config` dosyanıza, pipeline parametreleri bölümünden önce aşağıdakileri ekleyin:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

Bu, tükettiğimiz bilgi işlem miktarını azaltmaya yardımcı olacaktır.

### 5.3. Belirli bir süreç için kaynak tahsisleri ayarlama

Aynı zamanda, `cowpy` sürecinin diğerlerinden daha fazla kaynak gerektirdiğini varsayacağız, sadece bireysel bir süreç için tahsisleri nasıl ayarlayacağımızı gösterebilmemiz için.

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

Bu yapılandırmayla, tüm süreçler 1GB bellek ve tek bir CPU (ima edilen varsayılan) talep edecek, `cowpy` süreci hariç; bu 2GB ve 2 CPU talep edecek.

!!! info

    Az CPU'lu bir makineniz varsa ve süreç başına yüksek bir sayı tahsis ederseniz, süreç çağrılarının birbirinin arkasında sıraya girdiğini görebilirsiniz.
    Bunun nedeni Nextflow'un mevcut olandan daha fazla CPU talep etmememizi sağlamasıdır.

### 5.4. Güncellenmiş yapılandırmayla iş akışını çalıştırma

Hadi deneyelim, yapılandırma değişikliklerinden önce ve sonra performansı karşılaştırabilmemiz için profilleme raporu için farklı bir dosya adı sağlayalım.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Bu çok küçük bir iş yükü olduğu için muhtemelen gerçek bir fark fark etmeyeceksiniz, ancak gerçek dünya iş akışının performansını ve kaynak gereksinimlerini analiz etmek için kullanacağınız yaklaşım budur.

Süreçlerinizin farklı kaynak gereksinimleri olduğunda çok kullanışlıdır. Tahminlere değil, gerçek verilere dayalı olarak her süreç için ayarladığınız kaynak tahsislerini doğru boyutlandırmanızı sağlar.

!!! tip

    Bu, kaynakları kullanımınızı optimize etmek için yapabileceklerinizin sadece küçük bir tadımlığıdır.
    Nextflow'un kendisi, kaynak sınırlamaları nedeniyle başarısız olan işleri yeniden denemek için yerleşik gerçekten düzgün [dinamik yeniden deneme mantığı](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) içerir.
    Ek olarak, Seqera Platform kaynak tahsislerinizi otomatik olarak optimize etmek için YZ destekli araçlar da sunar.

### 5.5. Kaynak sınırları ekleme

Hangi bilgi işlem yürütücüsünü ve bilgi işlem altyapısını kullandığınıza bağlı olarak, tahsis edebileceğiniz (veya tahsis etmeniz gereken) bazı kısıtlamalar olabilir.
Örneğin, kümeniz belirli sınırlar içinde kalmanızı gerektirebilir.

İlgili sınırlamaları ayarlamak için `resourceLimits` yönergesini kullanabilirsiniz. Bir süreç bloğunda tek başına olduğunda sözdizimi şöyle görünür:

```groovy title="Sözdizimi örneği"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow bu değerleri belirttiğiniz yürütücüye bağlı olarak uygun talimatlara çevirecektir.

Eğitim ortamında ilgili altyapıya erişimimiz olmadığı için bunu çalıştırmayacağız.
Ancak, bu sınırları aşan kaynak tahsisleriyle iş akışını çalıştırmayı denerseniz, ardından `.command.run` betik dosyasındaki `sbatch` komutuna bakarsanız, yürütücüye gerçekten gönderilen isteklerin `resourceLimits` tarafından belirtilen değerlerle sınırlandırıldığını görürsünüz.

??? info "Kurumsal referans yapılandırmaları"

    nf-core projesi, dünya çapındaki çeşitli kurumlar tarafından paylaşılan, çok çeşitli HPC ve bulut yürütücülerini kapsayan bir [yapılandırma dosyaları koleksiyonu](https://nf-co.re/configs/) derlemiştir.

    Bu paylaşılan yapılandırmalar hem orada çalışan ve dolayısıyla kurumlarının yapılandırmasını kutudan çıkar çıkmaz kullanabilen insanlar hem de kendi altyapıları için bir yapılandırma geliştirmek isteyen insanlar için bir model olarak değerlidir.

### Özet

Kaynak kullanımını değerlendirmek için bir profilleme raporu nasıl oluşturacağınızı ve tüm süreçler ve/veya bireysel süreçler için kaynak tahsislerini nasıl değiştireceğinizi, ayrıca HPC'de çalıştırmak için kaynak sınırlamalarını nasıl ayarlayacağınızı biliyorsunuz.

### Sırada ne var?

Önceden ayarlanmış yapılandırma profillerini nasıl kuracağınızı ve çalışma zamanında bunlar arasında nasıl geçiş yapacağınızı öğrenin.

---

## 6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanma

??? example "Senaryo"

    Geliştirme için dizüstü bilgisayarınızda ve üretim çalıştırmaları için kurumunuzun HPC'sinde iş akışlarını düzenli olarak çalıştırıyorsunuz.
    Her ortam değiştirdiğinizde yapılandırma ayarlarını manuel olarak değiştirmekten bıktınız.

Size, üzerinde çalıştığınız projeye veya kullandığınız bilgi işlem ortamına bağlı olarak iş akışı yapılandırmanızı özelleştirebileceğiniz bir dizi yol gösterdik.

Hangi bilgi işlem altyapısını kullandığınıza bağlı olarak alternatif ayarlar arasında geçiş yapmak isteyebilirsiniz. Örneğin, dizüstü bilgisayarınızda küçük ölçekli testler geliştirmek ve çalıştırmak, ardından HPC veya bulutta tam ölçekli iş yüklerini çalıştırmak isteyebilirsiniz.

Nextflow, farklı yapılandırmaları tanımlayan istediğiniz sayıda [**profil**](https://nextflow.io/docs/latest/config.html#profiles) kurmanıza olanak tanır; bunları daha sonra yapılandırma dosyasının kendisini değiştirmek zorunda kalmadan bir komut satırı argümanı kullanarak çalışma zamanında seçebilirsiniz.

### 6.1. Yerel geliştirme ve HPC'de yürütme arasında geçiş yapmak için profiller oluşturma

İki alternatif profil kuralım; biri normal bir bilgisayarda küçük ölçekli yükler çalıştırmak için, burada Docker konteynerlerini kullanacağız, ve biri Slurm zamanlayıcılı bir üniversite HPC'sinde çalıştırmak için, burada Conda paketlerini kullanacağız.

#### 6.1.1. Profilleri ayarlama

`nextflow.config` dosyanıza, pipeline parametreleri bölümünden sonra ancak çıktı ayarlarından önce aşağıdakileri ekleyin:

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

Üniversite HPC'si için kaynak sınırlamaları da belirttiğimizi görüyorsunuz.

#### 6.1.2. Bir profille iş akışını çalıştırma

Nextflow komut satırımızda bir profil belirtmek için `-profile` argümanını kullanırız.

`my_laptop` yapılandırmasıyla iş akışını çalıştırmayı deneyelim.

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

Gördüğünüz gibi, bu çalışma zamanında yapılandırmalar arasında çok rahat bir şekilde geçiş yapmamızı sağlar.

!!! warning

    `univ_hpc` profili eğitim ortamında düzgün çalışmayacaktır çünkü bir Slurm zamanlayıcısına erişimimiz yok.

Gelecekte bu yapılandırmaların her zaman birlikte ortaya çıkan başka öğelerini bulursak, bunları ilgili profil(ler)e basitçe ekleyebiliriz.
Birlikte gruplamak istediğimiz başka yapılandırma öğeleri varsa ek profiller de oluşturabiliriz.

### 6.2. Test parametrelerinin bir profilini oluşturma

??? example "Senaryo"

    Başkalarının kendi girdi verilerini toplamak zorunda kalmadan iş akışınızı hızlıca denemesini istiyorsunuz.

Profiller sadece altyapı yapılandırması için değildir.
Bunları iş akışı parametreleri için varsayılan değerler ayarlamak için de kullanabiliriz; böylece başkalarının uygun girdi değerlerini kendileri toplamak zorunda kalmadan iş akışını denemelerini kolaylaştırabiliriz.
Bunu bir parametre dosyası kullanmaya alternatif olarak düşünebilirsiniz.

#### 6.2.1. Profili ayarlama

Bu bağlamda varsayılan değerleri ifade etme sözdizimi, `test` olarak adlandırdığımız bir profil için şöyle görünür:

```groovy title="Sözdizimi örneği"
    test {
        params.<parametre1>
        params.<parametre2>
        ...
    }
```

İş akışımız için bir test profili eklersek, `profiles` bloğu şöyle olur:

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

Tıpkı teknik yapılandırma profilleri gibi, istediğiniz herhangi bir isim altında parametreleri belirten birden fazla farklı profil kurabilirsiniz.

#### 6.2.2. Test profiliyle iş akışını yerel olarak çalıştırma

Uygun bir şekilde, profiller birbirini dışlamaz, bu nedenle komut satırımızda aşağıdaki sözdizimini kullanarak birden fazla profil belirtebiliriz `-profile <profil1>,<profil2>` (herhangi bir sayıda profil için).

Aynı yapılandırma öğeleri için değerler ayarlayan ve aynı yapılandırma dosyasında tanımlanan profilleri birleştirirseniz, Nextflow çatışmayı en son okuduğu değeri kullanarak çözecektir (_yani_ dosyada daha sonra gelen her neyse).
Çakışan ayarlar farklı yapılandırma kaynaklarında ayarlanmışsa, varsayılan [öncelik sırası](https://www.nextflow.io/docs/latest/config.html) geçerlidir.

Önceki komutumuza test profilini eklemeyi deneyelim:

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

Bu, mümkün olduğunda Docker kullanacak ve `results_config/test` altında çıktılar üretecek, ve bu sefer karakter komedi ikilisi `dragonandcow`.

??? abstract "Dosya içeriği"

    ```console title="results_config/test/"
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

Bu, iş akışı koduyla herhangi bir test veri dosyasını dağıttığımız sürece, herkesin komut satırı veya bir parametre dosyası aracılığıyla kendi girdilerini sağlamak zorunda kalmadan iş akışını hızlıca deneyebileceği anlamına gelir.

!!! tip

    Harici olarak depolanan daha büyük dosyalar için URL'lere işaret edebiliriz.
    Açık bir bağlantı olduğu sürece Nextflow bunları otomatik olarak indirecektir.

    Daha fazla ayrıntı için, Yan Görev [Working with Files](../side_quests/working_with_files.md)'a bakın

### 6.3. Çözümlenmiş yapılandırmayı görmek için `nextflow config` kullanma

Yukarıda belirtildiği gibi, bazen aynı parametre birleştirmek istediğiniz profillerde farklı değerlere ayarlanabilir.
Ve daha genel olarak, yapılandırma öğelerinin saklanabileceği çok sayıda yer vardır ve bazen aynı özellikler farklı yerlerde farklı değerlere ayarlanabilir.

Nextflow herhangi bir çatışmayı çözmek için bir [öncelik sırası](https://nextflow.io/docs/latest/config.html#configuration-file) uygular, ancak bunu kendiniz belirlemek zor olabilir.
Ve hiçbir şey çakışmasa bile, şeylerin yapılandırılabileceği tüm olası yerlere bakmak sıkıcı olabilir.

Neyse ki, Nextflow tüm bu süreci sizin için otomatikleştirebilecek `config` adlı kullanışlı bir yardımcı araç içerir.

`config` aracı mevcut çalışma dizininizdeki tüm içerikleri keşfedecek, herhangi bir yapılandırma dosyasını toplayacak ve Nextflow'un iş akışını çalıştırmak için kullanacağı tam olarak çözümlenmiş yapılandırmayı üretecektir.
Bu, herhangi bir şey başlatmak zorunda kalmadan hangi ayarların kullanılacağını öğrenmenizi sağlar.

#### 6.3.1. Varsayılan yapılandırmayı çözme

Varsayılan olarak uygulanacak yapılandırmayı çözmek için bu komutu çalıştırın.

```bash
nextflow config
```

??? success "Komut çıktısı"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

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

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Bu, komut satırında ekstra bir şey belirtmezseniz aldığınız temel yapılandırmayı gösterir.

#### 6.3.2. Belirli ayarlar etkinleştirilmiş yapılandırmayı çözme

Komut satırı parametreleri sağlarsanız, örneğin bir veya daha fazla profili etkinleştirerek veya bir parametre dosyası yükleyerek, komut bunları da ek olarak dikkate alacaktır.

```bash
nextflow config -profile my_laptop,test
```

??? success "Komut çıktısı"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

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

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Bu, birden fazla yapılandırma katmanı içeren karmaşık projeler için özellikle kullanışlı hale gelir.

### Özet

Çalışma zamanında minimum güçlükle önceden ayarlanmış bir yapılandırmayı seçmek için profilleri nasıl kullanacağınızı biliyorsunuz.
Daha genel olarak, farklı bilgi işlem platformlarına uyacak ve analizlerinizin tekrarlanabilirliğini artıracak şekilde iş akışı yürütmelerinizi nasıl yapılandıracağınızı biliyorsunuz.

### Sırada ne var?

GitHub gibi uzak depolardan doğrudan iş akışlarını nasıl çalıştıracağınızı öğrenin.

---

## 7. Uzak depolardan iş akışlarını çalıştırma

??? example "Senaryo"

    nf-core'dan olanlar gibi köklü bir iş akışını kodu kendiniz indirip yönetmek zorunda kalmadan çalıştırmak istiyorsunuz.

Şimdiye kadar mevcut dizinde bulunan iş akışı betiklerini çalıştırıyorduk.
Pratikte, genellikle GitHub gibi uzak depolarda saklanan iş akışlarını çalıştırmak isteyeceksiniz.

Nextflow bunu basitleştirir: herhangi bir iş akışını manuel olarak indirmeden önce doğrudan bir Git deposu URL'sinden çalıştırabilirsiniz.

### 7.1. GitHub'dan bir iş akışı çalıştırma

Uzak bir iş akışını çalıştırmanın temel sözdizimi `nextflow run <depo>`'dur; burada `<depo>` `nextflow-io/hello` gibi bir GitHub depo yolu, tam bir URL veya GitLab, Bitbucket veya diğer Git barındırma hizmetlerine bir yol olabilir.

Resmi Nextflow "hello" demo iş akışını çalıştırmayı deneyin:

```bash
nextflow run nextflow-io/hello
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Uzak bir iş akışını ilk kez çalıştırdığınızda, Nextflow onu indirir ve yerel olarak önbelleğe alır.
Sonraki çalıştırmalar, açıkça bir güncelleme talep etmediğiniz sürece önbelleğe alınmış sürümü kullanır.

### 7.2. Tekrarlanabilirlik için bir sürüm belirtme

Varsayılan olarak, Nextflow varsayılan daldan en son sürümü çalıştırır.
`-r` bayrağını kullanarak belirli bir sürümü (etiket), dalı veya commit'i belirtebilirsiniz:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Tam sürümleri belirtmek tekrarlanabilirlik için esastır.

### Özet

GitHub ve diğer uzak depolardan doğrudan iş akışlarını nasıl çalıştıracağınızı ve tekrarlanabilirlik için sürümleri nasıl belirteceğinizi biliyorsunuz.

### Sırada ne var?

Kendinize büyük bir alkış verin!
Nextflow iş akışlarını çalıştırmaya ve yönetmeye başlamak için bilmeniz gereken her şeyi biliyorsunuz.

Bu kurs burada sona eriyor, ancak öğrenmeye devam etmeye hevesliyseniz, iki ana önerimiz var:

- Kendi iş akışlarınızı geliştirmeye daha derinlemesine inmek istiyorsanız, kanallar ve operatörler hakkında aynı genel ilerlemeyi kapsayan ancak çok daha ayrıntılı bir yeni başlayanlar kursu olan [Hello Nextflow](../hello_nextflow/index.md)'a göz atın.
- Koda daha derinlemesine inmeden Nextflow iş akışlarını çalıştırmayı öğrenmeye devam etmek istiyorsanız, son derece popüler [nf-core](https://nf-co.re/) projesinden iş akışlarını bulmak ve çalıştırmak için araçları tanıtan [Hello nf-core](../hello_nf-core/index.md)'un ilk bölümüne göz atın.

İyi eğlenceler!

---

## Quiz

<quiz>
Parametre değerleri hem iş akışı dosyasında hem de `nextflow.config`'de ayarlandığında, hangisi önceliklidir?
- [ ] İş akışı dosyası değeri
- [x] Yapılandırma dosyası değeri
- [ ] Karşılaşılan ilk değer
- [ ] Bir hataya neden olur

Daha fazla bilgi: [1.1. `nextflow.config` dosyasında değerleri ayarlama](#11-nextflowconfig-dosyasında-değerleri-ayarlama)
</quiz>

<quiz>
Bir iş akışı dosyasında parametre varsayılanı ayarlama ile bir yapılandırma dosyasında ayarlama arasındaki sözdizimi farkı nedir?
- [ ] Aynı sözdizimini kullanırlar
- [x] İş akışı türü belirtilmiş bildirim (`#!groovy param: Type = value`) kullanır, yapılandırma atama (`#!groovy param = value`) kullanır
- [ ] Yapılandırma türü belirtilmiş bildirim kullanır, iş akışı atama kullanır
- [ ] Sadece yapılandırma dosyaları varsayılan değerler ayarlayabilir

Daha fazla bilgi: [1.1. `nextflow.config` dosyasında değerleri ayarlama](#11-nextflowconfig-dosyasında-değerleri-ayarlama)
</quiz>

<quiz>
Bir iş akışını çalıştırırken bir parametre dosyasını nasıl belirtirsiniz?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Daha fazla bilgi: [1.3. Bir parametre dosyası kullanma](#13-bir-parametre-dosyası-kullanma)
</quiz>

<quiz>
`outputDir` yapılandırma seçeneği neyi kontrol eder?
- [ ] Çalışma dizininin konumu
- [x] İş akışı çıktılarının yayınlandığı temel yol
- [ ] Günlük dosyaları için dizin
- [ ] Modül dosyalarının konumu

Daha fazla bilgi: [2.1. `outputDir` dizin adını özelleştirme](#21-outputdir-dizin-adını-özelleştirme)
</quiz>

<quiz>
Çıktı yolu yapılandırmasında bir süreç adına dinamik olarak nasıl referans verirsiniz?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

Daha fazla bilgi: [2.2. Çıktıları sürece göre düzenleme](#22-çıktıları-sürece-göre-düzenleme)
</quiz>

<quiz>
Hem Docker hem de Conda etkinleştirilmişse ve bir sürecin her iki yönergesi de varsa, hangisine öncelik verilir?
- [x] Docker (konteynerler)
- [ ] Conda
- [ ] Süreçte tanımlanan ilki
- [ ] Bir hataya neden olur

Daha fazla bilgi: [3. Bir yazılım paketleme teknolojisi seçme](#3-bir-yazılım-paketleme-teknolojisi-seçme)
</quiz>

<quiz>
Nextflow'da varsayılan yürütücü nedir?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Daha fazla bilgi: [4. Bir yürütme platformu seçme](#4-bir-yürütme-platformu-seçme)
</quiz>

<quiz>
Hangi komut bir kaynak kullanım raporu oluşturur?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Daha fazla bilgi: [5.1. Bir kaynak kullanım raporu oluşturmak için iş akışını çalıştırma](#51-bir-kaynak-kullanım-raporu-oluşturmak-için-iş-akışını-çalıştırma)
</quiz>

<quiz>
Yapılandırma dosyasında `cowpy` adlı belirli bir süreç için kaynak gereksinimlerini nasıl ayarlarsınız?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Daha fazla bilgi: [5.3. Belirli bir süreç için kaynak tahsisleri ayarlama](#53-belirli-bir-süreç-için-kaynak-tahsisleri-ayarlama)
</quiz>

<quiz>
`resourceLimits` yönergesi ne yapar?
- [ ] Minimum kaynak gereksinimlerini ayarlar
- [ ] Süreçlere kaynak tahsis eder
- [x] Talep edilebilecek maksimum kaynakları sınırlar
- [ ] Kaynak kullanımını gerçek zamanlı olarak izler

Daha fazla bilgi: [5.5. Kaynak sınırları ekleme](#55-kaynak-sınırları-ekleme)
</quiz>

<quiz>
Tek bir komutta birden fazla profili nasıl belirtirsiniz?
- [ ] `-profile profil1 -profile profil2`
- [ ] `-profiles profil1,profil2`
- [x] `-profile profil1,profil2`
- [ ] `--profile profil1 --profile profil2`

Daha fazla bilgi: [6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanma](#6-önceden-ayarlanmış-yapılandırmalar-arasında-geçiş-yapmak-için-profilleri-kullanma)
</quiz>

<quiz>
Nextflow'un kullanacağı tam olarak çözümlenmiş yapılandırmayı hangi komut gösterir?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Daha fazla bilgi: [6.3. Çözümlenmiş yapılandırmayı görmek için `nextflow config` kullanma](#63-çözümlenmiş-yapılandırmayı-görmek-için-nextflow-config-kullanma)
</quiz>

<quiz>
Profiller ne için kullanılabilir? (Geçerli olanların tümünü seçin)
- [x] Altyapıya özgü ayarları tanımlama (yürütücüler, konteynerler)
- [x] Farklı ortamlar için kaynak sınırları ayarlama
- [x] Kolay iş akışı testi için test parametreleri sağlama
- [ ] Yeni süreçler tanımlama

Daha fazla bilgi: [6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanma](#6-önceden-ayarlanmış-yapılandırmalar-arasında-geçiş-yapmak-için-profilleri-kullanma)
</quiz>
