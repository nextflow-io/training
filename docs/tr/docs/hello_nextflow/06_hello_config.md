# Bölüm 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube kanalında görün.

:green_book: Video metni [burada](./transcripts/06_hello_config.md) mevcuttur.
///

Bu bölümde, Nextflow pipeline'ınızın yapılandırmasını nasıl kuracağınızı ve yöneteceğinizi keşfedeceğiz; böylece _iş akışı kodunun kendisinde tek bir satır bile değiştirmeden_ davranışını özelleştirebilecek, farklı ortamlara uyarlayabilecek ve kaynak kullanımını optimize edebileceksiniz.

Bunu yapmanın birden fazla yolu vardır; bunlar birlikte kullanılabilir ve yapılandırma belgelerinde açıklanan [öncelik sırasına](https://nextflow.io/docs/latest/config.html) göre yorumlanır.

Kursun bu bölümünde, size en basit ve en yaygın yapılandırma dosyası mekanizmasını göstereceğiz: Bölüm 5: Hello Containers'da zaten karşılaştığınız [`nextflow.config`](https://nextflow.io/docs/latest/config.html) dosyası.

Süreç yönergeleri, yürütücüler, profiller ve parametre dosyaları gibi Nextflow yapılandırmasının temel bileşenlerini ele alacağız.
Bu yapılandırma seçeneklerini etkili bir şekilde kullanmayı öğrenerek, pipeline'larınızın esnekliğini, ölçeklenebilirliğini ve performansını artırabilirsiniz.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Hello Nextflow](./index.md) kursunun 1-5. Bölümlerini tamamladığınızı ve eksiksiz çalışan bir pipeline'a sahip olduğunuzu varsayar.

    Kursa bu noktadan başlıyorsanız, `modules` dizinini ve `nextflow.config` dosyasını çözümlerden kopyalamanız gerekecektir:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    `nextflow.config` dosyası, Docker konteynerlerinin kullanımını etkinleştiren `docker.enabled = true` satırını içerir.

    Hello pipeline'ına aşina değilseniz veya hatırlatmaya ihtiyacınız varsa, [bu bilgi sayfasına](../info/hello_pipeline.md) bakın.

---

## 0. Isınma: `hello-config.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-config.nf` iş akışı betiğini kullanacağız.
Bu eğitim kursunun 5. Bölümünde üretilen betiğe eşdeğerdir, ancak çıktı hedeflerini değiştirdik:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Her şeyin çalıştığından emin olmak için, herhangi bir değişiklik yapmadan önce betiği bir kez çalıştırın:

```bash
nextflow run hello-config.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Daha önce olduğu gibi, çıktı dosyalarını `output` bloğunda belirtilen dizinde (`results/hello_config/`) bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

Son ASCII sanat çıktısı, `results/hello_config/` dizininde `cowpy-COLLECTED-batch-output.txt` adı altındadır.

??? abstract "Dosya içeriği"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Bu sizin için çalıştıysa, pipeline'larınızı nasıl yapılandıracağınızı öğrenmeye hazırsınız.

---

## 1. İş akışı girdi parametrelerini yönetin

İş akışı betiğinin kendisinde varsayılan değerler ayarlanmış olarak, komut satırı üzerinden birkaç parametre değerini kabul edecek şekilde kurulmuş olan bir yapılandırma unsuruyla başlayacağız.
Ancak, bu varsayılanları komut satırında parametreler belirtmek veya orijinal betik dosyasını değiştirmek zorunda kalmadan geçersiz kılmak isteyebilirsiniz.

Bunu yapmanın birden fazla yolu vardır; size çok yaygın olarak kullanılan üç temel yolu göstereceğiz.

### 1.1. Varsayılan değerleri `nextflow.config` dosyasına taşıyın

Bu en basit yaklaşımdır, ancak ana `nextflow.config` dosyası her çalıştırma için düzenlemek isteyeceğiniz bir şey olmadığından muhtemelen en az esnektir.
Ancak, iş akışında parametreleri _bildirme_ (ki bu kesinlikle oraya aittir) ile _varsayılan değerler_ sağlama (bunlar bir yapılandırma dosyasında daha uygun) kaygılarını ayırma avantajına sahiptir.

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

`params` bloğunu iş akışından yapılandırma dosyasına basitçe kopyalamadığımıza dikkat edin.
Sözdizimi biraz farklıdır.
İş akışı dosyasında bunlar tipli bildirimlerdir.
Yapılandırmada bunlar değer atamalarıdır.

Teknik olarak, bu iş akışı dosyasında hala belirtilen varsayılan değerleri geçersiz kılmak için yeterlidir.
Örneğin karakteri değiştirebilir ve yapılandırma dosyasında ayarlanan değerin iş akışı dosyasında ayarlanan değeri geçersiz kıldığından emin olmak için iş akışını çalıştırabilirsiniz.

Ancak yapılandırmayı tamamen yapılandırma dosyasına taşıma ruhuyla, bu değerleri iş akışı dosyasından tamamen kaldıralım.

#### 1.1.2. İş akışı dosyasındaki `params` bloğundan değerleri kaldırın

`hello-config.nf` iş akışı dosyasında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
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

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parametreleri
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Artık iş akışı dosyasının kendisi bu parametreler için herhangi bir varsayılan değer ayarlamıyor.

#### 1.1.3. Pipeline'ı çalıştırın

Doğru çalıştığını test edelim.

```bash
nextflow run hello-config.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretir.

Son ASCII sanat çıktısı, daha önce olduğu gibi `results/hello_config/` dizininde `cowpy-COLLECTED-batch-output.txt` adı altındadır.

??? abstract "Dosya içeriği"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

İşlevsel olarak, bu taşıma hiçbir şeyi değiştirmedi, ancak kavramsal olarak varsayılan değerlerin yapılandırma dosyasında ayarlanması biraz daha temiz.

### 1.2. Çalıştırmaya özgü bir yapılandırma dosyası kullanın

Bu harika, ancak bazen ana yapılandırma dosyasını karıştırmadan farklı varsayılan değerlerle bazı geçici deneyler yapmak isteyebilirsiniz.
Bunu, deneyler için çalışma dizini olarak kullanacağınız bir alt dizinde yeni bir `nextflow.config` dosyası oluşturarak yapabilirsiniz.

#### 1.2.1. Boş bir yapılandırmayla çalışma dizini oluşturun

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

#### 1.2.2. Deneysel yapılandırmayı kurun

Şimdi yeni dosyayı açın ve özelleştirmek istediğiniz parametreleri ekleyin:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Girdi dosyasının yolunun dizin yapısını yansıtması gerektiğine dikkat edin.

#### 1.2.3. Pipeline'ı çalıştırın

Artık pipeline'ımızı yeni çalışma dizinimizin içinden çalıştırabiliriz.
Yolu buna göre uyarladığınızdan emin olun!

```bash
nextflow run ../hello-config.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Bu, `tux-run/work/` ve `tux-run/results/` dahil olmak üzere `tux-run/` altında yeni bir dizin seti oluşturacaktır.

Bu çalıştırmada, Nextflow mevcut dizinimizdeki `nextflow.config` dosyasını pipeline'ın kök dizinindeki `nextflow.config` ile birleştirir ve böylece varsayılan karakteri (turkey) tux karakteriyle geçersiz kılar.

Son çıktı dosyası, selamlamaları söyleyen tux karakterini içermelidir.

??? abstract "Dosya içeriği"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

İşte bu kadar; artık 'normal' yapılandırmanızı değiştirmeden deney yapmak için bir alanınız var.

!!! warning "Uyarı"

    Bir sonraki bölüme geçmeden önce önceki dizine geri döndüğünüzden emin olun!

    ```bash
    cd ..
    ```

Şimdi parametre değerlerini ayarlamanın başka bir kullanışlı yoluna bakalım.

### 1.3. Bir parametre dosyası kullanın

Alt dizin yaklaşımı deney yapmak için harika çalışır, ancak biraz kurulum gerektirir ve yolları buna göre uyarlamanızı gerektirir.
Pipeline'ınızı belirli bir değer seti ile çalıştırmak istediğinizde veya başka birinin bunu minimum çabayla yapmasını sağlamak istediğinizde daha basit bir yaklaşım vardır.

Nextflow, YAML veya JSON formatında bir [parametre dosyası](https://nextflow.io/docs/latest/config.html#params-file) aracılığıyla parametreleri belirtmemize olanak tanır; bu da alternatif varsayılan değer setlerini yönetmeyi ve dağıtmayı, örneğin çalıştırmaya özgü parametre değerlerini çok kullanışlı hale getirir.

#### 1.3.1. Örnek parametre dosyasını inceleyin

Bunu göstermek için, mevcut dizinde `test-params.yaml` adlı bir örnek parametre dosyası sağlıyoruz:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Bu parametre dosyası, belirtmek istediğimiz her girdi için bir anahtar-değer çifti içerir.
Sözdizimini yapılandırma dosyasıyla karşılaştırırsanız, eşittir işaretleri (`=`) yerine iki nokta üst üste (`:`) kullanımına dikkat edin.
Yapılandırma dosyası Groovy'de yazılırken, parametre dosyası YAML'de yazılır.

!!! info "Bilgi"

    Ayrıca örnek olarak parametre dosyasının JSON sürümünü de sağlıyoruz ancak burada bununla çalıştırmayacağız.
    Bunu kendi başınıza denemekten çekinmeyin.

#### 1.3.2. Pipeline'ı çalıştırın

İş akışını bu parametre dosyasıyla çalıştırmak için, temel komuta basitçe `-params-file <dosyaadı>` ekleyin.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Son çıktı dosyası, selamlamaları söyleyen stegosaurus karakterini içermelidir.

??? abstract "Dosya içeriği"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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

Yalnızca birkaç parametre belirtmeniz gerektiğinde bir parametre dosyası kullanmak aşırıya kaçmış gibi görünebilir, ancak bazı pipeline'lar düzinelerce parametre bekler.
Bu durumlarda, bir parametre dosyası kullanmak, çalışma zamanında devasa komut satırları yazmak ve iş akışı betiğini değiştirmek zorunda kalmadan parametre değerleri sağlamamıza olanak tanır.

Ayrıca parametre setlerini işbirlikçilere dağıtmayı veya örneğin bir yayın için destekleyici bilgi olarak sunmayı kolaylaştırır.
Bu, çalışmanızı başkaları tarafından daha tekrarlanabilir hale getirir.

### Özet

İş akışı girdilerini yönetmek için temel yapılandırma seçeneklerinden nasıl yararlanacağınızı biliyorsunuz.

### Sırada ne var?

İş akışı çıktılarınızın nerede ve nasıl yayınlanacağını nasıl yöneteceğinizi öğrenin.

---

## 2. İş akışı çıktılarını yönetin

Şimdiye kadar iş akışı düzeyindeki çıktı bildirimleri için tüm yolları sabit kodluyorduk ve birden fazla çıktı eklemeye başladığımızda belirttiğimiz gibi, biraz tekrar söz konusu olabilir.

Bunu daha esnek olacak şekilde yapılandırabileceğiniz birkaç yaygın yola bakalım.

### 2.1. Çıktı dizinini `-output-dir` ile özelleştirin

'Yayınlanan' çıktılarımızın nasıl organize edildiğini kontrol ederken iki farklı önceliğimiz var:

- Üst düzey çıktı dizini
- Bu dizin içinde dosyaların nasıl organize edildiği

Şimdiye kadar varsayılan üst düzey dizini kullanıyorduk: `results`.
`-output-dir` CLI seçeneğini kullanarak bunu özelleştirerek başlayalım.

#### 2.1.1. Pipeline'ı `-output-dir` ile çalıştırın

`-output-dir` seçeneği (kısaltma: `-o`) tüm iş akışı çıktıları için varsayılan çıktı dizinini (`results/`) geçersiz kılar.
Bu, çıktıların yayınlandığı kök yolu kontrol etmenin önerilen yoludur.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

Bu, çıktıları `results/` yerine `custom-outdir-cli/` dizinine yayınlar:

??? abstract "Dizin içeriği"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Çıktı bloğundaki `path` bildirimlerinden hala `hello_config` alt dizinine sahip olduğumuza dikkat edin.
Bunu temizleyelim.

#### 2.1.2. Çıktı bloğundan sabit kodlanmış yolları kaldırın

`hello_config/` öneki önceki bölümlerde sabit kodlanmıştı, ancak artık çıktı yollarını esnek bir şekilde yapılandırmayı öğrendiğimize göre, bu sabit kodlamayı kaldırabiliriz.
Alt dizine ihtiyaç duymayan çıktılar için `path` yönergesini boş bir dizeye ayarlayabilir veya tamamen kaldırabiliriz.

İş akışı dosyasında aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Pipeline'ı tekrar çalıştırın:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Artık çıktılar `hello_config` alt dizini olmadan doğrudan `custom-outdir-cli-2/` altında yayınlanıyor:

??? abstract "Dizin içeriği"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip "İpucu"

    `-output-dir` seçeneği çıktıların _nereye_ gideceğini kontrol etmek için kullanılırken, çıktı bloğundaki `path` yönergesi _alt dizin yapısını_ kontrol eder.

### 2.2. Dinamik çıktı yolları

CLI üzerinden çıktı dizinini değiştirmenin yanı sıra, `outputDir` kullanarak yapılandırma dosyasında özel bir varsayılan değer de ayarlayabiliriz.
Bu, dizin yolunu dinamik olarak ayarlamamıza olanak tanır - sadece statik dizeler kullanmakla kalmaz.

#### 2.2.1. Yapılandırma dosyasında `outputDir` ayarlayın

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
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/${params.batch}"
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

Bu, çıktı dizinini `custom-outdir-config/` artı `batch` parametresinin değerini alt dizin olarak ayarlar.
Artık `--batch` parametresini ayarlayarak çıktı konumunu değiştirebilirsiniz:

```bash
nextflow run hello-config.nf --batch my_run
```

Bu, çıktıları `custom-outdir-config/my_run/` dizinine yayınlar.

!!! note "Not"

    `-output-dir` CLI seçeneği, `outputDir` yapılandırma ayarından önceliklidir.
    Ayarlanırsa, yapılandırma seçeneği tamamen göz ardı edilecektir.

#### 2.2.2. Batch ve süreç adlarıyla alt dizinler

Ayrıca çıktı başına alt dizin çıktı `path` bildirimlerini dinamik olarak ayarlayabiliriz.

Örneğin, çıktı yolu bildiriminde `<process>.name` referansı vererek çıktılarımızı sürece göre organize edebiliriz:

=== "Sonra"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Daha da ileri gidebilir ve daha karmaşık alt dizin yolları oluşturabiliriz.

Yukarıdaki düzenlemede `intermediates` ile üst düzeydeki son çıktılar arasındaki ayrımı sildik.
Bunu geri koyalım ve ayrıca dosyaları bir `params.batch` alt dizinine koyalım.

!!! tip "İpucu"

    `params.batch` değerini `outputDir` yapılandırması yerine çıktı bloğu `path` içine dahil etmek, CLI'da `-output-dir` ile üzerine yazılmayacağı anlamına gelir.

İlk olarak, `outputDir` değerinden `${params.batch}` değerini kaldırmak için yapılandırma dosyasını güncelleyin (çünkü bunu yol bildirimlerine taşıyoruz):

=== "Sonra"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Ardından, iş akışı dosyasında aşağıdaki değişiklikleri yapın:

=== "Sonra"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. Pipeline'ı çalıştırın

Bunun pratikte nasıl çalıştığını görelim, hem `-output-dir` (veya kısaca `-o`) değerini `custom-outdir-config-2` olarak hem de batch adını komut satırından `rep2` olarak ayarlayalım:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

Bu, çıktıları belirtilen temel yol _ve_ batch adı alt dizini _ve_ sürece göre gruplandırılmış sonuçlarla `custom-outdir-config-2/rep2/` dizinine yayınlar:

??? abstract "Dizin içeriği"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Yayınlama modunu iş akışı düzeyinde ayarlayın

Son olarak, tekrarlayan kod miktarını azaltma ruhuyla, çıktı başına `mode` bildirimlerini yapılandırmada tek bir satırla değiştirebiliriz.

#### 2.3.1. Yapılandırma dosyasına `workflow.output.mode` ekleyin

`nextflow.config` dosyasına aşağıdaki kodu ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/"
    ```

Yapılandırma dosyasında `workflow.output.mode` ayarlamak, iş akışı dosyasında ayarlananı geçersiz kılmak için yeterlidir, ancak gereksiz kodu yine de kaldıralım.

#### 2.3.2. İş akışı dosyasından çıktı modunu kaldırın

İş akışı dosyasında aşağıdaki değişiklikleri yapın:

=== "Sonra"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

Bu daha özlü, değil mi?

#### 2.3.3. Pipeline'ı çalıştırın

Doğru çalıştığını test edelim:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

Bu, çıktıları `config-output-mode/` dizinine yayınlar ve hepsi hala sembolik bağlantılar değil, uygun kopyalardır.

??? abstract "Dizin içeriği"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

Çıktı başına mod ayarlama yöntemini kullanmak isteyebileceğiniz ana neden, aynı iş akışı içinde karıştırmak ve eşleştirmek istemenizdir, _yani_ bazı çıktıların kopyalanmasını ve bazılarının sembolik bağlantı olmasını sağlamak.

Bu şekilde özelleştirebileceğiniz başka birçok seçenek vardır, ancak umarım bu size seçenek yelpazesi ve tercihlerinize uygun şekilde bunları etkili bir şekilde nasıl kullanacağınız konusunda bir fikir verir.

### Özet

Çıktılarınızın yayınlandığı dizinlerin adlandırmasını ve yapısını, ayrıca iş akışı çıktı yayınlama modunu nasıl kontrol edeceğinizi biliyorsunuz.

### Sırada ne var?

İş akışı yapılandırmanızı bilgi işlem ortamınıza nasıl uyarlayacağınızı, yazılım paketleme teknolojisiyle başlayarak öğrenin.

---

## 3. Bir yazılım paketleme teknolojisi seçin

Şimdiye kadar girdilerin nasıl girdiğini ve girdilerin nereye çıktığını kontrol eden yapılandırma öğelerine bakıyorduk. Şimdi iş akışı yapılandırmanızı bilgi işlem ortamınıza uyarlamaya daha spesifik olarak odaklanma zamanı.

Bu yoldaki ilk adım, her adımda çalıştırılacak yazılım paketlerinin nereden geleceğini belirtmektir.
Bunlar yerel bilgi işlem ortamında zaten yüklü mü?
Görüntüleri almamız ve bunları bir konteyner sistemi aracılığıyla çalıştırmamız mı gerekiyor?
Yoksa Conda paketlerini almamız ve yerel bir Conda ortamı oluşturmamız mı gerekiyor?

Bu eğitim kursunun ilk bölümünde (Bölüm 1-4) iş akışımızda sadece yerel olarak yüklenmiş yazılımı kullandık.
Ardından Bölüm 5'te Docker konteynerlerini ve Docker konteynerlerinin kullanımını etkinleştirmek için kullandığımız `nextflow.config` dosyasını tanıttık.

Şimdi `nextflow.config` dosyası aracılığıyla alternatif bir yazılım paketleme seçeneğini nasıl yapılandırabileceğimizi görelim.

### 3.1. Yapılandırma dosyasında Docker'ı devre dışı bırakın ve Conda'yı etkinleştirin

Bir HPC kümesinde çalıştığımızı ve yöneticinin güvenlik nedenleriyle Docker kullanımına izin vermediğini varsayalım.
Neyse ki bizim için Nextflow, Singularity (HPC'de daha yaygın olarak kullanılır) dahil olmak üzere birden fazla başka konteyner teknolojisini ve Conda gibi yazılım paket yöneticilerini destekler.

Yapılandırma dosyamızı Docker yerine [Conda](https://nextflow.io/docs/latest/conda.html) kullanacak şekilde değiştirebiliriz.
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
Bu da artık `cowpy` sürecimize bunlardan birini eklememiz gerektiği anlamına gelir!

### 3.2. Süreç tanımında bir Conda paketi belirtin

`cowpy` aracını içeren bir Conda paketi için URI'yi zaten aldık: `conda-forge::cowpy==1.1.5`

Şimdi URI'yi `conda` yönergesini kullanarak `cowpy` süreç tanımına ekleyin:

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

!!! tip "İpucu"

    Belirli bir conda paketi için URI'yi almanın birkaç farklı yolu vardır.
    Bir konteyner oluşturmayı planlamasanız bile kopyalayıp yapıştırabileceğiniz bir URI verecek olan [Seqera Containers](https://seqera.io/containers/) arama sorgusunu kullanmanızı öneririz.

### 3.3. Conda kullanabileceğini doğrulamak için iş akışını çalıştırın

Hadi deneyelim.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Komut çıktısı"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

Bu sorunsuz çalışmalı ve `custom-outdir-config/conda` altında daha önce olduğu gibi aynı çıktıları üretmelidir.

Perde arkasında, Nextflow Conda paketlerini almış ve ortamı oluşturmuştur; bu normalde biraz iş gerektirir; bu yüzden bunların hiçbirini kendimiz yapmak zorunda kalmamamız güzel!

!!! note "Not"

    `cowpy` paketi oldukça küçük olduğu için bu hızlı çalışır, ancak büyük paketlerle çalışıyorsanız, ilk seferinde normalden biraz daha uzun sürebilir ve konsol çıktısının tamamlanmadan önce bir dakika kadar 'takılı' kaldığını görebilirsiniz.
    Bu normaldir ve Nextflow'un yeni bir paketi ilk kez kullandığınızda yaptığı ekstra işten kaynaklanır.

Bizim açımızdan, arka planda mekanikler biraz farklı olsa da Docker ile çalıştırmakla tamamen aynı şekilde çalışıyor gibi görünüyor.

Bu, gerekirse Conda ortamlarıyla çalıştırmaya hazır olduğumuz anlamına gelir.

??? info "Docker ve Conda'yı karıştırma ve eşleştirme"

    Bu yönergeler süreç başına atandığından, örneğin kullandığınız bilgi işlem altyapısı her ikisini de destekliyorsa, iş akışınızdaki bazı süreçleri Docker ile ve diğerlerini Conda ile çalıştıracak şekilde yapılandırmak, yani 'karıştırmak ve eşleştirmek' mümkündür.
    Bu durumda, yapılandırma dosyanızda hem Docker'ı hem de Conda'yı etkinleştirirsiniz.
    Belirli bir süreç için her ikisi de mevcutsa, Nextflow konteynerlere öncelik verecektir.

    Ve daha önce belirtildiği gibi, Nextflow birden fazla başka yazılım paketleme ve konteyner teknolojisini destekler, bu nedenle sadece bu ikisiyle sınırlı değilsiniz.

### Özet

Her sürecin hangi yazılım paketini kullanması gerektiğini nasıl yapılandıracağınızı ve teknolojiler arasında nasıl geçiş yapacağınızı biliyorsunuz.

### Sırada ne var?

Nextflow tarafından işi gerçekten yapmak için kullanılan yürütme platformunu nasıl değiştireceğinizi öğrenin.

---

## 4. Bir yürütme platformu seçin

Şimdiye kadar pipeline'ımızı yerel yürütücü ile çalıştırıyorduk.
Bu, her görevi Nextflow'un üzerinde çalıştığı makinede yürütür.
Nextflow başladığında, mevcut CPU'lara ve belleğe bakar.
Çalıştırılmaya hazır görevlerin kaynakları mevcut kaynakları aşarsa, Nextflow son görevleri, önceki görevlerden biri veya daha fazlası bitene ve gerekli kaynakları serbest bırakana kadar yürütmeden geri tutar.

Yerel yürütücü kullanışlı ve verimlidir, ancak tek bir makineyle sınırlıdır. Çok büyük iş yükleri için, yerel makinenizin bir darboğaz olduğunu keşfedebilirsiniz; ya sahip olduğunuzdan daha fazla kaynak gerektiren tek bir göreviniz olduğu için ya da tek bir makinenin bunları çalıştırmasını beklemenin çok uzun süreceği kadar çok göreviniz olduğu için.

Nextflow, HPC zamanlayıcıları (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor ve diğerleri) ve bulut yürütme arka uçları (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes ve daha fazlası) dahil olmak üzere [birçok farklı yürütücüyü](https://nextflow.io/docs/latest/executor.html) destekler.

### 4.1. Farklı bir arka ucu hedefleme

Yürütücü seçimi, `executor` adlı bir süreç yönergesi tarafından ayarlanır.
Varsayılan olarak `local` olarak ayarlanmıştır, bu nedenle aşağıdaki yapılandırma ima edilir:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Yürütücüyü farklı bir arka ucu hedefleyecek şekilde ayarlamak için, kaynak tahsisleri için yukarıda açıklanan benzer sözdizimini kullanarak istediğiniz yürütücüyü belirtirsiniz (tüm seçenekler için [yürütücü belgelerine](https://nextflow.io/docs/latest/executor.html) bakın).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Uyarı"

    Eğitim ortamı bir HPC'ye bağlanacak şekilde kurulmadığı için bunu gerçekten test edemeyiz.

### 4.2. Yürütme parametreleri için arka uca özgü sözdizimi ile başa çıkma

Çoğu yüksek performanslı bilgi işlem platformu, CPU sayısı ve bellek gibi kaynak tahsis isteklerini ve sınırlamalarını ve kullanılacak iş kuyruğunun adını belirtmenize izin verir (ve bazen gerektirir).

Ne yazık ki, bu sistemlerin her biri bir işin nasıl tanımlanması ve ilgili zamanlayıcıya gönderilmesi gerektiğini tanımlamak için farklı teknolojiler, sözdizimi ve yapılandırmalar kullanır.

??? abstract "Örnekler"

    Örneğin, "my-science-work" kuyruğunda yürütülmek üzere 8 CPU ve 4GB RAM gerektiren aynı işin arka uca bağlı olarak aşağıdaki farklı şekillerde ifade edilmesi gerekir.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Neyse ki, Nextflow tüm bunları basitleştirir.
[`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) ve [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) gibi ilgili özellikleri belirtebilmeniz için standartlaştırılmış bir sözdizimi sağlar (diğer özellikler için [süreç yönergelerine](https://nextflow.io/docs/latest/reference/process.html#process-directives) bakın) sadece bir kez.
Ardından, çalışma zamanında Nextflow, yürütücü ayarına göre uygun arka uca özgü betikleri oluşturmak için bu ayarları kullanacaktır.

Bu standartlaştırılmış sözdizimini bir sonraki bölümde ele alacağız.

### Özet

Artık farklı türde bilgi işlem altyapılarını kullanmak için yürütücüyü nasıl değiştireceğinizi biliyorsunuz.

### Sırada ne var?

Nextflow'da kaynak tahsislerini ve sınırlamalarını nasıl değerlendireceğinizi ve ifade edeceğinizi öğrenin.

---

## 5. Bilgi işlem kaynak tahsislerini kontrol edin

Çoğu yüksek performanslı bilgi işlem platformu, CPU sayısı ve bellek gibi belirli kaynak tahsis parametrelerini belirtmenize izin verir (ve bazen gerektirir).

Varsayılan olarak, Nextflow her süreç için tek bir CPU ve 2GB bellek kullanacaktır.
Karşılık gelen süreç yönergeleri `cpus` ve `memory` olarak adlandırılır, bu nedenle aşağıdaki yapılandırma ima edilir:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Bu değerleri, tüm süreçler için veya belirli adlandırılmış süreçler için, yapılandırma dosyanızda ek süreç yönergeleri kullanarak değiştirebilirsiniz.
Nextflow bunları seçilen yürütücü için uygun talimatlara çevirecektir.

Ancak hangi değerleri kullanacağınızı nasıl bilirsiniz?

### 5.1. Bir kaynak kullanım raporu oluşturmak için iş akışını çalıştırın

Süreçlerinizin ne kadar CPU ve belleğe ihtiyaç duyacağını önceden bilmiyorsanız, biraz kaynak profilleme yapabilirsiniz; yani iş akışını bazı varsayılan tahsislerle çalıştırır, her sürecin ne kadar kullandığını kaydeder ve oradan temel tahsisleri nasıl ayarlayacağınızı tahmin edersiniz.

Uygun bir şekilde, Nextflow bunu yapmak için yerleşik araçlar içerir ve talep üzerine sizin için memnuniyetle bir rapor oluşturacaktır.

Bunu yapmak için komut satırınıza `-with-report <dosyaadı>.html` ekleyin.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Rapor bir html dosyasıdır, indirebilir ve tarayıcınızda açabilirsiniz. Ayrıca soldaki dosya gezgininde sağ tıklayıp eğitim ortamında görüntülemek için `Show preview` seçeneğine tıklayabilirsiniz.

Raporu incelemek ve kaynakları ayarlamak için bazı fırsatları belirleyip belirleyemeyeceğinizi görmek için birkaç dakikanızı ayırın.
Kullanımı tahsis edilenin yüzdesi olarak gösteren sekmelere tıkladığınızdan emin olun.

Mevcut tüm özellikler hakkında belgeler için [Raporlar](https://nextflow.io/docs/latest/reports.html) bölümüne bakın.

### 5.2. Tüm süreçler için kaynak tahsisleri ayarlayın

Profilleme, eğitim iş akışımızdaki süreçlerin çok hafif olduğunu gösteriyor, bu yüzden varsayılan bellek tahsisini süreç başına 1GB'a düşürelim.

`nextflow.config` dosyanıza, pipeline parametreleri bölümünden önce aşağıdakileri ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Süreç ayarları
    */
    process {
        memory = 1.GB
    }

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
    docker.enabled = false
    conda.enabled = true

    /*
    * Pipeline parametreleri
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Bu, tükettiğimiz bilgi işlem miktarını azaltmaya yardımcı olacaktır.

### 5.3. Belirli bir süreç için kaynak tahsisleri ayarlayın

Aynı zamanda, `cowpy` sürecinin diğerlerinden daha fazla kaynak gerektirdiğini varsayacağız, sadece bireysel bir süreç için tahsisleri nasıl ayarlayacağımızı gösterebilmemiz için.

=== "Sonra"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Süreç ayarları
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
    * Süreç ayarları
    */
    process {
        memory = 1.GB
    }
    ```

Bu yapılandırmayla, tüm süreçler 1GB bellek ve tek bir CPU (ima edilen varsayılan) talep edecek, `cowpy` süreci hariç, bu 2GB ve 2 CPU talep edecektir.

!!! tip "İpucu"

    Az CPU'ya sahip bir makineniz varsa ve süreç başına yüksek bir sayı tahsis ederseniz, süreç çağrılarının birbirinin arkasında sıraya girdiğini görebilirsiniz.
    Bunun nedeni Nextflow'un mevcut olandan daha fazla CPU talep etmememizi sağlamasıdır.

### 5.4. Güncellenmiş yapılandırmayla iş akışını çalıştırın

Hadi deneyelim, yapılandırma değişikliklerinden önce ve sonra performansı karşılaştırabilmemiz için profilleme raporu için farklı bir dosya adı sağlayalım.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Bu çok küçük bir iş yükü olduğu için muhtemelen gerçek bir fark fark etmeyeceksiniz, ancak bu gerçek dünya iş akışının performansını ve kaynak gereksinimlerini analiz etmek için kullanacağınız yaklaşımdır.

Süreçlerinizin farklı kaynak gereksinimleri olduğunda çok kullanışlıdır. Tahminlere değil, gerçek verilere dayalı olarak her süreç için ayarladığınız kaynak tahsislerini doğru boyutlandırmanızı sağlar.

!!! tip "İpucu"

    Bu, kaynakların kullanımını optimize etmek için yapabileceklerinizin sadece küçük bir tadımlığıdır.
    Nextflow'un kendisi, kaynak sınırlamaları nedeniyle başarısız olan işleri yeniden denemek için yerleşik gerçekten düzgün [dinamik yeniden deneme mantığına](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) sahiptir.
    Ek olarak, Seqera Platform, kaynak tahsislerinizi otomatik olarak optimize etmek için yapay zeka destekli araçlar da sunar.

### 5.5. Kaynak sınırları ekleyin

Hangi bilgi işlem yürütücüsünü ve bilgi işlem altyapısını kullandığınıza bağlı olarak, tahsis edebileceğiniz (veya tahsis etmeniz gereken) şeyler üzerinde bazı kısıtlamalar olabilir.
Örneğin, kümeniz belirli sınırlar içinde kalmanızı gerektirebilir.

İlgili sınırlamaları ayarlamak için `resourceLimits` yönergesini kullanabilirsiniz. Sözdizimi, bir süreç bloğunda tek başına olduğunda şöyle görünür:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow bu değerleri, belirttiğiniz yürütücüye bağlı olarak uygun talimatlara çevirecektir.

Eğitim ortamında ilgili altyapıya erişimimiz olmadığı için bunu çalıştırmayacağız.
Ancak, bu sınırları aşan kaynak tahsisleriyle iş akışını çalıştırmayı denerseniz, ardından `.command.run` betik dosyasındaki `sbatch` komutuna bakarsanız, yürütücüye gerçekten gönderilen isteklerin `resourceLimits` tarafından belirtilen değerlerle sınırlandırıldığını görürsünüz.

??? info "Kurumsal referans yapılandırmaları"

    nf-core projesi, dünya çapındaki çeşitli kurumlar tarafından paylaşılan, çok çeşitli HPC ve bulut yürütücülerini kapsayan bir [yapılandırma dosyaları koleksiyonu](https://nf-co.re/configs/) derlemiştir.

    Bu paylaşılan yapılandırmalar hem orada çalışan ve dolayısıyla kurumlarının yapılandırmasını kutudan çıkar çıkmaz kullanabilen insanlar hem de kendi altyapıları için bir yapılandırma geliştirmek isteyen insanlar için bir model olarak değerlidir.

### Özet

Kaynak kullanımını değerlendirmek için bir profilleme raporu oluşturmayı ve tüm süreçler ve/veya bireysel süreçler için kaynak tahsislerini nasıl değiştireceğinizi, ayrıca HPC'de çalıştırmak için kaynak sınırlamalarını nasıl ayarlayacağınızı biliyorsunuz.

### Sırada ne var?

Önceden ayarlanmış yapılandırma profillerini nasıl kuracağınızı ve çalışma zamanında bunlar arasında nasıl geçiş yapacağınızı öğrenin.

---

## 6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın

Size, üzerinde çalıştığınız projeye veya kullandığınız bilgi işlem ortamına bağlı olarak pipeline yapılandırmanızı özelleştirebileceğiniz bir dizi yol gösterdik.

Hangi bilgi işlem altyapısını kullandığınıza bağlı olarak alternatif ayarlar arasında geçiş yapmak isteyebilirsiniz. Örneğin, dizüstü bilgisayarınızda yerel olarak küçük ölçekli testler geliştirmek ve çalıştırmak, ardından HPC veya bulutta tam ölçekli iş yüklerini çalıştırmak isteyebilirsiniz.

Nextflow, farklı yapılandırmaları tanımlayan herhangi bir sayıda [profil](https://nextflow.io/docs/latest/config.html#config-profiles) kurmanıza olanak tanır; bunları daha sonra yapılandırma dosyasının kendisini değiştirmek yerine bir komut satırı argümanı kullanarak çalışma zamanında seçebilirsiniz.

### 6.1. Yerel geliştirme ve HPC'de yürütme arasında geçiş yapmak için profiller oluşturun

İki alternatif profil kuralım; biri normal bir bilgisayarda küçük ölçekli yükler çalıştırmak için, burada Docker konteynerlerini kullanacağız, ve biri Slurm zamanlayıcısı olan bir üniversite HPC'sinde çalıştırmak için, burada Conda paketlerini kullanacağız.

#### 6.1.1. Profilleri kurun

`nextflow.config` dosyanıza, pipeline parametreleri bölümünden sonra ancak çıktı ayarlarından önce aşağıdakileri ekleyin:

=== "Sonra"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Pipeline parametreleri
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Profiller
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

    /*
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Önce"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * Pipeline parametreleri
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Çıktı ayarları
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Üniversite HPC'si için kaynak sınırlamalarını da belirttiğimizi görüyorsunuz.

#### 6.1.2. İş akışını bir profille çalıştırın

Nextflow komut satırımızda bir profil belirtmek için `-profile` argümanını kullanırız.

İş akışını `my_laptop` yapılandırmasıyla çalıştırmayı deneyelim.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Gördüğünüz gibi, bu çalışma zamanında yapılandırmalar arasında çok rahat bir şekilde geçiş yapmamızı sağlar.

!!! warning "Uyarı"

    `univ_hpc` profili, bir Slurm zamanlayıcısına erişimimiz olmadığı için eğitim ortamında düzgün çalışmayacaktır.

Gelecekte bunlarla her zaman birlikte ortaya çıkan başka yapılandırma öğeleri bulursak, bunları basitçe ilgili profil(ler)e ekleyebiliriz.
Birlikte gruplandırmak istediğimiz başka yapılandırma öğeleri varsa ek profiller de oluşturabiliriz.

### 6.2. Test parametrelerinin bir profilini oluşturun

Profiller sadece altyapı yapılandırması için değildir.
Bunları iş akışı parametreleri için varsayılan değerler ayarlamak için de kullanabiliriz; bu, başkalarının uygun girdi değerlerini kendileri toplamak zorunda kalmadan iş akışını denemelerini kolaylaştırır.
Bunu bir parametre dosyası kullanmaya alternatif olarak düşünebilirsiniz.

#### 6.2.1. Profili kurun

Bu bağlamda varsayılan değerleri ifade etme sözdizimi şöyle görünür, `test` olarak adlandırdığımız bir profil için:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

İş akışımız için bir test profili eklersek, `profiles` bloğu şöyle olur:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* Profiller
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

Teknik yapılandırma profilleri gibi, istediğiniz herhangi bir ad altında parametreleri belirten birden fazla farklı profil kurabilirsiniz.

#### 6.2.2. İş akışını test profiliyle yerel olarak çalıştırın

Uygun bir şekilde, profiller birbirini dışlamaz, bu nedenle aşağıdaki sözdizimini kullanarak komut satırımızda birden fazla profil belirtebiliriz `-profile <profil1>,<profil2>` (herhangi bir sayıda profil için).

Aynı yapılandırma öğeleri için değerler ayarlayan ve aynı yapılandırma dosyasında açıklanan profilleri birleştirirseniz, Nextflow çatışmayı en son okuduğu değeri kullanarak çözecektir (_yani_ dosyada daha sonra gelen her neyse).
Çakışan ayarlar farklı yapılandırma kaynaklarında ayarlanmışsa, varsayılan [öncelik sırası](https://nextflow.io/docs/latest/config.html) geçerlidir.

Önceki komutumuza test profilini eklemeyi deneyelim:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

Bu, mümkün olduğunda Docker kullanacak ve çıktıları `custom-outdir-config/test` altında üretecek ve bu sefer karakter komedi ikilisi `dragonandcow`.

??? abstract "Dosya içeriği"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

Bu, iş akışı koduyla herhangi bir test veri dosyası dağıttığımız sürece, herkesin komut satırı veya bir parametre dosyası aracılığıyla kendi girdilerini sağlamak zorunda kalmadan iş akışını hızlı bir şekilde deneyebileceği anlamına gelir.

!!! tip "İpucu"

    Harici olarak depolanan daha büyük dosyalar için URL'lere işaret edebiliriz.
    Açık bir bağlantı olduğu sürece Nextflow bunları otomatik olarak indirecektir.

    Daha fazla ayrıntı için, Yan Görev [Dosyalarla Çalışma](../side_quests/working_with_files.md) bölümüne bakın

### 6.3. Çözümlenmiş yapılandırmayı görmek için `nextflow config` kullanın

Yukarıda belirtildiği gibi, bazen aynı parametre birleştirmek istediğiniz profillerde farklı değerlere ayarlanabilir.
Ve daha genel olarak, yapılandırma öğelerinin depolanabileceği çok sayıda yer vardır ve bazen aynı özellikler farklı yerlerde farklı değerlere ayarlanabilir.

Nextflow, herhangi bir çatışmayı çözmek için bir set [öncelik sırası](https://nextflow.io/docs/latest/config.html) uygular, ancak bunu kendiniz belirlemek zor olabilir.
Ve hiçbir şey çakışmasa bile, şeylerin yapılandırılabileceği tüm olası yerlere bakmak sıkıcı olabilir.

Neyse ki, Nextflow, tüm bu süreci sizin için otomatikleştirebilecek `config` adlı kullanışlı bir yardımcı program aracı içerir.

`config` aracı, mevcut çalışma dizininizdeki tüm içerikleri keşfedecek, herhangi bir yapılandırma dosyasını toplayacak ve Nextflow'un iş akışını çalıştırmak için kullanacağı tamamen çözümlenmiş yapılandırmayı üretecektir.
Bu, herhangi bir şey başlatmak zorunda kalmadan hangi ayarların kullanılacağını öğrenmenizi sağlar.

#### 6.3.1. Varsayılan yapılandırmayı çözümleyin

Varsayılan olarak uygulanacak yapılandırmayı çözümlemek için bu komutu çalıştırın.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Bu, komut satırında ekstra bir şey belirtmezseniz aldığınız temel yapılandırmayı gösterir.

#### 6.3.2. Belirli ayarlar etkinleştirilmiş yapılandırmayı çözümleyin

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Bu, birden fazla yapılandırma katmanı içeren karmaşık projeler için özellikle kullanışlı hale gelir.

### Özet

Minimum güçlükle çalışma zamanında önceden ayarlanmış bir yapılandırmayı seçmek için profilleri nasıl kullanacağınızı biliyorsunuz.
Daha genel olarak, iş akışı yürütmelerinizi farklı bilgi işlem platformlarına uyacak şekilde nasıl yapılandıracağınızı ve analizlerinizin tekrarlanabilirliğini nasıl artıracağınızı biliyorsunuz.

### Sırada ne var?

Kutlayın ve kendinize büyük bir alkış verin! İlk Nextflow geliştirici kursunuzu tamamladınız.

Ne öğrendiğinizi gözden geçirmek ve sırada ne olduğunu öğrenmek için son [kurs özetine](./next_steps.md) gidin.

---

## Quiz

<quiz>
Nextflow'un otomatik olarak yüklediği yapılandırma dosyasının adı nedir?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
Aynı parametre hem yapılandırma dosyasında hem de komut satırında ayarlandığında hangisi önceliklidir?
- [ ] Yapılandırma dosyası değeri
- [x] Komut satırı değeri
- [ ] Karşılaşılan ilk değer
- [ ] Hiçbiri; bir hataya neden olur

Daha fazla bilgi: [1.1. Varsayılan değerleri `nextflow.config` dosyasına taşıyın](#11-move-default-values-to-nextflowconfig)
</quiz>

<quiz>
Aynı yapılandırmada hem Docker hem de Conda etkinleştirilebilir mi?
- [x] Evet, Nextflow süreç yönergelerine bağlı olarak her ikisini de kullanabilir
- [ ] Hayır, aynı anda yalnızca biri etkinleştirilebilir
- [ ] Evet, ancak yalnızca profillerde
- [ ] Hayır, birbirini dışlarlar
</quiz>

<quiz>
Hem Docker hem de Conda etkinleştirilmişse ve bir sürecin her iki yönergesi de varsa, hangisine öncelik verilir?
- [x] Docker (konteynerler)
- [ ] Conda
- [ ] Tanımlanan ilk
- [ ] Bir hataya neden olur

Daha fazla bilgi: [3. Bir yazılım paketleme teknolojisi seçin](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Nextflow süreçleri için varsayılan bellek tahsisi nedir?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Sınır yok
</quiz>

<quiz>
Yapılandırma dosyasında belirli bir süreç için kaynak gereksinimlerini nasıl ayarlarsınız?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Daha fazla bilgi: [5.3. Belirli bir süreç için kaynak tahsisleri ayarlayın](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
Bir kaynak kullanım raporu oluşturan komut satırı seçeneği nedir?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Daha fazla bilgi: [5.1. Bir kaynak kullanım raporu oluşturmak için iş akışını çalıştırın](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
`resourceLimits` yönergesi ne yapar?
- [ ] Minimum kaynak gereksinimlerini ayarlar
- [ ] Süreçlere kaynak tahsis eder
- [x] Talep edilebilecek maksimum kaynakları sınırlar
- [ ] Kaynak kullanımını izler

Daha fazla bilgi: [5.5. Kaynak sınırları ekleyin](#55-add-resource-limits)
</quiz>

<quiz>
Nextflow'da varsayılan yürütücü nedir?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Daha fazla bilgi: [4. Bir yürütme platformu seçin](#4-select-an-execution-platform)
</quiz>

<quiz>
Nextflow çalıştırırken bir parametre dosyasını nasıl belirtirsiniz?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Daha fazla bilgi: [1.3. Bir parametre dosyası kullanın](#13-use-a-parameter-file)
</quiz>

<quiz>
Profiller ne için kullanılabilir? (Tümünü seçin)
- [x] Altyapıya özgü ayarları tanımlama
- [x] Farklı ortamlar için kaynak sınırları ayarlama
- [x] Test parametreleri sağlama
- [ ] Yeni süreçler tanımlama

Daha fazla bilgi: [6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
Tek bir komutta birden fazla profili nasıl belirtirsiniz?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Daha fazla bilgi: [6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
