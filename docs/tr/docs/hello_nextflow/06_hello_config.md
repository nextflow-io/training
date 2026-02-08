# Bölüm 6: Hello Config

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=tr" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube kanalında [tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) izleyin.

:green_book: Video transkripti [burada](./transcripts/06_hello_config.md) mevcuttur.
///

Bu bölüm, _iş akışı kodunun tek bir satırını değiştirmeden_ davranışını özelleştirebilmeniz, farklı ortamlara uyarlayabilmeniz ve kaynak kullanımını optimize edebilmeniz için Nextflow pipeline'ınızın yapılandırmasını nasıl kuracağınızı ve yöneteceğinizi keşfedecektir.

Bunu yapmanın birden fazla yolu vardır; bunlar birlikte kullanılabilir ve [burada](https://nextflow.io/docs/latest/config.html) açıklanan öncelik sırasına göre yorumlanır.

Bu kursun bu bölümünde, Bölüm 5: Hello Containers'da zaten karşılaştığınız en basit ve en yaygın yapılandırma dosyası mekanizması olan [`nextflow.config`](https://nextflow.io/docs/latest/config.html) dosyasını göstereceğiz.

Süreç yönergeleri, yürütücüler, profiller ve parametre dosyaları gibi Nextflow yapılandırmasının temel bileşenlerini gözden geçireceğiz.
Bu yapılandırma seçeneklerini etkin bir şekilde kullanmayı öğrenerek, pipeline'larınızın esnekliğini, ölçeklenebilirliğini ve performansını artırabilirsiniz.

??? info "Bu bölümden nasıl başlanır"

    Bu kursun bu bölümü, [Hello Nextflow](./index.md) kursunun 1-5. Bölümlerini tamamladığınızı ve eksiksiz çalışan bir pipeline'ınız olduğunu varsayar.

    Kursa bu noktadan başlıyorsanız, `modules` dizinini ve `nextflow.config` dosyasını çözümlerden kopyalamanız gerekecek:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    `nextflow.config` dosyası, Docker konteynerlerinin kullanımını etkinleştiren `docker.enabled = true` satırını içerir.

    Hello pipeline'ına aşina değilseniz veya bir hatırlatmaya ihtiyacınız varsa, [bu bilgi sayfasına](../info/hello_pipeline.md) bakın.

---

## 0. Isınma: `hello-config.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-config.nf` iş akışı betiğini kullanacağız.
Bu betik, bu eğitim kursunun 5. Bölümünde üretilen betiğe eşdeğerdir; ancak çıktı hedeflerini değiştirdik:

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

Daha önce olduğu gibi, çıktı dosyalarını `output` bloğunda belirtilen dizinde bulacaksınız (`results/hello_config/`).

??? abstract "Dizin içerikleri"

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

Son ASCII sanat çıktısı `results/hello_config/` dizininde, `cowpy-COLLECTED-batch-output.txt` adı altında.

??? abstract "Dosya içerikleri"

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

Şimdiye kadar üzerinde çalıştığımız şeyin basit bir uzantısı olan bir yapılandırma yönüyle başlayacağız: girdi parametrelerinin yönetimi.

Şu anda, iş akışımız komut satırı aracılığıyla çeşitli parametre değerlerini kabul edecek şekilde ayarlanmış, varsayılan değerler iş akışı betiğinin kendisindeki bir `params` bloğunda ayarlanmış.
Ancak, bu varsayılanları, komut satırında parametreler belirtmek veya orijinal betik dosyasını değiştirmek zorunda kalmadan geçersiz kılmak isteyebilirsiniz.

Bunu yapmanın birden fazla yolu vardır; size çok yaygın olarak kullanılan üç temel yolu göstereceğiz.

### 1.1. Varsayılan değerleri `nextflow.config`'e taşıyın

Bu en basit yaklaşımdır, ancak ana `nextflow.config` dosyası her çalıştırma için düzenlemek isteyeceğiniz bir şey olmadığından muhtemelen en az esnek olanıdır.
Ancak, parametreleri iş akışında _tanımlamak_ (kesinlikle oraya ait) ile _varsayılan değerler_ sağlamak (bir yapılandırma dosyasında daha uygun) arasındaki endişeleri ayırma avantajına sahiptir.

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
Sözdizimi biraz farklı.
İş akışı dosyasında, bunlar tip tanımlı bildirimlerdir.
Yapılandırmada, bunlar değer atamalarıdır.

Teknik olarak, bu hâlâ iş akışı dosyasında belirtilen varsayılan değerleri geçersiz kılmak için yeterlidir.
Karakteri değiştirebilir, örneğin, ve yapılandırma dosyasında ayarlanan değerin iş akışı dosyasında ayarlanan değeri geçersiz kıldığından emin olmak için iş akışını çalıştırabilirsiniz.

Ancak, yapılandırmayı tamamen yapılandırma dosyasına taşıma ruhuna uygun olarak, bu değerleri iş akışı dosyasından tamamen kaldıralım.

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

Doğru çalışıp çalışmadığını test edelim.

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

Bu hâlâ daha önce olduğu gibi aynı çıktıyı üretiyor.

Son ASCII sanat çıktısı `results/hello_config/` dizininde, daha önce olduğu gibi `cowpy-COLLECTED-batch-output.txt` adı altında.

??? abstract "Dosya içerikleri"

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

Harika, ama bazen ana yapılandırma dosyasıyla uğraşmadan farklı varsayılan değerlerle bazı geçici deneyler yapmak isteyebilirsiniz.
Bunu, deneyeleriniz için çalışma dizini olarak kullanacağınız bir alt dizinde yeni bir `nextflow.config` dosyası oluşturarak yapabilirsiniz.

#### 1.2.1. Boş bir yapılandırma ile çalışma dizini oluşturun

Yeni bir dizin oluşturarak ve içine girerek başlayalım:

```bash
mkdir -p tux-run
cd tux-run
```

Ardından, o dizinde boş bir yapılandırma dosyası oluşturun:

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

Bu, `tux-run/work/` ve `tux-run/results/` dahil olmak üzere `tux-run/` altında yeni bir dizi dizin oluşturacaktır.

Bu çalıştırmada, Nextflow mevcut dizinimizdeki `nextflow.config`'i pipeline'ın kök dizinindeki `nextflow.config` ile birleştirir ve böylece varsayılan karakteri (turkey) tux karakteriyle geçersiz kılar.

Son çıktı dosyası, tux karakterinin selamlamaları söylediğini içermelidir.

??? abstract "Dosya içerikleri"

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

Şimdi parametre değerlerini ayarlamanın başka kullanışlı bir yoluna bakalım.

### 1.3. Bir parametre dosyası kullanın

Alt dizin yaklaşımı deney yapmak için harika çalışır, ancak biraz kurulum gerektirir ve yolları buna göre uyarlamanızı gerektirir.
Pipeline'ınızı belirli bir değer kümesiyle çalıştırmak veya başkasının minimum çabayla yapmasını sağlamak istediğinizde daha basit bir yaklaşım vardır.

Nextflow, parametreleri YAML veya JSON formatında bir [parametre dosyası](https://nextflow.io/docs/latest/config.html#params-file) aracılığıyla belirtmemize olanak tanır; bu, örneğin alternatif varsayılan değer kümelerini ve çalıştırmaya özgü parametre değerlerini yönetmeyi ve dağıtmayı çok kullanışlı hale getirir.

#### 1.3.1. Örnek parametre dosyasını inceleyin

Bunu göstermek için, mevcut dizinde `test-params.yaml` adında bir örnek parametre dosyası sağlıyoruz:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Bu parametre dosyası, belirtmek istediğimiz girdilerin her biri için bir anahtar-değer çifti içerir.
Sözdizimini yapılandırma dosyasıyla karşılaştırırsanız, eşit işaretleri (`=`) yerine iki nokta üst üste (`:`) kullanıldığına dikkat edin.
Yapılandırma dosyası Groovy'de yazılmışken, parametre dosyası YAML'de yazılmıştır.

!!! info "Bilgi"

    Ayrıca parametre dosyasının bir JSON versiyonunu örnek olarak sağlıyoruz ancak burada onunla çalıştırmayacağız.
    Onu kendi başınıza denemekten çekinmeyin.

#### 1.3.2. Pipeline'ı çalıştırın

İş akışını bu parametre dosyasıyla çalıştırmak için, temel komuta `-params-file <dosyaadı>` eklemeniz yeterli.

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

Son çıktı dosyası, stegosaurus karakterinin selamlamaları söylediğini içermelidir.

??? abstract "Dosya içerikleri"

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

Yalnızca birkaç parametreniz olduğunda bir parametre dosyası kullanmak aşırı gibi görünebilir, ancak bazı pipeline'lar düzinelerce parametre bekler.
Bu durumlarda, bir parametre dosyası kullanmak, büyük komut satırları yazmak zorunda kalmadan ve iş akışı betiğini değiştirmeden çalışma zamanında parametre değerleri sağlamamıza olanak tanıyacaktır.

Ayrıca parametre kümelerini işbirlikçilere veya örneğin bir yayın için destekleyici bilgi olarak dağıtmayı kolaylaştırır.
Bu, çalışmanızı başkaları tarafından daha tekrarlanabilir hale getirir.

### Özet

İş akışı girdilerini yönetmek için temel yapılandırma seçeneklerinden nasıl yararlanacağınızı biliyorsunuz.

### Sırada ne var?

İş akışı çıktılarınızın nerede ve nasıl yayınlanacağını nasıl yöneteceğinizi öğrenin.

---

## 2. İş akışı çıktılarını yönetin

Şimdiye kadar iş akışı düzeyindeki çıktı tanımlamaları için tüm yolları sabit kodluyorduk ve birden fazla çıktı eklemeye başladığımızda belirttiğimiz gibi, biraz tekrar olabilir.

Bunu daha esnek yapılandırmak için birkaç yaygın yola bakalım.

### 2.1. `-output-dir` ile çıktı dizinini özelleştirin

'Yayınlanmış' çıktılarımızın nasıl organize edileceğini kontrol ederken iki farklı önceliğimiz var:

- Üst düzey çıktı dizini
- Dosyaların bu dizin içinde nasıl organize edileceği

Şimdiye kadar varsayılan üst düzey dizini kullanıyorduk: `results`.
Bunu `-output-dir` CLI seçeneğini kullanarak özelleştirmeyle başlayalım.

#### 2.1.1. Pipeline'ı `-output-dir` ile çalıştırın

`-output-dir` seçeneği (kısa hali: `-o`) tüm iş akışı çıktıları için varsayılan çıktı dizinini (`results/`) geçersiz kılar.
Bu, çıktıların yayınlandığı kök yolu kontrol etmek için önerilen yoldur.

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

Bu, çıktıları `results/` yerine `custom-outdir-cli/` altında yayınlar:

??? abstract "Dizin içerikleri"

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

Çıktı bloğundaki `path` bildirimlerinden gelen `hello_config` alt dizinine hâlâ sahip olduğumuza dikkat edin.
Bunu temizleyelim.

#### 2.1.2. Çıktı bloğundan sabit kodlanmış yolları kaldırın

`hello_config/` öneki önceki bölümlerde sabit kodlanmıştı, ancak artık çıktı yollarını esnek bir şekilde yapılandırmayı öğrendiğimize göre, bu sabit kodlamayı kaldırabiliriz.
Alt dizine ihtiyaç duymayan çıktılar için `path` yönergesini boş bir dizeye ayarlayabiliriz veya tamamen kaldırabiliriz.

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

Şimdi çıktılar `hello_config` alt dizini olmadan doğrudan `custom-outdir-cli-2/` altında yayınlanıyor:

??? abstract "Dizin içerikleri"

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

Bu, çıktı dizinini `custom-outdir-config/` artı alt dizin olarak `batch` parametresinin değerine ayarlar.
Şimdi `--batch` parametresini ayarlayarak çıktı konumunu değiştirebilirsiniz:

```bash
nextflow run hello-config.nf --batch my_run
```

Bu, çıktıları `custom-outdir-config/my_run/` altında yayınlar.

!!! note "Not"

    `-output-dir` CLI seçeneği, `outputDir` yapılandırma ayarından önceliklidir.
    Eğer ayarlanmışsa, yapılandırma seçeneği tamamen yok sayılır.

#### 2.2.2. Grup ve süreç adlarıyla alt dizinler

Ayrıca, çıktı başına alt dizin çıktı `path` bildirimlerini dinamik olarak ayarlayabiliriz.

Örneğin, çıktı yolu bildiriminde `<süreç>.name`'e referans vererek çıktılarımızı sürece göre organize edebiliriz:

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

Yukarıdaki düzenlemede `intermediates` ile son çıktıların üst düzeyde olması arasındaki ayrımı sildik.
Bunu geri koyalım ve ayrıca dosyaları bir `params.batch` alt dizinine koyalım.

!!! tip "İpucu"

    `params.batch`'i çıktı bloğu `path`'ine dahil etmek, `outputDir` yapılandırması yerine, CLI'da `-output-dir` ile üzerine yazılmayacağı anlamına gelir.

Önce, yapılandırma dosyasını `outputDir`'den `${params.batch}`'i kaldıracak şekilde güncelleyin (çünkü onu yol bildirimlerine taşıyoruz):

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

Bunun pratikte nasıl çalıştığını görelim, hem `-output-dir` (veya kısaca `-o`) değerini `custom-outdir-config-2` hem de grup adını komut satırından `rep2` olarak ayarlayarak:

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

Bu, çıktıları `custom-outdir-config-2/rep2/` altında yayınlar, belirtilen temel yol _ve_ grup adı alt dizini _ve_ sürece göre gruplanmış sonuçlarla:

??? abstract "Dizin içerikleri"

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

Son olarak, tekrarlayan kod miktarını azaltma ruhuna uygun olarak, çıktı başına `mode` tanımlarını yapılandırmada tek bir satırla değiştirebiliriz.

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

Tıpkı `outputDir` seçeneği gibi, yapılandırma dosyasında `workflow.output.mode`'a bir değer vermek, iş akışı dosyasında ayarlananı geçersiz kılmak için yeterli olurdu, ancak yine de gereksiz kodu kaldıralım.

#### 2.3.2. Çıktı modunu iş akışı dosyasından kaldırın

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

Bu daha kısa, değil mi?

#### 2.3.3. Pipeline'ı çalıştırın

Doğru çalışıp çalışmadığını test edelim:

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

Bu, çıktıları `config-output-mode/` altında yayınlar ve hepsi hâlâ düzgün kopyalar, sembolik bağlantılar değil.

??? abstract "Dizin içerikleri"

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

Çıktı başına mod ayarlama yolunu hâlâ kullanmak istemenizin ana nedeni, aynı iş akışı içinde karıştırıp eşleştirmek istemeniz olabilir, _yani_ bazı çıktıların kopyalanması ve bazılarının sembolik bağlantı olması.

Bu şekilde özelleştirebileceğiniz birçok başka seçenek var, ancak umarız bu size seçenek yelpazesi ve tercihlerinize uygun olarak bunları etkili bir şekilde nasıl kullanacağınız hakkında bir fikir verir.

### Özet

Çıktılarınızın yayınlandığı dizinlerin adlandırılmasını ve yapısını ve iş akışı çıktısı yayınlama modunu nasıl kontrol edeceğinizi biliyorsunuz.

### Sırada ne var?

İş akışı yapılandırmanızı yazılım paketleme teknolojisinden başlayarak hesaplama ortamınıza nasıl uyarlayacağınızı öğrenin.

---

## 3. Bir yazılım paketleme teknolojisi seçin

Şimdiye kadar girdilerin nasıl girdiğini ve çıktıların nereden çıktığını kontrol eden yapılandırma öğelerine bakıyorduk. Şimdi iş akışı yapılandırmanızı hesaplama ortamınıza uyarlamaya daha özel olarak odaklanma zamanı.

Bu yoldaki ilk adım, her adımda çalıştırılacak yazılım paketlerinin nereden geleceğini belirtmektir.
Yerel hesaplama ortamında zaten yüklü mü?
İmajları alıp bir konteyner sistemi aracılığıyla çalıştırmamız mı gerekiyor?
Yoksa Conda paketlerini alıp yerel bir Conda ortamı mı oluşturmamız gerekiyor?

Bu eğitim kursunun en başında (Bölüm 1-4) iş akışımızda sadece yerel olarak yüklenmiş yazılımı kullandık.
Ardından Bölüm 5'te, Docker konteynerlerini ve Docker konteynerlerinin kullanımını etkinleştirmek için kullandığımız `nextflow.config` dosyasını tanıttık.

Şimdi `nextflow.config` dosyası aracılığıyla alternatif bir yazılım paketleme seçeneğini nasıl yapılandırabileceğimizi görelim.

### 3.1. Yapılandırma dosyasında Docker'ı devre dışı bırakın ve Conda'yı etkinleştirin

Bir HPC kümesinde çalıştığımızı ve yöneticinin güvenlik nedeniyle Docker kullanımına izin vermediğini varsayalım.
Neyse ki bizim için, Nextflow, Singularity (HPC'de daha yaygın olarak kullanılır) gibi birden fazla başka konteyner teknolojisini ve Conda gibi yazılım paket yöneticilerini destekler.

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

Bu, Nextflow'un Conda paketleri belirtilmiş süreçler için Conda ortamları oluşturmasına ve kullanmasına olanak tanıyacaktır.
Bu, şimdi `cowpy` sürecimize bunlardan birini eklememiz gerektiği anlamına gelir!

### 3.2. Süreç tanımında bir Conda paketi belirtin

`cowpy` aracını içeren bir Conda paketi için URI'yi zaten aldık: `conda-forge::cowpy==1.1.5`

Şimdi `conda` yönergesini kullanarak URI'yi `cowpy` süreç tanımına ekliyoruz:

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

Açık olmak gerekirse, `container` yönergesini _değiştirmiyoruz_, alternatif bir seçenek _ekliyoruz_.

!!! tip "İpucu"

    Belirli bir conda paketi için URI almanın birkaç farklı yolu var.
    [Seqera Containers](https://seqera.io/containers/) arama sorgusunu kullanmanızı öneriyoruz; bu, ondan bir konteyner oluşturmayı planlamasanız bile kopyalayıp yapıştırabileceğiniz bir URI verecektir.

### 3.3. Conda kullanabildiğini doğrulamak için iş akışını çalıştırın

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

Perde arkasında, Nextflow Conda paketlerini aldı ve ortamı oluşturdu, bu normalde biraz iş gerektirir; bu yüzden bunların hiçbirini kendimiz yapmak zorunda kalmamamız güzel!

!!! note "Not"

    Bu hızlı çalışır çünkü `cowpy` paketi oldukça küçüktür, ancak büyük paketlerle çalışıyorsanız, ilk seferde normalden biraz daha uzun sürebilir ve konsol çıktısının tamamlanmadan önce bir dakika kadar 'takılı' kaldığını görebilirsiniz.
    Bu normaldir ve Nextflow'un yeni bir paketi ilk kullandığınızda yaptığı ekstra işten kaynaklanır.

Bizim açımızdan, arka uçta mekanikler biraz farklı olsa da Docker ile çalıştırmakla tamamen aynı görünüyor.

Bu, gerektiğinde Conda ortamlarıyla çalıştırmaya hazır olduğumuz anlamına gelir.

??? info "Docker ve Conda'yı karıştırıp eşleştirme"

    Bu yönergeler süreç başına atandığından, 'karıştırıp eşleştirmek' mümkündür, _yani_ kullandığınız hesaplama altyapısı her ikisini de destekliyorsa, iş akışınızdaki bazı süreçleri Docker ile ve diğerlerini Conda ile çalıştıracak şekilde yapılandırabilirsiniz.
    Bu durumda, yapılandırma dosyanızda hem Docker hem de Conda'yı etkinleştirirsiniz.
    Her ikisi de belirli bir süreç için mevcutsa, Nextflow konteynerlere öncelik verecektir.

    Ve daha önce belirtildiği gibi, Nextflow birden fazla başka yazılım paketleme ve konteyner teknolojisini destekler, bu nedenle yalnızca bu ikisiyle sınırlı değilsiniz.

### Özet

Her sürecin hangi yazılım paketini kullanması gerektiğini nasıl yapılandıracağınızı ve teknolojiler arasında nasıl geçiş yapacağınızı biliyorsunuz.

### Sırada ne var?

Nextflow'un işi gerçekten yapmak için kullandığı yürütme platformunu nasıl değiştireceğinizi öğrenin.

---

## 4. Bir yürütme platformu seçin

Şimdiye kadar, pipeline'ımızı local yürütücüyle çalıştırıyorduk.
Bu, her görevi Nextflow'un çalıştığı makinede yürütür.
Nextflow başladığında, mevcut CPU'lara ve belleğe bakar.
Çalıştırılmaya hazır görevlerin kaynakları mevcut kaynakları aşarsa, Nextflow son görevleri, önceki görevlerden biri veya daha fazlası tamamlanıp gerekli kaynakları serbest bırakana kadar yürütmeden alıkoyacaktır.

Local yürütücü kullanışlı ve verimlidir, ancak o tek makineyle sınırlıdır. Çok büyük iş yükleri için, yerel makinenizin bir darboğaz olduğunu keşfedebilirsiniz; ya mevcut olandan daha fazla kaynak gerektiren tek bir göreviniz olduğu için ya da tek bir makinenin bunları çalıştırmasını beklemenin çok uzun süreceği kadar çok göreviniz olduğu için.

Nextflow, HPC zamanlayıcıları (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor ve diğerleri) dahil olmak üzere [birçok farklı yürütme arka ucunu](https://nextflow.io/docs/latest/executor.html) ve ayrıca bulut yürütme arka uçlarını (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes ve daha fazlası) destekler.

### 4.1. Farklı bir arka ucu hedefleme

Yürütücü seçimi, `executor` adlı bir süreç yönergesiyle belirlenir.
Varsayılan olarak `local` olarak ayarlanmıştır, bu nedenle aşağıdaki yapılandırma ima edilir:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Yürütücüyü farklı bir arka ucu hedefleyecek şekilde ayarlamak için, kaynak tahsisleri için yukarıda açıklandığı gibi benzer sözdizimi kullanarak istediğiniz yürütücüyü belirtmeniz yeterlidir (tüm seçenekler için [yürütücü dokümantasyonuna](https://nextflow.io/docs/latest/executor.html) bakın).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Uyarı"

    Bunu eğitim ortamında gerçekten test edemeyiz çünkü bir HPC'ye bağlanacak şekilde ayarlanmamış.

### 4.2. Yürütme parametreleri için arka uca özgü sözdizimi ile başa çıkma

Çoğu yüksek performanslı hesaplama platformu, kaynak tahsisi istekleri ve sınırlamaları (örn. CPU sayısı ve bellek) ve kullanılacak iş kuyruğunun adı gibi belirli parametreleri belirtmenize izin verir (ve bazen gerektirir).

Ne yazık ki, bu sistemlerin her biri, bir işin nasıl tanımlanması ve ilgili zamanlayıcıya nasıl gönderilmesi gerektiğini belirlemek için farklı teknolojiler, sözdizimler ve yapılandırmalar kullanır.

??? abstract "Örnekler"

    Örneğin, 8 CPU ve 4GB RAM gerektiren ve "my-science-work" kuyruğunda yürütülecek aynı iş, arka uca bağlı olarak aşağıdaki farklı şekillerde ifade edilmelidir.

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
[`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) ve [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) gibi ilgili özellikleri (diğer özellikler için [süreç yönergelerine](https://nextflow.io/docs/latest/reference/process.html#process-directives) bakın) yalnızca bir kez belirtebilmeniz için standartlaştırılmış bir sözdizimi sağlar.
Ardından, çalışma zamanında, Nextflow bu ayarları kullanarak yürütücü ayarına dayalı olarak uygun arka uca özgü betikleri oluşturacaktır.

Bu standartlaştırılmış sözdizimini bir sonraki bölümde ele alacağız.

### Özet

Artık farklı türde hesaplama altyapısı kullanmak için yürütücüyü nasıl değiştireceğinizi biliyorsunuz.

### Sırada ne var?

Nextflow'da kaynak tahsislerini ve sınırlamalarını nasıl değerlendirip ifade edeceğinizi öğrenin.

---

## 5. Hesaplama kaynak tahsislerini kontrol edin

Çoğu yüksek performanslı hesaplama platformu, CPU sayısı ve bellek gibi belirli kaynak tahsis parametrelerini belirtmenize izin verir (ve bazen gerektirir).

Varsayılan olarak, Nextflow her süreç için tek bir CPU ve 2GB bellek kullanacaktır.
İlgili süreç yönergeleri `cpus` ve `memory` olarak adlandırılır, bu nedenle aşağıdaki yapılandırma ima edilir:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Bu değerleri, yapılandırma dosyanızdaki ek süreç yönergelerini kullanarak tüm süreçler veya belirli adlandırılmış süreçler için değiştirebilirsiniz.
Nextflow bunları seçilen yürütücü için uygun talimatlara çevirecektir.

Ancak hangi değerleri kullanacağınızı nasıl bilirsiniz?

### 5.1. Kaynak kullanım raporu oluşturmak için iş akışını çalıştırın

Süreçlerinizin ne kadar CPU ve belleğe ihtiyaç duyacağını önceden bilmiyorsanız, bazı kaynak profilleme yapabilirsiniz; yani iş akışını bazı varsayılan tahsislerle çalıştırırsınız, her sürecin ne kadar kullandığını kaydedersiniz ve oradan temel tahsisleri nasıl ayarlayacağınızı tahmin edersiniz.

Kullanışlı bir şekilde, Nextflow bunu yapmak için yerleşik araçlar içerir ve istek üzerine sizin için bir rapor oluşturmaktan mutluluk duyar.

Bunu yapmak için, komut satırınıza `-with-report <dosyaadı>.html` ekleyin.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Rapor, tarayıcınızda indirip açabileceğiniz bir html dosyasıdır. Eğitim ortamında görüntülemek için soldaki dosya gezgininde sağ tıklayıp `Show preview`'a da tıklayabilirsiniz.

Raporu incelemek ve kaynakları ayarlama fırsatlarını belirleyip belirleyemeyeceğinizi görmek için birkaç dakika ayırın.
Kullanım sonuçlarını tahsis edilenin yüzdesi olarak gösteren sekmelere tıkladığınızdan emin olun.

Mevcut tüm özellikleri açıklayan [Raporlar](https://nextflow.io/docs/latest/reports.html) dokümantasyonuna bakın.

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

Bu, tükettiğimiz hesaplama miktarını azaltmaya yardımcı olacaktır.

### 5.3. Belirli bir süreç için kaynak tahsisleri ayarlayın

Aynı zamanda, `cowpy` sürecinin diğerlerinden daha fazla kaynak gerektirdiğini varsayacağız, böylece bireysel bir süreç için tahsisleri nasıl ayarlayacağımızı gösterebiliriz.

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

Bu yapılandırmayla, `cowpy` süreci dışında tüm süreçler 1GB bellek ve tek bir CPU (ima edilen varsayılan) isteyecektir; `cowpy` süreci 2GB ve 2 CPU isteyecektir.

!!! tip "İpucu"

    Az sayıda CPU'ya sahip bir makineniz varsa ve süreç başına yüksek sayıda tahsis ederseniz, süreç çağrılarının birbiri ardına sıraya girdiğini görebilirsiniz.
    Bunun nedeni, Nextflow'un mevcut olandan daha fazla CPU talep etmememizi sağlamasıdır.

### 5.4. Güncellenmiş yapılandırma ile iş akışını çalıştırın

Bunu deneyelim, yapılandırma değişikliklerinden önce ve sonra performansı karşılaştırabilmemiz için profilleme raporu için farklı bir dosya adı sağlayarak.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Bu kadar küçük bir iş yükü olduğundan muhtemelen gerçek bir fark fark etmeyeceksiniz, ancak bu, gerçek dünya bir iş akışının performansını ve kaynak gereksinimlerini analiz etmek için kullanacağınız yaklaşımdır.

Süreçlerinizin farklı kaynak gereksinimleri olduğunda çok faydalıdır. Tahminde bulunmak yerine gerçek verilere dayalı olarak her süreç için ayarladığınız kaynak tahsislerini doğru boyutlandırmanızı sağlar.

!!! tip "İpucu"

    Bu, kaynak kullanımınızı optimize etmek için yapabileceklerinizin sadece küçük bir tadımıdır.
    Nextflow'un kendisi, kaynak sınırlamaları nedeniyle başarısız olan işleri yeniden denemek için gerçekten zarif [dinamik yeniden deneme mantığı](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) yerleşik olarak içerir.
    Ek olarak, Seqera Platform kaynak tahsislerinizi otomatik olarak optimize etmek için yapay zeka destekli araçlar da sunmaktadır.

### 5.5. Kaynak sınırları ekleyin

Hangi hesaplama yürütücüsü ve hesaplama altyapısı kullandığınıza bağlı olarak, tahsis edebileceğiniz (veya etmeniz gereken) konusunda bazı kısıtlamalar olabilir.
Örneğin, kümeniz belirli sınırlar içinde kalmanızı gerektirebilir.

İlgili sınırlamaları ayarlamak için `resourceLimits` yönergesini kullanabilirsiniz. Sözdizimi, bir process bloğunda tek başına olduğunda şöyle görünür:

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

Bunu çalıştırmayacağız, çünkü eğitim ortamında ilgili altyapıya erişimimiz yok.
Ancak, bu sınırları aşan kaynak tahsisleriyle iş akışını çalıştırmayı deneyecek olsaydınız, ardından `.command.run` betik dosyasındaki `sbatch` komutunu ararsanız, yürütücüye gönderilen isteklerin `resourceLimits` tarafından belirtilen değerlerde sınırlandığını görürdünüz.

??? info "Kurumsal referans yapılandırmaları"

    nf-core projesi, çok çeşitli HPC ve bulut yürütücülerini kapsayan, dünya çapındaki çeşitli kurumlar tarafından paylaşılan bir [yapılandırma dosyaları koleksiyonu](https://nf-co.re/configs/) derlemiştir.

    Bu paylaşılan yapılandırmalar, hem orada çalışan ve dolayısıyla kurumlarının yapılandırmasını kutudan çıktığı gibi kullanabilen insanlar için hem de kendi altyapıları için bir yapılandırma geliştirmek isteyen insanlar için model olarak değerlidir.

### Özet

Kaynak kullanımını değerlendirmek için profilleme raporu oluşturmayı ve tüm süreçler ve/veya bireysel süreçler için kaynak tahsislerini nasıl değiştireceğinizi ve HPC'de çalıştırma için kaynak sınırlamalarını nasıl ayarlayacağınızı biliyorsunuz.

### Sırada ne var?

Önceden ayarlanmış yapılandırma profillerini nasıl kuracağınızı ve çalışma zamanında aralarında nasıl geçiş yapacağınızı öğrenin.

---

## 6. Önceden ayarlanmış yapılandırmalar arasında geçiş yapmak için profilleri kullanın

Size, üzerinde çalıştığınız projeye veya kullandığınız hesaplama ortamına bağlı olarak pipeline yapılandırmanızı özelleştirebileceğiniz birkaç yol gösterdik.

Hangi hesaplama altyapısını kullandığınıza bağlı olarak alternatif ayarlar arasında geçiş yapmak isteyebilirsiniz. Örneğin, dizüstü bilgisayarınızda yerel olarak geliştirip küçük ölçekli testler yapmak, ardından tam ölçekli iş yüklerini HPC veya bulutta çalıştırmak isteyebilirsiniz.

Nextflow, farklı yapılandırmaları tanımlayan herhangi bir sayıda [profil](https://nextflow.io/docs/latest/config.html#config-profiles) ayarlamanıza olanak tanır; bunları yapılandırma dosyasını değiştirmek yerine bir komut satırı argümanı kullanarak çalışma zamanında seçebilirsiniz.

### 6.1. Yerel geliştirme ve HPC'de yürütme arasında geçiş yapmak için profiller oluşturun

İki alternatif profil kuralım; biri normal bir bilgisayarda küçük ölçekli yükler çalıştırmak için, burada Docker konteynerlerini kullanacağız, ve biri Slurm zamanlayıcısına sahip bir üniversite HPC'sinde çalıştırmak için, burada Conda paketlerini kullanacağız.

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

Nextflow komut satırımızda bir profil belirtmek için `-profile` argümanını kullanıyoruz.

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

Gördüğünüz gibi, bu yapılandırmalar arasında çalışma zamanında çok kullanışlı bir şekilde geçiş yapmamıza olanak tanır.

!!! warning "Uyarı"

    `univ_hpc` profili eğitim ortamında düzgün çalışmayacaktır çünkü bir Slurm zamanlayıcısına erişimimiz yok.

Gelecekte bunlarla her zaman birlikte olan başka yapılandırma öğeleri bulursak, bunları ilgili profile/profillere ekleyebiliriz.
Birlikte gruplandırmak istediğimiz başka yapılandırma öğeleri varsa ek profiller de oluşturabiliriz.

### 6.2. Test parametreleri profili oluşturun

Profiller yalnızca altyapı yapılandırması için değildir.
Başkalarının uygun girdi değerlerini kendileri toplamak zorunda kalmadan iş akışını denemesini kolaylaştırmak için iş akışı parametreleri için varsayılan değerler ayarlamak amacıyla da kullanabiliriz.
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

Teknik yapılandırma profilleri gibi, istediğiniz herhangi bir rastgele ad altında parametreleri belirten birden fazla farklı profil ayarlayabilirsiniz.

#### 6.2.2. İş akışını test profiliyle yerel olarak çalıştırın

Kullanışlı bir şekilde, profiller birbirini dışlamaz, bu nedenle `-profile <profil1>,<profil2>` sözdizimini (herhangi bir sayıda profil için) kullanarak komut satırımızda birden fazla profil belirtebiliriz.

Aynı yapılandırma öğeleri için değerler ayarlayan ve aynı yapılandırma dosyasında tanımlanan profilleri birleştirirseniz, Nextflow çatışmayı en son okuduğu değeri kullanarak çözecektir (_yani_ dosyada daha sonra gelen).
Çakışan ayarlar farklı yapılandırma kaynaklarında ayarlanmışsa, varsayılan [öncelik sırası](https://nextflow.io/docs/latest/config.html) geçerlidir.

Önceki komutumuzda test profilini eklemeyi deneyelim:

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

Bu, mümkün olduğunda Docker kullanacak ve çıktıları `custom-outdir-config/test` altında üretecektir ve bu sefer karakter komik ikili `dragonandcow`'dur.

??? abstract "Dosya içerikleri"

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

Bu, test veri dosyalarını iş akışı koduyla birlikte dağıttığımız sürece, herkesin komut satırı veya parametre dosyası aracılığıyla kendi girdilerini sağlamak zorunda kalmadan iş akışını hızlıca deneyebileceği anlamına gelir.

!!! tip "İpucu"

    Harici olarak depolanan daha büyük dosyalar için URL'lere işaret edebiliriz.
    Açık bir bağlantı olduğu sürece Nextflow bunları otomatik olarak indirecektir.

    Daha fazla ayrıntı için, [Dosyalarla Çalışma](../side_quests/working_with_files.md) Yan Görevi'ne bakın

### 6.3. Çözümlenmiş yapılandırmayı görmek için `nextflow config` kullanın

Yukarıda belirtildiği gibi, bazen aynı parametre birleştirmek istediğiniz profillerde farklı değerlere ayarlanabilir.
Ve daha genel olarak, yapılandırma öğelerinin depolanabileceği çok sayıda yer vardır ve bazen aynı özellikler farklı yerlerde farklı değerlere ayarlanabilir.

Nextflow, herhangi bir çatışmayı çözmek için belirlenmiş bir [öncelik sırası](https://nextflow.io/docs/latest/config.html) uygular, ancak bunu kendiniz belirlemek zor olabilir.
Ve hiçbir şey çakışmıyor olsa bile, şeylerin yapılandırılabileceği tüm olası yerlere bakmak sıkıcı olabilir.

Neyse ki, Nextflow bu tüm süreci sizin için otomatikleştirebilen `config` adlı kullanışlı bir yardımcı araç içerir.

`config` aracı mevcut çalışma dizininizdeki tüm içerikleri keşfedecek, tüm yapılandırma dosyalarını toplayacak ve Nextflow'un iş akışını çalıştırmak için kullanacağı tamamen çözümlenmiş yapılandırmayı üretecektir.
Bu, hiçbir şey başlatmak zorunda kalmadan hangi ayarların kullanılacağını öğrenmenizi sağlar.

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

Bu, komut satırında fazladan bir şey belirtmezseniz elde ettiğiniz temel yapılandırmayı gösterir.

#### 6.3.2. Belirli ayarlar etkinleştirilmiş yapılandırmayı çözümleyin

Komut satırı parametreleri sağlarsanız, örn. bir veya daha fazla profili etkinleştirme veya bir parametre dosyası yükleme, komut bunları da dikkate alacaktır.

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

Bu, birden fazla yapılandırma katmanı içeren karmaşık projeler için özellikle faydalıdır.

### Özet

Çalışma zamanında minimum güçlükle önceden ayarlanmış bir yapılandırma seçmek için profilleri nasıl kullanacağınızı biliyorsunuz.
Daha genel olarak, iş akışı yürütmelerinizi farklı hesaplama platformlarına uyacak şekilde yapılandırmayı ve analizlerinizin tekrarlanabilirliğini artırmayı biliyorsunuz.

### Sırada ne var?

Kutlayın ve kendinize büyük bir övgü verin! İlk Nextflow geliştirici kursunuzu tamamladınız.

Ne öğrendiğinizi gözden geçirmek ve sırada ne olduğunu öğrenmek için son [kurs özeti](./next_steps.md)'ne gidin.

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
Aynı parametre hem yapılandırma dosyasında hem de komut satırında ayarlandığında hangisi öncelik alır?
- [ ] Yapılandırma dosyası değeri
- [x] Komut satırı değeri
- [ ] İlk karşılaşılan değer
- [ ] Hiçbiri; bir hataya neden olur

Daha fazla bilgi: [1.1. Varsayılan değerleri `nextflow.config`'e taşıyın](#11-move-default-values-to-nextflowconfig)
</quiz>

<quiz>
Aynı yapılandırmada hem Docker hem de Conda etkinleştirilebilir mi?
- [x] Evet, Nextflow süreç yönergelerine bağlı olarak her ikisini de kullanabilir
- [ ] Hayır, aynı anda yalnızca biri etkinleştirilebilir
- [ ] Evet, ama sadece profillerde
- [ ] Hayır, birbirini dışlarlar
</quiz>

<quiz>
Hem Docker hem de Conda etkinse ve bir süreçte her iki yönerge de varsa, hangisi öncelikli?
- [x] Docker (konteynerler)
- [ ] Conda
- [ ] İlk tanımlanan
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
Hangi komut satırı seçeneği kaynak kullanım raporu oluşturur?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Daha fazla bilgi: [5.1. Kaynak kullanım raporu oluşturmak için iş akışını çalıştırın](#51-run-the-workflow-to-generate-a-resource-utilization-report)
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
Nextflow'daki varsayılan yürütücü nedir?
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
Profiller ne için kullanılabilir? (Uygulanabilen tümünü seçin)
- [x] Altyapıya özgü ayarları tanımlamak
- [x] Farklı ortamlar için kaynak sınırları ayarlamak
- [x] Test parametreleri sağlamak
- [ ] Yeni süreçler tanımlamak

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
