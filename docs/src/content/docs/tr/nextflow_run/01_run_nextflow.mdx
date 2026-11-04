# Bölüm 1: Nextflow'u Çalıştırma

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu bölümde, Nextflow pipeline'larını çalıştırmanın temel kavramlarını tanıtıyoruz.
Basit bir Hello World iş akışıyla başlayıp, konteynerler kullanarak birden fazla girdiyi paralel olarak işleyen eksiksiz, çok adımlı bir pipeline'a kadar ilerliyoruz.

---

## 1. Hello World

`1-hello.nf` iş akışı, bir komut satırı argümanı aracılığıyla bir selamlama alır ve bunu bir dosyaya yazar.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. İş akışını başlatma

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

Çıktıdaki en önemli satır, süreç durum satırıdır:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Bu satır, `sayHello` sürecinin bir kez başarıyla çalıştığını gösterir.
`[6d/740edd]` öneki, görevin çalışma dizinine giden kısaltılmış bir yoldur — bununla ilgili daha fazla bilgi aşağıda verilmektedir.
Ardından gelen `Outputs:` bloğu, pipeline'ın yayımladığı her dosyayı listeler; etiketler, aşağıdaki [1.4](#14-optional-code-walkthrough) bölümünde ele alınan `output` bloğuna göre belirlenir.

### 1.2. Çıktıyı bulma

Bu iş akışı, çıktısını bir `results` dizinine yayımlamak üzere yapılandırılmıştır.
Çalıştırmanın ardından çıktıyı orada bulabilirsiniz:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Dosyayı açarak `Hello World!` içerdiğini doğrulayın.

### 1.3. `work/` dizinini keşfetme

Arka planda Nextflow, `work/` adlı bir dizin içinde her süreç çağrısı için benzersiz bir görev dizini oluşturur.
Konsol çıktısında gösterilen hash (`[6d/740edd]`), bu dizine giden yoldur.

```bash
ls work/6d/740edd*
```

İçeride çıktı dosyasının yanı sıra birkaç gizli log dosyası bulacaksınız:

- **`.command.sh`**: Nextflow'un çalıştırdığı tam komut
- **`.command.out`** / **`.command.err`**: süreçten gelen stdout ve stderr
- **`.command.log`**: birleşik log çıktısı
- **`.exitcode`**: sürecin çıkış kodu

`.command.sh` dosyası, hata ayıklama sırasında özellikle kullanışlıdır — tam olarak neyin çalıştırıldığını gösterir.

### 1.4. İsteğe bağlı: Kod incelemesi

Yalnızca pipeline çalıştırmak istiyorsanız kodu anlamak zorunlu değildir; ancak merak ediyorsanız incelemeye değer.

??? optional "Bu alıştırmayla ilişkili kodu keşfetmek için tıklayın"

    `1-hello.nf` dosyasını açalım ve ana bileşenlerine bakalım.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

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

    Şunları görüyoruz:

    - bir `process` modülüne işaret eden `include` ifadesi
    - pipeline parametrelerini tanımlayan `params` bloğu
    - yapılacak işi tanımlayan `workflow` bloğu
    - çıktılarla ne yapılacağını belirten `output` bloğu

    Her birine sırayla bakalım.

    ### `process` modülü

    `include` ifadesi, Nextflow'a ayrı bir kod dosyasından `sayHello` adlı bir şeyi yüklemesini söyler.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    O dosyada, `sayHello` adlı bir sürecin tanımını buluruz:

    ```groovy title="modules/sayHello.nf" linenums="4"
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

    Bir **process**, pipeline'daki tek bir adımı tanımlar.
    Girdilerini, çıktılarını ve çalıştırılacak betiği bildirir.
    `val` niteleyicisi, girdinin düz bir değer (string, sayı vb.) olduğu anlamına gelir.
    `path` niteleyicisi ise çıktının bir dosya yolu olduğu anlamına gelir.

    Süreç tanımını ana iş akışı dosyasına yazmak mümkündür; ancak bunları ayrı modül dosyalarında tutmak yeniden kullanılabilirliklerini sağlar: aynı modül birden fazla iş akışı betiği tarafından içe aktarılabilir.

    ### `params` bloğu

    `params` bloğu, iş akışının kabul ettiği komut satırı parametrelerini bildirir:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Burada bildirilen her parametre, komut satırında çift tire (`--input`) ile kullanılabilir hale gelir.
    Desteklenen türler şunlardır: `String`, `Integer`, `Float`, `Boolean` ve `Path`.

    !!! tip "İpucu"

        İş akışı parametreleri, Nextflow'un kendi CLI bayraklarından (tek tire kullanan, örn. `-resume`) ayırt edilmek için her zaman iki tire (`--input`) kullanır.

    ### `workflow` bloğu

    **workflow** bloğu, veri akışı mantığını tanımlar: hangi süreçlerin çalıştırılacağını ve hangi sırayla.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // bir selamlama yayınla
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Burada yalnızca bir süreç çağrıldığından oldukça basittir; daha gerçekçi örnekleri ilerleyen bölümlerde ele alacağız.

    `main:` bölümü, `sayHello` sürecini `--input` değeriyle çağırır.
    `publish:` bölümü ise hangi çıktıların results dizinine kopyalanacağını listeler.

    ### `output` bloğu

    Dosyanın alt kısmındaki `output` bloğu, hedef yolu ve kopyalama modunu belirtir.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Her adlandırılmış giriş, iş akışındaki bir `publish:` etiketine karşılık gelir ve onu `results/` altındaki bir alt dizine eşler.

### Özetle

Bir Nextflow pipeline'ını nasıl çalıştıracağınızı ve çıktılarını nasıl bulacağınızı öğrendiniz; ayrıca işin `work/` altındaki görev dizinlerinde yürütüldüğünü öğrendiniz.

### Sırada ne var?

Nextflow'un birden fazla girdiyi nasıl verimli bir şekilde işlediğini keşfedin.

---

## 2. Birden fazla girdiyi işleme

Gerçek dünya pipeline'ları genellikle yalnızca bir değil, pek çok veri parçasını işler.
`2-inputs.nf` iş akışı bir CSV dosyasından okur ve her satır için `sayHello`'yu paralel olarak bir kez çalıştırır.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Önce iş akışını çalıştıralım, ardından Nextflow'un bu birden fazla girdiyi işlemek için hangi mekanizmayı kullandığına bakalım.

### 2.1. İş akışını çalıştırma

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

`3 of 3` ifadesi, `sayHello` sürecinin CSV'deki her satır için bir kez olmak üzere üç kez çağrıldığını gösterir.

`results` dizininde artık her selamlama için bir tane olmak üzere üç çıktı dosyası görmelisiniz:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Her birinin bir selamlama içerdiğini doğrulamak için çıktı dosyalarından herhangi birini açın.

Yukarıdaki özet çıktı, `sayHello` için tek bir özet satırı gösterse de Nextflow aslında arka planda üç ayrı görev yürütmesi başlattı; CSV'deki her satır için bir tane ve makinenizin kaynakları elverdiği anda bunları paralel olarak çalıştırdı.

[1.3](#13-explore-the-work-directory) bölümünde incelediğiniz tek görevde olduğu gibi, bu üç yürütmenin her biri `work/` altında kendi görev dizinini alır ve diğerlerinden tamamen yalıtılmış olur:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Her `.command.sh` yalnızca o tek selamlama için komutu içerir:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Bu yalıtım, paralel yürütmeyi güvenli kılan şeydir: aynı anda çalışan üç görev hiçbir zaman bir çalışma dizinini paylaşmaz; dolayısıyla bir görevin yazdığı şey, aynı ada sahip dosyalar üretilse bile başka bir görevin yazdığıyla çakışamaz veya onu üzerine yazamaz.
Aynı zamanda `-resume`'un (bir sonraki konuda ele alınacak) bireysel görevleri bağımsız olarak önbelleğe alıp yeniden kullanabilmesinin de nedeni budur: her görevin girdileri, çıktıları ve logları tamamen kendi dizininde bulunur; görevler arasında senkronizasyonu bozabilecek hiçbir şey paylaşılmaz.

### 2.2. İş akışını `-ansi-log false` ile yeniden çalıştırma

Nextflow varsayılan olarak çıktıyı her süreç için tek bir özet satırına sıkıştırır.
Her süreç çağrısının ayrı ayrı listelenmesini görmek için `-ansi-log false` ekleyin:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Bu, üç süreç çağrısının tamamını ve her biri için oluşturulan benzersiz work alt dizinini gösterir.

### 2.3. Tamamlanan işi atlamak için `-resume` kullanma

Şimdi iki selamlama daha ekleyen genişletilmiş girdi dosyasına geçin ve komut satırına `-resume` ekleyin:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Komut çıktısı"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow yalnızca iki yeni girdiyi çalıştırdı.
Önceki çalıştırmada işlenen üç selamlama önbelleğe alındı ve otomatik olarak yeniden kullanıldı.

Bu özellik, çok adımlı bir pipeline'da zaten başarıyla çalıştırılmış adımların yeniden yürütülmesini atlamak için de işe yarar.
Örneğin, bir pipeline çalıştırması sistem hatası nedeniyle kesintiye uğradıysa ya da geliştirme aşamasındaki bir pipeline'a yeni adımlar eklediyseniz.

`-resume` özelliği, uzun pipeline'larda özellikle değerlidir; hatadan kurtarma, kritik zaman ve kaynak tasarrufu sağlayabilir.

### 2.4. İsteğe bağlı: Kod incelemesi

Yalnızca pipeline çalıştırmak istiyorsanız kodu anlamak zorunlu değildir; ancak merak ediyorsanız incelemeye değer.

??? optional "Bu alıştırmayla ilişkili kodu keşfetmek için tıklayın"

    `2-inputs.nf` dosyasındaki temel değişiklik, iş akışının `main:` bölümündedir:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // bir selamlama yayınla
        sayHello(greeting_ch)
    ```

    Burada gördüğünüz şeye **kanal** denir: işlemleri paralelleştirmeyi kolaylaştıracak şekilde girdi verilerini işleyen bir kuyruk yapısıdır.

    - `channel.fromPath(params.input)`, `--input` ile verilen dosya yolundan bir kanal oluşturur
    - `.splitCsv()`, CSV'yi satırlara ayrıştırır
    - `#!groovy .map { line -> line[0] }`, her satırdan ilk sütunu çıkarır

    Sonuç, `Hello`, `Bonjour` ve `Hola` değerlerini içeren bir kanaldır.
    `sayHello(greeting_ch)`'ye geçirildiğinde, Nextflow her öğe için süreci otomatik olarak bir kez çağırır ve kaynaklar elverdiğinde bunları paralel olarak çalıştırır.

### Özetle

Bir CSV dosyasından birden fazla girdiyi paralel olarak nasıl işleyeceğinizi ve tamamlanan işi tekrarlamamak için `-resume`'u nasıl kullanacağınızı öğrendiniz.

### Sırada ne var?

Eksiksiz bir çok adımlı pipeline'ın süreçleri kanallar aracılığıyla nasıl birbirine zincirlediğini ve analiz araçlarıyla bağımlılıklarını yönetmek için konteynerleri nasıl kullanacağınızı öğrenin.

---

## 3. Çok adımlı bir pipeline çalıştırma

Şimdiye kadar tek bir süreci çalıştırdınız, ardından bir dizi girdi üzerinde paralel olarak birden fazla kez çalıştırdınız.
Gerçek pipeline'lar genellikle daha ileri gider: birkaç süreci birbirine zincirler, birinin çıktısını bir sonrakine besler ve çoğunlukla birden fazla yazılım parçasına dayanır.
`main.nf` iş akışı, bunların ikisini de eksiksiz bir pipeline'da bir araya getirir.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Her girdi selamlama dört adımın tamamından geçer: `sayHello` onu bir dosyaya yazar, `convertToUpper` metni büyük harfe dönüştürür, `collectGreetings` tüm sonuçları tek bir dosyada birleştirir ve `cowpy` birleştirilmiş çıktıdan konteynerize bir araç kullanarak ASCII sanatı oluşturur.
Nextflow bu adımları kanallarla birbirine bağlar: bir sürecin çıktısı bir sonrakinin girdisi olur; böylece veriler hazır hale geldikçe tüm zincir otomatik olarak çalışır ve her adımı elle düzenlemenize gerek kalmaz.

Bu iş akışının modüller kullandığını unutmayın: her süreç `modules/` altında kendi dosyasında tanımlanmıştır ve `main.nf` bunları satır içi tanımlamak yerine `include` ifadeleriyle içe aktarır.
Bu, her süreci kod tekrarı olmadan birden fazla iş akışında yeniden kullanılabilir kılar. Daha fazla bilgi için aşağıdaki kod inceleme bölümüne bakın.

### 3.1. İş akışını çalıştırma

Terminalinizde aşağıdaki komutu çalıştırın.

```bash
nextflow run main.nf --input data/greetings.csv
```

`character` parametresi `nextflow.config` dosyasında varsayılan olarak `turkey` değerine ayarlanmıştır; bu nedenle siz geçersiz kılmadıkça ASCII sanatı bir hindi kullanır (deneyin: `--character tux` ekleyin).

??? success "Komut çıktısı"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Dört süreç çalıştı; ancak hepsi aynı sayıda değil.
`sayHello` ve `convertToUpper` her girdi için bir kez çalıştı (3 of 3): her selamlama ayrı ayrı yazılıp büyük harfe dönüştürülmesi gerekir.
`collectGreetings` ve `cowpy` yalnızca bir kez çalıştı (1 of 1): selamlamaları birleştirmek ve ASCII sanatı oluşturmak, ancak tüm bireysel sonuçlar hazır olduğunda anlam ifade eder.
Birkaç paralel görevin daha az sayıda aşağı akış görevine beslenmesi şeklindeki bu yelpaze açılıp kapanma biçimi, gerçek pipeline'larda yaygındır.

Nextflow, bir sonraki adımı başlatmadan önce tüm adımın bitmesini beklemez.
Bir `sayHello` çıktısı hazır olur olmaz, eşleşen `convertToUpper` görevi başlayabilir; böylece farklı süreçlerdeki görevler katı gruplar halinde değil, eş zamanlı olarak çalışır.
`collectGreetings` ve `cowpy` ise beklemek zorundadır; zira her ikisi de tüm yukarı akış sonuçlarının önce mevcut olmasına bağlıdır.

`results` dizini bu yelpaze kapanmasını ve pipeline yazarının neyi nereye yayımlamayı seçtiğini yansıtır: bu yapıyı tanımlayan şeyin 1.4 bölümündeki kod incelemesindeki `output` bloğu olduğunu hatırlayın.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

Üst düzey dizin, varsayılan olarak `batch` olan `batch` parametresinden adını alır; ilerleyen alıştırmalarda değiştiğini göreceksiniz.

ASCII sanat dosyası için `cowpy-COLLECTED-batch-output.txt` dosyasını kontrol edin.

??? abstract "Dosya içeriği"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
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

[2.1](#21-run-the-workflow) bölümünde olduğu gibi, dört süreçteki bu 8 görev yürütmesinin her biri `work/` altında kendi dizinini alır ve diğerlerinden tamamen yalıtılmış olur.
`collectGreetings`, bunun neden önemli olduğunun iyi bir örneğidir: üç farklı görev dizininde bulunan `convertToUpper` görevlerinin üçünün de çıktılarına bağlıdır; bu nedenle Nextflow, `collectGreetings`'in yukarı akış görevlerinin dizinlerinden doğrudan okuması yerine, bu dosyalara sembolik bağlantıları `collectGreetings`'in kendi dizinine yerleştirir:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Her görev yalnızca ihtiyaç duyduğu belirli dosyaları görür; nereden geldiklerinden bağımsız olarak, başka bir görevin dizininin iç içeriğini asla görmez.
Tüm bir pipeline boyunca, [2.1](#21-run-the-workflow) bölümünde tek bir süreçle gördüğünüz aynı yalıtım, Nextflow'un her süreçteki her görevi eş zamanlı ve güvenli bir şekilde çalıştırmasını sağlar.

!!! note "Not"

    `cowpy` adımı, yerel olarak yüklü yazılıma güvenmek yerine bir Docker konteyneri içinde çalışır.
    Bir konteyner, bir uygulamayı çalışması için gereken her şeyle birlikte paketler; böylece bağımlılıkları kendiniz yükleyip yönetmek zorunda kalmazsınız ve pipeline, konteyneri çalıştırabilen herhangi bir makinede aynı şekilde davranır.
    Nextflow, konteynerler için alternatif olarak Conda'yı da destekler; aralarında nasıl geçiş yapılacağını öğrenmek için [Bölüm 2](./02_configure_pipeline.md)'ye bakın.

### 3.2. İsteğe bağlı: Kod incelemesi

Yalnızca pipeline çalıştırmak istiyorsanız kodu anlamak zorunlu değildir; ancak merak ediyorsanız incelemeye değer.

??? optional "Bu alıştırmayla ilişkili kodu keşfetmek için tıklayın"

    ### Verinin bir adımdan diğerine nasıl aktığı

    Her süreç, çıktı kanalını bir sonrakine geçirir:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // CSV dosyasından girdiler için bir kanal oluştur
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    `processName.out` kalıbı, bir sürecin çıktı kanalına atıfta bulunur.

    `.collect()` operatörü, `convertToUpper`'dan gelen tüm bireysel çıktıları `collectGreetings`'e geçirmeden önce tek bir kanal öğesinde toplar.

    ### Süreç modüllerini kullanma

    `main.nf` doğrudan herhangi bir süreç kodu tanımlamaz.
    Bunun yerine, her süreci `modules/` altındaki kendi dosyasından içe aktarır:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Her modül dosyası, [1.4](#14-optional-code-walkthrough) bölümündeki `sayHello` modülüyle aynı şekilde yapılandırılmış tek bir süreç tanımı içerir.
    Süreçleri ayrı dosyalarda tutmak, kod tekrarı olmadan birden fazla iş akışında yeniden kullanılabilir olmalarını sağlar.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Konteynerize yazılım kullanma

    `cowpy` süreci, modül dosyasında belirtilen bir Docker konteyneri içinde çalışır:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
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

    Nextflow imajı otomatik olarak çeker, betiği konteyner içinde çalıştırır ve ardından temizlik yapar.
    Docker, bu proje için `nextflow.config` dosyasında etkinleştirilmiştir:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Bu tek satır, pipeline'daki konteyner belirtilmiş herhangi bir süreç için Docker'ı etkinleştirir.

### Özetle

Konteynerize bir araç kullanarak birden fazla girdiyi paralel olarak işleyen eksiksiz bir çok adımlı pipeline çalıştırdınız.

### Sırada ne var?

`nextflow.config` kullanarak pipeline davranışını nasıl yapılandıracağınızı öğreneceğiniz [Bölüm 2](./02_configure_pipeline.md)'ye geçin.

---

## Özet

Bu bölümde şunları öğrendiniz:

- Bir Nextflow iş akışını çalıştırmak ve çıktılarını bulmak
- `work/` dizinini ve log dosyalarını keşfetmek
- Bir CSV dosyasından birden fazla girdiyi paralel olarak işlemek
- Yeni girdiler eklerken tamamlanan işi atlamak için `-resume` kullanmak
- Konteynerize bir araç kullanan çok adımlı bir pipeline çalıştırmak
