# Bölüm 5: Merhaba Konteynerlar

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube kanalında görün.

:green_book: Video metni [burada](./transcripts/05_hello_containers.md) mevcuttur.
///

Bu eğitim kursunun 1-4. Bölümlerinde, bazı metinleri işleyebilen, birden fazla girdi varsa yürütmeyi paralel hale getirebilen ve sonuçları daha fazla işleme için toplayabilen basit bir iş akışı oluşturmak için Nextflow'un temel yapı taşlarını nasıl kullanacağınızı öğrendiniz.

Ancak, ortamınızda varsayılan olarak bulunan temel UNIX araçlarıyla sınırlıydınız.
Gerçek dünya görevleri genellikle varsayılan olarak dahil edilmeyen çeşitli araçlar ve paketler gerektirir.
Tipik olarak, bu araçları kurmanız, bağımlılıklarını yönetmeniz ve herhangi bir çakışmayı çözmeniz gerekir.

Bunların hepsi çok sıkıcı ve can sıkıcıdır, bu yüzden size bu sorunu çok daha rahat bir şekilde çözmek için **konteynerları** nasıl kullanacağınızı göstereceğiz.

**Konteyner**, kod, sistem kütüphaneleri ve ayarlar dahil olmak üzere bir uygulamayı çalıştırmak için gereken her şeyi içeren bir konteyner **imajından** oluşturulan hafif, bağımsız, çalıştırılabilir bir yazılım birimidir.
Tahmin edebileceğiniz gibi, bu boru hatlarınızı daha tekrarlanabilir hale getirmek için çok yardımcı olacaktır.

Bu eğitimi [Docker](https://www.docker.com/get-started/) kullanarak vereceğimizi unutmayın, ancak Nextflow'un [diğer birçok konteyner teknolojisini](https://nextflow.io/docs/latest/container.html) de desteklediğini aklınızda bulundurun.

??? info "Bu bölümden nasıl başlanır"

    Bu kursun bu bölümü, [Hello Nextflow](./index.md) kursunun 1-4. Bölümlerini tamamladığınızı ve eksiksiz çalışan bir boru hattınız olduğunu varsayar.

    Kursa bu noktadan başlıyorsanız, `modules` dizinini çözümlerden kopyalamanız gerekecektir:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Isınma: `hello-containers.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-containers.nf` iş akışı betiğini kullanacağız.
Bu betik, bu eğitim kursunun 4. Bölümünde üretilen betiğe eşdeğerdir, ancak çıktı hedeflerini değiştirdik:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Her şeyin çalıştığından emin olmak için, herhangi bir değişiklik yapmadan önce betiği bir kez çalıştırın:

```bash
nextflow run hello-containers.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Daha önce olduğu gibi, çıktı dosyalarını `output` bloğunda belirtilen dizinde (`results/hello_containers/`) bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Bu sizin için çalıştıysa, konteynerları nasıl kullanacağınızı öğrenmeye hazırsınız.

---

## 1. Bir konteyneri 'manuel olarak' kullanın

Yapmak istediğimiz şey, iş akışımıza yürütme için bir konteyner kullanacak bir adım eklemektir.

Ancak, önce Nextflow'da kullanmaya başlamadan önce konteynerların ne olduğuna dair anlayışınızı sağlamlaştırmak için bazı temel kavramlar ve işlemler üzerinden geçeceğiz.

### 1.1. Konteyner imajını çekin

Bir konteyner kullanmak için, genellikle bir konteyner kayıt defterinden bir konteyner imajı indirirsiniz veya _çekersiniz_ ve ardından bir konteyner örneği oluşturmak için konteyner imajını çalıştırırsınız.

Genel sözdizimi şu şekildedir:

```bash title="Sözdizimi"
docker pull '<container>'
```

`docker pull` kısmı, konteyner sistemine bir depodan bir konteyner imajı çekmesi talimatıdır.

`'<container>'` kısmı, konteyner imajının URI adresidir.

Örnek olarak, rastgele metin girdilerini eğlenceli bir şekilde görüntülemek için ASCII sanatı üreten `cowsay` adlı bir aracın python uygulaması olan [cowpy](https://github.com/jeffbuttars/cowpy)'yi içeren bir konteyner imajı çekelim.

```txt title="Örnek"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

Yayınlanmış konteynerleri bulabileceğiniz çeşitli depolar vardır.
Bu Docker konteyner imajını `cowpy` Conda paketinden oluşturmak için [Seqera Containers](https://seqera.io/containers/) hizmetini kullandık: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Tam çekme komutunu çalıştırın:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Komut çıktısı"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
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
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

İmajı daha önce hiç indirmediyseniz, bu işlemin tamamlanması bir dakika sürebilir.
İşlem tamamlandığında, konteyner imajının yerel bir kopyasına sahip olursunuz.

### 1.2. Tek seferlik bir komut olarak `cowpy` çalıştırmak için konteyneri kullanın

İnsanların konteynerleri kullanmasının çok yaygın bir yolu, onları doğrudan, _yani_ etkileşimsiz olarak çalıştırmaktır.
Bu, tek seferlik komutları çalıştırmak için harikadır.

Genel sözdizimi şu şekildedir:

```bash title="Sözdizimi"
docker run --rm '<container>' [araç komutu]
```

`docker run --rm '<container>'` kısmı, konteyner sistemine bir konteyner imajından bir konteyner örneği başlatması ve içinde bir komut yürütmesi talimatıdır.
`--rm` bayrağı, sisteme komut tamamlandıktan sonra konteyner örneğini kapatmasını söyler.

`[araç komutu]` sözdizimi, kullandığınız araca ve konteynerin nasıl kurulduğuna bağlıdır.
Sadece `cowpy` ile başlayalım.

Tamamen birleştirilmiş halde, konteyner yürütme komutu şöyle görünür; devam edin ve çalıştırın.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Komut çıktısı"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Sistem konteyneri başlattı, `cowpy` komutunu parametreleriyle çalıştırdı, çıktıyı konsola gönderdi ve son olarak konteyner örneğini kapattı.

### 1.3. `cowpy`'yi etkileşimli olarak çalıştırmak için konteyneri kullanın

Ayrıca bir konteyneri etkileşimli olarak çalıştırabilirsiniz, bu size konteyner içinde bir kabuk istemi verir ve komutla oynamanıza olanak tanır.

#### 1.3.1. Konteyneri başlatın

Etkileşimli olarak çalıştırmak için, `docker run` komutuna `-it` ekliyoruz.
İsteğe bağlı olarak, komutun sonuna _örneğin_ `/bin/bash` ekleyerek konteyner içinde kullanmak istediğimiz kabuğu belirtebiliriz.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

İsteminizin `(base) root@b645838b3314:/tmp#` gibi bir şeye değiştiğini fark edin, bu da artık konteyner içinde olduğunuzu gösterir.

Dosya sisteminin kökünden dizin içeriğini listelemek için `ls /` komutunu çalıştırarak bunu doğrulayabilirsiniz:

```bash
ls /
```

??? abstract "Komut çıktısı"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Burada `tree` yerine `ls` kullanıyoruz çünkü `tree` yardımcı programı bu konteynerde mevcut değil.
Konteyner içindeki dosya sisteminin ana sistem dosya sisteminizden farklı olduğunu görebilirsiniz.

Az önce yaptığımız şeyin bir sınırlaması, konteynerin varsayılan olarak ana sistemden tamamen izole olmasıdır.
Bu, açıkça izin vermediğiniz sürece konteynerin ana sistemdeki herhangi bir dosyaya erişemeyeceği anlamına gelir.

Bunu nasıl yapacağınızı birazdan göstereceğiz.

#### 1.3.2. İstenen araç komut(lar)ını çalıştırın

Artık konteyner içinde olduğunuza göre, `cowpy` komutunu doğrudan çalıştırabilir ve ona bazı parametreler verebilirsiniz.
Örneğin, araç belgeleri `-c` ile karakteri ('cowacter') değiştirebileceğimizi söylüyor.

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

Şimdi çıktı, `-c tux` parametresini belirttiğimiz için varsayılan inek yerine Linux pengueni Tux'u gösteriyor.

Konteyner içinde olduğunuz için, Docker komutlarıyla uğraşmak zorunda kalmadan, girdi parametrelerini değiştirerek `cowpy` komutunu istediğiniz kadar çalıştırabilirsiniz.

!!! Tip

    Farklı bir karakter seçmek için '-c' bayrağını kullanın, bunlar arasında:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Bu güzel. Daha da güzel olan, `greetings.csv` dosyamızı buna girdi olarak besleyebilmek olurdu.
Ancak dosya sistemine erişimimiz olmadığı için yapamıyoruz.

Hadi bunu düzeltelim.

#### 1.3.3. Konteynerden çıkın

Konteynerden çıkmak için, istemde `exit` yazabilir veya ++ctrl+d++ klavye kısayolunu kullanabilirsiniz.

```bash
exit
```

İsteminiz artık konteyneri başlatmadan önceki haline dönmüş olmalıdır.

#### 1.3.4. Verileri konteynere bağlayın

Daha önce belirtildiği gibi, konteyner varsayılan olarak ana sistemden izole edilmiştir.

Konteynerin ana dosya sistemine erişmesine izin vermek için, aşağıdaki sözdizimini kullanarak ana sistemden konteynere bir **birim** **bağlayabilirsiniz**:

```bash title="Sözdizimi"
-v <dış_yol>:<iç_yol>
```

Bizim durumumuzda `<dış_yol>` mevcut çalışma dizini olacak, bu yüzden sadece bir nokta (`.`) kullanabiliriz ve `<iç_yol>` uydurduğumuz bir takma addır; buna `/my_project` diyelim (iç yol mutlak olmalıdır).

Bir birim bağlamak için, yolları değiştirir ve birim bağlama argümanını docker run komutuna şu şekilde ekleriz:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Bu, mevcut çalışma dizinini konteyner içinde `/my_project` altında erişilebilir olacak bir birim olarak bağlar.

`/my_project` içeriğini listeleyerek çalıştığını kontrol edebilirsiniz:

```bash
ls /my_project
```

??? success "Komut çıktısı"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Artık `data/` altındaki `greetings.csv` dosyası dahil olmak üzere çalışma dizininin içeriğini konteyner içinden görebilirsiniz.

Bu, dosya sisteminizin o kısmına erişmek için kullanabileceğiniz konteyner duvarından bir tünel oluşturdu.

#### 1.3.5. Bağlanan verileri kullanın

Artık çalışma dizinini konteynere bağladığımıza göre, `greetings.csv` dosyasının içeriğini görüntülemek için `cowpy` komutunu kullanabiliriz.

Bunu yapmak için, CSV dosyasının içeriğini `cowpy` komutuna aktarmak için `cat /my_project/data/greetings.csv | ` kullanacağız.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Komut çıktısı"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
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

Bu, örnek selamlamalarımızı tekrarlayan bir hindinin istenen ASCII sanatını üretir!
Ancak burada hindi sadece selamlamalar yerine tam satırları tekrarlıyor.
Nextflow iş akışımızın daha iyi bir iş çıkaracağını zaten biliyoruz!

Bu komutla oynamaktan çekinmeyin.
İşiniz bittiğinde, daha önce olduğu gibi konteynerden çıkın:

```bash
exit
```

Kendinizi normal kabuğunuzda bulacaksınız.

### Özet

Bir konteyneri nasıl çekeceğinizi ve tek seferlik veya etkileşimli olarak nasıl çalıştıracağınızı biliyorsunuz. Ayrıca verilerinizi konteynerinizin içinden nasıl erişilebilir hale getireceğinizi de biliyorsunuz, bu da sisteminize herhangi bir yazılım yüklemek zorunda kalmadan ilgilendiğiniz herhangi bir aracı gerçek veriler üzerinde denemenize olanak tanır.

### Sırada ne var?

Nextflow süreçlerinin yürütülmesi için konteynerleri nasıl kullanacağınızı öğrenin.

---

## 2. Nextflow'da konteynerleri kullanın

Nextflow, bilgi işlem ortamınızda yüklü olmayan araçları çalıştırmanıza olanak tanımak için süreçleri konteynerler içinde çalıştırmak için yerleşik desteğe sahiptir.
Bu, süreçlerinizi çalıştırmak için istediğiniz herhangi bir konteyner imajını kullanabileceğiniz ve Nextflow'un imajı çekme, verileri bağlama ve süreci içinde çalıştırma işlemlerini halledeceği anlamına gelir.

Bunu göstermek için, geliştirdiğimiz boru hattına `collectGreetings` adımından sonra bir `cowpy` adımı ekleyeceğiz.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Bir `cowpy` modülü yazın

İlk olarak, `cowpy` süreç modülünü oluşturalım.

#### 2.1.1. Yeni modül için bir dosya taslağı oluşturun

Modül için `cowpy.nf` adında boş bir dosya oluşturun.

```bash
touch modules/cowpy.nf
```

Bu bize süreç kodunu koyacak bir yer verir.

#### 2.1.2. `cowpy` süreç kodunu modül dosyasına kopyalayın

`cowpy` sürecimizi daha önce yazdığımız diğer süreçlere göre modelleyebiliriz.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// cowpy ile ASCII sanatı oluştur
process cowpy {

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

Süreç, selamlamaları içeren bir `input_file` ve bir `character` değeri bekler.

Çıktı, `cowpy` aracı tarafından oluşturulan ASCII sanatını içeren yeni bir metin dosyası olacaktır.

### 2.2. İş akışına cowpy ekleyin

Şimdi modülü içe aktarmamız ve süreci çağırmamız gerekiyor.

#### 2.2.1. `cowpy` sürecini `hello-containers.nf` dosyasına aktarın

İçe aktarma bildirimini iş akışı bloğunun üstüne ekleyin ve uygun şekilde doldurun.

=== "Sonra"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Önce"

    ```groovy title="hello-containers.nf" linenums="3"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Artık `cowpy` modülü iş akışında kullanılabilir.

#### 2.2.2. İş akışına `cowpy` sürecine bir çağrı ekleyin

`cowpy()` sürecini, hatırlayabileceğiniz gibi iki çıktı üreten `collectGreetings()` sürecinin çıktısına bağlayalım:

- `collectGreetings.out.outfile` çıktı dosyasını içerir <--_istediğimiz şey_
- `collectGreetings.out.report` toplu iş başına selamlama sayısını içeren rapor dosyasını içerir

İş akışı bloğunda, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
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
    ```

=== "Önce"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // bir CSV dosyasından girdiler için bir kanal oluştur
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

Hangi karakterin selamlamaları söylemesini istediğimizi belirtmek için yeni bir CLI parametresi olan `params.character` bildirdiğimize dikkat edin.

#### 2.2.3. `character` parametresini `params` bloğuna ekleyin

Bu teknik olarak isteğe bağlıdır ancak önerilen uygulamadır ve bu arada karakter için varsayılan bir değer belirleme fırsatıdır.

=== "Sonra"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Boru hattı parametreleri
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Önce"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Boru hattı parametreleri
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Artık tembel olabilir ve komut satırlarımızda karakter parametresini yazmayı atlayabiliriz.

#### 2.2.4. İş akışı çıktılarını güncelleyin

`cowpy` sürecinin çıktısını yayınlamak için iş akışı çıktılarını güncellememiz gerekiyor.

##### 2.2.4.1. `publish:` bölümünü güncelleyin

`workflow bloğunda`, aşağıdaki kod değişikliğini yapın:

=== "Sonra"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Önce"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

`cowpy` süreci yalnızca bir çıktı üretir, bu nedenle `.out` ekleyerek her zamanki şekilde ona başvurabiliriz.

Ancak şimdilik, iş akışı düzeyindeki çıktıları güncellemeyi bitirelim.

##### 2.2.4.2. `output` bloğunu güncelleyin

Son `cowpy_art` çıktısını `output` bloğuna eklememiz gerekiyor. Bu arada, artık boru hattımız tamamlandığı ve gerçekten önemsediğimiz çıktıların ne olduğunu bildiğimiz için yayınlama hedeflerini de düzenleyelim.

`output` bloğunda, aşağıdaki kod değişikliklerini yapın:

=== "Sonra"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Önce"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Artık yayınlanan çıktılar biraz daha düzenli olacak.

#### 2.2.5. İş akışını çalıştırın

Özetlemek gerekirse, hedeflediğimiz şey bu:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Sizce işe yarayacak mı?

Temiz bir sayfa için önceki yayınlanmış çıktıları silelim ve iş akışını `-resume` bayrağıyla çalıştıralım.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Komut çıktısı (netlik için düzenlenmiş)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

Ah hayır, bir hata var!
`error exit status (127)` tarafından verilen hata kodu, istediğimiz çalıştırılabilir dosyanın bulunamadığı anlamına gelir.

Bu mantıklı, çünkü `cowpy` aracını çağırıyoruz ancak henüz bir konteyner belirtmedik (hata).

### 2.3. `cowpy` sürecini çalıştırmak için bir konteyner kullanın

Bir konteyner belirtmemiz ve Nextflow'a `cowpy()` süreci için kullanmasını söylememiz gerekiyor.

#### 2.3.1. `cowpy` için bir konteyner belirtin

Bu eğitimin ilk bölümünde doğrudan kullandığımız aynı imajı kullanabiliriz.

`cowpy.nf` modülünü düzenleyerek süreç tanımına `container` yönergesini şu şekilde ekleyin:

=== "Sonra"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
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

=== "Önce"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

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

Bu, Nextflow'a _Docker kullanımı etkinleştirilmişse_, süreci yürütmek için burada belirtilen konteyner imajını kullanması gerektiğini söyler.

#### 2.3.2. `nextflow.config` dosyası aracılığıyla Docker kullanımını etkinleştirin

_'Docker kullanımı etkinleştirilmişse'_ dediğimize dikkat edin. Varsayılan olarak etkin değildir, bu nedenle Nextflow'a Docker kullanmasına izin verildiğini söylememiz gerekiyor.
Bu amaçla, bu kursun bir sonraki ve son bölümünün (Bölüm 6) konusunu biraz önceden ele alacağız, bu bölüm yapılandırmayı kapsar.

Nextflow'un iş akışı yürütmesini yapılandırmak için sunduğu ana yollardan biri, bir `nextflow.config` dosyası kullanmaktır.
Böyle bir dosya mevcut dizinde bulunduğunda, Nextflow onu otomatik olarak yükleyecek ve içerdiği herhangi bir yapılandırmayı uygulayacaktır.

Docker'ı açıkça devre dışı bırakan tek satırlık bir kodla bir `nextflow.config` dosyası sağladık: `docker.enabled = false`.

Şimdi, Docker'ı etkinleştirmek için bunu `true` olarak değiştirelim:

=== "Sonra"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Önce"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip

    `-with-docker <container>` parametresini kullanarak komut satırından, çalıştırma başına Docker yürütmesini etkinleştirmek mümkündür.
    Ancak, bu bize tüm iş akışı için yalnızca bir konteyner belirtmemize izin verirken, az önce size gösterdiğimiz yaklaşım süreç başına farklı bir konteyner belirtmemize olanak tanır.
    Bu, modülerlik, kod bakımı ve tekrarlanabilirlik için daha iyidir.

#### 2.3.3. Docker etkinken iş akışını çalıştırın

İş akışını `-resume` bayrağıyla çalıştırın:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

Bu sefer gerçekten işe yarıyor!
Her zamanki gibi iş akışı çıktılarını ilgili sonuçlar dizininde bulabilirsiniz, ancak bu sefer sadece rapor ve son çıktı üst düzeyde ve tüm ara dosyalar bir alt dizine itilmiş olarak biraz daha düzenli bir şekilde organize edilmişlerdir.

??? abstract "Dizin içeriği"

    ```console
    results/hello_containers/
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

Son ASCII sanat çıktısı, `cowpy-COLLECTED-batch-output.txt` adı altında `results/hello_containers/` dizinindedir.

??? abstract "Dosya içeriği"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
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

Ve işte burada, selamlamaları istediğimiz gibi söyleyen güzel hindimiz.

#### 2.3.4. Nextflow'un konteynerleştirilmiş görevi nasıl başlattığını inceleyin

Bu bölümün son bir eki olarak, Nextflow'un konteynerlerle perde arkasında nasıl çalıştığına dair biraz daha içgörü elde etmek için `cowpy` süreç çağrılarından birinin çalışma alt dizinine bir göz atalım.

`cowpy` süreci için çalışma alt dizininin yolunu bulmak için `nextflow run` komutunuzdan çıktıyı kontrol edin.
Yukarıda gösterilen çalıştırma için aldığımız şeye bakarsak, `cowpy` süreci için konsol günlüğü satırı `[98/656c6c]` ile başlar.
Bu, şu kısaltılmış dizin yoluna karşılık gelir: `work/98/656c6c`.

Bu dizinde, Nextflow'un boru hattını yürütme sürecinde sizin adınıza çalıştırdığı tüm komutları içeren `.command.run` dosyasını bulacaksınız.

??? abstract "Dosya içeriği"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
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
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # girdi dosyalarını hazırla
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
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
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
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
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Gördüğünüz gibi, Nextflow süreç çağrısını başlatmak için `docker run` komutunu kullanıyor.
Ayrıca ilgili çalışma alt dizinini konteynere bağlar, konteyner içindeki çalışma dizinini buna göre ayarlar ve `.command.sh` dosyasındaki şablonlu bash betiğimizi çalıştırır.

İlk bölümde manuel olarak yapmak zorunda kaldığımız tüm zor işleri Nextflow perde arkasında bizim için yapıyor!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Özet

Nextflow'da süreçleri çalıştırmak için konteynerleri nasıl kullanacağınızı biliyorsunuz.

### Sırada ne var?

Bir mola verin!

Hazır olduğunuzda, boru hattınızın yürütülmesini altyapınıza uyacak şekilde nasıl yapılandıracağınızı ve girdilerin ve parametrelerin yapılandırmasını nasıl yöneteceğinizi öğrenmek için [**Bölüm 6: Merhaba Config**](./06_hello_config.md) bölümüne geçin.

Bu en son bölüm ve ardından bu kursla işiniz bitmiş olacak!

---

## Sınav

<quiz>
Konteyner nedir?
- [ ] Bir tür sanal makine
- [ ] Bir dosya sıkıştırma formatı
- [x] Bir uygulamayı çalıştırmak için gereken her şeyi içeren hafif, bağımsız, çalıştırılabilir bir birim
- [ ] Bir ağ protokolü
</quiz>

<quiz>
Konteyner imajı ile konteyner örneği arasındaki fark nedir?
- [ ] Aynı şeylerdir
- [x] İmaj bir şablondur; örnek, o imajdan oluşturulan çalışan bir konteynerdir
- [ ] Örnek bir şablondur; imaj çalışan bir konteynerdir
- [ ] İmajlar Docker içindir; örnekler Singularity içindir
</quiz>

<quiz>
`docker run` komutunda `-v` bayrağı ne yapar?
- [ ] Ayrıntılı çıktıyı etkinleştirir
- [ ] Konteyneri doğrular
- [x] Ana sistemden konteynere bir birim bağlar
- [ ] Konteynerin sürümünü belirtir

Daha fazla bilgi: [1.3.4. Verileri konteynere bağlayın](#134-mount-data-into-the-container)
</quiz>

<quiz>
Konteynerler kullanırken neden birimleri bağlamanız gerekir?
- [ ] Konteyner performansını artırmak için
- [ ] Disk alanından tasarruf etmek için
- [x] Çünkü konteynerler varsayılan olarak ana dosya sisteminden izole edilmiştir
- [ ] Ağı etkinleştirmek için

Daha fazla bilgi: [1.3.4. Verileri konteynere bağlayın](#134-mount-data-into-the-container)
</quiz>

<quiz>
Bir Nextflow süreci için nasıl bir konteyner belirtirsiniz?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Daha fazla bilgi: [2.3.1. cowpy için bir konteyner belirtin](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
İş akışınız için Docker'ı hangi `nextflow.config` ayarı etkinleştirir?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Daha fazla bilgi: [2.3.2. `nextflow.config` dosyası aracılığıyla Docker kullanımını etkinleştirin](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
Nextflow bir süreci bir konteynerde çalıştırırken otomatik olarak ne halleder? (Tümünü seçin)
- [x] Gerekirse konteyner imajını çekme
- [x] Çalışma dizinini bağlama
- [x] Süreç betiğini konteyner içinde çalıştırma
- [x] Yürütmeden sonra konteyner örneğini temizleme

Daha fazla bilgi: [2.3.4. Nextflow'un konteynerleştirilmiş görevi nasıl başlattığını inceleyin](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
