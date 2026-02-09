# Bölüm 4: Merhaba Modüller

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi edinin ve iyileştirmeler önerin](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [Tüm oynatma listesini](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube kanalında izleyin.

:green_book: Video transkripti [burada](./transcripts/04_hello_modules.md) mevcuttur.
///

Bu bölüm, iş akışı kodunuzu nasıl organize edeceğinizi ve böylece boru hattınızın geliştirilmesini ve bakımını daha verimli ve sürdürülebilir hale getireceğinizi kapsamaktadır.
Özellikle, [**modüllerin**](https://nextflow.io/docs/latest/module.html) nasıl kullanılacağını göstereceğiz.

Nextflow'da bir **modül**, genellikle tek bir süreç tanımını kapsayan bağımsız bir kod dosyasıdır.
Bir iş akışında modül kullanmak için, iş akışı kod dosyanıza tek satırlık bir `include` ifadesi eklemeniz yeterlidir; ardından süreci normalde yaptığınız gibi iş akışına entegre edebilirsiniz.
Bu, kodun birden fazla kopyasını oluşturmadan süreç tanımlarını birden fazla iş akışında yeniden kullanmayı mümkün kılar.

İş akışımızı geliştirmeye başladığımızda, her şeyi tek bir kod dosyasında yazdık.
Şimdi süreçleri ayrı modüllere taşıyacağız.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Bu, kodumuzun daha paylaşılabilir, esnek ve sürdürülebilir olmasını sağlayacaktır.

??? info "Bu bölümden nasıl başlanır"

    Kursun bu bölümü, [Hello Nextflow](./index.md) kursunun 1-3. Bölümlerini tamamladığınızı varsayar, ancak bu bölümlerde ele alınan temel konularda rahat hissediyorsanız, özel bir şey yapmadan buradan başlayabilirsiniz.

---

## 0. Isınma: `hello-modules.nf` dosyasını çalıştırın

Başlangıç noktası olarak `hello-modules.nf` iş akışı betiğini kullanacağız.
Bu betik, bu eğitim kursunun 3. Bölümünde üretilen betiğe eşdeğerdir, ancak çıktı hedeflerini değiştirdik:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Her şeyin çalıştığından emin olmak için, herhangi bir değişiklik yapmadan önce betiği bir kez çalıştırın:

```bash
nextflow run hello-modules.nf
```

??? success "Komut çıktısı"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Daha önce olduğu gibi, çıktı dosyalarını `output` bloğunda belirtilen dizinde (burada `results/hello_modules/`) bulacaksınız.

??? abstract "Dizin içeriği"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Bu sizin için çalıştıysa, iş akışı kodunuzu nasıl modülerleştireceğinizi öğrenmeye hazırsınız.

---

## 1. Modülleri depolamak için bir dizin oluşturun

Modüllerinizi belirli bir dizinde saklamak en iyi uygulamadır.
Bu dizine istediğiniz adı verebilirsiniz, ancak geleneksel olarak `modules/` olarak adlandırılır.

```bash
mkdir modules
```

---

## 2. `sayHello()` için bir modül oluşturun

En basit haliyle, mevcut bir süreci modüle dönüştürmek kopyala-yapıştır işleminden biraz daha fazlasıdır.
Modül için bir dosya taslağı oluşturacağız, ilgili kodu kopyalayacağız ve ardından ana iş akışı dosyasından sileceğiz.

Ardından tek yapmamız gereken, Nextflow'un çalışma zamanında ilgili kodu çekmesini sağlamak için bir `include` ifadesi eklemek olacaktır.

### 2.1. Yeni modül için bir dosya taslağı oluşturun

`sayHello.nf` adında modül için boş bir dosya oluşturalım.

```bash
touch modules/sayHello.nf
```

Bu bize süreç kodunu koyacak bir yer verir.

### 2.2. `sayHello` süreç kodunu modül dosyasına taşıyın

Tüm süreç tanımını iş akışı dosyasından modül dosyasına kopyalayın.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Bir dosyaya 'Hello World!' yazdırmak için echo kullan
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

Bu işlem tamamlandıktan sonra, süreç tanımını iş akışı dosyasından silin.

### 2.3. Workflow bloğundan önce bir include bildirimi ekleyin

Bir modülden süreç dahil etme sözdizimi oldukça basittir:

```groovy title="Sözdizimi: include bildirimi"
include { <PROCESS_NAME> } from '<modül_yolu>'
```

Bunu `params` bloğunun üstüne ekleyelim ve uygun şekilde dolduralım.

=== "Sonra"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Önce"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Pipeline parametreleri
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Süreç adını `sayHello` ve modül kodunu içeren dosyanın yolunu `./modules/sayHello.nf` olarak doldurduğumuzu görüyorsunuz.

### 2.4. İş akışını çalıştırın

İş akışını daha önce olduğu gibi temelde aynı kod ve girdilerle çalıştırıyoruz, bu yüzden `-resume` bayrağıyla çalıştıralım ve ne olduğunu görelim.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Bu çok hızlı çalışmalıdır çünkü her şey önbellekte.
Yayınlanan çıktıları kontrol etmekten çekinmeyin.

Nextflow, kod birden fazla dosyaya bölünmüş olsa bile, hala yapılacak aynı işin olduğunu fark etti.

### Özet

Bir süreci yerel bir modüle nasıl çıkaracağınızı biliyorsunuz ve bunu yapmanın iş akışının devam ettirilebilirliğini bozmadığını biliyorsunuz.

### Sırada ne var?

Daha fazla modül oluşturma pratiği yapın.
Birini yaptıktan sonra, bir milyon tane daha yapabilirsiniz...
Ama şimdilik sadece iki tane daha yapalım.

---

## 3. `convertToUpper()` sürecini modülerleştirin

### 3.1. Yeni modül için bir dosya taslağı oluşturun

`convertToUpper.nf` adında modül için boş bir dosya oluşturun.

```bash
touch modules/convertToUpper.nf
```

### 3.2. `convertToUpper` süreç kodunu modül dosyasına taşıyın

Tüm süreç tanımını iş akışı dosyasından modül dosyasına kopyalayın.

```groovy title="modules/convertToUpper.nf" linenums="1"
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
```

Bu işlem tamamlandıktan sonra, süreç tanımını iş akışı dosyasından silin.

### 3.3. `params` bloğundan önce bir include bildirimi ekleyin

Include bildirimini `params` bloğunun üstüne ekleyin ve uygun şekilde doldurun.

=== "Sonra"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Önce"

    ```groovy title="hello-modules.nf" linenums="23"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Bu size çok tanıdık gelmeye başlamalı.

### 3.4. İş akışını tekrar çalıştırın

Bunu `-resume` bayrağıyla çalıştırın.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretmelidir.

İki tane bitti, bir tane daha kaldı!

---

## 4. `collectGreetings()` sürecini modülerleştirin

### 4.1. Yeni modül için bir dosya taslağı oluşturun

`collectGreetings.nf` adında modül için boş bir dosya oluşturun.

```bash
touch modules/collectGreetings.nf
```

### 4.2. `collectGreetings` süreç kodunu modül dosyasına taşıyın

Tüm süreç tanımını iş akışı dosyasından modül dosyasına kopyalayın.

```groovy title="modules/collectGreetings.nf" linenums="1"
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
```

Bu işlem tamamlandıktan sonra, süreç tanımını iş akışı dosyasından silin.

### 4.3. `params` bloğundan önce bir include bildirimi ekleyin

Include bildirimini `params` bloğunun üstüne ekleyin ve uygun şekilde doldurun.

=== "Sonra"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Önce"

    ```groovy title="hello-modules.nf" linenums="3"
    // Modülleri dahil et
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Pipeline parametreleri
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Son bir tane!

### 4.4. İş akışını çalıştırın

Bunu `-resume` bayrağıyla çalıştırın.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Bu hala daha önce olduğu gibi aynı çıktıyı üretmelidir.

### Özet

Bir iş akışında birden fazla süreci nasıl modülerleştireceğinizi biliyorsunuz.

Tebrikler, tüm bu işi yaptınız ve boru hattının çalışma şeklinde kesinlikle hiçbir şey değişmedi!

Şaka bir yana, artık kodunuz daha modüler ve bu süreçlerden birini çağıran başka bir boru hattı yazmaya karar verirseniz, ilgili modülü kullanmak için sadece kısa bir `include` ifadesi yazmanız yeterli.
Bu, kodu kopyala-yapıştır yapmaktan daha iyidir, çünkü daha sonra modülü geliştirmeye karar verirseniz, tüm boru hatlarınız iyileştirmeleri miras alacaktır.

### Sırada ne var?

İsterseniz kısa bir mola verin.

Hazır olduğunuzda, yazılım bağımlılıklarını daha rahat ve tekrarlanabilir bir şekilde yönetmek için konteynırların nasıl kullanılacağını öğrenmek için [**Bölüm 5: Merhaba Konteynerler**](./05_hello_containers.md) bölümüne geçin.

---

## Test

<quiz>
Nextflow'da modül nedir?
- [ ] Bir yapılandırma dosyası
- [x] Süreç tanımları içerebilen bağımsız bir dosya
- [ ] Bir iş akışı tanımı
- [ ] Bir kanal operatörü

Daha fazla bilgi: [2. `sayHello()` için bir modül oluşturun](#2-sayhello-için-bir-modül-oluşturun)
</quiz>

<quiz>
Modül dosyalarını depolamak için genellikle hangi kural kullanılır?
- [ ] İş akışıyla aynı dizinde
- [ ] Bir `bin/` dizininde
- [x] Bir `modules/` dizininde
- [ ] Bir `lib/` dizininde

Daha fazla bilgi: [1. Modülleri depolamak için bir dizin oluşturun](#1-modülleri-depolamak-için-bir-dizin-oluşturun)
</quiz>

<quiz>
Bir modül kullanmak için doğru sözdizimi nedir?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Daha fazla bilgi: [2.3. Bir include bildirimi ekleyin](#23-workflow-bloğundan-önce-bir-include-bildirimi-ekleyin)
</quiz>

<quiz>
Modüller kullanıldığında `-resume` işlevselliğine ne olur?
- [ ] Artık çalışmaz
- [ ] Ek yapılandırma gerektirir
- [x] Daha önce olduğu gibi çalışır
- [ ] Sadece yerel modüller için çalışır
</quiz>

<quiz>
Modül kullanmanın faydaları nelerdir? (Tümünü seçin)
- [x] İş akışları arasında kod yeniden kullanılabilirliği
- [x] Daha kolay bakım
- [x] İş akışı kodunun daha iyi organizasyonu
- [ ] Daha hızlı yürütme hızı
</quiz>
