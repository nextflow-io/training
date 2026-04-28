# İş Akışlarının İş Akışları

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bir pipeline geliştirirken, farklı veri türleri veya analiz adımları için benzer süreç dizileri oluşturduğunuzu sık sık fark edersiniz. Bu süreç dizilerini kopyalayıp yapıştırarak bakımı zor olan tekrarlı kodlara yol açabilirsiniz; ya da anlaşılması ve değiştirilmesi güç olan tek bir devasa iş akışı oluşturabilirsiniz.

Nextflow'un en güçlü özelliklerinden biri, karmaşık pipeline'ları daha küçük, yeniden kullanılabilir iş akışı modüllerinden oluşturabilmesidir. Bu modüler yaklaşım, pipeline'ların geliştirilmesini, test edilmesini ve bakımını kolaylaştırır.

### Öğrenme hedefleri

Bu yan görevde, bağımsız olarak test edilip kullanılabilen iş akışı modüllerinin nasıl geliştirileceğini, bu modüllerin daha büyük bir pipeline'da nasıl bir araya getirileceğini ve modüller arasındaki veri akışının nasıl yönetileceğini inceleyeceğiz.

Bu yan görevin sonunda şunları yapabileceksiniz:

- Karmaşık pipeline'ları mantıksal, yeniden kullanılabilir birimlere ayırmak
- Her iş akışı modülünü bağımsız olarak test etmek
- Yeni pipeline'lar oluşturmak için iş akışlarını bir araya getirmek
- Ortak iş akışı modüllerini farklı pipeline'lar arasında paylaşmak
- Kodunuzu daha bakımı kolay ve anlaşılır hale getirmek

Bu beceriler, temiz ve bakımı kolay bir kod yapısını korurken karmaşık pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce şunları yapmış olmalısınız:

- [Hello Nextflow](../../hello_nextflow/index.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, modüller) rahatça kullanabilmek.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için gerekli dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/workflows_of_workflows
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

Düzenleyici, proje dizinine odaklanmış şekilde açılır.

#### Materyalleri inceleyin

Süreç tanımlarını içeren bir `modules` dizini, önceden yazılmış iki iş akışı betiğini içeren bir `workflows` dizini ve aşamalı olarak güncelleyeceğiniz bir `main.nf` dosyası bulacaksınız:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

`modules/` dizini tek tek süreç tanımlarını, `workflows/` dizini ise bu yan görevde çalışacağınız önceden yazılmış iki iş akışı betiğini içerir.

#### Görevi inceleyin

Göreviniz, bu modülleri daha sonra bir ana iş akışında bir araya getireceğimiz iki ayrı iş akışında toplamaktır:

- İsimleri doğrulayan, selamlamalar oluşturan ve zaman damgası ekleyen bir `GREETING_WORKFLOW`
- Metni büyük harfe dönüştüren ve tersine çeviren bir `TRANSFORM_WORKFLOW`

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Hazırlık kontrol listesi

Başlamaya hazır mısınız?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Kod alanım çalışıyor
- [ ] Çalışma dizinimi uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebildiyseniz, başlayabilirsiniz.

---

## 1. Greeting workflow'unu pipeline'a ekleyin

Greeting iş akışı, isimleri doğrular ve zaman damgalı selamlamalar üretir.

### 1.1. Greeting workflow'unu inceleyin ve çalıştırın

`workflows/greeting.nf` dosyasını açın ve kodu inceleyin:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Süreçleri zincirleyin: doğrula -> selamlama oluştur -> zaman damgası ekle
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

Bu, 'Hello Nextflow' eğitiminde gördüklerinizle aynı yapıya sahip, eksiksiz ve bağımsız bir iş akışıdır. Girdi isimlerini sabit olarak kodlar, üç süreci zincirler ve iki çıktı yayımlar.

Her şeyin çalıştığını doğrulamak için çalıştırın:

```bash
nextflow run workflows/greeting.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Diğer iş akışlarıyla birleştirilebilir hale getirmek için birkaç değişiklik yapmamız gerekiyor.

### 1.2. İş akışını birleştirilebilir hale getirin

Bir iş akışını birleştirilebilir hale getirmek için dört şeyin değişmesi gerekir: iş akışına bir ad verilir, girdiler `take:` bloğuna taşınır, çıktılar `emit:` bloğuna taşınır ve bağımsız `publish:`/`output {}` blokları kaldırılır (bunlar giriş iş akışına aittir).

Bu değişiklikleri tek tek inceleyelim.

#### 1.2.1. İş akışını adlandırın

İş akışına bir ad verin; böylece üst iş akışından içe aktarılabilir hale gelir.

=== "Sonra"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Önce"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

Bir adla birlikte iş akışı, diğer betiklere aktarılabilir.

#### 1.2.2. `take:` ile girdileri tanımlayın

Sabit kodlanmış kanal tanımını, iş akışının beklediği girdileri bildiren bir `take:` bloğuyla değiştirin. `take:` bloğu `main:`'den önce gelir ve `names_ch = channel.of(...)` satırı kaldırılır.

=== "Sonra"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch // İsimlerle dolu girdi kanalı

        main:
        // Süreçleri zincirleyin: doğrula -> selamlama oluştur -> zaman damgası ekle
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Önce"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Süreçleri zincirleyin: doğrula -> selamlama oluştur -> zaman damgası ekle
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

`take:` bloğu kanalı yalnızca adıyla tanımlar; içine ne gireceğinin ayrıntıları üst iş akışı tarafından belirlenecektir.

#### 1.2.3. `emit:` ile çıktıları tanımlayın

`publish:` bölümünü ve `output {}` bloğunu kaldırarak bunların yerine çıktıları adlandıran bir `emit:` bloğu ekleyin.

=== "Sonra"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch // Orijinal selamlamalar
        timestamped = timestamped_ch // Zaman damgalı selamlamalar
    }
    ```

=== "Önce"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

`emit:` bloğu, üst iş akışlarının `GREETING_WORKFLOW.out.greetings` ve `GREETING_WORKFLOW.out.timestamped` aracılığıyla erişebileceği adlandırılmış çıktıları ortaya koyar.

#### 1.2.4. Sonucu doğrulayın ve test edin

Üç değişikliğin tamamından sonra dosyanın tamamı şöyle görünmelidir:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch // İsimlerle dolu girdi kanalı

    main:
    // Süreçleri zincirleyin: doğrula -> selamlama oluştur -> zaman damgası ekle
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch // Orijinal selamlamalar
    timestamped = timestamped_ch // Zaman damgalı selamlamalar
}
```

Şimdi doğrudan çalıştırmayı deneyin:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Bu, önemli bir kavramı ortaya koymaktadır: **giriş iş akışı**.
Nextflow, bir betiği doğrudan çalıştırdığınızda adsız bir `workflow {}` bloğunu giriş noktası olarak kullanır.
`GREETING_WORKFLOW` adlandırılmış olduğundan Nextflow, onu kendi başına nasıl çalıştıracağını bilmez.

Bu kasıtlıdır; birleştirilebilir iş akışları doğrudan çalıştırılmak için değil, bir giriş iş akışından çağrılmak üzere tasarlanmıştır. Çözüm, `GREETING_WORKFLOW`'u içe aktarıp çağıran `main.nf` dosyasında bir giriş iş akışı oluşturmaktır.

### 1.3. Ana iş akışını güncelleyin ve test edin

Şimdi greeting iş akışını çağırmak için ana iş akışını güncelleyelim.

#### 1.3.1. Greeting workflow'unu dahil edin ve çağırın

`include` ifadesini ekleyin, iş akışı gövdesini `GREETING_WORKFLOW`'u çağıracak şekilde güncelleyin ve `publish:` bölümündeki `channel.empty()` yer tutucusunu değiştirin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="1 7 8 11"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Greeting iş akışını çalıştırın
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }
    ```

Nextflow'un bunu pipeline giriş noktası olarak kullanabilmesi için giriş iş akışı adsız kalır.

#### 1.3.2. Output bloğunu güncelleyin

Yayımlanan selamlamaları bir `greetings/` alt dizinine yönlendirmek için `path` yönergesini ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="14" hl_lines="2 3"
    output {
        greetings {
        }
    }
    ```

#### 1.3.3. İş akışını çalıştırın

Çalıştığını test etmek için iş akışını çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    ```

??? abstract "Dizin içeriği"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "Dosya içeriği"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

Selamlama dosyaları `results/greetings/` dizinine yayımlanır. Ana iş akışı `GREETING_WORKFLOW`'u çağırır ve çıktısını doğrudan `publish:` bölümüne bağlar.

### Özetle

Bu bölümde birkaç önemli kavramı öğrendiniz:

- **Adlandırılmış İş Akışları**: İçe aktarılıp yeniden kullanılabilen adlandırılmış bir iş akışı (`GREETING_WORKFLOW`) oluşturma
- **İş Akışı Arayüzleri**: Birleştirilebilir bir iş akışı oluşturmak için `take:` ile açık girdiler ve `emit:` ile açık çıktılar tanımlama
- **Giriş Noktaları**: Nextflow'un bir betiği çalıştırmak için adsız bir giriş iş akışına ihtiyaç duyduğunu anlama
- **İş Akışı Birleştirme**: Adlandırılmış bir iş akışını başka bir iş akışı içinde içe aktarma ve kullanma
- **İş Akışı Ad Alanları**: `.out` ad alanını kullanarak iş akışı çıktılarına erişme (`GREETING_WORKFLOW.out.greetings`)

Artık çalışan bir greeting iş akışınız var. Bu iş akışı:

- Girdi olarak bir isimler kanalı alır
- Her ismi doğrular
- Her geçerli isim için bir selamlama oluşturur
- Selamlamalara zaman damgası ekler
- Hem orijinal hem de zaman damgalı selamlamaları çıktı olarak sunar

Bu modüler yaklaşım, greeting iş akışını bağımsız olarak test etmenize veya daha büyük pipeline'larda bir bileşen olarak kullanmanıza olanak tanır.

---

## 2. Dönüşüm workflow'unu pipeline'a ekleyin

Transform iş akışı, zaman damgalı selamlamalara metin dönüşümleri uygular.

### 2.1. İş akışını inceleyin ve çalıştırın

`workflows/transform.nf` dosyasını açın ve kodu inceleyin:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped_*.txt')

    // Dönüşümleri sırayla uygulayın
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

Bu bağımsız iş akışı, `greeting.nf` tarafından üretilen `results/` dizinindeki zaman damgalı selamlama dosyalarını okur, bunları büyük harfe dönüştürür ve ardından metni tersine çevirir.

1.1. bölümündeki selamlama sonuçlarıyla çalıştığını doğrulamak için çalıştırın:

```bash
nextflow run workflows/transform.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/transform.nf` [blissful_curie] DSL2 - revision: 4e7b1c9f02
    executor >  local (6)
    [3e/a14c29] process > SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [c8/51b9e3] process > REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

`GREETING_WORKFLOW` ile birleştirilebilir hale getirmek için 1.2. bölümündeki üç değişikliğin aynısı geçerlidir.

### 2.2. Birleştirilebilir hale getirin

1.2. bölümündeki üç değişikliğin aynısını uygulayın: iş akışını adlandırın, sabit kodlanmış girdiyi `take:` ile değiştirin ve `publish:`/`output {}` bloklarını `emit:` ile değiştirin.

Tamamlanmış dosya şöyle görünmelidir:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch // Mesajlarla dolu girdi kanalı

    main:
    // Dönüşümleri sırayla uygulayın
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch // Büyük harfli selamlamalar
    reversed = reversed_ch // Tersine çevrilmiş büyük harfli selamlamalar
}
```

Transform iş akışı artık birleştirilebilir durumdadır ve ana iş akışına aktarılmaya hazırdır.

### 2.3. Ana iş akışını güncelleyin ve test edin

Şimdi dönüşüm iş akışını çağırmak için ana iş akışını güncelleyelim.

#### 2.3.1. Dönüşüm workflow'unu dahil edin ve çağırın

Include ifadesini, zaman damgalı selamlamalar üzerinde zincirlenen `TRANSFORM_WORKFLOW` çağrısını ve iki yeni `publish:` girişini ekleyin:

=== "Sonra"

    ```groovy title="main.nf" linenums="1" hl_lines="2 11 12 16 17"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Greeting iş akışını çalıştırın
        GREETING_WORKFLOW(names)

        // Transform iş akışını çalıştırın
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Greeting iş akışını çalıştırın
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }
    ```

Bu, dönüşüm iş akışını zaman damgalı selamlamalar üzerinde çalıştıracaktır.

#### 2.3.2. Output bloğunu güncelleyin

`output {}` bloğuna `upper` ve `reversed` girişlerini ekleyin; her biri kendi alt dizini için bir `path` yönergesiyle birlikte:

=== "Sonra"

    ```groovy title="main.nf" linenums="20" hl_lines="5 6 7 8 9 10"
    output {
        greetings {
            path 'greetings'
        }
        upper {
            path 'upper'
        }
        reversed {
            path 'reversed'
        }
    }
    ```

=== "Önce"

    ```groovy title="main.nf" linenums="20" hl_lines="2 3 4 5"
    output {
        greetings {
            path 'greetings'
        }
    }
    ```

Bu, nihai çıktıları uygun dizinlere yayımlayacaktır.

#### 2.3.3. Tam pipeline'ı çalıştırın

Her şeyin çalıştığını test etmek için pipeline'ı çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Dizin içeriği"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "Dosya içeriği"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

Pipeline uçtan uca çalışıyor: selamlama büyük harfe dönüştürülmüş ve tersine çevrilmiştir.

### Özetle

Artık şunları yapan eksiksiz bir pipeline'ınız olmalı:

- İsimleri greeting iş akışı aracılığıyla işler
- Zaman damgalı selamlamaları transform iş akışına iletir
- Selamlamaların hem büyük harfli hem de tersine çevrilmiş sürümlerini üretir

---

## Özet

Bu yan görevde, Nextflow'da iş akışı birleştirmenin güçlü kavramını inceledik; bu kavram, karmaşık pipeline'ları daha küçük, yeniden kullanılabilir bileşenlerden oluşturmamıza olanak tanır.

Bu modüler yaklaşım, tek parça pipeline'lara kıyasla çeşitli avantajlar sunar:

- Her iş akışı bağımsız olarak geliştirilebilir, test edilebilir ve hata ayıklanabilir
- İş akışları farklı pipeline'larda yeniden kullanılabilir
- Genel pipeline yapısı daha okunabilir ve bakımı kolay hale gelir
- Arayüzler tutarlı kaldığı sürece bir iş akışındaki değişiklikler diğerlerini etkilemez
- Giriş noktaları, pipeline'ınızın farklı bölümlerini çalıştıracak şekilde yapılandırılabilir

_Ancak şunu belirtmek önemlidir: İş akışlarını çağırmak süreçleri çağırmaya biraz benzese de aslında aynı şey değildir. Örneğin, N boyutunda bir kanalla çağırarak bir iş akışını N kez çalıştıramazsınız; N boyutunda bir kanalı iş akışına iletmeniz ve dahili olarak yinelemeniz gerekir._

Bu teknikleri kendi çalışmalarınızda uygulamak, bakımı kolay ve ölçeklenebilir kalırken karmaşık biyoinformatik görevleri yerine getirebilen daha gelişmiş Nextflow pipeline'ları oluşturmanızı sağlayacaktır.

### Temel kalıplar

1.  **İş akışı yapısı**: Her iş akışı için `take:` ve `emit:` sözdizimini kullanarak açık girdiler ve çıktılar tanımladık; bileşenler arasında iyi tanımlanmış arayüzler oluşturduk ve iş akışı mantığını `main:` bloğunun içine yerleştirdik.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Girdi kanalları burada tanımlanır
            input_ch

        main:
            // İş akışı mantığı buraya gelir
            // Süreçlerin çağrıldığı ve kanalların işlendiği yer burasıdır
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Çıktı kanalları burada tanımlanır
            output_ch = result_ch
    }
    ```

2.  **İş akışı içe aktarmaları:** İki bağımsız iş akışı modülü oluşturduk ve bunları include ifadeleriyle bir ana pipeline'a aktardık.

    - Tek bir iş akışını içe aktarma

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Birden fazla iş akışını içe aktarma

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Ad çakışmalarını önlemek için takma adla içe aktarma

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Giriş noktaları**: Nextflow, yürütmenin nereden başlayacağını bilmek için adsız bir giriş iş akışı gerektirir. Bu giriş iş akışı, adlandırılmış iş akışlarınızı çağırır.

    - Adsız iş akışı (giriş noktası)

    ```groovy
    workflow {
        // Betik çalıştırıldığında burası giriş noktasıdır
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Adlandırılmış iş akışı (giriş iş akışından çağrılır)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Giriş iş akışından çağrılmalıdır
    }
    ```

4.  **Veri akışını yönetme:** Ad alanı gösterimini (`WORKFLOW_NAME.out.channel_name`) kullanarak iş akışı çıktılarına nasıl erişileceğini ve bunların diğer iş akışlarına nasıl iletileceğini öğrendik.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Ek kaynaklar

- [Nextflow İş Akışı Belgeleri](https://www.nextflow.io/docs/latest/workflow.html)
- [Kanal Operatörleri Referansı](https://www.nextflow.io/docs/latest/operator.html)
- [Hata Stratejisi Belgeleri](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Sırada ne var?

[Yan Görevler menüsüne](../index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
