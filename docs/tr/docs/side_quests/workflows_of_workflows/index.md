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

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmak.
- Temel Nextflow kavramları ve mekanizmalarını (süreçler, kanallar, operatörler, modüller) rahatça kullanabilmek.

---

## 0. Başlarken

#### Eğitim kod alanını açın

Henüz yapmadıysanız, eğitim ortamını [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı şekilde açtığınızdan emin olun.

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

#### Materyalleri inceleyin

'Hello Nextflow' kursunda öğrendiklerinizi temel alan çeşitli süreç tanımlarını içeren bir `modules` dizini bulacaksınız:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

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

## 1. Greeting Workflow'u Oluşturun

İsimleri doğrulayan ve zaman damgalı selamlamalar üreten bir iş akışı oluşturarak başlayalım.

### 1.1. İş akışı yapısını oluşturun

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. İlk (alt) iş akışı kodunu ekleyin

Bu kodu `workflows/greeting.nf` dosyasına ekleyin:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Süreçleri zincirleyin: doğrula -> selamlama oluştur -> zaman damgası ekle
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Bu, 'Hello Nextflow' eğitiminde gördüklerinize benzer bir yapıya sahip, bağımsız olarak test edebileceğimiz eksiksiz bir iş akışıdır. Şimdi deneyelim:

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

Bu beklendiği gibi çalışıyor; ancak birleştirilebilir hale getirmek için birkaç değişiklik yapmamız gerekiyor.

### 1.3. İş akışını birleştirilebilir hale getirin

Birleştirilebilir iş akışları, 'Hello Nextflow' eğitiminde gördüklerinizden bazı farklılıklar taşır:

- İş akışı bloğunun adlandırılması gerekir
- Girdiler `take:` anahtar sözcüğü kullanılarak tanımlanır
- İş akışı içeriği `main:` bloğunun içine yerleştirilir
- Çıktılar `emit:` anahtar sözcüğü kullanılarak tanımlanır

Greeting iş akışını bu yapıya uyacak şekilde güncelleyelim. Kodu aşağıdaki şekilde değiştirin:

<!-- TODO: switch to before/after tabs -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // İsimlerle dolu girdi kanalı

    main:
        // Süreçleri zincirleyin: doğrula -> selamlama oluştur -> zaman damgası ekle
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Orijinal selamlamalar
        timestamped = timestamped_ch  // Zaman damgalı selamlamalar
}
```

İş akışının artık adlandırıldığını ve bir `take:` ile `emit:` bloğuna sahip olduğunu görebilirsiniz; bunlar, daha üst düzey bir iş akışı oluşturmak için kullanacağımız bağlantı noktalarıdır.
İş akışı içeriği de `main:` bloğunun içine yerleştirilmiştir. Ayrıca `names_ch` girdi kanalı tanımını kaldırdığımıza dikkat edin; çünkü bu kanal artık iş akışına argüman olarak iletilmektedir.

İş akışının beklendiği gibi çalışıp çalışmadığını görmek için tekrar test edelim:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Bu, 'giriş iş akışı' adı verilen yeni bir kavramı ortaya koymaktadır. Giriş iş akışı, bir Nextflow betiği çalıştırıldığında çağrılan iş akışıdır. Nextflow, varsayılan olarak mevcut olduğunda adsız bir iş akışını giriş iş akışı olarak kullanır; şimdiye kadar yaptığınız da buydu; iş akışı blokları şu şekilde başlıyordu:

```groovy title="hello.nf" linenums="1"
workflow {
```

Ancak greeting iş akışımızda adsız bir iş akışı yoktur; bunun yerine adlandırılmış bir iş akışımız vardır:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Nextflow'un hata verip istediğimizi yapmamasının nedeni budur.

`take:`/`emit:` sözdizimini iş akışını doğrudan çağırabilmek için eklemedik; bunu, diğer iş akışlarıyla birleştirebilmek için ekledik. Çözüm, adlandırılmış iş akışımızı içe aktarıp çağıran, adsız bir giriş iş akışına sahip bir ana betik oluşturmaktır.

### 1.4. Ana iş akışını oluşturun ve test edin

Şimdi `greeting` iş akışını içe aktarıp kullanan bir ana iş akışı oluşturacağız.

`main.nf` dosyasını oluşturun:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Bu dosyadaki iş akışı girişinin adsız olduğuna dikkat edin; bunun nedeni, bunu bir giriş iş akışı olarak kullanacak olmamızdır.

Çalıştırın ve çıktıyı inceleyin:

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
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Çalışıyor! Adlandırılmış greeting iş akışını, adsız bir giriş `workflow` bloğuna sahip bir ana iş akışıyla sardık. Ana iş akışı, `GREETING_WORKFLOW` iş akışını neredeyse (tam olarak değil) bir süreç gibi kullanıyor ve `names` kanalını argüman olarak iletiyor.

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

## 2. Transform Workflow'u Ekleyin

Şimdi selamlamalara metin dönüşümleri uygulayan bir iş akışı oluşturalım.

### 2.1. İş akışı dosyasını oluşturun

```bash
touch workflows/transform.nf
```

### 2.2. İş akışı kodunu ekleyin

Bu kodu `workflows/transform.nf` dosyasına ekleyin:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Mesajlarla dolu girdi kanalı

    main:
        // Dönüşümleri sırayla uygulayın
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Büyük harfli selamlamalar
        reversed = reversed_ch  // Tersine çevrilmiş büyük harfli selamlamalar
}
```

Birleştirilebilir sözdiziminin açıklamasını burada tekrarlamayacağız; ancak adlandırılmış iş akışının yine bir `take:` ve `emit:` bloğuyla tanımlandığına ve iş akışı içeriğinin `main:` bloğunun içine yerleştirildiğine dikkat edin.

### 2.3. Ana iş akışını güncelleyin

Her iki iş akışını da kullanmak için `main.nf` dosyasını güncelleyin:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Greeting iş akışını çalıştırın
    GREETING_WORKFLOW(names)

    // Transform iş akışını çalıştırın
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Sonuçları görüntüleyin
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Tam pipeline'ı çalıştırın:

```bash
nextflow run main.nf
```

??? success "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Tersine çevrilmiş dosyalardan birine bakarsanız, selamlamanın büyük harfli sürümünün tersine çevrildiğini göreceksiniz:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

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

[Yan Görevler menüsüne](../) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
