# İç İçe Workflow'lar

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay Zeka Destekli Çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bir pipeline geliştirirken, kendinizi farklı veri türleri veya analiz adımları için benzer process dizileri oluştururken bulursunuz. Bu process dizilerini kopyalayıp yapıştırarak, bakımı zor çoğaltılmış kod oluşturabilirsiniz; veya anlaşılması ve değiştirilmesi zor tek bir devasa iş akışı oluşturabilirsiniz.

Nextflow'un en güçlü özelliklerinden biri, karmaşık pipeline'ları daha küçük, yeniden kullanılabilir workflow modüllerinden oluşturma yeteneğidir. Bu modüler yaklaşım, pipeline'ların geliştirilmesini, test edilmesini ve bakımını kolaylaştırır.

### Öğrenme hedefleri

Bu yan görevde, ayrı ayrı test edilebilen ve kullanılabilen workflow modülleri geliştirmeyi, bu modülleri daha büyük bir pipeline'da birleştirmeyi ve modüller arasında veri akışını yönetmeyi keşfedeceğiz.

Bu yan görevin sonunda şunları yapabileceksiniz:

- Karmaşık pipeline'ları mantıksal, yeniden kullanılabilir birimlere ayırma
- Her workflow modülünü bağımsız olarak test etme
- Yeni pipeline'lar oluşturmak için workflow'ları birleştirme
- Ortak workflow modüllerini farklı pipeline'lar arasında paylaşma
- Kodunuzu daha bakımı kolay ve anlaşılır hale getirme

Bu beceriler, temiz ve bakımı kolay kod yapısını korurken karmaşık pipeline'lar oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan göreve başlamadan önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramları ve mekanizmalarını (process'ler, channel'lar, operatörler, modüller) rahatça kullanabiliyor olmalısınız

---

## 0. Başlangıç

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md)'nda açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Proje dizinine geçin

Bu eğitim için dosyaların bulunduğu dizine geçelim.

```bash
cd side-quests/workflows_of_workflows
```

VSCode'u bu dizine odaklanacak şekilde ayarlayabilirsiniz:

```bash
code .
```

#### Materyalleri inceleyin

'Hello Nextflow'da öğrendiklerinizin üzerine inşa eden birkaç process tanımı içeren bir `modules` dizini bulacaksınız:

```console title="Dizin içeriği"
modules/
├── say_hello.nf             # Bir selamlama oluşturur (Hello Nextflow'dan)
├── say_hello_upper.nf       # Büyük harfe dönüştürür (Hello Nextflow'dan)
├── timestamp_greeting.nf    # Selamlamalara zaman damgası ekler
├── validate_name.nf         # Girdi isimlerini doğrular
└── reverse_text.nf          # Metin içeriğini tersine çevirir
```

#### Görevi inceleyin

Göreviniz, bu modülleri iki ayrı workflow'da birleştirmek ve ardından bunları ana bir workflow'da compose etmektir:

- İsimleri doğrulayan, selamlamalar oluşturan ve zaman damgaları ekleyen bir `GREETING_WORKFLOW`
- Metni büyük harfe dönüştüren ve tersine çeviren bir `TRANSFORM_WORKFLOW`

<!-- TODO: Metadata yan görevinde yapıldığına benzer şekilde biraz daha ayrıntı verin -->

#### Hazırlık kontrol listesi

Başlamaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Greeting Workflow'u Oluşturma

İsimleri doğrulayan ve zaman damgalı selamlamalar oluşturan bir workflow oluşturarak başlayalım.

### 1.1. Workflow yapısını oluşturma

```bash title="Workflow dizini ve dosyası oluşturma"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. İlk (alt)workflow kodunu ekleme

Bu kodu `workflows/greeting.nf` dosyasına ekleyin:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Process'leri zincirle: doğrula -> selamlama oluştur -> zaman damgası ekle
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Bu, 'Hello Nextflow' eğitiminde gördüklerinize benzer bir yapıya sahip, bağımsız olarak test edebileceğimiz tam bir workflow. Şimdi deneyelim:

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

Bu beklendiği gibi çalışıyor, ancak composable yapmak için değiştirmemiz gereken birkaç şey var.

### 1.3. Workflow'u composable yapma

Composable workflow'ların 'Hello Nextflow' eğitiminde gördüklerinizden bazı farklılıkları vardır:

- Workflow bloğunun adlandırılması gerekir
- Girdiler `take:` anahtar kelimesi kullanılarak bildirilir
- Workflow içeriği `main:` bloğunun içine yerleştirilir
- Çıktılar `emit:` anahtar kelimesi kullanılarak bildirilir

Greeting workflow'u bu yapıya uygun hale getirelim. Kodu aşağıdaki gibi değiştirin:

<!-- TODO: önce/sonra sekmelerine geçin -->

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 8 14 15 16"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // İsimler içeren girdi channel'ı

    main:
        // Process'leri zincirle: doğrula -> selamlama oluştur -> zaman damgası ekle
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Orijinal selamlamalar
        timestamped = timestamped_ch  // Zaman damgalı selamlamalar
}
```

Workflow'un artık adlandırıldığını ve `take:` ile `emit:` bloğuna sahip olduğunu görebilirsiniz; bunlar daha üst düzey bir workflow oluşturmak için kullanacağımız bağlantılardır.
Workflow içeriği de `main:` bloğunun içine yerleştirilmiştir. Ayrıca `names_ch` girdi channel bildirimini kaldırdığımızı not edin, çünkü artık workflow'a argüman olarak geçiriliyor.

Workflow'un beklendiği gibi çalışıp çalışmadığını görmek için tekrar test edelim:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Komut çıktısı"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Bu size başka bir yeni kavram hakkında bilgi veriyor, bir 'giriş workflow'u'. Giriş workflow'u, bir Nextflow script'i çalıştırdığınızda çağrılan workflow'dur. Varsayılan olarak, Nextflow mevcut olduğunda adlandırılmamış bir workflow'u giriş workflow'u olarak kullanır ve şimdiye kadar yaptığınız şey de buydu, şöyle başlayan workflow bloklarıyla:

```groovy title="hello.nf" linenums="1"
workflow {
```

Ancak greeting workflow'umuzda adlandırılmamış bir workflow yok, bunun yerine adlandırılmış bir workflow var:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Bu yüzden Nextflow hata verdi ve istediğimizi yapmadı.

`take:`/`emit:` sözdizimini workflow'u doğrudan çağırabilmek için eklemedik - diğer workflow'larla compose edebilmek için ekledik. Çözüm, adlandırılmış workflow'umuzu içe aktaran ve çağıran adlandırılmamış bir giriş workflow'una sahip bir ana script oluşturmaktır.

### 1.4. Ana workflow'u oluşturma ve test etme

Şimdi `greeting` workflow'unu içe aktaran ve kullanan bir ana workflow oluşturacağız.

`main.nf` oluşturun:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Orijinal: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Zaman damgalı: $it" }
}

```

Bu dosyadaki workflow girişimizin adlandırılmamış olduğunu not edin, çünkü onu bir giriş workflow'u olarak kullanacağız.

Bunu çalıştırın ve çıktıyı görün:

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
    Orijinal: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Orijinal: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Orijinal: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Zaman damgalı: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Zaman damgalı: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Zaman damgalı: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Çalışıyor! Adlandırılmış greeting workflow'unu, adlandırılmamış bir giriş `workflow` bloğuna sahip bir ana workflow'a sardık. Ana workflow, `GREETING_WORKFLOW` workflow'unu neredeyse (tam olarak değil) bir process gibi kullanıyor ve `names` channel'ını argüman olarak geçiriyor.

### Özet

Bu bölümde birkaç önemli kavram öğrendiniz:

- **Adlandırılmış Workflow'lar**: İçe aktarılabilen ve yeniden kullanılabilen adlandırılmış bir workflow (`GREETING_WORKFLOW`) oluşturma
- **Workflow Arayüzleri**: Composable bir workflow oluşturmak için `take:` ile net girdiler ve `emit:` ile çıktılar tanımlama
- **Giriş Noktaları**: Nextflow'un bir script'i çalıştırmak için adlandırılmamış bir giriş workflow'una ihtiyaç duyduğunu anlama
- **Workflow Kompozisyonu**: Başka bir workflow içinde adlandırılmış bir workflow'u içe aktarma ve kullanma
- **Workflow Ad Alanları**: `.out` ad alanını kullanarak workflow çıktılarına erişme (`GREETING_WORKFLOW.out.greetings`)

Artık şunları yapan çalışan bir greeting workflow'unuz var:

- Girdi olarak bir isimler channel'ı alır
- Her ismi doğrular
- Her geçerli isim için bir selamlama oluşturur
- Selamlamalara zaman damgaları ekler
- Hem orijinal hem de zaman damgalı selamlamaları çıktı olarak sunar

Bu modüler yaklaşım, greeting workflow'unu bağımsız olarak test etmenize veya daha büyük pipeline'larda bir bileşen olarak kullanmanıza olanak tanır.

---

## 2. Transform Workflow'unu Ekleme

Şimdi selamlamalara metin dönüşümleri uygulayan bir workflow oluşturalım.

### 2.1. Workflow dosyasını oluşturma

```bash
touch workflows/transform.nf
```

### 2.2. Workflow kodunu ekleme

Bu kodu `workflows/transform.nf` dosyasına ekleyin:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Mesajlar içeren girdi channel'ı

    main:
        // Dönüşümleri sırayla uygula
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Büyük harfli selamlamalar
        reversed = reversed_ch  // Tersine çevrilmiş büyük harfli selamlamalar
}
```

Composable sözdiziminin açıklamasını burada tekrarlamayacağız, ancak adlandırılmış workflow'un yine `take:` ve `emit:` bloğuyla bildirildiğini ve workflow içeriğinin `main:` bloğunun içine yerleştirildiğini not edin.

### 2.3. Ana workflow'u güncelleme

Her iki workflow'u da kullanmak için `main.nf`'i güncelleyin:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Greeting workflow'unu çalıştır
    GREETING_WORKFLOW(names)

    // Transform workflow'unu çalıştır
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Sonuçları görüntüle
    TRANSFORM_WORKFLOW.out.upper.view { "Büyük harf: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Tersine çevrilmiş: $it" }
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
    Büyük harf: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Büyük harf: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Büyük harf: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Tersine çevrilmiş: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Tersine çevrilmiş: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Tersine çevrilmiş: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Tersine çevrilmiş dosyalardan birine bakarsanız, bunun selamlamanın büyük harfli versiyonunun tersine çevrilmiş hali olduğunu göreceksiniz:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Tersine çevrilmiş dosya içeriği"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Özet

Artık şunları yapan tam bir pipeline'ınız olmalı:

- İsimleri greeting workflow'u aracılığıyla işler
- Zaman damgalı selamlamaları transform workflow'una besler
- Selamlamaların hem büyük harfli hem de tersine çevrilmiş versiyonlarını üretir

---

## Özet

Bu yan görevde, daha küçük, yeniden kullanılabilir bileşenlerden karmaşık pipeline'lar oluşturmamıza olanak tanıyan Nextflow'daki güçlü workflow kompozisyonu kavramını keşfettik.

Bu modüler yaklaşım, monolitik pipeline'lara göre çeşitli avantajlar sunar:

- Her workflow bağımsız olarak geliştirilebilir, test edilebilir ve hata ayıklanabilir
- Workflow'lar farklı pipeline'lar arasında yeniden kullanılabilir
- Genel pipeline yapısı daha okunabilir ve bakımı kolay hale gelir
- Arayüzler tutarlı kaldığı sürece bir workflow'daki değişiklikler diğerlerini etkilemez
- Giriş noktaları, pipeline'ınızın farklı bölümlerini gerektiği gibi çalıştırmak için yapılandırılabilir

_Ancak, workflow'ları çağırmanın process'leri çağırmak gibi olduğunu ancak aslında aynı şey olmadığını not etmek önemlidir. Örneğin, N boyutunda bir channel ile çağırarak bir workflow'u N kez çalıştıramazsınız - workflow'a N boyutunda bir channel geçirmeniz ve dahili olarak yinelemeniz gerekir._

Bu teknikleri kendi çalışmanıza uygulamak, bakımı ve ölçeklenebilir kalırken karmaşık biyoinformatik görevleri işleyebilen daha sofistike Nextflow pipeline'ları oluşturmanıza olanak tanıyacaktır.

### Temel desenler

1.  **Workflow yapısı**: `take:` ve `emit:` sözdizimini kullanarak her workflow için net girdiler ve çıktılar tanımladık, bileşenler arasında iyi tanımlanmış arayüzler oluşturduk ve workflow mantığını `main:` bloğu içine sardık.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Girdi channel'ları burada bildirilir
            input_ch

        main:
            // Workflow mantığı buraya gelir
            // Process'ler çağrılır ve channel'lar manipüle edilir
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Çıktı channel'ları burada bildirilir
            output_ch = result_ch
    }
    ```

2.  **Workflow içe aktarmaları:** İki bağımsız workflow modülü oluşturduk ve bunları include ifadeleriyle ana pipeline'a aktardık.

    - Tek bir workflow'u içe aktarma

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Birden fazla workflow'u içe aktarma

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Ad çakışmalarını önlemek için takma adla içe aktarma

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Giriş noktaları**: Nextflow, yürütmeye nereden başlayacağını bilmek için adlandırılmamış bir giriş workflow'u gerektirir. Bu giriş workflow'u, adlandırılmış workflow'larınızı çağırır.

    - Adlandırılmamış workflow (giriş noktası)

    ```groovy
    workflow {
        // Script çalıştırıldığında burası giriş noktasıdır
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Adlandırılmış workflow (giriş workflow'undan çağrılır)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Giriş workflow'undan çağrılmalıdır
    }
    ```

4.  **Veri akışını yönetme:** Ad alanı gösterimini (`WORKFLOW_NAME.out.channel_name`) kullanarak workflow çıktılarına erişmeyi ve bunları diğer workflow'lara geçirmeyi öğrendik.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Ek kaynaklar

- [Nextflow Workflow Belgeleri](https://www.nextflow.io/docs/latest/workflow.html)
- [Channel Operatörleri Referansı](https://www.nextflow.io/docs/latest/operator.html)
- [Error Strategy Belgeleri](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Sırada ne var?

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
