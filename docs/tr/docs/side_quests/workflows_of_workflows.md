# İş Akışlarının İş Akışları

Bir boru hattı geliştirirken, genellikle farklı veri türleri veya analiz adımları için benzer süreç dizileri oluşturduğunuzu fark edersiniz. Bu süreç dizilerini kopyalayıp yapıştırarak, bakımı zor olan yinelenen kodlara yol açabilirsiniz; veya anlaşılması ve değiştirilmesi zor olan devasa bir iş akışı oluşturabilirsiniz.

Nextflow'un en güçlü özelliklerinden biri, karmaşık boru hatlarını daha küçük, yeniden kullanılabilir iş akışı modüllerinden oluşturabilme yeteneğidir. Bu modüler yaklaşım, boru hatlarının geliştirilmesini, test edilmesini ve bakımını kolaylaştırır.

### Öğrenme hedefleri

Bu yan görevde, ayrı olarak test edilebilen ve kullanılabilen iş akışı modüllerinin nasıl geliştirileceğini, bu modüllerin daha büyük bir boru hattında nasıl birleştirileceğini ve modüller arasında veri akışının nasıl yönetileceğini keşfedeceğiz.

Bu yan görevin sonunda şunları yapabileceksiniz:

- Karmaşık boru hatlarını mantıksal, yeniden kullanılabilir birimlere ayırmak
- Her iş akışı modülünü bağımsız olarak test etmek
- Yeni boru hatları oluşturmak için iş akışlarını karıştırıp eşleştirmek
- Ortak iş akışı modüllerini farklı boru hatları arasında paylaşmak
- Kodunuzu daha sürdürülebilir ve anlaşılır hale getirmek

Bu beceriler, temiz ve sürdürülebilir kod yapısını korurken karmaşık boru hatları oluşturmanıza yardımcı olacaktır.

### Ön koşullar

Bu yan görevi üstlenmeden önce:

- [Hello Nextflow](../hello_nextflow/README.md) eğitimini veya eşdeğer bir başlangıç kursunu tamamlamış olmalısınız.
- Temel Nextflow kavramlarını ve mekanizmalarını (süreçler, kanallar, operatörler, modüller) rahatça kullanabiliyor olmalısınız

---

## 0. Başlayın

#### Eğitim codespace'ini açın

Henüz yapmadıysanız, [Ortam Kurulumu](../envsetup/index.md) bölümünde açıklandığı gibi eğitim ortamını açtığınızdan emin olun.

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

#### Materyalleri gözden geçirin

'Hello Nextflow'da öğrendiklerinizin üzerine inşa edilen birkaç süreç tanımı içeren bir `modules` dizini bulacaksınız:

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

#### Görevi gözden geçirin

Göreviniz, bu modülleri daha sonra ana iş akışında birleştireceğimiz iki ayrı iş akışında bir araya getirmektir:

- İsimleri doğrulayan, selamlamalar oluşturan ve zaman damgaları ekleyen bir `GREETING_WORKFLOW`
- Metni büyük harfe dönüştüren ve tersine çeviren bir `TRANSFORM_WORKFLOW`

#### Hazırlık kontrol listesi

Dalmaya hazır olduğunuzu düşünüyor musunuz?

- [ ] Bu kursun amacını ve ön koşullarını anlıyorum
- [ ] Codespace'im çalışıyor
- [ ] Çalışma dizinini uygun şekilde ayarladım
- [ ] Görevi anlıyorum

Tüm kutuları işaretleyebiliyorsanız, başlamaya hazırsınız.

---

## 1. Selamlama İş Akışını Oluşturun

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

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Bu, 'Hello Nextflow' eğitiminde gördüklerinize benzer bir yapıya sahip, bağımsız olarak test edebileceğimiz eksiksiz bir iş akışıdır. Şimdi bunu deneyelim:

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

Bu beklendiği gibi çalışıyor, ancak birleştirilebilir hale getirmek için değiştirmemiz gereken birkaç şey var.

### 1.3. İş akışını birleştirilebilir hale getirin

Birleştirilebilir iş akışlarının 'Hello Nextflow' eğitiminde gördüklerinizden bazı farkları vardır:

- İş akışı bloğunun adlandırılması gerekir
- Girdiler `take:` anahtar kelimesi kullanılarak bildirilir
- İş akışı içeriği `main:` bloğunun içine yerleştirilir
- Çıktılar `emit:` anahtar kelimesi kullanılarak bildirilir

Selamlama iş akışını bu yapıya uyacak şekilde güncelleyelim. Kodu aşağıdaki gibi değiştirin:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Input channel with names

    main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
}
```

İş akışının artık adlandırıldığını ve bir `take:` ve `emit:` bloğuna sahip olduğunu görebilirsiniz, bunlar daha üst seviye bir iş akışı oluşturmak için kullanacağımız bağlantılardır.
İş akışı içeriği de `main:` bloğunun içine yerleştirilmiştir. Ayrıca `names_ch` girdi kanalı bildirimini kaldırdığımızı unutmayın, çünkü artık iş akışına argüman olarak geçiriliyor.

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

Bu size başka bir yeni kavramdan, 'giriş iş akışı'ndan bahsediyor. Giriş iş akışı, bir Nextflow betiğini çalıştırdığınızda çağrılan iş akışıdır. Varsayılan olarak, Nextflow mevcut olduğunda adsız bir iş akışını giriş iş akışı olarak kullanacaktır ve şimdiye kadar yaptığınız şey budur, şu şekilde başlayan iş akışı blokları ile:

```groovy title="hello.nf" linenums="1"
workflow {
```

Ancak selamlama iş akışımızda adsız bir iş akışı yok, bunun yerine adlandırılmış bir iş akışımız var:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Bu yüzden Nextflow bir hata verdi ve istediğimizi yapmadı.

`take:`/`emit:` sözdizimini iş akışını doğrudan çağırabilmek için eklemedik - onu diğer iş akışlarıyla birleştirebilmek için ekledik. Çözüm, adlandırılmış iş akışımızı içe aktaran ve çağıran adsız bir giriş iş akışına sahip bir ana betik oluşturmaktır.

### 1.4. Ana iş akışını oluşturun ve test edin

Şimdi `greeting` iş akışını içe aktaran ve kullanan bir ana iş akışı oluşturacağız.

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

Bu dosyadaki iş akışı girişimizin adsız olduğuna dikkat edin, bunun nedeni onu giriş iş akışı olarak kullanacak olmamızdır.

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
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Çalışıyor! Adlandırılmış selamlama iş akışını, adsız bir giriş `workflow` bloğuna sahip bir ana iş akışında sarmaladık. Ana iş akışı, `GREETING_WORKFLOW` iş akışını neredeyse (tam olarak değil) bir süreç gibi kullanıyor ve `names` kanalını argüman olarak geçiriyor.

### Özet

Bu bölümde, birkaç önemli kavram öğrendiniz:

- **Adlandırılmış İş Akışları**: İçe aktarılabilen ve yeniden kullanılabilen adlandırılmış bir iş akışı (`GREETING_WORKFLOW`) oluşturma
- **İş Akışı Arayüzleri**: Birleştirilebilir bir iş akışı oluşturmak için `take:` ile net girdiler ve `emit:` ile çıktılar tanımlama
- **Giriş Noktaları**: Nextflow'un bir betiği çalıştırmak için adsız bir giriş iş akışına ihtiyaç duyduğunu anlama
- **İş Akışı Birleştirme**: Adlandırılmış bir iş akışını başka bir iş akışı içinde içe aktarma ve kullanma
- **İş Akışı Ad Alanları**: `.out` ad alanını kullanarak iş akışı çıktılarına erişme (`GREETING_WORKFLOW.out.greetings`)

Artık şunları yapan çalışan bir selamlama iş akışınız var:

- Girdi olarak bir isim kanalı alır
- Her ismi doğrular
- Her geçerli isim için bir selamlama oluşturur
- Selamlamalara zaman damgaları ekler
- Hem orijinal hem de zaman damgalı selamlamaları çıktı olarak sunar

Bu modüler yaklaşım, selamlama iş akışını bağımsız olarak test etmenize veya daha büyük boru hatlarında bir bileşen olarak kullanmanıza olanak tanır.

---

## 2. Dönüştürme İş Akışını Ekleyin

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
        input_ch         // Input channel with messages

    main:
        // Apply transformations in sequence
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Uppercase greetings
        reversed = reversed_ch  // Reversed uppercase greetings
}
```

Birleştirilebilir sözdiziminin açıklamasını burada tekrarlamayacağız, ancak adlandırılmış iş akışının yine `take:` ve `emit:` bloğu ile bildirildiğini ve iş akışı içeriğinin `main:` bloğunun içine yerleştirildiğini unutmayın.

### 2.3. Ana iş akışını güncelleyin

Her iki iş akışını da kullanmak için `main.nf` dosyasını güncelleyin:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Tam boru hattını çalıştırın:

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

Bu tersine çevrilmiş dosyalardan birine bakarsanız, selamlamanın büyük harfli versiyonunun tersine çevrildiğini göreceksiniz:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Özet

Artık şunları yapan eksiksiz bir boru hattınız olmalı:

- İsimleri selamlama iş akışından geçirir
- Zaman damgalı selamlamaları dönüştürme iş akışına besler
- Selamlamaların hem büyük harfli hem de tersine çevrilmiş versiyonlarını üretir

---

## Özet

Bu yan görevde, Nextflow'da iş akışı birleştirme kavramını keşfettik; bu, karmaşık boru hatlarını daha küçük, yeniden kullanılabilir bileşenlerden oluşturmamıza olanak tanır.

Bu modüler yaklaşım, monolitik boru hatlarına göre çeşitli avantajlar sunar:

- Her iş akışı bağımsız olarak geliştirilebilir, test edilebilir ve hata ayıklanabilir
- İş akışları farklı boru hatları arasında yeniden kullanılabilir
- Genel boru hattı yapısı daha okunabilir ve sürdürülebilir hale gelir
- Arayüzler tutarlı kaldığı sürece bir iş akışındaki değişiklikler diğerlerini etkilemez
- Giriş noktaları, boru hattınızın farklı bölümlerini gerektiği gibi çalıştıracak şekilde yapılandırılabilir

_Ancak iş akışlarını çağırmanın süreçleri çağırmaya biraz benzese de, aslında aynı şey olmadığını belirtmek önemlidir. Örneğin, bir iş akışını N boyutunda bir kanalla çağırarak N kez çalıştıramazsınız - iş akışına N boyutunda bir kanal geçirmeniz ve dahili olarak yinelemeniz gerekir._

Bu teknikleri kendi çalışmanızda uygulamak, sürdürülebilir ve ölçeklenebilir kalırken karmaşık biyoinformatik görevleri işleyebilen daha sofistike Nextflow boru hatları oluşturmanıza olanak tanıyacaktır.

### Temel desenler

1.  **İş akışı yapısı**: Her iş akışı için `take:` ve `emit:` sözdizimini kullanarak net girdiler ve çıktılar tanımladık, bileşenler arasında iyi tanımlanmış arayüzler oluşturduk ve iş akışı mantığını `main:` bloğu içinde sarmaladık.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **İş akışı içe aktarmaları:** İki bağımsız iş akışı modülü oluşturduk ve bunları include ifadeleriyle ana boru hattına aktardık.

    - Tek bir iş akışını içe aktarma

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Birden fazla iş akışını içe aktarma

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - İsim çakışmalarını önlemek için takma ad ile içe aktarma

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Giriş noktaları**: Nextflow, yürütmeye nereden başlayacağını bilmek için adsız bir giriş iş akışına ihtiyaç duyar. Bu giriş iş akışı, adlandırılmış iş akışlarınızı çağırır.

    - Adsız iş akışı (giriş noktası)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Adlandırılmış iş akışı (giriş iş akışından çağrılır)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **Veri akışını yönetme:** İş akışı çıktılarına ad alanı notasyonu (`WORKFLOW_NAME.out.channel_name`) kullanarak nasıl erişileceğini ve bunların diğer iş akışlarına nasıl geçirileceğini öğrendik.

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

[Yan Görevler menüsüne](./index.md) dönün veya listedeki bir sonraki konuya geçmek için sayfanın sağ alt köşesindeki düğmeye tıklayın.
