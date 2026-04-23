# Bölüm 3: Özel Fonksiyonlar

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bu bölümün sonunda, eklentinizde özel fonksiyonlar oluşturmuş, bunları yerel olarak derleyip kurmuş ve gerçek bir iş akışında çalıştırmış olacaksınız.

!!! tip "Buradan mı başlıyorsunuz?"

    Bu bölüme doğrudan katılıyorsanız, başlangıç noktası olarak kullanmak üzere Bölüm 2'nin çözümünü kopyalayın:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Şablonun oluşturduklarını inceleyin

Kendi fonksiyonlarınızı yazmadan önce, kalıbı anlamak için şablonun oluşturduğu örnek fonksiyona bakın.

Eklenti dizinine geçin:

```bash
cd nf-greeting
```

Şablon, eklenti fonksiyonlarının tanımlandığı `GreetingExtension.groovy` adlı bir dosya oluşturdu.
Başlangıç noktasını görmek için dosyayı açın:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

```groovy title="Output" hl_lines="29 40-43"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Nextflow betikleri tarafından içe aktarılabilen
 * özel bir fonksiyon uygular.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Verilen hedefe merhaba de.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. Eklentinizin üzerine inşa edildiği sınıf. Nextflow'un fonksiyonlarınızı tanıması için bu gereklidir.
2. Eklenti yüklendiğinde çağrılır; başlatma işlemleri için kullanın
3. Bu metodu iş akışlarından `include` aracılığıyla çağrılabilir hale getirir

Şablon, örnek bir `sayHello` fonksiyonu içermektedir.
`@Function` anotasyonu, bir metodun Nextflow iş akışlarından çağrılabilmesini sağlayan unsurdur.
Bu anotasyon olmadan, metot yalnızca eklenti kodu içinde var olur.

Groovy'de (ve Java'da) metodlar, döndürdükleri türü ve parametrelerinin türlerini açıkça belirtir.
Örneğin, `String reverseGreeting(String greeting)`, bir `String` parametresi alan ve `String` döndüren bir metot tanımlar.
`void` anahtar sözcüğü, metodun yukarıdaki `sayHello` örneğinde olduğu gibi hiçbir şey döndürmediği anlamına gelir.
Bu durum, türlerin açıkça belirtilmesine gerek olmadığı Python veya R'dan farklıdır.

---

## 2. sayHello'yu reverseGreeting ile değiştirin

Şablonun `sayHello` fonksiyonu bir yer tutucudur.
Bir eklenti fonksiyonu yazma, derleme ve kullanma döngüsünün tamamını görmek için onu kendi fonksiyonunuzla değiştirin.

`sayHello` metodunu değiştirmek için `src/main/groovy/training/plugin/GreetingExtension.groovy` dosyasını düzenleyin:

=== "Sonra"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Bir selamlama dizesini tersine çevirir
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Metodu Nextflow iş akışlarından çağrılabilir hale getirir
    2. Bir String alır, bir String döndürür
    3. Groovy'nin yerleşik dize tersine çevirme metodu

=== "Önce"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Nextflow betikleri tarafından içe aktarılabilen
     * özel bir fonksiyon uygular.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Verilen hedefe merhaba de.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Bu fonksiyonun temel bileşenleri:

- **`@Function`**: Metodu Nextflow iş akışlarından çağrılabilir hale getirir
- **`String reverseGreeting(String greeting)`**: Bir String alır, bir String döndürür
- **`greeting.reverse()`**: Groovy'nin yerleşik dize tersine çevirme metodu

!!! tip "Public ve private metodlar"

    `@Function` anotasyonu olmayan metodlar Nextflow iş akışlarına açık değildir.
    İş akışı isim alanına sızma konusunda endişelenmeden sınıfınıza yardımcı metodlar ekleyebilirsiniz.

---

## 3. Eklentinizi derleyin ve kurun

Eklentiyi derleyip kurun:

```bash
make install
```

!!! tip "Derleme başarısız olursa"

    Hata mesajını dikkatlice okuyun; genellikle bir satır numarası içerir ve sorunu açıklar.
    Yaygın nedenler arasında sözdizimi hataları (eksik parantez veya tırnak işareti), yanlış yazılmış sınıf adları ve tür uyumsuzlukları sayılabilir.
    Takılıp kalırsanız, kodunuzu örneklerle karakter karakter karşılaştırın.

---

## 4. Fonksiyonunuzu bir iş akışında kullanın

Eklenti derlenmiş ve kurulmuştur.
Bir sonraki adım, uçtan uca çalıştığını doğrulamak için `reverseGreeting`'i bir iş akışında kullanmaktır.

Pipeline dizinine geri dönün:

```bash
cd ..
```

`reverseGreeting`'i içe aktarmak ve kullanmak için `greet.nf` dosyasını düzenleyin:

=== "Sonra"

    ```groovy title="greet.nf" hl_lines="4 23-25" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Önce"

    ```groovy title="greet.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Pipeline'ı çalıştırın:

```bash
nextflow run greet.nf
```

??? example "Çıktı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Output: Hello
    Output: Bonjour
    Output: Holà
    Output: Ciao
    Output: Hallo
    Pipeline complete! 👋
    ```

İlk özel eklenti fonksiyonunuz gerçek bir iş akışında çalışıyor.
Bölüm 1'de nf-hello ve nf-schema ile kullandığınız `include { ... } from 'plugin/...'` kalıbı, kendi eklentinizle de aynı şekilde çalışır.

---

## 5. decorateGreeting'i ekleyin

Bir eklenti birden fazla fonksiyon sağlayabilir.
Bir selamlamayı dekoratif işaretlerle saran ikinci bir fonksiyon ekleyin; bunu Bölüm 6'da yapılandırılabilir hale getireceksiniz.

Sınıfın kapanış parantezinden önce, `reverseGreeting`'in ardına `decorateGreeting` eklemek için `GreetingExtension.groovy` dosyasını düzenleyin:

=== "Sonra"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Bir selamlama dizesini tersine çevirir
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Bir selamlamayı kutlama işaretleriyle süsler
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Groovy dize enterpolasyonu: `#!groovy ${...}` değişkenin değerini dizeye ekler

=== "Önce"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Bir selamlama dizesini tersine çevirir
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Bu fonksiyon, selamlama değişkenini bir dizeye gömmek için Groovy dize enterpolasyonunu (`"*** ${greeting} ***"`) kullanır.

Derleyin, kurun ve iş akışını güncelleyin:

```bash
cd nf-greeting && make install && cd ..
```

`decorateGreeting`'i de içe aktarmak ve kullanmak için `greet.nf` dosyasını güncelleyin:

=== "Sonra"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Eklentimizden özel fonksiyonları içe aktar
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Selamlamayı süslemek için özel eklenti fonksiyonumuzu kullan
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // reverseGreeting fonksiyonunun kullanımını göster
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Aynı eklentiden birden fazla fonksiyon için ayrı `include` ifadeleri gerekir
    2. Eklenti fonksiyonları süreç `script:` bloklarında da çalışır

=== "Önce"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 12 15 28"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

```bash
nextflow run greet.nf
```

??? example "Çıktı"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    ```

Eklenti fonksiyonları hem süreç betiklerinde (`SAY_HELLO` içindeki `decorateGreeting` gibi) hem de iş akışı operasyonlarında (`map` içindeki `reverseGreeting` gibi) çalışır.

---

## Özetle

Şunları öğrendiniz:

- Fonksiyonlar, `PluginExtensionPoint` alt sınıflarında `@Function` anotasyonu ile tanımlanır
- `include` ile içe aktarılan eklenti fonksiyonları, kendi eklentinizden veya mevcut bir eklentiden gelip gelmediğine bakılmaksızın aynı şekilde çalışır
- Eklenti fonksiyonları hem süreç betiklerinde hem de iş akışı operasyonlarında kullanılabilir

---

## Sırada ne var?

Fonksiyonlarınız çalışıyor; ancak şimdiye kadar bunu yalnızca pipeline'ın tamamını çalıştırıp çıktıyı gözle kontrol ederek doğruladınız.
Bu yaklaşım ölçeklenmez: daha fazla fonksiyon ekledikçe, özellikle değişiklik yaptıktan sonra her birinin doğru davrandığını denetlemenin daha hızlı bir yoluna ihtiyaç duyarsınız.
Bir sonraki bölüm, pipeline çalıştırmadan tek tek fonksiyonları otomatik olarak doğrulamanızı sağlayan birim testlerini tanıtmaktadır.

[Bölüm 4'e geçin :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
