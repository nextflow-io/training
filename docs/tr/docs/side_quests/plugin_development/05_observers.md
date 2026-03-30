# Bölüm 5: Trace Observer'lar

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Yapay zeka destekli çeviri - [daha fazla bilgi ve iyileştirme önerileri](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Trace observer'lar, eklentinizin bir görevin tamamlanması, bir dosyanın yayımlanması veya pipeline'ın bitmesi gibi iş akışı olaylarına yanıt vermesini sağlar.
Bu sayede özel raporlar, Slack bildirimleri, metrik toplama veya harici izleme sistemleriyle entegrasyon gibi kullanım senaryoları mümkün olur.
Bu bölümde, tamamlanan görevleri sayan ve bir özet yazdıran bir observer oluşturacaksınız.

!!! tip "Buradan mı başlıyorsunuz?"

    Bu bölüme doğrudan katılıyorsanız, başlangıç noktası olarak Bölüm 4'ün çözümünü kopyalayın:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Mevcut trace observer'ı anlamak

Pipeline'ı çalıştırdığınızda görünen "Pipeline is starting!" mesajı, eklentinizdeki `GreetingObserver` sınıfından geliyordu.

Observer koduna bakın:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingObserver.groovy
```

```groovy title="Output" hl_lines="30 32-34 37-39"
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
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.trace.TraceObserver

/**
 * Özel mantığın Nextflow yürütme olaylarına bağlanmasını
 * sağlayan bir observer uygular.
 */
@Slf4j
@CompileStatic
class GreetingObserver implements TraceObserver {    // (1)!

    @Override
    void onFlowCreate(Session session) {            // (2)!
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {                         // (3)!
        println "Pipeline complete! 👋"
    }
}
```

1. İş akışı yaşam döngüsü olaylarına bağlanmak için kullanılan arayüz
2. İş akışı başladığında çağrılır; yapılandırmaya erişmek için oturumu alır
3. İş akışı başarıyla tamamlandığında çağrılır

Burada dikkat edilmesi gereken iki nokta vardır:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver`, Nextflow tarafından tanımlanmış bir arayüzdür. Sınıfınız bu arayüzü uygularsa, Nextflow olaylar gerçekleştiğinde metodlarınızı çağırabilir.
2. **`@Override`**: `TraceObserver` arayüzü, `onFlowCreate` ve `onFlowComplete` gibi metodlar tanımlar. Bu isimlerle metodlar yazıp `@Override` anotasyonunu eklediğinizde, Nextflow bunları uygun zamanda çağırır. Override etmediğiniz metodlar yok sayılır.

Yazım sırasında bağlanabileceğiniz yaşam döngüsü olaylarının tam listesi şöyledir:

| Metod               | Ne zaman çağrılır                  |
| ------------------- | ---------------------------------- |
| `onFlowCreate`      | İş akışı başladığında              |
| `onFlowComplete`    | İş akışı tamamlandığında           |
| `onProcessStart`    | Bir görev yürütmeye başladığında   |
| `onProcessComplete` | Bir görev tamamlandığında          |
| `onProcessCached`   | Önbelleğe alınmış görev kullanıldığında |
| `onFilePublish`     | Bir dosya yayımlandığında          |

Tam liste için Nextflow kaynak kodundaki [TraceObserver arayüzüne](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) bakın.

---

## 2. Görev sayacı observer'ı ekleme

Amaç, tamamlanan görevleri sayan ve sonunda bir özet yazdıran bir observer oluşturmaktır.
Bir eklentiye yeni bir observer eklemek iki şey gerektirir: observer sınıfını yazmak ve Nextflow'un onu yükleyebilmesi için factory'ye kaydetmek.

### 2.1. Minimal bir observer oluşturma

Yeni bir dosya oluşturun:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Herhangi bir görev tamamlandığında mesaj yazdıran mümkün olan en basit observer ile başlayın:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Görev tamamlanma olaylarına yanıt veren observer
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Gerekli sınıfları içe aktarın: `TraceObserver`, `TaskHandler` ve `TraceRecord`
2. `TraceObserver`'ı uygulayan (`implements`) bir sınıf oluşturun
3. Bir görev tamamlandığında kod çalıştırmak için `onProcessComplete`'i override edin

Bu, gereken asgari yapıdır:

- Gerekli sınıfları içe aktarın (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- `TraceObserver`'ı uygulayan bir sınıf oluşturun
- Bir görev tamamlandığında bir şeyler yapmak için `onProcessComplete`'i override edin

### 2.2. Observer'ı kaydetme

`GreetingFactory`, observer'ları oluşturur.
Bir göz atın:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingFactory.groovy
```

```groovy title="Output" hl_lines="25 27-29"
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
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }

}
```

Yeni observer'ı eklemek için `GreetingFactory.groovy`'yi düzenleyin:

=== "Sonra"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Önce"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Groovy liste sözdizimi"

    Java tarzı `List.<TraceObserver>of(...)` ifadesini Groovy'nin daha sade liste değişmezi `[...]` ile değiştirdik.
    Her ikisi de bir `Collection` döndürür; ancak birden fazla öğe eklerken Groovy sözdizimi daha okunabilirdir.

### 2.3. Derleme, kurulum ve test

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Neden `-ansi-log false`?"

    Nextflow, varsayılan olarak ANSI ilerleme görüntüsünü kullanır ve önceki satırların üzerine yazarak ilerlemenin temiz ve güncel bir görünümünü sunar.
    Bu, yalnızca *son* görev sayısını göreceğiniz anlamına gelir; ara mesajları göremezsiniz.

    `-ansi-log false` kullanmak bu davranışı devre dışı bırakır ve tüm çıktıyı sıralı olarak gösterir; bu da yürütme sırasında mesaj yazdıran observer'ları test ederken önemlidir.

"✓ Task completed!" ifadesinin beş kez (her görev için bir kez) mevcut pipeline çıktısıyla iç içe geçmiş şekilde yazdırıldığını görmelisiniz:

```console title="Output (partial)"
...
[9b/df7630] Submitted process > SAY_HELLO (4)
Decorated: *** Hello ***
✓ Task completed!
✓ Task completed!
Decorated: *** Holà ***
✓ Task completed!
...
Pipeline complete! 👋
```

Observer çalışıyor.
Bir görev her tamamlandığında, Nextflow `onProcessComplete`'i çağırır ve bizim uygulamamız bir mesaj yazdırır.

??? exercise "Mesajı özelleştirin"

    `onProcessComplete` içindeki mesajı kendi seçtiğiniz bir şeyle değiştirmeyi deneyin, yeniden derleyin ve çalıştırın.
    Bu, observer'lar için tam düzenleme-derleme-çalıştırma döngüsünün çalıştığını doğrular.

### 2.4. Sayma mantığı ekleme

Minimal observer, kancaların çalıştığını kanıtlar; ancak herhangi bir şeyi takip etmez.

Bir sınıf, nesnenin ömrü boyunca kalıcı olan değişkenler (alan veya örnek değişkeni olarak adlandırılır) tutabilir.
Bu, bir observer'ın pipeline çalışması sırasında birden fazla olay boyunca durum biriktirebileceği anlamına gelir.

Bir sonraki sürüm, sıfırdan başlayan bir sayaç değişkeni (`taskCount`) ekler.
Bir görev her tamamlandığında sayaç bir artar.
Tüm iş akışı tamamlandığında, observer son toplamı yazdırır.

`TaskCounterObserver.groovy`'yi vurgulanan değişikliklerle güncelleyin:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Tamamlanan görevleri sayan observer
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0                // (1)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++                          // (2)!
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFlowComplete() {                  // (3)!
        println "📈 Final task count: ${taskCount}"
    }
}
```

1. `taskCount`, observer nesnesine ait bir değişkendir. Metod çağrıları arasında değerini korur; bu sayede tüm iş akışı çalışması boyunca bir sayım biriktirebilir. `private`, yalnızca bu sınıfın ona erişebileceği anlamına gelir.
2. `taskCount++` sayaca bir ekler. Bu satır her görev tamamlandığında çalışır; böylece iş akışı ilerledikçe sayım artar.
3. `onFlowComplete`, ikinci bir yaşam döngüsü kancasıdır. İş akışı tamamlandığında bir kez çalışır ve özet yazdırmak için uygun bir yerdir.

Özetle:

- `taskCount`, metod çağrıları arasında kalıcıdır ve tüm çalışma boyunca bir sayım biriktirir
- `onProcessComplete`, bir görev her tamamlandığında sayacı artırır ve anlık toplamı yazdırır
- `onFlowComplete`, sonunda bir kez çalışır ve son sayımı yazdırır

Yeniden derleyin ve test edin:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Çıktı"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `greet.nf` [pensive_engelbart] DSL2 - revision: 85fefd90d0
    Pipeline is starting! 🚀
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    [be/bd8e72] Submitted process > SAY_HELLO (2)
    [5b/d24c2b] Submitted process > SAY_HELLO (1)
    [14/1f9dbe] Submitted process > SAY_HELLO (3)
    Decorated: *** Bonjour ***
    Decorated: *** Hello ***
    [85/a6b3ad] Submitted process > SAY_HELLO (4)
    📊 Tasks completed so far: 1
    📊 Tasks completed so far: 2
    Decorated: *** Holà ***
    📊 Tasks completed so far: 3
    Decorated: *** Ciao ***
    [3c/be6686] Submitted process > SAY_HELLO (5)
    📊 Tasks completed so far: 4
    Decorated: *** Hallo ***
    📊 Tasks completed so far: 5
    Pipeline complete! 👋
    📈 Final task count: 5
    ```

    Sayaç mesajları, görev gönderimleriyle iç içe geçmiştir; çünkü observer'lar görevler tamamlandıkça çalışır.

---

## 3. Yayımlanan dosyaları takip etme

Observer, dosyalar yayımlandığında da yanıt verebilir.
`onFilePublish` metodu, hedef ve kaynak yollarını alır; bunları yayımlanan çıktıları günlüğe kaydetmek, doğrulamak veya işlemek için kullanabilirsiniz.

### 3.1. Yayımlama dizini ekleme

Önce, `SAY_HELLO` sürecinin çıktı dosyalarını yayımlaması için `greet.nf`'yi güncelleyin:

=== "Sonra"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Selamlamayı süslemek için özel eklenti fonksiyonumuzu kullan
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Önce"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Selamlamayı süslemek için özel eklenti fonksiyonumuzu kullan
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. onFilePublish metodunu ekleme

`TaskCounterObserver.groovy`'ye bir `onFilePublish` metodu ve gerekli içe aktarmayı ekleyin:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Tamamlanan görevleri sayan observer
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "📁 Published: ${destination.fileName}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
```

### 3.3. Derleme ve test

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Her çıktı dosyası için görev sayacı çıktısıyla birlikte "Published:" mesajlarını görmelisiniz:

```console title="Output (partial)"
...
📊 Tasks completed so far: 1
📁 Published: greeting.txt
📊 Tasks completed so far: 2
📁 Published: greeting.txt
...
📈 Final task count: 5
Pipeline complete! 👋
```

`onFilePublish` metodu, Nextflow her `results` dizinine bir dosya yayımladığında tetiklenir.
Bu kalıp, denetim günlükleri oluşturmak, aşağı akış eylemlerini tetiklemek veya çıktılar üretildikçe doğrulamak için kullanışlıdır.

---

## Özetle

Şunları öğrendiniz:

- Trace observer'lar, `onFlowCreate`, `onProcessComplete`, `onFilePublish` ve `onFlowComplete` gibi iş akışı yaşam döngüsü olaylarına bağlanır
- `TraceObserver`'ı uygulayarak ve bir Factory'ye kaydederek observer oluşturulur
- Observer'lar, olaylar arasında durum biriktirmek için örnek değişkenleri tutabilir
- Observer'lar özel günlükleme, metrik toplama, bildirimler ve raporlama için kullanışlıdır

---

## Sırada ne var?

Görev sayacı çalışıyor; ancak her zaman açık durumda.
Gerçek bir eklentide, kullanıcıların eklenti kaynak kodunu düzenlemeden `nextflow.config` üzerinden özellikleri etkinleştirip devre dışı bırakabilmesi veya davranışı ayarlayabilmesi gerekir.
Bir sonraki bölüm, observer'ınızı yapılandırılabilir hale getirmeyi ve tamamlanmış eklentinizi başkalarıyla paylaşmayı göstermektedir.

[Bölüm 6'ya devam edin :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
