# भाग 5: Trace Observers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Trace observers तुम्हारे plugin को workflow events पर प्रतिक्रिया देने देते हैं, जैसे कि कोई कार्य पूरा होना, कोई फ़ाइल publish होना, या पाइपलाइन समाप्त होना।
इससे custom reports, Slack notifications, metrics collection, या बाहरी monitoring systems के साथ integration जैसे use cases संभव होते हैं।
इस section में, तुम एक ऐसा observer बनाओगे जो पूरे हुए कार्यों की गिनती करे और एक सारांश print करे।

!!! tip "यहाँ से शुरू कर रहे हो?"

    अगर तुम इस भाग से जुड़ रहे हो, तो शुरुआती बिंदु के रूप में Part 4 का समाधान copy करो:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. मौजूदा trace observer को समझना

जब तुमने पाइपलाइन चलाई थी, तब "Pipeline is starting!" संदेश तुम्हारे plugin की `GreetingObserver` class से आया था।

observer का कोड देखो:

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
 * एक observer implement करता है जो nextflow execution events पर
 * custom logic implement करने देता है।
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

1. workflow lifecycle events में hook करने के लिए interface
2. workflow शुरू होने पर call होता है; config access करने के लिए session प्राप्त करता है
3. workflow सफलतापूर्वक समाप्त होने पर call होता है

यहाँ दो बातें ध्यान देने योग्य हैं:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` Nextflow द्वारा परिभाषित एक interface है। अगर तुम्हारी class इस interface को implement करती है, तो Nextflow इसमें hook कर सकता है और events होने पर तुम्हारे methods को call कर सकता है।
2. **`@Override`**: `TraceObserver` interface `onFlowCreate` और `onFlowComplete` जैसे methods परिभाषित करता है। जब तुम इन नामों के methods लिखते हो और `@Override` annotation जोड़ते हो, तो Nextflow उन्हें उचित समय पर call करता है। जो methods तुम override नहीं करते, उन्हें ignore किया जाता है।

लिखते समय जिन lifecycle events में तुम hook कर सकते हो, उनका पूरा सेट है:

| Method              | कब call होता है              |
| ------------------- | ---------------------------- |
| `onFlowCreate`      | Workflow शुरू होता है        |
| `onFlowComplete`    | Workflow समाप्त होता है      |
| `onProcessStart`    | कोई कार्य execution शुरू करता है |
| `onProcessComplete` | कोई कार्य समाप्त होता है     |
| `onProcessCached`   | cached कार्य पुनः उपयोग होता है |
| `onFilePublish`     | कोई फ़ाइल publish होती है    |

पूरी सूची के लिए, Nextflow source में [TraceObserver interface](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) देखो।

---

## 2. एक task counter observer जोड़ना

लक्ष्य है एक ऐसा observer बनाना जो पूरे हुए कार्यों की गिनती करे और अंत में एक सारांश print करे।
किसी plugin में नया observer जोड़ने के लिए दो चीज़ें चाहिए: observer class लिखना, और उसे factory में register करना ताकि Nextflow उसे load करे।

### 2.1. एक minimal observer बनाना

एक नई फ़ाइल बनाओ:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

सबसे सरल observer से शुरू करो जो किसी भी कार्य के पूरा होने पर एक संदेश print करे:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * कार्य पूर्णता पर प्रतिक्रिया देने वाला observer
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. आवश्यक classes import करो: `TraceObserver`, `TaskHandler`, और `TraceRecord`
2. एक class बनाओ जो `TraceObserver` को `implements` करे
3. कार्य समाप्त होने पर कोड चलाने के लिए `onProcessComplete` को override करो

यह न्यूनतम आवश्यकता है:

- आवश्यक classes import करो (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- एक class बनाओ जो `TraceObserver` को `implements` करे
- कार्य समाप्त होने पर कुछ करने के लिए `onProcessComplete` को override करो

### 2.2. Observer को register करना

`GreetingFactory` observers बनाता है।
इसे देखो:

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

नया observer जोड़ने के लिए `GreetingFactory.groovy` को edit करो:

=== "बाद में"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "पहले"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Groovy list syntax"

    हमने Java-style `List.<TraceObserver>of(...)` को Groovy के सरल list literal `[...]` से बदल दिया है।
    दोनों एक `Collection` return करते हैं, लेकिन Groovy syntax कई items जोड़ते समय अधिक पठनीय है।

### 2.3. Build, install, और test करना

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "`-ansi-log false` क्यों?"

    डिफ़ॉल्ट रूप से, Nextflow का ANSI progress display progress का एक साफ़, अपडेट होता दृश्य दिखाने के लिए पिछली lines को overwrite करता है।
    इसका मतलब है कि तुम केवल *अंतिम* task count देखोगे, intermediate संदेश नहीं।

    `-ansi-log false` का उपयोग इस व्यवहार को disable करता है और सभी आउटपुट क्रमिक रूप से दिखाता है, जो execution के दौरान संदेश print करने वाले observers का परीक्षण करते समय आवश्यक है।

तुम्हें "✓ Task completed!" पाँच बार print होते दिखना चाहिए (प्रत्येक कार्य के लिए एक बार), मौजूदा पाइपलाइन आउटपुट के साथ मिलाकर:

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

observer काम कर रहा है।
हर बार जब कोई कार्य समाप्त होता है, Nextflow `onProcessComplete` को call करता है, और हमारी implementation एक संदेश print करती है।

??? exercise "संदेश को customize करो"

    `onProcessComplete` में संदेश को अपनी पसंद के किसी संदेश से बदलने की कोशिश करो, rebuild करो, और फिर से चलाओ।
    यह पुष्टि करता है कि observers के लिए पूरा edit-build-run cycle काम करता है।

### 2.4. गिनती का logic जोड़ना

minimal observer यह साबित करता है कि hook काम करता है, लेकिन यह कुछ track नहीं करता।

एक class में variables (जिन्हें fields या instance variables कहते हैं) हो सकते हैं जो object के जीवनकाल तक बने रहते हैं।
इसका मतलब है कि एक observer पाइपलाइन run के दौरान कई events में state accumulate कर सकता है।

अगला version एक counter variable (`taskCount`) जोड़ता है जो शून्य से शुरू होता है।
हर बार जब कोई कार्य पूरा होता है, counter एक बढ़ जाता है।
जब पूरा workflow समाप्त होता है, observer अंतिम कुल print करता है।

highlighted changes के साथ `TaskCounterObserver.groovy` को update करो:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * पूरे हुए कार्यों की गिनती करने वाला observer
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

1. `taskCount` एक variable है जो observer object से संबंधित है। यह method calls के बीच अपना मान बनाए रखता है, इसलिए यह पूरे workflow run में एक count accumulate कर सकता है। `private` का मतलब है कि केवल यह class इसे access कर सकती है।
2. `taskCount++` counter में एक जोड़ता है। यह line हर बार चलती है जब कोई कार्य पूरा होता है, इसलिए workflow आगे बढ़ने के साथ count बढ़ता है।
3. `onFlowComplete` एक दूसरा lifecycle hook है। यह workflow समाप्त होने पर एक बार चलता है, जो इसे सारांश print करने के लिए एक अच्छी जगह बनाता है।

संक्षेप में:

- `taskCount` method calls में बना रहता है, पूरे run में एक count accumulate करता है
- `onProcessComplete` counter को बढ़ाता है और हर बार कार्य समाप्त होने पर running total print करता है
- `onFlowComplete` अंत में एक बार चलता है, अंतिम count print करता है

Rebuild और test करो:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "आउटपुट"

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

    Counter संदेश कार्य submissions के साथ मिले हुए हैं क्योंकि observers कार्य पूरे होने पर चलते हैं।

---

## 3. Published फ़ाइलों को track करना

Observer फ़ाइलें publish होने पर भी प्रतिक्रिया दे सकता है।
`onFilePublish` method destination और source paths प्राप्त करता है, जिनका उपयोग तुम published आउटपुट को log, validate, या process करने के लिए कर सकते हो।

### 3.1. एक publish डायरेक्टरी जोड़ना

पहले, `greet.nf` को update करो ताकि `SAY_HELLO` प्रोसेस अपनी आउटपुट फ़ाइलें publish करे:

=== "बाद में"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // हमारे custom plugin function का उपयोग करके greeting को decorate करो
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "पहले"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // हमारे custom plugin function का उपयोग करके greeting को decorate करो
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. onFilePublish method जोड़ना

`TaskCounterObserver.groovy` में एक `onFilePublish` method और आवश्यक import जोड़ो:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * पूरे हुए कार्यों की गिनती करने वाला observer
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

### 3.3. Build और test करना

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

तुम्हें task counter आउटपुट के साथ प्रत्येक आउटपुट फ़ाइल के लिए "Published:" संदेश दिखने चाहिए:

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

`onFilePublish` method हर बार fire होता है जब Nextflow `results` डायरेक्टरी में कोई फ़ाइल publish करता है।
यह pattern audit logs बनाने, downstream actions trigger करने, या produce होते समय आउटपुट validate करने के लिए उपयोगी है।

---

## सारांश

तुमने सीखा कि:

- Trace observers workflow lifecycle events जैसे `onFlowCreate`, `onProcessComplete`, `onFilePublish`, और `onFlowComplete` में hook करते हैं
- `TraceObserver` implement करके और उन्हें Factory में register करके observers बनाओ
- Observers events में state accumulate करने के लिए instance variables रख सकते हैं
- Observers custom logging, metrics collection, notifications, और reporting के लिए उपयोगी हैं

---

## आगे क्या है?

Task counter काम करता है, लेकिन यह हमेशा चालू रहता है।
एक real plugin में, users को plugin source code edit किए बिना `nextflow.config` से features enable या disable करने, या व्यवहार adjust करने में सक्षम होना चाहिए।
अगला section दिखाता है कि अपने observer को configurable कैसे बनाएं और अपना तैयार plugin दूसरों के साथ कैसे share करें।

[भाग 6 पर जारी रखो :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
