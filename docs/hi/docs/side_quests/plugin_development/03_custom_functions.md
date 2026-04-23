# भाग 3: कस्टम फ़ंक्शन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस सेक्शन के अंत तक, तुम्हारे plugin में कस्टम फ़ंक्शन होंगे, जो लोकल में build और install होंगे और एक असली वर्कफ़्लो में चलेंगे।

!!! tip "यहाँ से शुरू कर रहे हो?"

    अगर तुम इस भाग से जुड़ रहे हो, तो शुरुआती बिंदु के रूप में Part 2 का समाधान कॉपी करो:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. देखो कि टेम्पलेट ने क्या बनाया

अपने खुद के फ़ंक्शन लिखने से पहले, टेम्पलेट द्वारा बनाए गए उदाहरण फ़ंक्शन को देखो ताकि pattern समझ सको।

plugin डायरेक्टरी में जाओ:

```bash
cd nf-greeting
```

टेम्पलेट ने `GreetingExtension.groovy` नाम की एक फ़ाइल बनाई है जहाँ plugin फ़ंक्शन परिभाषित होते हैं।
शुरुआती बिंदु देखने के लिए इसे खोलो:

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
 * एक कस्टम फ़ंक्शन implement करता है जिसे
 * Nextflow scripts import कर सकती हैं।
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * दिए गए target को hello कहो।
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. वह class जिस पर तुम्हारा extension बनता है। Nextflow को तुम्हारे फ़ंक्शन पहचानने के लिए यह ज़रूरी है।
2. plugin लोड होने पर call होता है; initialization के लिए उपयोग करो
3. इस method को `include` के ज़रिए वर्कफ़्लो से callable बनाता है

टेम्पलेट में एक sample `sayHello` फ़ंक्शन शामिल है।
`@Function` annotation ही किसी method को Nextflow वर्कफ़्लो से callable बनाता है।
इसके बिना, method केवल plugin code के अंदर ही मौजूद रहती है।

Groovy (और Java) में, methods यह declare करती हैं कि वे क्या return करती हैं और उनके पैरामीटर किस type के हैं।
उदाहरण के लिए, `String reverseGreeting(String greeting)` एक ऐसी method declare करता है जो एक `String` पैरामीटर लेती है और एक `String` return करती है।
`void` keyword का मतलब है कि method कुछ भी return नहीं करती, जैसे ऊपर `sayHello` में।
यह Python या R से अलग है, जहाँ types को explicitly declare करने की ज़रूरत नहीं होती।

---

## 2. sayHello को reverseGreeting से बदलो

टेम्पलेट का `sayHello` फ़ंक्शन एक placeholder है।
इसे अपने खुद के फ़ंक्शन से बदलो ताकि plugin फ़ंक्शन लिखने, build करने और उपयोग करने का पूरा cycle देख सको।

`sayHello` method को बदलने के लिए `src/main/groovy/training/plugin/GreetingExtension.groovy` edit करो:

=== "बाद में"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * एक greeting string को उल्टा करो
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Method को Nextflow वर्कफ़्लो से callable बनाता है
    2. एक String लेता है, एक String return करता है
    3. Groovy की built-in string reversal method

=== "पहले"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * एक कस्टम फ़ंक्शन implement करता है जिसे
     * Nextflow scripts import कर सकती हैं।
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * दिए गए target को hello कहो।
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

इस फ़ंक्शन के मुख्य भाग:

- **`@Function`**: Method को Nextflow वर्कफ़्लो से callable बनाता है
- **`String reverseGreeting(String greeting)`**: एक String लेता है, एक String return करता है
- **`greeting.reverse()`**: Groovy की built-in string reversal method

!!! tip "Public और private methods"

    `@Function` के बिना methods Nextflow वर्कफ़्लो के लिए expose नहीं होतीं।
    तुम अपनी class में helper methods जोड़ सकते हो बिना इस चिंता के कि वे वर्कफ़्लो namespace में leak हो जाएंगी।

---

## 3. अपना plugin build और install करो

Plugin को build और install करो:

```bash
make install
```

!!! tip "अगर build fail हो जाए"

    error message को ध्यान से पढ़ो; इसमें आमतौर पर एक line number होता है और समस्या का वर्णन होता है।
    सामान्य कारण हैं syntax errors (missing bracket या quote), गलत spelled class names, और type mismatches।
    अगर तुम फँस गए हो, तो अपने code को examples के साथ character-by-character compare करो।

---

## 4. अपने फ़ंक्शन को एक वर्कफ़्लो में उपयोग करो

Plugin build और install हो गया है।
अगला कदम है `reverseGreeting` को एक वर्कफ़्लो में उपयोग करना ताकि verify हो सके कि यह end-to-end काम करता है।

pipeline डायरेक्टरी में वापस जाओ:

```bash
cd ..
```

`reverseGreeting` को import और उपयोग करने के लिए `greet.nf` edit करो:

=== "बाद में"

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

=== "पहले"

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

Pipeline चलाओ:

```bash
nextflow run greet.nf
```

??? example "आउटपुट"

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

तुम्हारा पहला कस्टम plugin फ़ंक्शन एक असली वर्कफ़्लो में काम कर रहा है।
वही `include { ... } from 'plugin/...'` pattern जो तुमने Part 1 में nf-hello और nf-schema के साथ उपयोग किया था, तुम्हारे अपने plugin के साथ भी काम करता है।

---

## 5. decorateGreeting जोड़ो

एक plugin कई फ़ंक्शन provide कर सकता है।
एक दूसरा फ़ंक्शन जोड़ो जो greeting को decorative markers से wrap करे; तुम इसे Part 6 में configurable बनाओगे।

`GreetingExtension.groovy` edit करो और `reverseGreeting` के बाद, class के closing brace से पहले `decorateGreeting` जोड़ो:

=== "बाद में"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * एक greeting string को उल्टा करो
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * greeting को celebratory markers से decorate करो
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Groovy string interpolation: `#!groovy ${...}` variable की value को string में insert करता है

=== "पहले"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * एक greeting string को उल्टा करो
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

यह फ़ंक्शन Groovy string interpolation (`"*** ${greeting} ***"`) का उपयोग करता है ताकि greeting variable को एक string के अंदर embed किया जा सके।

Build, install करो और वर्कफ़्लो update करो:

```bash
cd nf-greeting && make install && cd ..
```

`decorateGreeting` को भी import और उपयोग करने के लिए `greet.nf` update करो:

=== "बाद में"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // हमारे plugin से कस्टम फ़ंक्शन import करो
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // greeting को decorate करने के लिए हमारा कस्टम plugin फ़ंक्शन उपयोग करो
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // reverseGreeting फ़ंक्शन का उपयोग demonstrate करो
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. एक ही plugin के कई फ़ंक्शन के लिए अलग-अलग `include` statements की ज़रूरत होती है
    2. Plugin फ़ंक्शन process `script:` blocks के अंदर भी काम करते हैं

=== "पहले"

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

??? example "आउटपुट"

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

Plugin फ़ंक्शन process scripts (जैसे `SAY_HELLO` के अंदर `decorateGreeting`) और वर्कफ़्लो operations (जैसे `map` में `reverseGreeting`) दोनों में काम करते हैं।

---

## सारांश

तुमने सीखा कि:

- फ़ंक्शन `PluginExtensionPoint` subclasses में `@Function` annotation के साथ परिभाषित होते हैं
- `include` से import किए गए plugin फ़ंक्शन एक जैसे काम करते हैं, चाहे वे तुम्हारे अपने plugin से हों या किसी existing plugin से
- Plugin फ़ंक्शन process scripts और वर्कफ़्लो operations दोनों में काम करते हैं

---

## आगे क्या है?

तुम्हारे फ़ंक्शन काम कर रहे हैं, लेकिन अभी तक तुमने इसे केवल पूरी pipeline चलाकर और आउटपुट को आँखों से देखकर verify किया है।
यह तरीका scale नहीं होता: जैसे-जैसे तुम और फ़ंक्शन जोड़ते हो, तुम्हें एक तेज़ तरीके की ज़रूरत होती है यह जाँचने के लिए कि हर एक सही तरह से behave करता है, खासकर बदलाव करने के बाद।
अगला सेक्शन unit tests से परिचित कराता है, जो तुम्हें pipeline चलाए बिना अलग-अलग फ़ंक्शन को automatically verify करने देते हैं।

[Part 4 पर जारी रखो :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
