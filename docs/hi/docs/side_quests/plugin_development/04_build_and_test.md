# भाग 4: परीक्षण

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Plugins स्वतंत्र सॉफ़्टवेयर हैं जिन पर पाइपलाइन डेवलपर्स को भरोसा करना होता है।
प्रत्येक फ़ीचर को स्वतंत्र रूप से, पाइपलाइन के बाहर, परखना यह सुनिश्चित करता है कि plugin किसी के वर्कफ़्लो में integrate करने से पहले सही तरीके से काम करे।
इस सेक्शन में, तुम Spock testing framework का उपयोग करके टेस्ट लिखोगे और चलाओगे।

!!! tip "यहाँ से शुरू कर रहे हो?"

    अगर तुम इस भाग से जुड़ रहे हो, तो अपने शुरुआती बिंदु के रूप में भाग 3 का समाधान कॉपी करो:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    फिर plugin डायरेक्टरी में जाओ:

    ```bash
    cd nf-greeting
    ```

सुनिश्चित करो कि तुम plugin डायरेक्टरी में हो:

```bash
cd nf-greeting
```

---

## 1. परीक्षण क्यों करें?

एक सफल build का मतलब है कि कोड compile हो गया, लेकिन यह नहीं जाँचता कि वह अपेक्षित रूप से काम करता है या नहीं।
Unit tests कोड के छोटे टुकड़े होते हैं जो स्वचालित रूप से जाँचते हैं कि तुम्हारे फ़ंक्शन किसी दिए गए इनपुट के लिए सही आउटपुट देते हैं या नहीं।
उदाहरण के लिए, एक टेस्ट यह जाँच सकता है कि `#!groovy reverseGreeting("Hello")` `"olleH"` लौटाता है।

टेस्ट मूल्यवान हैं क्योंकि:

- ये bugs को users से पहले पकड़ते हैं
- ये तुम्हें बिना कुछ तोड़े बदलाव करने का आत्मविश्वास देते हैं
- ये दस्तावेज़ीकरण के रूप में काम करते हैं जो दिखाते हैं कि फ़ंक्शन कैसे उपयोग किए जाने चाहिए

---

## 2. Spock टेस्ट को समझना

Plugin template [Spock](https://spockframework.org/) का उपयोग करता है, जो Groovy के लिए एक testing framework है।
Spock पहले से ही प्रोजेक्ट में कॉन्फ़िगर है (`build.gradle` के ज़रिए), इसलिए तुम्हें कुछ भी जोड़ने की ज़रूरत नहीं है।

अगर तुमने पहले testing tools का उपयोग किया है (जैसे Python में `pytest` या R में `testthat`), तो Spock वही भूमिका निभाता है: तुम छोटे फ़ंक्शन लिखते हो जो तुम्हारे कोड को ज्ञात इनपुट के साथ कॉल करते हैं और आउटपुट की जाँच करते हैं।
अंतर यह है कि Spock लेबल किए गए blocks (`given:`, `expect:`, `when:`, `then:`) का उपयोग करता है जो एक Nextflow process या workflow के समान हैं।

यहाँ बुनियादी संरचना है:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **उद्धरण चिह्नों में टेस्ट का नाम**: बताता है कि टेस्ट क्या जाँचता है। सरल अंग्रेज़ी का उपयोग करो।
2. **`given:` block**: टेस्ट के लिए जो चाहिए वह सेट अप करो (objects बनाओ, डेटा तैयार करो)
3. **`expect:` block**: वास्तविक जाँच। टेस्ट पास होने के लिए प्रत्येक पंक्ति `true` होनी चाहिए

यह संरचना टेस्ट को पठनीय बनाती है: "एक extension object दिया गया है, अपेक्षा है कि `reverseGreeting('Hello')` `'olleH'` के बराबर हो।"

---

## 3. टेस्ट लिखो

भाग 3 में बनाए गए दो फ़ंक्शन के लिए टेस्ट लिखो: `reverseGreeting` और `decorateGreeting`।

### 3.1. टेस्ट class बनाओ

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

इसे अपने editor में खोलो और खाली टेस्ट class का skeleton जोड़ो:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * greeting extension फ़ंक्शन के लिए टेस्ट
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. सभी Spock टेस्ट classes `Specification` को extend करती हैं। यह किसी भी Spock टेस्ट फ़ाइल का शुरुआती बिंदु है।

### 3.2. reverseGreeting का परीक्षण करो

class body के अंदर एक टेस्ट method जोड़ो।
`given:` block एक `GreetingExtension` instance बनाता है, और `expect:` block जाँचता है कि `reverseGreeting` दो अलग-अलग इनपुट को सही तरीके से उलटता है।
यह पाइपलाइन चलाए बिना, सीधे फ़ंक्शन का परीक्षण करता है।

=== "बाद में"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension फ़ंक्शन के लिए टेस्ट
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()            // (1)!

            expect:
            ext.reverseGreeting('Hello') == 'olleH'     // (2)!
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

    1. पाइपलाइन चलाए बिना, सीधे परीक्षण करने के लिए अपने extension का एक instance बनाओ
    2. `expect:` में प्रत्येक पंक्ति एक assertion है; टेस्ट तभी पास होता है जब सभी `true` हों

=== "पहले"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension फ़ंक्शन के लिए टेस्ट
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. decorateGreeting का परीक्षण करो

पहले के बाद एक दूसरी टेस्ट method जोड़ो।
यह जाँचता है कि `decorateGreeting` इनपुट string को दोनों तरफ `***` से लपेटता है।

=== "बाद में"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension फ़ंक्शन के लिए टेस्ट
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }

        def 'should decorate a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.decorateGreeting('Hello') == '*** Hello ***'
        }
    }
    ```

=== "पहले"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension फ़ंक्शन के लिए टेस्ट
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

---

## 4. टेस्ट चलाओ

```bash
make test
```

??? example "टेस्ट आउटपुट"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **टेस्ट के परिणाम कहाँ हैं?** जब सभी टेस्ट पास हो जाते हैं तो Gradle विस्तृत आउटपुट छुपा देता है।
    "BUILD SUCCESSFUL" का मतलब है कि सब कुछ काम किया।
    अगर कोई टेस्ट विफल होता है, तो तुम्हें विस्तृत error messages दिखेंगे।

??? exercise "एक edge case टेस्ट जोड़ो"

    एक टेस्ट जोड़ो जो जाँचे कि `reverseGreeting` एक खाली string को कैसे संभालता है।
    `reverseGreeting('')` क्या लौटाना चाहिए?
    टेस्ट जोड़ो, `make test` चलाओ, और सत्यापित करो कि यह पास होता है।

    ??? solution "समाधान"

        `GreetingExtensionTest.groovy` में यह टेस्ट method जोड़ो:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        एक खाली string को उलटने पर भी खाली string ही मिलती है।

---

## 5. टेस्ट रिपोर्ट देखो

Gradle प्रत्येक टेस्ट के विस्तृत परिणामों के साथ एक HTML टेस्ट रिपोर्ट बनाता है।
रिपोर्ट डायरेक्टरी में एक web server शुरू करो:

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code तुम्हें browser में application खोलने का संकेत देगा।
अपनी टेस्ट class तक क्लिक करके जाओ और व्यक्तिगत टेस्ट परिणाम देखो:

![टेस्ट रिपोर्ट जो सभी टेस्ट पास दिखा रही है](./img/test_report.png)

रिपोर्ट प्रत्येक टेस्ट method और उसके पास या फेल होने की स्थिति दिखाती है।

Server बंद करने के लिए ++ctrl+c++ दबाओ, फिर पिछली डायरेक्टरी में वापस जाओ:

```bash
popd
```

मुख्य प्रोजेक्ट डायरेक्टरी में वापस जाओ:

```bash
cd ..
```

---

## सारांश

तुमने सीखा कि:

- Spock टेस्ट एक पठनीय `given:`/`expect:` संरचना का उपयोग करते हैं
- टेस्ट चलाने के लिए `make test` और HTML रिपोर्ट के लिए `build/reports/tests/test/` का उपयोग करो
- टेस्ट व्यवहार की जाँच करते हैं और दस्तावेज़ीकरण के रूप में काम करते हैं जो दिखाते हैं कि फ़ंक्शन कैसे उपयोग किए जाने चाहिए

---

## आगे क्या है?

अब तक, तुम्हारा plugin कस्टम फ़ंक्शन जोड़ता है जिन्हें पाइपलाइन कॉल कर सकती हैं।
Plugins trace observers का उपयोग करके workflow events (एक कार्य पूरा होना, एक फ़ाइल publish होना, पाइपलाइन समाप्त होना) पर भी प्रतिक्रिया दे सकते हैं।
अगले सेक्शन में, तुम एक observer बनाओगे जो पूर्ण हुए कार्यों की गिनती करता है और पाइपलाइन समाप्त होने पर एक सारांश प्रिंट करता है।

[भाग 5 पर जारी रखो :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
