# भाग 3: रन कॉन्फ़िगरेशन

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

यह सेक्शन यह एक्सप्लोर करेगा कि Nextflow pipeline की कॉन्फ़िगरेशन को कैसे मैनेज करें ताकि इसके व्यवहार को कस्टमाइज़ किया जा सके, इसे विभिन्न वातावरणों के अनुकूल बनाया जा सके, और रिसोर्स उपयोग को ऑप्टिमाइज़ किया जा सके _बिना workflow कोड की एक भी लाइन बदले_।

ऐसा करने के कई तरीके हैं, जिन्हें संयोजन में उपयोग किया जा सकता है और [Configuration](https://nextflow.io/docs/latest/config.html) डॉक्यूमेंटेशन में वर्णित प्राथमिकता के क्रम के अनुसार व्याख्या की जाती है।

कोर्स के इस भाग में, हम तुम्हें सबसे सरल और सबसे सामान्य कॉन्फ़िगरेशन फ़ाइल मैकेनिज़्म दिखाने जा रहे हैं, `nextflow.config` फ़ाइल, जिसका तुम पहले ही भाग 2 में containers के सेक्शन में सामना कर चुके हो।

हम Nextflow कॉन्फ़िगरेशन के आवश्यक घटकों जैसे process directives, executors, profiles, और parameter फ़ाइलों पर चर्चा करेंगे।
इन कॉन्फ़िगरेशन विकल्पों का प्रभावी ढंग से उपयोग करना सीखकर, तुम Nextflow pipelines की लचीलापन, स्केलेबिलिटी और परफ़ॉर्मेंस का पूरा लाभ उठा सकते हो।

कॉन्फ़िगरेशन के इन तत्वों का अभ्यास करने के लिए, हम उस workflow की एक नई कॉपी चलाने जा रहे हैं जिसे हमने इस ट्रेनिंग कोर्स के भाग 2 के अंत में चलाया था, जिसका नाम बदलकर `3-main.nf` कर दिया गया है।

यदि तुम Hello pipeline से परिचित नहीं हो या तुम्हें याद दिलाने की ज़रूरत है, तो [यह जानकारी पेज](../info/hello_pipeline.md) देखो।

---

## 1. Workflow इनपुट पैरामीटर्स को मैनेज करना

??? example "परिदृश्य"

    तुमने एक pipeline डाउनलोड की है और इसे समान इनपुट फ़ाइलों और सेटिंग्स के साथ बार-बार चलाना चाहते हो, लेकिन तुम हर बार सभी पैरामीटर्स टाइप नहीं करना चाहते। या शायद तुम किसी सहकर्मी के लिए pipeline सेट अप कर रहे हो जो command-line arguments के साथ सहज नहीं है।

हम कॉन्फ़िगरेशन के एक पहलू से शुरुआत करने जा रहे हैं जो केवल उसका विस्तार है जिसके साथ हम अब तक काम कर रहे हैं: इनपुट पैरामीटर्स का प्रबंधन।

वर्तमान में, हमारी workflow command-line के माध्यम से कई पैरामीटर वैल्यूज़ स्वीकार करने के लिए सेट अप है, जो workflow स्क्रिप्ट में ही एक `params` ब्लॉक में घोषित हैं।
एक में इसकी घोषणा के हिस्से के रूप में एक डिफ़ॉल्ट वैल्यू सेट है।

हालांकि, तुम उन सभी के लिए डिफ़ॉल्ट सेट करना चाह सकते हो, या मौजूदा डिफ़ॉल्ट को ओवरराइड करना चाह सकते हो बिना command line पर पैरामीटर्स निर्दिष्ट किए, या मूल स्क्रिप्ट फ़ाइल को संशोधित किए।

ऐसा करने के कई तरीके हैं; हम तुम्हें तीन बुनियादी तरीके दिखाने जा रहे हैं जो बहुत आम तौर पर उपयोग किए जाते हैं।

### 1.1. `nextflow.config` में वैल्यूज़ सेट अप करना

यह सबसे सरल दृष्टिकोण है, हालांकि यह संभवतः सबसे कम लचीला है क्योंकि मुख्य `nextflow.config` फ़ाइल कुछ ऐसी नहीं है जिसे तुम हर रन के लिए एडिट करना चाहोगे।
लेकिन इसमें workflow में पैरामीटर्स को _घोषित_ करने (जो निश्चित रूप से वहां होना चाहिए) बनाम _डिफ़ॉल्ट वैल्यूज़_ प्रदान करने की चिंताओं को अलग करने का फायदा है, जो कॉन्फ़िगरेशन फ़ाइल में अधिक उपयुक्त हैं।

चलो इसे दो चरणों में करते हैं।

#### 1.1.1. कॉन्फ़िगरेशन फ़ाइल में एक `params` ब्लॉक बनाना

`nextflow.config` फ़ाइल में निम्नलिखित कोड परिवर्तन करो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

ध्यान दो कि हमने workflow से `params` ब्लॉक को कॉन्फ़िगरेशन फ़ाइल में केवल कॉपी नहीं किया।
`batch` पैरामीटर के लिए जिसमें पहले से ही एक डिफ़ॉल्ट वैल्यू घोषित थी, सिंटैक्स थोड़ा अलग है।
workflow फ़ाइल में, वह एक typed declaration है।
कॉन्फ़िगरेशन में, वे वैल्यू असाइनमेंट हैं।

तकनीकी रूप से, यह workflow फ़ाइल में अभी भी निर्दिष्ट डिफ़ॉल्ट वैल्यूज़ को ओवरराइड करने के लिए पर्याप्त है।
तुम `batch` के लिए डिफ़ॉल्ट वैल्यू को संशोधित कर सकते हो और workflow चला सकते हो यह संतुष्ट करने के लिए कि कॉन्फ़िगरेशन फ़ाइल में सेट की गई वैल्यू workflow फ़ाइल में सेट की गई वैल्यू को ओवरराइड करती है।

लेकिन कॉन्फ़िगरेशन को पूरी तरह से कॉन्फ़िगरेशन फ़ाइल में ले जाने की भावना में, चलो workflow फ़ाइल से उस डिफ़ॉल्ट वैल्यू को पूरी तरह से हटा देते हैं।

#### 1.1.2. workflow फ़ाइल में `batch` के लिए डिफ़ॉल्ट वैल्यू हटाना

`3-main.nf` workflow फ़ाइल में निम्नलिखित कोड परिवर्तन करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

अब workflow फ़ाइल स्वयं इन पैरामीटर्स के लिए कोई डिफ़ॉल्ट वैल्यूज़ सेट नहीं करती है।

#### 1.1.3. Pipeline चलाना

चलो टेस्ट करते हैं कि यह command line में कोई पैरामीटर निर्दिष्ट किए बिना सही तरीके से काम करता है।

```bash
nextflow run 3-main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करता है।

अंतिम ASCII art आउटपुट `results/3-main/` डायरेक्टरी में है, `cowpy-COLLECTED-batch-output.txt` नाम के तहत, पहले जैसा ही।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

कार्यात्मक रूप से, इस कदम ने कुछ भी नहीं बदला है, लेकिन वैचारिक रूप से कॉन्फ़िगरेशन फ़ाइल में डिफ़ॉल्ट वैल्यूज़ सेट करना थोड़ा साफ है।

### 1.2. रन-विशिष्ट कॉन्फ़िगरेशन फ़ाइल का उपयोग करना

??? example "परिदृश्य"

    तुम अपनी मुख्य कॉन्फ़िगरेशन फ़ाइल को संशोधित किए बिना विभिन्न सेटिंग्स के साथ प्रयोग करना चाहते हो।

तुम एक सबडायरेक्टरी में एक नई `nextflow.config` फ़ाइल बनाकर ऐसा कर सकते हो जिसे तुम अपने प्रयोगों के लिए working directory के रूप में उपयोग करोगे।

#### 1.2.1. खाली कॉन्फ़िगरेशन के साथ working directory बनाना

चलो एक नई डायरेक्टरी बनाकर और उसमें जाकर शुरू करते हैं:

```bash
mkdir -p tux-run
cd tux-run
```

फिर, उस डायरेक्टरी में एक खाली कॉन्फ़िगरेशन फ़ाइल बनाओ:

```bash
touch nextflow.config
```

यह एक खाली फ़ाइल उत्पन्न करता है।

#### 1.2.2. प्रयोगात्मक कॉन्फ़िगरेशन सेट अप करना

अब नई फ़ाइल खोलो और उन पैरामीटर्स को जोड़ो जिन्हें तुम कस्टमाइज़ करना चाहते हो:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

ध्यान दो कि इनपुट फ़ाइल का पाथ डायरेक्टरी संरचना को प्रतिबिंबित करना चाहिए।

#### 1.2.3. Pipeline चलाना

अब हम अपनी नई working directory के भीतर से अपनी pipeline चला सकते हैं।
पाथ को तदनुसार अनुकूलित करना सुनिश्चित करो!

```bash
nextflow run ../3-main.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

यह `tux-run/` के तहत डायरेक्टरियों का एक नया सेट बनाएगा जिसमें `tux-run/work/` और `tux-run/results/` शामिल हैं।

इस रन में, Nextflow हमारी वर्तमान डायरेक्टरी में `nextflow.config` को pipeline की रूट डायरेक्टरी में `nextflow.config` के साथ जोड़ता है, और इस प्रकार डिफ़ॉल्ट character (turkey) को tux character के साथ ओवरराइड करता है।

अंतिम आउटपुट फ़ाइल में tux character को greetings कहते हुए होना चाहिए।

??? abstract "फ़ाइल सामग्री"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/

    ```

बस इतना ही; अब तुम्हारे पास अपनी 'सामान्य' कॉन्फ़िगरेशन को संशोधित किए बिना प्रयोग करने के लिए एक स्थान है।

!!! warning

    अगले सेक्शन में जाने से पहले पिछली डायरेक्टरी में वापस जाना सुनिश्चित करो!

    ```bash
    cd ..
    ```

अब चलो पैरामीटर वैल्यूज़ सेट करने के एक और उपयोगी तरीके को देखते हैं।

### 1.3. पैरामीटर फ़ाइल का उपयोग करना

??? example "परिदृश्य"

    तुम्हें किसी सहयोगी के साथ सटीक रन पैरामीटर्स साझा करने की ज़रूरत है, या उन्हें किसी प्रकाशन के लिए रिकॉर्ड करने की ज़रूरत है।

सबडायरेक्टरी दृष्टिकोण प्रयोग के लिए बहुत अच्छा काम करता है, लेकिन इसमें थोड़ा सेटअप शामिल है और आवश्यक है कि तुम पाथ को तदनुसार अनुकूलित करो।
जब तुम अपनी pipeline को वैल्यूज़ के एक विशिष्ट सेट के साथ चलाना चाहते हो, या किसी और को न्यूनतम प्रयास के साथ ऐसा करने में सक्षम बनाना चाहते हो, तो एक सरल दृष्टिकोण है।

Nextflow हमें YAML या JSON फ़ॉर्मेट में [parameter file](https://nextflow.io/docs/latest/config.html#parameter-file) के माध्यम से पैरामीटर्स निर्दिष्ट करने की अनुमति देता है, जो डिफ़ॉल्ट वैल्यूज़ के वैकल्पिक सेट को मैनेज और वितरित करना बहुत सुविधाजनक बनाता है, उदाहरण के लिए, साथ ही रन-विशिष्ट पैरामीटर वैल्यूज़।

#### 1.3.1. उदाहरण पैरामीटर फ़ाइल की जांच करना

इसे प्रदर्शित करने के लिए, हम वर्तमान डायरेक्टरी में एक उदाहरण पैरामीटर फ़ाइल प्रदान करते हैं, जिसे `test-params.yaml` कहा जाता है:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

इस पैरामीटर फ़ाइल में प्रत्येक इनपुट के लिए एक key-value pair है जिसे हम निर्दिष्ट करना चाहते हैं।
यदि तुम सिंटैक्स की तुलना कॉन्फ़िगरेशन फ़ाइल से करते हो तो colons (`:`) के उपयोग पर ध्यान दो equal signs (`=`) के बजाय।
config फ़ाइल Groovy में लिखी गई है, जबकि पैरामीटर फ़ाइल YAML में लिखी गई है।

!!! info

    हम एक उदाहरण के रूप में पैरामीटर फ़ाइल का JSON संस्करण भी प्रदान करते हैं लेकिन हम इसे यहां चलाने नहीं जा रहे हैं।
    अपने दम पर उसे आज़माने के लिए स्वतंत्र महसूस करो।

#### 1.3.2. Pipeline चलाना

इस पैरामीटर फ़ाइल के साथ workflow चलाने के लिए, बस बेस कमांड में `-params-file <filename>` जोड़ो।

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

अंतिम आउटपुट फ़ाइल में stegosaurus character को greetings कहते हुए होना चाहिए।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

पैरामीटर फ़ाइल का उपयोग करना ओवरकिल लग सकता है जब तुम्हारे पास निर्दिष्ट करने के लिए केवल कुछ पैरामीटर्स हैं, लेकिन कुछ pipelines दर्जनों पैरामीटर्स की अपेक्षा करती हैं।
उन मामलों में, पैरामीटर फ़ाइल का उपयोग करने से हमें बड़े command lines टाइप किए बिना और workflow स्क्रिप्ट को संशोधित किए बिना रनटाइम पर पैरामीटर वैल्यूज़ प्रदान करने की अनुमति मिलेगी।

यह सहयोगियों को पैरामीटर्स के सेट वितरित करना भी आसान बनाता है, या उदाहरण के लिए, किसी प्रकाशन के लिए सहायक जानकारी के रूप में।
यह तुम्हारे काम को दूसरों द्वारा अधिक पुनरुत्पादनीय बनाता है।

### सारांश

तुम workflow इनपुट्स को मैनेज करने के लिए प्रमुख कॉन्फ़िगरेशन विकल्पों का लाभ उठाना जानते हो।

### आगे क्या है?

सीखो कि अपने workflow आउटपुट्स को कहां और कैसे प्रकाशित किया जाए, इसे कैसे मैनेज करें।

---

## 2. Workflow आउटपुट्स को मैनेज करना

??? example "परिदृश्य"

    तुम्हारी pipeline आउटपुट्स को एक हार्डकोडेड डायरेक्टरी में प्रकाशित करती है, लेकिन तुम हर बार workflow कोड एडिट किए बिना प्रोजेक्ट या प्रयोग नाम के अनुसार परिणामों को व्यवस्थित करना चाहते हो।

हमें जो workflow विरासत में मिली है वह workflow-level आउटपुट घोषणाओं के लिए पाथ का उपयोग करती है, जो बहुत लचीली नहीं है और बहुत सारी पुनरावृत्ति शामिल है।

चलो कुछ सामान्य तरीकों को देखते हैं जिनसे तुम इसे अधिक लचीला बनाने के लिए कॉन्फ़िगर कर सकते हो।

### 2.1. `outputDir` डायरेक्टरी नाम को कस्टमाइज़ करना

हमने अब तक चलाई गई workflow के प्रत्येक संस्करण ने अपने आउटपुट्स को आउटपुट परिभाषाओं में हार्डकोडेड एक अलग सबडायरेक्टरी में प्रकाशित किया है।

हमने भाग 1 में `-output-dir` CLI flag का उपयोग करके बदल दिया कि वह सबडायरेक्टरी कहां थी, लेकिन वह अभी भी केवल एक स्टैटिक स्ट्रिंग है।
चलो इसके बजाय इसे एक config फ़ाइल में कॉन्फ़िगर करते हैं, जहां हम अधिक जटिल डायनामिक पाथ परिभाषित कर सकते हैं।
हम इसके लिए एक पूरा नया पैरामीटर बना सकते थे, लेकिन चलो `batch` पैरामीटर का उपयोग करते हैं क्योंकि यह वहीं है।

#### 2.1.1. कॉन्फ़िगरेशन फ़ाइल में `outputDir` के लिए एक वैल्यू सेट करना

Nextflow जो पाथ आउटपुट्स प्रकाशित करने के लिए उपयोग करता है वह `outputDir` विकल्प द्वारा नियंत्रित होता है।
सभी आउटपुट्स के लिए पाथ बदलने के लिए, तुम `nextflow.config` कॉन्फ़िगरेशन फ़ाइल में इस विकल्प के लिए एक वैल्यू सेट कर सकते हो।

`nextflow.config` फ़ाइल में निम्नलिखित कोड जोड़ो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

यह बिल्ट-इन डिफ़ॉल्ट पाथ, `results/`, को `results_config/` प्लस `batch` पैरामीटर की वैल्यू को सबडायरेक्टरी के रूप में बदल देगा।

याद रखो कि तुम अपनी कमांड में `-output-dir` पैरामीटर का उपयोग करके command-line से भी इस विकल्प को सेट कर सकते हो (संक्षेप में `-o`), लेकिन फिर तुम `batch` पैरामीटर वैल्यू का उपयोग नहीं कर सकते।
CLI flag का उपयोग करने से config में `outputDir` को ओवरराइट कर दिया जाएगा यदि यह सेट है।

#### 2.1.2. हार्डकोडेड पाथ के दोहराए गए हिस्से को हटाना

हमारे पास अभी भी आउटपुट विकल्पों में एक सबडायरेक्टरी हार्डकोडेड है, तो चलो अब उससे छुटकारा पाते हैं।

workflow फ़ाइल में निम्नलिखित कोड परिवर्तन करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

हम प्रत्येक पाथ में `${params.batch}` भी जोड़ सकते थे `outputDir` डिफ़ॉल्ट को संशोधित करने के बजाय, लेकिन यह अधिक संक्षिप्त है।

#### 2.1.3. Pipeline चलाना

चलो टेस्ट करते हैं कि यह सही तरीके से काम करता है, command line से batch नाम को `outdir` पर सेट करते हुए।

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करता है, सिवाय इस बार हम अपने आउटपुट्स `results_config/outdir/` के तहत पाते हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

तुम इस दृष्टिकोण को कस्टम पाथ परिभाषाओं के साथ जोड़ सकते हो ताकि तुम्हें पसंद की कोई भी डायरेक्टरी पदानुक्रम बना सको।

### 2.2. Process द्वारा आउटपुट्स को व्यवस्थित करना

आउटपुट्स को आगे व्यवस्थित करने का एक लोकप्रिय तरीका process द्वारा करना है, _यानी_ pipeline में चलाए गए प्रत्येक process के लिए सबडायरेक्टरियां बनाना।

#### 2.2.1. आउटपुट पाथ को process नामों के संदर्भ से बदलना

तुम्हें बस आउटपुट पाथ घोषणा में `<process>.name` के रूप में process के नाम का संदर्भ देना है।

workflow फ़ाइल में निम्नलिखित परिवर्तन करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

यह आउटपुट पाथ कॉन्फ़िगरेशन से शेष हार्डकोडेड तत्वों को हटा देता है।

#### 2.2.2. Pipeline चलाना

चलो टेस्ट करते हैं कि यह सही तरीके से काम करता है, command line से batch नाम को `pnames` पर सेट करते हुए।

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करता है, सिवाय इस बार हम अपने आउटपुट्स `results_config/pnames/` के तहत पाते हैं, और वे process द्वारा समूहीकृत हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note

    ध्यान दो कि यहां हमने `intermediates` बनाम शीर्ष स्तर पर अंतिम आउटपुट्स के बीच के अंतर को मिटा दिया है।
    तुम इन दृष्टिकोणों को मिला और मिलान कर सकते हो और यहां तक कि कई वेरिएबल्स भी शामिल कर सकते हो, उदाहरण के लिए पहले आउटपुट के पाथ को `#!groovy "${params.batch}/intermediates/${sayHello.name}"` के रूप में सेट करके

### 2.3. Workflow स्तर पर publish mode सेट करना

अंत में, दोहराए जाने वाले कोड की मात्रा को कम करने की भावना में, हम प्रति-आउटपुट `mode` घोषणाओं को कॉन्फ़िगरेशन में एक एकल लाइन से बदल सकते हैं।

#### 2.3.1. कॉन्फ़िगरेशन फ़ाइल में `workflow.output.mode` जोड़ना

`nextflow.config` फ़ाइल में निम्नलिखित कोड जोड़ो:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results_config/${params.batch}"
    ```

`outputDir` विकल्प की तरह, कॉन्फ़िगरेशन फ़ाइल में `workflow.output.mode` को एक वैल्यू देना workflow फ़ाइल में सेट की गई चीज़ को ओवरराइड करने के लिए पर्याप्त होगा, लेकिन चलो वैसे भी अनावश्यक कोड को हटा देते हैं।

#### 2.3.2. Workflow फ़ाइल से आउटपुट mode हटाना

workflow फ़ाइल में निम्नलिखित परिवर्तन करो:

=== "बाद में"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "पहले"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

यह अधिक संक्षिप्त है, है ना?

#### 2.3.3. Pipeline चलाना

चलो टेस्ट करते हैं कि यह सही तरीके से काम करता है, command line से batch नाम को `outmode` पर सेट करते हुए।

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करता है, सिवाय इस बार हम अपने आउटपुट्स `results_config/outmode/` के तहत पाते हैं।
वे अभी भी सभी उचित कॉपियां हैं, symlinks नहीं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

मुख्य कारण जिससे तुम अभी भी प्रति-आउटपुट तरीके से mode सेट करना चाह सकते हो वह है यदि तुम एक ही workflow के भीतर मिश्रण और मिलान करना चाहते हो, _यानी_ कुछ आउटपुट्स को कॉपी किया जाए और कुछ को symlinked किया जाए।

बहुत सारे अन्य विकल्प हैं जिन्हें तुम इस तरह से कस्टमाइज़ कर सकते हो, लेकिन उम्मीद है कि यह तुम्हें विकल्पों की रेंज और अपनी प्राथमिकताओं के अनुरूप उन्हें प्रभावी ढंग से उपयोग करने का एहसास देता है।

### सारांश

तुम जानते हो कि उन डायरेक्टरियों के नामकरण और संरचना को कैसे नियंत्रित करें जहां तुम्हारे आउटपुट्स प्रकाशित होते हैं, साथ ही workflow आउटपुट publishing mode भी।

### आगे क्या है?

सीखो कि अपने workflow कॉन्फ़िगरेशन को अपने compute वातावरण के अनुकूल कैसे बनाएं, software packaging technology से शुरू करते हुए।

---

## 3. Software packaging technology का चयन करना

अब तक हम कॉन्फ़िगरेशन तत्वों को देख रहे हैं जो नियंत्रित करते हैं कि इनपुट्स कैसे अंदर जाते हैं और इनपुट्स कहां बाहर आते हैं। अब समय है अपने workflow कॉन्फ़िगरेशन को अपने compute वातावरण के अनुकूप बनाने पर अधिक विशेष रूप से ध्यान केंद्रित करने का।

उस पथ पर पहला कदम यह निर्दिष्ट करना है कि प्रत्येक चरण में चलाए जाने वाले software packages कहां से आने वाले हैं।
क्या वे पहले से ही local compute वातावरण में इंस्टॉल हैं?
क्या हमें images प्राप्त करने और उन्हें container system के माध्यम से चलाने की ज़रूरत है?
या क्या हमें Conda packages प्राप्त करने और एक local Conda वातावरण बनाने की ज़रूरत है?

इस ट्रेनिंग कोर्स के पहले भाग (भाग 1-4) में हमने अपनी workflow में केवल locally installed software का उपयोग किया।
फिर भाग 5 में, हमने Docker containers और `nextflow.config` फ़ाइल पेश की, जिसका उपयोग हमने Docker containers के उपयोग को सक्षम करने के लिए किया।

अब चलो देखते हैं कि हम `nextflow.config` फ़ाइल के माध्यम से एक वैकल्पिक software packaging विकल्प को कैसे कॉन्फ़िगर कर सकते हैं।

### 3.1. Config फ़ाइल में Docker को अक्षम करना और Conda को सक्षम करना

??? example "परिदृश्य"

    तुम अपनी pipeline को एक HPC cluster में ले जा रहे हो जहां सुरक्षा कारणों से Docker की अनुमति नहीं है।
    Cluster Singularity और Conda का समर्थन करता है, इसलिए तुम्हें तदनुसार अपनी कॉन्फ़िगरेशन को स्विच करने की ज़रूरत है।

जैसा कि पहले उल्लेख किया गया है, Nextflow कई container technologies का समर्थन करता है जिसमें Singularity (जो HPC पर अधिक व्यापक रूप से उपयोग किया जाता है) शामिल है, साथ ही software package managers जैसे Conda भी।

हम अपनी कॉन्फ़िगरेशन फ़ाइल को Docker के बजाय Conda का उपयोग करने के लिए बदल सकते हैं।
ऐसा करने के लिए, चलो `docker.enabled` की वैल्यू को `false` में स्विच करते हैं, और Conda के उपयोग को सक्षम करने वाला एक directive जोड़ते हैं:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

यह Nextflow को उन processes के लिए Conda वातावरण बनाने और उपयोग करने की अनुमति देगा जिनमें Conda packages निर्दिष्ट हैं।
जिसका मतलब है कि अब हमें अपनी `cowpy` process में उनमें से एक जोड़ने की ज़रूरत है!

### 3.2. Process परिभाषा में एक Conda package निर्दिष्ट करना

हमने पहले ही `cowpy` tool वाले Conda package के लिए URI प्राप्त कर लिया है: `conda-forge::cowpy==1.1.5`

अब हम `conda` directive का उपयोग करके `cowpy` process परिभाषा में URI जोड़ते हैं:

=== "बाद में"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "पहले"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

स्पष्ट होने के लिए, हम `docker` directive को _बदल_ नहीं रहे हैं, हम एक वैकल्पिक विकल्प _जोड़_ रहे हैं।

!!! tip

    किसी दिए गए conda package के लिए URI प्राप्त करने के कुछ अलग तरीके हैं।
    हम [Seqera Containers](https://seqera.io/containers/) search query का उपयोग करने की सलाह देते हैं, जो तुम्हें एक URI देगा जिसे तुम कॉपी और पेस्ट कर सकते हो, भले ही तुम इससे container बनाने की योजना नहीं बना रहे हो।

### 3.3. Workflow चलाकर सत्यापित करना कि यह Conda का उपयोग कर सकती है

चलो इसे आज़माते हैं।

```bash
nextflow run 3-main.nf --batch conda
```

??? success "कमांड आउटपुट"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

यह बिना किसी समस्या के काम करना चाहिए और `results_config/conda` के तहत पहले जैसे ही आउटपुट्स उत्पन्न करना चाहिए।

पर्दे के पीछे, Nextflow ने Conda packages प्राप्त किए हैं और वातावरण बनाया है, जो सामान्य रूप से थोड़ा काम लेता है; तो यह अच्छा है कि हमें वह सब कुछ खुद नहीं करना पड़ता!

!!! info

    यह जल्दी चलता है क्योंकि `cowpy` package काफी छोटा है, लेकिन यदि तुम बड़े packages के साथ काम कर रहे हो, तो पहली बार सामान्य से थोड़ा अधिक समय लग सकता है, और तुम console output को पूरा होने से पहले एक मिनट या उससे अधिक समय तक 'अटका हुआ' देख सकते हो।
    यह सामान्य है और पहली बार जब तुम एक नए package का उपयोग करते हो तो Nextflow द्वारा किए गए अतिरिक्त काम के कारण है।

हमारे दृष्टिकोण से, ऐसा लगता है कि यह Docker के साथ चलने के समान ही काम करता है, भले ही backend पर mechanics थोड़ा अलग हैं।

इसका मतलब है कि हम ज़रूरत पड़ने पर Conda वातावरण के साथ चलाने के लिए तैयार हैं।

??? info "Docker और Conda को मिलाना और मिलान करना"

    चूंकि ये directives प्रति process असाइन किए जाते हैं, इसलिए 'मिश्रण और मिलान' करना संभव है, _यानी_ अपनी workflow में कुछ processes को Docker के साथ चलाने के लिए कॉन्फ़िगर करना और अन्य को Conda के साथ, उदाहरण के लिए, यदि तुम जो compute infrastructure उपयोग कर रहे हो वह दोनों का समर्थन करता है।
    उस स्थिति में, तुम अपनी कॉन्फ़िगरेशन फ़ाइल में Docker और Conda दोनों को सक्षम करोगे।
    यदि किसी दिए गए process के लिए दोनों उपलब्ध हैं, तो Nextflow containers को प्राथमिकता देगा।

    और जैसा कि पहले उल्लेख किया गया है, Nextflow कई अन्य software packaging और container technologies का समर्थन करता है, इसलिए तुम केवल उन दोनों तक सीमित नहीं हो।

### सारांश

तुम जानते हो कि प्रत्येक process को कौन सा software package उपयोग करना चाहिए, और technologies के बीच कैसे स्विच करना है, यह कैसे कॉन्फ़िगर करें।

### आगे क्या है?

सीखो कि Nextflow द्वारा वास्तव में काम करने के लिए उपयोग किए जाने वाले execution platform को कैसे बदलें।

---

## 4. Execution platform का चयन करना

??? example "परिदृश्य"

    तुम अपने laptop पर अपनी pipeline विकसित और परीक्षण कर रहे हो, लेकिन अब तुम्हें इसे हजारों samples पर चलाने की ज़रूरत है।
    तुम्हारी संस्था के पास Slurm scheduler के साथ एक HPC cluster है जिसे तुम इसके बजाय उपयोग करना चाहते हो।

अब तक, हम अपनी pipeline को local executor के साथ चला रहे हैं।
यह प्रत्येक task को उस machine पर execute करता है जिस पर Nextflow चल रहा है।
जब Nextflow शुरू होता है, तो यह उपलब्ध CPUs और memory को देखता है।
यदि चलाने के लिए तैयार tasks के resources उपलब्ध resources से अधिक हैं, तो Nextflow अंतिम tasks को execution से रोक देगा जब तक कि पहले के एक या अधिक tasks समाप्त नहीं हो जाते, आवश्यक resources को मुक्त करते हुए।

Local executor सुविधाजनक और कुशल है, लेकिन यह उस एकल machine तक सीमित है। बहुत बड़े workloads के लिए, तुम पा सकते हो कि तुम्हारी local machine एक bottleneck है, या तो क्योंकि तुम्हारे पास एक एकल task है जिसे तुम्हारे पास उपलब्ध से अधिक resources की आवश्यकता है, या क्योंकि तुम्हारे पास इतने सारे tasks हैं कि एक एकल machine के उन्हें चलाने की प्रतीक्षा करने में बहुत अधिक समय लगेगा।

Nextflow [कई अलग-अलग execution backends](https://nextflow.io/docs/latest/executor.html) का समर्थन करता है, जिसमें HPC schedulers (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor और अन्य) के साथ-साथ cloud execution backends जैसे (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes और अधिक) शामिल हैं।

### 4.1. एक अलग backend को लक्षित करना

Executor का चुनाव `executor` नामक process directive द्वारा सेट किया जाता है।
डिफ़ॉल्ट रूप से यह `local` पर सेट है, इसलिए निम्नलिखित कॉन्फ़िगरेशन निहित है:

```groovy title="बिल्ट-इन कॉन्फ़िगरेशन"
process {
    executor = 'local'
}
```

Executor को एक अलग backend को लक्षित करने के लिए सेट करने के लिए, तुम बस उस executor को निर्दिष्ट करोगे जो तुम चाहते हो resource allocations के लिए ऊपर वर्णित समान सिंटैक्स का उपयोग करके (सभी विकल्पों के लिए [Executors](https://nextflow.io/docs/latest/executor.html) देखें)।

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    हम वास्तव में इसे ट्रेनिंग वातावरण में परीक्षण नहीं कर सकते क्योंकि यह HPC से कनेक्ट करने के लिए सेट अप नहीं है।

### 4.2. Execution पैरामीटर्स के लिए backend-विशिष्ट सिंटैक्स से निपटना

अधिकांश high-performance computing platforms आपको कुछ पैरामीटर्स निर्दिष्ट करने की अनुमति देते हैं (और कभी-कभी आवश्यकता होती है) जैसे resource allocation requests और limitations (उदाहरण के लिए CPUs और memory की संख्या) और उपयोग करने के लिए job queue का नाम।

दुर्भाग्य से, इनमें से प्रत्येक system विभिन्न technologies, syntaxes और configurations का उपयोग करता है यह परिभाषित करने के लिए कि एक job को कैसे परिभाषित और relevant scheduler को submit किया जाना चाहिए।

??? abstract "उदाहरण"

    उदाहरण के लिए, 8 CPUs और 4GB RAM की आवश्यकता वाली एक ही job को queue "my-science-work" पर execute करने के लिए backend के आधार पर निम्नलिखित विभिन्न तरीकों से व्यक्त करने की आवश्यकता है।

    ```bash title="SLURM के लिए Config / sbatch का उपयोग करके submit करें"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="PBS के लिए Config / qsub का उपयोग करके submit करें"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="SGE के लिए Config / qsub का उपयोग करके submit करें"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

सौभाग्य से, Nextflow यह सब सरल बनाता है।
यह एक मानकीकृत सिंटैक्स प्रदान करता है ताकि तुम `cpus`, `memory` और `queue` जैसी प्रासंगिक properties को केवल एक बार निर्दिष्ट कर सको (सभी उपलब्ध विकल्पों के लिए [Process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) देखें)।
फिर, runtime पर, Nextflow executor setting के आधार पर उपयुक्त backend-विशिष्ट scripts उत्पन्न करने के लिए उन settings का उपयोग करेगा।

हम अगले सेक्शन में उस मानकीकृत सिंटैक्स को कवर करेंगे।

### सारांश

तुम अब जानते हो कि विभिन्न प्रकार के computing infrastructure का उपयोग करने के लिए executor को कैसे बदलें।

### आगे क्या है?

सीखो कि Nextflow में resource allocations और limitations का मूल्यांकन और व्यक्त कैसे करें।

---

## 5. Compute resource allocations को नियंत्रित करना

??? example "परिदृश्य"

    तुम्हारी pipeline cluster पर विफल होती रहती है क्योंकि tasks memory limits से अधिक होने के लिए मारे जा रहे हैं।
    या शायद तुम्हें उन resources के लिए शुल्क लिया जा रहा है जिनका तुम उपयोग नहीं कर रहे हो और लागत को अनुकूलित करना चाहते हो।

अधिकांश high-performance computing platforms आपको कुछ resource allocation पैरामीटर्स निर्दिष्ट करने की अनुमति देते हैं (और कभी-कभी आवश्यकता होती है) जैसे CPUs और memory की संख्या।

डिफ़ॉल्ट रूप से, Nextflow प्रत्येक process के लिए एक एकल CPU और 2GB memory का उपयोग करेगा।
संबंधित process directives को `cpus` और `memory` कहा जाता है, इसलिए निम्नलिखित कॉन्फ़िगरेशन निहित है:

```groovy title="बिल्ट-इन कॉन्फ़िगरेशन" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

तुम इन वैल्यूज़ को संशोधित कर सकते हो, या तो सभी processes के लिए या विशिष्ट नामित processes के लिए, अपनी कॉन्फ़िगरेशन फ़ाइल में अतिरिक्त process directives का उपयोग करके।
Nextflow उन्हें चुने गए executor के लिए उपयुक्त निर्देशों में अनुवादित करेगा।

लेकिन तुम कैसे जानते हो कि किन वैल्यूज़ का उपयोग करना है?

### 5.1. Resource utilization report उत्पन्न करने के लिए workflow चलाना

??? example "परिदृश्य"

    तुम नहीं जानते कि तुम्हारी processes को कितनी memory या CPU की ज़रूरत है और resources बर्बाद करने या jobs को मारे जाने से बचना चाहते हो।

यदि तुम पहले से नहीं जानते कि तुम्हारी processes को कितनी CPU और memory की आवश्यकता होने की संभावना है, तो तुम कुछ resource profiling कर सकते हो, जिसका अर्थ है कि तुम कुछ डिफ़ॉल्ट allocations के साथ workflow चलाते हो, रिकॉर्ड करते हो कि प्रत्येक process ने कितना उपयोग किया, और वहां से, बेस allocations को कैसे समायोजित करना है इसका अनुमान लगाते हो।

सुविधाजनक रूप से, Nextflow में इसे करने के लिए बिल्ट-इन tools शामिल हैं, और अनुरोध पर तुम्हारे लिए खुशी से एक report उत्पन्न करेगा।

ऐसा करने के लिए, अपनी command line में `-with-report <filename>.html` जोड़ो।

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Report एक html फ़ाइल है, जिसे तुम डाउनलोड कर सकते हो और अपने browser में खोल सकते हो। तुम इसे बाईं ओर file explorer में right click भी कर सकते हो और ट्रेनिंग वातावरण में इसे देखने के लिए `Show preview` पर click कर सकते हो।

Report को देखने में कुछ मिनट लगाओ और देखो कि क्या तुम resources को समायोजित करने के लिए कुछ अवसरों की पहचान कर सकते हो।
उन tabs पर click करना सुनिश्चित करो जो utilization परिणामों को आवंटित किए गए प्रतिशत के रूप में दिखाते हैं।

सभी उपलब्ध सुविधाओं पर डॉक्यूमेंटेशन के लिए [Reports](https://nextflow.io/docs/latest/reports.html) देखें।

### 5.2. सभी processes के लिए resource allocations सेट करना

Profiling दिखाती है कि हमारी ट्रेनिंग workflow में processes बहुत हल्की हैं, तो चलो प्रति process डिफ़ॉल्ट memory allocation को 1GB तक कम करते हैं।

अपनी `nextflow.config` फ़ाइल में निम्नलिखित जोड़ो, pipeline parameters सेक्शन से पहले:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

यह हमारे द्वारा उपभोग किए जाने वाले compute की मात्रा को कम करने में मदद करेगा।

### 5.3. एक विशिष्ट process के लिए resource allocations सेट करना

साथ ही, हम यह दिखावा करने जा रहे हैं कि `cowpy` process को अन्य की तुलना में अधिक resources की आवश्यकता है, बस इसलिए कि हम प्रदर्शित कर सकें कि एक व्यक्तिगत process के लिए allocations को कैसे समायोजित करें।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

इस कॉन्फ़िगरेशन के साथ, सभी processes 1GB memory और एक एकल CPU (निहित डिफ़ॉल्ट) का अनुरोध करेंगी, सिवाय `cowpy` process के, जो 2GB और 2 CPUs का अनुरोध करेगी।

!!! info

    यदि तुम्हारे पास कम CPUs वाली machine है और तुम प्रति process एक उच्च संख्या आवंटित करते हो, तो तुम process calls को एक दूसरे के पीछे queued होते हुए देख सकते हो।
    ऐसा इसलिए है क्योंकि Nextflow सुनिश्चित करता है कि हम उपलब्ध से अधिक CPUs का अनुरोध नहीं करते।

### 5.4. अपडेट की गई कॉन्फ़िगरेशन के साथ workflow चलाना

चलो इसे आज़माते हैं, profiling report के लिए एक अलग filename प्रदान करते हुए ताकि हम कॉन्फ़िगरेशन परिवर्तनों से पहले और बाद में performance की तुलना कर सकें।

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

तुम शायद कोई वास्तविक अंतर नहीं देखोगे क्योंकि यह इतना छोटा workload है, लेकिन यह वह दृष्टिकोण है जिसका उपयोग तुम वास्तविक दुनिया की workflow के performance और resource requirements का विश्लेषण करने के लिए करोगे।

यह बहुत उपयोगी है जब तुम्हारी processes की अलग-अलग resource requirements हैं। यह तुम्हें वास्तविक डेटा के आधार पर प्रत्येक process के लिए सेट किए गए resource allocations को सही आकार देने का अधिकार देता है, अनुमान नहीं।

!!! tip

    यह केवल एक छोटा सा स्वाद है जो तुम resources के अपने उपयोग को अनुकूलित करने के लिए कर सकते हो।
    Nextflow में स्वयं कुछ वास्तव में साफ [dynamic retry logic](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) बिल्ट-इन है जो resource limitations के कारण विफल होने वाली jobs को retry करने के लिए है।
    इसके अतिरिक्त, Seqera Platform तुम्हारे resource allocations को स्वचालित रूप से अनुकूलित करने के लिए AI-driven tooling भी प्रदान करता है।

### 5.5. Resource limits जोड़ना

तुम जो computing executor और compute infrastructure उपयोग कर रहे हो उसके आधार पर, तुम क्या आवंटित कर सकते हो (या आवंटित करना चाहिए) पर कुछ बाधाएं हो सकती हैं।
उदाहरण के लिए, तुम्हारा cluster तुम्हें कुछ सीमाओं के भीतर रहने की आवश्यकता कर सकता है।

तुम प्रासंगिक limitations सेट करने के लिए `resourceLimits` directive का उपयोग कर सकते हो। जब यह process block में अकेला होता है तो सिंटैक्स इस तरह दिखता है:

```groovy title="सिंटैक्स उदाहरण"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow इन वैल्यूज़ को तुम्हारे द्वारा निर्दिष्ट executor के आधार पर उपयुक्त निर्देशों में अनुवादित करेगा।

हम इसे नहीं चलाने जा रहे हैं, क्योंकि हमारे पास ट्रेनिंग वातावरण में प्रासंगिक infrastructure तक पहुंच नहीं है।
हालांकि, यदि तुम इन सीमाओं से अधिक resource allocations के साथ workflow चलाने की कोशिश करते हो, तो `.command.run` स्क्रिप्ट फ़ाइल में `sbatch` कमांड को देखो, तुम देखोगे कि वास्तव में executor को भेजे जाने वाले अनुरोध `resourceLimits` द्वारा निर्दिष्ट वैल्यूज़ पर capped हैं।

??? info "संस्थागत संदर्भ कॉन्फ़िगरेशन"

    nf-core प्रोजेक्ट ने दुनिया भर की विभिन्न संस्थाओं द्वारा साझा की गई [कॉन्फ़िगरेशन फ़ाइलों का संग्रह](https://nf-co.re/configs/) संकलित किया है, जो HPC और cloud executors की एक विस्तृत श्रृंखला को कवर करता है।

    वे साझा configs उन लोगों के लिए मूल्यवान हैं जो वहां काम करते हैं और इसलिए बस अपनी संस्था की कॉन्फ़िगरेशन को out of the box उपयोग कर सकते हैं, और उन लोगों के लिए एक मॉडल के रूप में जो अपने स्वयं के infrastructure के लिए कॉन्फ़िगरेशन विकसित करना चाहते हैं।

### सारांश

तुम जानते हो कि resource utilization का आकलन करने के लिए profiling report कैसे उत्पन्न करें और सभी processes और/या व्यक्तिगत processes के लिए resource allocations को कैसे संशोधित करें, साथ ही HPC पर चलाने के लिए resource limitations कैसे सेट करें।

### आगे क्या है?

सीखो कि preset कॉन्फ़िगरेशन profiles कैसे सेट अप करें और runtime पर उनके बीच कैसे स्विच करें।

---

## 6. Preset कॉन्फ़िगरेशन के बीच स्विच करने के लिए profiles का उपयोग करना

??? example "परिदृश्य"

    तुम नियमित रूप से development के लिए अपने laptop पर और production runs के लिए अपनी संस्था के HPC पर pipelines चलाने के बीच स्विच करते हो।
    तुम हर बार वातावरण बदलने पर manually कॉन्फ़िगरेशन सेटिंग्स बदलने से थक गए हो।

हमने तुम्हें कई तरीके दिखाए हैं जिनसे तुम अपनी pipeline कॉन्फ़िगरेशन को कस्टमाइज़ कर सकते हो जो तुम जिस प्रोजेक्ट पर काम कर रहे हो या तुम जो compute वातावरण उपयोग कर रहे हो उसके आधार पर।

तुम वैकल्पिक सेटिंग्स के बीच स्विच करना चाह सकते हो जो तुम जो computing infrastructure उपयोग कर रहे हो उसके आधार पर। उदाहरण के लिए, तुम अपने laptop पर locally छोटे पैमाने के परीक्षण विकसित और चलाना चाह सकते हो, फिर HPC या cloud पर पूर्ण पैमाने के workloads चला सकते हो।

Nextflow तुम्हें किसी भी संख्या में [**profiles**](https://nextflow.io/docs/latest/config.html#profiles) सेट अप करने देता है जो विभिन्न कॉन्फ़िगरेशन का वर्णन करते हैं, जिन्हें तुम फिर command-line argument का उपयोग करके runtime पर चुन सकते हो, बजाय कॉन्फ़िगरेशन फ़ाइल को स्वयं संशोधित करने के।

### 6.1. Local development और HPC पर execution के बीच स्विच करने के लिए profiles बनाना

चलो दो वैकल्पिक profiles सेट अप करते हैं; एक नियमित computer पर छोटे पैमाने के loads चलाने के लिए, जहां हम Docker containers का उपयोग करेंगे, और एक Slurm scheduler के साथ university HPC पर चलाने के लिए, जहां हम Conda packages का उपयोग करेंगे।

#### 6.1.1. Profiles सेट अप करना

अपनी `nextflow.config` फ़ाइल में निम्नलिखित जोड़ो, pipeline parameters सेक्शन के बाद लेकिन output settings से पहले:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

तुम देखते हो कि university HPC के लिए, हम resource limitations भी निर्दिष्ट कर रहे हैं।

#### 6.1.2. एक profile के साथ workflow चलाना

अपनी Nextflow command line में एक profile निर्दिष्ट करने के लिए, हम `-profile` argument का उपयोग करते हैं।

चलो `my_laptop` कॉन्फ़िगरेशन के साथ workflow चलाने की कोशिश करते हैं।

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

जैसा कि तुम देख सकते हो, यह हमें runtime पर कॉन्फ़िगरेशन के बीच बहुत सुविधाजनक रूप से toggle करने की अनुमति देता है।

!!! warning

    `univ_hpc` profile ट्रेनिंग वातावरण में ठीक से नहीं चलेगा क्योंकि हमारे पास Slurm scheduler तक पहुंच नहीं है।

यदि भविष्य में हम कॉन्फ़िगरेशन के अन्य तत्व पाते हैं जो हमेशा इनके साथ सह-घटित होते हैं, तो हम बस उन्हें संबंधित profile(s) में जोड़ सकते हैं।
यदि कॉन्फ़िगरेशन के अन्य तत्व हैं जिन्हें हम एक साथ समूहित करना चाहते हैं तो हम अतिरिक्त profiles भी बना सकते हैं।

### 6.2. Test पैरामीटर्स का एक profile बनाना

??? example "परिदृश्य"

    तुम चाहते हो कि अन्य लोग अपने स्वयं के इनपुट डेटा एकत्र किए बिना जल्दी से तुम्हारी pipeline आज़मा सकें।

Profiles केवल infrastructure कॉन्फ़िगरेशन के लिए नहीं हैं।
हम workflow पैरामीटर्स के लिए डिफ़ॉल्ट वैल्यूज़ सेट करने के लिए भी उनका उपयोग कर सकते हैं, ताकि दूसरों के लिए उपयुक्त इनपुट वैल्यूज़ स्वयं एकत्र किए बिना workflow को आज़माना आसान हो जाए।
तुम इसे पैरामीटर फ़ाइल का उपयोग करने के विकल्प के रूप में मान सकते हो।

#### 6.2.1. Profile सेट अप करना

इस संदर्भ में डिफ़ॉल्ट वैल्यूज़ व्यक्त करने के लिए सिंटैक्स इस तरह दिखता है, एक profile के लिए जिसे हम `test` नाम देते हैं:

```groovy title="सिंटैक्स उदाहरण"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

यदि हम अपनी workflow के लिए एक test profile जोड़ते हैं, तो `profiles` block बन जाता है:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

तकनीकी कॉन्फ़िगरेशन profiles की तरह, तुम किसी भी मनमाने नाम के तहत पैरामीटर्स निर्दिष्ट करने वाले कई अलग-अलग profiles सेट अप कर सकते हो जो तुम्हें पसंद हो।

#### 6.2.2. Test profile के साथ workflow को locally चलाना

सुविधाजनक रूप से, profiles परस्पर अनन्य नहीं हैं, इसलिए हम निम्नलिखित सिंटैक्स `-profile <profile1>,<profile2>` (किसी भी संख्या में profiles के लिए) का उपयोग करके अपनी command line में कई profiles निर्दिष्ट कर सकते हैं।

यदि तुम profiles को जोड़ते हो जो कॉन्फ़िगरेशन के समान तत्वों के लिए वैल्यूज़ सेट करते हैं और एक ही कॉन्फ़िगरेशन फ़ाइल में वर्णित हैं, तो Nextflow जो भी वैल्यू इसने अंतिम में पढ़ी (_यानी_ फ़ाइल में बाद में जो कुछ भी आता है) उसका उपयोग करके संघर्ष को हल करेगा।
यदि विरोधाभासी सेटिंग्स विभिन्न कॉन्फ़िगरेशन स्रोतों में सेट हैं, तो डिफ़ॉल्ट [order of precedence](https://www.nextflow.io/docs/latest/config.html) लागू होता है।

चलो अपनी पिछली कमांड में test profile जोड़ने की कोशिश करते हैं:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

यह जहां संभव हो Docker का उपयोग करेगा और `results_config/test` के तहत आउटपुट्स उत्पन्न करेगा, और इस बार character हास्य जोड़ी `dragonandcow` है।

??? abstract "फ़ाइल सामग्री"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

इसका मतलब है कि जब तक हम workflow कोड के साथ कोई भी test data फ़ाइलें वितरित करते हैं, कोई भी command line या पैरामीटर फ़ाइल के माध्यम से अपने स्वयं के इनपुट्स प्रदान किए बिना जल्दी से workflow को आज़मा सकता है।

!!! tip

    हम बाहरी रूप से संग्रहीत बड़ी फ़ाइलों के लिए URLs की ओर इशारा कर सकते हैं।
    Nextflow उन्हें स्वचालित रूप से डाउनलोड करेगा जब तक कि एक open connection है।

    अधिक विवरण के लिए, Side Quest [Working with Files](../side_quests/working_with_files.md) देखें

### 6.3. Resolved कॉन्फ़िगरेशन देखने के लिए `nextflow config` का उपयोग करना

जैसा कि ऊपर उल्लेख किया गया है, कभी-कभी एक ही पैरामीटर को profiles में अलग-अलग वैल्यूज़ पर सेट किया जा सकता है जिन्हें तुम जोड़ना चाहते हो।
और अधिक सामान्य रूप से, कई स्थान हैं जहां कॉन्फ़िगरेशन के तत्वों को संग्रहीत किया जा सकता है, और कभी-कभी समान properties को विभिन्न स्थानों पर अलग-अलग वैल्यूज़ पर सेट किया जा सकता है।

Nextflow किसी भी संघर्ष को हल करने के लिए एक सेट [order of precedence](https://nextflow.io/docs/latest/config.html#configuration-file) लागू करता है, लेकिन वह स्वयं निर्धारित करना मुश्किल हो सकता है।
और भले ही कुछ भी विरोधाभासी न हो, सभी संभावित स्थानों को देखना थकाऊ हो सकता है जहां चीजों को कॉन्फ़िगर किया जा सकता है।

सौभाग्य से, Nextflow में `config` नामक एक सुविधाजनक utility tool शामिल है जो तुम्हारे लिए उस पूरी प्रक्रिया को स्वचालित कर सकता है।

`config` tool तुम्हारी वर्तमान working directory में सभी सामग्री का पता लगाएगा, किसी भी कॉन्फ़िगरेशन फ़ाइलों को इकट्ठा करेगा, और पूरी तरह से resolved कॉन्फ़िगरेशन उत्पन्न करेगा जिसे Nextflow workflow चलाने के लिए उपयोग करेगा।
यह तुम्हें यह पता लगाने की अनुमति देता है कि कुछ भी लॉन्च किए बिना कौन सी सेटिंग्स का उपयोग किया जाएगा।

#### 6.3.1. डिफ़ॉल्ट कॉन्फ़िगरेशन को resolve करना

कॉन्फ़िगरेशन को resolve करने के लिए यह कमांड चलाओ जो डिफ़ॉल्ट रूप से लागू होगी।

```bash
nextflow config
```

??? success "कमांड आउटपुट"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

यह तुम्हें बेस कॉन्फ़िगरेशन दिखाता है जो तुम्हें मिलती है यदि तुम command line में कुछ अतिरिक्त निर्दिष्ट नहीं करते।

#### 6.3.2. विशिष्ट सेटिंग्स सक्रिय के साथ कॉन्फ़िगरेशन को resolve करना

यदि तुम command-line पैरामीटर्स प्रदान करते हो, उदाहरण के लिए एक या अधिक profiles को सक्षम करना या पैरामीटर फ़ाइल लोड करना, तो कमांड अतिरिक्त रूप से उन्हें ध्यान में रखेगी।

```bash
nextflow config -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

यह जटिल प्रोजेक्ट्स के लिए विशेष रूप से उपयोगी हो जाता है जिसमें कॉन्फ़िगरेशन की कई परतें शामिल हैं।

### सारांश

तुम जानते हो कि न्यूनतम परेशानी के साथ runtime पर एक preset कॉन्फ़िगरेशन चुनने के लिए profiles का उपयोग कैसे करें।
अधिक सामान्य रूप से, तुम जानते हो कि विभिन्न compute platforms के अनुरूप अपने workflow executions को कैसे कॉन्फ़िगर करें और अपने विश्लेषणों की पुनरुत्पादनीयता को कैसे बढ़ाएं।

### आगे क्या है?

सीखो कि GitHub जैसे remote repositories से सीधे pipelines कैसे चलाएं।

---

## 7. Remote repositories से pipelines चलाना

??? example "परिदृश्य"

    तुम nf-core जैसी एक अच्छी तरह से स्थापित pipeline चलाना चाहते हो बिना कोड को स्वयं डाउनलोड और मैनेज किए।

अब तक हम वर्तमान डायरेक्टरी में स्थित workflow scripts चला रहे हैं।
व्यवहार में, तुम अक्सर remote repositories में संग्रहीत pipelines चलाना चाहोगे, जैसे GitHub।

Nextflow इसे सीधा बनाता है: तुम पहले manually डाउनलोड किए बिना सीधे Git repository URL से कोई भी pipeline चला सकते हो।

### 7.1. GitHub से एक pipeline चलाना

Remote pipeline चलाने के लिए बुनियादी सिंटैक्स `nextflow run <repository>` है, जहां `<repository>` एक GitHub repository path हो सकता है जैसे `nextflow-io/hello`, एक पूर्ण URL, या GitLab, Bitbucket, या अन्य Git hosting services का पाथ।

आधिकारिक Nextflow "hello" demo pipeline चलाने की कोशिश करो:

```bash
nextflow run nextflow-io/hello
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

पहली बार जब तुम एक remote pipeline चलाते हो, तो Nextflow इसे डाउनलोड करता है और locally cache करता है।
बाद के runs cached संस्करण का उपयोग करते हैं जब तक कि तुम स्पष्ट रूप से अपडेट का अनुरोध नहीं करते।

### 7.2. पुनरुत्पादनीयता के लिए एक संस्करण निर्दिष्ट करना

डिफ़ॉल्ट रूप से, Nextflow डिफ़ॉल्ट branch से नवीनतम संस्करण चलाता है।
तुम `-r` flag का उपयोग करके एक विशेष संस्करण (tag), branch, या commit निर्दिष्ट कर सकते हो:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

सटीक संस्करण निर्दिष्ट करना पुनरुत्पादनीयता के लिए आवश्यक है।

### सारांश

तुम जानते हो कि GitHub और अन्य remote repositories से सीधे pipelines कैसे चलाएं, और पुनरुत्पादनीयता के लिए संस्करण कैसे निर्दिष्ट करें।

### आगे क्या है?

अपने आप को एक बड़ी थपथपाहट दो!
तुम Nextflow pipelines को चलाने और मैनेज करने के लिए शुरू करने के लिए जो कुछ भी जानने की ज़रूरत है वह जानते हो।

यह इस कोर्स का समापन करता है, लेकिन यदि तुम सीखना जारी रखने के लिए उत्सुक हो, तो हमारे पास दो मुख्य सिफारिशें हैं:

- यदि तुम अपनी स्वयं की pipelines विकसित करने में गहराई से जाना चाहते हो, तो [Hello Nextflow](../hello_nextflow/index.md) पर एक नज़र डालो, शुरुआती लोगों के लिए एक कोर्स जो इस एक के समान सामान्य प्रगति को कवर करता है लेकिन channels और operators के बारे में अधिक विस्तार से जाता है।
- यदि तुम कोड में गहराई से गए बिना Nextflow pipelines चलाना सीखना जारी रखना चाहते हो, तो [Hello nf-core](../hello_nf-core/index.md) के पहले भाग पर एक नज़र डालो, जो बेहद लोकप्रिय [nf-core](https://nf-co.re/) प्रोजेक्ट से pipelines खोजने और चलाने के लिए tooling का परिचय देता है।

मज़े करो!

---

## क्विज़

<quiz>
जब पैरामीटर वैल्यूज़ workflow फ़ाइल और `nextflow.config` दोनों में सेट होती हैं, तो किसे प्राथमिकता मिलती है?
- [ ] Workflow फ़ाइल वैल्यू
- [x] कॉन्फ़िगरेशन फ़ाइल वैल्यू
- [ ] पहली वैल्यू जो मिली
- [ ] यह एक error का कारण बनता है

अधिक जानें: [1.1. `nextflow.config` में वैल्यूज़ सेट अप करना](#11-nextflowconfig-में-वैल्यूज़-सेट-अप-करना)
</quiz>

<quiz>
Workflow फ़ाइल बनाम config फ़ाइल में पैरामीटर डिफ़ॉल्ट सेट करने के बीच सिंटैक्स अंतर क्या है?
- [ ] वे समान सिंटैक्स का उपयोग करते हैं
- [x] Workflow typed declaration (`#!groovy param: Type = value`) का उपयोग करता है, config assignment (`#!groovy param = value`) का उपयोग करता है
- [ ] Config typed declaration का उपयोग करता है, workflow assignment का उपयोग करता है
- [ ] केवल config फ़ाइलें डिफ़ॉल्ट वैल्यूज़ सेट कर सकती हैं

अधिक जानें: [1.1. `nextflow.config` में वैल्यूज़ सेट अप करना](#11-nextflowconfig-में-वैल्यूज़-सेट-अप-करना)
</quiz>

<quiz>
Workflow चलाते समय पैरामीटर फ़ाइल कैसे निर्दिष्ट करते हो?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

अधिक जानें: [1.3. पैरामीटर फ़ाइल का उपयोग करना](#13-पैरामीटर-फ़ाइल-का-उपयोग-करना)
</quiz>

<quiz>
`outputDir` कॉन्फ़िगरेशन विकल्प क्या नियंत्रित करता है?
- [ ] Work directory का स्थान
- [x] बेस पाथ जहां workflow आउटपुट्स प्रकाशित होते हैं
- [ ] Log फ़ाइलों के लिए डायरेक्टरी
- [ ] Module फ़ाइलों का स्थान

अधिक जानें: [2.1. `outputDir` डायरेक्टरी नाम को कस्टमाइज़ करना](#21-outputdir-डायरेक्टरी-नाम-को-कस्टमाइज़-करना)
</quiz>

<quiz>
आउटपुट पाथ कॉन्फ़िगरेशन में process नाम को dynamically कैसे संदर्भित करते हो?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

अधिक जानें: [2.2. Process द्वारा आउटपुट्स को व्यवस्थित करना](#22-process-द्वारा-आउटपुट्स-को-व्यवस्थित-करना)
</quiz>

<quiz>
यदि Docker और Conda दोनों सक्षम हैं और एक process में दोनों directives हैं, तो किसे प्राथमिकता दी जाती है?
- [x] Docker (containers)
- [ ] Conda
- [ ] Process में परिभाषित पहला
- [ ] यह एक error का कारण बनता है

अधिक जानें: [3. Software packaging technology का चयन करना](#3-software-packaging-technology-का-चयन-करना)
</quiz>

<quiz>
Nextflow में डिफ़ॉल्ट executor क्या है?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

अधिक जानें: [4. Execution platform का चयन करना](#4-execution-platform-का-चयन-करना)
</quiz>

<quiz>
कौन सी कमांड resource utilization report उत्पन्न करती है?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

अधिक जानें: [5.1. Resource utilization report उत्पन्न करने के लिए workflow चलाना](#51-resource-utilization-report-उत्पन्न-करने-के-लिए-workflow-चलाना)
</quiz>

<quiz>
Config फ़ाइल में `cowpy` नामक एक विशिष्ट process के लिए resource requirements कैसे सेट करते हो?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

अधिक जानें: [5.3. एक विशिष्ट process के लिए resource allocations सेट करना](#53-एक-विशिष्ट-process-के-लिए-resource-allocations-सेट-करना)
</quiz>

<quiz>
`resourceLimits` directive क्या करता है?
- [ ] न्यूनतम resource requirements सेट करता है
- [ ] Processes को resources आवंटित करता है
- [x] अधिकतम resources को cap करता है जिनका अनुरोध किया जा सकता है
- [ ] Real-time में resource उपयोग की निगरानी करता है

अधिक जानें: [5.5. Resource limits जोड़ना](#55-resource-limits-जोड़ना)
</quiz>

<quiz>
एक एकल कमांड में कई profiles कैसे निर्दिष्ट करते हो?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

अधिक जानें: [6. Preset कॉन्फ़िगरेशन के बीच स्विच करने के लिए profiles का उपयोग करना](#6-preset-कॉन्फ़िगरेशन-के-बीच-स्विच-करने-के-लिए-profiles-का-उपयोग-करना)
</quiz>

<quiz>
कौन सी कमांड पूरी तरह से resolved कॉन्फ़िगरेशन दिखाती है जिसे Nextflow उपयोग करेगा?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

अधिक जानें: [6.3. Resolved कॉन्फ़िगरेशन देखने के लिए `nextflow config` का उपयोग करना](#63-resolved-कॉन्फ़िगरेशन-देखने-के-लिए-nextflow-config-का-उपयोग-करना)
</quiz>

<quiz>
Profiles का उपयोग किसके लिए किया जा सकता है? (सभी लागू का चयन करो)
- [x] Infrastructure-विशिष्ट सेटिंग्स परिभाषित करना (executors, containers)
- [x] विभिन्न वातावरणों के लिए resource limits सेट करना
- [x] आसान workflow परीक्षण के लिए test पैरामीटर्स प्रदान करना
- [ ] नई processes परिभाषित करना

अधिक जानें: [6. Preset कॉन्फ़िगरेशन के बीच स्विच करने के लिए profiles का उपयोग करना](#6-preset-कॉन्फ़िगरेशन-के-बीच-स्विच-करने-के-लिए-profiles-का-उपयोग-करना)
</quiz>
