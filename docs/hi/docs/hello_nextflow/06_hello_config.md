# भाग 6: Hello Config

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } [पूरी प्लेलिस्ट](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) Nextflow YouTube चैनल पर देखें।

:green_book: वीडियो ट्रांसक्रिप्ट [यहाँ](./transcripts/06_hello_config.md) उपलब्ध है।
///

यह सेक्शन आपके Nextflow पाइपलाइन के कॉन्फ़िगरेशन को सेट अप और मैनेज करने का तरीका बताएगा ताकि तुम इसके व्यवहार को कस्टमाइज़ कर सको, इसे अलग-अलग वातावरण में अनुकूलित कर सको, और रिसोर्स उपयोग को ऑप्टिमाइज़ कर सको _बिना वर्कफ़्लो कोड की एक भी लाइन बदले_।

ऐसा करने के कई तरीके हैं, जिन्हें संयोजन में उपयोग किया जा सकता है और कॉन्फ़िगरेशन डॉक्यूमेंटेशन में वर्णित [प्राथमिकता के क्रम](https://nextflow.io/docs/latest/config.html) के अनुसार व्याख्या किया जाता है।

कोर्स के इस भाग में, हम तुम्हें सबसे सरल और सबसे सामान्य कॉन्फ़िगरेशन फ़ाइल मैकेनिज़्म दिखाएंगे, [`nextflow.config`](https://nextflow.io/docs/latest/config.html) फ़ाइल, जिसका तुम पहले ही भाग 5: Hello Containers में सामना कर चुके हो।

हम Nextflow कॉन्फ़िगरेशन के आवश्यक घटकों जैसे प्रोसेस डायरेक्टिव, एक्ज़ीक्यूटर, प्रोफ़ाइल और पैरामीटर फ़ाइलों पर चर्चा करेंगे।
इन कॉन्फ़िगरेशन विकल्पों का प्रभावी ढंग से उपयोग करना सीखकर, तुम अपनी पाइपलाइन की लचीलापन, स्केलेबिलिटी और परफ़ॉर्मेंस बढ़ा सकते हो।

??? info "इस सेक्शन से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [Hello Nextflow](./index.md) कोर्स के भाग 1-5 पूरे कर लिए हैं और तुम्हारे पास एक पूर्ण कार्यशील पाइपलाइन है।

    यदि तुम कोर्स को इस बिंदु से शुरू कर रहे हो, तो तुम्हें `modules` डायरेक्टरी और `nextflow.config` फ़ाइल को सॉल्यूशन से कॉपी करना होगा:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    `nextflow.config` फ़ाइल में `docker.enabled = true` लाइन है जो Docker कंटेनर के उपयोग को सक्षम करती है।

    यदि तुम Hello पाइपलाइन से परिचित नहीं हो या तुम्हें रिमाइंडर चाहिए, तो [यह इन्फो पेज](../info/hello_pipeline.md) देखें।

---

## 0. वार्मअप: `hello-config.nf` चलाएं

हम वर्कफ़्लो स्क्रिप्ट `hello-config.nf` को शुरुआती बिंदु के रूप में उपयोग करेंगे।
यह इस ट्रेनिंग कोर्स के भाग 5 को पूरा करके बनाई गई स्क्रिप्ट के बराबर है, सिवाय इसके कि हमने आउटपुट डेस्टिनेशन बदल दिए हैं:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

बस यह सुनिश्चित करने के लिए कि सब कुछ काम कर रहा है, कोई भी बदलाव करने से पहले स्क्रिप्ट को एक बार चलाएं:

```bash
nextflow run hello-config.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

पहले की तरह, तुम्हें आउटपुट फ़ाइलें `output` ब्लॉक में निर्दिष्ट डायरेक्टरी (`results/hello_config/`) में मिलेंगी।

??? abstract "डायरेक्टरी कंटेंट"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

अंतिम ASCII आर्ट आउटपुट `results/hello_config/` डायरेक्टरी में है, `cowpy-COLLECTED-batch-output.txt` नाम के तहत।

??? abstract "फ़ाइल कंटेंट"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

यदि यह तुम्हारे लिए काम कर गया, तो तुम अपनी पाइपलाइन को कॉन्फ़िगर करना सीखने के लिए तैयार हो।

---

## 1. वर्कफ़्लो इनपुट पैरामीटर मैनेज करें

हम कॉन्फ़िगरेशन के एक पहलू से शुरू करेंगे जो बस अब तक हम जो काम कर रहे हैं उसका विस्तार है: इनपुट पैरामीटर का प्रबंधन।

वर्तमान में, हमारा वर्कफ़्लो कमांड-लाइन के माध्यम से कई पैरामीटर वैल्यू स्वीकार करने के लिए सेट अप है, जिसमें वर्कफ़्लो स्क्रिप्ट में ही एक `params` ब्लॉक में डिफ़ॉल्ट वैल्यू सेट हैं।
हालाँकि, तुम उन डिफ़ॉल्ट को ओवरराइड करना चाह सकते हो बिना कमांड लाइन पर पैरामीटर निर्दिष्ट किए, या मूल स्क्रिप्ट फ़ाइल को संशोधित किए।

ऐसा करने के कई तरीके हैं; हम तुम्हें तीन बुनियादी तरीके दिखाएंगे जो बहुत आम तौर पर उपयोग किए जाते हैं।

### 1.1. डिफ़ॉल्ट वैल्यू को `nextflow.config` में मूव करें

यह सबसे सरल दृष्टिकोण है, हालांकि यह संभवतः सबसे कम लचीला है क्योंकि मुख्य `nextflow.config` फ़ाइल कुछ ऐसी नहीं है जिसे तुम हर रन के लिए एडिट करना चाहोगे।
लेकिन इसका फायदा यह है कि यह वर्कफ़्लो में पैरामीटर _घोषित_ करने (जो निश्चित रूप से वहाँ होना चाहिए) बनाम _डिफ़ॉल्ट वैल्यू_ प्रदान करने की चिंताओं को अलग करता है, जो कॉन्फ़िगरेशन फ़ाइल में अधिक उपयुक्त हैं।

चलो इसे दो चरणों में करते हैं।

#### 1.1.1. कॉन्फ़िगरेशन फ़ाइल में एक `params` ब्लॉक बनाएं

`nextflow.config` फ़ाइल में निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline पैरामीटर
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

ध्यान दें कि हमने वर्कफ़्लो से कॉन्फ़िगरेशन फ़ाइल में `params` ब्लॉक को सिर्फ कॉपी नहीं किया।
सिंटैक्स थोड़ा अलग है।
वर्कफ़्लो फ़ाइल में, वे टाइप की गई घोषणाएं हैं।
कॉन्फ़िगरेशन में, वे वैल्यू असाइनमेंट हैं।

तकनीकी रूप से, यह वर्कफ़्लो फ़ाइल में अभी भी निर्दिष्ट डिफ़ॉल्ट वैल्यू को ओवरराइड करने के लिए पर्याप्त है।
तुम कैरेक्टर को संशोधित कर सकते हो, उदाहरण के लिए, और वर्कफ़्लो चला सकते हो यह संतुष्ट करने के लिए कि कॉन्फ़िगरेशन फ़ाइल में सेट की गई वैल्यू वर्कफ़्लो फ़ाइल में सेट की गई वैल्यू को ओवरराइड करती है।

लेकिन कॉन्फ़िगरेशन को पूरी तरह से कॉन्फ़िगरेशन फ़ाइल में मूव करने की भावना में, चलो वर्कफ़्लो फ़ाइल से उन वैल्यू को पूरी तरह से हटा दें।

#### 1.1.2. वर्कफ़्लो फ़ाइल में `params` ब्लॉक से वैल्यू हटाएं

`hello-config.nf` वर्कफ़्लो फ़ाइल में निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "पहले"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

अब वर्कफ़्लो फ़ाइल स्वयं इन पैरामीटर के लिए कोई डिफ़ॉल्ट वैल्यू सेट नहीं करती है।

#### 1.1.3. पाइपलाइन चलाएं

चलो टेस्ट करते हैं कि यह सही तरीके से काम करता है।

```bash
nextflow run hello-config.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

यह अभी भी पहले जैसा ही आउटपुट उत्पन्न करता है।

अंतिम ASCII आर्ट आउटपुट `results/hello_config/` डायरेक्टरी में है, `cowpy-COLLECTED-batch-output.txt` नाम के तहत, पहले की तरह।

??? abstract "फ़ाइल कंटेंट"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

कार्यात्मक रूप से, इस मूव ने कुछ भी नहीं बदला है, लेकिन वैचारिक रूप से कॉन्फ़िगरेशन फ़ाइल में डिफ़ॉल्ट वैल्यू सेट करना थोड़ा साफ है।

### 1.2. रन-स्पेसिफिक कॉन्फ़िगरेशन फ़ाइल का उपयोग करें

यह बढ़िया है, लेकिन कभी-कभी तुम मुख्य कॉन्फ़िगरेशन फ़ाइल के साथ गड़बड़ किए बिना अलग-अलग डिफ़ॉल्ट वैल्यू के साथ कुछ अस्थायी प्रयोग चलाना चाह सकते हो।
तुम एक सबडायरेक्टरी में एक नई `nextflow.config` फ़ाइल बनाकर ऐसा कर सकते हो जिसे तुम अपने प्रयोगों के लिए वर्किंग डायरेक्टरी के रूप में उपयोग करोगे।

#### 1.2.1. ब्लैंक कॉन्फ़िगरेशन के साथ वर्किंग डायरेक्टरी बनाएं

चलो एक नई डायरेक्टरी बनाकर और उसमें जाकर शुरू करते हैं:

```bash
mkdir -p tux-run
cd tux-run
```

फिर, उस डायरेक्टरी में एक ब्लैंक कॉन्फ़िगरेशन फ़ाइल बनाएं:

```bash
touch nextflow.config
```

यह एक खाली फ़ाइल बनाता है।

#### 1.2.2. एक्सपेरिमेंटल कॉन्फ़िगरेशन सेट अप करें

अब नई फ़ाइल खोलें और उन पैरामीटर को जोड़ें जिन्हें तुम कस्टमाइज़ करना चाहते हो:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

ध्यान दें कि इनपुट फ़ाइल का पाथ डायरेक्टरी स्ट्रक्चर को प्रतिबिंबित करना चाहिए।

#### 1.2.3. पाइपलाइन चलाएं

अब हम अपनी नई वर्किंग डायरेक्टरी के भीतर से अपनी पाइपलाइन चला सकते हैं।
पाथ को तदनुसार अनुकूलित करना सुनिश्चित करें!

```bash
nextflow run ../hello-config.nf
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

यह `tux-run/` के तहत डायरेक्टरी का एक नया सेट बनाएगा जिसमें `tux-run/work/` और `tux-run/results/` शामिल हैं।

इस रन में, Nextflow हमारी वर्तमान डायरेक्टरी में `nextflow.config` को पाइपलाइन की रूट डायरेक्टरी में `nextflow.config` के साथ जोड़ता है, और इस तरह डिफ़ॉल्ट कैरेक्टर (turkey) को tux कैरेक्टर के साथ ओवरराइड करता है।

अंतिम आउटपुट फ़ाइल में tux कैरेक्टर को अभिवादन कहते हुए होना चाहिए।

??? abstract "फ़ाइल कंटेंट"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

बस इतना ही; अब तुम्हारे पास अपने 'सामान्य' कॉन्फ़िगरेशन को संशोधित किए बिना प्रयोग करने के लिए एक स्पेस है।

!!! warning

    अगले सेक्शन में जाने से पहले पिछली डायरेक्टरी में वापस जाना सुनिश्चित करें!

    ```bash
    cd ..
    ```

अब चलो पैरामीटर वैल्यू सेट करने का एक और उपयोगी तरीका देखते हैं।

### 1.3. पैरामीटर फ़ाइल का उपयोग करें

सबडायरेक्टरी दृष्टिकोण प्रयोग के लिए बहुत अच्छा काम करता है, लेकिन इसमें थोड़ा सेटअप शामिल है और आवश्यक है कि तुम पाथ को तदनुसार अनुकूलित करो।
जब तुम अपनी पाइपलाइन को वैल्यू के एक विशिष्ट सेट के साथ चलाना चाहते हो, या किसी और को न्यूनतम प्रयास के साथ ऐसा करने में सक्षम बनाना चाहते हो, तो एक सरल दृष्टिकोण है।

Nextflow हमें YAML या JSON फॉर्मेट में एक [पैरामीटर फ़ाइल](https://nextflow.io/docs/latest/config.html#params-file) के माध्यम से पैरामीटर निर्दिष्ट करने की अनुमति देता है, जो डिफ़ॉल्ट वैल्यू के वैकल्पिक सेट को मैनेज और वितरित करना बहुत सुविधाजनक बनाता है, उदाहरण के लिए, साथ ही रन-स्पेसिफिक पैरामीटर वैल्यू।

#### 1.3.1. उदाहरण पैरामीटर फ़ाइल की जांच करें

इसे प्रदर्शित करने के लिए, हम वर्तमान डायरेक्टरी में एक उदाहरण पैरामीटर फ़ाइल प्रदान करते हैं, जिसे `test-params.yaml` कहा जाता है:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

इस पैरामीटर फ़ाइल में प्रत्येक इनपुट के लिए एक की-वैल्यू पेअर है जिसे हम निर्दिष्ट करना चाहते हैं।
यदि तुम सिंटैक्स की तुलना कॉन्फ़िगरेशन फ़ाइल से करते हो तो बराबर चिह्न (`=`) के बजाय कोलन (`:`) के उपयोग पर ध्यान दें।
कॉन्फ़िग फ़ाइल Groovy में लिखी गई है, जबकि पैरामीटर फ़ाइल YAML में लिखी गई है।

!!! info

    हम एक उदाहरण के रूप में पैरामीटर फ़ाइल का JSON संस्करण भी प्रदान करते हैं लेकिन हम यहाँ इसके साथ रन नहीं करेंगे।
    उसे अपने दम पर आज़माने के लिए स्वतंत्र महसूस करो।

#### 1.3.2. पाइपलाइन चलाएं

इस पैरामीटर फ़ाइल के साथ वर्कफ़्लो चलाने के लिए, बस बेस कमांड में `-params-file <filename>` जोड़ें।

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

अंतिम आउटपुट फ़ाइल में stegosaurus कैरेक्टर को अभिवादन कहते हुए होना चाहिए।

??? abstract "फ़ाइल कंटेंट"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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

पैरामीटर फ़ाइल का उपयोग करना ओवरकिल लग सकता है जब तुम्हारे पास निर्दिष्ट करने के लिए केवल कुछ पैरामीटर हैं, लेकिन कुछ पाइपलाइन दर्जनों पैरामीटर की अपेक्षा करती हैं।
उन मामलों में, पैरामीटर फ़ाइल का उपयोग करने से हमें बड़े पैमाने पर कमांड लाइन टाइप किए बिना और वर्कफ़्लो स्क्रिप्ट को संशोधित किए बिना रनटाइम पर पैरामीटर वैल्यू प्रदान करने की अनुमति मिलेगी।

यह सहयोगियों को पैरामीटर के सेट वितरित करना, या उदाहरण के लिए, एक प्रकाशन के लिए सहायक जानकारी के रूप में आसान बनाता है।
यह दूसरों द्वारा तुम्हारे काम को अधिक पुनरुत्पादनीय बनाता है।

### सारांश

तुम वर्कफ़्लो इनपुट को मैनेज करने के लिए प्रमुख कॉन्फ़िगरेशन विकल्पों का लाभ उठाना जानते हो।

### आगे क्या है?

सीखें कि अपने वर्कफ़्लो आउटपुट को कहाँ और कैसे प्रकाशित किया जाए, इसे कैसे मैनेज करें।

---

## 2. वर्कफ़्लो आउटपुट मैनेज करें

अब तक हम वर्कफ़्लो-लेवल आउटपुट घोषणाओं के लिए सभी पाथ को हार्डकोड कर रहे हैं, और जैसा कि हमने कई आउटपुट जोड़ना शुरू किया तब नोट किया, थोड़ी पुनरावृत्ति शामिल हो सकती है।

चलो कुछ सामान्य तरीकों को देखते हैं जिनसे तुम इसे अधिक लचीला बनाने के लिए कॉन्फ़िगर कर सकते हो।

### 2.1. `-output-dir` के साथ आउटपुट डायरेक्टरी को कस्टमाइज़ करें

जब हम नियंत्रित कर रहे हैं कि हमारे 'प्रकाशित' आउटपुट कैसे व्यवस्थित हैं तो हमारी दो अलग प्राथमिकताएं हैं:

- टॉप-लेवल आउटपुट डायरेक्टरी
- इस डायरेक्टरी के भीतर फ़ाइलें कैसे व्यवस्थित हैं

हम अब तक डिफ़ॉल्ट टॉप-लेवल डायरेक्टरी का उपयोग कर रहे हैं: `results`।
चलो `-output-dir` CLI ऑप्शन का उपयोग करके उसे कस्टमाइज़ करके शुरू करते हैं।

#### 2.1.1. `-output-dir` के साथ पाइपलाइन चलाएं

`-output-dir` ऑप्शन (शॉर्टहैंड: `-o`) सभी वर्कफ़्लो आउटपुट के लिए डिफ़ॉल्ट आउटपुट डायरेक्टरी (`results/`) को ओवरराइड करता है।
यह रूट पाथ को नियंत्रित करने का अनुशंसित तरीका है जहाँ आउटपुट प्रकाशित होते हैं।

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

यह `results/` के बजाय `custom-outdir-cli/` में आउटपुट प्रकाशित करता है:

??? abstract "डायरेक्टरी कंटेंट"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

ध्यान दें कि हमारे पास अभी भी output ब्लॉक में `path` घोषणाओं से `hello_config` सबडायरेक्टरी है।
चलो इसे साफ करते हैं।

#### 2.1.2. output ब्लॉक से हार्डकोडेड पाथ हटाएं

`hello_config/` प्रीफ़िक्स पहले के अध्यायों में हार्डकोड किया गया था, लेकिन चूंकि अब हम आउटपुट पाथ को लचीले ढंग से कॉन्फ़िगर करना सीख रहे हैं, हम इस हार्डकोडिंग को हटा सकते हैं।
जिन आउटपुट को सबडायरेक्टरी की आवश्यकता नहीं है, उनके लिए हम `path` डायरेक्टिव को एक खाली स्ट्रिंग पर सेट कर सकते हैं, या इसे पूरी तरह से हटा सकते हैं।

वर्कफ़्लो फ़ाइल में निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

पाइपलाइन फिर से चलाएं:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

अब आउटपुट सीधे `custom-outdir-cli-2/` के तहत प्रकाशित होते हैं, `hello_config` सबडायरेक्टरी के बिना:

??? abstract "डायरेक्टरी कंटेंट"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip

    `-output-dir` ऑप्शन का उपयोग _कहाँ_ आउटपुट जाते हैं यह नियंत्रित करने के लिए किया जाता है, जबकि output ब्लॉक में `path` डायरेक्टिव _सबडायरेक्टरी स्ट्रक्चर_ को नियंत्रित करता है।

### 2.2. डायनामिक आउटपुट पाथ

CLI के माध्यम से आउटपुट डायरेक्टरी बदलने के अलावा, हम `outputDir` का उपयोग करके कॉन्फ़िग फ़ाइल में एक कस्टम डिफ़ॉल्ट वैल्यू भी सेट कर सकते हैं।
यह हमें डायरेक्टरी पाथ को डायनामिक रूप से सेट करने की अनुमति देता है - न केवल स्टैटिक स्ट्रिंग का उपयोग करके।

#### 2.2.1. कॉन्फ़िगरेशन फ़ाइल में `outputDir` सेट करें

`nextflow.config` फ़ाइल में निम्नलिखित कोड जोड़ें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

यह आउटपुट डायरेक्टरी को `custom-outdir-config/` प्लस `batch` पैरामीटर की वैल्यू को सबडायरेक्टरी के रूप में सेट करता है।
अब तुम `--batch` पैरामीटर सेट करके आउटपुट लोकेशन बदल सकते हो:

```bash
nextflow run hello-config.nf --batch my_run
```

यह `custom-outdir-config/my_run/` में आउटपुट प्रकाशित करता है।

!!! note

    `-output-dir` CLI ऑप्शन `outputDir` कॉन्फ़िगरेशन सेटिंग पर प्राथमिकता लेता है।
    यदि यह सेट है, तो कॉन्फ़िग ऑप्शन को पूरी तरह से अनदेखा कर दिया जाएगा।

#### 2.2.2. batch और process नाम के साथ सबडायरेक्टरी

हम प्रति-आउटपुट आधार पर सबडायरेक्टरी आउटपुट `path` घोषणाओं को भी डायनामिक रूप से सेट कर सकते हैं।

उदाहरण के लिए, हम आउटपुट पाथ घोषणा में `<process>.name` का संदर्भ देकर अपने आउटपुट को प्रोसेस द्वारा व्यवस्थित कर सकते हैं:

=== "बाद में"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

हम आगे जा सकते हैं और अधिक जटिल सबडायरेक्टरी पाथ बना सकते हैं।

उपरोक्त एडिट में हमने `intermediates` बनाम टॉप लेवल पर अंतिम आउटपुट के बीच के अंतर को मिटा दिया।
चलो इसे वापस रखते हैं, और फ़ाइलों को `params.batch` सबडायरेक्टरी में भी रखते हैं।

!!! tip

    output ब्लॉक `path` में `params.batch` को शामिल करना, `outputDir` कॉन्फ़िग के बजाय, इसका मतलब है कि यह CLI पर `-output-dir` के साथ ओवरराइट नहीं होगा।

पहले, कॉन्फ़िग फ़ाइल को अपडेट करें ताकि `outputDir` से `${params.batch}` हटा सकें (क्योंकि हम इसे पाथ घोषणाओं में मूव कर रहे हैं):

=== "बाद में"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/"
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

फिर, वर्कफ़्लो फ़ाइल में निम्नलिखित परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "पहले"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. पाइपलाइन चलाएं

चलो देखते हैं कि यह व्यवहार में कैसे काम करता है, कमांड लाइन से `-output-dir` (या संक्षेप में `-o`) को `custom-outdir-config-2` और batch नाम को `rep2` दोनों सेट करते हुए:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

यह `custom-outdir-config-2/rep2/` में आउटपुट प्रकाशित करता है, निर्दिष्ट बेस पाथ _और_ batch नाम सबडायरेक्टरी _और_ प्रोसेस द्वारा समूहीकृत परिणामों के साथ:

??? abstract "डायरेक्टरी कंटेंट"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. वर्कफ़्लो लेवल पर publish mode सेट करें

अंत में, दोहराए जाने वाले कोड की मात्रा को कम करने की भावना में, हम प्रति-आउटपुट `mode` घोषणाओं को कॉन्फ़िगरेशन में एक सिंगल लाइन से बदल सकते हैं।

#### 2.3.1. कॉन्फ़िगरेशन फ़ाइल में `workflow.output.mode` जोड़ें

`nextflow.config` फ़ाइल में निम्नलिखित कोड जोड़ें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/"
    ```

कॉन्फ़िगरेशन फ़ाइल में `workflow.output.mode` सेट करना वर्कफ़्लो फ़ाइल में सेट किए गए को ओवरराइड करने के लिए पर्याप्त है, लेकिन चलो अनावश्यक कोड को वैसे भी हटा दें।

#### 2.3.2. वर्कफ़्लो फ़ाइल से output mode हटाएं

वर्कफ़्लो फ़ाइल में निम्नलिखित परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "पहले"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

यह अधिक संक्षिप्त है, है ना?

#### 2.3.3. पाइपलाइन चलाएं

चलो टेस्ट करते हैं कि यह सही तरीके से काम करता है:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

यह `config-output-mode/` में आउटपुट प्रकाशित करता है, और वे सभी अभी भी उचित कॉपी हैं, symlinks नहीं।

??? abstract "डायरेक्टरी कंटेंट"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

मुख्य कारण जिससे तुम अभी भी प्रति-आउटपुट तरीके से mode सेट करना चाह सकते हो, यदि तुम एक ही वर्कफ़्लो के भीतर मिक्स और मैच करना चाहते हो, _यानी_ कुछ आउटपुट कॉपी किए जाएं और कुछ symlinked हों।

ऐसे बहुत सारे अन्य ऑप्शन हैं जिन्हें तुम इस तरह से कस्टमाइज़ कर सकते हो, लेकिन उम्मीद है कि यह तुम्हें ऑप्शन की रेंज और अपनी प्राथमिकताओं के अनुरूप उन्हें प्रभावी ढंग से उपयोग करने का एहसास देता है।

### सारांश

तुम जानते हो कि उन डायरेक्टरी के नामकरण और संरचना को कैसे नियंत्रित करें जहाँ तुम्हारे आउटपुट प्रकाशित होते हैं, साथ ही वर्कफ़्लो आउटपुट publishing mode भी।

### आगे क्या है?

सीखें कि अपने वर्कफ़्लो कॉन्फ़िगरेशन को अपने कंप्यूट वातावरण में कैसे अनुकूलित करें, सॉफ़्टवेयर पैकेजिंग तकनीक से शुरू करते हुए।

---

## 3. सॉफ़्टवेयर पैकेजिंग तकनीक चुनें

अब तक हम कॉन्फ़िगरेशन तत्वों को देख रहे हैं जो नियंत्रित करते हैं कि इनपुट कैसे अंदर जाते हैं और इनपुट कहाँ बाहर आते हैं। अब यह समय है कि हम अपने कंप्यूट वातावरण में अपने वर्कफ़्लो कॉन्फ़िगरेशन को अनुकूलित करने पर अधिक विशेष रूप से ध्यान केंद्रित करें।

उस पाथ पर पहला कदम यह निर्दिष्ट करना है कि प्रत्येक चरण में चलाए जाने वाले सॉफ़्टवेयर पैकेज कहाँ से आएंगे।
क्या वे पहले से ही लोकल कंप्यूट वातावरण में इंस्टॉल हैं?
क्या हमें इमेज प्राप्त करने और उन्हें कंटेनर सिस्टम के माध्यम से चलाने की आवश्यकता है?
या क्या हमें Conda पैकेज प्राप्त करने और एक लोकल Conda वातावरण बनाने की आवश्यकता है?

ट्रेनिंग कोर्स के पहले भाग (भाग 1-4) में हमने अपने वर्कफ़्लो में बस लोकली इंस्टॉल किए गए सॉफ़्टवेयर का उपयोग किया।
फिर भाग 5 में, हमने Docker कंटेनर पेश किए और `nextflow.config` फ़ाइल, जिसका उपयोग हमने Docker कंटेनर के उपयोग को सक्षम करने के लिए किया।

अब चलो देखते हैं कि हम `nextflow.config` फ़ाइल के माध्यम से एक वैकल्पिक सॉफ़्टवेयर पैकेजिंग ऑप्शन को कैसे कॉन्फ़िगर कर सकते हैं।

### 3.1. कॉन्फ़िग फ़ाइल में Docker को डिसेबल करें और Conda को एनेबल करें

चलो मान लेते हैं कि हम एक HPC क्लस्टर पर काम कर रहे हैं और एडमिन सुरक्षा कारणों से Docker के उपयोग की अनुमति नहीं देता है।
सौभाग्य से हमारे लिए, Nextflow कई अन्य कंटेनर तकनीकों जैसे Singularity (जो HPC पर अधिक व्यापक रूप से उपयोग की जाती है) और सॉफ़्टवेयर पैकेज मैनेजर जैसे Conda का समर्थन करता है।

हम अपनी कॉन्फ़िगरेशन फ़ाइल को Docker के बजाय [Conda](https://nextflow.io/docs/latest/conda.html) का उपयोग करने के लिए बदल सकते हैं।
ऐसा करने के लिए, चलो `docker.enabled` की वैल्यू को `false` में बदलें, और Conda के उपयोग को सक्षम करने वाला एक डायरेक्टिव जोड़ें:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

यह Nextflow को उन प्रोसेस के लिए Conda वातावरण बनाने और उपयोग करने की अनुमति देगा जिनमें Conda पैकेज निर्दिष्ट हैं।
जिसका मतलब है कि अब हमें अपनी `cowpy` प्रोसेस में उनमें से एक जोड़ने की आवश्यकता है!

### 3.2. प्रोसेस डेफिनिशन में Conda पैकेज निर्दिष्ट करें

हमने पहले ही `cowpy` टूल वाले Conda पैकेज के लिए URI प्राप्त कर लिया है: `conda-forge::cowpy==1.1.5`

अब हम `conda` डायरेक्टिव का उपयोग करके `cowpy` प्रोसेस डेफिनिशन में URI जोड़ते हैं:

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

स्पष्ट होने के लिए, हम `docker` डायरेक्टिव को _बदल_ नहीं रहे हैं, हम एक वैकल्पिक ऑप्शन _जोड़_ रहे हैं।

!!! tip

    किसी दिए गए conda पैकेज के लिए URI प्राप्त करने के कुछ अलग तरीके हैं।
    हम [Seqera Containers](https://seqera.io/containers/) सर्च क्वेरी का उपयोग करने की सलाह देते हैं, जो तुम्हें एक URI देगा जिसे तुम कॉपी और पेस्ट कर सकते हो, भले ही तुम इससे कंटेनर बनाने की योजना नहीं बना रहे हो।

### 3.3. वर्कफ़्लो चलाएं यह सत्यापित करने के लिए कि यह Conda का उपयोग कर सकता है

चलो इसे आज़माते हैं।

```bash
nextflow run hello-config.nf --batch conda
```

??? success "कमांड आउटपुट"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

यह बिना किसी समस्या के काम करना चाहिए और `custom-outdir-config/conda` के तहत पहले जैसे ही आउटपुट उत्पन्न करना चाहिए।

पर्दे के पीछे, Nextflow ने Conda पैकेज प्राप्त किए हैं और वातावरण बनाया है, जो सामान्य रूप से थोड़ा काम लेता है; तो यह अच्छा है कि हमें वह सब कुछ खुद नहीं करना पड़ता!

!!! note

    यह जल्दी चलता है क्योंकि `cowpy` पैकेज काफी छोटा है, लेकिन यदि तुम बड़े पैकेज के साथ काम कर रहे हो, तो पहली बार सामान्य से थोड़ा अधिक समय लग सकता है, और तुम कंसोल आउटपुट को पूरा होने से पहले एक मिनट या उससे अधिक समय तक 'अटका हुआ' देख सकते हो।
    यह सामान्य है और पहली बार जब तुम एक नया पैकेज उपयोग करते हो तो Nextflow द्वारा किए गए अतिरिक्त काम के कारण है।

हमारे दृष्टिकोण से, ऐसा लगता है कि यह Docker के साथ चलने के बिल्कुल समान काम करता है, भले ही बैकएंड पर मैकेनिक्स थोड़े अलग हों।

इसका मतलब है कि हम आवश्यकता पड़ने पर Conda वातावरण के साथ चलाने के लिए तैयार हैं।

??? info "Docker और Conda को मिक्स और मैच करना"

    चूंकि ये डायरेक्टिव प्रति प्रोसेस असाइन किए जाते हैं, इसलिए 'मिक्स और मैच' करना संभव है, _यानी_ अपने वर्कफ़्लो में कुछ प्रोसेस को Docker के साथ और अन्य को Conda के साथ चलाने के लिए कॉन्फ़िगर करना, उदाहरण के लिए, यदि तुम जिस कंप्यूट इंफ्रास्ट्रक्चर का उपयोग कर रहे हो वह दोनों का समर्थन करता है।
    उस स्थिति में, तुम अपनी कॉन्फ़िगरेशन फ़ाइल में Docker और Conda दोनों को सक्षम करोगे।
    यदि किसी दिए गए प्रोसेस के लिए दोनों उपलब्ध हैं, तो Nextflow कंटेनर को प्राथमिकता देगा।

    और जैसा कि पहले उल्लेख किया गया है, Nextflow कई अन्य सॉफ़्टवेयर पैकेजिंग और कंटेनर तकनीकों का समर्थन करता है, इसलिए तुम केवल उन दोनों तक सीमित नहीं हो।

### सारांश

तुम जानते हो कि प्रत्येक प्रोसेस को किस सॉफ़्टवेयर पैकेज का उपयोग करना चाहिए यह कैसे कॉन्फ़िगर करें, और तकनीकों के बीच कैसे स्विच करें।

### आगे क्या है?

सीखें कि Nextflow द्वारा वास्तव में काम करने के लिए उपयोग किए जाने वाले execution प्लेटफ़ॉर्म को कैसे बदलें।

---

## 4. execution प्लेटफ़ॉर्म चुनें

अब तक, हम अपनी पाइपलाइन को लोकल एक्ज़ीक्यूटर के साथ चला रहे हैं।
यह प्रत्येक कार्य को उस मशीन पर निष्पादित करता है जिस पर Nextflow चल रहा है।
जब Nextflow शुरू होता है, तो यह उपलब्ध CPU और मेमोरी को देखता है।
यदि चलाने के लिए तैयार कार्यों के रिसोर्स उपलब्ध रिसोर्स से अधिक हैं, तो Nextflow अंतिम कार्यों को execution से रोक देगा जब तक कि एक या अधिक पहले के कार्य समाप्त नहीं हो जाते, आवश्यक रिसोर्स को मुक्त करते हुए।

लोकल एक्ज़ीक्यूटर सुविधाजनक और कुशल है, लेकिन यह उस सिंगल मशीन तक सीमित है। बहुत बड़े वर्कलोड के लिए, तुम पा सकते हो कि तुम्हारी लोकल मशीन एक बाधा है, या तो क्योंकि तुम्हारे पास एक सिंगल कार्य है जिसे तुम्हारे पास उपलब्ध से अधिक रिसोर्स की आवश्यकता है, या क्योंकि तुम्हारे पास इतने सारे कार्य हैं कि एक सिंगल मशीन के लिए उन्हें चलाने की प्रतीक्षा करना बहुत लंबा समय लेगा।

Nextflow [कई अलग-अलग एक्ज़ीक्यूटर](https://nextflow.io/docs/latest/executor.html) का समर्थन करता है, जिसमें HPC शेड्यूलर (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor और अन्य) के साथ-साथ क्लाउड execution बैकएंड (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes और अधिक) शामिल हैं।

### 4.1. एक अलग बैकएंड को टार्गेट करना

एक्ज़ीक्यूटर की पसंद एक प्रोसेस डायरेक्टिव द्वारा सेट की जाती है जिसे `executor` कहा जाता है।
डिफ़ॉल्ट रूप से यह `local` पर सेट है, इसलिए निम्नलिखित कॉन्फ़िगरेशन निहित है:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

एक अलग बैकएंड को टार्गेट करने के लिए एक्ज़ीक्यूटर सेट करने के लिए, तुम बस उसी सिंटैक्स का उपयोग करके जो एक्ज़ीक्यूटर चाहते हो उसे निर्दिष्ट करोगे जैसा कि रिसोर्स आवंटन के लिए ऊपर वर्णित है (सभी ऑप्शन के लिए [executor डॉक्यूमेंटेशन](https://nextflow.io/docs/latest/executor.html) देखें)।

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    हम वास्तव में ट्रेनिंग वातावरण में इसे टेस्ट नहीं कर सकते क्योंकि यह HPC से कनेक्ट करने के लिए सेट अप नहीं है।

### 4.2. execution पैरामीटर के लिए बैकएंड-स्पेसिफिक सिंटैक्स से निपटना

अधिकांश हाई-परफ़ॉर्मेंस कंप्यूटिंग प्लेटफ़ॉर्म अनुमति देते हैं (और कभी-कभी आवश्यक होते हैं) कि तुम कुछ पैरामीटर निर्दिष्ट करो जैसे रिसोर्स आवंटन अनुरोध और सीमाएं (उदाहरण के लिए CPU की संख्या और मेमोरी) और उपयोग करने के लिए job queue का नाम।

दुर्भाग्य से, इनमें से प्रत्येक सिस्टम एक job को कैसे परिभाषित और संबंधित शेड्यूलर को सबमिट किया जाना चाहिए, इसके लिए अलग-अलग तकनीकों, सिंटैक्स और कॉन्फ़िगरेशन का उपयोग करता है।

??? abstract "उदाहरण"

    उदाहरण के लिए, 8 CPU और 4GB RAM की आवश्यकता वाली एक ही job को queue "my-science-work" पर निष्पादित करने के लिए बैकएंड के आधार पर निम्नलिखित अलग-अलग तरीकों से व्यक्त करने की आवश्यकता है।

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

सौभाग्य से, Nextflow यह सब सरल बनाता है।
यह एक मानकीकृत सिंटैक्स प्रदान करता है ताकि तुम संबंधित गुणों जैसे [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) और [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) को बस एक बार निर्दिष्ट कर सको (अन्य गुणों के लिए [process directives](https://nextflow.io/docs/latest/reference/process.html#process-directives) देखें)।
फिर, रनटाइम पर, Nextflow उन सेटिंग्स का उपयोग executor सेटिंग के आधार पर उपयुक्त बैकएंड-स्पेसिफिक स्क्रिप्ट उत्पन्न करने के लिए करेगा।

हम अगले सेक्शन में उस मानकीकृत सिंटैक्स को कवर करेंगे।

### सारांश

तुम अब जानते हो कि विभिन्न प्रकार के कंप्यूटिंग इंफ्रास्ट्रक्चर का उपयोग करने के लिए एक्ज़ीक्यूटर को कैसे बदलें।

### आगे क्या है?

सीखें कि Nextflow में रिसोर्स आवंटन और सीमाओं का मूल्यांकन और व्यक्त कैसे करें।

---

## 5. कंप्यूट रिसोर्स आवंटन को नियंत्रित करें

अधिकांश हाई-परफ़ॉर्मेंस कंप्यूटिंग प्लेटफ़ॉर्म अनुमति देते हैं (और कभी-कभी आवश्यक होते हैं) कि तुम कुछ रिसोर्स आवंटन पैरामीटर निर्दिष्ट करो जैसे CPU की संख्या और मेमोरी।

डिफ़ॉल्ट रूप से, Nextflow प्रत्येक प्रोसेस के लिए एक सिंगल CPU और 2GB मेमोरी का उपयोग करेगा।
संबंधित प्रोसेस डायरेक्टिव को `cpus` और `memory` कहा जाता है, इसलिए निम्नलिखित कॉन्फ़िगरेशन निहित है:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

तुम इन वैल्यू को संशोधित कर सकते हो, या तो सभी प्रोसेस के लिए या विशिष्ट नामित प्रोसेस के लिए, अपनी कॉन्फ़िगरेशन फ़ाइल में अतिरिक्त प्रोसेस डायरेक्टिव का उपयोग करके।
Nextflow उन्हें चुने गए एक्ज़ीक्यूटर के लिए उपयुक्त निर्देशों में अनुवाद करेगा।

लेकिन तुम कैसे जानते हो कि किन वैल्यू का उपयोग करना है?

### 5.1. रिसोर्स उपयोग रिपोर्ट उत्पन्न करने के लिए वर्कफ़्लो चलाएं

यदि तुम पहले से नहीं जानते कि तुम्हारी प्रोसेस को कितनी CPU और मेमोरी की आवश्यकता होने की संभावना है, तो तुम कुछ रिसोर्स प्रोफाइलिंग कर सकते हो, जिसका अर्थ है कि तुम कुछ डिफ़ॉल्ट आवंटन के साथ वर्कफ़्लो चलाते हो, रिकॉर्ड करते हो कि प्रत्येक प्रोसेस ने कितना उपयोग किया, और वहाँ से, बेस आवंटन को कैसे समायोजित करना है इसका अनुमान लगाते हो।

सुविधाजनक रूप से, Nextflow में ऐसा करने के लिए बिल्ट-इन टूल शामिल हैं, और अनुरोध पर तुम्हारे लिए खुशी से एक रिपोर्ट उत्पन्न करेगा।

ऐसा करने के लिए, अपनी कमांड लाइन में `-with-report <filename>.html` जोड़ें।

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

रिपोर्ट एक html फ़ाइल है, जिसे तुम डाउनलोड कर सकते हो और अपने ब्राउज़र में खोल सकते हो। तुम इसे बाईं ओर फ़ाइल एक्सप्लोरर में राइट क्लिक भी कर सकते हो और ट्रेनिंग वातावरण में इसे देखने के लिए `Show preview` पर क्लिक कर सकते हो।

रिपोर्ट को देखने में कुछ मिनट लें और देखें कि क्या तुम रिसोर्स को समायोजित करने के लिए कुछ अवसरों की पहचान कर सकते हो।
उन टैब पर क्लिक करना सुनिश्चित करें जो आवंटित की गई चीज़ के प्रतिशत के रूप में उपयोग परिणाम दिखाते हैं।

सभी उपलब्ध सुविधाओं पर डॉक्यूमेंटेशन के लिए [Reports](https://nextflow.io/docs/latest/reports.html) देखें।

### 5.2. सभी प्रोसेस के लिए रिसोर्स आवंटन सेट करें

प्रोफाइलिंग दिखाती है कि हमारे ट्रेनिंग वर्कफ़्लो में प्रोसेस बहुत हल्के हैं, तो चलो प्रति प्रोसेस डिफ़ॉल्ट मेमोरी आवंटन को 1GB तक कम करें।

अपनी `nextflow.config` फ़ाइल में निम्नलिखित जोड़ें, पाइपलाइन पैरामीटर सेक्शन से पहले:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * प्रोसेस सेटिंग्स
    */
    process {
        memory = 1.GB
    }

    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    conda.enabled = true

    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

यह हमारे द्वारा उपभोग किए जाने वाले कंप्यूट की मात्रा को कम करने में मदद करेगा।

### 5.3. एक विशिष्ट प्रोसेस के लिए रिसोर्स आवंटन सेट करें

उसी समय, हम मान लेंगे कि `cowpy` प्रोसेस को अन्य की तुलना में अधिक रिसोर्स की आवश्यकता है, बस इसलिए कि हम प्रदर्शित कर सकें कि एक व्यक्तिगत प्रोसेस के लिए आवंटन को कैसे समायोजित करें।

=== "बाद में"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * प्रोसेस सेटिंग्स
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
    * प्रोसेस सेटिंग्स
    */
    process {
        memory = 1.GB
    }
    ```

इस कॉन्फ़िगरेशन के साथ, सभी प्रोसेस 1GB मेमोरी और एक सिंगल CPU (निहित डिफ़ॉल्ट) का अनुरोध करेंगी, सिवाय `cowpy` प्रोसेस के, जो 2GB और 2 CPU का अनुरोध करेगी।

!!! tip

    यदि तुम्हारे पास कम CPU वाली मशीन है और तुम प्रति प्रोसेस एक उच्च संख्या आवंटित करते हो, तो तुम प्रोसेस कॉल को एक दूसरे के पीछे कतारबद्ध होते हुए देख सकते हो।
    ऐसा इसलिए है क्योंकि Nextflow सुनिश्चित करता है कि हम उपलब्ध से अधिक CPU का अनुरोध नहीं करते हैं।

### 5.4. अपडेट किए गए कॉन्फ़िगरेशन के साथ वर्कफ़्लो चलाएं

चलो इसे आज़माते हैं, प्रोफाइलिंग रिपोर्ट के लिए एक अलग फ़ाइलनाम प्रदान करते हुए ताकि हम कॉन्फ़िगरेशन परिवर्तनों से पहले और बाद में परफ़ॉर्मेंस की तुलना कर सकें।

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

तुम शायद कोई वास्तविक अंतर नोटिस नहीं करोगे क्योंकि यह इतना छोटा वर्कलोड है, लेकिन यह वह दृष्टिकोण है जिसका उपयोग तुम वास्तविक दुनिया के वर्कफ़्लो के परफ़ॉर्मेंस और रिसोर्स आवश्यकताओं का विश्लेषण करने के लिए करोगे।

यह बहुत उपयोगी है जब तुम्हारी प्रोसेस की अलग-अलग रिसोर्स आवश्यकताएं हों। यह तुम्हें वास्तविक डेटा के आधार पर प्रत्येक प्रोसेस के लिए सेट किए गए रिसोर्स आवंटन को सही आकार देने में सक्षम बनाता है, अनुमान नहीं।

!!! tip

    यह बस एक छोटा सा स्वाद है कि तुम रिसोर्स के अपने उपयोग को ऑप्टिमाइज़ करने के लिए क्या कर सकते हो।
    Nextflow स्वयं में कुछ वास्तव में साफ [डायनामिक रिट्राई लॉजिक](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) बिल्ट-इन है जो रिसोर्स सीमाओं के कारण विफल होने वाली job को रिट्राई करने के लिए है।
    इसके अतिरिक्त, Seqera Platform तुम्हारे रिसोर्स आवंटन को स्वचालित रूप से ऑप्टिमाइज़ करने के लिए AI-संचालित टूलिंग भी प्रदान करता है।

### 5.5. रिसोर्स सीमाएं जोड़ें

तुम किस कंप्यूटिंग एक्ज़ीक्यूटर और कंप्यूट इंफ्रास्ट्रक्चर का उपयोग कर रहे हो, इस पर निर्भर करते हुए, कुछ बाधाएं हो सकती हैं कि तुम क्या आवंटित कर सकते हो (या करना चाहिए)।
उदाहरण के लिए, तुम्हारा क्लस्टर तुम्हें कुछ सीमाओं के भीतर रहने की आवश्यकता कर सकता है।

तुम संबंधित सीमाओं को सेट करने के लिए `resourceLimits` डायरेक्टिव का उपयोग कर सकते हो। सिंटैक्स इस तरह दिखता है जब यह प्रोसेस ब्लॉक में अकेला होता है:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow इन वैल्यू को तुम्हारे द्वारा निर्दिष्ट एक्ज़ीक्यूटर के आधार पर उपयुक्त निर्देशों में अनुवाद करेगा।

हम इसे नहीं चलाएंगे, क्योंकि हमारे पास ट्रेनिंग वातावरण में संबंधित इंफ्रास्ट्रक्चर तक पहुंच नहीं है।
हालाँकि, यदि तुम इन सीमाओं से अधिक रिसोर्स आवंटन के साथ वर्कफ़्लो चलाने की कोशिश करते हो, फिर `.command.run` स्क्रिप्ट फ़ाइल में `sbatch` कमांड को देखते हो, तो तुम देखोगे कि वास्तव में एक्ज़ीक्यूटर को भेजे जाने वाले अनुरोध `resourceLimits` द्वारा निर्दिष्ट वैल्यू पर कैप किए गए हैं।

??? info "संस्थागत संदर्भ कॉन्फ़िगरेशन"

    nf-core प्रोजेक्ट ने दुनिया भर के विभिन्न संस्थानों द्वारा साझा की गई [कॉन्फ़िगरेशन फ़ाइलों का संग्रह](https://nf-co.re/configs/) संकलित किया है, जो HPC और क्लाउड एक्ज़ीक्यूटर की एक विस्तृत श्रृंखला को कवर करता है।

    वे साझा कॉन्फ़िग उन लोगों के लिए मूल्यवान हैं जो वहाँ काम करते हैं और इसलिए बस अपने संस्थान के कॉन्फ़िगरेशन का उपयोग बॉक्स से बाहर कर सकते हैं, और उन लोगों के लिए एक मॉडल के रूप में जो अपने स्वयं के इंफ्रास्ट्रक्चर के लिए कॉन्फ़िगरेशन विकसित करना चाहते हैं।

### सारांश

तुम जानते हो कि रिसोर्स उपयोग का आकलन करने के लिए प्रोफाइलिंग रिपोर्ट कैसे उत्पन्न करें और सभी प्रोसेस और/या व्यक्तिगत प्रोसेस के लिए रिसोर्स आवंटन को कैसे संशोधित करें, साथ ही HPC पर चलाने के लिए रिसोर्स सीमाएं सेट करें।

### आगे क्या है?

सीखें कि प्रीसेट कॉन्फ़िगरेशन प्रोफ़ाइल कैसे सेट अप करें और रनटाइम पर उनके बीच कैसे स्विच करें।

---

## 6. प्रीसेट कॉन्फ़िगरेशन के बीच स्विच करने के लिए प्रोफ़ाइल का उपयोग करें

हमने तुम्हें कई तरीके दिखाए हैं जिनसे तुम अपनी पाइपलाइन कॉन्फ़िगरेशन को कस्टमाइज़ कर सकते हो, यह इस बात पर निर्भर करता है कि तुम किस प्रोजेक्ट पर काम कर रहे हो या तुम किस कंप्यूट वातावरण का उपयोग कर रहे हो।

तुम वैकल्पिक सेटिंग्स के बीच स्विच करना चाह सकते हो, यह इस बात पर निर्भर करता है कि तुम किस कंप्यूटिंग इंफ्रास्ट्रक्चर का उपयोग कर रहे हो। उदाहरण के लिए, तुम अपने लैपटॉप पर लोकली छोटे पैमाने के टेस्ट विकसित और चलाना चाह सकते हो, फिर HPC या क्लाउड पर पूर्ण पैमाने के वर्कलोड चला सकते हो।

Nextflow तुम्हें किसी भी संख्या में [प्रोफ़ाइल](https://nextflow.io/docs/latest/config.html#config-profiles) सेट अप करने देता है जो विभिन्न कॉन्फ़िगरेशन का वर्णन करते हैं, जिन्हें तुम फिर कमांड-लाइन आर्गुमेंट का उपयोग करके रनटाइम पर चुन सकते हो, बजाय कॉन्फ़िगरेशन फ़ाइल को स्वयं संशोधित करने के।

### 6.1. लोकल डेवलपमेंट और HPC पर execution के बीच स्विच करने के लिए प्रोफ़ाइल बनाएं

चलो दो वैकल्पिक प्रोफ़ाइल सेट अप करते हैं; एक नियमित कंप्यूटर पर छोटे पैमाने के लोड चलाने के लिए, जहाँ हम Docker कंटेनर का उपयोग करेंगे, और एक Slurm शेड्यूलर के साथ यूनिवर्सिटी HPC पर चलाने के लिए, जहाँ हम Conda पैकेज का उपयोग करेंगे।

#### 6.1.1. प्रोफ़ाइल सेट अप करें

अपनी `nextflow.config` फ़ाइल में निम्नलिखित जोड़ें, पाइपलाइन पैरामीटर सेक्शन के बाद लेकिन आउटपुट सेटिंग्स से पहले:

=== "बाद में"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * प्रोफ़ाइल
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

    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "पहले"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * आउटपुट सेटिंग्स
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

तुम देखते हो कि यूनिवर्सिटी HPC के लिए, हम रिसोर्स सीमाएं भी निर्दिष्ट कर रहे हैं।

#### 6.1.2. प्रोफ़ाइल के साथ वर्कफ़्लो चलाएं

हमारी Nextflow कमांड लाइन में प्रोफ़ाइल निर्दिष्ट करने के लिए, हम `-profile` आर्गुमेंट का उपयोग करते हैं।

चलो `my_laptop` कॉन्फ़िगरेशन के साथ वर्कफ़्लो चलाने की कोशिश करते हैं।

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

जैसा कि तुम देख सकते हो, यह हमें रनटाइम पर कॉन्फ़िगरेशन के बीच बहुत सुविधाजनक रूप से टॉगल करने की अनुमति देता है।

!!! warning

    `univ_hpc` प्रोफ़ाइल ट्रेनिंग वातावरण में ठीक से नहीं चलेगी क्योंकि हमारे पास Slurm शेड्यूलर तक पहुंच नहीं है।

यदि भविष्य में हम कॉन्फ़िगरेशन के अन्य तत्व पाते हैं जो हमेशा इनके साथ सह-घटित होते हैं, तो हम बस उन्हें संबंधित प्रोफ़ाइल(लों) में जोड़ सकते हैं।
यदि कॉन्फ़िगरेशन के अन्य तत्व हैं जिन्हें हम एक साथ समूहित करना चाहते हैं तो हम अतिरिक्त प्रोफ़ाइल भी बना सकते हैं।

### 6.2. टेस्ट पैरामीटर की प्रोफ़ाइल बनाएं

प्रोफ़ाइल केवल इंफ्रास्ट्रक्चर कॉन्फ़िगरेशन के लिए नहीं हैं।
हम उनका उपयोग वर्कफ़्लो पैरामीटर के लिए डिफ़ॉल्ट वैल्यू सेट करने के लिए भी कर सकते हैं, ताकि दूसरों के लिए खुद उपयुक्त इनपुट वैल्यू इकट्ठा किए बिना वर्कफ़्लो को आज़माना आसान हो सके।
तुम इसे पैरामीटर फ़ाइल का उपयोग करने के विकल्प के रूप में मान सकते हो।

#### 6.2.1. प्रोफ़ाइल सेट अप करें

इस संदर्भ में डिफ़ॉल्ट वैल्यू व्यक्त करने के लिए सिंटैक्स इस तरह दिखता है, एक प्रोफ़ाइल के लिए जिसे हम `test` नाम देते हैं:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

यदि हम अपने वर्कफ़्लो के लिए एक टेस्ट प्रोफ़ाइल जोड़ते हैं, तो `profiles` ब्लॉक बन जाता है:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* प्रोफ़ाइल
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

तकनीकी कॉन्फ़िगरेशन प्रोफ़ाइल की तरह, तुम किसी भी मनमाने नाम के तहत पैरामीटर निर्दिष्ट करने वाली कई अलग-अलग प्रोफ़ाइल सेट अप कर सकते हो।

#### 6.2.2. टेस्ट प्रोफ़ाइल के साथ वर्कफ़्लो लोकली चलाएं

सुविधाजनक रूप से, प्रोफ़ाइल परस्पर अनन्य नहीं हैं, इसलिए हम अपनी कमांड लाइन में निम्नलिखित सिंटैक्स का उपयोग करके कई प्रोफ़ाइल निर्दिष्ट कर सकते हैं `-profile <profile1>,<profile2>` (किसी भी संख्या में प्रोफ़ाइल के लिए)।

यदि तुम ऐसी प्रोफ़ाइल को संयोजित करते हो जो कॉन्फ़िगरेशन के समान तत्वों के लिए वैल्यू सेट करती हैं और एक ही कॉन्फ़िगरेशन फ़ाइल में वर्णित हैं, तो Nextflow जो भी वैल्यू उसने अंतिम में पढ़ी (_यानी_ फ़ाइल में जो बाद में आती है) उसका उपयोग करके संघर्ष को हल करेगा।
यदि विरोधाभासी सेटिंग्स विभिन्न कॉन्फ़िगरेशन स्रोतों में सेट हैं, तो डिफ़ॉल्ट [प्राथमिकता का क्रम](https://nextflow.io/docs/latest/config.html) लागू होता है।

चलो हमारी पिछली कमांड में टेस्ट प्रोफ़ाइल जोड़ने की कोशिश करते हैं:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

यह जहाँ संभव हो Docker का उपयोग करेगा और `custom-outdir-config/test` के तहत आउटपुट उत्पन्न करेगा, और इस बार कैरेक्टर कॉमेडिक जोड़ी `dragonandcow` है।

??? abstract "फ़ाइल कंटेंट"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

इसका मतलब है कि जब तक हम वर्कफ़्लो कोड के साथ किसी भी टेस्ट डेटा फ़ाइलों को वितरित करते हैं, कोई भी कमांड लाइन या पैरामीटर फ़ाइल के माध्यम से अपने स्वयं के इनपुट प्रदान किए बिना जल्दी से वर्कफ़्लो को आज़मा सकता है।

!!! tip

    हम बाहरी रूप से संग्रहीत बड़ी फ़ाइलों के लिए URL की ओर इशारा कर सकते हैं।
    Nextflow उन्हें स्वचालित रूप से डाउनलोड करेगा जब तक कि एक खुला कनेक्शन हो।

    अधिक विवरण के लिए, साइड क्वेस्ट [Working with Files](../side_quests/working_with_files.md) देखें

### 6.3. रिज़ॉल्व्ड कॉन्फ़िगरेशन देखने के लिए `nextflow config` का उपयोग करें

जैसा कि ऊपर उल्लेख किया गया है, कभी-कभी एक ही पैरामीटर को उन प्रोफ़ाइल में अलग-अलग वैल्यू पर सेट किया जा सकता है जिन्हें तुम संयोजित करना चाहते हो।
और अधिक सामान्य रूप से, कई स्थान हैं जहाँ कॉन्फ़िगरेशन के तत्व संग्रहीत किए जा सकते हैं, और कभी-कभी समान गुणों को अलग-अलग स्थानों पर अलग-अलग वैल्यू पर सेट किया जा सकता है।

Nextflow किसी भी संघर्ष को हल करने के लिए एक सेट [प्राथमिकता का क्रम](https://nextflow.io/docs/latest/config.html) लागू करता है, लेकिन वह खुद निर्धारित करना मुश्किल हो सकता है।
और भले ही कुछ भी विरोधाभासी न हो, सभी संभावित स्थानों को देखना थकाऊ हो सकता है जहाँ चीजें कॉन्फ़िगर की जा सकती हैं।

सौभाग्य से, Nextflow में `config` नामक एक सुविधाजनक उपयोगिता टूल शामिल है जो तुम्हारे लिए उस पूरी प्रक्रिया को स्वचालित कर सकता है।

`config` टूल तुम्हारी वर्तमान वर्किंग डायरेक्टरी में सभी कंटेंट का पता लगाएगा, किसी भी कॉन्फ़िगरेशन फ़ाइलों को इकट्ठा करेगा, और पूरी तरह से रिज़ॉल्व्ड कॉन्फ़िगरेशन उत्पन्न करेगा जिसका उपयोग Nextflow वर्कफ़्लो चलाने के लिए करेगा।
यह तुम्हें यह पता लगाने की अनुमति देता है कि कुछ भी लॉन्च किए बिना किन सेटिंग्स का उपयोग किया जाएगा।

#### 6.3.1. डिफ़ॉल्ट कॉन्फ़िगरेशन को रिज़ॉल्व करें

उस कॉन्फ़िगरेशन को रिज़ॉल्व करने के लिए यह कमांड चलाएं जो डिफ़ॉल्ट रूप से लागू होगा।

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

यह तुम्हें बेस कॉन्फ़िगरेशन दिखाता है जो तुम्हें मिलता है यदि तुम कमांड लाइन में कुछ अतिरिक्त निर्दिष्ट नहीं करते हो।

#### 6.3.2. विशिष्ट सेटिंग्स सक्रिय के साथ कॉन्फ़िगरेशन को रिज़ॉल्व करें

यदि तुम कमांड-लाइन पैरामीटर प्रदान करते हो, उदाहरण के लिए एक या अधिक प्रोफ़ाइल को सक्षम करना या पैरामीटर फ़ाइल लोड करना, तो कमांड अतिरिक्त रूप से उन्हें ध्यान में रखेगी।

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

यह जटिल प्रोजेक्ट के लिए विशेष रूप से उपयोगी हो जाता है जिसमें कॉन्फ़िगरेशन की कई परतें शामिल हैं।

### सारांश

तुम जानते हो कि न्यूनतम परेशानी के साथ रनटाइम पर प्रीसेट कॉन्फ़िगरेशन चुनने के लिए प्रोफ़ाइल का उपयोग कैसे करें।
अधिक सामान्य रूप से, तुम जानते हो कि विभिन्न कंप्यूट प्लेटफ़ॉर्म के अनुरूप अपने वर्कफ़्लो execution को कैसे कॉन्फ़िगर करें और अपने विश्लेषण की पुनरुत्पादनीयता को कैसे बढ़ाएं।

### आगे क्या है?

जश्न मनाओ और अपनी पीठ थपथपाओ! तुमने अपना पहला Nextflow डेवलपर कोर्स पूरा कर लिया है।

अंतिम [कोर्स सारांश](./next_steps.md) पर जाएं यह समीक्षा करने के लिए कि तुमने क्या सीखा और आगे क्या आता है यह जानने के लिए।

---

## क्विज़

<quiz>
कॉन्फ़िगरेशन फ़ाइल का नाम क्या है जिसे Nextflow स्वचालित रूप से लोड करता है?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
जब एक ही पैरामीटर कॉन्फ़िग फ़ाइल और कमांड लाइन दोनों में सेट हो तो क्या प्राथमिकता लेता है?
- [ ] कॉन्फ़िग फ़ाइल वैल्यू
- [x] कमांड लाइन वैल्यू
- [ ] पहली वैल्यू जो मिली
- [ ] न तो; यह एक एरर का कारण बनता है

अधिक जानें: [1.1. डिफ़ॉल्ट वैल्यू को `nextflow.config` में मूव करें](#11-move-default-values-to-nextflowconfig)
</quiz>

<quiz>
क्या तुम एक ही कॉन्फ़िगरेशन में Docker और Conda दोनों को सक्षम कर सकते हो?
- [x] हाँ, Nextflow प्रोसेस डायरेक्टिव के आधार पर दोनों का उपयोग कर सकता है
- [ ] नहीं, एक समय में केवल एक ही सक्षम किया जा सकता है
- [ ] हाँ, लेकिन केवल प्रोफ़ाइल में
- [ ] नहीं, वे परस्पर अनन्य हैं
</quiz>

<quiz>
यदि Docker और Conda दोनों सक्षम हैं और एक प्रोसेस में दोनों डायरेक्टिव हैं, तो किसे प्राथमिकता दी जाती है?
- [x] Docker (कंटेनर)
- [ ] Conda
- [ ] पहला जो परिभाषित किया गया
- [ ] यह एक एरर का कारण बनता है

अधिक जानें: [3. सॉफ़्टवेयर पैकेजिंग तकनीक चुनें](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Nextflow प्रोसेस के लिए डिफ़ॉल्ट मेमोरी आवंटन क्या है?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] कोई सीमा नहीं
</quiz>

<quiz>
कॉन्फ़िग फ़ाइल में एक विशिष्ट प्रोसेस के लिए रिसोर्स आवश्यकताएं कैसे सेट करें?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

अधिक जानें: [5.3. एक विशिष्ट प्रोसेस के लिए रिसोर्स आवंटन सेट करें](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
कौन सा कमांड लाइन ऑप्शन रिसोर्स उपयोग रिपोर्ट उत्पन्न करता है?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

अधिक जानें: [5.1. रिसोर्स उपयोग रिपोर्ट उत्पन्न करने के लिए वर्कफ़्लो चलाएं](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
`resourceLimits` डायरेक्टिव क्या करता है?
- [ ] न्यूनतम रिसोर्स आवश्यकताएं सेट करता है
- [ ] प्रोसेस को रिसोर्स आवंटित करता है
- [x] अधिकतम रिसोर्स को कैप करता है जिनका अनुरोध किया जा सकता है
- [ ] रिसोर्स उपयोग की निगरानी करता है

अधिक जानें: [5.5. रिसोर्स सीमाएं जोड़ें](#55-add-resource-limits)
</quiz>

<quiz>
Nextflow में डिफ़ॉल्ट एक्ज़ीक्यूटर क्या है?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

अधिक जानें: [4. execution प्लेटफ़ॉर्म चुनें](#4-select-an-execution-platform)
</quiz>

<quiz>
Nextflow चलाते समय पैरामीटर फ़ाइल कैसे निर्दिष्ट करें?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

अधिक जानें: [1.3. पैरामीटर फ़ाइल का उपयोग करें](#13-use-a-parameter-file)
</quiz>

<quiz>
प्रोफ़ाइल का उपयोग किसके लिए किया जा सकता है? (सभी लागू चुनें)
- [x] इंफ्रास्ट्रक्चर-स्पेसिफिक सेटिंग्स परिभाषित करना
- [x] विभिन्न वातावरण के लिए रिसोर्स सीमाएं सेट करना
- [x] टेस्ट पैरामीटर प्रदान करना
- [ ] नई प्रोसेस परिभाषित करना

अधिक जानें: [6. प्रीसेट कॉन्फ़िगरेशन के बीच स्विच करने के लिए प्रोफ़ाइल का उपयोग करें](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
एक सिंगल कमांड में कई प्रोफ़ाइल कैसे निर्दिष्ट करें?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

अधिक जानें: [6. प्रीसेट कॉन्फ़िगरेशन के बीच स्विच करने के लिए प्रोफ़ाइल का उपयोग करें](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
