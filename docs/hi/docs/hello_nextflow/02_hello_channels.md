# पार्ट 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Nextflow YouTube चैनल पर [पूरी प्लेलिस्ट](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) देखें।

:green_book: वीडियो ट्रांसक्रिप्ट [यहाँ](./transcripts/02_hello_channels.md) उपलब्ध है।
///

इस कोर्स के पार्ट 1 (Hello World) में, हमने तुम्हें दिखाया कि प्रोसेस कॉल में सीधे इनपुट देकर एक प्रोसेस को वेरिएबल इनपुट कैसे प्रदान करें: `sayHello(params.input)`।
यह जानबूझकर सरलीकृत दृष्टिकोण था।
व्यावहारिक रूप से, इस दृष्टिकोण की बड़ी सीमाएं हैं; अर्थात् यह केवल बहुत सरल मामलों के लिए काम करता है जहां हम प्रोसेस को केवल एक बार, एक ही वैल्यू पर चलाना चाहते हैं।
अधिकांश यथार्थवादी वर्कफ़्लो उपयोग के मामलों में, हम कई वैल्यू को प्रोसेस करना चाहते हैं (उदाहरण के लिए, कई नमूनों के लिए प्रयोगात्मक डेटा), इसलिए हमें इनपुट को संभालने के लिए अधिक परिष्कृत तरीके की आवश्यकता है।

यही वह है जिसके लिए Nextflow [**channels**](https://nextflow.io/docs/latest/channel.html) हैं।
चैनल ऐसी कतारें हैं जो इनपुट को कुशलता से संभालने और बहु-चरणीय वर्कफ़्लो में एक चरण से दूसरे चरण में उन्हें शटल करने के लिए डिज़ाइन की गई हैं, जबकि अंतर्निहित समानांतरता और कई अतिरिक्त लाभ प्रदान करती हैं।

इस कोर्स के इस भाग में, तुम सीखोगे कि विभिन्न स्रोतों से कई इनपुट को संभालने के लिए चैनल का उपयोग कैसे करें।
तुम यह भी सीखोगे कि चैनल की सामग्री को आवश्यकतानुसार रूपांतरित करने के लिए [**operators**](https://nextflow.io/docs/latest/reference/operator.html) का उपयोग कैसे करें।

??? info "इस सेक्शन से कैसे शुरू करें"

    कोर्स का यह सेक्शन मानता है कि तुमने [Hello Nextflow](./index.md) कोर्स का पार्ट 1 पूरा कर लिया है, लेकिन यदि तुम उस सेक्शन में शामिल मूल बातों से सहज हो, तो तुम बिना कुछ विशेष किए यहां से शुरू कर सकते हो।

---

## 0. वार्मअप: `hello-channels.nf` चलाएं

हम वर्कफ़्लो स्क्रिप्ट `hello-channels.nf` को शुरुआती बिंदु के रूप में उपयोग करने जा रहे हैं।
यह इस प्रशिक्षण कोर्स के पार्ट 1 के माध्यम से काम करके तैयार की गई स्क्रिप्ट के बराबर है, सिवाय इसके कि हमने आउटपुट गंतव्य बदल दिया है:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

बस यह सुनिश्चित करने के लिए कि सब कुछ काम कर रहा है, कोई भी बदलाव करने से पहले स्क्रिप्ट को एक बार चलाएं:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

पहले की तरह, तुम्हें `results/hello_channels` डायरेक्टरी में `output.txt` नाम की आउटपुट फ़ाइल मिलेगी (जैसा कि वर्कफ़्लो स्क्रिप्ट के `output` ब्लॉक में निर्दिष्ट है, ऊपर दिखाया गया है)।

??? abstract "डायरेक्टरी सामग्री"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

यदि यह तुम्हारे लिए काम किया, तो तुम चैनल के बारे में सीखने के लिए तैयार हो।

---

## 1. चैनल के माध्यम से स्पष्ट रूप से वेरिएबल इनपुट प्रदान करें

हम `sayHello()` प्रोसेस को वेरिएबल इनपुट पास करने के लिए एक **channel** बनाने जा रहे हैं, बजाय इसके कि अंतर्निहित हैंडलिंग पर निर्भर रहें, जिसकी कुछ सीमाएं हैं।

### 1.1. एक इनपुट चैनल बनाएं

विभिन्न प्रकार के [**channel factories**](https://nextflow.io/docs/latest/reference/channel.html) हैं जिनका उपयोग हम चैनल सेट अप करने के लिए कर सकते हैं।
अभी के लिए चीजों को सरल रखने के लिए, हम सबसे बुनियादी चैनल फैक्ट्री का उपयोग करने जा रहे हैं, जिसे [`channel.of`](https://nextflow.io/docs/latest/reference/channel.html#of) कहा जाता है, जो एक ही वैल्यू वाला चैनल बनाएगा।
कार्यात्मक रूप से यह पहले की तरह होगा, लेकिन Nextflow को अंतर्निहित रूप से चैनल बनाने के बजाय, हम अब इसे स्पष्ट रूप से कर रहे हैं।

यह वह कोड की लाइन है जिसका हम उपयोग करने जा रहे हैं:

```console title="सिंटैक्स"
greeting_ch = channel.of('Hello Channels!')
```

यह `channel.of()` चैनल फैक्ट्री का उपयोग करके `greeting_ch` नामक एक चैनल बनाता है, जो एक सरल queue चैनल सेट अप करता है, और अभिवादन वैल्यू के रूप में उपयोग करने के लिए स्ट्रिंग `'Hello Channels!'` लोड करता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note

    हम पठनीयता के लिए CLI पैरामीटर का उपयोग करने के बजाय अस्थायी रूप से हार्डकोडेड स्ट्रिंग पर वापस स्विच कर रहे हैं। एक बार जब हम चैनल के स्तर पर क्या हो रहा है, यह कवर कर लेंगे, तो हम CLI पैरामीटर का उपयोग करने पर वापस जाएंगे।

वर्कफ़्लो ब्लॉक में, चैनल फैक्ट्री कोड जोड़ें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

यह अभी तक कार्यात्मक नहीं है क्योंकि हमने अभी तक प्रोसेस कॉल में इनपुट को स्विच नहीं किया है।

### 1.2. प्रोसेस कॉल में इनपुट के रूप में चैनल जोड़ें

अब हमें वास्तव में अपने नए बनाए गए चैनल को `sayHello()` प्रोसेस कॉल में प्लग करना होगा, CLI पैरामीटर को बदलते हुए जिसे हम पहले सीधे प्रदान कर रहे थे।

वर्कफ़्लो ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

यह Nextflow को `greeting_ch` चैनल की सामग्री पर `sayHello` प्रोसेस चलाने के लिए कहता है।

अब हमारा वर्कफ़्लो ठीक से कार्यात्मक है; यह `sayHello('Hello Channels!')` लिखने के बराबर है।

### 1.3. वर्कफ़्लो चलाएं

चलो इसे चलाते हैं!

```bash
nextflow run hello-channels.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

यदि तुमने दोनों संपादन सही ढंग से किए, तो तुम्हें एक सफल निष्पादन मिलना चाहिए।
तुम परिणाम डायरेक्टरी की जांच कर सकते हो कि परिणाम अभी भी पहले जैसा ही है।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

तो हमने अपने वर्कफ़्लो की लचीलापन बढ़ा दी है जबकि वही अंतिम परिणाम प्राप्त किया है।
यह ऐसा लग सकता है कि हम बिना किसी ठोस लाभ के अधिक कोड लिख रहे हैं, लेकिन जैसे ही हम अधिक इनपुट को संभालना शुरू करेंगे, मूल्य स्पष्ट हो जाएगा।

इसके पूर्वावलोकन के रूप में, आगे बढ़ने से पहले एक और चीज़ देखते हैं: डेटा इनपुट को प्रबंधित करने के लिए एक स्पष्ट चैनल का उपयोग करने का एक छोटा लेकिन सुविधाजनक लाभ।

### 1.4. चैनल सामग्री का निरीक्षण करने के लिए `view()` का उपयोग करें

Nextflow चैनल इस तरह से बनाए गए हैं कि हम ऑपरेटरों का उपयोग करके उनकी सामग्री पर काम कर सकते हैं, जिसे हम इस अध्याय में बाद में विस्तार से कवर करेंगे।

अभी के लिए, हम तुम्हें केवल यह दिखाने जा रहे हैं कि चैनल की सामग्री का निरीक्षण करने के लिए [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) नामक एक सुपर सरल ऑपरेटर का उपयोग कैसे करें।
तुम `view()` को डिबगिंग टूल के रूप में सोच सकते हो, जैसे Python में `print()` स्टेटमेंट, या अन्य भाषाओं में इसके समकक्ष।

वर्कफ़्लो ब्लॉक में यह छोटी लाइन जोड़ें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello Channels!')
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

स्पेस की सटीक मात्रा मायने नहीं रखती जब तक यह 4 का गुणक है; हम बस `.view()` स्टेटमेंट की शुरुआत को चैनल निर्माण के `.of()` भाग के साथ संरेखित करने का लक्ष्य रख रहे हैं।

अब वर्कफ़्लो को फिर से चलाएं:

```bash
nextflow run hello-channels.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

जैसा कि तुम देख सकते हो, यह चैनल सामग्री को कंसोल पर आउटपुट करता है।
यहां हमारे पास केवल एक तत्व है, लेकिन जब हम अगले सेक्शन में चैनल में कई वैल्यू लोड करना शुरू करेंगे, तो तुम देखोगे कि यह प्रति पंक्ति एक तत्व आउटपुट करने के लिए सेट है।

### सारांश

तुम जानते हो कि एक प्रोसेस को इनपुट प्रदान करने के लिए एक बुनियादी चैनल फैक्ट्री का उपयोग कैसे करें।

### आगे क्या है?

सीखें कि कई इनपुट वैल्यू पर वर्कफ़्लो को पुनरावृत्त करने के लिए चैनल का उपयोग कैसे करें।

---

## 2. कई इनपुट वैल्यू पर चलने के लिए वर्कफ़्लो को संशोधित करें

वर्कफ़्लो आमतौर पर इनपुट के बैच पर चलते हैं जिन्हें थोक में संसाधित किया जाना है, इसलिए हम कई इनपुट वैल्यू स्वीकार करने के लिए वर्कफ़्लो को अपग्रेड करना चाहते हैं।

### 2.1. इनपुट चैनल में कई अभिवादन लोड करें

सुविधाजनक रूप से, `channel.of()` चैनल फैक्ट्री जिसका हम उपयोग कर रहे हैं, एक से अधिक वैल्यू स्वीकार करने में काफी खुश है, इसलिए हमें इसे बिल्कुल भी संशोधित करने की आवश्यकता नहीं है।
हम बस चैनल में कई वैल्यू लोड कर सकते हैं।

चलो उन्हें `'Hello'`, `'Bonjour'` और `'Holà'` बनाते हैं।

#### 2.1.1. अधिक अभिवादन जोड़ें

वर्कफ़्लो ब्लॉक से पहले, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // इनपुट के लिए एक चैनल बनाएं
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // इनपुट के लिए एक चैनल बनाएं
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

दस्तावेज़ हमें बताता है कि यह काम करना चाहिए। क्या यह वास्तव में इतना सरल हो सकता है?

#### 2.1.2. कमांड चलाएं और लॉग आउटपुट देखें

चलो इसे आज़माते हैं।

```bash
nextflow run hello-channels.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

यह निश्चित रूप से ठीक चला है।
निष्पादन मॉनिटर दिखाता है कि `sayHello` प्रोसेस के लिए `3 of 3` कॉल किए गए, और हम `view()` स्टेटमेंट द्वारा गिने गए तीन अभिवादन देखते हैं, प्रति पंक्ति एक जैसा वादा किया गया था।

हालांकि, परिणाम डायरेक्टरी में अभी भी केवल एक आउटपुट है:

??? abstract "डायरेक्टरी सामग्री"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

तुम्हें वहां तीन अभिवादनों में से एक दिखना चाहिए, हालांकि तुम्हें जो मिला वह यहां दिखाए गए से अलग हो सकता है।
क्या तुम सोच सकते हो कि ऐसा क्यों हो सकता है?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_आरेख में, चैनल को हरे रंग में दर्शाया गया है, और तत्वों का क्रम पाइप में मार्बल्स की तरह दर्शाया गया है: पहला लोड किया गया दाईं ओर है, फिर दूसरा बीच में, फिर तीसरा बाईं ओर है।_

निष्पादन मॉनिटर को देखते हुए, इसने हमें केवल एक सबडायरेक्टरी पथ (`f4/c9962c`) दिया।
चलो वहां देखते हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

यह वह अभिवादन भी नहीं है जो हमें परिणाम डायरेक्टरी में मिला! क्या हो रहा है?

इस बिंदु पर, हमें तुम्हें बताना होगा कि डिफ़ॉल्ट रूप से, ANSI लॉगिंग सिस्टम एक ही प्रोसेस के कई कॉल से लॉगिंग को एक ही लाइन पर लिखता है।
तो sayHello() प्रोसेस के सभी तीन कॉल की स्थिति एक ही स्थान पर आ रही है।

सौभाग्य से, हम प्रोसेस कॉल की पूरी सूची देखने के लिए उस व्यवहार को अक्षम कर सकते हैं।

#### 2.1.3. `-ansi-log false` विकल्प के साथ कमांड फिर से चलाएं

प्रति प्रोसेस कॉल एक लाइन प्रदर्शित करने के लिए लॉगिंग का विस्तार करने के लिए, कमांड में `-ansi-log false` जोड़ें।

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

इस बार हम सभी तीन प्रोसेस रन और उनकी संबंधित कार्य सबडायरेक्टरी आउटपुट में सूचीबद्ध देखते हैं।

यह बहुत बेहतर है, कम से कम एक सरल वर्कफ़्लो के लिए।
एक जटिल वर्कफ़्लो के लिए, या बड़ी संख्या में इनपुट के लिए, टर्मिनल पर आउटपुट की पूरी सूची थोड़ी भारी हो जाएगी।
यही कारण है कि `-ansi-log false` डिफ़ॉल्ट व्यवहार नहीं है।

!!! tip

    दोनों लॉगिंग मोड के बीच स्थिति की रिपोर्ट करने का तरीका थोड़ा अलग है।
    संघनित मोड में, Nextflow रिपोर्ट करता है कि कॉल सफलतापूर्वक पूर्ण हुए या नहीं।
    इस विस्तारित मोड में, यह केवल रिपोर्ट करता है कि वे सबमिट किए गए थे।

वैसे भी, अब जब हमारे पास प्रत्येक प्रोसेस कॉल की सबडायरेक्टरी है, तो हम उनके लॉग और आउटपुट की तलाश कर सकते हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "फ़ाइल सामग्री"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

यह दिखाता है कि सभी तीन प्रोसेस सफलतापूर्वक चले (यय)।

उस ने कहा, हमारे पास अभी भी समस्या है कि परिणाम डायरेक्टरी में केवल एक आउटपुट फ़ाइल है।

तुम्हें याद होगा कि हमने `sayHello` प्रोसेस के लिए आउटपुट फ़ाइल नाम हार्डकोड किया था, इसलिए सभी तीन कॉल ने `output.txt` नामक एक फ़ाइल तैयार की।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-task-dirs.svg"
</figure>

जब तक आउटपुट फ़ाइलें कार्य सबडायरेक्टरी में रहती हैं, अन्य प्रोसेस से अलग, यह ठीक है।
लेकिन जब उन्हें एक ही परिणाम डायरेक्टरी में प्रकाशित किया जाता है, तो जो भी पहले वहां कॉपी हुआ वह अगले द्वारा अधिलेखित हो जाता है, और इसी तरह।

### 2.2. सुनिश्चित करें कि आउटपुट फ़ाइल नाम अद्वितीय होंगे

हम सभी आउटपुट को एक ही परिणाम डायरेक्टरी में प्रकाशित करना जारी रख सकते हैं, लेकिन हमें यह सुनिश्चित करने की आवश्यकता है कि उनके अद्वितीय नाम होंगे।
विशेष रूप से, हमें पहली प्रोसेस को गतिशील रूप से एक फ़ाइल नाम उत्पन्न करने के लिए संशोधित करने की आवश्यकता है ताकि अंतिम फ़ाइल नाम अद्वितीय हों।

तो हम फ़ाइल नामों को अद्वितीय कैसे बनाएं?
ऐसा करने का एक सामान्य तरीका आउटपुट फ़ाइल नाम के हिस्से के रूप में इनपुट (इनपुट चैनल से प्राप्त) से कुछ अद्वितीय मेटाडेटा का उपयोग करना है।
यहां, सुविधा के लिए, हम केवल अभिवादन का उपयोग करेंगे क्योंकि यह सिर्फ एक छोटी स्ट्रिंग है, और इसे बेस आउटपुट फ़ाइलनाम में जोड़ेंगे।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi-unique.svg"
</figure>

#### 2.2.1. एक गतिशील आउटपुट फ़ाइल नाम बनाएं

प्रोसेस ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

आउटपुट परिभाषा और `script:` कमांड ब्लॉक दोनों में `output.txt` को बदलना सुनिश्चित करें।

!!! tip

    आउटपुट परिभाषा में, तुम्हें आउटपुट फ़ाइलनाम एक्सप्रेशन के चारों ओर डबल कोट्स का उपयोग करना चाहिए (सिंगल कोट्स नहीं), अन्यथा यह विफल हो जाएगा।

यह हर बार प्रोसेस को कॉल करने पर एक अद्वितीय आउटपुट फ़ाइल नाम उत्पन्न करना चाहिए, ताकि इसे आउटपुट डायरेक्टरी में एक ही प्रोसेस के अन्य कॉल के आउटपुट से अलग किया जा सके।

#### 2.2.2. वर्कफ़्लो चलाएं

चलो इसे चलाते हैं। ध्यान दें कि हम डिफ़ॉल्ट ANSI लॉग सेटिंग्स के साथ चलाने के लिए वापस आ गए हैं।

```bash
nextflow run hello-channels.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

सारांश दृश्य पर वापस लौटते हुए, आउटपुट फिर से एक पंक्ति पर संक्षेपित है।
यह देखने के लिए `results` डायरेक्टरी देखें कि क्या सभी आउटपुट अभिवादन वहां हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

हाँ! और उनमें से प्रत्येक में अपेक्षित सामग्री है।

??? abstract "फ़ाइल सामग्री"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

सफलता! अब हम आउटपुट फ़ाइलों के अधिलेखित होने की चिंता किए बिना जितने चाहें उतने अभिवादन जोड़ सकते हैं।

!!! tip

    व्यवहार में, इनपुट डेटा के आधार पर फ़ाइलों का नामकरण लगभग हमेशा अव्यावहारिक होता है।
    गतिशील फ़ाइलनाम उत्पन्न करने का बेहतर तरीका इनपुट फ़ाइलों के साथ एक प्रोसेस को मेटाडेटा पास करना है।
    मेटाडेटा आमतौर पर 'sample sheet' या समकक्षों के माध्यम से प्रदान किया जाता है।
    तुम अपने Nextflow प्रशिक्षण में बाद में सीखोगे कि यह कैसे करें (देखें [Metadata side quest](../side_quests/metadata.md))।

### सारांश

तुम जानते हो कि एक चैनल के माध्यम से कई इनपुट तत्वों को कैसे फीड करें।

### आगे क्या है?

चैनल की सामग्री को रूपांतरित करने के लिए एक ऑपरेटर का उपयोग करना सीखें।

---

## 3. एक ऐरे के माध्यम से कई इनपुट प्रदान करें

हमने अभी तुम्हें दिखाया कि कई इनपुट तत्वों को कैसे संभालें जो सीधे चैनल फैक्ट्री में हार्डकोड किए गए थे।
यदि हम उन कई इनपुट को एक अलग तरीके से प्रदान करना चाहते हैं तो क्या होगा?

उदाहरण के लिए, कल्पना करें कि हम इस तरह तत्वों की एक ऐरे वाला एक इनपुट वेरिएबल सेट अप करते हैं:

`greetings_array = ['Hello','Bonjour','Holà']`

क्या हम इसे अपने आउटपुट चैनल में लोड कर सकते हैं और इसके काम करने की उम्मीद कर सकते हैं?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

चलो पता लगाते हैं।

### 3.1. चैनल के इनपुट के रूप में वैल्यू की एक ऐरे प्रदान करें

सामान्य ज्ञान बताता है कि हमें एक ही वैल्यू के बजाय वैल्यू की एक ऐरे पास करने में सक्षम होना चाहिए।
चलो इसे आज़माते हैं; हमें इनपुट वेरिएबल सेट अप करना होगा और इसे चैनल फैक्ट्री में लोड करना होगा।

#### 3.1.1. इनपुट वेरिएबल सेट अप करें

चलो `greetings_array` वेरिएबल लेते हैं जिसकी हमने अभी कल्पना की थी और इसे वर्कफ़्लो ब्लॉक में जोड़कर वास्तविकता बनाते हैं:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यह अभी तक कार्यात्मक नहीं है, हमने बस ऐरे के लिए एक घोषणा जोड़ी है।

#### 3.1.2. चैनल फैक्ट्री के इनपुट के रूप में अभिवादनों की ऐरे सेट करें

अब हम चैनल फैक्ट्री में वर्तमान में हार्डकोड की गई वैल्यू `'Hello','Bonjour','Holà'` को अभी बनाए गए `greetings_array` से बदलने जा रहे हैं।

वर्कफ़्लो ब्लॉक में, निम्नलिखित परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यह अब कार्यात्मक होना चाहिए।

#### 3.1.3. वर्कफ़्लो चलाएं

चलो इसे चलाने की कोशिश करते हैं:

```bash
nextflow run hello-channels.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

ओह नहीं! एक त्रुटि है!

`view()` के आउटपुट और त्रुटि संदेशों को देखें।

ऐसा लगता है कि Nextflow ने `[Hello, Bonjour, Holà]` को एक ही स्ट्रिंग वैल्यू के रूप में उपयोग करते हुए एक ही प्रोसेस कॉल चलाने की कोशिश की, बजाय ऐरे में तीन स्ट्रिंग को अलग वैल्यू के रूप में उपयोग करने के।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-fail.svg"
</figure>

तो यह 'पैकेजिंग' है जो समस्या पैदा कर रही है।
हम Nextflow को ऐरे को अनपैक करने और व्यक्तिगत स्ट्रिंग को चैनल में लोड करने के लिए कैसे कहें?

### 3.2. चैनल सामग्री को रूपांतरित करने के लिए एक ऑपरेटर का उपयोग करें

यहीं पर [**operators**](https://nextflow.io/docs/latest/reference/operator.html) काम में आते हैं।
तुमने पहले ही `.view()` ऑपरेटर का उपयोग किया है, जो बस देखता है कि वहां क्या है।
अब हम ऐसे ऑपरेटरों को देखने जा रहे हैं जो हमें चैनल की सामग्री पर कार्य करने की अनुमति देते हैं।

यदि तुम Nextflow दस्तावेज़ में [ऑपरेटरों की सूची](https://nextflow.io/docs/latest/reference/operator.html) को स्किम करते हो, तो तुम्हें [`flatten()`](https://nextflow.io/docs/latest/reference/operator.html#flatten) मिलेगा, जो बिल्कुल वही करता है जो हमें चाहिए: एक ऐरे की सामग्री को अनपैक करें और उन्हें व्यक्तिगत आइटम के रूप में emit करें।

#### 3.2.1. `flatten()` ऑपरेटर जोड़ें

हमारे इनपुट चैनल पर `flatten()` ऑपरेटर लागू करने के लिए, हम इसे पहले की तरह चैनल फैक्ट्री घोषणा में जोड़ते हैं।

वर्कफ़्लो ब्लॉक में, `flatten()` को `splitcsv()` (अनकमेंटेड) से बदलने के लिए निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

यहां हमने पठनीयता के लिए अगली पंक्ति पर ऑपरेटर जोड़ा, लेकिन यदि तुम चाहो तो चैनल फैक्ट्री के रूप में एक ही पंक्ति पर ऑपरेटर जोड़ सकते हो, इस तरह:
`greeting_ch = channel.of(greetings_array).view().flatten()`

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-array-success.svg"
</figure>

#### 3.2.2. `view()` स्टेटमेंट को परिष्कृत करें

हम यह परीक्षण करने के लिए इसे तुरंत चला सकते हैं कि यह काम करता है या नहीं, लेकिन जब हम इस पर हैं, तो हम यह परिष्कृत करने जा रहे हैं कि हम चैनल सामग्री का निरीक्षण कैसे करते हैं।

हम `flatten()` ऑपरेटर लागू होने से पहले और बाद में सामग्री कैसी दिखती है, इसकी तुलना करने में सक्षम होना चाहते हैं, इसलिए हम दूसरा जोड़ने जा रहे हैं, और हम आउटपुट में उन्हें अधिक स्पष्ट रूप से लेबल करने के लिए थोड़ा कोड जोड़ने जा रहे हैं।

वर्कफ़्लो ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "flatten से पहले: $greeting" }
                             .flatten()
                             .view { greeting -> "flatten के बाद: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

तुम देखते हो कि हमने दूसरा `.view` स्टेटमेंट जोड़ा है, और उनमें से प्रत्येक के लिए, हमने खाली कोष्ठक (`()`) को कुछ कोड वाले घुंघराले ब्रेसिज़ से बदल दिया है, जैसे `{ greeting -> "flatten से पहले: $greeting" }`।

इन्हें _closures_ कहा जाता है। उनमें मौजूद कोड चैनल में प्रत्येक आइटम के लिए निष्पादित किया जाएगा।
हम आंतरिक वैल्यू के लिए एक अस्थायी वेरिएबल परिभाषित करते हैं, यहां `greeting` कहा जाता है (लेकिन यह कोई भी मनमाना नाम हो सकता है), जिसका उपयोग केवल उस closure के दायरे में किया जाता है।

इस उदाहरण में, `$greeting` चैनल में लोड किए गए प्रत्येक व्यक्तिगत आइटम का प्रतिनिधित्व करता है।
इसके परिणामस्वरूप साफ-सुथरा लेबल वाला कंसोल आउटपुट होगा।

!!! info

    कुछ पाइपलाइनों में तुम ऑपरेटर closures के अंदर उपयोग किए जाने वाले `$it` नामक एक विशेष वेरिएबल देख सकते हो।
    यह एक _implicit_ वेरिएबल है जो आंतरिक वेरिएबल तक शॉर्ट-हैंड एक्सेस की अनुमति देता है,
    बिना इसे `->` के साथ परिभाषित करने की आवश्यकता के।

    हम कोड स्पष्टता में सहायता के लिए स्पष्ट होना पसंद करते हैं, इसलिए `$it` सिंटैक्स को हतोत्साहित किया जाता है और धीरे-धीरे Nextflow भाषा से चरणबद्ध किया जाएगा।

#### 3.2.3. वर्कफ़्लो चलाएं

अंत में, तुम वर्कफ़्लो को फिर से चलाने की कोशिश कर सकते हो!

```bash
nextflow run hello-channels.nf
```

??? success "कमांड आउटपुट"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    flatten से पहले: [Hello, Bonjour, Holà]
    flatten के बाद: Hello
    flatten के बाद: Bonjour
    flatten के बाद: Holà
    ```

इस बार यह काम करता है और हमें `flatten()` ऑपरेटर चलाने से पहले और बाद में चैनल की सामग्री कैसी दिखती है, इसमें अतिरिक्त अंतर्दृष्टि देता है।

- एक ही `flatten से पहले:` स्टेटमेंट क्योंकि उस बिंदु पर चैनल में एक आइटम होता है, मूल ऐरे।
- तीन अलग `flatten के बाद:` स्टेटमेंट, प्रत्येक अभिवादन के लिए एक, जो अब चैनल में व्यक्तिगत आइटम हैं।

महत्वपूर्ण रूप से, इसका मतलब है कि अब प्रत्येक आइटम को वर्कफ़्लो द्वारा अलग से संसाधित किया जा सकता है।

!!! tip

    तकनीकी रूप से एक अलग चैनल फैक्ट्री, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist) का उपयोग करके समान परिणाम प्राप्त करना संभव है, जिसमें इसके संचालन में एक अंतर्निहित मैपिंग चरण शामिल है।
    यहां हमने एक सरल उपयोग के मामले पर एक ऑपरेटर के उपयोग को प्रदर्शित करने के लिए इसका उपयोग नहीं करना चुना।

### सारांश

तुम जानते हो कि चैनल की सामग्री को रूपांतरित करने के लिए `flatten()` जैसे ऑपरेटर का उपयोग कैसे करें, और ऑपरेटर लागू करने से पहले और बाद में चैनल सामग्री का निरीक्षण करने के लिए `view()` ऑपरेटर का उपयोग कैसे करें।

### आगे क्या है?

सीखें कि वर्कफ़्लो को इनपुट वैल्यू के स्रोत के रूप में एक फ़ाइल कैसे लें।

---

## 4. CSV फ़ाइल से इनपुट वैल्यू पढ़ें

यथार्थवादी रूप से, हम शायद ही कभी वैल्यू की एक ऐरे से शुरू करने जा रहे हैं।
सबसे अधिक संभावना है, हमारे पास एक या अधिक फ़ाइलें होंगी जिनमें डेटा होगा जिसे संसाधित करने की आवश्यकता है, किसी प्रकार के संरचित प्रारूप में।

हमने `greetings.csv` नामक एक CSV फ़ाइल तैयार की है जिसमें कई इनपुट अभिवादन हैं, जो उस प्रकार के स्तंभीय डेटा की नकल करते हैं जिसे तुम वास्तविक डेटा विश्लेषण में संसाधित करना चाहते हो, `data/` के तहत संग्रहीत।
(संख्याएं सार्थक नहीं हैं, वे केवल उदाहरणात्मक उद्देश्यों के लिए हैं।)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

हमारा अगला कार्य इस फ़ाइल से वैल्यू पढ़ने के लिए अपने वर्कफ़्लो को अनुकूलित करना है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

चलो देखते हैं कि हम ऐसा कैसे कर सकते हैं।

### 4.1. अभिवादनों के स्रोत के रूप में CSV फ़ाइल की अपेक्षा करने के लिए स्क्रिप्ट को संशोधित करें

शुरू करने के लिए, हमें स्क्रिप्ट में दो प्रमुख परिवर्तन करने होंगे:

- CSV फ़ाइल को इंगित करने के लिए इनपुट पैरामीटर को स्विच करें
- एक फ़ाइल को संभालने के लिए डिज़ाइन किए गए चैनल फैक्ट्री पर स्विच करें

#### 4.1.1. CSV फ़ाइल को इंगित करने के लिए इनपुट पैरामीटर को स्विच करें

याद रखें `params.input` पैरामीटर जिसे हमने पार्ट 1 में सेट अप किया था?
हम इसे हमारे अभिवादनों वाली CSV फ़ाइल को इंगित करने के लिए अपडेट करने जा रहे हैं।

पैरामीटर घोषणा में निम्नलिखित संपादन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline पैरामीटर
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline पैरामीटर
     */
    input: String = 'Holà mundo!'
    ```

यह मानता है कि फ़ाइल वर्कफ़्लो कोड के साथ सह-स्थित है।
तुम बाद में अपनी Nextflow यात्रा में अन्य डेटा स्थानों से निपटने का तरीका सीखोगे।

#### 4.1.2. एक फ़ाइल को संभालने के लिए डिज़ाइन किए गए चैनल फैक्ट्री पर स्विच करें

चूंकि अब हम इनपुट के रूप में सरल स्ट्रिंग के बजाय एक फ़ाइल का उपयोग करना चाहते हैं, हम पहले से `channel.of()` चैनल फैक्ट्री का उपयोग नहीं कर सकते।
हमें एक नए चैनल फैक्ट्री, [`channel.fromPath()`](https://nextflow.io/docs/latest/reference/channel.html#frompath) का उपयोग करने के लिए स्विच करने की आवश्यकता है, जिसमें फ़ाइल पथों को संभालने के लिए कुछ अंतर्निहित कार्यक्षमता है।

वर्कफ़्लो ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "flatten से पहले: $greeting" }
                             // .flatten()
                             // .view { greeting -> "flatten के बाद: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // इनपुट अभिवादनों की एक ऐरे घोषित करें
        greetings_array = ['Hello','Bonjour','Holà']
        // इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "flatten से पहले: $greeting" }
                             .flatten()
                             .view { greeting -> "flatten के बाद: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

तुम देखोगे कि हमने चैनल इनपुट को वापस `param.input` पर स्विच किया है, और `greetings_array` घोषणा को हटा दिया है क्योंकि हमें अब इसकी आवश्यकता नहीं होगी।
हमने `flatten()` और दूसरे `view()` स्टेटमेंट को भी कमेंट आउट कर दिया है।

#### 4.1.3. वर्कफ़्लो चलाएं

चलो नए चैनल फैक्ट्री और इनपुट फ़ाइल के साथ वर्कफ़्लो चलाने की कोशिश करते हैं।

```bash
nextflow run hello-channels.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    flatten से पहले: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

ओह नहीं, यह काम नहीं करता। कंसोल आउटपुट और त्रुटि संदेश की शुरुआत देखें।
`Command executed:` बिट यहां विशेष रूप से सहायक है।

यह थोड़ा परिचित लग सकता है।
ऐसा लगता है कि Nextflow ने फ़ाइल पथ को ही एक स्ट्रिंग वैल्यू के रूप में उपयोग करते हुए एक ही प्रोसेस कॉल चलाने की कोशिश की।
तो इसने फ़ाइल पथ को सही ढंग से हल किया है, लेकिन इसने वास्तव में इसकी सामग्री को पार्स नहीं किया, जो हम चाहते थे।

हम Nextflow को फ़ाइल खोलने और इसकी सामग्री को चैनल में लोड करने के लिए कैसे कहें?

ऐसा लगता है कि हमें एक और [ऑपरेटर](https://nextflow.io/docs/latest/reference/operator.html) की आवश्यकता है!

### 4.2. फ़ाइल को पार्स करने के लिए `splitCsv()` ऑपरेटर का उपयोग करें

ऑपरेटरों की सूची को फिर से देखते हुए, हमें [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) मिलता है, जो CSV-स्वरूपित टेक्स्ट को पार्स और विभाजित करने के लिए डिज़ाइन किया गया है।

#### 4.2.1. चैनल पर `splitCsv()` लागू करें

ऑपरेटर लागू करने के लिए, हम इसे पहले की तरह चैनल फैक्ट्री लाइन में जोड़ते हैं।

वर्कफ़्लो ब्लॉक में, `flatten()` को `splitcsv()` (अनकमेंटेड) से बदलने के लिए निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "splitCsv से पहले: $csv" }
                             .splitCsv()
                             .view { csv -> "splitCsv के बाद: $csv" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "flatten से पहले: $greeting" }
                             // .flatten()
                             // .view { greeting -> "flatten के बाद: $greeting" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

जैसा कि तुम देख सकते हो, हमने पहले/बाद के `view()` स्टेटमेंट को भी अपडेट किया है।
तकनीकी रूप से हम वही वेरिएबल नाम (`greeting`) का उपयोग कर सकते थे लेकिन हमने इसे दूसरों द्वारा कोड को अधिक पठनीय बनाने के लिए कुछ अधिक उपयुक्त (`csv`) में अपडेट किया।

#### 4.2.2. वर्कफ़्लो को फिर से चलाएं

चलो जोड़े गए CSV-पार्सिंग लॉजिक के साथ वर्कफ़्लो चलाने की कोशिश करते हैं।

```bash
nextflow run hello-channels.nf
```

??? failure "कमांड आउटपुट"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    splitCsv से पहले: /workspaces/training/hello-nextflow/data/greetings.csv
    splitCsv के बाद: [Hello, English, 123]
    splitCsv के बाद: [Bonjour, French, 456]
    splitCsv के बाद: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

दिलचस्प बात यह है कि यह भी विफल हो जाता है, लेकिन एक अलग त्रुटि के साथ।
इस बार Nextflow ने फ़ाइल की सामग्री को पार्स किया है (यय!) लेकिन इसने प्रत्येक पंक्ति को एक ऐरे के रूप में लोड किया है, और प्रत्येक ऐरे चैनल में एक तत्व है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-fail.svg"
</figure>

हमें इसे बताना होगा कि प्रत्येक पंक्ति में केवल पहला कॉलम लें।
तो हम इसे कैसे अनपैक करें?

हमने पहले चैनल की सामग्री को अनपैक करने के लिए `flatten()` का उपयोग किया है, लेकिन यह यहां काम नहीं करेगा क्योंकि flatten _सब कुछ_ अनपैक करता है (यदि तुम चाहो तो इसे स्वयं देखने के लिए इसे आज़माने के लिए स्वतंत्र महसूस करो)।

इसके बजाय, हम `map()` नामक एक और ऑपरेटर का उपयोग करेंगे जो वास्तव में उपयोगी है और Nextflow पाइपलाइनों में बहुत कुछ पॉप अप होता है।

### 4.3. अभिवादनों को निकालने के लिए `map()` ऑपरेटर का उपयोग करें

[`map()`](https://nextflow.io/docs/latest/reference/operator.html#map) ऑपरेटर एक बहुत ही आसान छोटा उपकरण है जो हमें चैनल की सामग्री के लिए सभी प्रकार की मैपिंग करने की अनुमति देता है।

इस मामले में, हम इसका उपयोग उस एक तत्व को निकालने के लिए करने जा रहे हैं जो हम अपनी डेटा फ़ाइल में प्रत्येक पंक्ति से चाहते हैं।
यह वह सिंटैक्स है जैसा दिखता है:

```groovy title="सिंटैक्स"
.map { row -> row[0] }
```

इसका मतलब है 'चैनल में प्रत्येक पंक्ति के लिए, 0वां (पहला) आइटम लें जो इसमें है'।

तो चलो इसे अपनी CSV पार्सिंग पर लागू करते हैं।

#### 4.3.1. चैनल पर `map()` लागू करें

वर्कफ़्लो ब्लॉक में, निम्नलिखित कोड परिवर्तन करें:

=== "बाद में"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "splitCsv से पहले: $csv" }
                             .splitCsv()
                             .view { csv -> "splitCsv के बाद: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "map के बाद: $csv" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "पहले"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // CSV फ़ाइल से इनपुट के लिए एक चैनल बनाएं
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "splitCsv से पहले: $csv" }
                             .splitCsv()
                             .view { csv -> "splitCsv के बाद: $csv" }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

तुम देखते हो कि हमने यह पुष्टि करने के लिए एक और `view()` कॉल जोड़ा कि ऑपरेटर वही करता है जो हम उम्मीद करते हैं।

#### 4.3.2. वर्कफ़्लो चलाएं

चलो इसे एक बार और चलाते हैं:

```bash
nextflow run hello-channels.nf
```

??? success "कमांड आउटपुट"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    splitCsv से पहले: /workspaces/training/hello-nextflow/data/greetings.csv
    splitCsv के बाद: [Hello, English, 123]
    splitCsv के बाद: [Bonjour, French, 456]
    splitCsv के बाद: [Holà, Spanish, 789]
    map के बाद: Hello
    map के बाद: Bonjour
    map के बाद: Holà
    ```

इस बार इसे बिना त्रुटि के चलना चाहिए।

`view()` स्टेटमेंट के आउटपुट को देखते हुए, तुम निम्नलिखित देखते हो:

- एक ही `splitCsv से पहले:` स्टेटमेंट: उस बिंदु पर चैनल में एक आइटम होता है, मूल फ़ाइल पथ।
- तीन अलग `splitCsv के बाद:` स्टेटमेंट: प्रत्येक अभिवादन के लिए एक, लेकिन प्रत्येक एक ऐरे के भीतर निहित है जो फ़ाइल में उस पंक्ति से मेल खाता है।
- तीन अलग `map के बाद:` स्टेटमेंट: प्रत्येक अभिवादन के लिए एक, जो अब चैनल में व्यक्तिगत तत्व हैं।

_ध्यान दें कि तुम्हारे आउटपुट में लाइनें एक अलग क्रम में दिखाई दे सकती हैं।_

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-split-and-map.svg"
</figure>

तुम यह सत्यापित करने के लिए आउटपुट फ़ाइलों को भी देख सकते हो कि प्रत्येक अभिवादन को सही ढंग से निकाला गया और वर्कफ़्लो के माध्यम से संसाधित किया गया।

हमने पहले जैसा ही परिणाम प्राप्त किया है, लेकिन अब हमारे पास किसी भी कोड को संशोधित किए बिना एक इनपुट फ़ाइल को संशोधित करके अभिवादनों के चैनल में अधिक तत्व जोड़ने के लिए बहुत अधिक लचीलापन है।
तुम बाद के प्रशिक्षण में जटिल इनपुट को संभालने के लिए अधिक परिष्कृत दृष्टिकोण सीखोगे।

### सारांश

तुम जानते हो कि इनपुट वैल्यू की एक फ़ाइल पढ़ने और उन्हें उचित रूप से संभालने के लिए `.fromPath()` चैनल कंस्ट्रक्टर और ऑपरेटर `splitCsv()` और `map()` का उपयोग कैसे करें।

अधिक सामान्य रूप से, तुम्हें इस बात की बुनियादी समझ है कि Nextflow प्रोसेस के इनपुट को प्रबंधित करने के लिए **channels** का उपयोग कैसे करता है और उनकी सामग्री को रूपांतरित करने के लिए **operators** का उपयोग कैसे करता है।
तुमने यह भी देखा है कि चैनल अंतर्निहित रूप से समानांतर निष्पादन को कैसे संभालते हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-channels-parallel.svg"
</figure>

### आगे क्या है?

एक बड़ा ब्रेक लो, तुमने इसमें कड़ी मेहनत की!

जब तुम तैयार हो, तो अधिक चरण जोड़ने और उन्हें एक उचित वर्कफ़्लो में एक साथ जोड़ने के बारे में जानने के लिए [**पार्ट 3: Hello Workflow**](./03_hello_workflow.md) पर जाएं।

---

## क्विज़

<quiz>
Nextflow में चैनल क्या है?
- [ ] एक फ़ाइल पथ विनिर्देश
- [ ] एक प्रोसेस परिभाषा
- [x] प्रोसेस के बीच डेटा पास करने के लिए एक queue-जैसी संरचना
- [ ] एक कॉन्फ़िगरेशन सेटिंग

अधिक जानें: [1.1. एक इनपुट चैनल बनाएं](#11-एक-इनपुट-चैनल-बनाएं)
</quiz>

<quiz>
यह कोड क्या आउटपुट करेगा?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (एक ही सूची)
- [x] प्रत्येक तत्व एक अलग पंक्ति पर: `Hello`, `Bonjour`, `Hola`
- [ ] कुछ नहीं (चैनल डिफ़ॉल्ट रूप से प्रिंट नहीं होते)
- [ ] एक त्रुटि (अमान्य सिंटैक्स)

अधिक जानें: [1.1. एक इनपुट चैनल बनाएं](#11-एक-इनपुट-चैनल-बनाएं)
</quiz>

<quiz>
जब एक चैनल में कई वैल्यू होती हैं, तो Nextflow प्रोसेस निष्पादन को कैसे संभालता है?
- [ ] प्रोसेस सभी वैल्यू के साथ एक बार चलता है
- [x] प्रोसेस चैनल में प्रत्येक वैल्यू के लिए एक बार चलता है
- [ ] प्रोसेस केवल पहली वैल्यू के साथ चलता है
- [ ] प्रोसेस केवल अंतिम वैल्यू के साथ चलता है

अधिक जानें: [2. कई इनपुट वैल्यू पर चलने के लिए वर्कफ़्लो को संशोधित करें](#2-कई-इनपुट-वैल्यू-पर-चलने-के-लिए-वर्कफ़्लो-को-संशोधित-करें)
</quiz>

<quiz>
`flatten()` ऑपरेटर क्या करता है?
- [ ] कई चैनलों को एक में जोड़ता है
- [ ] चैनल तत्वों को सॉर्ट करता है
- [x] ऐरे को व्यक्तिगत तत्वों में अनपैक करता है
- [ ] डुप्लिकेट तत्वों को हटाता है

अधिक जानें: [3.2.1. `flatten()` ऑपरेटर जोड़ें](#321-flatten-ऑपरेटर-जोड़ें)
</quiz>

<quiz>
`view()` ऑपरेटर का उद्देश्य क्या है?
- [ ] चैनल सामग्री को फ़िल्टर करना
- [ ] चैनल तत्वों को रूपांतरित करना
- [x] चैनल सामग्री का निरीक्षण और डिबग करना
- [ ] चैनल सामग्री को एक फ़ाइल में सहेजना

अधिक जानें: [1.4. चैनल सामग्री का निरीक्षण करने के लिए `view()` का उपयोग करें](#14-चैनल-सामग्री-का-निरीक्षण-करने-के-लिए-view-का-उपयोग-करें)
</quiz>

<quiz>
`splitCsv()` क्या करता है?
- [ ] चैनल सामग्री से एक CSV फ़ाइल बनाता है
- [ ] एक स्ट्रिंग को कॉमा द्वारा विभाजित करता है
- [x] एक CSV फ़ाइल को प्रत्येक पंक्ति का प्रतिनिधित्व करने वाले ऐरे में पार्स करता है
- [ ] कई CSV फ़ाइलों को मर्ज करता है

अधिक जानें: [4.2. फ़ाइल को पार्स करने के लिए `splitCsv()` ऑपरेटर का उपयोग करें](#42-फ़ाइल-को-पार्स-करने-के-लिए-splitcsv-ऑपरेटर-का-उपयोग-करें)
</quiz>

<quiz>
`map()` ऑपरेटर का उद्देश्य क्या है?
- [ ] चैनल से तत्वों को फ़िल्टर करना
- [ ] कई चैनलों को संयोजित करना
- [x] चैनल में प्रत्येक तत्व को रूपांतरित करना
- [ ] चैनल में तत्वों को गिनना

अधिक जानें: [4.3. अभिवादनों को निकालने के लिए `map()` ऑपरेटर का उपयोग करें](#43-अभिवादनों-को-निकालने-के-लिए-map-ऑपरेटर-का-उपयोग-करें)
</quiz>

<quiz>
कई इनपुट को संसाधित करते समय गतिशील आउटपुट फ़ाइलनाम का उपयोग करना क्यों महत्वपूर्ण है?
- [ ] प्रदर्शन में सुधार के लिए
- [ ] डिस्क स्पेस कम करने के लिए
- [x] आउटपुट फ़ाइलों को एक दूसरे को अधिलेखित करने से रोकने के लिए
- [ ] resume कार्यक्षमता को सक्षम करने के लिए

अधिक जानें: [2.2. सुनिश्चित करें कि आउटपुट फ़ाइल नाम अद्वितीय होंगे](#22-सुनिश्चित-करें-कि-आउटपुट-फ़ाइल-नाम-अद्वितीय-होंगे)
</quiz>
