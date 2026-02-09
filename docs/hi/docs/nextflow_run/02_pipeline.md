# भाग 2: वास्तविक पाइपलाइन चलाना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

इस कोर्स के भाग 1 (बुनियादी ऑपरेशन चलाना) में, हमने एक उदाहरण workflow से शुरुआत की जिसमें कोड की जटिलता को कम रखने के लिए केवल न्यूनतम फीचर्स थे।
उदाहरण के लिए, `1-hello.nf` ने एक बार में एक ही वैल्यू प्रदान करने के लिए एक command-line पैरामीटर (`--input`) का उपयोग किया।

हालांकि, अधिकांश वास्तविक दुनिया की पाइपलाइन बड़े पैमाने पर बड़ी मात्रा में डेटा की कुशल प्रोसेसिंग को सक्षम करने के लिए अधिक परिष्कृत फीचर्स का उपयोग करती हैं, और कभी-कभी जटिल लॉजिक द्वारा एक साथ जुड़े कई प्रोसेसिंग स्टेप्स लागू करती हैं।

प्रशिक्षण के इस भाग में, हम मूल Hello World पाइपलाइन के विस्तारित संस्करणों को आज़माकर वास्तविक दुनिया की पाइपलाइनों की प्रमुख विशेषताओं का प्रदर्शन करते हैं।

## 1. फ़ाइल से इनपुट डेटा प्रोसेस करना

वास्तविक दुनिया की पाइपलाइन में, हम आमतौर पर एक या अधिक इनपुट फ़ाइलों में निहित कई डेटा पॉइंट्स (या डेटा सीरीज़) को प्रोसेस करना चाहते हैं।
और जहां भी संभव हो, हम विश्लेषण की प्रतीक्षा में बिताए गए समय को कम करने के लिए स्वतंत्र डेटा की प्रोसेसिंग को समानांतर में चलाना चाहते हैं।

यह प्रदर्शित करने के लिए कि Nextflow यह कैसे करता है, हमने `greetings.csv` नामक एक CSV फ़ाइल तैयार की है जिसमें कई इनपुट अभिवादन हैं, जो उस प्रकार के स्तंभीय डेटा की नकल करते हैं जिसे तुम वास्तविक डेटा विश्लेषण में प्रोसेस करना चाहोगे।
ध्यान दो कि संख्याएं सार्थक नहीं हैं, वे केवल उदाहरण के उद्देश्यों के लिए हैं।

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

हमने मूल workflow का एक बेहतर संस्करण भी लिखा है, जिसे अब `2a-inputs.nf` कहा जाता है, जो CSV फ़ाइल को पढ़ेगा, अभिवादनों को निकालेगा और उनमें से प्रत्येक को एक अलग फ़ाइल में लिखेगा।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

चलो पहले workflow चलाते हैं, और फिर हम प्रासंगिक Nextflow कोड पर एक नज़र डालेंगे।

### 1.1. Workflow चलाना

अपने टर्मिनल में निम्नलिखित कमांड चलाओ।

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

रोमांचक रूप से, यह संकेत देता है कि प्रोसेस के लिए '3 of 3' कॉल किए गए, जो उत्साहजनक है, क्योंकि हमने इनपुट के रूप में प्रदान की गई CSV में डेटा की तीन पंक्तियां थीं।
यह सुझाव देता है कि `sayHello()` प्रोसेस को तीन बार कॉल किया गया, प्रत्येक इनपुट पंक्ति पर एक बार।

### 1.2. `results` डायरेक्टरी में प्रकाशित आउटपुट खोजना

चलो 'results' डायरेक्टरी को देखते हैं कि क्या हमारा workflow अभी भी वहां हमारे आउटपुट की एक कॉपी लिख रहा है।

??? abstract "डायरेक्टरी सामग्री"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

हां! हम `2a-inputs` नामक एक नई डायरेक्टरी देखते हैं जिसमें विभिन्न नामों वाली तीन आउटपुट फ़ाइलें हैं, काफी सुविधाजनक रूप से।

तुम उनमें से प्रत्येक को खोल सकते हो यह सुनिश्चित करने के लिए कि उनमें उपयुक्त अभिवादन स्ट्रिंग है।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

यह पुष्टि करता है कि इनपुट फ़ाइल में प्रत्येक अभिवादन को उचित रूप से प्रोसेस किया गया है।

### 1.3. मूल आउटपुट और लॉग खोजना

तुमने देखा होगा कि ऊपर दिए गए console आउटपुट में केवल एक task डायरेक्टरी का उल्लेख था।
क्या इसका मतलब यह है कि `sayHello()` के सभी तीन कॉल उस एक task डायरेक्टरी के भीतर निष्पादित किए गए थे?

#### 1.3.1. टर्मिनल में दी गई task डायरेक्टरी की जांच करना

चलो उस `8e/0eb066` task डायरेक्टरी के अंदर देखते हैं।

??? abstract "डायरेक्टरी सामग्री"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

हमें केवल एक अभिवादन के अनुरूप आउटपुट मिलता है (साथ ही सहायक फ़ाइलें यदि हम छिपी हुई फ़ाइलों का प्रदर्शन सक्षम करते हैं)।

तो यहां क्या हो रहा है?

डिफ़ॉल्ट रूप से, ANSI लॉगिंग सिस्टम एक ही लाइन पर एक ही प्रोसेस के सभी कॉल के लिए स्टेटस जानकारी लिखता है।
परिणामस्वरूप, इसने हमें console आउटपुट में केवल तीन task डायरेक्टरी पथों (`8e/0eb066`) में से एक दिखाया।
दो अन्य हैं जो वहां सूचीबद्ध नहीं हैं।

#### 1.3.2. टर्मिनल को अधिक विवरण दिखाना

हम प्रोसेस कॉल की पूरी सूची देखने के लिए लॉगिंग व्यवहार को संशोधित कर सकते हैं, कमांड में `-ansi-log false` जोड़कर निम्नानुसार:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

इस बार हम सभी तीन प्रोसेस रन और उनसे जुड़ी work subdirectories को आउटपुट में सूचीबद्ध देखते हैं।
ANSI लॉगिंग को अक्षम करने से Nextflow को टर्मिनल आउटपुट में रंगों का उपयोग करने से भी रोका गया।

ध्यान दो कि दोनों लॉगिंग मोड के बीच स्टेटस रिपोर्ट करने का तरीका थोड़ा अलग है।
संघनित मोड में, Nextflow रिपोर्ट करता है कि कॉल सफलतापूर्वक पूर्ण हुए या नहीं।
इस विस्तारित मोड में, यह केवल रिपोर्ट करता है कि वे सबमिट किए गए थे।

यह पुष्टि करता है कि `sayHello()` प्रोसेस को तीन बार कॉल किया जाता है, और प्रत्येक के लिए एक अलग task डायरेक्टरी बनाई जाती है।

यदि हम वहां सूचीबद्ध प्रत्येक task डायरेक्टरी के अंदर देखें, तो हम सत्यापित कर सकते हैं कि प्रत्येक एक अभिवादन से मेल खाती है।

??? abstract "डायरेक्टरी सामग्री"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

यह पुष्टि करता है कि प्रत्येक प्रोसेस कॉल अन्य सभी से अलगाव में निष्पादित किया जाता है।
इसके कई फायदे हैं, जिसमें टकराव से बचना शामिल है यदि प्रोसेस गैर-अद्वितीय नामों के साथ कोई मध्यवर्ती फ़ाइलें उत्पन्न करती है।

!!! tip

    एक जटिल workflow के लिए, या बड़ी संख्या में इनपुट के लिए, टर्मिनल पर पूरी सूची का आउटपुट थोड़ा भारी हो सकता है, इसलिए लोग आमतौर पर नियमित उपयोग में `-ansi-log false` का उपयोग नहीं करते हैं।

### 1.4. Workflow कोड की जांच करना

तो workflow का यह संस्करण इनपुट की एक CSV फ़ाइल पढ़ने, इनपुट को अलग से प्रोसेस करने और आउटपुट को विशिष्ट रूप से नामकरण करने में सक्षम है।

चलो देखते हैं कि workflow कोड में यह क्या संभव बनाता है।

??? full-code "पूर्ण कोड फ़ाइल"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

एक बार फिर, तुम्हें कोड सिंटैक्स याद रखने की आवश्यकता नहीं है, लेकिन महत्वपूर्ण कार्यक्षमता प्रदान करने वाले workflow के प्रमुख घटकों को पहचानना सीखना अच्छा है।

#### 1.4.1. CSV से इनपुट डेटा लोड करना

यह सबसे दिलचस्प हिस्सा है: हम command-line से एक ही वैल्यू लेने से CSV फ़ाइल लेने, इसे पार्स करने और इसमें निहित व्यक्तिगत अभिवादनों को प्रोसेस करने में कैसे स्विच हुए?

Nextflow में, हम इसे एक [**channel**](https://nextflow.io/docs/latest/channel.html) के साथ करते हैं: एक queue construct जो इनपुट को कुशलता से संभालने और बहु-चरणीय workflows में एक चरण से दूसरे चरण में शटल करने के लिए डिज़ाइन किया गया है, जबकि अंतर्निहित समानांतरता और कई अतिरिक्त लाभ प्रदान करता है।

चलो इसे तोड़ते हैं।

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // एक अभिवादन emit करें
    sayHello(greeting_ch)
```

यह कोड `greeting_ch` नामक एक channel बनाता है जो CSV फ़ाइल को पढ़ता है, इसे पार्स करता है, और प्रत्येक पंक्ति से पहले कॉलम को निकालता है।
परिणाम `Hello`, `Bonjour`, और `Holà` युक्त एक channel है।

??? tip "यह कैसे काम करता है?"

    यहां बताया गया है कि उस लाइन का सादे अंग्रेजी में क्या अर्थ है:

    - `channel.fromPath` एक **channel factory** है जो फ़ाइल पथ(ओं) से एक channel बनाता है
    - `(params.input)` निर्दिष्ट करता है कि filepath command line पर `--input` द्वारा प्रदान किया गया है

    दूसरे शब्दों में, वह लाइन Nextflow को बताती है: `--input` के साथ दिए गए filepath को लो और इसकी सामग्री को इनपुट डेटा के रूप में मानने के लिए तैयार हो जाओ।

    फिर अगली दो लाइनें **operators** लागू करती हैं जो फ़ाइल की वास्तविक पार्सिंग और उपयुक्त डेटा संरचना में डेटा की लोडिंग करती हैं:

    - `.splitCsv()` Nextflow को CSV फ़ाइल को पंक्तियों और कॉलमों का प्रतिनिधित्व करने वाले एक array में पार्स करने के लिए कहता है
    - `.map { line -> line[0] }` Nextflow को प्रत्येक पंक्ति से केवल पहले कॉलम में तत्व लेने के लिए कहता है

    तो व्यवहार में, निम्नलिखित CSV फ़ाइल से शुरू करते हुए:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    हमने इसे एक array में बदल दिया है जो इस तरह दिखता है:

    ```txt title="Array सामग्री"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    और फिर हमने तीनों पंक्तियों में से प्रत्येक से पहला तत्व लिया है और उन्हें एक Nextflow channel में लोड किया है जिसमें अब शामिल है: `Hello`, `Bonjour`, और `Holà`।

    यदि तुम channels और operators को गहराई से समझना चाहते हो, जिसमें उन्हें स्वयं कैसे लिखना है, देखो [Hello Nextflow Part 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file)।

#### 1.4.2. प्रत्येक अभिवादन पर प्रोसेस को कॉल करना

इसके बाद, workflow के `main:` ब्लॉक की अंतिम पंक्ति में, हम लोड किए गए `greeting_ch` channel को `sayHello()` प्रोसेस के इनपुट के रूप में प्रदान करते हैं।

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // एक अभिवादन emit करें
    sayHello(greeting_ch)
```

यह Nextflow को channel में प्रत्येक तत्व पर व्यक्तिगत रूप से प्रोसेस चलाने के लिए कहता है, _यानी_ प्रत्येक अभिवादन पर।
और क्योंकि Nextflow इस तरह स्मार्ट है, यह उपलब्ध कंप्यूटिंग इंफ्रास्ट्रक्चर के आधार पर, यदि संभव हो तो इन प्रोसेस कॉल को समानांतर में चलाएगा।

यही वह तरीका है जिससे तुम तुलनात्मक रूप से बहुत कम कोड के साथ बहुत सारे डेटा (कई नमूने, या डेटा पॉइंट्स, जो भी तुम्हारी अनुसंधान की इकाई है) की कुशल और स्केलेबल प्रोसेसिंग प्राप्त कर सकते हो।

#### 1.4.3. आउटपुट कैसे नामित किए जाते हैं

अंत में, यह देखने लायक है कि हम आउटपुट फ़ाइलों को विशिष्ट रूप से नामित कैसे करते हैं, प्रोसेस कोड पर एक त्वरित नज़र डालना।

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
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

तुम देखते हो कि, `1-hello.nf` में इस प्रोसेस के संस्करण की तुलना में, आउटपुट घोषणा और कमांड का प्रासंगिक हिस्सा आउटपुट फ़ाइल नाम में अभिवादन वैल्यू शामिल करने के लिए बदल गया है।

यह सुनिश्चित करने का एक तरीका है कि आउटपुट फ़ाइल नाम टकराएंगे नहीं जब वे सामान्य results डायरेक्टरी में प्रकाशित हों।

और यह एकमात्र परिवर्तन है जो हमें प्रोसेस घोषणा के अंदर करना पड़ा है!

### सारांश

तुम बुनियादी स्तर पर समझते हो कि channels और operators हमें कई इनपुट को कुशलता से प्रोसेस करने में कैसे सक्षम बनाते हैं।

### आगे क्या है?

जानो कि बहु-चरणीय workflows कैसे बनाए जाते हैं और वे कैसे संचालित होते हैं।

---

## 2. बहु-चरणीय workflows चलाना

अधिकांश वास्तविक दुनिया की workflows में एक से अधिक चरण शामिल होते हैं।
चलो channels के बारे में हमने जो सीखा है उस पर निर्माण करते हैं, और देखते हैं कि Nextflow बहु-चरणीय workflow में प्रोसेस को एक साथ जोड़ने के लिए channels और operators का उपयोग कैसे करता है।

इसके लिए, हम तुम्हें एक उदाहरण workflow प्रदान करते हैं जो तीन अलग-अलग चरणों को एक साथ जोड़ता है और निम्नलिखित प्रदर्शित करता है:

1. एक प्रोसेस से दूसरे में डेटा प्रवाह बनाना
2. कई प्रोसेस कॉल से आउटपुट को एक ही प्रोसेस कॉल में एकत्र करना

विशेष रूप से, हमने workflow का एक विस्तारित संस्करण बनाया जिसे `2b-multistep.nf` कहा जाता है जो प्रत्येक इनपुट अभिवादन लेता है, इसे uppercase में परिवर्तित करता है, फिर सभी uppercased अभिवादनों को एक ही आउटपुट फ़ाइल में एकत्र करता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

पहले की तरह, हम पहले workflow चलाएंगे फिर कोड को देखेंगे कि क्या नया है।

### 2.1. Workflow चलाना

अपने टर्मिनल में निम्नलिखित कमांड चलाओ:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

तुम देखते हो कि जैसा वादा किया गया था, workflow के हिस्से के रूप में कई चरण चलाए गए; पहले दो (`sayHello` और `convertToUpper`) संभवतः प्रत्येक व्यक्तिगत अभिवादन पर चलाए गए, और तीसरा (`collectGreetings`) केवल एक बार चलाया गया होगा, सभी तीन `convertToUpper` कॉल के आउटपुट पर।

### 2.2. आउटपुट खोजना

चलो `results` डायरेक्टरी में देखकर सत्यापित करते हैं कि वास्तव में ऐसा ही हुआ।

??? abstract "डायरेक्टरी सामग्री"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

जैसा कि तुम देख सकते हो, हमारे पास `2b-multistep` नामक एक नई डायरेक्टरी है, और इसमें पहले की तुलना में काफी अधिक फ़ाइलें हैं।
कुछ फ़ाइलों को `intermediates` नामक एक subdirectory में समूहीकृत किया गया है, जबकि दो फ़ाइलें शीर्ष स्तर पर स्थित हैं।

वे दो बहु-चरणीय workflow के अंतिम परिणाम हैं।
फ़ाइल नामों को देखने और उनकी सामग्री की जांच करने के लिए एक मिनट लो यह पुष्टि करने के लिए कि वे वही हैं जो तुम उम्मीद करते हो।

??? abstract "फ़ाइल सामग्री"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

पहले में हमारे तीन अभिवादन हैं, uppercased और वादे के अनुसार एक ही फ़ाइल में वापस एकत्र किए गए।
दूसरा एक रिपोर्ट फ़ाइल है जो रन के बारे में कुछ जानकारी को सारांशित करती है।

### 2.3. कोड की जांच करना

चलो कोड को देखते हैं और बहु-चरणीय workflows के लिए प्रमुख पैटर्न की पहचान करते हैं।

??? full-code "पूर्ण कोड फ़ाइल"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

वहां बहुत कुछ चल रहा है, लेकिन workflow के पिछले संस्करण की तुलना में सबसे स्पष्ट अंतर यह है कि अब कई प्रोसेस परिभाषाएं हैं, और तदनुसार, workflow ब्लॉक में कई प्रोसेस कॉल हैं।

चलो करीब से देखते हैं और देखते हैं कि क्या हम सबसे दिलचस्प टुकड़ों की पहचान कर सकते हैं।

#### 2.3.1. Workflow संरचना को विज़ुअलाइज़ करना

यदि तुम Nextflow एक्सटेंशन के साथ VSCode का उपयोग कर रहे हो, तो तुम किसी भी Nextflow स्क्रिप्ट में workflow ब्लॉक के ठीक ऊपर प्रदर्शित छोटे `DAG preview` लिंक पर क्लिक करके प्रोसेस कैसे जुड़े हैं इसका एक सहायक आरेख प्राप्त कर सकते हो।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

यह तुम्हें एक अच्छा अवलोकन देता है कि प्रोसेस कैसे जुड़े हैं और वे क्या उत्पन्न करते हैं।

तुम देखते हो कि मूल `sayHello` प्रोसेस के अलावा, हमारे पास अब `convertToUpper` और `collectGreetings` भी हैं, जो उन प्रोसेस के नामों से मेल खाते हैं जो हमने console आउटपुट में देखे थे।
दो नई प्रोसेस परिभाषाएं `sayHello` प्रोसेस के समान तरीके से संरचित हैं, सिवाय इसके कि `collectGreetings` `batch` नामक एक अतिरिक्त इनपुट पैरामीटर लेता है और दो आउटपुट उत्पन्न करता है।

हम प्रत्येक के लिए कोड में विस्तार से नहीं जाएंगे, लेकिन यदि तुम उत्सुक हो, तो तुम [Hello Nextflow के भाग 2](../hello_nextflow/03_hello_workflow.md) में विवरण देख सकते हो।

अभी के लिए, चलो देखते हैं कि प्रोसेस एक दूसरे से कैसे जुड़े हैं।

#### 2.3.2. प्रोसेस कैसे जुड़े हैं

यहां देखने के लिए वास्तव में दिलचस्प बात यह है कि workflow के `main:` ब्लॉक में प्रोसेस कॉल कैसे एक साथ जुड़े हैं।

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // एक अभिवादन emit करें
    sayHello(greeting_ch)
    // अभिवादन को uppercase में बदलें
    convertToUpper(sayHello.out)
    // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

तुम देख सकते हो कि पहला प्रोसेस कॉल, `sayHello(greeting_ch)`, अपरिवर्तित है।
फिर अगला प्रोसेस कॉल, `convertToUpper` के लिए, `sayHello` के आउटपुट को `sayHello.out` के रूप में संदर्भित करता है।

पैटर्न सरल है: `processName.out` एक प्रोसेस के आउटपुट channel को संदर्भित करता है, जिसे सीधे अगले प्रोसेस को पास किया जा सकता है।
यही वह तरीका है जिससे हम Nextflow में एक चरण से दूसरे चरण में डेटा शटल करते हैं।

#### 2.3.3. एक प्रोसेस कई इनपुट ले सकता है

तीसरा प्रोसेस कॉल, `collectGreetings` के लिए, थोड़ा अलग है।

```groovy title="2b-multistep.nf" linenums="77"
    // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

तुम देखते हो कि इस कॉल को दो इनपुट दिए गए हैं, `convertToUpper.out.collect()` और `params.batch`।
अभी के लिए `.collect()` बिट को अनदेखा करते हुए, हम इसे `collectGreetings(input1, input2)` के रूप में सामान्यीकृत कर सकते हैं।

यह प्रोसेस मॉड्यूल में दो इनपुट घोषणाओं से मेल खाता है:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

जब Nextflow इसे पार्स करता है, तो यह कॉल में पहले इनपुट को `path input_files` को असाइन करेगा, और दूसरे को `val batch_name` को।

तो अब तुम जानते हो कि एक प्रोसेस कई इनपुट ले सकता है, और workflow ब्लॉक में कॉल कैसा दिखता है।

अब चलो उस पहले इनपुट, `convertToUpper.out.collect()` पर करीब से नज़र डालते हैं।

#### 2.3.4. `collectGreetings` कॉल में `collect()` क्या करता है

`sayHello` के आउटपुट को `convertToUpper` में पास करने के लिए, हमने बस `sayHello` के आउटपुट channel को `sayHello.out` के रूप में संदर्भित किया। लेकिन अगले चरण के लिए, हम `convertToUpper.out.collect()` का संदर्भ देख रहे हैं।

यह `collect()` बिट क्या है और यह क्या करता है?

यह एक operator है, निश्चित रूप से। बिल्कुल `splitCsv` और `map` operators की तरह जिनका हमने पहले सामना किया।
इस बार operator को `collect` कहा जाता है, और `convertToUpper` द्वारा उत्पादित आउटपुट channel पर लागू किया जाता है।

`collect` operator का उपयोग एक ही प्रोसेस के कई कॉल से आउटपुट एकत्र करने और उन्हें एक ही channel तत्व में पैकेज करने के लिए किया जाता है।

इस workflow के संदर्भ में, यह `convertToUpper.out` channel में तीन uppercased अभिवादनों को ले रहा है (जो तीन अलग-अलग channel आइटम हैं, और आमतौर पर अगले प्रोसेस द्वारा अलग-अलग कॉल में संभाले जाएंगे) और उन्हें एक ही आइटम में पैकेज कर रहा है।
यही वह तरीका है जिससे हम सभी अभिवादनों को एक ही फ़ाइल में वापस लाते हैं।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

इसके विपरीत, यदि हम `collectGreetings()` को खिलाने से पहले `convertToUpper()` के आउटपुट पर `collect()` लागू नहीं करते, तो Nextflow बस प्रत्येक अभिवादन पर स्वतंत्र रूप से `collectGreetings()` चलाएगा, जो हमारे लक्ष्य को प्राप्त नहीं करेगा।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

प्रोसेस कॉल के बीच channels की सामग्री में परिवर्तन लागू करने के लिए कई अन्य [operators](https://nextflow.io/docs/latest/reference/operator.html) उपलब्ध हैं।

यह पाइपलाइन डेवलपर्स को उनकी पाइपलाइन के flow logic को अनुकूलित करने के लिए बहुत लचीलापन देता है।
नकारात्मक पक्ष यह है कि यह कभी-कभी यह समझना कठिन बना सकता है कि पाइपलाइन क्या कर रही है।

#### 2.3.5. एक इनपुट पैरामीटर का डिफ़ॉल्ट मान हो सकता है

तुमने देखा होगा कि `collectGreetings` एक दूसरा इनपुट लेता है, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

यह workflow को `--batch` नामक एक CLI पैरामीटर पास करता है।
हालांकि, जब हमने पहले workflow लॉन्च किया, तो हमने `--batch` पैरामीटर निर्दिष्ट नहीं किया।

वहां क्या हो रहा है?
`params` ब्लॉक पर एक नज़र डालो:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

workflow में एक डिफ़ॉल्ट मान कॉन्फ़िगर किया गया है, इसलिए हमें इसे प्रदान करने की आवश्यकता नहीं है।
लेकिन यदि हम command line पर एक प्रदान करते हैं, तो हम जो मान निर्दिष्ट करते हैं वह डिफ़ॉल्ट के बजाय उपयोग किया जाएगा।

इसे आज़माओ:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "कमांड आउटपुट"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

तुम्हें अपने कस्टम batch नाम के साथ नामित नए अंतिम आउटपुट दिखाई देने चाहिए।

??? abstract "डायरेक्टरी सामग्री"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

यह इनपुट कॉन्फ़िगरेशन का एक पहलू है, जिसे हम भाग 3 में अधिक विस्तार से कवर करेंगे, लेकिन अभी के लिए महत्वपूर्ण बात यह जानना है कि इनपुट पैरामीटर को डिफ़ॉल्ट मान दिए जा सकते हैं।

#### 2.3.6. एक प्रोसेस कई आउटपुट उत्पन्न कर सकता है

`collectGreetings` प्रोसेस परिभाषा में, हम निम्नलिखित आउटपुट घोषणाएं देखते हैं:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

जिन्हें फिर `publish:` ब्लॉक में `emit:` के साथ दिए गए नाम से संदर्भित किया जाता है:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

यह workflow में विभिन्न operators के संयोजन में विशिष्ट आउटपुट को व्यक्तिगत रूप से अन्य प्रोसेस में पास करना आसान बनाता है।

#### 2.3.7. प्रकाशित आउटपुट को व्यवस्थित किया जा सकता है

`output` ब्लॉक में, हमने workflow के केवल अंतिम आउटपुट को चुनना आसान बनाने के लिए मध्यवर्ती परिणामों को समूहित करने के लिए कस्टम पथों का उपयोग किया है।

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

प्रकाशित आउटपुट को व्यवस्थित करने के अधिक परिष्कृत तरीके हैं; हम कॉन्फ़िगरेशन पर भाग में कुछ को छूएंगे।

!!! tip "Workflows बनाने के बारे में अधिक जानना चाहते हो?"

    बहु-चरणीय workflows बनाने के विस्तृत कवरेज के लिए, देखो [Hello Nextflow Part 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md)।

### सारांश

तुम बुनियादी स्तर पर समझते हो कि channels और operators का उपयोग करके बहु-चरणीय workflows कैसे बनाए जाते हैं और वे कैसे संचालित होते हैं।
तुमने यह भी देखा है कि प्रोसेस कई इनपुट ले सकते हैं और कई आउटपुट उत्पन्न कर सकते हैं, और इन्हें संरचित तरीके से प्रकाशित किया जा सकता है।

### आगे क्या है?

जानो कि कोड पुन: उपयोग और रखरखाव को बढ़ावा देने के लिए Nextflow पाइपलाइनों को कैसे मॉड्यूलराइज़ किया जा सकता है।

---

## 3. मॉड्यूलराइज़्ड पाइपलाइन चलाना

अब तक, हमने जो सभी workflows देखे हैं वे एक ही workflow फ़ाइल से मिलकर बने हैं जिसमें सभी प्रासंगिक कोड हैं।

हालांकि, वास्तविक दुनिया की पाइपलाइन आमतौर पर _मॉड्यूलराइज़्ड_ होने से लाभान्वित होती हैं, जिसका अर्थ है कि कोड को विभिन्न फ़ाइलों में विभाजित किया जाता है।
यह उनके विकास और रखरखाव को अधिक कुशल और टिकाऊ बना सकता है।

यहां हम Nextflow में कोड मॉड्यूलरिटी के सबसे सामान्य रूप का प्रदर्शन करने जा रहे हैं, जो **modules** का उपयोग है।

Nextflow में, एक [**module**](https://nextflow.io/docs/latest/module.html) एक एकल प्रोसेस परिभाषा है जो एक स्टैंडअलोन कोड फ़ाइल में स्वयं द्वारा समाहित है।
एक workflow में एक module का उपयोग करने के लिए, तुम बस अपनी workflow कोड फ़ाइल में एक single-line import statement जोड़ते हो; फिर तुम प्रोसेस को workflow में उसी तरह एकीकृत कर सकते हो जैसे तुम सामान्य रूप से करोगे।
यह कोड की कई प्रतियां बनाए बिना कई workflows में प्रोसेस परिभाषाओं का पुन: उपयोग करना संभव बनाता है।

अब तक हम workflows चला रहे थे जिनमें उनकी सभी प्रोसेस एक मोनोलिथिक कोड फ़ाइल में शामिल थीं।
अब हम देखने जा रहे हैं कि यह कैसा दिखता है जब प्रोसेस व्यक्तिगत modules में संग्रहीत होते हैं।

हमने निश्चित रूप से एक बार फिर प्रदर्शन उद्देश्यों के लिए एक उपयुक्त workflow तैयार किया है, जिसे `2c-modules.nf` कहा जाता है, साथ ही `modules/` डायरेक्टरी में स्थित modules का एक सेट।

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "डायरेक्टरी सामग्री"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

तुम देखते हो कि चार Nextflow फ़ाइलें हैं, प्रत्येक एक प्रोसेस के नाम पर। 
अभी के लिए `cowpy.nf` फ़ाइल को अनदेखा कर सकते हो; हम उस पर बाद में आएंगे।

### 3.1. कोड की जांच करना

इस बार हम पहले कोड को देखने जा रहे हैं।
`2c-modules.nf` workflow फ़ाइल खोलकर शुरू करो।

??? full-code "पूर्ण कोड फ़ाइल"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // modules शामिल करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

तुम देखते हो कि workflow लॉजिक workflow के पिछले संस्करण के बिल्कुल समान है।
हालांकि, प्रोसेस कोड workflow फ़ाइल से गायब है, और इसके बजाय `modules` के तहत अलग-अलग फ़ाइलों की ओर इशारा करते हुए `include` statements हैं।

```groovy title="hello-modules.nf" linenums="3"
// modules शामिल करें
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

उन फ़ाइलों में से एक को खोलो और तुम्हें संबंधित प्रोसेस के लिए कोड मिलेगा।

??? full-code "पूर्ण कोड फ़ाइल"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
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

जैसा कि तुम देख सकते हो, प्रोसेस कोड नहीं बदला है; इसे बस मुख्य workflow फ़ाइल में होने के बजाय एक व्यक्तिगत module फ़ाइल में कॉपी किया गया है।
यही बात अन्य दो प्रोसेस पर भी लागू होती है।

तो चलो देखते हैं कि इस नए संस्करण को चलाना कैसा दिखता है।

### 3.2. Workflow चलाना

अपने टर्मिनल में यह कमांड चलाओ, `-resume` फ्लैग के साथ:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

तुम देखोगे कि प्रोसेस निष्पादन सभी सफलतापूर्वक cached हुए, जिसका अर्थ है कि Nextflow ने पहचाना कि उसने पहले ही अनुरोधित कार्य कर लिया है, भले ही कोड को विभाजित किया गया हो और मुख्य workflow फ़ाइल का नाम बदल दिया गया हो।

इनमें से कुछ भी Nextflow के लिए मायने नहीं रखता; जो मायने रखता है वह job script है जो सभी कोड को एक साथ खींचने और मूल्यांकन करने के बाद उत्पन्न होती है।

!!! tip

    एक workflow के एक खंड को 'subworkflow' के रूप में समाहित करना भी संभव है जिसे एक बड़ी पाइपलाइन में आयात किया जा सकता है, लेकिन यह इस कोर्स के दायरे से बाहर है।

    तुम [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/) पर Side Quest में composable workflows विकसित करने के बारे में अधिक जान सकते हो।

### सारांश

तुम जानते हो कि कोड पुन: उपयोग को बढ़ावा देने और रखरखाव में सुधार के लिए प्रोसेस को स्टैंडअलोन modules में कैसे संग्रहीत किया जा सकता है।

### आगे क्या है?

सॉफ़्टवेयर निर्भरताओं के प्रबंधन के लिए containers का उपयोग करना सीखो।

---

## 4. Containerized सॉफ़्टवेयर का उपयोग करना

अब तक हम जो workflows उदाहरण के रूप में उपयोग कर रहे थे, उन्हें बस हमारे वातावरण में उपलब्ध UNIX टूल का उपयोग करके बहुत बुनियादी टेक्स्ट प्रोसेसिंग ऑपरेशन चलाने की आवश्यकता थी।

हालांकि, वास्तविक दुनिया की पाइपलाइनों को आमतौर पर विशेष टूल और पैकेज की आवश्यकता होती है जो अधिकांश वातावरणों में डिफ़ॉल्ट रूप से शामिल नहीं होते हैं।
आमतौर पर, तुम्हें इन टूल को इंस्टॉल करने, उनकी निर्भरताओं को प्रबंधित करने और किसी भी संघर्ष को हल करने की आवश्यकता होगी।

यह सब बहुत थकाऊ और कष्टप्रद है।
इस समस्या को हल करने का एक बहुत बेहतर तरीका **containers** का उपयोग करना है।

एक **container** एक हल्की, स्टैंडअलोन, निष्पादन योग्य सॉफ़्टवेयर इकाई है जो एक container **image** से बनाई गई है जिसमें कोड, सिस्टम लाइब्रेरी और सेटिंग्स सहित एक एप्लिकेशन चलाने के लिए आवश्यक सब कुछ शामिल है।

!!! Tip

    हम इसे [Docker](https://www.docker.com/get-started/) तकनीक का उपयोग करके सिखाते हैं, लेकिन Nextflow कई अन्य container तकनीकों का भी समर्थन करता है।
    तुम Nextflow के containers के लिए समर्थन के बारे में [यहां](https://nextflow.io/docs/latest/container.html) अधिक जान सकते हो।

### 4.1. सीधे एक container का उपयोग करना

पहले, चलो सीधे एक container के साथ इंटरैक्ट करने की कोशिश करते हैं।
यह Nextflow में उनका उपयोग करना शुरू करने से पहले containers क्या हैं इसकी तुम्हारी समझ को मजबूत करने में मदद करेगा।

#### 4.1.1. Container image को pull करना

एक container का उपयोग करने के लिए, तुम आमतौर पर एक container registry से एक container image डाउनलोड या "pull" करते हो, और फिर एक container instance बनाने के लिए container image चलाते हो।

सामान्य सिंटैक्स इस प्रकार है:

```bash title="सिंटैक्स"
docker pull '<container>'
```

- `docker pull` container सिस्टम को एक repository से एक container image pull करने का निर्देश है।
- `'<container>'` container image का URI पता है।

एक उदाहरण के रूप में, चलो एक container image pull करते हैं जिसमें [cowpy](https://github.com/jeffbuttars/cowpy) है, एक टूल का python implementation जिसे `cowsay` कहा जाता है जो मज़ेदार तरीके से मनमाने टेक्स्ट इनपुट प्रदर्शित करने के लिए ASCII art उत्पन्न करता है।

विभिन्न repositories हैं जहां तुम प्रकाशित containers पा सकते हो।
हमने `cowpy` Conda पैकेज से इस Docker container image को उत्पन्न करने के लिए [Seqera Containers](https://seqera.io/containers/) सेवा का उपयोग किया: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`।

पूर्ण pull कमांड चलाओ:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "कमांड आउटपुट"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

यह सिस्टम को निर्दिष्ट image डाउनलोड करने के लिए कहता है।
एक बार डाउनलोड पूरा हो जाने के बाद, तुम्हारे पास container image की एक local कॉपी है।

#### 4.1.2. Container को spin up करना

Containers को एक one-off कमांड के रूप में चलाया जा सकता है, लेकिन तुम उन्हें interactively भी उपयोग कर सकते हो, जो तुम्हें container के अंदर एक shell prompt देता है और तुम्हें कमांड के साथ खेलने की अनुमति देता है।

सामान्य सिंटैक्स इस प्रकार है:

```bash title="सिंटैक्स"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` container सिस्टम को एक container image से एक container instance spin up करने और उसमें एक कमांड निष्पादित करने का निर्देश है।
- `--rm` सिस्टम को कमांड पूर्ण होने के बाद container instance को shut down करने के लिए कहता है।

पूरी तरह से असेंबल किया गया, container execution कमांड इस तरह दिखता है:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

वह कमांड चलाओ, और तुम्हें अपना prompt `(base) root@b645838b3314:/tmp#` जैसा कुछ बदलता हुआ दिखाई देना चाहिए, जो इंगित करता है कि तुम अब container के अंदर हो।

तुम directory सामग्री को सूचीबद्ध करने के लिए `ls` चलाकर इसे सत्यापित कर सकते हो:

```bash
ls /
```

??? success "कमांड आउटपुट"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

तुम देखते हो कि container के अंदर filesystem तुम्हारे host सिस्टम पर filesystem से अलग है।

!!! Tip

    जब तुम एक container चलाते हो, तो यह डिफ़ॉल्ट रूप से host सिस्टम से अलग होता है।
    इसका मतलब है कि container host सिस्टम पर किसी भी फ़ाइल तक पहुंच नहीं सकता जब तक कि तुम स्पष्ट रूप से इसे निम्नलिखित सिंटैक्स का उपयोग करके `docker run` कमांड के हिस्से के रूप में एक volume mount करने की अनुमति नहीं देते:

    ```bash title="सिंटैक्स"
    -v <outside_path>:<inside_path>
    ```

    यह प्रभावी रूप से container की दीवार के माध्यम से एक सुरंग स्थापित करता है जिसका उपयोग तुम अपने filesystem के उस हिस्से तक पहुंचने के लिए कर सकते हो।

    यह [Hello Nextflow के भाग 5](../hello_nextflow/05_hello_containers.md) में अधिक विस्तार से कवर किया गया है।

#### 4.1.3. `cowpy` टूल चलाना

Container के अंदर से, तुम `cowpy` कमांड सीधे चला सकते हो।

```bash
cowpy "Hello Containers"
```

??? success "कमांड आउटपुट"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

यह हमारे द्वारा निर्दिष्ट टेक्स्ट युक्त एक speech bubble के साथ डिफ़ॉल्ट cow character (या 'cowacter') का ASCII art उत्पन्न करता है।

अब जब तुमने बुनियादी उपयोग का परीक्षण कर लिया है, तो तुम इसे कुछ पैरामीटर देने की कोशिश कर सकते हो।
उदाहरण के लिए, टूल documentation कहता है कि हम `-c` के साथ character सेट कर सकते हैं।

```bash
cowpy "Hello Containers" -c tux
```

??? success "कमांड आउटपुट"

    ```console
    __________________
    < Hello Containers >
    ------------------
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

इस बार ASCII art आउटपुट Linux penguin, Tux को दिखाता है, क्योंकि हमने `-c tux` पैरामीटर निर्दिष्ट किया।

चूंकि तुम container के अंदर हो, तुम cowpy कमांड को जितनी बार चाहो चला सकते हो, इनपुट पैरामीटर को बदलते हुए, बिना अपने सिस्टम पर किसी भी लाइब्रेरी को इंस्टॉल करने की चिंता किए।

??? tip "अन्य उपलब्ध characters"

    एक अलग character चुनने के लिए '-c' फ्लैग का उपयोग करो, जिसमें शामिल हैं:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

इसके साथ खेलने के लिए स्वतंत्र महसूस करो।
जब तुम कर लो, तो `exit` कमांड का उपयोग करके container से बाहर निकलो:

```bash
exit
```

तुम अपने सामान्य shell में वापस आ जाओगे।

### 4.2. Workflow में एक container का उपयोग करना

जब हम एक पाइपलाइन चलाते हैं, तो हम Nextflow को बताना चाहते हैं कि प्रत्येक चरण में किस container का उपयोग करना है, और महत्वपूर्ण रूप से, हम चाहते हैं कि यह वह सारा काम संभाले जो हमने अभी किया: container को pull करो, इसे spin up करो, कमांड चलाओ और जब यह हो जाए तो container को tear down करो।

अच्छी खबर: यही वह है जो Nextflow हमारे लिए करने जा रहा है।
हमें बस प्रत्येक प्रोसेस के लिए एक container निर्दिष्ट करने की आवश्यकता है।

यह प्रदर्शित करने के लिए कि यह कैसे काम करता है, हमने अपने workflow का एक और संस्करण बनाया जो तीसरे चरण में उत्पादित एकत्रित अभिवादनों की फ़ाइल पर `cowpy` चलाता है।

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

यह speech bubble में तीन अभिवादनों के साथ ASCII art युक्त एक फ़ाइल आउटपुट करना चाहिए।

#### 4.2.1. कोड की जांच करना

Workflow पिछले वाले के समान है, साथ ही `cowpy` चलाने के लिए अतिरिक्त चरण।

??? full-code "पूर्ण कोड फ़ाइल"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // modules शामिल करें
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // इनपुट के लिए CSV फ़ाइल से एक channel बनाएं
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // एक अभिवादन emit करें
        sayHello(greeting_ch)
        // अभिवादन को uppercase में बदलें
        convertToUpper(sayHello.out)
        // सभी अभिवादनों को एक फ़ाइल में एकत्र करें
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // cowpy के साथ अभिवादनों का ASCII art उत्पन्न करें
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

तुम देखते हो कि यह workflow एक module फ़ाइल से एक `cowpy` प्रोसेस import करता है, और इसे `collectGreetings()` कॉल के आउटपुट पर कॉल करता है, साथ ही `params.character` नामक एक इनपुट पैरामीटर।

```groovy title="2d-container.nf" linenums="31"
// cowpy के साथ अभिवादनों का ASCII art उत्पन्न करें
cowpy(collectGreetings.out.outfile, params.character)
```

`cowpy` प्रोसेस, जो ASCII art उत्पन्न करने के लिए cowpy कमांड को wrap करता है, `cowpy.nf` module में परिभाषित है।

??? full-code "पूर्ण कोड फ़ाइल"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // cowpy के साथ ASCII art उत्पन्न करें (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

`cowpy` प्रोसेस को दो इनपुट की आवश्यकता है: speech bubble में डालने के लिए टेक्स्ट युक्त एक इनपुट फ़ाइल का पथ (`input_file`), और character variable के लिए एक मान।

महत्वपूर्ण रूप से, इसमें `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'` लाइन भी शामिल है, जो उस container URI की ओर इशारा करती है जिसका हमने पहले उपयोग किया था।

#### 4.2.2. जांचें कि कॉन्फ़िगरेशन में Docker सक्षम है

हम `nextflow.config` कॉन्फ़िगरेशन फ़ाइल को पेश करके प्रशिक्षण के भाग 3 को थोड़ा anticipate करने जा रहे हैं, जो workflow निष्पादन को कॉन्फ़िगर करने के लिए Nextflow द्वारा प्रदान किए जाने वाले मुख्य तरीकों में से एक है।
जब वर्तमान डायरेक्टरी में `nextflow.config` नामक एक फ़ाइल मौजूद होती है, तो Nextflow स्वचालित रूप से इसे लोड करेगा और इसमें निहित किसी भी कॉन्फ़िगरेशन को लागू करेगा।

इसके लिए, हमने एक `nextflow.config` फ़ाइल शामिल की जो Docker को सक्षम करने वाली कोड की एक ही लाइन के साथ है।

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

यह कॉन्फ़िगरेशन Nextflow को किसी भी प्रोसेस के लिए Docker का उपयोग करने के लिए कहता है जो एक compatible container निर्दिष्ट करता है।

!!! tip

    तकनीकी रूप से command-line से, प्रति-रन आधार पर, `-with-docker <container>` पैरामीटर का उपयोग करके Docker निष्पादन को सक्षम करना संभव है।
    हालांकि, यह हमें केवल पूरे workflow के लिए एक container निर्दिष्ट करने की अनुमति देता है, जबकि हमने अभी जो दृष्टिकोण दिखाया है वह हमें प्रति प्रोसेस एक अलग container निर्दिष्ट करने की अनुमति देता है।
    बाद वाला मॉड्यूलरिटी, कोड रखरखाव और reproducibility के लिए बहुत बेहतर है।

#### 4.2.3. Workflow चलाना

बस recap करने के लिए, यह वह है जो हम चलाने वाले हैं:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

क्या तुम्हें लगता है कि यह काम करेगा?

चलो `-resume` फ्लैग के साथ workflow चलाते हैं, और निर्दिष्ट करते हैं कि हम character को turkey चाहते हैं।

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

पहले तीन चरण cached हुए क्योंकि हमने उन्हें पहले ही चला लिया था, लेकिन `cowpy` प्रोसेस नया है इसलिए वह वास्तव में चलता है।

तुम `results` डायरेक्टरी में `cowpy` चरण का आउटपुट पा सकते हो।

??? abstract "फ़ाइल सामग्री"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
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

तुम देखते हो कि character सभी अभिवादन कह रहा है, क्योंकि यह एकत्रित uppercased अभिवादनों की फ़ाइल पर चला।

अधिक महत्वपूर्ण बात, हम cowpy और इसकी सभी निर्भरताओं की उचित स्थापना किए बिना अपनी पाइपलाइन के हिस्से के रूप में इसे चलाने में सक्षम थे।
और अब हम सहयोगियों के साथ पाइपलाइन साझा कर सकते हैं और उन्हें अपने infrastructure पर चलाने के लिए कह सकते हैं बिना उन्हें Docker या इसके विकल्पों में से एक (जैसे Singularity/Apptainer) के अलावा कुछ भी इंस्टॉल करने की आवश्यकता के।

#### 4.2.4. जांचें कि Nextflow ने containerized task को कैसे लॉन्च किया

इस खंड के अंतिम coda के रूप में, चलो `cowpy` प्रोसेस कॉल में से एक के लिए work subdirectory के अंदर देखते हैं ताकि यह थोड़ा और अंतर्दृष्टि प्राप्त हो सके कि Nextflow containers के साथ कैसे काम करता है।

`cowpy` प्रोसेस के लिए work subdirectory का पथ खोजने के लिए अपने `nextflow run` कमांड से आउटपुट की जांच करो।
ऊपर दिखाए गए रन के लिए हमें जो मिला उसे देखते हुए, `cowpy` प्रोसेस के लिए console log line `[7f/caf718]` से शुरू होती है।
यह निम्नलिखित truncated डायरेक्टरी पथ से मेल खाती है: `work/7f/caf718`।

उस डायरेक्टरी में, तुम्हें `.command.run` फ़ाइल मिलेगी जिसमें सभी कमांड हैं जो Nextflow ने पाइपलाइन को निष्पादित करने के दौरान तुम्हारी ओर से चलाए।

??? abstract "फ़ाइल सामग्री"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

यदि तुम इस फ़ाइल में `nxf_launch` खोजते हो, तो तुम्हें कुछ ऐसा दिखाई देना चाहिए:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

यह launch कमांड दिखाता है कि Nextflow प्रोसेस कॉल को लॉन्च करने के लिए बहुत समान `docker run` कमांड का उपयोग कर रहा है जैसा हमने किया जब हमने इसे manually चलाया।
यह संबंधित work subdirectory को container में mount भी करता है, container के अंदर working डायरेक्टरी को तदनुसार सेट करता है, और `.command.sh` फ़ाइल में हमारी templated bash script चलाता है।

यह पुष्टि करता है कि पिछले खंड में हमें manually करना पड़ा सारा कठिन काम अब Nextflow द्वारा हमारे लिए किया जाता है!

### सारांश

तुम समझते हो कि सॉफ़्टवेयर टूल संस्करणों के प्रबंधन और reproducibility सुनिश्चित करने में containers क्या भूमिका निभाते हैं।

अधिक सामान्य रूप से, तुम्हें बुनियादी समझ है कि वास्तविक दुनिया की Nextflow पाइपलाइनों के मुख्य घटक क्या हैं और वे कैसे व्यवस्थित हैं।
तुम जानते हो कि Nextflow कई इनपुट को कुशलता से कैसे प्रोसेस कर सकता है, एक साथ जुड़े कई चरणों से बने workflows चला सकता है, modular कोड घटकों का लाभ उठा सकता है, और अधिक reproducibility और portability के लिए containers का उपयोग कर सकता है।

### आगे क्या है?

एक और ब्रेक लो! यह Nextflow पाइपलाइनों के काम करने के बारे में जानकारी का एक बड़ा ढेर था।

इस प्रशिक्षण के अंतिम खंड में, हम कॉन्फ़िगरेशन के विषय में गहराई से जाने वाले हैं।
तुम सीखोगे कि अपनी पाइपलाइन के निष्पादन को अपने infrastructure के अनुरूप कैसे कॉन्फ़िगर करें साथ ही इनपुट और पैरामीटर के कॉन्फ़िगरेशन को कैसे प्रबंधित करें।

---

## क्विज़

<quiz>
Nextflow प्रत्येक प्रोसेस कॉल के लिए एक अलग task डायरेक्टरी क्यों बनाता है?
- [ ] निष्पादन गति में सुधार के लिए
- [ ] मेमोरी उपयोग को कम करने के लिए
- [x] निष्पादन को अलग करने और आउटपुट के बीच टकराव से बचने के लिए
- [ ] समानांतर फ़ाइल संपीड़न को सक्षम करने के लिए

अधिक जानें: [1.3. मूल आउटपुट और लॉग खोजना](#13-find-the-original-outputs-and-logs)
</quiz>

<quiz>
Workflow चलाते समय `-ansi-log false` विकल्प क्या करता है?
- [ ] सभी console आउटपुट को अक्षम करता है
- [x] आउटपुट से रंग हटाता है
- [x] उन्हें एक लाइन पर संघनित करने के बजाय सभी task डायरेक्टरी पथ दिखाता है
- [ ] verbose debugging मोड सक्षम करता है

अधिक जानें: [1.3.2. टर्मिनल को अधिक विवरण दिखाना](#132-make-the-terminal-show-more-details)

तुम निम्नलिखित environment variables में से किसी एक का भी उपयोग कर सकते हो यदि तुम इस शैली को पसंद करते हो:

```bash
export NXF_ANSI_LOG=0
# या
export NO_COLOR=1
```

</quiz>

<quiz>
कोड `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }` में, `#!groovy .map { line -> line[0] }` क्या करता है?
- [ ] खाली लाइनों को फ़िल्टर करता है
- [ ] लाइनों को वर्णानुक्रम में सॉर्ट करता है
- [x] प्रत्येक CSV पंक्ति से पहला कॉलम निकालता है
- [ ] लाइनों की संख्या गिनता है

अधिक जानें: [1.4.1. CSV से इनपुट डेटा लोड करना](#141-loading-the-input-data-from-the-csv)
</quiz>

<quiz>
आउटपुट फ़ाइल नामों में इनपुट मान शामिल करना क्यों महत्वपूर्ण है (जैसे, `#!groovy "${greeting}-output.txt"`)?
- [ ] प्रोसेसिंग गति में सुधार के लिए
- [ ] resume कार्यक्षमता को सक्षम करने के लिए
- [x] कई इनपुट को प्रोसेस करते समय आउटपुट फ़ाइलों को एक दूसरे को overwrite करने से रोकने के लिए
- [ ] फ़ाइलों को संपीड़ित करना आसान बनाने के लिए

अधिक जानें: [1.4.3. आउटपुट कैसे नामित किए जाते हैं](#143-how-the-outputs-are-named)
</quiz>

<quiz>
मॉड्यूलराइज़्ड workflow में `include` statement का उद्देश्य क्या है?
- [ ] प्रोसेस कोड को workflow फ़ाइल में कॉपी करना
- [x] एक बाहरी module फ़ाइल से एक प्रोसेस परिभाषा import करना
- [ ] कॉन्फ़िगरेशन सेटिंग्स शामिल करना
- [ ] documentation comments जोड़ना

अधिक जानें: [3. मॉड्यूलराइज़्ड पाइपलाइन चलाना](#3-running-modularized-pipelines)
</quiz>

<quiz>
जब तुम एक workflow को modularize करते हो और इसे `-resume` के साथ चलाते हो, तो क्या होता है?
- [ ] modular प्रोसेस के लिए caching अक्षम है
- [ ] सभी tasks को फिर से निष्पादित किया जाना चाहिए
- [x] Caching उत्पन्न job scripts के आधार पर सामान्य रूप से काम करती है
- [ ] केवल मुख्य workflow फ़ाइल cached है

अधिक जानें: [3.2. Workflow चलाना](#32-run-the-workflow)
</quiz>

<quiz>
प्रोसेस परिभाषा में `container` directive क्या निर्दिष्ट करता है?
- [ ] प्रोसेस के लिए working डायरेक्टरी
- [ ] अधिकतम मेमोरी आवंटन
- [x] प्रोसेस चलाने के लिए उपयोग करने के लिए container image URI
- [ ] आउटपुट फ़ाइल प्रारूप

अधिक जानें: [4.2. Workflow में एक container का उपयोग करना](#42-use-a-container-in-a-workflow)
</quiz>

<quiz>
`.command.run` फ़ाइल में, `nxf_launch` फ़ंक्शन में क्या है?
- [ ] Nextflow संस्करण जानकारी
- [ ] Workflow पैरामीटर
- [x] Volume mounts और container सेटिंग्स के साथ `docker run` कमांड
- [ ] प्रोसेस इनपुट घोषणाएं

अधिक जानें: [4.2.4. जांचें कि Nextflow ने containerized task को कैसे लॉन्च किया](#424-inspect-how-nextflow-launched-the-containerized-task)
</quiz>

<quiz>
Nextflow एक containerized प्रोसेस चलाते समय स्वचालित रूप से क्या संभालता है? (सभी लागू का चयन करें)
- [x] यदि आवश्यक हो तो container image को pull करना
- [x] Work डायरेक्टरी को container में mount करना
- [x] Container के अंदर प्रोसेस script चलाना
- [x] निष्पादन के बाद container instance को साफ करना

अधिक जानें: [4. Containerized सॉफ़्टवेयर का उपयोग करना](#4-using-containerized-software)
</quiz>
