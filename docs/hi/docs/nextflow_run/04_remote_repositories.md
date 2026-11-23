# भाग 4: रिमोट पाइपलाइन चलाना

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI-सहायता प्राप्त अनुवाद - [अधिक जानें और सुधार सुझाएं](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

अब तक, तुमने लोकल रूप से स्टोर किए गए वर्कफ़्लो स्क्रिप्ट चलाए हैं।
व्यवहार में, तुम अक्सर रिमोट रिपॉजिटरी में प्रकाशित पाइपलाइन चलाना चाहोगे, जैसे GitHub पर, बिना उन्हें खुद डाउनलोड किए।

Nextflow इसे सरल बनाता है: तुम किसी भी पाइपलाइन को सीधे Git रिपॉजिटरी URL से चला सकते हो।

---

## 1. GitHub से पाइपलाइन चलाना

रिमोट पाइपलाइन चलाने का बुनियादी सिंटैक्स है `nextflow run <repository>`, जहाँ `<repository>` एक GitHub रिपॉजिटरी पाथ हो सकता है जैसे `nextflow-io/hello`, एक पूरा URL, या GitLab, Bitbucket, या किसी अन्य Git होस्टिंग सेवा का पाथ।

### 1.1. पाइपलाइन लॉन्च करना

आधिकारिक Nextflow "hello" डेमो पाइपलाइन चलाओ।
यह इस कोर्स में तुम जो पाइपलाइन चला रहे हो उससे अलग, बहुत सरल पाइपलाइन है: यह इस प्रशिक्षण में उपयोग की गई "Hello" पाइपलाइन से पहले की है, और बस कुछ हार्डकोड की गई भाषाओं में अभिवादन प्रिंट करती है, इसलिए CSV इनपुट या ASCII आर्ट की उम्मीद मत करो।

```bash
nextflow run nextflow-io/hello
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. पाइपलाइन कहाँ cache होती है यह पता करना

पहली बार जब तुम कोई रिमोट पाइपलाइन चलाते हो, Nextflow उसे डाउनलोड करके लोकल रूप से cache करता है।
बाद के रन cache किए गए वर्शन का उपयोग करते हैं जब तक तुम स्पष्ट रूप से अपडेट का अनुरोध नहीं करते।

डिफ़ॉल्ट रूप से, Nextflow pull की गई पाइपलाइन को `$NXF_HOME/assets` के अंतर्गत सेव करता है।
यह पता करने के लिए कि कोई विशेष पाइपलाइन कहाँ गई, और कौन से revisions उपलब्ध हैं, Nextflow से सीधे पूछो:

```bash
nextflow info nextflow-io/hello
```

??? success "कमांड आउटपुट"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow हर उस revision को `>` से चिह्नित करता है जिसे तुमने पहले से लोकल रूप से checkout किया है; बाकी उपलब्ध हैं लेकिन अभी तक working copy में fetch नहीं किए गए।

तुम `nextflow list` से अब तक pull की गई हर पाइपलाइन की सूची भी देख सकते हो:

```bash
nextflow list
```

??? success "कमांड आउटपुट"

    ```console
    nextflow-io/hello
    ```

[Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) कोर्स इस caching तंत्र को और गहराई से कवर करता है, जिसमें pull की गई पाइपलाइन का सोर्स कोड कैसे देखें यह भी शामिल है।

### सारांश

तुम जानते हो कि GitHub रिपॉजिटरी से सीधे पाइपलाइन कैसे चलाएं बिना खुद डाउनलोड किए, और बाद में इसे लोकल रूप से कहाँ खोजें।

### आगे क्या है?

जानो कि reproducibility के लिए रिमोट पाइपलाइन का एक विशेष वर्शन कैसे pin करें।

---

## 2. Reproducibility के लिए वर्शन निर्दिष्ट करना

डिफ़ॉल्ट रूप से, Nextflow default branch से नवीनतम revision चलाता है।
तुम `-r` फ़्लैग का उपयोग करके एक विशेष वर्शन (tag), branch, या commit pin कर सकते हो।

### 2.1. एक विशेष revision pin करना

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow पहली बार इस revision को fetch करता है जब तुम इसका अनुरोध करते हो, इसलिए `Pulling` और `downloaded from` लाइनें दिखती हैं; बाद में उसी revision का अनुरोध करने पर सीधे `Launching` पर जाता है।
एक सटीक revision pin करना reproducibility के लिए आवश्यक है।
यह गारंटी देता है कि तुम और तुम्हारे सहयोगी बिल्कुल एक ही पाइपलाइन कोड चलाते हैं, चाहे रिपॉजिटरी में तब से कुछ भी बदल गया हो।

### 2.2. Revisions केवल प्रति invocation लागू होते हैं

`-r` के साथ revision pin करना केवल उस रन को प्रभावित करता है जहाँ तुम इसे निर्दिष्ट करते हो: यह नहीं बदलता कि बाद का, सामान्य `nextflow run` क्या उपयोग करता है।
पाइपलाइन को `-r` के बिना फिर से चलाकर देखो:

```bash
nextflow run nextflow-io/hello
```

??? success "कमांड आउटपुट"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

भले ही पिछले रन ने स्पष्ट रूप से `v1.3` pin किया था, यह रन सीधे default branch (`master`) पर वापस जाता है।
Nextflow हर उस revision के लिए एक अलग लोकल working copy रखता है जिसे तुमने उपयोग किया है, जो `nextflow info` में `>` मार्कर दिखाता है, लेकिन यह कभी याद नहीं रखता कि तुमने आखिरी बार कौन सा चलाया था।
तुम `nextflow info <pipeline>` चलाकर पाइपलाइन की default branch का नाम पता कर सकते हो; यह वही है जो `(default)` से चिह्नित है।
Reproducibility पूरी तरह तुम पर निर्भर है: जब भी यह मायने रखे, हमेशा `-r` स्पष्ट रूप से पास करो, बजाय यह मानने के कि पिछले रन में pin किया गया revision अभी भी लागू होता है।

### सारांश

तुम जानते हो कि reproducible execution के लिए रिमोट पाइपलाइन को एक विशेष वर्शन, branch, या commit पर कैसे pin करें, और यह कि pin केवल उस एक invocation पर लागू होता है, बाद के रन पर नहीं।

### आगे क्या है?

तुमने Nextflow पाइपलाइन चलाने और प्रबंधित करने की बुनियादी बातें सीख ली हैं।
आगे कहाँ जाना है इसके लिए [कोर्स सारांश](next_steps.md) देखो।

---

## सारांश

इस भाग में तुमने सीखा:

- GitHub रिपॉजिटरी से सीधे पाइपलाइन चलाना बिना उसे डाउनलोड किए
- Reproducibility के लिए रिमोट पाइपलाइन को एक विशेष revision पर pin करना, और यह समझना कि pin केवल उस एक invocation पर लागू होता है
